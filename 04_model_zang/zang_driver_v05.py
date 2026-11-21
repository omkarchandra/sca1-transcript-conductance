# Zang, Dieudonne & De Schutter 2018 Purkinje cell (ModelDB 243446) driver P4 v05.
# v05 = v04 (a363c9d3...) with the beta4 shift inverted to the plan's definition (pre-specified):
#   r = Kcnmb4 ratio = fraction of abBK that remains; in each abBK section gabkbar x r, gbar_mslo + (1 - r) x gabkbar.
# v04 = v03 (e9434afd...) plus an optional "aux" object in the factor file (Zang P9, review 26 Oct ), applied
#   after the five-gene scaling to every segment of every section carrying the mechanism; RECORD.aux_applied lists
#   each value. Allowed keys only (AUX): newCaP.pcabar x, Kv4.vshift +, Kv4s.vshift + (vshift enters hinf only:
#   V1/2 = vh - vshift), and beta4_shift_r (definition as in the v05 line above).
#   Without aux, v04 does exactly what v03 does.
# v03 = v02 (961524e5...) with the soma spike level moved from 0 mV to -30 mV: in this model the somatic spike
#   peaks lie between about -12 and -1 mV (AIS +4 to +12 mV), so 0 mV missed all but the start-up spike.
# v02 = v01 (a81bcb95...) plus --hold-nA: a DC current added to the OU noise in the same somaA(0.5) IClamp
#   (pre-specified), written to RECORD/SMOKE as stimulus.hold_nA. The noise series and its sha are unchanged.
# Set-up: the published purkinje_demo.hoc block from after load_file("nrngui.hoc") up to xopen("dendv_arnd21.ses")
# is executed byte-for-byte (stdrun.hoc replaces nrngui.hoc). The 45-s protocol replaces the published stimulation:
# dendv_arnd21.ses (graphs), electrode.hoc (stim1, amp 0 -> cfdistrib.hoc: 500 syn2 climbing-fibre synapses at
# 400 ms) and distri.hoc (shape plot) are not loaded. celsius 34 and v_init -65-3 as published; dt 0.02 ms.
# Five-gene factors scale every segment of every section carrying the mapped mechanism(s); RECORD lists each value.
# Usage: python -I -B <this file> <MODEL_DIR> <OUT_DIR> <NRNMECH_LIB> <FACTOR_JSON> <RUN_S>
#        [--noise-sigma-nA S --noise-tau-ms T --noise-seed N] [--hold-nA I] [--inspect | --smoke-ms MS]
# Single thread only: the mod files are not THREADSAFE (NEURON: "cdp4Nsp is not thread safe"; multisplit segfaults).
import argparse
import hashlib
import json
import math
import os
import platform
import re
import socket
import sys
import time
from pathlib import Path

import numpy as np

DT_MS = 0.02
REC_DT_MS = 0.1
CELSIUS = 34.0
V_INIT = -68.0
SPIKE_MV = -30.0
REFRACT_MS = 2.0
PEAK_MS = 5.0
DEND_MV = -30.0
DEND_REFRACT_MS = 20.0
DEND_PEAK_MS = 20.0
SOMA = ("somaA", 0.5)
DISTAL = ("dendA1_001011110110010110", 0.5)  # farthest dendritic site plotted in the published dendv_arnd21.ses
MAIN = ("dendA1_0010111101", 0.5)  # farthest smooth (maindend) site plotted there
PINS = {
    "purkinje_demo.hoc": "4890a8f2c12e68d9258d138222e6b51515e9a3738dc50e28603a7276a1052740",
    "Purkinje19b972-1.nrn": "e6ff8733eaae8f1e5312f42c70526463e21c2b72fc5c9feb02d2eb782c94366c",
    "dendv_arnd21.ses": "66d513162e78e9eac92ec2feaf465ca493cd244331b59acd54b22c5528fc5189",
    "electrode.hoc": None,
    "cfdistrib.hoc": None,
}
HEAD_MARK = 'load_file("nrngui.hoc")\n'
SES_MARK = 'xopen("dendv_arnd21.ses")'
PUBLISHED_TAIL = ("\ncelsius = 34\n", "\ndt = 0.02\n", "\nv_init = -65-3\n", "\nstim1.amp =1.06*0\n")
MECHS = {  # gene -> (mechanism, density parameter); mapping, review 26 Oct
    "Kcnc3": (("Kv3", "gbar_Kv3"),),
    "Kcnd3": (("Kv4", "gbar_Kv4"), ("Kv4s", "gbar_Kv4s")),
    "Kcnma1": (("mslo", "gbar_mslo"), ("abBK", "gabkbar_abBK")),
    "Cacna1a": (("newCaP", "pcabar_newCaP"),),
    "Cacna1g": (("CaT3_1", "pcabar_CaT3_1"),),
}
DENSITY = {"Kv3": "gbar_Kv3", "Kv4": "gbar_Kv4", "Kv4s": "gbar_Kv4s", "Kv1": "gbar_Kv1", "mslo": "gbar_mslo",
           "abBK": "gabkbar_abBK", "SK2": "gkbar_SK2", "newCaP": "pcabar_newCaP", "CaT3_1": "pcabar_CaT3_1",
           "naRsg": "gbar_naRsg", "nap": "gbar_nap", "hpkj": "ghbar_hpkj", "pas": "g_pas"}
AUX = {  # key -> (mechanism, attribute, operation); beta4_shift_r is handled on its own
    "newCaP.pcabar": ("newCaP", "pcabar_newCaP", "multiply"),
    "Kv4.vshift": ("Kv4", "vshift_Kv4", "add"),
    "Kv4s.vshift": ("Kv4s", "vshift_Kv4s", "add"),
    "beta4_shift_r": ("abBK", "gabkbar_abBK", "beta4"),
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def crossings(t, v, level, refract):
    idx = np.where((v[:-1] < level) & (v[1:] >= level))[0] + 1
    out, last = [], -1e9
    for i in idx:
        if t[i] - last >= refract:
            out.append(i)
            last = t[i]
    return np.array(out, dtype=int)


def make_noise(sigma, tau_ms, seed, n):
    """Ornstein-Uhlenbeck series in nA, one sample per dt step, x[0] = 0 (Masoli driver v09 generator)."""
    a = math.exp(-DT_MS / tau_ms)
    rng = np.random.default_rng(seed)
    z = rng.standard_normal(n - 1)
    x = np.zeros(n)
    k = sigma * math.sqrt(1.0 - a * a)
    for i in range(n - 1):
        x[i + 1] = x[i] * a + k * z[i]
    return x


def region_of(h, sec, lists):
    name = sec.name()
    if name == "somaA":
        return "soma"
    if name == "AIS":
        return "AIS"
    if name in lists["maindend"]:
        return "maindend"
    if name in lists["spinydend"]:
        return "spinydend"
    return "other"


def seg_values(sec, attr):
    vals = [getattr(seg, attr) for seg in sec]
    return vals[0] if len(set(vals)) == 1 else vals


def apply_factors(h, factors, lists):
    applied, counts = [], {}
    for gene, pairs in MECHS.items():
        gf = float(factors[gene])
        for mech, attr in pairs:
            n = 0
            for sec in h.allsec():
                if not sec.has_membrane(mech):
                    continue
                published = seg_values(sec, attr)
                for seg in sec:
                    setattr(seg, attr, getattr(seg, attr) * gf)
                applied.append({"gene": gene, "mechanism": mech, "attr": attr, "section": sec.name(),
                                "region": region_of(h, sec, lists), "nseg": sec.nseg, "published": published,
                                "factor": gf, "applied": seg_values(sec, attr)})
                n += 1
            counts[f"{gene}|{mech}"] = n
    zero = [k for k, c in counts.items() if c == 0]
    if zero:
        raise SystemExit(f"zero sections carry: {zero}")
    return applied, counts


def apply_aux(h, aux, lists):
    applied, counts = [], {}
    for key, val in aux.items():
        mech, attr, op = AUX[key]
        val, n = float(val), 0
        for sec in h.allsec():
            if not sec.has_membrane(mech):
                continue
            row = {"key": key, "mechanism": mech, "attr": attr, "op": op, "value": val, "section": sec.name(),
                   "region": region_of(h, sec, lists), "nseg": sec.nseg, "before": seg_values(sec, attr)}
            if op == "beta4":
                if not sec.has_membrane("mslo"):
                    raise SystemExit(f"beta4 shift: {sec.name()} carries abBK but not mslo")
                row["mslo_before"] = seg_values(sec, "gbar_mslo")
                for seg in sec:
                    g = seg.gabkbar_abBK
                    seg.gabkbar_abBK = g * val
                    seg.gbar_mslo = seg.gbar_mslo + (1.0 - val) * g
                row["mslo_after"] = seg_values(sec, "gbar_mslo")
            else:
                for seg in sec:
                    x = getattr(seg, attr)
                    setattr(seg, attr, x * val if op == "multiply" else x + val)
            row["after"] = seg_values(sec, attr)
            applied.append(row)
            n += 1
        counts[key] = n
        if n == 0:
            raise SystemExit(f"zero sections carry aux mechanism: {key}")
    return applied, counts


def inspect(h, lists, text_by_file):
    secs = list(h.allsec())
    mechs = {}
    for sec in secs:
        reg = region_of(h, sec, lists)
        for m in sec.psection()["density_mechs"]:
            d = mechs.setdefault(m, {"regions": {}, "values": {}})
            d["regions"][reg] = d["regions"].get(reg, 0) + 1
            if m in DENSITY:
                vals = d["values"].setdefault(reg, set())
                v = seg_values(sec, DENSITY[m])
                vals.update(v if isinstance(v, list) else [v])
    for m in mechs:
        mechs[m]["values"] = {r: sorted(v) for r, v in mechs[m]["values"].items()}
        if m in DENSITY:
            mechs[m]["density_attr"] = DENSITY[m]
    ses = text_by_file["dendv_arnd21.ses"]
    sites = []
    for name, x in re.findall(r"(dendA1_[0-9]+)\.v\(\s*([0-9.]+)\s*\)", ses):
        sec = getattr(h, name)
        sites.append({"section": name, "x": float(x), "region": region_of(h, sec, lists),
                      "path_distance_um_from_somaA1": h.distance(h.somaA(1), sec(float(x))),
                      "diam_um": sec.diam, "mechanisms": sorted(sec.psection()["density_mechs"])})
    drive = {f: [ln.strip() for ln in text_by_file[f].splitlines()
                 if re.search(r"stim1\.|synlist\.object\(i\)\.|NUMSYN =|cfstart =|tstop =", ln)]
             for f in ("purkinje_demo.hoc", "electrode.hoc", "cfdistrib.hoc")}
    return {"sections": len(secs), "segments": sum(s.nseg for s in secs),
            "section_lists": {k: len(v) for k, v in lists.items()},
            "non_dendritic_sections": [s.name() for s in secs if not s.name().startswith("dend")],
            "maindend_sections": sorted(lists["maindend"]),
            "mechanisms": mechs, "published_plot_sites": sites, "published_drive_lines": drive}


def main():
    ap = argparse.ArgumentParser()
    for name in ("model_dir", "out_dir", "lib", "factor_json", "run_s"):
        ap.add_argument(name)
    ap.add_argument("--noise-sigma-nA", type=float, default=0.0)
    ap.add_argument("--noise-tau-ms", type=float, default=5.0)
    ap.add_argument("--noise-seed", type=int, default=None)
    ap.add_argument("--hold-nA", type=float, default=0.0)
    ap.add_argument("--inspect", action="store_true")
    ap.add_argument("--smoke-ms", type=float, default=None)
    a = ap.parse_args()
    model_dir, out_dir, lib, factor_json = Path(a.model_dir), Path(a.out_dir), Path(a.lib), Path(a.factor_json)
    run_s = float(a.run_s)
    if out_dir.exists():
        raise SystemExit(f"output exists; refusing to overwrite: {out_dir}")
    if a.noise_sigma_nA > 0 and a.noise_seed is None:
        raise SystemExit("--noise-seed is required when sigma > 0")
    if (model_dir / "x86_64").exists():
        raise SystemExit("model dir contains x86_64; NEURON would auto-load it")
    for f, pin in PINS.items():
        if pin is not None and sha256(model_dir / f) != pin:
            raise SystemExit(f"published file hash mismatch: {f}")
    text_by_file = {f: (model_dir / f).read_text() for f in PINS}
    demo = text_by_file["purkinje_demo.hoc"]
    if demo.count(HEAD_MARK) != 1 or demo.count(SES_MARK) != 1:
        raise SystemExit("set-up markers not found exactly once in purkinje_demo.hoc")
    i0, i1 = demo.index(HEAD_MARK) + len(HEAD_MARK), demo.index(SES_MARK)
    setup, tail = demo[i0:i1], demo[i1:]
    missing = [s.strip() for s in PUBLISHED_TAIL if s not in tail]
    if missing:
        raise SystemExit(f"published settings not found: {missing}")
    fdoc = json.loads(factor_json.read_text())
    if not (isinstance(fdoc, dict) and isinstance(fdoc.get("scales"), dict)):
        raise SystemExit("factor file must be JSON with a scales object")
    factors = fdoc["scales"]
    lacking = [g for g in MECHS if g not in factors]
    if lacking:
        raise SystemExit(f"factor file lacks mapped genes: {lacking}")
    unmapped = {g: v for g, v in factors.items() if g not in MECHS and float(v) != 1.0}
    if unmapped:
        raise SystemExit(f"unmapped genes must be 1.0: {unmapped}")
    aux = fdoc.get("aux", {})
    if not isinstance(aux, dict) or any(k not in AUX for k in aux):
        raise SystemExit(f"aux keys must be among {sorted(AUX)}: {aux}")
    if "beta4_shift_r" in aux and not 0.0 <= float(aux["beta4_shift_r"]) <= 1.0:
        raise SystemExit("beta4_shift_r must lie in [0, 1]")

    os.chdir(model_dir)
    from neuron import h
    h.nrn_load_dll(str(lib))
    h.load_file("stdrun.hoc")
    t_setup = time.perf_counter()
    if not h(setup):
        raise SystemExit("published set-up block failed in hoc")
    lists = {n: set(s.name() for s in getattr(h, n)) for n in ("alldend", "maindend", "spinydend", "cf", "bs")}
    h.celsius = CELSIUS
    h.dt = DT_MS
    h.steps_per_ms = 1.0 / DT_MS
    h.v_init = V_INIT
    h.CVode().active(0)

    if a.inspect:
        doc = inspect(h, lists, text_by_file)
        doc.update({"celsius": h.celsius, "dt_ms": h.dt, "v_init": h.v_init, "published_dt_ms": 0.02,
                    "pins": {f: sha256(model_dir / f) for f in PINS}, "nrnmech_lib_sha256": sha256(lib),
                    "neuron_version": h.nrnversion(0), "hostname": socket.gethostname()})
        out_dir.mkdir()
        with open(out_dir / "INSPECT.json", "x", encoding="utf-8", newline="") as f:
            f.write(json.dumps(doc, indent=1) + "\n")
        print("inspect ok", doc["sections"], doc["segments"])
        return

    applied, counts = apply_factors(h, factors, lists)
    aux_applied, aux_counts = apply_aux(h, aux, lists)
    tstop_ms = a.smoke_ms if a.smoke_ms is not None else run_s * 1000.0
    noise_info = {"sigma_nA": a.noise_sigma_nA, "tau_ms": a.noise_tau_ms, "seed": a.noise_seed, "sha256": None,
                  "site": f"{SOMA[0]}({SOMA[1]})"}
    stimulus = {"site": f"{SOMA[0]}({SOMA[1]})", "hold_nA": a.hold_nA, "noise_sigma_nA": a.noise_sigma_nA,
                "amp": "hold_nA + OU noise sample per dt step" if a.noise_sigma_nA > 0 else "hold_nA constant"}
    if a.noise_sigma_nA > 0 or a.hold_nA != 0.0:
        stim = h.IClamp(getattr(h, SOMA[0])(SOMA[1]))
        stim.delay = 0
        stim.dur = 1e9
        stim.amp = a.hold_nA
    if a.noise_sigma_nA > 0:
        noise = make_noise(a.noise_sigma_nA, a.noise_tau_ms, a.noise_seed, int(round(tstop_ms / DT_MS)) + 1)
        noise_info.update({"sha256": hashlib.sha256(noise.tobytes()).hexdigest(), "length": int(len(noise)),
                           "mean": float(noise.mean()), "sd": float(noise.std())})
        amp_vec = h.Vector(noise + a.hold_nA)
        amp_vec.play(stim, stim._ref_amp, DT_MS)
    threads = int(h.ParallelContext().nthread())
    if threads != 1:
        raise SystemExit(f"expected one thread, got {threads}")
    recs = {}
    for key, (name, x) in (("soma", SOMA), ("distal", DISTAL), ("main", MAIN)):
        recs[key] = h.Vector()
        recs[key].record(getattr(h, name)(x)._ref_v, REC_DT_MS)
    setup_s = time.perf_counter() - t_setup
    t_run = time.perf_counter()
    h.finitialize(V_INIT)
    h.continuerun(tstop_ms)
    wall_s = time.perf_counter() - t_run

    v = {k: np.array(r) for k, r in recs.items()}
    t = np.arange(len(v["soma"])) * REC_DT_MS
    i_s = crossings(t, v["soma"], SPIKE_MV, REFRACT_MS)
    i_d = crossings(t, v["distal"], DEND_MV, DEND_REFRACT_MS)
    i_m = crossings(t, v["main"], DEND_MV, DEND_REFRACT_MS)
    n_pk, n_dpk = int(PEAK_MS / REC_DT_MS), int(DEND_PEAK_MS / REC_DT_MS)
    out_dir.mkdir()
    if a.smoke_ms is not None:
        smoke = {"stage": "Zang 2018 driver v05 SMOKE", "smoke_ms": tstop_ms, "threads": threads,
                 "aux": aux, "aux_section_counts": aux_counts, "aux_applied": aux_applied,
                 "stimulus": stimulus,
                 "wall_s": wall_s, "setup_s": setup_s, "factor_file_sha256": sha256(factor_json),
                 "noise": noise_info, "soma_spike_times_ms": [float(x) for x in t[i_s]],
                 "distal_crossings_ms": [float(x) for x in t[i_d]], "main_crossings_ms": [float(x) for x in t[i_m]],
                 "v_range": {k: [float(x.min()), float(x.max())] for k, x in v.items()},
                 "all_finite": bool(all(np.all(np.isfinite(x)) for x in v.values())),
                 "nrnmech_lib_sha256": sha256(lib), "neuron_version": h.nrnversion(0),
                 "hostname": socket.gethostname()}
        with open(out_dir / "SMOKE.json", "x", encoding="utf-8", newline="") as f:
            f.write(json.dumps(smoke, indent=1) + "\n")
        print("smoke ok", f"wall_s={wall_s:.1f}", f"spikes={len(i_s)}")
        return
    with open(out_dir / "SPIKES.csv", "x", encoding="utf-8", newline="") as f:
        f.write("time_ms,peak_mV\n")
        f.write("".join(f"{t[i]:.4f},{v['soma'][i:i + n_pk].max():.4f}\n" for i in i_s))
    with open(out_dir / "DSPIKES.csv", "x", encoding="utf-8", newline="") as f:
        f.write("time_ms,peak_mV\n")
        f.write("".join(f"{t[i]:.4f},{v['distal'][i:i + n_dpk].max():.4f}\n" for i in i_d))
    names = ["SPIKES.csv", "DSPIKES.csv"]
    for key, fname in (("soma", "V_SOMA.npy"), ("distal", "V_DISTAL.npy"), ("main", "V_MAIN.npy")):
        with open(out_dir / fname, "xb") as f:
            np.save(f, v[key].astype(np.float32))
        names.append(fname)
    record = {
        "stage": "Zang 2018 (ModelDB 243446) P4 driver v05",
        "condition": fdoc.get("condition"),
        "factor_file": str(factor_json),
        "factor_file_sha256": sha256(factor_json),
        "scales": {g: float(factors[g]) for g in MECHS},
        "mapping": {g: [list(p) for p in pairs] for g, pairs in MECHS.items()},
        "mechanism_section_counts": counts,
        "run_s": run_s, "tstop_ms": tstop_ms, "dt_ms": DT_MS, "published_dt_ms": 0.02,
        "record_dt_ms": REC_DT_MS, "celsius": h.celsius, "v_init": V_INIT,
        "time_axis": "index x record_dt_ms from 0",
        "n_samples": int(len(t)),
        "setup": "purkinje_demo.hoc between load_file(\"nrngui.hoc\") and xopen(\"dendv_arnd21.ses\"), "
                 "executed byte-for-byte after stdrun.hoc",
        "omitted_published": ["dendv_arnd21.ses (graphs)", "electrode.hoc (stim1 IClamp, amp 1.06*0 in demo)",
                              "cfdistrib.hoc (500 syn2 climbing-fibre synapses, onset 400 ms)",
                              "distri.hoc (shape plot)"],
        "pins": {f: sha256(model_dir / f) for f in PINS},
        "nrnmech_lib_sha256": sha256(lib),
        "neuron_version": h.nrnversion(0),
        "python": sys.version.split()[0], "numpy": np.__version__, "platform": platform.platform(),
        "hostname": socket.gethostname(),
        "task": {k: os.environ.get(k) for k in ("JOB_ID", "TASK_ID")},
        "threads": threads,
        "noise": noise_info,
        "stimulus": stimulus,
        "sites": {"soma": f"{SOMA[0]}({SOMA[1]})", "distal": f"{DISTAL[0]}({DISTAL[1]})",
                  "main": f"{MAIN[0]}({MAIN[1]})"},
        "soma_spike_definition": {"level_mV": SPIKE_MV, "refractory_ms": REFRACT_MS, "peak_window_ms": PEAK_MS},
        "dendritic_spike_definition": {"site": "distal", "level_mV": DEND_MV, "direction": "upward",
                                       "refractory_ms": DEND_REFRACT_MS, "peak_window_ms": DEND_PEAK_MS},
        "n_soma_spikes": int(len(i_s)), "n_distal_dspikes": int(len(i_d)), "n_main_crossings": int(len(i_m)),
        "wall_s": wall_s, "setup_s": setup_s,
        "applied": applied,
        "aux": aux, "aux_section_counts": aux_counts, "aux_applied": aux_applied,
        "saved": names,
    }
    with open(out_dir / "RECORD.json", "x", encoding="utf-8", newline="") as f:
        f.write(json.dumps(record, indent=1) + "\n")
    names.append("RECORD.json")
    with open(out_dir / "SHA256SUMS", "x", encoding="utf-8", newline="") as f:
        f.write("".join(f"{sha256(out_dir / n)}  {n}\n" for n in names))
    print("done", f"spikes={len(i_s)}", f"dspikes={len(i_d)}", f"wall_s={wall_s:.0f}")


if __name__ == "__main__":
    main()

#%% debug
