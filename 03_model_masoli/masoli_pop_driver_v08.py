# Masoli 2015 (ModelDB 229585) population driver S2 DRAFT.
# v08 = v07 (c1d9bc7c...) plus an optional "cdp5_scales" object: per-parameter multipliers
#   for Buffnull2_cdp5, CBnull_cdp5 and PVnull_cdp5. "cdp5_scale" (one multiplier for all three)
#   still works. Five-gene behaviour is unchanged.
# Derived from driver v05 (d1fb515c...). Changes: member scale x gene factor applied to the listed
# mechanisms; run length is a parameter; V_TRACE only with --save-v; RECORD lists every applied value.
# Usage: python -I -B <this file> <MODEL_DIR> <OUT_DIR> <NRNMECH_LIB> <MEMBER_ID> <MEMBERS_TSV> <FACTOR_JSON> <RUN_S> [--save-v]
import argparse
import csv
import hashlib
import json
import math
import os
import platform
import sys
from pathlib import Path

import numpy as np

DT_MS = 0.025
REC_DT_MS = 0.1
THREADS = 8
SPIKE_MV = 0.0
REPORT_MV = -20.0
REFRACT_MS = 2.0
WINDOW_MS = 2000.0
TAIL_MS = 10000.0
PEAK_MS = 5.0
MODEL_HEAD = "7b1b3b99d5719bf320d8d249fbcb83a3da9e36dd"
MECH_OF_GENE = {  # gene -> (population group, RANGE prefix, mech suffix)
    "Kcnc3": ("Kv3_3", "gbar", "Kv3_3", ["dend"]),
    "Kcnma1": ("Kca1_1", "gbar", "Kca1_1", ["soma", "dend"]),
    "Kcnd3": ("Kv4_3", "gkbar", "Kv4_3", ["dend"]),
    "Cacna1a": ("Cav2_1", "pcabar", "Cav2_1", ["soma", "dend", "AIS"]),
    "Cacna1g": ("Cav3_1", "pcabar", "Cav3_1", ["soma", "dend", "AIS"]),
    "Scn8a": ("Nav1_6", "gbar", "Nav1_6", ["soma", "dend", "AIS"]),
    "Hcn1": ("HCN1", "gbar", "HCN1", ["soma", "dend"]),
    "Kcnn2": ("Kca2_2", "gkbar", "Kca2_2", ["soma", "dend"]),
    "Kcna1": (None, "gbar", "Kv1_1", ["soma", "dend"]),  # held at 1.0 member scale
}
NO_MECH = {"Scn4b", "Fgf14"}
FIVE = ("Kcnc3", "Kcnma1", "Kcnd3", "Cacna1a", "Cacna1g")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def crossings(t, v, level):
    idx = np.where((v[:-1] < level) & (v[1:] >= level))[0] + 1
    out, last = [], -1e9
    for i in idx:
        if t[i] - last >= REFRACT_MS:
            out.append(i)
            last = t[i]
    return np.array(out, dtype=int)


def load_member(tsv, member_id):
    with open(tsv, encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    row = [r for r in rows if int(r["member_id"]) == member_id]
    if len(row) != 1:
        raise SystemExit(f"member {member_id} not found once")
    return {k: float(v) for k, v in row[0].items() if k != "member_id"}


def apply_factors(cell, member, factors, applied):
    regions = {"soma": [cell.soma], "dend": list(cell.dend), "AIS": [cell.axonAIS]}
    counts = {}
    for gene in MECH_OF_GENE:
        if gene not in factors:
            raise SystemExit(f"factor file lacks mapped gene: {gene}")
        group, prefix, suffix, where = MECH_OF_GENE[gene]
        gf = factors[gene]
        ms = member[group] if group else 1.0
        for reg in where:
            n = 0
            for sec in regions[reg]:
                if not sec.has_membrane(suffix):
                    continue
                attr = f"{prefix}_{suffix}"
                published = getattr(sec, attr)
                new = published * ms * gf
                setattr(sec, attr, new)
                n += 1
                applied.append({"gene": gene, "section": reg, "name": sec.name(), "attr": attr,
                                "published": published, "member_scale": ms, "gene_factor": gf,
                                "applied": new})
            counts[f"{gene}|{reg}"] = n
    for gene in factors:
        if gene not in MECH_OF_GENE and gene not in NO_MECH:
            raise SystemExit(f"gene key not in mapping: {gene}")
    return counts


def apply_cdp5(cell, scale_doc, applied):
    """Apply one common or per-parameter cdp5 calcium-buffer multipliers."""
    allowed = ("Buffnull2_cdp5", "CBnull_cdp5", "PVnull_cdp5")
    if set(scale_doc) - set(allowed):
        raise SystemExit(f"unknown cdp5 buffer parameter(s): {sorted(set(scale_doc) - set(allowed))}")
    for attr, val in scale_doc.items():
        if not (0.1 <= float(val) <= 10.0):
            raise SystemExit(f"cdp5 scale outside [0.1, 10] for {attr}: {val}")
    n = 0
    for sec in cell.allsec():
        if not sec.has_membrane("cdp5"):
            continue
        for attr in scale_doc:
            published = getattr(sec, attr)
            new = published * float(scale_doc[attr])
            setattr(sec, attr, new)
            applied.append({"mechanism": "cdp5", "section": sec.name(), "attr": attr,
                            "published": published, "buffer_scale": float(scale_doc[attr]),
                            "applied": new})
        n += 1
    if n == 0:
        raise SystemExit("no cdp5 sections found")
    return n


def make_noise(sigma, tau_ms, seed, n):
    """Ornstein-Uhlenbeck series in nA, one sample per dt step, x[0] = 0."""
    a = math.exp(-DT_MS / tau_ms)
    rng = np.random.default_rng(seed)
    z = rng.standard_normal(n - 1)
    x = np.zeros(n)
    k = sigma * math.sqrt(1.0 - a * a)
    for i in range(n - 1):
        x[i + 1] = x[i] * a + k * z[i]
    return x


def main():
    ap = argparse.ArgumentParser()
    for name in ("model_dir", "out_dir", "lib", "member_id", "members_tsv", "factor_json", "run_s"):
        ap.add_argument(name)
    ap.add_argument("--save-v", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--published", action="store_true")
    ap.add_argument("--noise-sigma-nA", type=float, default=0.0)
    ap.add_argument("--noise-tau-ms", type=float, default=5.0)
    ap.add_argument("--noise-seed", type=int, default=None)
    ap.add_argument("--smoke-ms", type=float, default=None)
    a = ap.parse_args()
    save_v, dry_run = a.save_v, a.dry_run
    model_dir, out_dir, lib = Path(a.model_dir), Path(a.out_dir), Path(a.lib)
    member_id, members_tsv, factor_json, run_s = int(a.member_id), Path(a.members_tsv), Path(a.factor_json), float(a.run_s)
    if a.noise_sigma_nA > 0 and a.noise_seed is None:
        raise SystemExit("--noise-seed is required when sigma > 0")
    member = load_member(members_tsv, member_id)
    if a.published:
        member = {k: 1.0 for k in member}
    fdoc = json.loads(factor_json.read_text())
    if not (isinstance(fdoc, dict) and isinstance(fdoc.get("scales"), dict)):
        raise SystemExit("factor file must be JSON with a scales object")
    factors = fdoc["scales"]
    missing = [g for g in FIVE if g not in factors]
    if missing:
        raise SystemExit(f"factor file lacks scaled genes: {missing}")
    cdp5_scale = fdoc.get("cdp5_scale")
    cdp5_scales = fdoc.get("cdp5_scales")
    if cdp5_scale is not None and not isinstance(cdp5_scale, (int, float)):
        raise SystemExit("cdp5_scale must be a number")
    if cdp5_scales is not None:
        if not isinstance(cdp5_scales, dict):
            raise SystemExit("cdp5_scales must be an object")
        if cdp5_scale is not None:
            raise SystemExit("use either cdp5_scale or cdp5_scales, not both")
    buffer_scales = cdp5_scales if cdp5_scales is not None else (
        {k: cdp5_scale for k in ("Buffnull2_cdp5", "CBnull_cdp5", "PVnull_cdp5")}
        if cdp5_scale is not None else None)
    os.chdir(model_dir)
    sys.path.insert(0, str(model_dir))

    from neuron import h
    h.nrn_load_dll(str(lib))
    h.load_file("stdrun.hoc")
    from Purkinje_py3 import Purkinje_py3

    h.dt = DT_MS
    h.celsius = 37
    h.v_init = -65
    h.CVode().active(0)
    cell = Purkinje_py3()
    applied = []
    counts = apply_factors(cell, member, factors, applied)
    cdp5_sections = 0
    if buffer_scales is not None:
        cdp5_sections = apply_cdp5(cell, buffer_scales, applied)
    zero = [k for k, c in counts.items() if c == 0]
    if zero:
        raise SystemExit(f"zero sections changed for: {zero}")
    tstop_ms = run_s * 1000.0
    noise = None
    noise_info = {"sigma_nA": a.noise_sigma_nA, "tau_ms": a.noise_tau_ms, "seed": a.noise_seed, "sha256": None}
    if a.noise_sigma_nA > 0:
        noise = make_noise(a.noise_sigma_nA, a.noise_tau_ms, a.noise_seed, int(round(tstop_ms / DT_MS)) + 1)
        noise_info.update({"sha256": hashlib.sha256(noise.tobytes()).hexdigest(), "length": int(len(noise)),
                           "mean": float(noise.mean()), "sd": float(noise.std())})
    if dry_run:
        out_dir.mkdir(exist_ok=False)
        dry = {"stage": "Masoli 2015 population v07 DRY RUN", "member_id": member_id,
               "members_tsv_sha256": sha256(members_tsv), "member_scales": member,
               "factor_file": str(factor_json), "factor_file_sha256": sha256(factor_json),
               "cdp5_scale": cdp5_scale, "cdp5_scales": cdp5_scales,
               "cdp5_sections": cdp5_sections,
               "region_counts": counts, "applied": applied, "noise": noise_info}
        with open(out_dir / "DRYRUN.json", "x", encoding="utf-8", newline="") as f:
            f.write(json.dumps(dry, indent=1) + "\n")
        print("dry-run ok", f"pairs={len(counts)}", f"applied={len(applied)}")
        return
    stim = None
    if noise is not None:
        stim = h.IClamp(cell.soma(0.5))
        stim.delay = 0
        stim.dur = 1e9
        amp_vec = h.Vector(noise)
        amp_vec.play(stim, stim._ref_amp, DT_MS)
    cores = THREADS
    h.load_file("parcom.hoc")
    p = h.ParallelComputeTool()
    p.change_nthread(cores,1)
    p.multisplit(1)
    vv = h.Vector()
    vv.record(cell.soma(0.5)._ref_v, REC_DT_MS)
    if a.smoke_ms is not None:
        h.tstop = a.smoke_ms
        h.finitialize(-65)
        h.continuerun(a.smoke_ms)
        vals = np.array(vv)
        out_dir.mkdir(exist_ok=False)
        smoke = {"stage": "Masoli 2015 population v06 SMOKE", "smoke_ms": a.smoke_ms,
                 "member_id": member_id, "published": a.published, "sigma_nA": a.noise_sigma_nA,
                 "seed": a.noise_seed, "noise_sha256": noise_info["sha256"],
                 "n_samples": int(len(vals)), "last10_soma_mV": [float(x) for x in vals[-10:]],
                 "all_finite": bool(np.all(np.isfinite(vals)))}
        with open(out_dir / "SMOKE.json", "x", encoding="utf-8", newline="") as f:
            f.write(json.dumps(smoke, indent=1) + "\n")
        print("smoke ok", f"n={len(vals)}", f"finite={smoke['all_finite']}")
        return
    h.tstop = tstop_ms
    h.finitialize(-65)
    h.run()
    out_dir.mkdir(exist_ok=False)

    v = np.array(vv)
    t = np.arange(len(v)) * REC_DT_MS
    i0 = crossings(t, v, SPIKE_MV)
    i20 = crossings(t, v, REPORT_MV)
    t0 = t[i0]
    n_peak = int(PEAK_MS / REC_DT_MS)
    peaks = np.array([v[i:i + n_peak].max() for i in i0])
    win_rates = [
        float(np.count_nonzero((t0 >= k * WINDOW_MS) & (t0 < (k + 1) * WINDOW_MS)) / (WINDOW_MS / 1000.0))
        for k in range(int(tstop_ms / WINDOW_MS))
    ]
    tail = t0[t0 >= tstop_ms - TAIL_MS]
    tail_rate = len(tail) / (TAIL_MS / 1000.0)
    tail_isi = np.diff(tail)
    tail_isi_cv = float(np.std(tail_isi) / np.mean(tail_isi)) if len(tail_isi) > 1 else None

    with open(out_dir / "SPIKES.csv", "x", encoding="utf-8", newline="") as f:
        f.write("time_ms,peak_mV\n")
        f.write("".join(f"{a:.4f},{b:.4f}\n" for a, b in zip(t0, peaks)))
    names = ["SPIKES.csv"]
    if save_v:
        with open(out_dir / "V_TRACE.csv", "x", encoding="utf-8", newline="") as f:
            f.write("time_ms,vsoma_mV\n")
            f.write("".join(f"{a:.4f},{b:.4f}\n" for a, b in zip(t, v)))
        names.append("V_TRACE.csv")
    record = {
        "stage": "Masoli 2015 population v07 (E7 exploratory)",
        "amendment": "S2,S3,S4",
        "member_id": member_id,
        "members_tsv_sha256": sha256(members_tsv),
        "member_scales": {k: member[k] for k in member},
        "factor_file": str(factor_json),
        "factor_file_sha256": sha256(factor_json),
        "condition": fdoc.get("condition"),
        "cdp5_scale": cdp5_scale,
        "cdp5_scales": cdp5_scales,
        "cdp5_sections": cdp5_sections,
        "run_s": run_s,
        "threads": THREADS,
        "expected_head": MODEL_HEAD,
        "py_sha256": {n: sha256(model_dir / n) for n in ("Purkinje_py3.py", "PC_param.py")},
        "nrnmech_lib_sha256": sha256(lib),
        "neuron_version": h.nrnversion(0),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "platform": platform.platform(),
        "dt_ms": DT_MS, "record_dt_ms": REC_DT_MS, "tstop_ms": tstop_ms,
        "time_axis": "index x record_dt_ms from 0 (t not recorded; multithread)",
        "n_samples": int(len(v)),
        "applied": applied,
        "noise": noise_info,
        "published": a.published,
        "tail_isi_cv": tail_isi_cv,
        "region_counts": counts,
        "n_spikes_0mV": int(len(i0)), "n_spikes_m20mV": int(len(i20)),
        "window_rates_hz_2s": win_rates,
        "tail_rate_hz_last10s": tail_rate,
        "saved": names,
    }
    with open(out_dir / "RECORD.json", "x", encoding="utf-8", newline="") as f:
        f.write(json.dumps(record, indent=1) + "\n")
    names.append("RECORD.json")
    sums = {n: sha256(out_dir / n) for n in names}
    with open(out_dir / "SHA256SUMS", "x", encoding="utf-8", newline="") as f:
        f.write("".join(f"{h_}  {n}\n" for n, h_ in sums.items()))
    print("done", f"spikes={len(i0)}", f"tail_rate_hz={tail_rate:.2f}")


if __name__ == "__main__":
    main()

#%% debug
