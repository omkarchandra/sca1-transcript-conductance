# Zang & De Schutter 2021 Purkinje cell (ModelDB 266864) driver with background synapses, P16 v01.
# P16 v01 = P14 v01 (883e987e...) plus the parameter-set path of the Zang 2018 P13 driver p13 v01 (36c91de7...) and
#   nothing else: an optional --popset PATH, one parameter set
#   {"set_id": int, "seed": int, "scales": {Kv3, Kv4, BK, newCaP, CaT3_1, naRsg, hpkj, SK2}}, each factor finite and in
#   [0.5, 1.5], validated before anything is built. Each factor multiplies the published density/permeability of its
#   group (Kv4 = Kv4 + Kv4s; BK = mslo + abBK; see POP) in every segment of every section carrying the mechanism. Order
#   per segment: (a) population scale (this set's WT), (b) five-gene factors, (c) aux block, all unchanged code.
#   The set factors are applied after the background synapses are built and before --inspect (INSPECT shows the scaled
#   values); with a set, "published" in each RECORD.applied row is the value before the gene factor.
#   RECORD/SMOKE/INSPECT carry popset_file, popset_sha256, popset_set_id (the set number is the file's set_id).
#   Not taken from P13: --prespec, --prespec-commit and the fields popset_seed, popset_scales, pop_applied,
#   pop_section_counts, prespec_*. Without --popset, behaviour is that of P14 v01 (the three new fields null); an
#   all-1.0 set gives identical values.
# Lineage (P14 driver):
# Zang & De Schutter 2021 Purkinje cell (ModelDB 266864) driver with background synapses P14 v01
# (plan P14_zang2021_synaptic_plan_20261110.md, sha256 036b0a78...).
# P14 v01 = P4 driver v05 (5847a7cb...) moved to the 2021 model:
#   Set-up, all from the published fig8.hoc (hash-pinned), executed byte-for-byte after stdrun.hoc:
#     (1) biophysics: from after load_file("nrngui.hoc") to just before "celsius = 34";
#     (2) parameters: "celsius = 34" through "Freq_bs = 14.4";
#     (3) xopen("background_syn_distrib.hoc"), as published;
#     (4) background synapses: from "// following are background synapses" to just before "forsec br1 {",
#         with randomseed0 = --syn-seed (the published per-trial number is replaced by ours).
#   Omitted: stim1 (-0.3 nA for a random 0-25 ms), the CCn Normal(0, 0.6) nA current (100 ms), the two clustered
#   syn2 parallel-fibre test inputs (onset 406 ms), electrode.hoc, distri.hoc, distri_synapse.hoc, the MPI farm.
#   v_init -70 as published in 2021. No climbing fibre exists in the 2021 scripts.
#   Each driven NetStim gets one extra NetCon without a target that only records its event times:
#   SYN_EVENTS.npy = float64 (n, 2): time_ms, id (PF i; 10000 + stellate i; 20000 + basket i), in delivery order.
#   --seedtest: no model is loaded; tests whether NetStim.seed() sets one global generator (SEEDTEST.json).
#   RECORD adds per-region section counts for the five-gene factors and for each aux change (abBK is also in
#   the spiny dendrites in 2021; the beta4 shift applies in every section that carries abBK, as before).
# Lineage (2018 driver):
# Zang, Dieudonne & De Schutter 2018 Purkinje cell (ModelDB 243446) driver P4 v05.
# v05 = v04 (a363c9d3...) with the beta4 shift inverted to the plan's definition:
#   r = Kcnmb4 ratio = fraction of abBK that remains; in each abBK section gabkbar x r, gbar_mslo + (1 - r) x gabkbar.
# v04 = v03 (91407de6...) plus an optional "aux" object in the factor file, applied
#   after the five-gene scaling to every segment of every section carrying the mechanism; RECORD.aux_applied lists
#   each value. Allowed keys only (AUX): newCaP.pcabar x, Kv4.vshift +, Kv4s.vshift + (vshift enters hinf only:
#   V1/2 = vh - vshift), and beta4_shift_r (definition as in the v05 line above).
#   Without aux, v04 does exactly what v03 does.
# v03 = v02 (961524e5...) with the soma spike level moved from 0 mV to -30 mV: in this model the somatic spike
#   peaks lie between about -12 and -1 mV (AIS +4 to +12 mV), so 0 mV missed all but the start-up spike.
# v02 = v01 (a81bcb95...) plus --hold-nA: a DC current added to the OU noise in the same somaA(0.5) IClamp
# written to RECORD/SMOKE as stimulus.hold_nA. The noise series and its sha are unchanged.
# Set-up: the published purkinje_demo.hoc block from after load_file("nrngui.hoc") up to xopen("dendv_arnd21.ses")
# is executed byte-for-byte (stdrun.hoc replaces nrngui.hoc). The 45-s protocol replaces the published stimulation:
# dendv_arnd21.ses (graphs), electrode.hoc (stim1, amp 0 -> cfdistrib.hoc: 500 syn2 climbing-fibre synapses at
# 400 ms) and distri.hoc (shape plot) are not loaded. celsius 34 and v_init -65-3 as published; dt 0.02 ms.
# Five-gene factors scale every segment of every section carrying the mapped mechanism(s); RECORD lists each value.
# Usage: python -I -B <this file> <MODEL_DIR> <OUT_DIR> <NRNMECH_LIB> <FACTOR_JSON> <RUN_S>
#        [--noise-sigma-nA S --noise-tau-ms T --noise-seed N] [--hold-nA I] [--inspect | --smoke-ms MS]
#        P14: --syn-seed N is required (except with --seedtest); --seedtest needs only OUT_DIR to be new.
#        P16: [--popset POPSET_JSON] (optional; the set number is the file's set_id)
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
V_INIT = -70.0
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
    "fig8.hoc": "ab94d32ed928bdf3c859e9087be5518783e4a06b15bdc2f387c5153a257726b2",
    "Purkinje19b972-1.nrn": "927f2740df707c1433928ba63559abf155e29fc562cf8e4d2b395fc7846657b3",
    "background_syn_distrib.hoc": "2c06c0b30245ac29d8d98f95071d4f533a42cc7b5045cbe40cee1e68ce2bc963",
    "dendv_arnd21.ses": "aad7517d7492de21f4b27f3ce5cc7a01d2c011e03ed2237351dd6224c611ef5c",
    "electrode.hoc": "1ada2187417612c1de70c0dcc727e26e917ce88dd6418569da97eacfeeb19184",
}
HEAD_MARK = 'load_file("nrngui.hoc")\n'
BIO_END_MARK = "\ncelsius = 34\n"  # the biophysics block ends at the newline before it
PARAM_END_MARK = "\nFreq_bs = 14.4\n"
BG_MARK = "// following are background synapses\n"
BG_END_MARK = "forsec br1 {\n"
PUBLISHED_LINES = ("\ncelsius = 34\n", "\ndt = 0.02\n", "\nv_init = -70\n",
                   '\nxopen("background_syn_distrib.hoc")\n')
SYN_POPS = (("pf", 0, "Ngrc", "ampsynprelist", "ampsynlist", "ampcon", "Freq_pf"),
            ("st", 10000, "NST", "gabastsynprelist", "gabastsynlist", "gabastcon", "Freq_st"),
            ("bs", 20000, "NBS", "gababssynprelist", "gababssynlist", "gababscon", "Freq_bs"))
MECHS = {  # gene -> (mechanism, density parameter); mapping
    "Kcnc3": (("Kv3", "gbar_Kv3"),),
    "Kcnd3": (("Kv4", "gbar_Kv4"), ("Kv4s", "gbar_Kv4s")),
    "Kcnma1": (("mslo", "gbar_mslo"), ("abBK", "gabkbar_abBK")),
    "Cacna1a": (("newCaP", "pcabar_newCaP"),),
    "Cacna1g": (("CaT3_1", "pcabar_CaT3_1"),),
}
DENSITY = {"Kv3": "gbar_Kv3", "Kv4": "gbar_Kv4", "Kv4s": "gbar_Kv4s", "Kv1": "gbar_Kv1", "mslo": "gbar_mslo",
           "abBK": "gabkbar_abBK", "SK2": "gkbar_SK2", "newCaP": "pcabar_newCaP", "CaT3_1": "pcabar_CaT3_1",
           "naRsg": "gbar_naRsg", "nap": "gbar_nap", "hpkj": "ghbar_hpkj", "pas": "g_pas"}
POP = {  # population group -> ((mechanism, density parameter), ...); same attributes as MECHS/DENSITY (P13 pre-spec 3)
    "Kv3": (("Kv3", "gbar_Kv3"),),
    "Kv4": (("Kv4", "gbar_Kv4"), ("Kv4s", "gbar_Kv4s")),
    "BK": (("mslo", "gbar_mslo"), ("abBK", "gabkbar_abBK")),
    "newCaP": (("newCaP", "pcabar_newCaP"),),
    "CaT3_1": (("CaT3_1", "pcabar_CaT3_1"),),
    "naRsg": (("naRsg", "gbar_naRsg"),),
    "hpkj": (("hpkj", "ghbar_hpkj"),),
    "SK2": (("SK2", "gkbar_SK2"),),
}
POP_MIN, POP_MAX = 0.5, 1.5
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


def validate_popset(doc):
    """Return (set_id, seed, {group: float}) or exit non-zero. Exactly the eight POP keys, each finite in [0.5, 1.5]."""
    def num(x):
        return isinstance(x, (int, float)) and not isinstance(x, bool)
    if not isinstance(doc, dict):
        raise SystemExit("popset must be a JSON object")
    for k in ("set_id", "seed"):
        if not isinstance(doc.get(k), int) or isinstance(doc.get(k), bool):
            raise SystemExit(f"popset {k} must be an integer")
    sc = doc.get("scales")
    if not isinstance(sc, dict):
        raise SystemExit("popset must contain a scales object")
    if set(sc) != set(POP):
        raise SystemExit(f"popset scales keys must be exactly {sorted(POP)}: missing {sorted(set(POP) - set(sc))}, "
                         f"extra {sorted(set(sc) - set(POP))}")
    out = {}
    for g in POP:
        x = sc[g]
        if not num(x) or not math.isfinite(float(x)) or not POP_MIN <= float(x) <= POP_MAX:
            raise SystemExit(f"popset factor {g}={x!r} must be finite and within [{POP_MIN}, {POP_MAX}]")
        out[g] = float(x)
    return doc["set_id"], doc["seed"], out


def apply_popset(h, scales, lists):
    """Population scaling of the published density of every segment of every section carrying the mechanism."""
    applied, counts = [], {}
    for group, pairs in POP.items():
        gf = scales[group]
        for mech, attr in pairs:
            n = 0
            for sec in h.allsec():
                if not sec.has_membrane(mech):
                    continue
                published = seg_values(sec, attr)
                for seg in sec:
                    setattr(seg, attr, getattr(seg, attr) * gf)
                applied.append({"group": group, "mechanism": mech, "attr": attr, "section": sec.name(),
                                "region": region_of(h, sec, lists), "nseg": sec.nseg, "published": published,
                                "factor": gf, "applied": seg_values(sec, attr)})
                n += 1
            counts[f"{group}|{mech}"] = n
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


def region_counts(rows, key_of):
    out = {}
    for r in rows:
        d = out.setdefault(key_of(r), {})
        d[r["region"]] = d.get(r["region"], 0) + 1
    return out


def build_background(h, background, syn_seed, lists):
    """Published background synapses (fig8.hoc block) with randomseed0 = syn_seed, plus event recorders."""
    if not h('xopen("background_syn_distrib.hoc")'):
        raise SystemExit("background_syn_distrib.hoc failed in hoc")
    if not h(f"randomseed0 = {int(syn_seed)}"):
        raise SystemExit("could not set randomseed0")
    if not h(background):
        raise SystemExit("published background-synapse block failed in hoc")
    ev_t, ev_id, recorders, pops = h.Vector(), h.Vector(), [], {}
    for key, offset, n_name, pre_name, syn_name, con_name, freq_name in SYN_POPS:
        n_driven, pre, syn = int(getattr(h, n_name)), getattr(h, pre_name), getattr(h, syn_name)
        for i in range(n_driven):
            nc = h.NetCon(pre.object(i), None)
            nc.record(ev_t, ev_id, offset + i)
            recorders.append(nc)
        sites = {}
        for j in range(int(syn.count())):
            reg = region_of(h, syn.object(j).get_segment().sec, lists)
            sites[reg] = sites.get(reg, 0) + 1
        s0, p0, c0 = syn.object(0), pre.object(0), getattr(h, con_name)[0]
        pops[key] = {"netstim_driven": n_driven, "netstim_created": int(pre.count()),
                     "synapses": int(syn.count()), "synapse_regions": sites, "event_id_offset": offset,
                     "rate_hz": float(getattr(h, freq_name)), "interval_ms": float(p0.interval),
                     "number": float(p0.number), "start_ms": float(p0.start), "noise": float(p0.noise),
                     "tau1_ms": float(s0.tau1), "tau2_ms": float(s0.tau2), "e_mV": float(s0.e),
                     "weight_uS": float(c0.weight[0])}
    info = {"source": "fig8.hoc background block, byte-for-byte; background_syn_distrib.hoc by xopen",
            "block_sha256": hashlib.sha256(background.encode("utf-8")).hexdigest(),
            "randomseed0": int(h.randomseed0), "SEED_0": int(h.SEED_0), "SEED_1": int(h.SEED_1),
            "SEED_2": int(h.SEED_2),
            "seed_calls": "NetStim i of each population: seed(i + SEED_k + randomseed0), k = 0 PF, 1 stellate, "
                          "2 basket; use_mcell_ran4(1) as published",
            "populations": pops, "n_recorders": len(recorders),
            "events_file": None, "events_sha256": None, "n_events": None}
    return info, ev_t, ev_id, recorders


def events_array(ev_t, ev_id, info):
    ev = np.column_stack([np.array(ev_t, dtype=np.float64), np.array(ev_id, dtype=np.float64)])
    ids = ev[:, 1].astype(int) if len(ev) else np.zeros(0, dtype=int)
    info["n_events"] = {"pf": int((ids < 10000).sum()), "st": int(((ids >= 10000) & (ids < 20000)).sum()),
                        "bs": int((ids >= 20000).sum()), "total": int(len(ids))}
    info["events_sha256"] = hashlib.sha256(ev.tobytes()).hexdigest()
    return ev


def seedtest(out_dir):
    """Does NetStim.seed() set one global generator? Three noisy NetStims, 2 s, no model."""
    from neuron import h
    h.load_file("stdrun.hoc")
    h("use_mcell_ran4(1)")
    sec = h.Section(name="seedtest")

    def events(seeds):
        stims, cons, tv, iv = [], [], h.Vector(), h.Vector()
        for i, s in enumerate(seeds):
            ns = h.NetStim()
            ns.interval, ns.number, ns.start, ns.noise = 10, 1e9, 0, 1
            ns.seed(s)
            nc = h.NetCon(ns, None)
            nc.record(tv, iv, i)
            stims.append(ns)
            cons.append(nc)
        h.finitialize(-65)
        h.continuerun(2000)
        out = [[float(x) for x in tv], [int(x) for x in iv]]
        del cons, stims
        return out

    a, a2, b, c = events([1, 2, 3]), events([1, 2, 3]), events([7, 8, 3]), events([1, 2, 4])
    doc = {"test": "three NetStim (interval 10 ms, noise 1), seeds given in creation order, 2000 ms, "
                   "use_mcell_ran4(1)",
           "n_events": [len(x[0]) for x in (a, a2, b, c)],
           "repeat_same_seeds_identical": a == a2,
           "only_earlier_seeds_changed_identical": a == b,
           "only_last_seed_changed_identical": a == c,
           "one_global_generator_last_seed_call_fixes_stream": a == a2 and a == b and a != c,
           "neuron_version": h.nrnversion(0), "hostname": socket.gethostname()}
    del sec
    out_dir.mkdir()
    with open(out_dir / "SEEDTEST.json", "x", encoding="utf-8", newline="") as f:
        f.write(json.dumps(doc, indent=1) + "\n")
    print("seedtest", doc["one_global_generator_last_seed_call_fixes_stream"], doc["n_events"])


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
                 if re.search(r"stim1\.|CCn|cfstart =|cfdistrib|tstop =|Freq_|SEED_|\.seed\(|distSyns\(", ln)]
             for f in ("fig8.hoc", "electrode.hoc", "background_syn_distrib.hoc")}
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
    ap.add_argument("--syn-seed", type=int, default=None)
    ap.add_argument("--seedtest", action="store_true")
    ap.add_argument("--inspect", action="store_true")
    ap.add_argument("--smoke-ms", type=float, default=None)
    ap.add_argument("--popset", default=None)
    a = ap.parse_args()
    model_dir, out_dir, lib, factor_json = Path(a.model_dir), Path(a.out_dir), Path(a.lib), Path(a.factor_json)
    run_s = float(a.run_s)
    if out_dir.exists():
        raise SystemExit(f"output exists; refusing to overwrite: {out_dir}")
    if a.seedtest:
        seedtest(out_dir)
        return
    if a.syn_seed is None or a.syn_seed < 0:
        raise SystemExit("--syn-seed (a non-negative integer) is required")
    if a.noise_sigma_nA > 0 and a.noise_seed is None:
        raise SystemExit("--noise-seed is required when sigma > 0")
    if (model_dir / "x86_64").exists():
        raise SystemExit("model dir contains x86_64; NEURON would auto-load it")
    for f, pin in PINS.items():
        if pin is not None and sha256(model_dir / f) != pin:
            raise SystemExit(f"published file hash mismatch: {f}")
    text_by_file = {f: (model_dir / f).read_text() for f in PINS}
    fig8 = text_by_file["fig8.hoc"]
    if any(fig8.count(m) != 1 for m in (HEAD_MARK, BIO_END_MARK, PARAM_END_MARK, BG_MARK, BG_END_MARK)):
        raise SystemExit("set-up markers not found exactly once in fig8.hoc")
    i0, i1 = fig8.index(HEAD_MARK) + len(HEAD_MARK), fig8.index(BIO_END_MARK) + 1
    i2 = fig8.index(PARAM_END_MARK) + len(PARAM_END_MARK)
    i3, i4 = fig8.index(BG_MARK), fig8.index(BG_END_MARK)
    if not i0 < i1 < i2 < i3 < i4:
        raise SystemExit("set-up markers out of order in fig8.hoc")
    setup, params, background = fig8[i0:i1], fig8[i1:i2], fig8[i3:i4]
    missing = [s.strip() for s in PUBLISHED_LINES if s not in fig8]
    if missing:
        raise SystemExit(f"published settings not found: {missing}")
    popset_path = Path(a.popset) if a.popset is not None else None
    pop_set_id = pop_seed = pop_scales = None
    if popset_path is not None:
        pop_set_id, pop_seed, pop_scales = validate_popset(json.loads(popset_path.read_text()))
    pop_info = {"popset_file": str(popset_path) if popset_path else None,
                "popset_sha256": sha256(popset_path) if popset_path else None,
                "popset_set_id": pop_set_id}
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
    if not h(params):
        raise SystemExit("published parameter block failed in hoc")
    h.celsius = CELSIUS
    h.dt = DT_MS
    h.steps_per_ms = 1.0 / DT_MS
    h.v_init = V_INIT
    h.CVode().active(0)
    syn_info, ev_t, ev_id, recorders = build_background(h, background, a.syn_seed, lists)

    pop_applied, pop_counts = apply_popset(h, pop_scales, lists) if pop_scales is not None else (None, None)
    if a.inspect:
        doc = inspect(h, lists, text_by_file)
        doc["background"] = syn_info
        doc.update(pop_info)
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
    ev = events_array(ev_t, ev_id, syn_info)
    n_pk, n_dpk = int(PEAK_MS / REC_DT_MS), int(DEND_PEAK_MS / REC_DT_MS)
    out_dir.mkdir()
    if a.smoke_ms is not None:
        smoke = {"stage": "Zang 2021 P16 driver v01 SMOKE", "smoke_ms": tstop_ms, "threads": threads,
                 "background": syn_info, "syn_event_times_ms_first20": [float(x) for x in ev[:20, 0]],
                 **pop_info,
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
    with open(out_dir / "SYN_EVENTS.npy", "xb") as f:
        np.save(f, ev)
    names.append("SYN_EVENTS.npy")
    syn_info["events_file"] = "SYN_EVENTS.npy (float64, columns time_ms, id; sha256 below is of the array bytes)"
    record = {
        "stage": "Zang & De Schutter 2021 (ModelDB 266864) P16 driver v01 (P14 v01 6511f2ec + popset)",
        "condition": fdoc.get("condition"),
        "factor_file": str(factor_json),
        "factor_file_sha256": sha256(factor_json),
        "scales": {g: float(factors[g]) for g in MECHS},
        "mapping": {g: [list(p) for p in pairs] for g, pairs in MECHS.items()},
        "mechanism_section_counts": counts,
        "mechanism_region_counts": region_counts(applied, lambda r: f"{r['gene']}|{r['mechanism']}"),
        "run_s": run_s, "tstop_ms": tstop_ms, "dt_ms": DT_MS, "published_dt_ms": 0.02,
        "record_dt_ms": REC_DT_MS, "celsius": h.celsius, "v_init": V_INIT,
        "time_axis": "index x record_dt_ms from 0",
        "n_samples": int(len(t)),
        "setup": "fig8.hoc: biophysics (after load_file(\"nrngui.hoc\") to before \"celsius = 34\"), parameters "
                 "(\"celsius = 34\" .. \"Freq_bs = 14.4\"), xopen background_syn_distrib.hoc, background block "
                 "(\"// following are background synapses\" to before \"forsec br1 {\"); byte-for-byte after "
                 "stdrun.hoc, with randomseed0 = syn_seed",
        "setup_sha256": {"biophysics": hashlib.sha256(setup.encode("utf-8")).hexdigest(),
                         "parameters": hashlib.sha256(params.encode("utf-8")).hexdigest(),
                         "background": hashlib.sha256(background.encode("utf-8")).hexdigest()},
        "omitted_published": ["stim1 IClamp (-0.3 nA at somaA for a random 0-25 ms)",
                              "CCn IClamp (amp played from Normal(0, 0.6) nA, 100 ms)",
                              "two clustered syn2 parallel-fibre test inputs (distSyns, onset 406 ms)",
                              "electrode.hoc", "distri.hoc (shape plot)", "distri_synapse.hoc",
                              "MPI trial farm (ParallelContext), tstop 600 ms",
                              "climbing fibre: none in the 2021 scripts (cfdistrib.hoc commented out, not shipped)"],
        "syn_seed": a.syn_seed,
        "background": syn_info,
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
        **pop_info,
        "applied": applied,
        "aux": aux, "aux_section_counts": aux_counts, "aux_applied": aux_applied,
        "aux_region_counts": region_counts(aux_applied, lambda r: r["key"]),
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
