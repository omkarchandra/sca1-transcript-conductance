"""Phase 4 Track 2 Stage B driver  — v03 DRAFT.

v02 -> v03 (per review , drafting only): (a) NCX = exchanger only
(<AUTHOR> 00:12) — pump variant removed; ncx_scope field kept, refuses anything
but "exchanger_only"; "h" access kind dropped (nothing uses it); mypumpsoma
density stays at its published value (it is not section-addressable — v02
dry probe; effect of h-level sets UNVERIFIED, hence exchanger-only);
(b) EXPECT_RANGED reverted to the v04 cross-section tuples
(KNa soma/Couple, tau soma/Couple, KK Couple/soma); (c) dendritic count from
detect_spikes(tw, vdw) on the window slice; (d) LABEL_MAP asserted as a
subset of the sorter LABELS at startup — LABELS L21, Bounds L23-35,
classify(t,v,bounds) L91 of p4_t2_regime.py (fdab1cb9…); Bounds(**json)
is strict by dataclass.

One condition per invocation; scales ONLY from the factors JSON (sha
recorded); v_init required CLI (gate -80..-50 mV, reference -65).
Run with the pinned venv python ONLY (3.14.4):
  python -I -B 05_model_forrest/p4_t2_stageB.py <model_dir> <mech_lib> <factors_json> <bounds_json> <new_out_dir> <v_init_mv>
NOT RUN until review clears these exact bytes AND <AUTHOR> gives the Stage B
yes. Exclusive-create outputs; no deletion; no retry (STOP + verbatim).
This file never reads the frozen prediction or the B3 tables.
"""
from pathlib import Path
import hashlib
import json
import sys
import time

T_SETTLE_MS = 60000.0
T_WINDOW_MS = 120000.0
DT_MS = 0.025
REC_MS = 0.1
V_INIT_MIN_MV = -80.0
V_INIT_MAX_MV = -50.0
SPIKE_THRESHOLD_MV = 0.0
GAP_SILENT_MS = 1000.0
V_PLATEAU_MEASURE_MV = -40.0
SORTER_SHA256 = "b3714ca27132582b01fbfcf31ab9d4e25e95764af8292ec82027cab1cca71885"
SIX_IDS = {"g_BK", "g_Kv3", "g_KA", "g_CaP", "g_CaT", "p_NCX"}
EXPECT_GLOBALS = {"alcohol": 0.0, "celsius": 36.0}
EXPECT_RANGED = {("KNa", "soma", "Couple", "pump", "km", 40.0),
                 ("tau", "soma", "Couple", "nadifl", "tau", 5000.0),
                 ("KK", "Couple", "soma", "newnew", "km", 2.245)}
# (section, mechanism, variable) — every target is RANGE accessed via
# sec(0.5).var_mech with refuse-neither-way; no getattr on a bare Section.
TARGETS_BASE = {
    "g_BK": [("soma", "bkpkj", "gkbar"), ("Couple", "kc3", "gkbar")],
    "g_Kv3": [("soma", "kpkj", "gkbar")],
    "g_KA": [("Couple", "ka", "gkbar")],
    "g_CaP": [("soma", "captain", "pcabar"), ("Couple", "cap2", "gcabar")],
    "g_CaT": [("Couple", "cat", "gcabar")],
}
# <AUTHOR> 00:12 via review 00:25 — exchanger only, pumps stay published.
NCX_SCOPE_ONLY = "exchanger_only"
NCX_TARGETS = [("soma", "myexchangersoma", "ImaxNax"),
               ("Couple", "myexchanger", "ImaxNax")]
LABEL_MAP = {"TONIC": "pacemaking", "BURSTING": "bursting",
             "QUIESCENT": "silent-hyperpolarized",
             "DEPOLARIZATION_BLOCK": "depolarized plateau"}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json_exclusive(path, obj):
    with Path(path).open("x", encoding="utf-8", newline="") as stream:
        stream.write(json.dumps(obj, indent=1) + "\n")


def resolve(h, sec_name, mech, var):
    name = f"{var}_{mech}"
    sec = getattr(h, sec_name)
    try:
        return (f"{sec_name}(0.5).{name}", float(getattr(sec(0.5), name)))
    except AttributeError:
        raise ValueError(f"target resolves neither way: {sec_name}(0.5).{name}")


def read_ranged(h, primary, secondary, mech, var):
    out = []
    for sec_name in (primary, secondary):
        sec = getattr(h, sec_name)
        try:
            out.append((f"{sec_name}.{var}_{mech}", float(getattr(sec(0.5), f"{var}_{mech}"))))
        except AttributeError:
            continue
    return out


def detect_spikes(t, v, threshold=SPIKE_THRESHOLD_MV):
    t = [float(x) for x in t]
    v = [float(x) for x in v]
    spikes = []
    for i in range(len(v) - 1):
        if v[i] < threshold <= v[i + 1]:
            frac = (threshold - v[i]) / (v[i + 1] - v[i])
            spikes.append(t[i] + frac * (t[i + 1] - t[i]))
    return spikes


def main(argv=None):
    from neuron import h  # noqa: PLC0415 — venv NEURON only

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    sorter_path = Path(__file__).resolve().parent / "p4_t2_regime.py"
    sorter_sha = sha(sorter_path)
    if sorter_sha != SORTER_SHA256:
        raise ValueError("sorter differs from its pin fdab1cb9")
    import p4_t2_regime as regime  # noqa: PLC0415 — reviewed sorter
    if not set(LABEL_MAP) <= set(regime.LABELS):
        raise ValueError("LABEL_MAP disagrees with sorter LABELS")

    parts = (argv or sys.argv[1:])[:6]
    model_s, lib_s, fact_s, bounds_s, out_s, v_init_s = parts
    model_dir, mech_lib = Path(model_s), Path(lib_s)
    factors_path, bounds_path, out_dir = Path(fact_s), Path(bounds_s), Path(out_s)
    factors = json.loads(factors_path.read_text())
    bounds = regime.Bounds(**json.loads(bounds_path.read_text()))
    if set(factors["scales"]) != SIX_IDS:
        raise ValueError(f"factors JSON must name exactly {sorted(SIX_IDS)}")
    scope = factors.get("ncx_scope")
    if scope != NCX_SCOPE_ONLY:
        raise ValueError('factors JSON needs ncx_scope "exchanger_only"')
    v_init_mv = float(v_init_s)
    if not (V_INIT_MIN_MV <= v_init_mv <= V_INIT_MAX_MV):
        raise ValueError(f"v_init {v_init_mv} outside the allowed range [{V_INIT_MIN_MV},{V_INIT_MAX_MV}]")
    targets = dict(TARGETS_BASE)
    targets["p_NCX"] = NCX_TARGETS
    out_dir.mkdir(parents=False, exist_ok=False)
    started = time.time()
    h.nrn_load_dll(str(mech_lib))
    h.load_file("stdrun.hoc")
    h.load_file(str(model_dir / "2_compartment.hoc"))
    got = {"alcohol": float(h.alcohol), "celsius": float(h.celsius)}
    for key, want in EXPECT_GLOBALS.items():
        if got[key] != want:
            raise ValueError(f"published parameter differs: {key}={got[key]} != {want}")
    ranged = {}
    for label, primary, secondary, mech, var, want in sorted(EXPECT_RANGED):
        found = read_ranged(h, primary, secondary, mech, var)
        if not found:
            raise ValueError(f"parameter not found in any section: {var}_{mech}")
        for where, val in found:
            if val != want:
                raise ValueError(f"published parameter differs: {where}={val} != {want}")
        ranged[label] = {"expected": want, "read": found}
    applied = {}
    for cid, factor in factors["scales"].items():
        for sec_name, mech, var in targets[cid]:
            where, current = resolve(h, sec_name, mech, var)
            sec = getattr(h, sec_name)
            setattr(sec(0.5), f"{var}_{mech}", current * float(factor))
            applied[where] = {"access": "seg", "before": current,
                              "after": float(getattr(sec(0.5), f"{var}_{mech}"))}
    stims = []
    for ns in h.List("NetStim"):
        before = {"number": float(ns.number), "interval": float(ns.interval),
                  "start": float(ns.start), "noise": float(ns.noise)}
        ns.number = 0
        stims.append({"before": before, "after": {"number": 0.0}})
    h.dt = DT_MS
    h.celsius = 36.0
    soma = h.soma
    couple = h.Couple
    t_vec = h.Vector()
    v_vec = h.Vector()
    vd_vec = h.Vector()
    t_vec.record(h._ref_t, REC_MS)
    v_vec.record(soma(0.5)._ref_v, REC_MS)
    vd_vec.record(couple(0.5)._ref_v, REC_MS)
    h.finitialize(v_init_mv)
    h.continuerun(T_SETTLE_MS + T_WINDOW_MS)
    t = list(t_vec)
    v = list(v_vec)
    vd = list(vd_vec)
    if len(t) < 2:
        raise ValueError("empty recording")
    win = [i for i, x in enumerate(t) if x >= T_SETTLE_MS]
    tw = [t[i] - T_SETTLE_MS for i in win]
    vw = [v[i] for i in win]
    vdw = [vd[i] for i in win]
    spikes_all = detect_spikes(t, v)
    spikes = [s for s in spikes_all if s >= T_SETTLE_MS]
    dend_cross = detect_spikes(tw, vdw)
    gaps = [b - a for a, b in zip(spikes, spikes[1:]) if b - a >= GAP_SILENT_MS]
    gap_onsets = [b for a, b in zip(spikes, spikes[1:]) if b - a >= GAP_SILENT_MS]
    verdict = regime.classify(tw, vw, bounds)
    mapped = LABEL_MAP.get(verdict["label"], verdict["label"] + " (unmapped)")
    frac_above = sum(1 for x in vw if x > V_PLATEAU_MEASURE_MV) / len(vw)
    summary = {"window_ms": [T_SETTLE_MS, T_SETTLE_MS + T_WINDOW_MS],
               "n_samples_full": len(t), "n_samples_window": len(tw),
               "v_min_window": min(vw), "v_max_window": max(vw),
               "n_spikes_window": len(spikes),
               "n_quiescent_gaps": len(gaps), "gap_durations_ms": gaps,
               "cycle_periods_ms": [b - a for a, b in zip(gap_onsets, gap_onsets[1:])],
               "n_dend_crossings_0mV_supplementary": len(dend_cross),
               "frac_time_above_m40": frac_above,
               "sorter_label": verdict["label"], "mapped_class": mapped,
               "sorter_reason": verdict.get("reason")}
    with (out_dir / "V_TRACE.tsv").open("x", encoding="utf-8", newline="") as f:
        f.write("t_ms\tv_mV\n")
        for ti, vi in zip(t, v):
            f.write(f"{ti!r}\t{vi!r}\n")
    with (out_dir / "V_DEND_TRACE.tsv").open("x", encoding="utf-8", newline="") as f:
        f.write("t_ms\tv_mV\n")
        for ti, vi in zip(t, vd):
            f.write(f"{ti!r}\t{vi!r}\n")
    with (out_dir / "SPIKES.tsv").open("x", encoding="utf-8", newline="") as f:
        f.write("spike_ms\n")
        for s in spikes_all:
            f.write(f"{s!r}\n")
    record = {"condition": factors["condition"], "factors_sha256": sha(factors_path),
              "bounds_sha256": sha(bounds_path), "ncx_scope": scope,
              "sorter_sha256": sorter_sha, "python_version": sys.version.split()[0],
              "published_params": got, "ranged_params": ranged,
              "applied": applied, "netstims": stims,
              "v_init_mv": v_init_mv, "v_init_source": "cli-required (-80..-50 mV gate)",
              "t_settle_ms": T_SETTLE_MS, "t_window_ms": T_WINDOW_MS,
              "dt_ms": DT_MS, "rec_ms": REC_MS,
              "wall_seconds": time.time() - started,
              "neuron_version": h.nrnversion(0),
              "mech_lib_sha256": sha(mech_lib),
              "driver_sha256": sha(__file__),
              "model": "ModelDB 180789 2-compartment surrogate @1a63182",
              "summary": summary}
    write_json_exclusive(out_dir / "RECORD.json", record)
    sums = []
    for name in ("V_TRACE.tsv", "V_DEND_TRACE.tsv", "SPIKES.tsv", "RECORD.json"):
        sums.append(f"{sha(out_dir / name)}  {name}")
    with (out_dir / "SHA256SUMS").open("x", encoding="utf-8", newline="") as f:
        f.write("\n".join(sums) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
