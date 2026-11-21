"""Phase 4 Track 2 Stage A run driver  — v03 DRAFT.

v02 -> v03 (per review , nothing else changes): parameter refusal
check reads RANGE variables via explicit sections — KNa via
h.soma(0.5).km_pump (hoc L68, inside `soma {` L38), tau via
h.soma(0.5).tau_nadifl (hoc L77), KK via h.Couple(0.5).km_newnew (hoc L151,
inside `Couple {` L99). If KNa/tau are also set in Couple, or KK in soma,
both values are read and recorded; refusal unless every value read equals
40 / 5000 / 2.245. (v02 failed because bare h.km_pump is illegal for RANGE.)

Reproduces Forrest 2012 Fig 1C: the 2-compartment surrogate at its published
WT parameters fires a spontaneous repeating trimodal pattern (tonic, burst,
silence; ~20 s repeat). No B3 factors, no KI parameters, no project data, no
frozen-prediction reference (not read, not named by path here).

Run with the pinned venv python ONLY:
  python -I -B 05_model_forrest/p4_t2_stageA_run.py <model_dir> <mech_lib> <new_out_dir>
NOT RUN until review clears these exact bytes AND <AUTHOR> gives the re-run yes.
One run, exclusive-create outputs, no deletion, no retry on failure (STOP +
verbatim report).

Published-parameter handling: NOTHING in the model is changed except the
stimulus definition. After loading 2_compartment.hoc the driver READS
alcohol, celsius (globals) plus the section-based km/tau values and REFUSES
unless they equal 0, 36, 40, 2.245, 5000. Rationale recorded here: the hoc
wires NetStim->synapse connections with nonzero weights, but Fig 1C is "in
the absence of synaptic input", so every NetStim is captured
(number/interval/start/noise) into RECORD.json and then silenced via
number=0 — the single declared condition change. GUI-only files
(lesbos.ses, simple_data_writer.hoc, mex.hoc which is commented out upstream)
are NOT loaded; headless run uses stdrun.hoc + 2_compartment.hoc only.
Sections `soma` and `Couple` confirmed from 2_compartment.hoc L8-9.
"""
from pathlib import Path
import hashlib
import json
import re
import sys
import time

T_TOTAL_MS = 90000.0
DT_MS = 0.025
REC_MS = 0.1
V_INIT_DEFAULT_MV = -70.0
SPIKE_THRESHOLD_MV = 0.0
GAP_SILENT_MS = 1000.0
EXPECT_GLOBALS = {"alcohol": 0.0, "celsius": 36.0}
EXPECT_RANGED = {("KNa", "soma", "Couple", "pump", "km", 40.0),
                 ("tau", "soma", "Couple", "nadifl", "tau", 5000.0),
                 ("KK", "Couple", "soma", "newnew", "km", 2.245)}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json_exclusive(path, obj):
    with Path(path).open("x", encoding="utf-8", newline="") as stream:
        stream.write(json.dumps(obj, indent=1) + "\n")


def read_ranged(h, primary, secondary, mech, var):
    """Read var_mech from primary section and, if also present, secondary;
    returns a list of (where, value). Absence in a section is not an error."""
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

    model_dir, mech_lib, out_dir = (Path(a) for a in (argv or sys.argv[1:])[:3])
    out_dir.mkdir(parents=False, exist_ok=False)
    started = time.time()
    h.nrn_load_dll(str(mech_lib))
    h.load_file("stdrun.hoc")
    hoc_text = (model_dir / "2_compartment.hoc").read_text()
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
    if re.search(r"(?m)^\s*v_init\s*=", hoc_text):
        v_init_mv, v_init_source = float(h.v_init), "hoc v_init"
    else:
        v_init_mv, v_init_source = V_INIT_DEFAULT_MV, "default -70 (hoc sets none)"
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
    h.continuerun(T_TOTAL_MS)
    t = list(t_vec)
    v = list(v_vec)
    vd = list(vd_vec)
    if len(t) < 2:
        raise ValueError("empty recording")
    spikes = detect_spikes(t, v)
    dend_cross = detect_spikes(t, vd)
    gaps = [b - a for a, b in zip(spikes, spikes[1:]) if b - a >= GAP_SILENT_MS]
    gap_onsets = [b for a, b in zip(spikes, spikes[1:]) if b - a >= GAP_SILENT_MS]
    summary = {"t_ms": [t[0], t[-1]], "n_samples": len(t),
               "v_min": min(v), "v_max": max(v),
               "n_spikes": len(spikes),
               "n_quiescent_gaps": len(gaps), "gap_durations_ms": gaps,
               "cycle_periods_ms": [b - a for a, b in zip(gap_onsets, gap_onsets[1:])],
               "n_dend_crossings_0mV_supplementary": len(dend_cross)}
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
        for s in spikes:
            f.write(f"{s!r}\n")
    record = {"published_params": got, "ranged_params": ranged, "netstims": stims,
              "condition": "no synaptic input (all NetStim number=0)",
              "v_init_mv": v_init_mv, "v_init_source": v_init_source,
              "t_total_ms": T_TOTAL_MS, "dt_ms": DT_MS, "rec_ms": REC_MS,
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
