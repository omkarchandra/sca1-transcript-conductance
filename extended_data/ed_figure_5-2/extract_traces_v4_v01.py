"""Extract a 1-s window of somatic voltage from seven V_SOMA.npy files (published parameter set, with synaptic input)
into one TSV, plus a small JSON record. Read-only on inputs; outputs by exclusive create; deletes nothing.
Usage: python3 -I extract_traces_v4_v01.py IN_DIR OUT_DIR
IN_DIR holds seven subfolders named as below, each with V_SOMA.npy (copied from the wave-8 runs)."""
import sys, os, json, hashlib
import numpy as np
IN, OUT = sys.argv[1], sys.argv[2]
NAMES = [  # subfolder (copy of wave8_syn_<cond>/set000_<...>) -> column label
 ("wave8_syn_WT", "WT"),
 ("wave8_syn_KI_18w_mid", "all_five_transcript_ratios"),
 ("wave8_syn_KI_1c_restore_Cacna1a", "Cav2.1_unchanged"),
 ("wave8_syn_KI_4_protein_Kcnma1_Cacna1g", "staining_ratios"),
 ("wave8_syn_KI_5_protein_Kcnma1_Cacna1g_Cacna1a_0p85", "staining_Cav2.1_0.85"),
 ("wave8_syn_KI_9_area_A0p8", "area_A0.8"),
 ("wave8_syn_KI_9_area_A0p75", "area_A0.75"),
]
N_EXPECT, DT = 100001, 0.1          # 10 s at 0.1 ms, as in the 15 Nov extraction
I0, I1 = 50000, 60000               # 5000 to 6000 ms inclusive (10,001 rows)
os.makedirs(OUT, exist_ok=False)
cols, rec = [], {"window_ms": [I0*DT, I1*DT], "dt_ms": DT, "inputs": {}}
for sub, lab in NAMES:
    p = os.path.join(IN, sub, "V_SOMA.npy")
    b = open(p, "rb").read()
    v = np.load(p, allow_pickle=False)
    if v.ndim != 1 or v.shape[0] != N_EXPECT:
        raise SystemExit(f"STOP: {p} has shape {v.shape}, expected ({N_EXPECT},)")
    up = lambda x: int(np.sum((x[1:] >= -20.0) & (x[:-1] < -20.0)))
    rec["inputs"][lab] = {"path": p, "sha256": hashlib.sha256(b).hexdigest(), "dtype": str(v.dtype),
                          "spikes_10s": up(v), "spikes_window": up(v[I0:I1+1]),
                          "min_mV": float(v.min()), "max_mV": float(v.max())}
    cols.append((lab, v[I0:I1+1].astype(float)))
fn = os.path.join(OUT, "TRACES_v4_with_synaptic_input.tsv")
with open(fn, "x") as f:
    f.write("time_ms\t" + "\t".join(l for l, _ in cols) + "\n")
    for i in range(I1 - I0 + 1):
        f.write(f"{(I0+i)*DT:.1f}\t" + "\t".join(f"{c[i]:.3f}" for _, c in cols) + "\n")
rec["tsv_sha256"] = hashlib.sha256(open(fn, "rb").read()).hexdigest()
with open(os.path.join(OUT, "RECORD.json"), "x") as f: json.dump(rec, f, indent=1)
print(json.dumps({k: (v["spikes_window"], v["spikes_10s"]) for k, v in rec["inputs"].items()}))
print("tsv", rec["tsv_sha256"])
