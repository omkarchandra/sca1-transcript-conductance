"""Final Zang P4/P8/P9 tables from the per-run METRICS.tsv files at full precision.
 Reads every score folder written by zang_collect_20261027.py (scorer v02 ed3c7859, each run scored alone, so its
own d_* columns are nan). Writes ONE new folder with:
  RUNS.tsv        one row per run: task, condition, seed, all scorer metrics unrounded, min peak, wall, and the d_* paired
                  differences against the stage-2 WT run of the same seed, filled from the unrounded values (WT rows: 0);
  CONDITIONS.tsv  per condition: n pairs, median and IQR (numpy percentiles 25/50/75) of rate, rate ratio, dCV, dCV2,
                  rise counts, seeds with a pause, dspikes, lowest peak, and the three-part definition (median dCV2 > 0 AND
                  median rate ratio in [0.90, 1.10] AND a pause in >= 2 of 5 seeds), computed from unrounded values;
  SHA256SUMS.
No nan anywhere: a field that does not apply is written as "NA" with the reason in the header comment.
Usage: python -I -B zang_final_table_20261027.py <OUT_DIR>
"""
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

PUB = Path("${LOCAL_ROOT}/results/publication_20261026")
SCORES = PUB / "zang_stage1_20261026"
ROOTS = {"E": "zang_E_25s_sig0.0065_20261026", "P8": "zang_P8_atrophy_25s_20261026",
         "P9": "zang_P9_auxiliary_25s_20261026"}
KEYS = ("rate_hz", "cv", "cv2", "pause_frac", "log_isi_mode_ms", "dspikes", "main_crossings")
ORDER = ["WT", "KI_18w_mid", "KI_24w_mid", "KI_30w_mid", "P8_atrophy_A0.8", "P8_atrophy_A0.6", "ZP9_iv_aux_only",
         "ZP9_v_five_aux"]


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def runs():
    out = []
    for tag, root in ROOTS.items():
        for d in sorted(p for p in (PUB / root).iterdir() if p.is_dir()):
            sc = SCORES / f"score_{tag}_{d.name}_v02_20261027"
            m = next(csv.DictReader(open(sc / "METRICS.tsv"), delimiter="\t"))
            rec = json.load(open(d / "RECORD.json"))
            r = list(csv.DictReader(open(d / "SPIKES.csv")))
            t = np.array([float(x["time_ms"]) for x in r])
            pk = np.array([float(x["peak_mV"]) for x in r])
            w = (t >= 10000) & (t <= 25000)
            out.append(dict(task=(lambda t: f"{t.get('JOB_ID')}_{t.get('TASK_ID')}" if t else "")(rec.get("task", {})),
                            condition=rec["condition"], run=d.name, seed=int(m["seed"]), k=int(m["seed"]) - 20261013,
                            n_spikes=int(m["n_spikes"]), **{k: float(m[k]) for k in KEYS},
                            min_peak_mV=float(pk[w].min()), wall_h=rec["wall_s"] / 3600,
                            sha256sums=sha256(d / "SHA256SUMS"), metrics_sha256=sha256(sc / "METRICS.tsv")))
    return out


def q(x):
    p = np.percentile(np.asarray(x, dtype=float), [25, 50, 75])
    return float(p[1]), float(p[0]), float(p[2])


def main():
    out = Path(sys.argv[1])
    if out.exists():
        raise SystemExit(f"output exists: {out}")
    rows = runs()
    wt = {r["k"]: r for r in rows if r["condition"] == "WT"}
    for r in rows:
        w = wt.get(r["k"])
        for k in KEYS:
            r[f"d_{k}"] = r[k] - w[k] if w else "NA"
        r["rate_ratio"] = r["rate_hz"] / w["rate_hz"] if w else "NA"
    out.mkdir()
    cols = ["task", "condition", "run", "seed", "k", "n_spikes", *KEYS, "min_peak_mV", "wall_h", "rate_ratio",
            *[f"d_{k}" for k in KEYS], "sha256sums", "metrics_sha256"]
    with open(out / "RUNS.tsv", "x", newline="") as f:
        f.write("# per-run values at full precision; d_* and rate_ratio = run minus / over the stage-2 WT run of the same "
                "seed k (WT rows: 0 and 1 by construction); NA = no WT run of that seed\n")
        w_ = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w_.writeheader()
        for r in sorted(rows, key=lambda r: (ORDER.index(r["condition"]), r["k"])):
            w_.writerow({c: (repr(r[c]) if isinstance(r[c], float) else r[c]) for c in cols})
    ccols = ["condition", "n_runs", "n_pairs", "rate_med", "rate_q1", "rate_q3", "ratio_med", "ratio_q1", "ratio_q3",
             "dcv_med", "dcv_q1", "dcv_q3", "dcv2_med", "dcv2_q1", "dcv2_q3", "cv_rises", "cv2_rises",
             "seeds_with_pause", "dspikes_total", "lowest_peak_mV", "dcv2_gt_0", "ratio_in_band", "pause_ge_2_of_5",
             "moves_toward_phenotype"]
    with open(out / "CONDITIONS.tsv", "x", newline="") as f:
        f.write("# medians and quartiles: numpy percentile (linear) over the paired seeds, from unrounded values; "
                "three-part definition only for conditions with 5 pairs (else NA); WT row: pair fields NA\n")
        w_ = csv.DictWriter(f, fieldnames=ccols, delimiter="\t")
        w_.writeheader()
        for c in ORDER:
            rr = [r for r in rows if r["condition"] == c]
            if not rr:
                continue
            pr = [r for r in rr if r["k"] in wt] if c != "WT" else []
            row = dict(condition=c, n_runs=len(rr), n_pairs=len(pr), dspikes_total=int(sum(r["dspikes"] for r in rr)),
                       lowest_peak_mV=repr(min(r["min_peak_mV"] for r in rr)),
                       seeds_with_pause=sum(r["pause_frac"] > 0 for r in (pr or rr)))
            row.update(zip(("rate_med", "rate_q1", "rate_q3"), map(repr, q([r["rate_hz"] for r in rr]))))
            if pr:
                for name, key in (("ratio", "rate_ratio"), ("dcv", "d_cv"), ("dcv2", "d_cv2")):
                    row.update(zip((f"{name}_med", f"{name}_q1", f"{name}_q3"), map(repr, q([r[key] for r in pr]))))
                row["cv_rises"] = sum(r["d_cv"] > 0 for r in pr)
                row["cv2_rises"] = sum(r["d_cv2"] > 0 for r in pr)
                if len(pr) == 5:
                    a = float(np.median([r["d_cv2"] for r in pr])) > 0
                    b = 0.90 <= float(np.median([r["rate_ratio"] for r in pr])) <= 1.10
                    p = row["seeds_with_pause"] >= 2
                    row.update(dcv2_gt_0=a, ratio_in_band=b, pause_ge_2_of_5=p,
                               moves_toward_phenotype="MET" if (a and b and p) else "NOT MET")
            w_.writerow({k: row.get(k, "NA") for k in ccols})
    with open(out / "SHA256SUMS", "x") as f:
        f.write("".join(f"{sha256(out / n)}  {n}\n" for n in ("RUNS.tsv", "CONDITIONS.tsv")))
    print(open(out / "CONDITIONS.tsv").read())


if __name__ == "__main__":
    main()

#%% debug
