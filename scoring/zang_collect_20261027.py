"""Collect finished Zang runs (stage 2, P8, P9) from local copies of the run outputs, verify, score each with scorer v02, and print one table.
27 Oct 2026. Reads each run folder under the local results root, checks its own
SHA256SUMS, scores it alone with zang_score_20261026_v02.py (ed3c7859) into a new score folder, and adds the minimum
soma spike peak and steadiness (1-s counts over 10-25 s). Pairs each non-WT run with the stage-2 WT run of the same seed.
Usage: python -I -B zang_collect_20261027.py
"""
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

PROJ = Path("${LOCAL_ROOT}")
PUB = PROJ / "results/publication_20261026"
SCORER = PROJ / "docs/20261026_cc/zang_score_20261026_v02.py"
SCORER_SHA = "ed3c78598deb3e8363d1192ee0fff704c09ea93ec4964d9e63cd30ecb03b6fbb"
SCORES = PUB / "zang_stage1_20261026"
SCRATCH = Path("/tmp/scratch")
RUNS = PUB
ROOTS = {"zang_E_25s_sig0.0065_20261026": 5493773, "zang_P8_atrophy_25s_20261026": 5501712,
         "zang_P9_auxiliary_25s_20261026": 5502003}


def sh(cmd):
    return subprocess.run(cmd, capture_output=True, text=True).stdout


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def local_done():
    out = {}
    for root in ROOTS:
        base = RUNS / root
        out[root] = sorted(d.name for d in base.iterdir() if (d / "SHA256SUMS").is_file()) if base.is_dir() else []
    return out


def extra(d):
    r = list(csv.DictReader(open(d / "SPIKES.csv")))
    t = np.array([float(x["time_ms"]) for x in r])
    p = np.array([float(x["peak_mV"]) for x in r])
    w = (t >= 10000) & (t <= 25000)
    c = np.histogram(t[t >= 10000], bins=np.arange(10000, 25001, 1000))[0]
    rec = json.load(open(d / "RECORD.json"))
    return dict(min_peak=float(p[w].min()) if w.any() else float("nan"), c_min=int(c.min()), c_max=int(c.max()),
                wall_h=rec["wall_s"] / 3600, task=(lambda t: f"{t.get('JOB_ID')}_{t.get('TASK_ID')}" if t else "")(rec.get("task", {})),
                condition=rec["condition"])


def score(root, name):
    out = SCORES / f"score_{root.split('_')[1]}_{name}_v02_20261027"
    if not out.exists():
        tmp = SCRATCH / f"root_{root}_{name}"
        tmp.mkdir(parents=True, exist_ok=True)
        link = tmp / name
        if not link.exists():
            link.symlink_to((RUNS / root / name).resolve())
        subprocess.run([sys.executable, "-I", "-B", str(SCORER), str(tmp), str(out)], capture_output=True, check=True)
    row = next(csv.DictReader(open(out / "METRICS.tsv"), delimiter="\t"))
    return out, row


def main():
    if sha256(SCORER) != SCORER_SHA:
        raise SystemExit("scorer hash differs")
    rows = []
    for root, names in local_done().items():
        for n in names:
            loc = RUNS / root / n
            for line in (loc / "SHA256SUMS").read_text().splitlines():
                h, f = line.split("  ", 1)
                if sha256(loc / f) != h:
                    raise SystemExit(f"checksum mismatch: {loc / f}")
            out, m = score(root, n)
            rows.append(dict(root=root, name=n, seed=int(m["seed"]), sums=sha256(loc / "SHA256SUMS")[:8],
                             metrics=sha256(out / "METRICS.tsv")[:8], **{k: float(m[k]) for k in
                             ("rate_hz", "cv", "cv2", "pause_frac", "log_isi_mode_ms", "dspikes")}, **extra(loc)))
    wt = {r["seed"]: r for r in rows if r["name"].startswith("WT_")}
    print("| task | condition | run | seed | rate Hz | CV | CV2 | pause | mode ms | dspk | min peak | 1-s n | wall h | d rate | rate ratio | d CV2 | sums | metrics |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in sorted(rows, key=lambda r: (r["root"], r["name"])):
        w = wt.get(r["seed"]) if not r["name"].startswith("WT_") else None
        d = (f"{r['rate_hz'] - w['rate_hz']:+.2f} | {r['rate_hz'] / w['rate_hz']:.3f} | {r['cv2'] - w['cv2']:+.4f}"
             if w else "- | - | -")
        print(f"| {r['task']} | {r['condition']} | {r['name']} | {r['seed'] - 20261013} | {r['rate_hz']:.2f} | {r['cv']:.4f} | {r['cv2']:.4f} | "
              f"{r['pause_frac']:.3f} | {r['log_isi_mode_ms']:.2f} | {int(r['dspikes'])} | {r['min_peak']:.2f} | "
              f"{r['c_min']}-{r['c_max']} | {r['wall_h']:.2f} | {d} | {r['sums']} | {r['metrics']} |")


if __name__ == "__main__":
    main()

#%% debug
