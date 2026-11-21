"""E7 scoring v04 (E1 definitions): rate, CV, CV2, pause fraction in the 30-45 s window.

Window [30000, 45000] ms, both ends inclusive; rate = spikes in window / 15 s; CV = sd(ISI)/mean(ISI)
with ddof 0 over the window ISIs; CV2 = mean(2|dISI|/(ISI_i+ISI_i+1)); pause fraction = share of
window ISIs > 3 x median ISI. Every run is paired with its member's PRIMARY WT run and the noise
seed must match, else stop. Read-only; writes ONE new output folder.

Usage:
  python -I -B e7_score_20261024_v04_e1.py <RUN_ROOT> [<RUN_ROOT> ...] <WT_ROOT> <OUT_DIR>
"""
import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

T0_MS, TSTOP_MS = 30000.0, 45000.0
WINDOW_S = (TSTOP_MS - T0_MS) / 1000.0


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def metrics(times_ms):
    t = np.asarray(times_ms, dtype=float)
    isi = np.diff(t)
    out = dict(n_spikes=int(t.size), rate_hz=float(t.size / WINDOW_S))
    if isi.size >= 2:
        out.update(cv=float(np.std(isi, ddof=0) / np.mean(isi)),
                   cv2=float(np.mean(2.0 * np.abs(np.diff(isi)) / (isi[1:] + isi[:-1]))),
                   mean_isi_ms=float(np.mean(isi)))
        med = float(np.median(isi))
        out["pause_frac"] = float(np.mean(isi > 3.0 * med)) if med > 0 else float("nan")
    else:
        out.update(cv=float("nan"), cv2=float("nan"), pause_frac=float("nan"),
                   mean_isi_ms=float("nan"))
    return out


def read_run(d):
    rec = json.loads((d / "RECORD.json").read_text())
    with open(d / "SPIKES.csv", encoding="utf-8") as f:
        t = np.array([float(x["time_ms"]) for x in csv.DictReader(f)])
    window = t[(t >= T0_MS) & (t <= TSTOP_MS)]
    m = metrics(window)
    m["n_spikes_total"] = int(t.size)
    return rec, m


def main():
    if len(sys.argv) < 4:
        raise SystemExit("need RUN_ROOT [RUN_ROOT...] WT_ROOT OUT_DIR")
    roots = [Path(a) for a in sys.argv[1:-2]]
    wt_root, out = Path(sys.argv[-2]), Path(sys.argv[-1])
    if out.exists():
        raise SystemExit(f"output exists; refusing to overwrite: {out}")
    rows, seed_checks = [], 0
    for root in roots:
        folders = sorted(p for p in root.iterdir() if p.is_dir() and p.name.startswith("m"))
        if not folders:
            raise SystemExit(f"no member x condition folders under {root}")
        for d in folders:
            rec, m = read_run(d)
            member = int(rec["member_id"])
            wt_dir = wt_root / f"m{member:03d}_c00"
            if not (wt_dir / "RECORD.json").is_file():
                raise SystemExit(f"PRIMARY WT run missing: {wt_dir}")
            wrec, wm = read_run(wt_dir)
            if wrec.get("condition") != "WT" or int(wrec["member_id"]) != member:
                raise SystemExit(f"WT reference mismatch: {wt_dir}")
            if int(rec["noise"]["seed"]) != int(wrec["noise"]["seed"]):
                raise SystemExit(f"seed mismatch: {d}")
            if float(rec["run_s"]) != float(wrec["run_s"]):
                raise SystemExit(f"run_s mismatch: {d}")
            seed_checks += 1
            row = dict(root=root.name, folder=d.name, condition=rec.get("condition"),
                       member_id=member, seed=int(rec["noise"]["seed"]),
                       cdp5_scale=rec.get("cdp5_scale"), window_s=WINDOW_S, **m)
            for k in ("rate_hz", "cv", "cv2", "pause_frac", "mean_isi_ms", "n_spikes"):
                row[f"wt_{k}"] = wm[k]
            for k in ("rate_hz", "cv", "cv2", "pause_frac"):
                row[f"d_{k}"] = m[k] - wm[k]
            rows.append(row)
    out.mkdir(parents=True)
    cols = ["root", "folder", "condition", "member_id", "seed", "cdp5_scale", "window_s",
            "n_spikes", "rate_hz", "cv", "cv2", "pause_frac", "mean_isi_ms",
            "wt_rate_hz", "wt_cv", "wt_cv2", "wt_pause_frac", "wt_mean_isi_ms", "wt_n_spikes",
            "d_rate_hz", "d_cv", "d_cv2", "d_pause_frac"]
    with open(out / "E7_METRICS.tsv", "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in cols})
    conds = sorted({r["condition"] for r in rows})
    L = ["# E7 metrics, E1 definitions (30-45 s inclusive), paired against each member's PRIMARY WT", "",
         "rate = window spikes / 15 s; CV = sd(ISI)/mean(ISI), ddof 0; CV2 = mean(2|dISI|/(ISI_i+ISI_i+1)); "
         "pause fraction = share of window ISIs > 3 x median ISI. d_* = condition minus that member's WT.", "",
         "| condition | n | rate Hz median (IQR) | d_rate Hz | CV median | d_CV | CV2 median | d_CV2 | pause median | d_pause |",
         "|---|---:|---|---|---|---|---|---|---|---|"]
    for c in conds:
        rr = [r for r in rows if r["condition"] == c]
        def iqr(key):
            x = np.array([r[key] for r in rr], dtype=float)
            x = x[np.isfinite(x)]
            if not x.size:
                return "n/a"
            q = np.percentile(x, [25, 50, 75])
            return f"{q[1]:.3g} ({q[0]:.3g} to {q[2]:.3g})"
        L.append(f"| {c} | {len(rr)} | {iqr('rate_hz')} | {iqr('d_rate_hz')} | {iqr('cv')} | "
                 f"{iqr('d_cv')} | {iqr('cv2')} | {iqr('d_cv2')} | {iqr('pause_frac')} | {iqr('d_pause_frac')} |")
    (out / "SUMMARY.md").write_text("\n".join(L) + "\n")
    rec = dict(stage="E7 scoring v04 (E1 30-45 s, paired)",
               created_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
               run_roots=[str(r) for r in roots], wt_root=str(wt_root),
               window=[T0_MS, TSTOP_MS], conditions=conds, folders=len(rows),
               seed_checks_passed=seed_checks, script_sha256=sha256(Path(__file__)),
               outputs={n: sha256(out / n) for n in ("E7_METRICS.tsv", "SUMMARY.md")})
    (out / "RECORD.json").write_text(json.dumps(rec, indent=1) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()

#%% debug
