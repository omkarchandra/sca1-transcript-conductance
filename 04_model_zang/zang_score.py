"""Zang P4 scoring v02 (E1 definitions, 25-s runs) for run folders written by the Zang driver.

v02 = v01 with the window moved to the last 15 s of a 25-s run (pre-specified); PRIMARY used 30-45 s.
Window [10000, 25000] ms, both ends inclusive; rate = window spikes / 15 s; CV = sd(ISI)/mean(ISI), ddof 0;
CV2 = mean(2|dISI|/(ISI_i+ISI_i+1)); pause fraction = share of window ISIs > 3 x median ISI; log-ISI mode =
centre of the fullest 0.05-decade bin of log10(ISI ms) on [-1, 5] (values clipped to the range, half-open bins
except the inclusive rightmost edge, lowest bin wins ties). Dendritic spikes =
the RECORD's definition (upward crossings of -30 mV, 20 ms refractory, distal site) counted in the same window;
the count is re-derived from V_DISTAL.npy and must equal DSPIKES.csv over the whole run, else stop. The same
definition is applied to V_MAIN.npy (smooth-dendrite site) as a secondary count.
If a WT condition and other conditions are present, each run is paired with the WT run of the same noise seed.
Read-only; writes ONE new output folder.

python -I -B 04_model_zang/zang_score.py <RUN_ROOT> <OUT_DIR>
"""
import csv
import hashlib
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

T0_MS, TSTOP_MS = 10000.0, 25000.0
WINDOW_S = (TSTOP_MS - T0_MS) / 1000.0
EDGES = np.arange(-1.0, 5.025, 0.05)
KEYS = ("rate_hz", "cv", "cv2", "pause_frac", "log_isi_mode_ms", "dspikes", "dspike_rate_hz", "main_crossings")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def log_isi_mode(isi):
    values = [max(-1.0, min(float(EDGES[-1]), math.log10(x))) for x in isi]
    bins = np.zeros(len(EDGES) - 1, dtype=int)
    for x in values:
        bins[min(len(bins) - 1, int(np.searchsorted(EDGES, x, side="right")) - 1)] += 1
    k = int(np.argmax(bins))  # first maximum = lowest bin on ties
    return float(10 ** ((EDGES[k] + EDGES[k + 1]) / 2))


def spike_metrics(times_ms):
    t = np.asarray(times_ms, dtype=float)
    isi = np.diff(t)
    out = dict(n_spikes=int(t.size), rate_hz=float(t.size / WINDOW_S), cv=float("nan"), cv2=float("nan"),
               pause_frac=float("nan"), log_isi_mode_ms=float("nan"))
    if isi.size >= 1:
        med = float(np.median(isi))
        out["pause_frac"] = float(np.mean(isi > 3.0 * med))
        out["log_isi_mode_ms"] = log_isi_mode(isi)
    if isi.size >= 2:
        out["cv"] = float(np.std(isi, ddof=0) / np.mean(isi))
        out["cv2"] = float(np.mean(2.0 * np.abs(np.diff(isi)) / (isi[1:] + isi[:-1])))
    return out


def crossings(t, v, level, refract):
    idx = np.where((v[:-1] < level) & (v[1:] >= level))[0] + 1
    out, last = [], -1e9
    for i in idx:
        if t[i] - last >= refract:
            out.append(t[i])
            last = t[i]
    return np.array(out)


def read_times(path):
    with open(path, encoding="utf-8") as f:
        return np.array([float(r["time_ms"]) for r in csv.DictReader(f)])


def read_run(d):
    for line in (d / "SHA256SUMS").read_text().splitlines():
        h_, name = line.split("  ", 1)
        if sha256(d / name) != h_:
            raise SystemExit(f"checksum mismatch: {d / name}")
    rec = json.loads((d / "RECORD.json").read_text())
    t = read_times(d / "SPIKES.csv")
    m = spike_metrics(t[(t >= T0_MS) & (t <= TSTOP_MS)])
    m["n_spikes_total"] = int(t.size)
    dd = rec["dendritic_spike_definition"]
    dt = float(rec["record_dt_ms"])
    vd = np.load(d / "V_DISTAL.npy").astype(float)
    vm = np.load(d / "V_MAIN.npy").astype(float)
    tt = np.arange(len(vd)) * dt
    dts = crossings(tt, vd, float(dd["level_mV"]), float(dd["refractory_ms"]))
    if dts.size != read_times(d / "DSPIKES.csv").size or dts.size != int(rec["n_distal_dspikes"]):
        raise SystemExit(f"dendritic spike count differs from DSPIKES.csv/RECORD: {d}")
    mts = crossings(tt, vm, float(dd["level_mV"]), float(dd["refractory_ms"]))
    inwin = lambda x: int(np.count_nonzero((x >= T0_MS) & (x <= TSTOP_MS)))
    m.update(dspikes=inwin(dts), dspike_rate_hz=inwin(dts) / WINDOW_S, main_crossings=inwin(mts),
             dspikes_total=int(dts.size), wall_s=float(rec["wall_s"]))
    return rec, m


def q(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if not x.size:
        return "n/a"
    p = np.percentile(x, [25, 50, 75])
    return f"{p[1]:.4g} ({p[0]:.4g} to {p[2]:.4g})"


def main():
    if len(sys.argv) != 3:
        raise SystemExit("need RUN_ROOT OUT_DIR")
    root, out = Path(sys.argv[1]), Path(sys.argv[2])
    if out.exists():
        raise SystemExit(f"output exists; refusing to overwrite: {out}")
    folders = sorted(p for p in root.iterdir() if p.is_dir() and (p / "RECORD.json").is_file())
    if not folders:
        raise SystemExit(f"no run folders under {root}")
    rows = []
    for d in folders:
        rec, m = read_run(d)
        if float(rec["tstop_ms"]) < TSTOP_MS:
            raise SystemExit(f"run shorter than the window: {d}")
        sigma = float(rec["noise"]["sigma_nA"])
        rows.append(dict(folder=d.name, condition=rec["condition"], seed=rec["noise"]["seed"], sigma_nA=sigma,
                         group=rec["condition"] if sigma > 0 else f"{rec['condition']} (no input)", **m))
    wt = {r["seed"]: r for r in rows if r["group"] == "WT"}
    paired = bool(wt) and any(r["group"] != "WT" for r in rows)
    for r in rows:
        w = wt.get(r["seed"]) if paired else None
        if paired and r["group"] != "WT" and w is None:
            raise SystemExit(f"no WT run with seed {r['seed']} for {r['folder']}")
        for k in KEYS:
            r[f"d_{k}"] = (r[k] - w[k]) if w is not None else float("nan")
    out.mkdir(parents=True)
    cols = ["folder", "condition", "group", "seed", "sigma_nA", "n_spikes", "n_spikes_total", *KEYS, "dspikes_total",
            "wall_s", *[f"d_{k}" for k in KEYS]]
    with open(out / "METRICS.tsv", "x", encoding="utf-8", newline="") as f:
        w_ = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w_.writeheader()
        for r in rows:
            w_.writerow({k: r.get(k) for k in cols})
    L = ["# Zang 2018 P4 metrics, E1 definitions on 10-25 s inclusive (25-s runs; PRIMARY used 30-45 s)", "",
         "rate = window spikes / 15 s; CV ddof 0; CV2 adjacent-ISI; pause = share of ISIs > 3 x median; "
         "log-ISI mode = fullest 0.05-decade bin on [-1, 5]; dspikes = distal upward -30 mV crossings, "
         "20 ms refractory; main = same rule at the smooth-dendrite site."
         + (" d_* = run minus the WT run of the same seed." if paired else ""), "",
         "| condition | n | rate Hz | CV | CV2 | pause | log-ISI mode ms | dspikes | main crossings | wall h |",
         "|---|---:|---|---|---|---|---|---|---|---|"]
    conds = sorted({r["group"] for r in rows}, key=lambda c: (c != "WT", c))
    for c in conds:
        rr = [r for r in rows if r["group"] == c]
        L.append(f"| {c} | {len(rr)} | " + " | ".join(q([r[k] for r in rr]) for k in
                 ("rate_hz", "cv", "cv2", "pause_frac", "log_isi_mode_ms", "dspikes", "main_crossings"))
                 + f" | {q([r['wall_s'] / 3600 for r in rr])} |")
    if paired:
        L += ["", "Paired change against the WT run of the same seed, median (IQR); counts of seeds above WT.", "",
              "| condition | n | d rate Hz | d CV | d CV2 | d pause | d dspikes | CV > WT | CV2 > WT |",
              "|---|---:|---|---|---|---|---|---:|---:|"]
        for c in conds:
            if c == "WT":
                continue
            rr = [r for r in rows if r["group"] == c]
            L.append(f"| {c} | {len(rr)} | " + " | ".join(q([r[f'd_{k}'] for r in rr]) for k in
                     ("rate_hz", "cv", "cv2", "pause_frac", "dspikes"))
                     + f" | {sum(r['d_cv'] > 0 for r in rr)} | {sum(r['d_cv2'] > 0 for r in rr)} |")
    L += ["", "Per run.", "", "| folder | group | seed | spikes | rate Hz | CV | CV2 | pause | log-ISI mode ms "
          "| dspikes | main crossings | wall h |", "|---|---|---|---:|---|---|---|---|---|---:|---:|---|"]
    for r in rows:
        L.append(f"| {r['folder']} | {r['group']} | {r['seed']} | {r['n_spikes']} | " + " | ".join(
            f"{r[k]:.4g}" for k in ("rate_hz", "cv", "cv2", "pause_frac", "log_isi_mode_ms"))
            + f" | {r['dspikes']} | {r['main_crossings']} | {r['wall_s'] / 3600:.2f} |")
    (out / "SUMMARY.md").write_text("\n".join(L) + "\n")
    rec = dict(stage="Zang P4 scoring v02 (E1 10-25 s)", created_utc=datetime.now(timezone.utc).isoformat(
               timespec="seconds"), run_root=str(root), window=[T0_MS, TSTOP_MS], folders=len(rows),
               conditions=conds, paired=paired, script_sha256=sha256(Path(__file__)),
               outputs={n: sha256(out / n) for n in ("METRICS.tsv", "SUMMARY.md")})
    (out / "RECORD.json").write_text(json.dumps(rec, indent=1) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()

#%% debug
