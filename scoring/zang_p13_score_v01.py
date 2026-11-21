"""Zang P13 scoring v01 = zang_score_20261026_v02.py (sha256 76c6cc89...) plus a window parameter and two new modes
(stage-1 screening, stage-2 population reading R1-R5 of P13_zang_population_prespecification_20261109.md).
Metric code (spike_metrics, log_isi_mode, crossings, read_run) is unchanged from v02; new code only wraps it.
Default mode (no --mode) reproduces v02: window [10000, 25000] ms, outputs METRICS.tsv, SUMMARY.md, RECORD.json.

Zang P4 scoring v02 (E1 definitions, 25-s runs) for run folders written by zang_driver_20261026_v01.py.

v02 = v01 with the window moved to the last 15 s of a 25-s run; PRIMARY used 30-45 s.
Window [10000, 25000] ms, both ends inclusive; rate = window spikes / 15 s; CV = sd(ISI)/mean(ISI), ddof 0;
CV2 = mean(2|dISI|/(ISI_i+ISI_i+1)); pause fraction = share of window ISIs > 3 x median ISI; log-ISI mode =
centre of the fullest 0.05-decade bin of log10(ISI ms) on [-1, 5] (values clipped to the range, half-open bins
except the inclusive rightmost edge, lowest bin wins ties). Dendritic spikes =
the RECORD's definition (upward crossings of -30 mV, 20 ms refractory, distal site) counted in the same window;
the count is re-derived from V_DISTAL.npy and must equal DSPIKES.csv over the whole run, else stop. The same
definition is applied to V_MAIN.npy (smooth-dendrite site) as a secondary count.
If a WT condition and other conditions are present, each run is paired with the WT run of the same noise seed.
Read-only; writes ONE new output folder.

Usage: python -I -B zang_p13_score_v01.py <RUN_ROOT> <OUT_DIR> [--mode legacy|screen|population]
           [--window-ms T0 TSTOP] [--n-sets N] [--kept-list SCREEN.tsv]
  legacy (default)  = v02 behaviour, window 10000-25000 ms.
  screen            = stage 1: set{NN}_WT runs, window default 5000-10000 ms; writes SCREEN.tsv, SCORE_RECORD.json.
  population        = stage 2: set{NN}_{COND} runs, window default 10000-25000 ms; writes RUNS.tsv, PER_SET.tsv,
                      READING.json, SCORE_RECORD.json.
"""
import argparse
import csv
import hashlib
import json
import math
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

T0_MS, TSTOP_MS = 10000.0, 25000.0
WINDOW_S = (TSTOP_MS - T0_MS) / 1000.0
SCRIPT_NAME = "zang_p13_score_v01.py"
BASE_SCRIPT_SHA256_PREFIX = "ed3c7859"  # zang_score_20261026_v02.py


def set_window(t0_ms, tstop_ms):
    """Rebind the three window globals together (spike_metrics / read_run read them at call time)."""
    global T0_MS, TSTOP_MS, WINDOW_S
    if not tstop_ms > t0_ms:
        raise SystemExit("window must have tstop > t0")
    T0_MS, TSTOP_MS = float(t0_ms), float(tstop_ms)
    WINDOW_S = (TSTOP_MS - T0_MS) / 1000.0


# --- P13 constants (pre-specification 9 Nov 2026) ---
CONDITIONS = ("WT", "KI_18w_mid", "KI_24w_mid", "KI_30w_mid", "KI_18w_mid_A0.6", "five_aux")
AGE_CONDITIONS = ("KI_18w_mid", "KI_24w_mid", "KI_30w_mid")
ATROPHY_CONDITION = "KI_18w_mid_A0.6"
AUX_CONDITION = "five_aux"
REQUIRED_CONDITIONS = ("WT",) + AGE_CONDITIONS  # the others may be dropped for cost (section 4)
SCREEN_RATE_HZ = (20.0, 120.0)  # inclusive
SCREEN_CV_MAX = 0.2  # strict
PLATEAU_LEVEL_MV = -30.0  # soma voltage strictly above this ...
PLATEAU_MIN_S = 1.0  # ... continuously for strictly more than this
SCREEN_MIN_KEPT = 5
SCREEN_WINDOW_MS = (5000.0, 10000.0)
POP_WINDOW_MS = (10000.0, 25000.0)
RATIO_PHENOTYPE = (0.90, 1.10)  # R1 inclusive band
R2_RATE_RATIO_GT = 1.10
R2_FRACTION = (2, 3)  # at least 2/3 of kept sets
R1_MIN_SETS_ABS, R1_MIN_SETS_FRAC = 3, 0.20
R4_PRODUCT_FACTOR = 1.42
BURST_CV = 0.5  # R5 exploratory proxy for "bursting" (window ISI CV >= this); NOT in the pre-spec
EXPECT_SIGMA_NA = 0.0065
SEED_BASE = 20261109  # seed = SEED_BASE + set number (warning only)
EPS = 1e-12
INPUT_SHA = {}  # path -> sha256 of every input file verified in read_run
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
        INPUT_SHA[str(d / name)] = h_
    INPUT_SHA[str(d / "SHA256SUMS")] = sha256(d / "SHA256SUMS")
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


def main_legacy(root, out):
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
    L = [f"# Zang 2018 P4 metrics, E1 definitions on {T0_MS / 1000:g}-{TSTOP_MS / 1000:g} s inclusive "
         "(25-s runs; PRIMARY used 30-45 s)", "",
         f"rate = window spikes / {WINDOW_S:g} s; CV ddof 0; CV2 adjacent-ISI; pause = share of ISIs > 3 x median; "
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
    rec = dict(stage=f"Zang P4 scoring v02 (E1 {T0_MS / 1000:g}-{TSTOP_MS / 1000:g} s), run by {SCRIPT_NAME}", created_utc=datetime.now(timezone.utc).isoformat(
               timespec="seconds"), run_root=str(root), window=[T0_MS, TSTOP_MS], folders=len(rows),
               conditions=conds, paired=paired, script_sha256=sha256(Path(__file__)),
               outputs={n: sha256(out / n) for n in ("METRICS.tsv", "SUMMARY.md")})
    (out / "RECORD.json").write_text(json.dumps(rec, indent=1) + "\n")
    print("\n".join(L))
    return rows


# ====================== P13 additions (new code only; v02 metric code above is untouched) ======================


def require_isolated():
    if not (sys.flags.isolated and sys.flags.dont_write_bytecode):
        raise SystemExit("screen/population modes must be run as: python -I -B " + SCRIPT_NAME + " ...")


def scan_folders(root):
    """-> ({(set_id, cond_suffix): path}, skipped_names). Every directory named set{NN}_{suffix} is returned,
    with or without RECORD.json (a missing RECORD.json is reported by the caller, never silently dropped)."""
    runs, skipped = {}, []
    for p in sorted(root.iterdir()):
        m = re.fullmatch(r"set(\d+)_(.+)", p.name)
        if not p.is_dir() or not m:
            skipped.append(p.name)
            continue
        key = (int(m.group(1)), m.group(2))
        if key in runs:
            raise SystemExit(f"two folders map to set {key[0]} condition {key[1]}: {p.name}")
        runs[key] = p
    return runs, skipped


def longest_above_s(d, dt_ms, level=PLATEAU_LEVEL_MV):
    """Longest continuous stretch with V_SOMA strictly above `level`, over the WHOLE run; L consecutive
    samples count as L*dt (each sample represents one record_dt)."""
    p = d / "V_SOMA.npy"
    if str(p) not in INPUT_SHA:
        raise SystemExit(f"V_SOMA.npy not covered by SHA256SUMS: {d}")
    v = np.load(p).astype(float)
    x = np.concatenate(([0], (v > level).astype(np.int8), [0]))
    dx = np.diff(x)
    starts, ends = np.where(dx == 1)[0], np.where(dx == -1)[0]
    n = int((ends - starts).max()) if starts.size else 0
    return n * dt_ms / 1000.0


def score_run(d, set_id, cond):
    """v02 metrics (read_run) for the current window plus the new V_SOMA plateau length and pause count."""
    rec, m = read_run(d)
    if float(rec["tstop_ms"]) < TSTOP_MS:
        raise SystemExit(f"run shorter than the window: {d}")
    n_isi = max(m["n_spikes"] - 1, 0)
    pause_count = int(round(m["pause_frac"] * n_isi)) if n_isi else 0  # share of window ISIs > 3 x median, as a count
    long_s = longest_above_s(d, float(rec["record_dt_ms"]))
    sigma = float(rec["noise"]["sigma_nA"])
    return dict(set_id=set_id, condition=cond, folder=d.name, record_condition=rec.get("condition"),
                seed=rec["noise"]["seed"], sigma_nA=sigma, n_isi=n_isi, n_pauses=pause_count,
                longest_above_m30_s=long_s, plateau=int(long_s > PLATEAU_MIN_S), **m)


def write_tsv(path, cols, rows):
    with open(path, "x", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def jnum(x):
    """JSON-safe number: NaN/inf -> None."""
    if x is None:
        return None
    x = float(x)
    return x if math.isfinite(x) else None


def write_score_record(out, root, mode, args, counts, warnings, extra=None):
    outputs = {p.name: sha256(p) for p in sorted(out.iterdir()) if p.is_file()}
    rec = dict(stage=f"Zang P13 scoring v01 ({mode})", created_utc=datetime.now(timezone.utc).isoformat(
               timespec="seconds"), run_root=str(root), argv=sys.argv[1:], window_ms=[T0_MS, TSTOP_MS],
               script=str(Path(__file__).resolve()), script_sha256=sha256(Path(__file__)),
               base_script=f"zang_score_20261026_v02.py (sha256 {BASE_SCRIPT_SHA256_PREFIX}...)",
               python=sys.version.split()[0], numpy=np.__version__,
               flags=dict(isolated=bool(sys.flags.isolated), dont_write_bytecode=bool(sys.flags.dont_write_bytecode)),
               counts=counts, warnings=warnings, input_sha256=dict(sorted(INPUT_SHA.items())),
               n_input_files=len(INPUT_SHA), outputs_sha256=outputs, **(extra or {}))
    with open(out / "SCORE_RECORD.json", "x", encoding="utf-8", newline="") as f:
        f.write(json.dumps(rec, indent=1) + "\n")


def main_screen(root, out, n_sets):
    require_isolated()
    if out.exists():
        raise SystemExit(f"output exists; refusing to overwrite: {out}")
    runs, skipped = scan_folders(root)
    warnings = [f"skipped (not a set{{NN}}_{{COND}} folder): {n}" for n in skipped]
    wt_runs = {k[0]: d for k, d in runs.items() if k[1] == "WT"}
    for k in sorted(runs):
        if k[1] != "WT":
            warnings.append(f"ignored non-WT folder in screening root: {runs[k].name}")
    ids = sorted(set(wt_runs) | (set(range(1, n_sets + 1)) if n_sets else set()))
    if not ids:
        raise SystemExit(f"no set*_WT folders under {root}")
    rows = []
    for sid in ids:
        d = wt_runs.get(sid)
        if d is None:
            rows.append(dict(set_id=sid, seed="NA", rate="NA", isi_cv="NA", longest_above_m30_s="NA", kept=0,
                             reason="no run folder"))
            continue
        if not (d / "RECORD.json").is_file():
            rows.append(dict(set_id=sid, seed="NA", rate="NA", isi_cv="NA", longest_above_m30_s="NA", kept=0,
                             reason="no RECORD.json (failed or incomplete run)"))
            continue
        r = score_run(d, sid, "WT")
        why = []
        if not r["rate_hz"] >= SCREEN_RATE_HZ[0]:
            why.append(f"rate<{SCREEN_RATE_HZ[0]:g}")
        if r["rate_hz"] > SCREEN_RATE_HZ[1]:
            why.append(f"rate>{SCREEN_RATE_HZ[1]:g}")
        if not math.isfinite(r["cv"]):
            why.append("cv undefined (<3 spikes in window)")
        elif not r["cv"] < SCREEN_CV_MAX:
            why.append(f"cv>={SCREEN_CV_MAX:g}")
        if r["plateau"]:
            why.append(f"plateau (>{PLATEAU_LEVEL_MV:g} mV for {r['longest_above_m30_s']:.4g} s > {PLATEAU_MIN_S:g} s)")
        if abs(r["sigma_nA"] - EXPECT_SIGMA_NA) > 1e-12:
            warnings.append(f"set {sid}: sigma_nA {r['sigma_nA']} != {EXPECT_SIGMA_NA}")
        if r["seed"] != SEED_BASE + sid:
            warnings.append(f"set {sid}: seed {r['seed']} != {SEED_BASE + sid}")
        rows.append(dict(set_id=sid, seed=r["seed"], rate=r["rate_hz"], isi_cv=r["cv"],
                         longest_above_m30_s=r["longest_above_m30_s"], kept=0 if why else 1,
                         reason="; ".join(why) if why else "ok"))
    out.mkdir(parents=True)
    write_tsv(out / "SCREEN.tsv", ["set_id", "seed", "rate", "isi_cv", "longest_above_m30_s", "kept", "reason"], rows)
    kept = [r["set_id"] for r in rows if r["kept"]]
    crit = dict(rate_below=0, rate_above=0, cv_undefined=0, cv_too_high=0, plateau=0, no_run=0)
    for r in rows:
        t = r["reason"]
        crit["rate_below"] += "rate<" in t
        crit["rate_above"] += "rate>" in t
        crit["cv_undefined"] += "cv undefined" in t
        crit["cv_too_high"] += "cv>=" in t
        crit["plateau"] += "plateau" in t
        crit["no_run"] += t.startswith("no ")
    counts = dict(
        n_sets=len(rows), n_kept=len(kept), n_rejected=len(rows) - len(kept), kept_set_ids=kept,
        criterion_hits_overlapping={k: int(v) for k, v in crit.items()},
        screen_rule=dict(rate_hz_inclusive=SCREEN_RATE_HZ, cv_strictly_below=SCREEN_CV_MAX,
                         plateau_above_mV=PLATEAU_LEVEL_MV, plateau_longer_than_s=PLATEAU_MIN_S),
        fewer_than_min_kept=len(kept) < SCREEN_MIN_KEPT, min_kept=SCREEN_MIN_KEPT)
    write_score_record(out, root, "screen", sys.argv[1:], counts, warnings)
    print(f"screen: {len(rows)} sets, kept {len(kept)}, rejected {len(rows) - len(kept)}; kept ids {kept}")
    for w in warnings:
        print("warning:", w)
    if len(kept) < SCREEN_MIN_KEPT:
        print(f"STOP: fewer than {SCREEN_MIN_KEPT} sets kept; report to the owner (pre-spec section 4).")
        sys.exit(3)


# ---------------------------------------------- stage 2 ----------------------------------------------
def med(x):
    a = np.asarray([v for v in x], dtype=float)
    a = a[np.isfinite(a)]
    return (float(np.median(a)) if a.size else None), int(a.size)


def cond_stats(cond, sets, by, kept):
    present = [s for s in sets if (s, cond) in by]
    rr = [by[(s, cond)] for s in present]
    ratio = [r["rate_ratio"] for r in rr]
    dcv2 = [r["d_cv2"] for r in rr]
    mr, nr = med(ratio)
    md, nd = med(dcv2)
    ids = lambda f: [r["set_id"] for r in rr if f(r)]
    pause_ids = ids(lambda r: r["n_pauses"] >= 1)
    dsp_ids = ids(lambda r: r["dspikes"] >= 1)
    pd_ids = ids(lambda r: r["n_pauses"] >= 1 or r["dspikes"] >= 1)
    gt110 = ids(lambda r: r["rate_ratio"] > R2_RATE_RATIO_GT + EPS)
    cv2pos = ids(lambda r: r["d_cv2"] > 0)
    slower = ids(lambda r: r["rate_ratio"] < 1.0 - EPS)
    thr = max(R1_MIN_SETS_ABS, (kept * 20 + 99) // 100)  # max(3, ceil(0.2 * kept)), exact integer arithmetic
    lo, hi = RATIO_PHENOTYPE
    r1 = dict(median_d_cv2=md, n_valid_d_cv2=nd, median_rate_ratio=mr, n_valid_rate_ratio=nr,
              n_pause_sets=len(pause_ids), pause_set_ids=pause_ids, pause_set_threshold=thr,
              c_median_d_cv2_gt_0=bool(md is not None and md > 0),
              c_median_ratio_in_band=bool(mr is not None and lo - EPS <= mr <= hi + EPS),
              c_pause_sets_ge_threshold=len(pause_ids) >= thr)
    r1["R1_moves_toward_phenotype"] = bool(r1["c_median_d_cv2_gt_0"] and r1["c_median_ratio_in_band"]
                                           and r1["c_pause_sets_ge_threshold"])
    pd_ok = len(pd_ids) < 3 and 5 * len(pd_ids) < kept  # fewer than 3 AND < 20 % of kept sets
    r2 = dict(n_rate_ratio_gt_1_10=len(gt110), rate_ratio_gt_1_10_set_ids=gt110,
              n_d_cv2_gt_0=len(cv2pos), d_cv2_gt_0_set_ids=cv2pos,
              n_pause_or_dspike_sets=len(pd_ids), pause_or_dspike_set_ids=pd_ids,
              n_dspike_sets=len(dsp_ids), dspike_set_ids=dsp_ids,
              c_rate_ratio_gt_1_10_ge_2_3=bool(R2_FRACTION[1] * len(gt110) >= R2_FRACTION[0] * kept),
              c_d_cv2_gt_0_ge_2_3=bool(R2_FRACTION[1] * len(cv2pos) >= R2_FRACTION[0] * kept),
              c_pause_or_dspike_fewer_than_3_and_lt_20pct=bool(pd_ok), withdrawal_trigger=bool(not pd_ok))
    r2["R2_representative"] = bool(r2["c_rate_ratio_gt_1_10_ge_2_3"] and r2["c_d_cv2_gt_0_ge_2_3"] and pd_ok)
    r5 = dict(silent_runs=ids(lambda r: r["n_spikes"] == 0), plateau_runs=ids(lambda r: r["plateau"] == 1),
              bursting_runs_cv_ge_proxy=ids(lambda r: math.isfinite(r["cv"]) and r["cv"] >= BURST_CV),
              pause_runs=pause_ids, dspike_runs=dsp_ids)
    r5.update({"n_" + k: len(v) for k, v in list(r5.items())})
    return dict(n_present=len(present), n_missing=kept - len(present), n_slower_than_wt=len(slower),
                slower_than_wt_set_ids=slower, R1=r1, R2_inputs=r2, R5=r5)


def main_population(root, out, kept_list):
    require_isolated()
    if out.exists():
        raise SystemExit(f"output exists; refusing to overwrite: {out}")
    runs, skipped = scan_folders(root)
    warnings = [f"skipped (not a set{{NN}}_{{COND}} folder): {n}" for n in skipped]
    unknown = sorted({c for _, c in runs if c not in CONDITIONS})
    if unknown:
        raise SystemExit(f"unknown condition suffix(es) {unknown}; allowed {list(CONDITIONS)}")
    incomplete = [d.name for k, d in sorted(runs.items()) if not (d / "RECORD.json").is_file()]
    for n in incomplete:
        warnings.append(f"no RECORD.json (failed or incomplete run), treated as missing: {n}")
    runs = {k: d for k, d in runs.items() if (d / "RECORD.json").is_file()}
    sets = sorted({k[0] for k in runs})
    for s_ in sets:
        if (s_, "WT") not in runs:
            raise SystemExit(f"set {s_} has runs but no WT run; cannot compare")
    if not sets:
        raise SystemExit(f"no set*_COND run folders under {root}")
    if kept_list is not None:
        with open(kept_list, encoding="utf-8") as f:
            ids = {int(r["set_id"]) for r in csv.DictReader(f, delimiter="\t") if r["kept"] == "1"}
        INPUT_SHA[str(kept_list)] = sha256(kept_list)
        if ids != set(sets):
            raise SystemExit(f"sets under {root} {sets} differ from kept sets in {kept_list}: {sorted(ids)}")
    kept = len(sets)
    conds_present = [c for c in CONDITIONS if any((s, c) in runs for s in sets)]
    for c in REQUIRED_CONDITIONS:
        miss = [s for s in sets if (s, c) not in runs]
        if miss:
            raise SystemExit(f"required condition {c} missing for sets {miss}")
    dropped = [c for c in CONDITIONS if c not in conds_present]
    for c in conds_present:
        miss = [s for s in sets if (s, c) not in runs]
        if miss:
            warnings.append(f"condition {c} missing for sets {miss}")
    rows = []
    for (sid, cond), d in sorted(runs.items()):
        r = score_run(d, sid, cond)
        if abs(r["sigma_nA"] - EXPECT_SIGMA_NA) > 1e-12:
            warnings.append(f"{d.name}: sigma_nA {r['sigma_nA']} != {EXPECT_SIGMA_NA}")
        if r["seed"] != SEED_BASE + sid:
            warnings.append(f"{d.name}: seed {r['seed']} != {SEED_BASE + sid}")
        r["silent"] = int(r["n_spikes"] == 0)
        r["bursting_proxy"] = int(math.isfinite(r["cv"]) and r["cv"] >= BURST_CV)
        rows.append(r)
    for sid in sets:
        seeds = {r["seed"] for r in rows if r["set_id"] == sid}
        if len(seeds) != 1:
            raise SystemExit(f"set {sid}: conditions do not share one noise seed: {sorted(seeds)}")
    by = {(r["set_id"], r["condition"]): r for r in rows}
    for r in rows:
        w = by[(r["set_id"], "WT")]
        r["wt_rate_hz"], r["wt_cv2"] = w["rate_hz"], w["cv2"]
        r["rate_ratio"] = r["rate_hz"] / w["rate_hz"] if w["rate_hz"] > 0 else float("nan")
        r["d_cv2"] = r["cv2"] - w["cv2"]
        r["pause_or_dspike"] = int(r["n_pauses"] >= 1 or r["dspikes"] >= 1)
        r["slower_than_wt"] = int(r["rate_ratio"] < 1.0 - EPS)
    out.mkdir(parents=True)
    write_tsv(out / "RUNS.tsv", ["set_id", "condition", "folder", "record_condition", "seed", "sigma_nA", "n_spikes",
              "n_spikes_total", "n_isi", *KEYS, "n_pauses", "dspikes_total", "longest_above_m30_s", "plateau",
              "silent", "bursting_proxy", "wall_s", "wt_rate_hz", "wt_cv2", "rate_ratio", "d_cv2"], rows)
    write_tsv(out / "PER_SET.tsv", ["set_id", "seed", "condition", "wt_rate_hz", "rate_hz", "rate_ratio", "wt_cv2",
              "cv2", "d_cv2", "n_pauses", "dspikes", "pause_or_dspike", "slower_than_wt"],
              [r for r in rows if r["condition"] != "WT"])
    ki = [c for c in conds_present if c != "WT"]
    stats = {c: cond_stats(c, sets, by, kept) for c in ki}
    # R2 (per KI age) and withdrawal (any KI condition)
    r2_age = {}
    for c in AGE_CONDITIONS:
        i = stats[c]["R2_inputs"]
        verdict = ("representative" if i["R2_representative"] else
                   "set-dependent (withdrawal trigger met under this condition)" if i["withdrawal_trigger"] else
                   "partly representative")
        r2_age[c] = dict(inputs=i, verdict=verdict)
    trig_all = [c for c in ki if stats[c]["R2_inputs"]["withdrawal_trigger"]]
    trig_ages = [c for c in AGE_CONDITIONS if stats[c]["R2_inputs"]["withdrawal_trigger"]]
    R2 = dict(per_age=r2_age, withdrawal_triggered_by_any_KI_condition=bool(trig_all),
              withdrawal_trigger_conditions_all_KI=trig_all, withdrawal_triggered_by_three_ages_only=bool(trig_ages),
              withdrawal_trigger_conditions_ages=trig_ages,
              statement_if_withdrawn='replace "no run paused and no dendritic calcium spike occurred in the Zang '
                                     'model" with "set-dependent", with the counts',
              verdict_keyed_to="any KI condition present (literal reading); three-ages-only value given for comparison",
              overall="statement withdrawn: set-dependent" if trig_all else "not withdrawn")
    # R3
    R3 = None
    if ATROPHY_CONDITION in ki:
        both = [s for s in sets if (s, ATROPHY_CONDITION) in by and (s, "KI_18w_mid") in by]
        ma, _ = med([by[(s, ATROPHY_CONDITION)]["rate_ratio"] for s in both])
        mk, _ = med([by[(s, "KI_18w_mid")]["rate_ratio"] for s in both])
        da, _ = med([by[(s, ATROPHY_CONDITION)]["d_cv2"] for s in both])
        dk, _ = med([by[(s, "KI_18w_mid")]["d_cv2"] for s in both])
        pr, _ = med([by[(s, ATROPHY_CONDITION)]["rate_ratio"] - by[(s, "KI_18w_mid")]["rate_ratio"] for s in both])
        pd_, _ = med([by[(s, ATROPHY_CONDITION)]["d_cv2"] - by[(s, "KI_18w_mid")]["d_cv2"] for s in both])
        R3 = dict(sets_with_both=both, n_sets_with_both=len(both),
                  median_rate_ratio_A0_6=ma, median_rate_ratio_KI_18w=mk, median_d_cv2_A0_6=da,
                  median_d_cv2_KI_18w=dk, median_paired_diff_rate_ratio=pr, median_paired_diff_d_cv2=pd_,
                  rate_ratio_moves_toward_1=bool(ma is not None and mk is not None and abs(ma - 1) < abs(mk - 1)),
                  d_cv2_moves_toward_0=bool(da is not None and dk is not None and abs(da) < abs(dk)),
                  n_sets_slower_than_WT_A0_6=stats[ATROPHY_CONDITION]["n_slower_than_wt"],
                  slower_than_WT_set_ids_A0_6=stats[ATROPHY_CONDITION]["slower_than_wt_set_ids"],
                  geometry_reverses_sign=bool(ma is not None and ma < 1.0))
    # R4
    R4 = None
    if AUX_CONDITION in ki:
        both = [s for s in sets if (s, AUX_CONDITION) in by and (s, "KI_18w_mid") in by]
        ma, _ = med([by[(s, AUX_CONDITION)]["rate_ratio"] for s in both])
        mk, _ = med([by[(s, "KI_18w_mid")]["rate_ratio"] for s in both])
        pp, npp = med([by[(s, AUX_CONDITION)]["rate_ratio"] / (by[(s, "KI_18w_mid")]["rate_ratio"] * R4_PRODUCT_FACTOR)
                       if by[(s, "KI_18w_mid")]["rate_ratio"] > 0 else float("nan") for s in both])
        pr, _ = med([by[(s, AUX_CONDITION)]["rate_ratio"] / by[(s, "KI_18w_mid")]["rate_ratio"]
                     if by[(s, "KI_18w_mid")]["rate_ratio"] > 0 else float("nan") for s in both])
        r1a, r1k = stats[AUX_CONDITION]["R1"], stats["KI_18w_mid"]["R1"]
        R4 = dict(R1_five_aux=r1a, R1_KI_18w_mid_moves_toward_phenotype=r1k["R1_moves_toward_phenotype"],
                  five_aux_moves_toward_phenotype_where_KI_18w_did_not=bool(
                      r1a["R1_moves_toward_phenotype"] and not r1k["R1_moves_toward_phenotype"]),
                  product_check_descriptive=dict(factor=R4_PRODUCT_FACTOR, sets_with_both=both,
                      median_rate_ratio_five_aux=ma, median_rate_ratio_KI_18w=mk,
                      median_rate_ratio_KI_18w_times_factor=(mk * R4_PRODUCT_FACTOR if mk is not None else None),
                      median_per_set_five_aux_over_KI_18w=pr,
                      median_per_set_five_aux_over_KI_18w_times_factor=pp, n_valid=npp))
    R5 = {"definitions": dict(silent="no soma spike in the scoring window", plateau=f"V_SOMA > {PLATEAU_LEVEL_MV:g} mV "
          f"continuously for > {PLATEAU_MIN_S:g} s anywhere in the run",
          bursting_proxy=f"window ISI CV >= {BURST_CV:g} (proxy; the pre-spec gives no definition)",
          pause="at least one window ISI > 3 x median", dspike="at least one distal dendritic spike in the window"),
          "WT": {k: v for k, v in cond_stats("WT", sets, by, kept)["R5"].items()}, **{c: stats[c]["R5"] for c in ki}}
    reading = dict(
        prespec="P13_zang_population_prespecification_20261109.md sections 2-5", window_ms=[T0_MS, TSTOP_MS],
        n_kept=kept, kept_set_ids=sets, conditions_present=conds_present, conditions_absent=dropped,
        constants=dict(RATIO_PHENOTYPE=RATIO_PHENOTYPE, R2_RATE_RATIO_GT=R2_RATE_RATIO_GT, R2_FRACTION=R2_FRACTION,
                       R1_MIN_SETS_ABS=R1_MIN_SETS_ABS, R1_MIN_SETS_FRAC=R1_MIN_SETS_FRAC,
                       R4_PRODUCT_FACTOR=R4_PRODUCT_FACTOR, BURST_CV=BURST_CV, EPS=EPS),
        interpretation=dict(
            R1="n_pause_sets >= max(3, ceil(0.2*kept)); pause set = >=1 window ISI > 3 x median in the KI run; "
               "medians over finite values (n_valid reported)",
            R2="representative: >=2/3 of kept sets rate ratio > 1.10 AND >=2/3 d_CV2 > 0 AND n_pause_or_dspike < 3 AND "
               "n/kept < 0.20 (both); withdrawal trigger = the complement (n >= 3 OR n/kept >= 0.20), evaluated per "
               "condition and keyed to any KI condition present; absent/NaN counts as not satisfying",
            R3="medians over sets holding both conditions; toward 1 / toward 0 = smaller |median - target| than KI_18w_mid",
            R4="R1 on five_aux; product check descriptive only"),
        per_condition=stats, R2=R2, R3=R3, R4=R4, R5=R5)
    with open(out / "READING.json", "x", encoding="utf-8", newline="") as f:
        f.write(json.dumps(reading, indent=1) + "\n")
    counts = dict(n_runs=len(rows), n_kept_sets=kept, conditions_present=conds_present, conditions_absent=dropped,
                  runs_per_condition={c: sum(r["condition"] == c for r in rows) for c in conds_present},
                  incomplete_folders=incomplete)
    write_score_record(out, root, "population", sys.argv[1:], counts, warnings)
    print(f"population: {len(sets)} kept sets, {len(rows)} runs, conditions {conds_present}")
    for c in ki:
        r1 = stats[c]["R1"]
        print(f"  {c}: median ratio {r1['median_rate_ratio']}, median dCV2 {r1['median_d_cv2']}, pause sets "
              f"{r1['n_pause_sets']} (thr {r1['pause_set_threshold']}), R1 {r1['R1_moves_toward_phenotype']}")
    print("  R2 overall:", R2["overall"], "|", {c: v["verdict"] for c, v in r2_age.items()})
    for w in warnings:
        print("warning:", w)


def main():
    ap = argparse.ArgumentParser(prog=SCRIPT_NAME)
    ap.add_argument("run_root")
    ap.add_argument("out_dir")
    ap.add_argument("--mode", choices=("legacy", "screen", "population"), default="legacy")
    ap.add_argument("--window-ms", nargs=2, type=float, metavar=("T0", "TSTOP"))
    ap.add_argument("--n-sets", type=int, default=0, help="screen: expect sets 1..N (missing ones are listed as rejected)")
    ap.add_argument("--kept-list", type=Path, help="population: SCREEN.tsv whose kept sets must equal the sets present")
    a = ap.parse_args()
    default = {"legacy": (T0_MS, TSTOP_MS), "screen": SCREEN_WINDOW_MS, "population": POP_WINDOW_MS}[a.mode]
    set_window(*(a.window_ms or default))
    root, out = Path(a.run_root), Path(a.out_dir)
    if a.mode == "legacy":
        main_legacy(root, out)
    elif a.mode == "screen":
        main_screen(root, out, a.n_sets)
    else:
        main_population(root, out, a.kept_list)


if __name__ == "__main__":
    main()
