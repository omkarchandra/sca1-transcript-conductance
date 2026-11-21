"""P15 variant scoring: one variant root against the WT root, parameter set by parameter set (14 Nov 2026).

Usage: python3 -I -B p15_variant_score_v01.py WT_ROOT KI_ROOT FACTOR_SHA256 OUT_TSV [--partial]

WT_ROOT holds folders set<NNN>_WT, KI_ROOT holds folders set<NNN>_KI_18w_mid (the folder name is the same in every
variant root; which variant a run belongs to is read from its RECORD.json). FACTOR_SHA256 is the sha256 of the factor
file the variant was run with. OUT_TSV must not exist; it is written by exclusive create with one row per set.
Nothing else is written; nothing is changed or removed. The script stops unless KI_ROOT holds exactly 130 sets
(the kept population); --partial lifts that one check, for tests.

Checks (the script stops on the first failure): RECORD.json, SPIKES.csv, DSPIKES.csv and V_SOMA.npy of every run
against the folder's SHA256SUMS; every KI record carries FACTOR_SHA256 and one and the same condition; run length
10 s, record step 0.1 ms, no injected current; KI and WT run of a set have the same noise seed and sigma;
V_SOMA.npy has 100,001 samples.

Measures, window 5-10 s (both ends included), as in the earlier P15 prints:
  rate        spikes in the window / 5 s
  CV          population SD of the window interspike intervals / their mean
  CV2         mean of 2|b - a| / (a + b) over consecutive intervals a, b
  pause       an interval longer than 3 x the run's median interval; share = time in pauses / 5 s
  between     number of non-pause intervals / their summed duration (rate between pauses)
  trough      median, over the non-pause intervals, of the lowest somatic voltage between the two spikes (mV)
  mean_v      mean somatic voltage of the window (mV)
  dspikes     dendritic spikes in the window (DSPIKES.csv; its time_ms column, or its first column if there is none)
State of a run: fires (3 or more window spikes, no pause), pauses (3 or more spikes, at least one pause),
block (fewer than 3 spikes and mean_v above -55 mV), rest (fewer than 3 spikes and mean_v at or below -55 mV).
Ratios and differences are KI relative to the WT run of the same set and are given when both runs have 3 or more spikes.
"""
import csv, hashlib, json, math, os, re, statistics, sys

import numpy as np

T0, T1, DT, N_SAMPLES, RUN_S = 5000.0, 10000.0, 0.1, 100001, 10.0
WIN_S = (T1 - T0) / 1000.0
NAN = float("nan")


def sums(folder):
    out = {}
    with open(os.path.join(folder, "SHA256SUMS")) as fh:
        for line in fh:
            h, n = line.rstrip("\n").split("  ", 1)
            out[n] = h
    return out


def checked(folder, name, table):
    path = os.path.join(folder, name)
    with open(path, "rb") as fh:
        got = hashlib.sha256(fh.read()).hexdigest()
    if table.get(name) != got:
        sys.exit("STOP: hash check failed for %s in %s" % (name, folder))
    return path


def times(path, strict):
    with open(path) as fh:
        rd = csv.DictReader(fh)
        names = rd.fieldnames or []
        if "time_ms" in names:
            key = "time_ms"
        elif strict or not names:
            if strict:
                sys.exit("STOP: no time_ms column in " + path)
            return []
        else:
            key = names[0]
        return [float(r[key]) for r in rd]


def one(folder):
    tab = sums(folder)
    with open(checked(folder, "RECORD.json", tab)) as fh:
        rec = json.load(fh)
    stim, noise = rec.get("stimulus", {}), rec.get("noise", {})
    if abs(float(rec.get("run_s", -1)) - RUN_S) > 1e-9 or abs(float(rec.get("record_dt_ms", -1)) - DT) > 1e-12:
        sys.exit("STOP: run length or record step not 10 s / 0.1 ms in " + folder)
    if abs(float(stim.get("hold_nA", NAN))) > 1e-12:
        sys.exit("STOP: injected current not zero in " + folder)
    t = times(checked(folder, "SPIKES.csv", tab), True)
    d = times(checked(folder, "DSPIKES.csv", tab), False)
    v = np.load(checked(folder, "V_SOMA.npy", tab), allow_pickle=False)
    if v.ndim != 1 or len(v) != N_SAMPLES:
        sys.exit("STOP: V_SOMA.npy in %s has shape %r, expected (%d,)" % (folder, v.shape, N_SAMPLES))
    i0, i1 = int(round(T0 / DT)), int(round(T1 / DT))
    w = [x for x in t if T0 <= x <= T1]
    r = dict(cond=rec.get("condition"), fsha=rec.get("factor_file_sha256"), seed=noise.get("seed"),
             sigma=noise.get("sigma_nA"), n=len(w), rate=len(w) / WIN_S, cv=NAN, cv2=NAN, pauses=0, share=0.0,
             between=NAN, trough=NAN, mean_v=float(np.mean(v[i0:i1 + 1], dtype=np.float64)),
             dsp=sum(1 for x in d if T0 <= x <= T1))
    isi = [b - a for a, b in zip(w, w[1:])]
    if len(isi) >= 2:
        m = statistics.median(isi)
        r["cv"] = statistics.pstdev(isi) / (sum(isi) / len(isi))
        r["cv2"] = sum(2.0 * abs(b - a) / (a + b) for a, b in zip(isi, isi[1:])) / (len(isi) - 1)
        keep = [x for x in isi if x <= 3 * m]
        r["pauses"] = len(isi) - len(keep)
        r["share"] = sum(x for x in isi if x > 3 * m) / (T1 - T0)
        r["between"] = len(keep) / (sum(keep) / 1000.0)
        tr = []
        for a, b in zip(w, w[1:]):
            if b - a > 3 * m:
                continue
            ia, ib = int(round(a / DT)), int(round(b / DT))
            if ib > ia:
                tr.append(float(v[ia:ib + 1].min()))
        if tr:
            r["trough"] = statistics.median(tr)
        r["state"] = "pauses" if r["pauses"] else "fires"
    else:
        r["state"] = "block" if r["mean_v"] > -55.0 else "rest"
    return r


def q(xs, p):
    xs = sorted(xs)
    if not xs:
        return NAN
    k = (len(xs) - 1) * p
    lo, hi = int(math.floor(k)), int(math.ceil(k))
    return xs[lo] + (xs[hi] - xs[lo]) * (k - lo)


def five(label, xs, fmt="%.3f"):
    xs = [x for x in xs if x == x]
    f = "%s: n %d; min " + fmt + ", quartiles " + fmt + " / " + fmt + " / " + fmt + ", max " + fmt
    print(f % (label, len(xs), min(xs) if xs else NAN, q(xs, 0.25), q(xs, 0.5), q(xs, 0.75), max(xs) if xs else NAN))


def main(wt_root, ki_root, fsha, out_tsv, partial=False):
    if not re.fullmatch(r"[0-9a-f]{64}", fsha):
        sys.exit("STOP: FACTOR_SHA256 is not 64 hex characters")
    if os.path.exists(out_tsv):
        sys.exit("STOP: OUT_TSV exists: " + out_tsv)
    ids = sorted(int(m.group(1)) for m in (re.fullmatch(r"set(\d+)_KI_18w_mid", n) for n in os.listdir(ki_root)) if m)
    if not ids:
        sys.exit("STOP: no set<NNN>_KI_18w_mid folders in " + ki_root)
    if len(ids) != 130 and not partial:
        sys.exit("STOP: %d sets in KI_ROOT, expected 130" % len(ids))
    rows, conds = [], {}
    for s in ids:
        wf, kf = os.path.join(wt_root, "set%03d_WT" % s), os.path.join(ki_root, "set%03d_KI_18w_mid" % s)
        if not os.path.isdir(wf):
            sys.exit("STOP: WT folder missing for set %d: %s" % (s, wf))
        if not os.path.isfile(os.path.join(kf, "SHA256SUMS")):
            sys.exit("STOP: run not finished (no SHA256SUMS) in " + kf)
        w, k = one(wf), one(kf)
        if k["fsha"] != fsha:
            sys.exit("STOP: set %d was run with factor file %s, not %s" % (s, k["fsha"], fsha))
        if k["seed"] != w["seed"] or k["sigma"] != w["sigma"]:
            sys.exit("STOP: set %d: noise seed or sigma differs between KI and WT" % s)
        conds[k["cond"]] = conds.get(k["cond"], 0) + 1
        both = w["n"] >= 3 and k["n"] >= 3
        rows.append(dict(set=s, seed=k["seed"], w=w, k=k,
                         rr=(k["rate"] / w["rate"]) if both else NAN, cvr=(k["cv"] / w["cv"]) if both and w["cv"] > 0 else NAN,
                         dcv2=(k["cv2"] - w["cv2"]) if both else NAN, dtr=(k["trough"] - w["trough"]) if both else NAN,
                         bwr=(k["between"] / w["between"]) if both else NAN, dmv=k["mean_v"] - w["mean_v"]))
    if len(conds) != 1:
        sys.exit("STOP: more than one condition in KI_ROOT: %r" % conds)
    cols = ["set", "noise_seed", "wt_spikes", "ki_spikes", "wt_rate_hz", "ki_rate_hz", "rate_ratio", "wt_cv", "ki_cv", "cv_ratio",
            "wt_cv2", "ki_cv2", "d_cv2", "wt_pauses", "ki_pauses", "ki_pause_share", "wt_between_hz", "ki_between_hz",
            "between_ratio", "wt_trough_mv", "ki_trough_mv", "d_trough_mv", "wt_mean_v", "ki_mean_v", "d_mean_v",
            "wt_dspikes", "ki_dspikes", "wt_state", "ki_state"]
    lines = ["\t".join(cols)]
    for r in rows:
        w, k = r["w"], r["k"]
        vals = [r["set"], r["seed"], w["n"], k["n"], "%.2f" % w["rate"], "%.2f" % k["rate"], "%.6f" % r["rr"], "%.6f" % w["cv"],
                "%.6f" % k["cv"], "%.6f" % r["cvr"], "%.6f" % w["cv2"], "%.6f" % k["cv2"], "%.6f" % r["dcv2"], w["pauses"],
                k["pauses"], "%.6f" % k["share"], "%.3f" % w["between"], "%.3f" % k["between"], "%.6f" % r["bwr"],
                "%.3f" % w["trough"], "%.3f" % k["trough"], "%.3f" % r["dtr"], "%.3f" % w["mean_v"], "%.3f" % k["mean_v"],
                "%.3f" % r["dmv"], w["dsp"], k["dsp"], w["state"], k["state"]]
        lines.append("\t".join(str(x) for x in vals))
    data = ("\n".join(lines) + "\n").encode("ascii")
    with open(out_tsv, "xb") as fh:
        fh.write(data)

    def ids_of(test):
        return [r["set"] for r in rows if test(r)]

    print("WT_ROOT %s\nKI_ROOT %s" % (wt_root, ki_root))
    print("condition %s; factor file sha256 %s; sets %d; noise sigma %s nA; every file hash-checked"
          % (list(conds)[0], fsha, len(rows), sorted(set(str(r["k"]["sigma"]) for r in rows))))
    print("noise seed equals 40261112 + set number in %d of %d sets" % (sum(r["seed"] == 40261112 + r["set"] for r in rows), len(rows)))
    for tag in ("w", "k"):
        st = [r[tag]["state"] for r in rows]
        print("%s states: fires %d, pauses %d, block %d, rest %d"
              % ("WT" if tag == "w" else "KI", st.count("fires"), st.count("pauses"), st.count("block"), st.count("rest")))
    fr = [r for r in rows if r["rr"] == r["rr"]]
    nop = [r for r in fr if r["k"]["pauses"] == 0 and r["w"]["pauses"] == 0]
    print("sets with 3 or more spikes in both runs: %d; of those without a pause in either run: %d" % (len(fr), len(nop)))
    rr = [r["rr"] for r in fr]
    five("firing rate ratio KI/WT", rr)
    print("  ratio below 0.85: %d; 0.85 to below 0.90: %d; 0.90-1.10: %d; above 1.10 to 1.33: %d; above 1.33 to 1.50: %d; above 1.50: %d"
          % (sum(x < 0.85 for x in rr), sum(0.85 <= x < 0.90 for x in rr), sum(0.90 <= x <= 1.10 for x in rr),
             sum(1.10 < x <= 1.33 for x in rr), sum(1.33 < x <= 1.50 for x in rr), sum(x > 1.50 for x in rr)))
    five("rate between pauses, KI/WT", [r["bwr"] for r in fr])
    five("WT CV", [r["w"]["cv"] for r in fr], "%.4f")
    five("KI CV", [r["k"]["cv"] for r in fr], "%.4f")
    cv = [r["cvr"] for r in fr if r["cvr"] == r["cvr"]]
    five("CV ratio KI/WT", cv)
    print("  CV ratio below 1.0: %d; 1.0 to below 1.5: %d; 1.5 to below 1.8: %d; 1.8 to below 2.8: %d; 2.8 or more: %d"
          % (sum(x < 1.0 for x in cv), sum(1.0 <= x < 1.5 for x in cv), sum(1.5 <= x < 1.8 for x in cv),
             sum(1.8 <= x < 2.8 for x in cv), sum(x >= 2.8 for x in cv)))
    five("CV ratio KI/WT, sets without pauses", [r["cvr"] for r in nop])
    five("CV2 difference KI - WT", [r["dcv2"] for r in fr], "%.4f")
    dt = [r["dtr"] for r in fr if r["dtr"] == r["dtr"]]
    five("trough difference KI - WT (mV)", dt)
    print("  KI trough above WT: %d; equal: %d; below: %d; median WT trough %.3f, median KI trough %.3f mV"
          % (sum(x > 0 for x in dt), sum(x == 0 for x in dt), sum(x < 0 for x in dt),
             q([r["w"]["trough"] for r in fr if r["dtr"] == r["dtr"]], 0.5), q([r["k"]["trough"] for r in fr if r["dtr"] == r["dtr"]], 0.5)))
    five("trough difference, sets without pauses (mV)", [r["dtr"] for r in nop])
    five("mean window voltage difference KI - WT (mV), all sets", [r["dmv"] for r in rows])
    both_ok = [r for r in fr if 0.85 <= r["rr"] <= 1.33]
    print("sets with rate ratio 0.85-1.33: %d; of those with CV ratio 1.5 or more: %d; with trough difference above 0: %d; with both: %d"
          % (len(both_ok), sum(r["cvr"] >= 1.5 for r in both_ok), sum(r["dtr"] > 0 for r in both_ok),
             sum(r["cvr"] >= 1.5 and r["dtr"] > 0 for r in both_ok)))
    print("KI sets with pauses: %s" % " ".join("%d(%d;%.1f%%;ratio %.3f)" % (r["set"], r["k"]["pauses"], 100 * r["k"]["share"], r["rr"])
                                              for r in rows if r["k"]["state"] == "pauses"))
    print("KI sets in block: %s" % " ".join("%d(%.1f mV)" % (r["set"], r["k"]["mean_v"]) for r in rows if r["k"]["state"] == "block"))
    print("KI sets at rest: %s" % " ".join("%d(%.1f mV)" % (r["set"], r["k"]["mean_v"]) for r in rows if r["k"]["state"] == "rest"))
    print("WT sets not firing without pauses: %s" % ids_of(lambda r: r["w"]["state"] != "fires"))
    print("sets with dendritic spikes in the window: WT %d, KI %d" % (sum(r["w"]["dsp"] > 0 for r in rows), sum(r["k"]["dsp"] > 0 for r in rows)))
    print("OUT_TSV %s: %d rows, sha256 %s" % (out_tsv, len(rows), hashlib.sha256(data).hexdigest()))


if __name__ == "__main__":
    if len(sys.argv) == 6 and sys.argv[5] == "--partial":
        main(*sys.argv[1:5], partial=True)
    elif len(sys.argv) == 5:
        main(*sys.argv[1:])
    else:
        sys.exit(__doc__)
