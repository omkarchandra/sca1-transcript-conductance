"""Runs with synaptic input (P16, wave 8) against the runs without synaptic input (P15), version by version and
parameter set by parameter set (14 Nov 2026). Read-only on its inputs; prints, and writes one table by exclusive create.

Usage: python3 -I -B p16_vs_p15_print_v01.py P15_DIR P16_DIR OUT_TSV

P15_DIR and P16_DIR hold the per-version tables written by the scoring scripts (<version>.tsv, one row per parameter
set; columns set, rate_ratio, cv_ratio, d_trough_mv, wt_rate_hz, wt_state, ki_state are used). Every version that has
a table in both folders is compared. Set 0 (the published values) exists only with synaptic input; it is printed apart
and is not part of the set-by-set comparison. A set enters a comparison of a measure only if the measure is a number in
both conditions (three or more spikes in the window in the KI and in the WT run of both conditions).

Printed and written per version:
  n            sets compared (numbers in both conditions)
  median       median over those sets, without and with synaptic input
  in range     sets with a firing rate ratio inside the recorded range, 23.78/27.81 to 42.04/31.65 (Bushart et al. 2021,
               Supplementary Table 1), among the n sets, in each condition, and in both
  side         where the median firing rate ratio lies: below, inside or above the recorded range
  rho          Spearman rank correlation across the n sets between the two conditions (ties get the mean rank)
  same sign    sets in which the measure minus its neutral value (1 for ratios, 0 for the difference) has the same sign
"""
import csv, hashlib, math, os, sys
import numpy as np

LO, HI = 23.78 / 27.81, 42.04 / 31.65
NAN = float("nan")


def ranks(x):
    x = np.asarray(x, dtype=float)
    order = np.argsort(x, kind="mergesort")
    r = np.empty(len(x), dtype=float)
    i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and x[order[j + 1]] == x[order[i]]:
            j += 1
        r[order[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return r


def spearman(a, b):
    if len(a) < 3:
        return NAN
    ra, rb = ranks(a), ranks(b)
    sa, sb = ra.std(), rb.std()
    if sa == 0 or sb == 0:
        return NAN
    return float(((ra - ra.mean()) * (rb - rb.mean())).mean() / (sa * sb))


def read(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    out = {}
    for r in rows:
        s = int(r["set"])
        if s in out:
            sys.exit("STOP: set %d twice in %s" % (s, path))
        out[s] = r
    return out


def num(r, key):
    v = float(r[key])
    return v if math.isfinite(v) else NAN


def side(x):
    if not math.isfinite(x):
        return "-"
    return "below" if x < LO else ("above" if x > HI else "inside")


def med(xs):
    return float(np.median(xs)) if len(xs) else NAN


def pair(t15, t16, key, sets):
    a, b, used = [], [], []
    for s in sets:
        x, y = num(t15[s], key), num(t16[s], key)
        if math.isfinite(x) and math.isfinite(y):
            a.append(x); b.append(y); used.append(s)
    return a, b, used


def main(d15, d16, out_tsv):
    if os.path.exists(out_tsv):
        sys.exit("STOP: OUT_TSV exists: " + out_tsv)
    names = sorted(f[:-4] for f in os.listdir(d16) if f.endswith(".tsv") and os.path.isfile(os.path.join(d15, f)))
    only16 = sorted(f[:-4] for f in os.listdir(d16) if f.endswith(".tsv") and not os.path.isfile(os.path.join(d15, f)))
    if not names:
        sys.exit("STOP: no version has a table in both folders")
    print("recorded firing rate ratio range: %.6f to %.6f, both ends included" % (LO, HI))
    print("versions with a table in both folders: %d; with synaptic input only: %s" % (len(names), " ".join(only16) if only16 else "none"))
    cols = ["version", "sets_common", "n_rate", "rate_med_noinput", "rate_med_syn", "side_noinput", "side_syn", "same_side",
            "n_in_range_noinput", "n_in_range_syn", "n_in_range_both", "rho_rate", "n_same_sign_rate",
            "n_cv", "cv_med_noinput", "cv_med_syn", "rho_cv", "n_ahp", "ahp_med_noinput", "ahp_med_syn", "rho_ahp", "n_same_sign_ahp",
            "set0_rate_syn", "set0_cv_syn", "set0_ahp_syn", "sha12_noinput", "sha12_syn"]
    lines = ["\t".join(cols)]
    wt_done = False
    for name in names:
        p15, p16 = os.path.join(d15, name + ".tsv"), os.path.join(d16, name + ".tsv")
        t15, t16 = read(p15), read(p16)
        sets = sorted(s for s in t16 if s in t15 and s != 0)
        if not wt_done:
            a, b, used = pair(t15, t16, "wt_rate_hz", sets)
            print("WT firing rate (Hz) in the %d common sets: median %.2f without and %.2f with synaptic input; rho %.3f"
                  % (len(used), med(a), med(b), spearman(a, b)))
            print("WT sets not in the state 'fires': %d without, %d with synaptic input"
                  % (sum(t15[s]["wt_state"] != "fires" for s in sets), sum(t16[s]["wt_state"] != "fires" for s in sets)))
            wt_done = True
        a, b, used = pair(t15, t16, "rate_ratio", sets)
        in15, in16 = [LO <= x <= HI for x in a], [LO <= x <= HI for x in b]
        c1, c2, _ = pair(t15, t16, "cv_ratio", sets)
        h1, h2, _ = pair(t15, t16, "d_trough_mv", sets)
        s0 = t16.get(0)
        vals = [name, len(sets), len(used), "%.3f" % med(a), "%.3f" % med(b), side(med(a)), side(med(b)),
                "yes" if side(med(a)) == side(med(b)) else "NO", sum(in15), sum(in16), sum(x and y for x, y in zip(in15, in16)),
                "%.3f" % spearman(a, b), sum((x > 1) == (y > 1) and (x < 1) == (y < 1) for x, y in zip(a, b)),
                len(c1), "%.3f" % med(c1), "%.3f" % med(c2), "%.3f" % spearman(c1, c2),
                len(h1), "%.3f" % med(h1), "%.3f" % med(h2), "%.3f" % spearman(h1, h2),
                sum((x > 0) == (y > 0) and (x < 0) == (y < 0) for x, y in zip(h1, h2)),
                ("%.3f" % num(s0, "rate_ratio")) if s0 else "-", ("%.3f" % num(s0, "cv_ratio")) if s0 else "-",
                ("%.3f" % num(s0, "d_trough_mv")) if s0 else "-",
                hashlib.sha256(open(p15, "rb").read()).hexdigest()[:12], hashlib.sha256(open(p16, "rb").read()).hexdigest()[:12]]
        lines.append("\t".join(str(v) for v in vals))
    data = ("\n".join(lines) + "\n").encode("ascii")
    with open(out_tsv, "xb") as fh:
        fh.write(data)
    sys.stdout.write(data.decode("ascii"))
    print("OUT_TSV %s: %d versions, sha256 %s" % (out_tsv, len(names), hashlib.sha256(data).hexdigest()))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
