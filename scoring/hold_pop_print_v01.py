"""Block-threshold print for a hold-current series run on the population of parameter sets (14 Nov 2026).

Usage: python3 -I -B hold_pop_print_v01.py HOLD_ROOT TASKS_TSV OUT_RUNS_TSV OUT_SETS_TSV [--partial]

Reads the 2-s runs set<NNN>_hold<hold_label>_<condition> under HOLD_ROOT that are listed in TASKS_TSV (columns task,
set_id, condition, hold_nA, hold_label). Writes one row per run to OUT_RUNS_TSV and one row per set and condition to
OUT_SETS_TSV, both by exclusive create, and prints a summary. Nothing else is written. Without --partial the script
stops unless every run of the table is finished (has a SHA256SUMS file); with it, only sets whose runs of a condition
are all finished enter the summary of that condition.

Checks (stop on the first failure): RECORD.json, SPIKES.csv and V_SOMA.npy against the folder's SHA256SUMS; from the
record the factor-file hash of the condition (table FACT below), the injected current of the task, noise seed
40261112 + set and sigma 0.0065 nA, run length 2 s, record step 0.1 ms, a parameter-set file named popset_k<NNN>.json
of that set; 20,000 or 20,001 voltage samples, equal to the record's n_samples when it has one.

Per run (sample i at i x 0.1 ms): spikes in the whole run and in 1.0-2.0 s; firing rate in 1.0-2.0 s; pauses in
1.0-2.0 s (interval > 3 x the median interval of that window); trough = median over non-pause intervals in 1.0-2.0 s
of the lowest voltage between consecutive spikes; time of the last spike; mean voltage in 1.75-2.0 s.
State of a run at its end (the rule of hold_print_v03.py and v04):
  block   no spike in 1.75-2.0 s and mean voltage there above -55 mV
  rest    no spike in 1.75-2.0 s and mean voltage there at or below -55 mV
  pauses  spikes in 1.75-2.0 s and at least one pause in 1.0-2.0 s
  fires   spikes in 1.75-2.0 s and no pause in 1.0-2.0 s
Per set and condition: the states by current; first_block = the lowest current in block; sustained_block = the
lowest current from which every tested current is in block ("none" when there is no such current up to the highest
one tested). Two conditions of a set are compared on each of the two currents: lower, equal or higher; a current of
"none" counts as above every tested current; "neither" when both are "none". Ratios are given for the sets in which
both currents exist.
"""
import csv, hashlib, json, math, os, statistics, sys

import numpy as np

FACT = {
    "WT": "bdbb0f458069e46aefc13f6dfe0dd488de397fcf1546821dca2c3b804157c204",
    "KI_4_protein_Kcnma1_Cacna1g": "367ed13bdaa0f07fa2067a76df2c53e96df998ed4c8e44f42da0b49187555b22",
    "KI_4_protein_Kcnma1_Cacna1g_plus_Kcnc3": "7227127e499155e4665f212896cbb18c886aa26da5bfa690e8de4bd494e92b2e",
    "KI_4_protein_Kcnma1_Cacna1g_plus_Kcnd3": "0fe043257c838d8df8adf3b12081fad2ccd8c85786d0c3828bd6905bda31bf1d",
    "KI_4_protein_Kcnma1_Cacna1g_plus_Kcnc3_Kcnd3": "1940e87b678313438b3fa6db72239c8297698169a3e6229426ae544ca55fe4e6",
}
RUN_S, DT, SIGMA, SEED0, NAN = 2.0, 0.1, 0.0065, 40261112, float("nan")
LETTER = {"fires": "F", "pauses": "P", "block": "B", "rest": "R"}


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


def times(path):
    with open(path) as fh:
        rd = csv.DictReader(fh)
        if not rd.fieldnames or "time_ms" not in rd.fieldnames:
            sys.exit("STOP: no time_ms column in " + path)
        return [float(r["time_ms"]) for r in rd]


def one(folder, task):
    k, cond, hold = int(task["set_id"]), task["condition"], float(task["hold_nA"])
    tab = sums(folder)
    with open(checked(folder, "RECORD.json", tab)) as fh:
        rec = json.load(fh)
    stim, noise = rec.get("stimulus", {}), rec.get("noise", {})
    problems = []
    if rec.get("factor_file_sha256") != FACT[cond]:
        problems.append("factor hash")
    if abs(float(stim.get("hold_nA", NAN)) - hold) > 1e-9:
        problems.append("hold_nA %r" % stim.get("hold_nA"))
    if noise.get("seed") != SEED0 + k or abs(float(noise.get("sigma_nA", -1)) - SIGMA) > 1e-12:
        problems.append("noise seed or sigma")
    if abs(float(rec.get("run_s", -1)) - RUN_S) > 1e-9 or abs(float(rec.get("record_dt_ms", -1)) - DT) > 1e-12:
        problems.append("run length or record step")
    pops = [v for key, v in rec.items() if "popset" in key and isinstance(v, str)]
    if not any(v.endswith("popset_k%03d.json" % k) for v in pops):
        problems.append("parameter-set file %r" % pops)
    if problems:
        sys.exit("STOP: RECORD.json of %s does not match its task: %s" % (folder, "; ".join(problems)))
    t = times(checked(folder, "SPIKES.csv", tab))
    v = np.load(checked(folder, "V_SOMA.npy", tab), allow_pickle=False)
    n_run = int(round(RUN_S * 1000.0 / DT))
    if v.ndim != 1 or len(v) not in (n_run, n_run + 1) or rec.get("n_samples") not in (None, len(v)):
        sys.exit("STOP: V_SOMA.npy in %s has shape %r; expected %d or %d samples, equal to the record's n_samples (%r)"
                 % (folder, v.shape, n_run, n_run + 1, rec.get("n_samples")))
    late = [x for x in t if 1000 <= x <= 2000]
    isi = [b - a for a, b in zip(late, late[1:])]
    pauses, trough = 0, NAN
    if isi:
        m = statistics.median(isi)
        tr = []
        for a, b in zip(late, late[1:]):
            if b - a > 3 * m:
                pauses += 1
                continue
            ia, ib = int(round(a / DT)), int(round(b / DT))
            if ib > ia:
                tr.append(float(v[ia:ib + 1].min()))
        if tr:
            trough = statistics.median(tr)
    mean_end = float(np.mean(v[int(round(1750.0 / DT)):], dtype=np.float64))
    n_end = sum(1 for x in t if x >= 1750)
    state = ("block" if mean_end > -55.0 else "rest") if n_end == 0 else ("pauses" if pauses else "fires")
    return dict(set=k, cond=cond, hold=hold, n=len(t), n_late=len(late), rate=len(late) / 1.0, pauses=pauses,
                last=(max(t) if t else NAN), mean_end=mean_end, trough=trough, state=state)


def q(xs, p):
    xs = sorted(xs)
    if not xs:
        return NAN
    k = (len(xs) - 1) * p
    lo, hi = int(math.floor(k)), int(math.ceil(k))
    return xs[lo] + (xs[hi] - xs[lo]) * (k - lo)


def quart(xs):
    return "n %d; min %.3f, quartiles %.3f / %.3f / %.3f, max %.3f" % (len(xs), min(xs), q(xs, 0.25), q(xs, 0.5), q(xs, 0.75), max(xs)) if xs else "n 0"


def compare(a, b, key):
    """a against b over the sets present in both: counts lower / equal / higher / neither, and the ratios a/b."""
    lo = eq = hi = nei = 0
    ratios = []
    for k in sorted(set(a) & set(b)):
        x, y = a[k][key], b[k][key]
        if x is None and y is None:
            nei += 1
        elif y is None or (x is not None and x < y):
            lo += 1
        elif x is None or x > y:
            hi += 1
        else:
            eq += 1
        if x is not None and y is not None:
            ratios.append(x / y)
    return lo, eq, hi, nei, ratios


def main(root, tasks_tsv, out_runs, out_sets, partial):
    for p in (out_runs, out_sets):
        if os.path.exists(p):
            sys.exit("STOP: output exists: " + p)
    with open(tasks_tsv, "rb") as fh:
        raw = fh.read()
    tasks = list(csv.DictReader(raw.decode("ascii").splitlines(), delimiter="\t"))
    conds, sets = [], []
    for k in tasks:
        if k["condition"] not in FACT:
            sys.exit("STOP: condition without a pinned factor hash: " + k["condition"])
        if k["condition"] not in conds:
            conds.append(k["condition"])
        if int(k["set_id"]) not in sets:
            sets.append(int(k["set_id"]))
    rows, missing = [], []
    for k in tasks:
        folder = os.path.join(root, "set%03d_hold%s_%s" % (int(k["set_id"]), k["hold_label"], k["condition"]))
        if not os.path.isfile(os.path.join(folder, "SHA256SUMS")):
            missing.append(os.path.basename(folder))
            continue
        rows.append(one(folder, k))
    print("HOLD_ROOT %s\nTASKS_TSV %s sha256 %s" % (root, tasks_tsv, hashlib.sha256(raw).hexdigest()))
    print("runs in the table: %d; finished: %d; sets: %d; conditions: %d" % (len(tasks), len(rows), len(sets), len(conds)))
    if missing and not partial:
        sys.exit("STOP: %d runs not finished, first %s (use --partial for a partial print)" % (len(missing), missing[:3]))
    cols = [("set", "set", "%d"), ("condition", "cond", "%s"), ("hold_nA", "hold", "%.3f"), ("spikes_total", "n", "%d"),
            ("spikes_1_2s", "n_late", "%d"), ("rate_1_2s_hz", "rate", "%.1f"), ("pauses_1_2s", "pauses", "%d"),
            ("last_spike_ms", "last", "%.1f"), ("mean_v_1p75_2s", "mean_end", "%.3f"), ("trough_mv", "trough", "%.3f"),
            ("state", "state", "%s")]
    data = ("\n".join(["\t".join(c[0] for c in cols)] + ["\t".join(c[2] % r[c[1]] for c in cols) for r in rows]) + "\n").encode("ascii")
    with open(out_runs, "xb") as fh:
        fh.write(data)
    want = sorted(set(float(k["hold_nA"]) for k in tasks))
    got = {}
    for r in rows:
        got.setdefault((r["cond"], r["set"]), {})[round(r["hold"], 6)] = r
    res = dict((c, {}) for c in conds)
    lines = ["\t".join(["set", "condition", "states", "first_block_nA", "sustained_block_nA", "highest_fires_nA",
                        "rate_at_lowest_current_hz", "finished"])]
    fmt = lambda x: "%.3f" % x if x is not None else "none"
    for k in sets:
        for c in conds:
            have = got.get((c, k), {})
            seq = "".join(LETTER[have[round(h, 6)]["state"]] if round(h, 6) in have else "-" for h in want)
            block = [h for h, s in zip(want, seq) if s == "B"]
            first = min(block) if block else None
            sus = None
            for i, h in enumerate(want):
                if all(s == "B" for s in seq[i:]):
                    sus = h
                    break
            fires = [h for h, s in zip(want, seq) if s == "F"]
            low = have.get(round(want[0], 6))
            complete = len(have) == len(want)
            lines.append("\t".join([str(k), c, seq, fmt(first), fmt(sus), fmt(max(fires) if fires else None),
                                    "%.1f" % low["rate"] if low else "na", "%d/%d" % (len(have), len(want))]))
            if complete:
                res[c][k] = dict(seq=seq, first=first, sus=sus, low=low)
    sdata = ("\n".join(lines) + "\n").encode("ascii")
    with open(out_sets, "xb") as fh:
        fh.write(sdata)
    print("currents (nA): %s" % " ".join("%.3f" % h for h in want))
    print("states by current: F fires, P pauses, B block, R rest")
    for c in conds:
        x = res[c]
        first = [v["first"] for v in x.values() if v["first"] is not None]
        sus = [v["sus"] for v in x.values() if v["sus"] is not None]
        low = [v["seq"][0] for v in x.values()]
        print("%s: sets with every run finished %d" % (c, len(x)))
        print("  first block current (nA): %s; no block up to %.3f nA in %d sets" % (quart(first), want[-1], len(x) - len(first)))
        print("  sustained block current (nA): %s; none in %d sets; first and sustained differ in %d sets"
              % (quart(sus), len(x) - len(sus), sum(1 for v in x.values() if v["first"] != v["sus"])))
        print("  at %.3f nA: fires %d, pauses %d, block %d, rest %d; firing rate there (Hz, sets that fire or pause): %s"
              % (want[0], low.count("F"), low.count("P"), low.count("B"), low.count("R"),
                 quart([v["low"]["rate"] for v in x.values() if v["seq"][0] in "FP"])))
    print("comparisons, set by set (first condition against second; a current of none counts as above every tested current):")
    for i, a in enumerate(conds):
        for b in conds[:i]:
            for key, name in (("first", "first block"), ("sus", "sustained block")):
                lo, eq, hi, nei, ratios = compare(res[a], res[b], key)
                print("  %s against %s, %s: lower %d, equal %d, higher %d, neither %d; ratio where both exist: %s"
                      % (a, b, name, lo, eq, hi, nei, quart(ratios)))
    print("OUT_RUNS_TSV %s: %d rows, sha256 %s" % (out_runs, len(rows), hashlib.sha256(data).hexdigest()))
    print("OUT_SETS_TSV %s: %d rows, sha256 %s" % (out_sets, len(lines) - 1, hashlib.sha256(sdata).hexdigest()))


if __name__ == "__main__":
    a = sys.argv[1:]
    if len(a) == 5 and a[4] == "--partial":
        main(a[0], a[1], a[2], a[3], True)
    elif len(a) == 4:
        main(a[0], a[1], a[2], a[3], False)
    else:
        sys.exit(__doc__)
