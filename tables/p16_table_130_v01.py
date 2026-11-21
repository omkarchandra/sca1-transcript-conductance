"""P16 (synaptic-input runs, wave 8): the one-line-per-version table over the 130 kept parameter sets, without set 0
(set 0 = the published values), with set 0 given apart (15 Nov 2026, v01).

Usage: python3 -I -B p16_table_130_v01.py P16_DIR P16_TABLE OUT_TSV

P16_DIR    the folder with the per-version tables <version>.tsv written by p16_variant_score_v02.py
P16_TABLE  the table written by p16_score_all_v01.py (P16_TABLE_v01_0900.tsv; its sha256 is checked below)
OUT_TSV    the new table; must not exist; written by exclusive create

Reads only; changes none of the existing files. The columns are those of P16_TABLE, recomputed from each version's
per-set table with the definitions of p16_variant_score_v02.py (its printed summary, from which p16_score_all_v01.py
takes the table entries) and of p16_score_all_v01.py (recorded_counts). The quartile function q and the recorded-range
bounds REC_LO and REC_HI are taken from p16_variant_score_v02.py itself: the file is loaded by path after its sha256 is
checked; its main part does not run (it runs only under __main__).

Self-test, run first: over all 131 rows, every cell of the 21 scored rows of P16_TABLE except "status" must be
reproduced exactly as printed, and the first 12 hex characters of each per-set table's sha256 must equal its
tsv_sha256_12. On the first difference the script stops and writes nothing.
Then the same recomputation over the 130 rows with set != 0. OUT_TSV holds the 21 scored versions, with the header of
P16_TABLE (status = "130 sets, set 0 apart"; sets = 130; tsv_sha256_12 = that of the per-set table read) and five more
columns for set 0 alone, copied as printed in the per-set table: rate ratio, CV ratio, trough difference, KI state,
WT state.
"""
import csv, hashlib, importlib.util, os, sys

SCORER = "${LOCAL_ROOT}/docs/20261026_cc/p16_score_20261114/p16_variant_score_v02.py"
SCORER_SHA = "1e0a200d7d6a6e066e6d48796b84fec47f527a5b37a4185274a950fd184d9a5e"
TABLE_SHA = "6620377a263cfb016845e6d1b6fff21101aa79543da608f4a1d5114eafa6c660"
SCORED = ("scored now", "scored earlier")
SET0_COLS = ["set0_rate_ratio", "set0_cv_ratio", "set0_d_trough_mv", "set0_ki_state", "set0_wt_state"]


def sha_of(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def load_scorer():
    if sha_of(SCORER) != SCORER_SHA:
        sys.exit("STOP: p16_variant_score_v02.py does not have the expected sha256")
    spec = importlib.util.spec_from_file_location("p16_variant_score_v02", SCORER)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_tsv(path):
    with open(path, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        head = next(rd)
        return head, [r for r in rd]


def entries(rows, idx, m):
    """The table entries from the given per-set rows, as p16_variant_score_v02.py prints them and p16_score_all_v01.py
    takes them: sets, both_fire, rate quartiles, n 0.85-1.33, CV medians, CV counts, trough median and count, n of all
    three, KI states, and the two recorded-range counts."""
    f = lambda r, c: float(r[idx[c]])
    fr = [r for r in rows if f(r, "rate_ratio") == f(r, "rate_ratio")]
    nop = [r for r in fr if int(r[idx["ki_pauses"]]) == 0 and int(r[idx["wt_pauses"]]) == 0]
    rr = [f(r, "rate_ratio") for r in fr]
    cv = [f(r, "cv_ratio") for r in fr if f(r, "cv_ratio") == f(r, "cv_ratio")]
    cvnop = [x for x in (f(r, "cv_ratio") for r in nop) if x == x]
    dt = [f(r, "d_trough_mv") for r in fr if f(r, "d_trough_mv") == f(r, "d_trough_mv")]
    band = [r for r in fr if 0.85 <= f(r, "rate_ratio") <= 1.33]
    st = [r[idx["ki_state"]] for r in rows]
    n_rec = n_rec_all = 0
    for r in rows:
        x, c, d = f(r, "rate_ratio"), f(r, "cv_ratio"), f(r, "d_trough_mv")
        if m.REC_LO <= x <= m.REC_HI:
            n_rec += 1
            n_rec_all += c >= 1.5 and d > 0
    return {
        "sets": "%d" % len(rows), "both_fire": "%d" % len(fr),
        "rate_q25": "%.3f" % m.q(rr, 0.25), "rate_median": "%.3f" % m.q(rr, 0.5), "rate_q75": "%.3f" % m.q(rr, 0.75),
        "n_rate_0.85-1.33": "%d" % len(band),
        "cv_ratio_median": "%.3f" % m.q(cv, 0.5), "cv_ratio_median_no_pause": "%.3f" % m.q(cvnop, 0.5),
        "n_cv_ge_1.5": "%d" % sum(x >= 1.5 for x in cv), "n_cv_ge_1.8": "%d" % sum(x >= 1.8 for x in cv),
        "d_trough_median_mv": "%.3f" % m.q(dt, 0.5), "n_trough_above": "%d" % sum(x > 0 for x in dt),
        "n_all_three": "%d" % sum(f(r, "cv_ratio") >= 1.5 and f(r, "d_trough_mv") > 0 for r in band),
        "ki_fires": "%d" % st.count("fires"), "ki_pauses": "%d" % st.count("pauses"),
        "ki_block": "%d" % st.count("block"), "ki_rest": "%d" % st.count("rest"),
        "n_rate_recorded": "%d" % n_rec, "n_all_three_recorded": "%d" % n_rec_all,
    }


def main(p16_dir, table, out_tsv):
    if os.path.exists(out_tsv):
        sys.exit("STOP: OUT_TSV exists: " + out_tsv)
    if sha_of(table) != TABLE_SHA:
        sys.exit("STOP: P16_TABLE does not have the expected sha256")
    m = load_scorer()
    head, trows = read_tsv(table)
    ti = {c: i for i, c in enumerate(head)}
    scored = [r for r in trows if len(r) == len(head) and r[ti["status"]] in SCORED]
    if len(scored) != 21:
        sys.exit("STOP: %d scored rows in P16_TABLE, expected 21" % len(scored))
    data, cells = [], 0
    for tr in scored:
        v = tr[ti["variant"]]
        path = os.path.join(p16_dir, v + ".tsv")
        tsha = sha_of(path)
        if tsha[:12] != tr[ti["tsv_sha256_12"]]:
            sys.exit("SELF-TEST FAIL: %s: per-set table sha256 %s, table says %s" % (v, tsha[:12], tr[ti["tsv_sha256_12"]]))
        phead, prows = read_tsv(path)
        pi = {c: i for i, c in enumerate(phead)}
        if len(prows) != 131 or sorted(int(r[pi["set"]]) for r in prows).count(0) != 1:
            sys.exit("SELF-TEST FAIL: %s: %d rows in the per-set table, or set 0 not there once" % (v, len(prows)))
        e = entries(prows, pi, m)
        for c in head:
            if c in ("variant", "status", "tsv_sha256_12"):
                continue
            got = e[c] if c in e else tr[ti[c]]
            if c in e:
                cells += 1
            if got != tr[ti[c]]:
                sys.exit("SELF-TEST FAIL: first differing cell: version %s, column %s: table %s, recomputed %s"
                         % (v, c, tr[ti[c]], got))
        data.append((tr, tsha, prows, pi))
    print("SELF-TEST PASS: 21 versions, %d recomputed cells over 131 sets equal P16_TABLE (sha256 %s) as printed; "
          "21 per-set table hashes equal tsv_sha256_12" % (cells, TABLE_SHA))
    lines = ["\t".join(head + SET0_COLS)]
    for tr, tsha, prows, pi in data:
        kept = [r for r in prows if int(r[pi["set"]]) != 0]
        s0 = [r for r in prows if int(r[pi["set"]]) == 0][0]
        e = entries(kept, pi, m)
        out = []
        for c in head:
            if c == "status":
                out.append("130 sets, set 0 apart")
            elif c == "tsv_sha256_12":
                out.append(tsha[:12])
            else:
                out.append(e[c] if c in e else tr[ti[c]])
        out += [s0[pi["rate_ratio"]], s0[pi["cv_ratio"]], s0[pi["d_trough_mv"]], s0[pi["ki_state"]], s0[pi["wt_state"]]]
        lines.append("\t".join(out))
    blob = ("\n".join(lines) + "\n").encode("ascii")
    with open(out_tsv, "xb") as fh:
        fh.write(blob)
    print("OUT_TSV %s: %d versions, sha256 %s" % (out_tsv, len(data), hashlib.sha256(blob).hexdigest()))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])

#%% debug
