"""Hold-current series print for any task table (14 Nov 2026).

Usage: python3 -I -B hold_print_v04.py HOLD_ROOT TASKS_TSV OUT_TSV [--partial]

Reads the 2-s runs hold<hold_label>_<condition> under HOLD_ROOT that are listed in TASKS_TSV (columns task, sigma_nA,
sigma_label, condition, noise_seed, hold_nA, hold_label). Writes one row per run to OUT_TSV by exclusive create and
prints a summary per condition. Nothing else is written. Without --partial the script stops unless every run of the
table is finished (has a SHA256SUMS file).

Checks (stop on the first failure): RECORD.json, SPIKES.csv, DSPIKES.csv and V_SOMA.npy against the folder's
SHA256SUMS; from the record the factor-file hash of the condition (table FACT below), the injected current, the noise
seed and sigma of the task table, run length 2 s, record step 0.1 ms, no popset; 20,000 or 20,001 voltage samples,
equal to the record's n_samples when it has one.

Per run, from SPIKES.csv and V_SOMA.npy (sample i at i x 0.1 ms):
  spikes in the whole run and in 1.0-2.0 s; firing rate in 1.0-2.0 s; time of the last spike;
  pauses in 1.0-2.0 s (interval > 3 x the median interval of that window);
  trough = median over non-pause intervals in 1.0-2.0 s of the lowest voltage between consecutive spikes;
  mean voltage in 1.75-2.0 s and in 1.5-2.0 s; lowest and highest voltage in 1.5-2.0 s;
  dendritic spikes in 1.0-2.0 s (DSPIKES.csv; -1 if that file has no readable time column).
State of a run at its end:
  block   no spike in 1.75-2.0 s and mean voltage there above -55 mV
  rest    no spike in 1.75-2.0 s and mean voltage there at or below -55 mV
  pauses  spikes in 1.75-2.0 s and at least one pause in 1.0-2.0 s
  fires   spikes in 1.75-2.0 s and no pause in 1.0-2.0 s
The column state_v02 gives the state by the rule of hold_print_v02.py (window 1.5-2.0 s), for comparison.
Per condition: the lowest current in block, whether every higher current is in block, that current relative to WT,
and the trough difference from WT at the same firing rate (WT trough interpolated against WT rate over the WT runs
in state fires; median over the condition's runs in state fires whose rate lies inside the WT range).
"""
import csv, hashlib, json, os, statistics, sys

import numpy as np

FACT = {
    "WT": "bdbb0f458069e46aefc13f6dfe0dd488de397fcf1546821dca2c3b804157c204",
    "KI_18w_mid": "660f8462aec5b6671a279c5041ee539909b8a8395de3c7eed30994ca283c937a",
    "KI_1c_restore_Kcnc3": "df44d4080056959e957a431436c19fc5f7978538259d6fbdc4a987479ca6d7f0",
    "KI_1c_restore_Kcnma1": "f1d8bcb702484db8cfd3087407c467a91e284e3e10c42601548c4513b4e733e7",
    "KI_1c_restore_Kcnd3": "d54bd6e4d307beecd363e16e065c83151702d64dec3e91ceae7fd08f92c881bf",
    "KI_1c_restore_Cacna1a": "5b11a8fcbf75fb4999bf739a22096b5216f50f240b78efc59af8b4fde6b51fbe",
    "KI_1c_restore_Cacna1g": "1b6b5bd5f3ff27a8ca557033b3568e488e649ac63a18f5213eaec5d75eb9fc1f",
    "KI_2b_only_Kcnma1_Cacna1g": "aa182bf2b36f35b758ae9235baa62ef89f560506db69c16d35de725c8ac322ee",
    "KI_3_only_Kcnc3": "120bd84515cc864a8efa0c3dd0d97ba37646557bbf1dc5d83cd01166605a9982",
    "KI_3_only_Kcnma1": "10610f9327d94f7b5563f743aee0980c6a8b9e13741fdc836f7ffff44f27a4e9",
    "KI_3_only_Kcnd3": "04b388690f34cd9ad2520e9ddcddc44c853cb1e4bc1a7984561e49d126ec9581",
    "KI_3_only_Cacna1a": "541a8b940179df4c80cf1ed17d5bc25e3abc74ff4e10a110df46055d8c9c83ea",
    "KI_3_only_Cacna1g": "65f5be13402bb0077714805dc8c90695ef80d4fb8b0a15bb8dd9cb9a526ccb20",
    "KI_3_only_Kcnc3_Kcnma1": "b9f105ce7cf168c75ea197c77550e93e85dc33d8635095c5aee6ced6d6834cee",
    "KI_3_only_Kcnc3_Kcnd3": "34083ac4e53a384b54a7629fb89c5a9af5effbde840a30770b52f702f439aa6e",
    "KI_3_only_Kcnc3_Cacna1a": "af845509b03ace32b3476e4ceeef7b650fc615a149712c2d55faae44838a9f23",
    "KI_3_only_Kcnc3_Cacna1g": "f058363e93a0bd0ccda509746210fa3d699604197807a9e5d9567ad1626763e0",
    "KI_3_only_Kcnma1_Kcnd3": "993b670e1fc5c4b145f83744ec40d35a89804aa65bab767fced5d4aab4560648",
    "KI_3_only_Kcnma1_Cacna1a": "13874e5eab4530cbcbbb1d8e1aae9fac696c5eb1720a6e6b4115a63e9ebdf5f1",
    "KI_3_only_Kcnd3_Cacna1a": "9563472a18a9a394e2b76e1c813c8bd9da72198746d5f118c75bf2a37de34ae4",
    "KI_3_only_Kcnd3_Cacna1g": "70130bc2a12bbf1fb43b5fb13a36394e80b1a65690a9612c00b99bbab4f5c2bf",
    "KI_3_only_Cacna1a_Cacna1g": "fef91b4a6ffbd7d824600721ac0573364c107ff3429848334265af1d10db5708",
    "KI_3_only_Kcnc3_Kcnma1_Kcnd3": "caff5fb8d0a17ecfdafee8c4e79c3c141431546b548eda1cd71830447f9ae3ba",
    "KI_3_only_Kcnc3_Kcnma1_Cacna1a": "0339b5e53d54557b1da973b48a2898f09e93501815a56d1fdccf6e8bb1dddb04",
    "KI_3_only_Kcnc3_Kcnma1_Cacna1g": "89453f3af78955f8c2d9aa48761748d82098b5dd75596c1c846bc7e6c9b0cdaf",
    "KI_3_only_Kcnc3_Kcnd3_Cacna1a": "ff3cef75f7af83eb9032224d13301e7b0afae539f50992def38d1eee2c3cb6c5",
    "KI_3_only_Kcnc3_Kcnd3_Cacna1g": "676f6ba8381a37a8263bc7d12f084e7d42059a73e986a698ac6a760b9de344cc",
    "KI_3_only_Kcnc3_Cacna1a_Cacna1g": "bfe8ff885528b299f78a96d3ba9564e408573c8214943618bff3ed35971bafdc",
    "KI_3_only_Kcnma1_Kcnd3_Cacna1a": "735afd2bd39f9e35668af2bfe16b21f65d8fb1ff66fbca3eec0b485ac7e7e503",
    "KI_3_only_Kcnma1_Kcnd3_Cacna1g": "35885d4c2124156a830676a005e0230070a4cc9cf124219e90a863661a3224ad",
    "KI_3_only_Kcnma1_Cacna1a_Cacna1g": "99eeed1006ef972048e87268c5955450401b4ebe109ab1911a510e780ec33a0f",
    "KI_3_only_Kcnd3_Cacna1a_Cacna1g": "6702c6624ce7b3a9d60e28a42c29760ce022db6effb9988eb88e9a994c234acd",
    "KI_4_protein_Kcnma1_Cacna1g": "367ed13bdaa0f07fa2067a76df2c53e96df998ed4c8e44f42da0b49187555b22",
    "KI_4_protein_Kcnma1_Cacna1g_plus_Kcnc3": "7227127e499155e4665f212896cbb18c886aa26da5bfa690e8de4bd494e92b2e",
    "KI_4_protein_Kcnma1_Cacna1g_plus_Kcnd3": "0fe043257c838d8df8adf3b12081fad2ccd8c85786d0c3828bd6905bda31bf1d",
    "KI_4_protein_Kcnma1_Cacna1g_plus_Kcnc3_Kcnd3": "1940e87b678313438b3fa6db72239c8297698169a3e6229426ae544ca55fe4e6",
    "KI_5_protein_Kcnma1_Cacna1g_Cacna1a_0p95": "106e921da57e838ee8dbeda719dae66c832e0d1eba256ac63314a6756a66da9f",
    "KI_5_protein_Kcnma1_Cacna1g_Cacna1a_0p90": "34489efdf2075a044b8912ce7911353f4b8a118d0b688dc2cbbf73af344f7fb1",
    "KI_5_protein_Kcnma1_Cacna1g_Cacna1a_0p85": "a657ad914d32cfb3d1705def32f16eb3abe817750200958f1f3a24a8db15268f",
    "KI_5_protein_Kcnma1_Cacna1g_Cacna1a_0p80": "7d5b1070f759ce90d9660c0c312ce62a5cc0be97beafd960c227025178a46d92",
    "KI_6_aux_beta4": "20acc7edcb87392da16fd3d254598751c2627ea311ee8713e223eae90adeaaac",
    "KI_6_aux_a2d2": "fe4245dafe16f1e1accf50e37e854c9f6bccc395508e77e1e90b95aabc8c2a4e",
    "KI_6_aux_beta4_a2d2": "78f2c7ceda57422f8e22eb0fcbe69ee85544cee30ed64eaac1b330b57cc50c0c",
    "KI_6_protein_Kcnma1_Cacna1g_aux_beta4": "caf7c99bb694273463e642f0b2b6913ee9aebac96166a75113716bc9a2f27c2e",
    "KI_6_protein_Kcnma1_Cacna1g_aux_beta4_a2d2": "43309ae41665bb8cffd48f990099ed42fd29e2f8cb9bed297df2054f6e01acf1",
    "KI_6_restore_Cacna1a_aux_beta4": "956f7830e72fb517454f1921b11bb3194ae0095843b68e9dea2822218101697b",
}
RUN_S, DT, NAN = 2.0, 0.1, float("nan")
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
        key = "time_ms" if rd.fieldnames and "time_ms" in rd.fieldnames else (rd.fieldnames[0] if rd.fieldnames else None)
        return [float(r[key]) for r in rd] if key else []


def one(folder, task):
    cond, hold = task["condition"], float(task["hold_nA"])
    tab = sums(folder)
    with open(checked(folder, "RECORD.json", tab)) as fh:
        rec = json.load(fh)
    stim, noise = rec.get("stimulus", {}), rec.get("noise", {})
    problems = []
    if rec.get("factor_file_sha256") != FACT[cond]:
        problems.append("factor hash")
    if abs(float(stim.get("hold_nA", NAN)) - hold) > 1e-9:
        problems.append("hold_nA %r" % stim.get("hold_nA"))
    if noise.get("seed") != int(task["noise_seed"]) or abs(float(noise.get("sigma_nA", -1)) - float(task["sigma_nA"])) > 1e-12:
        problems.append("noise seed or sigma")
    if abs(float(rec.get("run_s", -1)) - RUN_S) > 1e-9 or abs(float(rec.get("record_dt_ms", -1)) - DT) > 1e-12:
        problems.append("run length or record step")
    if rec.get("popset_file") not in (None, ""):
        problems.append("popset present")
    if problems:
        sys.exit("STOP: RECORD.json of %s does not match its task: %s" % (folder, "; ".join(problems)))
    t = times(checked(folder, "SPIKES.csv", tab))
    v = np.load(checked(folder, "V_SOMA.npy", tab), allow_pickle=False)
    n_run = int(round(RUN_S * 1000.0 / DT))
    if v.ndim != 1 or len(v) not in (n_run, n_run + 1) or rec.get("n_samples") not in (None, len(v)):
        sys.exit("STOP: V_SOMA.npy in %s has shape %r; expected %d or %d samples, equal to the record's n_samples (%r)"
                 % (folder, v.shape, n_run, n_run + 1, rec.get("n_samples")))
    try:
        dsp = times(checked(folder, "DSPIKES.csv", tab)) if "DSPIKES.csv" in tab else []
        n_dsp = sum(1 for x in dsp if 1000 <= x <= 2000)
    except (ValueError, KeyError):
        n_dsp = -1
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
    i15, i175 = int(round(1500.0 / DT)), int(round(1750.0 / DT))
    seg = v[i15:]
    mean_end = float(np.mean(v[i175:], dtype=np.float64))
    mean_15 = float(np.mean(seg, dtype=np.float64))
    n_end = sum(1 for x in t if x >= 1750)
    n_15 = sum(1 for x in t if x >= 1500)
    state = ("block" if mean_end > -55.0 else "rest") if n_end == 0 else ("pauses" if pauses else "fires")
    state2 = ("block" if mean_15 > -55.0 else "silent") if n_15 == 0 else ("pauses" if pauses else "fires")
    return dict(cond=cond, hold=hold, n=len(t), n_late=len(late), rate=len(late) / 1.0, pauses=pauses,
                last=(max(t) if t else NAN), mean_end=mean_end, mean_15=mean_15, min_v=float(seg.min()),
                max_v=float(seg.max()), trough=trough, dsp=n_dsp, state=state, state2=state2)


def main(root, tasks_tsv, out_tsv, partial):
    if os.path.exists(out_tsv):
        sys.exit("STOP: OUT_TSV exists: " + out_tsv)
    with open(tasks_tsv, "rb") as fh:
        raw = fh.read()
    tasks = list(csv.DictReader(raw.decode("ascii").splitlines(), delimiter="\t"))
    conds = []
    for k in tasks:
        if k["condition"] not in FACT:
            sys.exit("STOP: condition without a pinned factor hash: " + k["condition"])
        if k["condition"] not in conds:
            conds.append(k["condition"])
    rows, missing = [], []
    for k in tasks:
        folder = os.path.join(root, "hold%s_%s" % (k["hold_label"], k["condition"]))
        if not os.path.isfile(os.path.join(folder, "SHA256SUMS")):
            missing.append(os.path.basename(folder))
            continue
        rows.append(one(folder, k))
    print("HOLD_ROOT %s\nTASKS_TSV %s sha256 %s" % (root, tasks_tsv, hashlib.sha256(raw).hexdigest()))
    print("runs in the table: %d; finished: %d; conditions: %d" % (len(tasks), len(rows), len(conds)))
    if missing and not partial:
        sys.exit("STOP: %d runs not finished, first %s (use --partial for a partial print)" % (len(missing), missing[:3]))
    cols = [("condition", "cond", "%s"), ("hold_nA", "hold", "%.3f"), ("spikes_total", "n", "%d"), ("spikes_1_2s", "n_late", "%d"),
            ("rate_1_2s_hz", "rate", "%.1f"), ("pauses_1_2s", "pauses", "%d"), ("last_spike_ms", "last", "%.1f"),
            ("mean_v_1p75_2s", "mean_end", "%.3f"), ("mean_v_1p5_2s", "mean_15", "%.3f"), ("min_v_1p5_2s", "min_v", "%.3f"),
            ("max_v_1p5_2s", "max_v", "%.3f"), ("trough_mv", "trough", "%.3f"), ("dendritic_spikes_1_2s", "dsp", "%d"),
            ("state", "state", "%s"), ("state_v02", "state2", "%s")]
    data = ("\n".join(["\t".join(c[0] for c in cols)] + ["\t".join(c[2] % r[c[1]] for c in cols) for r in rows]) + "\n").encode("ascii")
    with open(out_tsv, "xb") as fh:
        fh.write(data)
    by = dict((c, sorted((r for r in rows if r["cond"] == c), key=lambda r: r["hold"])) for c in conds)
    want = dict((c, sorted(float(k["hold_nA"]) for k in tasks if k["condition"] == c)) for c in conds)
    wt = [(r["rate"], r["trough"]) for r in by.get("WT", []) if r["state"] == "fires" and r["trough"] == r["trough"]]
    wt.sort()

    def wt_trough(rate):
        for (r0, t0), (r1, t1) in zip(wt, wt[1:]):
            if r0 <= rate <= r1 and r1 > r0:
                return t0 + (t1 - t0) * (rate - r0) / (r1 - r0)
        return None

    res = {}
    for c in conds:
        rs = by[c]
        have = dict((round(r["hold"], 6), r) for r in rs)
        seq = "".join(LETTER[have[round(h, 6)]["state"]] if round(h, 6) in have else "-" for h in want[c])
        block = [r["hold"] for r in rs if r["state"] == "block"]
        first = min(block) if block else None
        stays = all(r["state"] == "block" for r in rs if first is not None and r["hold"] >= first)
        fires = [r["hold"] for r in rs if r["state"] == "fires"]
        dif = []
        for r in rs:
            if r["state"] == "fires" and r["trough"] == r["trough"] and c != "WT":
                w = wt_trough(r["rate"])
                if w is not None:
                    dif.append(r["trough"] - w)
        res[c] = dict(seq=seq, first=first, stays=stays, top=(max(fires) if fires else None), complete=len(rs) == len(want[c]),
                      dif=dif, n=len(rs), low=rs[0] if rs else None)
    w = res.get("WT")
    w_first = w["first"] if w and w["complete"] else None
    print("currents (nA): %s" % " ".join("%.3f" % h for h in want[conds[0]]))
    print("states by current: F fires, P pauses, B block, R rest, - not finished")
    print("\t".join(["condition", "states", "lowest_block_nA", "block_at_every_higher_current", "relative_to_WT", "highest_fires_nA",
                     "trough_minus_WT_at_same_rate_mV(median;n)", "rate_at_lowest_current_hz", "trough_at_lowest_current_mV", "finished"]))
    for c in conds:
        x = res[c]
        rel = "%.3f" % (x["first"] / w_first) if (x["first"] is not None and w_first) else "na"
        dif = "%.2f;%d" % (statistics.median(x["dif"]), len(x["dif"])) if x["dif"] else "na;0"
        low = x["low"]
        print("\t".join([c, x["seq"], "%.3f" % x["first"] if x["first"] is not None else "none",
                         str(x["stays"]) if x["first"] is not None else "na", rel,
                         "%.3f" % x["top"] if x["top"] is not None else "none", dif,
                         "%.0f" % low["rate"] if low else "na", "%.2f" % low["trough"] if low else "na",
                         "%d/%d" % (x["n"], len(want[c]))]))
    odd = [(r["cond"], r["hold"], r["state"], r["state2"]) for r in rows if r["state"] != r["state2"] and not (r["state"] == "rest" and r["state2"] == "silent")]
    print("runs whose state differs from the v02 rule: %d %s" % (len(odd), " ".join("%s@%.3f:%s/%s" % o for o in odd[:20])))
    print("runs in state block with a last spike after 1,000 ms: %s"
          % " ".join("%s@%.3f:%.0fms" % (r["cond"], r["hold"], r["last"]) for r in rows if r["state"] == "block" and r["last"] == r["last"] and r["last"] > 1000))
    print("OUT_TSV %s: %d rows, sha256 %s" % (out_tsv, len(rows), hashlib.sha256(data).hexdigest()))


if __name__ == "__main__":
    a = sys.argv[1:]
    if len(a) == 4 and a[3] == "--partial":
        main(a[0], a[1], a[2], True)
    elif len(a) == 3:
        main(a[0], a[1], a[2], False)
    else:
        sys.exit(__doc__)
