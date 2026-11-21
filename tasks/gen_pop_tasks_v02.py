"""Task table for a block-threshold series on the kept parameter sets (14 Nov 2026; v02 of gen_pop1.py).

Usage: python3 -I -B gen_pop_tasks_v02.py WAVE2_KI_TASKS_TSV ITS_SHA256_FILE OUT_TSV CONDITION [CONDITION ...]

Same checks, currents and row format as gen_pop1.py; the conditions come from the command line, in the order given.
OUT_TSV is written by exclusive create. With the five conditions of the first series the output is byte-identical to
that of gen_pop1.py.
"""
import hashlib, re, sys
if len(sys.argv) < 5:
    sys.exit(__doc__)
src, pin, out, conds = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
if len(set(conds)) != len(conds) or any(not re.fullmatch(r"[A-Za-z0-9_.]+", c) for c in conds):
    sys.exit("STOP: a condition is given twice or has a character outside letters, digits, underscore and point")
raw = open(src, "rb").read()
if hashlib.sha256(raw).hexdigest() != open(pin).read().split()[0]:
    sys.exit("STOP: WAVE2_KI_TASKS.tsv does not have the hash of its .sha256 file")
lines = raw.decode("ascii").rstrip("\n").split("\n")
if lines[0].split("\t")[:3] != ["task", "set_id", "condition"]:
    sys.exit("STOP: unexpected header")
sets = [int(l.split("\t")[1]) for l in lines[1:]]
if len(sets) != 130 or len(set(sets)) != 130 or any(l.split("\t")[2] != "KI_18w_mid" for l in lines[1:]):
    sys.exit("STOP: expected 130 different sets, all with condition KI_18w_mid")
holds = ["0.250", "0.500", "0.750", "1.000", "1.250", "1.500", "1.750", "2.000", "2.500", "3.000"]
rows = ["task\tset_id\tcondition\thold_nA\thold_label"]
for c in conds:
    for k in sets:
        for h in holds:
            rows.append("%d\t%d\t%s\t%s\t%s" % (len(rows) - 1, k, c, h, h.replace(".", "p")))
data = ("\n".join(rows) + "\n").encode("ascii")
open(out, "xb").write(data)
print(len(rows) - 1, "tasks;", len(sets), "sets;", len(conds), "conditions;", len(holds), "currents; first set", sets[0], "last set", sets[-1])
print("conditions:", " ".join(conds))
print("sha256", hashlib.sha256(data).hexdigest())
