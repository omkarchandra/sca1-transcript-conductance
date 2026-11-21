"""P8 atrophy factors (pre-specified): the five-gene KI_18w_mid factors divided by A (0.8, 0.6).
Other genes stay as in the source (all 1.0). Writes two new JSON files and a SHA256SUMS; refuses to overwrite.
python -I -B 03_model_masoli/p8_make_factors.py <KI_18w_mid.json> <OUT_DIR>
"""
import hashlib, json, sys
from pathlib import Path

FIVE = ("Kcnc3", "Kcnma1", "Kcnd3", "Cacna1a", "Cacna1g")
SRC_SHA = "660f8462aec5b6671a279c5041ee539909b8a8395de3c7eed30994ca283c937a"
src, out = Path(sys.argv[1]), Path(sys.argv[2])
if hashlib.sha256(src.read_bytes()).hexdigest() != SRC_SHA:
    raise SystemExit("source factor file hash mismatch")
doc = json.loads(src.read_text())
others = {g: v for g, v in doc["scales"].items() if g not in FIVE}
if any(float(v) != 1.0 for v in others.values()):
    raise SystemExit(f"non-five genes not 1.0: {others}")
names = []
for a in (0.8, 0.6):
    scales = {g: (doc["scales"][g] / a if g in FIVE else doc["scales"][g]) for g in doc["scales"]}
    new = {"condition": f"P8_atrophy_A{a}", "age_weeks": doc["age_weeks"], "atrophy_A": a,
           "source_file": "KI_18w_mid.json", "source_sha256": SRC_SHA,
           "rule": "five-gene factor / A; other genes unchanged", "scales": scales}
    name = f"P8_atrophy_A{a}.json"
    with open(out / name, "x", encoding="utf-8", newline="") as f:
        f.write(json.dumps(new, indent=2) + "\n")
    names.append(name)
    print(name, {g: round(scales[g], 6) for g in FIVE})
with open(out / "SHA256SUMS", "x", encoding="utf-8", newline="") as f:
    f.write("".join(f"{hashlib.sha256((out / n).read_bytes()).hexdigest()}  {n}\n" for n in names))

#%% debug
