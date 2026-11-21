# Masoli 2015 population generator . Local Python only: no NEURON, no model run, no cluster.
# Draws N members; each channel group gets one scale per member, log-uniform in [0.5, 1.5] (+/-50 %).
# Published values and section keys come from PC_param.py (read by path; plain dict, no code run beyond it).
# Usage: python -I -B <this file> <PC_param.py path> <out TSV> <seed> <N>   (writes a NEW file; refuses if it exists)
import csv
import hashlib
import importlib.util
import json
import math
import os
import sys
from pathlib import Path

import numpy as np

# channel group -> PC_param keys (section) that the group's scale multiplies
GROUPS = {
    "Kv3_3": ["Kv3.3Dend"],
    "Kca1_1": ["Kca1.1Soma", "Kca1.1Dend"],
    "Kv4_3": ["Kv4.3Dend"],
    "Cav2_1": ["Cav2.1Soma", "Cav2.1Dend", "Cav2.1AIS"],
    "Cav3_1": ["Cav3.1Soma", "Cav3.1Dend", "Cav3.1Ais"],
    "Nav1_6": ["Nav1.6Soma", "Nav1.6Dend", "Nav1.6AIS"],
    "HCN1": ["HCNSoma", "HCNDend"],
    "Kca2_2": ["Kca2.2Soma", "Kca2.2Dend"],
}
LO, HI = 0.5, 1.5


def main():
    pc_path, out_path, seed, n = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    spec = importlib.util.spec_from_file_location("pc_param_mod", pc_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    pc = mod.pc_param
    rng = np.random.default_rng(seed)
    groups = list(GROUPS)
    rows = []
    for m in range(n):
        u = rng.uniform(math.log(LO), math.log(HI), size=len(groups))
        scales = np.exp(u)
        rows.append([m] + [f"{s:.6f}" for s in scales])
    with open(out_path, "x", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["member_id"] + groups)
        w.writerows(rows)
    rec = {
        "seed": seed, "N": n, "distribution": "log-uniform in [0.5, 1.5] per group per member",
        "pc_param_sha256": hashlib.sha256(pc_path.read_bytes()).hexdigest(),
        "groups": {g: {k: pc[k] for k in keys} for g, keys in GROUPS.items()},
        "members_tsv_sha256": hashlib.sha256(out_path.read_bytes()).hexdigest(),
    }
    rec_path = out_path.with_suffix(".json")
    with open(rec_path, "x", encoding="utf-8", newline="") as f:
        f.write(json.dumps(rec, indent=1) + "\n")
    print(rec["members_tsv_sha256"], out_path, rec_path)


if __name__ == "__main__":
    main()

#%% debug
