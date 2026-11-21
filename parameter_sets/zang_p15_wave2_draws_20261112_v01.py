#!/usr/bin/env python3
"""P15 wave 2 parameter-set draws. Latin hypercube, numpy only.

N = 200 sets over the eight P13 conductance groups (same order and same popset JSON format as P13, so the P13 driver's
--popset path reads them unchanged). Each factor is spread evenly (one value per stratum, stratum order shuffled
independently per group) over [1/1.5, 1.5] on a log scale, i.e. a box centred on the published values.
RNG: numpy.random.default_rng(20261112), only Generator.random() is used (its stream is the 53-bit double of the
PCG64 output and is the same in every numpy version; no permutation or uniform helpers). For group j (fixed order below):
r = random(200); stratum rank = rank of r (stable argsort); jitter = random(200); u = (rank + jitter) / 200;
factor = exp(ln(1/1.5) + u * (ln 1.5 - ln(1/1.5))). Set k (1..200) gets set_id k and noise seed 40261112 + k, stored as the
popset "seed" (the launchers pass the same number as --noise-seed).
Writes into a NEW exclusive-create folder: popset_kNNN.json (200), POPSETS.tsv, SHA256SUMS. Checks before writing: every
group has exactly one value in each of the 200 strata; every factor is inside [1/1.5, 1.5] and inside the driver's [0.5, 1.5].
Usage: python3 -I -B zang_p15_wave2_draws_20261112_v01.py OUT_DIR
"""
import hashlib, json, math, sys
from pathlib import Path

import numpy as np

N_SETS = 200
LHS_SEED = 20261112
NOISE_BASE = 40261112
LO, HI = 1.0 / 1.5, 1.5
GROUPS = ["Kv3", "Kv4", "BK", "newCaP", "CaT3_1", "naRsg", "hpkj", "SK2"]   # P13 order: five disease + three pacemaking
RULE = ("Latin hypercube, 200 strata per group, default_rng(20261112).random() only; per group: ranks of random(200) as strata, "
        "u=(rank+random(200))/200, factor=exp(ln(1/1.5)+u*2*ln(1.5))")


def draw():
    rng = np.random.default_rng(LHS_SEED)
    u = np.empty((N_SETS, len(GROUPS)))
    for j in range(len(GROUPS)):
        r = rng.random(N_SETS)
        rank = np.empty(N_SETS, dtype=np.int64)
        rank[np.argsort(r, kind="stable")] = np.arange(N_SETS)
        u[:, j] = (rank + rng.random(N_SETS)) / N_SETS
    x = np.exp(math.log(LO) + u * (math.log(HI) - math.log(LO)))
    for j in range(len(GROUPS)):
        assert sorted(np.floor(u[:, j] * N_SETS).astype(int).tolist()) == list(range(N_SETS)), ("strata", GROUPS[j])
    assert bool(np.all(x >= 0.5) and np.all(x <= 1.5) and np.all(x >= LO * (1 - 1e-12)) and np.all(x <= HI * (1 + 1e-12)))
    return x


def main(out):
    out = Path(out)
    x = draw()
    out.mkdir(parents=False, exist_ok=False)
    rows = []
    for k in range(1, N_SETS + 1):
        scales = {g: float(x[k - 1, j]) for j, g in enumerate(GROUPS)}
        rec = {"set_id": k, "seed": NOISE_BASE + k, "scales": scales, "wave": 2, "lhs_seed": LHS_SEED, "rule": RULE, "numpy": np.__version__}
        with (out / f"popset_k{k:03d}.json").open("x") as f:
            f.write(json.dumps(rec, indent=1, sort_keys=True) + "\n")
        rows.append([k, NOISE_BASE + k] + [f"{scales[g]:.17g}" for g in GROUPS])
    with (out / "POPSETS.tsv").open("x") as f:
        f.write("\t".join(["set_id", "seed"] + GROUPS) + "\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    with (out / "SHA256SUMS").open("x") as f:
        for p in sorted(out.iterdir()):
            if p.name != "SHA256SUMS":
                f.write(hashlib.sha256(p.read_bytes()).hexdigest() + "  " + p.name + "\n")
    print(f"wrote {N_SETS} sets to {out}")
    print("POPSETS.tsv sha256", hashlib.sha256((out / "POPSETS.tsv").read_bytes()).hexdigest())
    print("SHA256SUMS  sha256", hashlib.sha256((out / "SHA256SUMS").read_bytes()).hexdigest())


if __name__ == "__main__":
    main(sys.argv[1])


#%% debug
