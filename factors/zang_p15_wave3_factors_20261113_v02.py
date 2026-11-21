#!/usr/bin/env python3
"""P15 wave 3 factor files: every remaining combination of the five reduced channel genes, and four files that use
measured protein ratios (13 Nov 2026).

Writes into OUT_DIR (must not exist) 28 factor files in the byte format of the files, plus WAVE3_VARIANTS.tsv
and SHA256SUMS. Nothing is overwritten (exclusive create). No random numbers. Standard library only.
  24 files KI_3_only_<genes>.json: the named genes at their KI 18-week transcript ratios, the other channels at 1.0.
     These are the 5 single genes, 9 pairs and 10 triples not run before (the pair Kcnma1 + Cacna1g is wave 2b,
     the five four-gene sets are the 1c files, all five genes is KI_18w_mid).
  4 files KI_4_protein_Kcnma1_Cacna1g[_plus_<genes>].json: Kcnma1 0.78 and Cacna1g 0.75, the protein ratios read from
     Figure 2 of Bushart et al. 2021 (mean of the anterior and nodular values); Cacna1a at 1.0 (protein unchanged in
     that figure); Kcnc3 and Kcnd3, whose proteins were not measured there, each at 1.0 or at its transcript ratio.
WAVE3_VARIANTS.tsv lists each variant with its mask (1 = gene changed; P = protein ratios for Kcnma1 and Cacna1g),
its submission group (1 = single genes and the protein files; 2 = pairs and triples with Cacna1a at WT; 3 = pairs
and triples with Cacna1a reduced), file, sha256 and output root.
Usage: python3 -I -B zang_p15_wave3_factors_20261113_v02.py OUT_DIR
"""
import hashlib, itertools, os, sys

GENES = ["Kcnc3", "Kcnma1", "Kcnd3", "Cacna1a", "Cacna1g"]
RATIO = {"Kcnc3": "0.6666093556941", "Kcnma1": "0.7190924725561185", "Kcnd3": "0.5876215270426435",
         "Cacna1a": "0.6722697088324129", "Cacna1g": "0.6491953709552294"}
PROTEIN = {"Kcnma1": "0.78", "Cacna1g": "0.75"}
OTHERS = ["Scn8a", "Scn4b", "Fgf14", "Kcnn2", "Kcna1", "Hcn1"]
SOURCE = "40593a55c8538dc7d687cc07eea13b8d640e4fe4d696794879d159e74b625298"


def text(condition, scales, variant, note):
    lines = ["{", '  "condition": "%s",' % condition, '  "age_weeks": 18,', '  "source_sha256": "%s",' % SOURCE, '  "scales": {']
    keys = GENES + OTHERS
    for i, k in enumerate(keys):
        lines.append('    "%s": %s%s' % (k, scales.get(k, "1.0"), "," if i < len(keys) - 1 else ""))
    lines += ["  },", '  "p15_variant": "%s",' % variant, '  "p15_note": "%s"' % note, "}"]
    return "\n".join(lines) + "\n"


def variants():
    """Rows: (mask, group, name, scales, condition, note, output root)."""
    out = []
    for size in (1, 2, 3):
        for combo in itertools.combinations(GENES, size):
            if combo == ("Kcnma1", "Cacna1g"):
                continue
            mask = "".join("1" if g in combo else "0" for g in GENES)
            tag = "_".join(combo)
            note = "KI_18w_mid factor file, only the named channel genes at their KI ratios, the other channels at WT (wave 3)"
            group = "1" if size == 1 else ("2" if mask[3] == "0" else "3")
            out.append((mask, group, "KI_3_only_" + tag, dict((g, RATIO[g]) for g in combo),
                        "KI_18w_mid_3_only_" + tag, note, "results/wave3_only_" + tag))
    note = ("Kcnma1 0.78 and Cacna1g 0.75 are protein ratios read from Figure 2 of Bushart et al. 2021 (mean of anterior and "
            "nodular); all other channels at WT; the age and source fields are kept only for the file format (wave 3)")
    out.append(("P01001", "1", "KI_4_protein_Kcnma1_Cacna1g", dict(PROTEIN), "KI_18w_mid_4_protein_Kcnma1_Cacna1g", note,
                "results/wave3_protein_Kcnma1_Cacna1g"))
    for extra in (("Kcnc3",), ("Kcnd3",), ("Kcnc3", "Kcnd3")):
        scales = dict(PROTEIN)
        scales.update((g, RATIO[g]) for g in extra)
        mask = "P" + "".join("1" if g in scales else "0" for g in GENES)
        tag = "protein_Kcnma1_Cacna1g_plus_" + "_".join(extra)
        note = ("Kcnma1 0.78 and Cacna1g 0.75 are protein ratios read from Figure 2 of Bushart et al. 2021 (mean of anterior and "
                "nodular); %s at the KI 18-week transcript ratio; all other channels at WT; the age and source fields are kept "
                "only for the file format (wave 3)" % " and ".join(extra))
        out.append((mask, "1", "KI_4_" + tag, scales, "KI_18w_mid_4_" + tag, note, "results/wave3_" + tag))
    return out


def main(out):
    v = variants()
    if len(v) != 28 or len(set(x[2] for x in v)) != 28 or len(set(x[0] for x in v)) != 28 or len(set(x[6] for x in v)) != 28:
        sys.exit("unexpected variant list: %d rows" % len(v))
    os.mkdir(out)
    rows = ["variant\tmask_Kcnc3_Kcnma1_Kcnd3_Cacna1a_Cacna1g\tgroup\tfile\tsha256\toutroot"]
    for mask, group, name, scales, condition, note, outroot in v:
        fn = name + ".json"
        data = text(condition, scales, name, note).encode("ascii")
        with open(os.path.join(out, fn), "xb") as f:
            f.write(data)
        rows.append("\t".join([name, mask, group, fn, hashlib.sha256(data).hexdigest(), outroot]))
    with open(os.path.join(out, "WAVE3_VARIANTS.tsv"), "x") as f:
        f.write("\n".join(rows) + "\n")
    with open(os.path.join(out, "SHA256SUMS"), "x") as f:
        for fn in sorted(os.listdir(out)):
            if fn != "SHA256SUMS":
                with open(os.path.join(out, fn), "rb") as g:
                    f.write(hashlib.sha256(g.read()).hexdigest() + "  " + fn + "\n")
    with open(os.path.join(out, "SHA256SUMS"), "rb") as g:
        print("wrote %d factor files, WAVE3_VARIANTS.tsv and SHA256SUMS into %s" % (len(v), out))
        print("SHA256SUMS sha256 " + hashlib.sha256(g.read()).hexdigest())


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
