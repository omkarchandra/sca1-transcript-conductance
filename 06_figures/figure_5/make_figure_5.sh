#!/usr/bin/env bash
# Draws Figure 5 from the tables in inputs/.
# Usage (from this folder):  bash make_figure_5.sh OUT_DIR      OUT_DIR must not exist yet.
set -euo pipefail
OUT="${1:?usage: bash make_figure_5.sh OUT_DIR (a folder that does not exist yet)}"
mkdir "$OUT"
REF=../../data/reference_values/bushart2021_reference_values.md

python3 figure5_panel_A.py inputs/traces "$OUT/panel_A"

python3 figure5_panel_B.py inputs/per_set_values.txt "$REF" "$OUT/panel_B" \
    "KI_18w_mid=All five|ratios" \
    "KI_1c_restore_Cacna1a=Cav2.1|unchanged" \
    "KI_4_protein_Kcnma1_Cacna1g=Staining|ratios" \
    "P8_atrophy_A0.8=Area,|A = 0.8" \
    "P8_atrophy_A0.75=Area,|A = 0.75"

python3 figure5_join.py "$OUT/panel_A.png" "$OUT/panel_B.png" "$OUT/panel_A.pdf" "$OUT/panel_B.pdf" "$OUT/Figure_5"
