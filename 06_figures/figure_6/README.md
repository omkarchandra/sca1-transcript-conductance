# Figure 6

Firing rate ratio with the Cav2.1 conductance unchanged and reduced (A), and CV ratio against firing rate ratio (B), across
the versions of the Purkinje cell model.

| Panel | Script | Reads |
|---|---|---|
| A, B | `figure6_panels.py` | `inputs/versions_without_synaptic_input.json`, `inputs/versions_with_synaptic_input.tsv`, the reference values |
| whole figure | `figure6_page.py` | the panel file |

The reference values are in `data/reference_values/bushart2021_reference_values.md`.

## Run

```
cd 06_figures/figure_6
bash make_figure_6.sh /path/to/a/new/folder
```

The folder must not exist yet; no script overwrites a file. The run writes the panel file, the plotted values
(`panels.tsv`) and `Figure_6.png` (1950 x 2613 px, 300 dpi) with `Figure_6_vector.pdf` (468 x 627 pt).
`expected_output/` holds the files of the figure in the manuscript, with their SHA-256 sums.

Needs Python 3 with matplotlib, numpy, pillow and pypdf, and the fonts Liberation Sans and DejaVu Sans (the second is
shipped with matplotlib). The figure in the manuscript was drawn with Python 3.13.18, matplotlib 3.11.2, numpy 2.5.3,
pillow 12.3.0, pypdf 6.17.0 and FreeType 2.14.3; with these versions the output is byte-identical to
`expected_output/`. With other versions the values file stays the same, and the images can differ in single pixels.

## Inputs

All inputs are small tables scored from the simulation runs; the runs themselves are not in this repository.

- `versions_without_synaptic_input.json`: one entry per version (58 versions), with the five conductance factors, the
  quartiles and the median of the firing rate ratio over the 130 parameter sets, the median CV ratio and the counts
  of parameter sets by state. The figure draws 35 of these versions: the 31 combinations of
  the five transcript ratios and the four staining-ratio versions.
- `versions_with_synaptic_input.tsv`: the same for the 21 versions that were run with synaptic input; the figure
  draws 14 of them: those without a membrane-area correction, apart from the staining version with Cav2.1 at 0.85.

## Names of the versions

The tables use the working names of the versions. The conductance factors of every version are in the tables.

| Name in the tables | Version in the manuscript |
|---|---|
| `KI_18w_mid` | all five conductances at their 18-week KI/WT transcript ratios (Kv3.3 0.667, BK 0.719, Kv4.3 0.588, Cav2.1 0.672, Cav3.1 0.649) |
| `KI_1c_restore_<gene>` | four conductances at their transcript ratios, the named one unchanged |
| `KI_3_only_<genes>`, `KI_2b_only_<genes>` | only the named conductances at their transcript ratios |
| `KI_4_protein_Kcnma1_Cacna1g` | staining ratios: BK 0.78 and Cav3.1 0.75, the others unchanged |
| `..._plus_Kcnc3`, `..._plus_Kcnd3` | staining ratios plus the transcript ratio of Kv3.3 or Kv4.3 |
| `KI_5_protein_..._Cacna1a_0p85` | staining ratios with the Cav2.1 conductance at 0.85 (also 0.95, 0.90, 0.80) |
| `KI_6_...`, `KI_7_...`, `KI_8_...` | versions with a change of an auxiliary subunit (column `aux`) |
| `P8_atrophy_A0.8` | membrane-area correction with A = 0.8 (also 0.9, 0.85, 0.75, 0.7, 0.6) |

Genes and channels: Kcnc3, Kv3.3; Kcnma1, BK; Kcnd3, Kv4.3; Cacna1a, Cav2.1; Cacna1g, Cav3.1.
