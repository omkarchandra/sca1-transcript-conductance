# Figure 5

Simulated neurons scaled by the KI/WT transcript ratios, compared with the firing recorded in KI Purkinje cells
(Bushart et al. 2021).

| Panel | Script | Reads |
|---|---|---|
| A | `figure5_panel_A.py` | `inputs/traces/` (four voltage traces of one second each) |
| B | `figure5_panel_B.py` | `inputs/per_set_values.txt`, the reference values |
| whole figure | `figure5_join.py` | the two panel files |

The reference values are in `data/reference_values/bushart2021_reference_values.md`.

## Run

```
cd 06_figures/figure_5
bash make_figure_5.sh /path/to/a/new/folder
```

The folder must not exist yet; no script overwrites a file. The run writes the two panel files, the plotted values
(`panel_B.tsv`) and `Figure_5.png` (1950 x 2341 px, 300 dpi) with `Figure_5_vector.pdf` (468 x 562 pt).
`expected_output/` holds the files of the figure in the manuscript, with their SHA-256 sums.

Needs Python 3 with matplotlib, numpy, pillow and pypdf, and the fonts Liberation Sans and DejaVu Sans (the second is
shipped with matplotlib). The figure in the manuscript was drawn with Python 3.13.18, matplotlib 3.11.2, numpy 2.5.3,
pillow 12.3.0, pypdf 6.17.0 and FreeType 2.14.3; with these versions the output is byte-identical to
`expected_output/`. With other versions the values file stays the same, and the images can differ in single pixels.

## Inputs

All inputs are small tables scored from the simulation runs; the runs themselves are not in this repository.

- `per_set_values.txt`: one value per parameter set for four versions and the six membrane-area versions, each
  without and with synaptic input: the firing rate ratio (RATE) and the CV ratio (CV) to three decimals, the
  afterhyperpolarization difference in mV (AHP), the WT firing rate in Hz (WTHZ) and the state of the KI-scaled and
  of the WT run as one letter (F fires, P fires with pauses, B depolarization block, R at rest). `x` marks a set
  without a value, and `= first` a line that equals the same line of the first version of its block. Set 0 is the
  published parameter set. Each version carries the SHA-256 of the per-set table it was taken from and the sums
  of its values, which the script checks.
- `traces/`: somatic voltage of the published parameter set from 5 to 6 s of the run, sampled every 0.1 ms, for WT
  and for the version with all five transcript ratios, without and with synaptic input.

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
