# Code for the SCA1 transcript-to-conductance study

This package holds the analysis code, small frozen inputs and result tables behind the manuscript. Large inputs are not included (see "Data").

## Contents

| Folder | What it holds |
|---|---|
| `01_snrnaseq/` | Single-nucleus RNA-seq analyses of mouse cerebellum: count splitting, the 12-week held-out check, far-set occupancy, RNA content, PCA, limma tests, lane checks, spread and state-space analyses |
| `02_replication/` | Replication in public data: Purkinje-resolution datasets (`a1_purkinje/`) and bulk cerebellum (`p2_bulk/`) |
| `03_model_masoli/` | Population runs of the Purkinje cell model of ModelDB 229585: generator, drivers and scorers |
| `04_model_zang/` | The Purkinje cell model of ModelDB 243446: build, drivers, scorer and table builder |
| `05_model_forrest/` | The two-compartment model of ModelDB 180789: gene-to-conductance tables, scale factors, Stage A and Stage B drivers and scorers |
| `06_figures/` | Figure scripts |
| `06_figures/figure_1/`, `figure_5/`, `figure_6/` | The builders of Figures 1, 5 and 6, with their inputs and expected outputs. `figure_1/LICENCES.md` lists the third-party drawings and font. |
| `extended_data/` | The builders of Extended Data Figures 2-3, 3-2, 5-2 and 6-1 (panel B), with their inputs and expected outputs, and the trace extractor of Extended Data Figure 5-2 |
| `data/reference_values/` | The published reference values read by the Figure 5 and 6 builders |
| `config/` | Path configuration (`paths.py`, `paths.R`, `local.yaml.example`) |
| `data/` | Small frozen inputs and records: sample table, unit metadata, limma results, model build records, scale-factor files, `INPUT_MANIFEST.tsv` |
| `results/` | Frozen result tables and run records used by the figures and the text |
| `drivers/`, `factors/`, `parameter_sets/`, `tasks/`, `scoring/`, `tables/`, `masoli/`, `models/`, `provenance/`, `environment/` | The later population and parameter-search runs of the ModelDB 243446 and 266864 models and the matching ModelDB 229585 scoring: drivers, scale factors, parameter sets, task tables, scorers, table builders, model checksum lists and the script that built the run environment |
| `tests/smoke.sh` | Syntax check of every Python, shell and R file |
| `PROVENANCE.tsv` | SHA-256 values for every file (see below) |

## Setup

1. Python 3 with PyYAML is needed to read the configuration.
2. Copy `config/local.yaml.example` to `config/local.yaml` and set every key. `config/paths.py` stops if any key is missing or empty. A key set as an environment variable takes precedence.
   - `DATA_ROOT`: the folder with the external inputs (see "Data").
   - `OUTPUT_ROOT`: where new outputs are written. Scripts create new folders and do not overwrite.
   - `MODEL_ROOT`: the compiled NEURON models.
   - `PYTHON_NEURON`: a Python interpreter with NEURON 8.2.6.
   - `WORK_ROOT`: a workspace for intermediate files.
   - `CONDA_BIN` and `CONDA_ENV`: the conda installation and environment of the state-space analysis.
   - The R scripts are run with `Rscript` from the shell.
   - `REPO_ROOT` is not set: `config/paths.py` takes it to be the package folder.
3. The scripts that launched the population and parameter-search runs (the launchers) are not included. Each run called its driver once per task, with the values of that task as command-line arguments. For `drivers/zang_p13_driver_v01.py`, `drivers/zang_p16_driver_v01.py` and the Masoli population drivers, these values were one line of a table, chosen by the task number. In the wave 2 WT runs, the task number was the parameter-set number. For the Zang drivers v03 and v05 (in `04_model_zang/` and `drivers/`) and the Masoli calibration runs, the values were written in the launcher itself. The Masoli atrophy runs took the population member from a table and the condition from the launcher. The packaged tables, parameter sets and factor files are in `tasks/`, `parameter_sets/`, `masoli/parameter_sets/`, `factors/` and `data/model_masoli/`. These values were set only in the launchers: the conditions and seeds of the Zang E, P8 and P9 runs, the noise levels of the Masoli calibration, the condition names of the Masoli atrophy runs, and the run length and noise level of each Masoli population run. The drivers in `04_model_zang/` and `drivers/` write the environment variables `TASK_ID` (the task number) and `JOB_ID`, if they are set, into each run record.
4. The later runs (the folders in the second-to-last row of the table above) were written for the cluster layout of the original runs. They name their folders as `${PROJECT_ROOT}`, `${SCRATCH_ROOT}`, `${LOCAL_ROOT}` and `${RESULTS_ROOT}`. These are not read from `config/local.yaml`. Replace them, or export them for the shell scripts, before rerunning.
5. `bash tests/smoke.sh` checks that every file parses. It runs no analysis.

## Data

Public data used:
- GEO GSE246183: the mouse cerebellar single-nucleus RNA-seq data (the sample table in `data/snrnaseq/sample_table.tsv` is taken from its Table S1) and its human samples, used in the replication.
- GEO GSE215336 and GSE269430: Purkinje-resolution replication. GSE269430 (SCA7) is exploratory.
- GEO GSE75778 and GSE122099: bulk cerebellum.
- ModelDB 229585, 243446, 266864 and 180789: the neuron models.

Not in the package, and needed to rerun:
- under `DATA_ROOT`: the annotated AnnData file, the derived feature, representation and pseudobulk folders, the 12-week data and the helper modules;
- under `MODEL_ROOT`: the compiled models (`data/model_builds/` and `models/` hold only their build records and checksums);
- the outputs of the simulation runs. The tables in `results/` are the frozen copies used in the paper.

The five-age limma-voom run that gives the 12-week values (`01_snrnaseq/s08_limma/e3_limma_12w.R`) reads these external inputs, listed with their as-run sha256. The first five are pinned in the script. The three 12-week files are not, and their values are those of the copies the run was given.
- `raw_counts.tsv.gz`: 3597cca1e4f690125b414ad48fe517c0f913d2f119f2741d118b43c451942778
- `unit_metadata.tsv`: 11f96a5eaac14dcd523f65935bf255aa18caf2a734597dafee55919b3e97abb4
- `gene_metadata.tsv`: 765cf0b8c4ae14711a87c3f64612b6144ca00451bfe6bc60d3fab7652ffa62a3
- `dictionary_feature_registry.tsv`: ea40dde47e744d9930df038458f4683264736d57e442d450cbc5f1d9707056b8
- `S08_S11_release_purkinje_bergmann.tsv.gz`: 40593a55c8538dc7d687cc07eea13b8d640e4fe4d696794879d159e74b625298 (in this package as `data/snrnaseq/S08_S11_release_purkinje_bergmann.tsv.gz`)
- `12w_unit_metadata.tsv`: 1aa64ede36764a53ca9ad0f0a34605c7e101031605e606c40631f1aaebfdd03a
- `12w_purkinje_raw_counts.tsv.gz`: 5e609ccceeb362b15ceda6f9e33072d22875189575398ebe5eb47b9c6205af60
- `12w_bergmann_raw_counts.tsv.gz`: f697b1c5d43258c7cf89faae301bb9aa5b232cca42001dde6f9582b53824af3b

This run used the earlier unit table (`unit_metadata.tsv`, 11f96a5e…), not the packaged `data/snrnaseq/unit_metadata.tsv` (d2a0b8b2…). The two differ in two columns and in the last digit of some values.

## Provenance

`PROVENANCE.tsv` lists every file of the package with three SHA-256 values:
- `sha256_in_package`: the file as distributed;
- `sha256_in_stage`: the file before the cleaning steps, which removed personal names, machine paths and internal notes, and the code that sent jobs to the cluster or copied files from it;
- `sha256_as_run`: the file as it was when the analyses were run, taken from the build records (for some files, the bytes retrieved from the cluster after the runs), or a dash where the file was written for the package and was never run.

Hashes inside the code that check files outside the package keep their as-run values, and all others refer to the packaged files. `data/INPUT_MANIFEST.tsv` records the as-run hashes of its files, and `PROVENANCE.tsv` maps them to the packaged files.

## Notes on particular scripts

- The original run records carry the scheduler's job and task numbers under another field, which the packaged collectors no longer read. The frozen tables in `results/` keep the original task labels.
- The four scripts in `01_snrnaseq/s12_lanes/` and `01_snrnaseq/s12_own_gene/` support only `--real`. The `--fixture` option remains in the code, but no fixture is included.
- `01_snrnaseq/rna_content/p7_purkinje_output_v01.py` is the record of the run whose outputs are in `01_snrnaseq/rna_content/real_v01/`. It reads the dated folder layout of the original project under `DATA_ROOT`, so it cannot be run with the folder names of this release.

## Software

These versions are given by the scripts' run records, or, for the ModelDB 229585 and 243446 simulations, by the build records and environment script in this package. The records of the scripts not named here give no versions.

- R analyses (`01_snrnaseq/s08_limma/`, `01_snrnaseq/count_splitting/a2_dispersion.R`, `02_replication/a1_purkinje/a1_limma.R`): R 4.4.3 with limma 3.62.1 and edgeR 4.4.0, and R 4.5.3 with limma 3.66.0 and edgeR 4.8.2. Each run writes its versions into its output.
- Python analyses (`01_snrnaseq/count_splitting/` Python scripts, `01_snrnaseq/occupancy/far_set.py`, `01_snrnaseq/rna_content/`, and the download, control-set and test scripts of `02_replication/`): Python 3.14.4 with numpy 2.5.3, pandas 3.0.5, scipy 1.18.1 and h5py 3.16.0.
- State-space analysis (`01_snrnaseq/state_space/`, `01_snrnaseq/s12_own_gene/s12_t5_decomp_r1.py`): Python 3.9.23 with scanpy 1.10.3, anndata 0.10.8, numpy 1.26.4, scipy 1.13.1, pandas 2.3.3, h5py 3.13.0 and decoupler 1.9.0, and R 4.3.1 with UCell 2.6.2.
- Simulations of the ModelDB 229585 and 243446 models (`03_model_masoli/`, `04_model_zang/` and the later-run folders): NEURON 8.2.6 under Python 3.10.1 with numpy 1.26.4. Mechanisms were compiled with gcc 8.5.0.
- Simulations of the two-compartment model (`05_model_forrest/` Stage A and Stage B drivers): NEURON 9.0.2 under Python 3.14.4.
- Scoring of the later runs (`scoring/zang_p13_score_v01.py`): Python 3.14.4 with numpy 2.3.5.
- Figures (`06_figures/`): matplotlib 3.11.1.
- Builders of Figures 1, 5 and 6 and of Extended Data Figures 3-2, 5-2 and 6-1 (`06_figures/figure_1/`, `figure_5/`, `figure_6/`, `extended_data/`): Python 3.13.18 with matplotlib 3.11.2, numpy 2.5.3, pillow 12.3.0 and pypdf 6.17.0. Figure 1 is rendered with Node.js 22.23.4, playwright 1.63.0 and Chromium 153.0.8010.12. Re-run with these versions, each builder reproduced its expected output byte for byte.

## Figures and tables

Each item of the manuscript is listed with the scripts in this package that draw it or build its table.
- "Numbers from" names the script that wrote the plotted or tabulated values.
- Where marked "inferred", that link comes from the scripts' descriptions, not from a run record.
- The Extended Data workbook is not written by any script in this package.

| Item | Scripts |
|---|---|
| **Figure 1** | `06_figures/figure_1/`: `Figure_1.svg`, rendered by `render_figure_1.js`. Third-party drawings and font: `06_figures/figure_1/LICENCES.md`. |
| **Figure 2** | `06_figures/fig1_purkinje_only.py`. Numbers from `01_snrnaseq/s12_lanes/s12_offroad_existing_r1.py` and `01_snrnaseq/heldout_12w/`. |
| **Figure 3** | `06_figures/fig3_channel_heatmap_v04.py` (labels as published). Numbers from `01_snrnaseq/s08_limma/s08_limma_real.R`. `06_figures/fig3_channel_heatmap.py` draws the same data with the earlier labels. |
| **Figure 4** | `06_figures/fig4_replication_v07b.py` with `06_figures/replication_values.py`. The asterisks use the q values of `data/snrnaseq/E3_LIMMA_RESULTS.tsv` (held-out row) and the stored adj.P.Val of the `a1_limma.R` runs (`results/a1/primary/NINE_GENES.tsv` and `results/a1/exploratory/NINE_GENES.tsv`). Numbers from `02_replication/a1_purkinje/a1_limma.R`, `a1_pseudobulk.py`, `02_replication/p2_bulk/p2_tests.py` and `01_snrnaseq/s08_limma/e3_limma_12w.R`. `06_figures/fig4_replication.py` draws the same bars with the earlier labels. |
| **Figures 5 and 6** | `06_figures/figure_5/` (`make_figure_5.sh`, `figure5_panel_A.py`, `figure5_panel_B.py`, `figure5_join.py`) and `06_figures/figure_6/` (`make_figure_6.sh`, `figure6_panels.py`, `figure6_page.py`), both with `data/reference_values/bushart2021_reference_values.md`. |
| **Table 2** | `tables/p16_table_130_v01.py` (with synaptic input) and `scoring/p15_score_all_v06.py` (without). |
| **Extended Data Figure 2-1** | `06_figures/s4_heldout_purkinje.py`. Numbers from `01_snrnaseq/heldout_12w/phase4_12w_heldout_check.py`. |
| **Extended Data Figure 2-2** | a sheet of the Extended Data workbook (not drawn). Numbers from `01_snrnaseq/s12_lanes/s12_t4_pc_r1.py` and `s12_offroad_existing_r1.py` (inferred). |
| **Extended Data Figure 2-3** | `extended_data/ed_figure_2-3/figS7_score_densities.py` (drawn as Supplementary Figure S7 of an earlier version), with its inputs in `extended_data/ed_figure_2-3/inputs/` and the packaged `data/snrnaseq/unit_metadata.tsv` and `results/far_set_ratios/real_v02/OCCUPANCY.tsv`. |
| **Extended Data Figure 3-1** | no drawing script. The four-age values and the gene-set test come from `01_snrnaseq/s08_limma/s08_limma_real.R` (inferred). The 12-week values come from the five-age run, `01_snrnaseq/s08_limma/e3_limma_12w.R` (`data/snrnaseq/E3_LIMMA_RESULTS.tsv`, with its inputs listed under "Data"). |
| **Extended Data Figure 3-2** | `extended_data/ed_figure_3-2/draw_ed_figure_3-2.py`. |
| **Extended Data Figure 4-1** | a sheet of the Extended Data workbook (not drawn): the values plotted in Figure 4 (`results/replication/PLOTTED_BARS.tsv`) and the q values behind its asterisks (`data/snrnaseq/E3_LIMMA_RESULTS.tsv`, `results/a1/primary/NINE_GENES.tsv` and `results/a1/exploratory/NINE_GENES.tsv`). Numbers from `02_replication/a1_purkinje/a1_limma.R`, `a1_pseudobulk.py`, `02_replication/p2_bulk/p2_tests.py` and `01_snrnaseq/s08_limma/e3_limma_12w.R`. |
| **Extended Data Figure 5-1** | a sheet of the Extended Data workbook (not drawn). Rows built by `04_model_zang/zang_final_table.py` and `tables/zang_final_table_20261027.py`. |
| **Extended Data Figure 5-2** | `extended_data/ed_figure_5-2/draw_ed_figure_5-2.py`, from the trace table in `extended_data/ed_figure_5-2/inputs/`, which `extended_data/ed_figure_5-2/extract_traces_v4_v01.py` extracted from runs of `drivers/zang_p16_driver_v01.py`. |
| **Extended Data Figure 6-1** | Panel A: `06_figures/fig5_masoli.py` with `fig5a_values.py`. Numbers from `03_model_masoli/p4_t2_population_scorer.py`. Panel B as published (recorded rate range 0.86–1.33): `extended_data/ed_figure_6-1/draw_ed_figure_6-1_panel_B.py`. `06_figures/fig5_masoli.py` with `fig5b_values.py` draws the earlier panel B (0.75–1.25). |
| **Extended Data Figure 6-2** | a sheet of the Extended Data workbook (not drawn). The block column from `scoring/hold_print_v04.py` (as run). |
| **Extended Data Figure 6-3** | a sheet of the Extended Data workbook (not drawn). The block counts from `scoring/hold_pop_print_v01.py` (as run). `scoring/hold_pop_print_v03.py` reproduces its two tables byte for byte. |
| **Extended Data Figure 6-4** | a sheet of the Extended Data workbook (not drawn). Rows built by `03_model_masoli/p4_t2_population_scorer.py`, `03_model_masoli/e7_score.py`, `masoli/scoring/e6b_far_set_cv2_v02.py` and `masoli/scoring/e7_score_20261024_v04_e1.py`. |
| **Extended Data Table 2-1** | `scoring/p15_variant_score_v01.py` and `scoring/p15_score_all_v06.py`. |
| **Extended Data Table 2-2** | `scoring/p16_variant_score_v02.py`, `scoring/p16_score_all_v01.py`, `tables/build_p16_score_all_v03.py` and `tables/p16_table_130_v01.py`. |
| **Extended Data Table 2-3** | `parameter_sets/zang_p15_wave2_draws_20261112_v01.py`, `scoring/zang_p13_score_v01.py`, `scoring/p15_variant_score_v01.py` and `scoring/p15_score_all_v06.py`. |
| **Extended Data Table 2-4** | `scoring/p15_score_all_v06.py` and `scoring/p16_score_all_v01.py`. |

`06_figures/s3_validity_purkinje.py` draws an age-ordering check of an earlier version (Spearman ρ between age and the WT mouse median, and the number of 12-week WT mice between the 5- and 18-week medians). It is not a figure of the manuscript. `06_figures/s1_noise_calibration.py` draws a calibration figure of an earlier version.

## Licence

MIT (see `LICENSE`).
