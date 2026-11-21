# Extended Data figures drawn by scripts in this folder

Each script writes into a new folder (or, for Extended Data Figure 3-2, into the current folder) and does not
overwrite files. `expected_output/` holds the image of the figure as it stands in the Extended Data, with its SHA-256.
The figures were drawn with Python 3.13.18, matplotlib 3.11.2, numpy 2.5.3 and pillow 12.3.0; with these versions the
images are byte-identical to `expected_output/`.

| Figure | Script | Inputs |
|---|---|---|
| Extended Data Figure 3-2 | `ed_figure_3-2/draw_ed_figure_3-2.py` | values written in the script: the KI/WT staining ratios read from the figure of Bushart et al. (2021), and the 18-week log2 fold changes of Extended Data Figure 3-1 |
| Extended Data Figure 5-2 | `ed_figure_5-2/draw_ed_figure_5-2.py` | `inputs/TRACES_v4_with_synaptic_input.tsv` (one second of somatic voltage per version, with synaptic input) and `inputs/RECORD.json` (the spike counts it is checked against) |
| Extended Data Figure 6-1 | `ed_figure_6-1/draw_ed_figure_6-1_panel_B.py` | `inputs/E_SUPPLEMENTARY_TABLE_S4.tsv` (per-run values of the population of the second model) for panel B; panel A is taken unchanged from the earlier rendering `inputs/ED_Figure_6-1_2026-11-18.png` |

## Run

```
cd extended_data/ed_figure_3-2 && mkdir /path/to/new && cd /path/to/new && python3 <this folder>/draw_ed_figure_3-2.py
cd extended_data/ed_figure_5-2 && python3 draw_ed_figure_5-2.py inputs/TRACES_v4_with_synaptic_input.tsv inputs/RECORD.json /path/to/new/folder
cd extended_data/ed_figure_6-1 && python3 draw_ed_figure_6-1_panel_B.py inputs/ED_Figure_6-1_2026-11-18.png inputs/E_SUPPLEMENTARY_TABLE_S4.tsv /path/to/new/folder
```
