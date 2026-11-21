# Figure 1

Overview of the study design. Figure 1 is a drawing, not a plot of data: `Figure_1.svg` is its source, with the text kept as
text. It shows what was done, not the results; the numbers it quotes (cohort, gene panel, parameter sets) are those of
the Methods. The font Source Sans 3 (Adobe; SIL Open Font License 1.1) is embedded
in the file.

## Render

```
cd 06_figures/figure_1
node render_figure_1.js Figure_1.svg /path/to/Figure_1
```

Writes `Figure_1.png` (2067 x 2766 px, that is 175 mm wide at 300 dpi) and `Figure_1.pdf` (vector, 175 mm x 234.167 mm).
Existing files are not overwritten. Needs Node.js with the package `playwright` and its Chromium; the figure in the
manuscript was rendered with Node.js 22.23.4, playwright 1.63.0 and Chromium 153.0.8010.12. `expected_output/` holds
the two files of the manuscript with their SHA-256 sums.
