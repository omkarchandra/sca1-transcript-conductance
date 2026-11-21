# Recorded values used as the reference phenotype

The firing phenotype that the simulated neurons are compared with is the one published for Purkinje cells of
Atxn1 154Q/2Q knock-in (KI) mice and wild-type (WT) littermates:

Bushart DD, Huang H, Man LJ, Morrison LM, Shakkottai VG (2021). A chlorzoxazone-baclofen combination improves
cerebellar impairment in spinocerebellar ataxia type 1. Mov Disord 36(3):622–631. doi:10.1002/mds.28355 (PMC7987844).

Mice were 14 weeks old; recordings were made in acute cerebellar slices at 33 °C, separately in the anterior
cerebellum (lobules IV and V) and in the nodular zone (lobule X).

## Group means (their Supplementary Table 1)

Mean ± SEM, two-tailed t-test, copied from Supplementary Table 1 of Bushart et al. (2021). The table gives no units;
the figure axes of the paper give spikes per second, mV and megaohm. The scripts of Figures 5 and 6 read the first two rows
of this table.

| Measure | WT anterior | KI anterior | P | WT nodular | KI nodular | P |
|---|---|---|---|---|---|---|
| Firing frequency | 31.65 ± 2.96 | 42.04 ± 4.47 | 0.1359 | 27.81 ± 1.76 | 23.78 ± 3.12 | 0.255 |
| CV of the interspike interval | 0.1171 ± 0.015 | 0.2082 ± 0.019 | <0.001 | 0.1145 ± 0.011 | 0.3274 ± 0.050 | <0.001 |
| Afterhyperpolarization minimum | −65.93 ± 1.30 | −63.88 ± 0.415 | 0.0073 | −66.65 ± 0.676 | −65.03 ± 0.967 | 0.0259 |
| Input resistance | 53.51 ± 6.66 | 57.21 ± 2.85 | 0.5529 | 73.97 ± 7.11 | 79.96 ± 5.11 | 0.5058 |

## KI relative to WT (computed from the table above)

| Measure | Anterior | Nodular |
|---|---|---|
| Firing frequency, KI / WT | 1.33 | 0.86 |
| CV of the interspike interval, KI / WT | 1.78 | 2.86 |
| Afterhyperpolarization minimum, KI − WT | +2.05 mV | +1.62 mV |
| Input resistance, KI / WT | 1.07 | 1.08 |

The firing-rate ratios 0.86 and 1.33 are the bounds of the recorded firing-rate range of the manuscript and of Figures 5 and 6.

## Values read from the published figures

The paper gives the threshold to depolarization block and the staining of channel proteins in figures only. These
values were measured on the figure images of the author manuscript (PMC7987844): the axis tick marks were located in
the image and used as the scale, and the bar ends were read by pixel row. The reading error is about one pixel
(13 pA for the block threshold, 0.004 for the staining).

| Measure | WT anterior | KI anterior | WT nodular | KI nodular |
|---|---|---|---|---|
| Threshold to depolarization block, pA (their Figure 1E) | 1,260 | 853 | 1,654 | 840 |

KI / WT: 0.68 (anterior) and 0.51 (nodular).

| Staining in the molecular layer, KI as a fraction of WT (their Figure 2) | Anterior | Nodular |
|---|---|---|
| BK | 0.79 | 0.77 |
| Cav3.1 | 0.77 | 0.73 |
| Cav2.1 | 1.01 | 0.93 |

The staining-ratio versions of the model use BK 0.78 and Cav3.1 0.75, the means of the two regions.

## Membrane capacitance of Purkinje cells in the same mouse line

Watase K, et al. (2002). A long CAG repeat in the mouse Sca1 locus replicates SCA1 features and reveals the impact of
protein solubility on selective neurodegeneration. Neuron 34:905–919. Their Table 1 (whole-cell recordings at
6–11 weeks, mean ± SEM):

| Parameter | WT (n = 14) | Sca1 154Q/2Q (n = 16) | KI / WT (computed here) |
|---|---|---|---|
| C1, soma and main proximal dendrites, pF | 184 ± 19 | 164 ± 16 | 0.89 |
| C2, dendritic tree, pF | 1038 ± 60 | 776 ± 50 | 0.75 |
| C1 + C2, whole cell, pF | 1222 | 940 | 0.77 |

The membrane-area versions of the model (A = 0.9 to 0.6) bracket these values.
