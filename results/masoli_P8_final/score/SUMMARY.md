# E7 metrics, E1 definitions (30-45 s inclusive), paired against each member's PRIMARY WT

rate = window spikes / 15 s; CV = sd(ISI)/mean(ISI), ddof 0; CV2 = mean(2|dISI|/(ISI_i+ISI_i+1)); pause fraction = share of window ISIs > 3 x median ISI. d_* = condition minus that member's WT.

| condition | n | rate Hz median (IQR) | d_rate Hz | CV median | d_CV | CV2 median | d_CV2 | pause median | d_pause |
|---|---:|---|---|---|---|---|---|---|---|
| P8_atrophy_A0.6 | 14 | 71.6 (64.9 to 94.4) | -3.73 (-4.38 to -3.17) | 0.0363 (0.0207 to 0.0448) | 0.00401 (0.00192 to 0.00551) | 0.0428 (0.0209 to 0.0561) | 0.00503 (0.0022 to 0.00766) | 0 (0 to 0) | 0 (0 to 0) |
| P8_atrophy_A0.8 | 14 | 81.6 (74.2 to 104) | 5.27 (4.75 to 5.9) | 0.027 (0.0166 to 0.0325) | -0.00424 (-0.00626 to -0.00218) | 0.0289 (0.0161 to 0.0371) | -0.00551 (-0.009 to -0.00263) | 0 (0 to 0) | 0 (0 to 0) |
