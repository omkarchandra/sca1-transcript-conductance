# E7 metrics, E1 definitions (30-45 s inclusive), paired against each member's PRIMARY WT

rate = window spikes / 15 s; CV = sd(ISI)/mean(ISI), ddof 0; CV2 = mean(2|dISI|/(ISI_i+ISI_i+1)); pause fraction = share of window ISIs > 3 x median ISI. d_* = condition minus that member's WT.

| condition | n | rate Hz median (IQR) | d_rate Hz | CV median | d_CV | CV2 median | d_CV2 | pause median | d_pause |
|---|---:|---|---|---|---|---|---|---|---|
| amp1.5_KI_24w_mid | 14 | 92.5 (80.1 to 112) | 11.9 (10.5 to 16) | 0.021 (0.0145 to 0.0263) | -0.00865 (-0.0144 to -0.00467) | 0.0214 (0.0136 to 0.0289) | -0.0106 (-0.0194 to -0.00529) | 0 (0 to 0) | 0 (0 to 0) |
| amp2_KI_24w_mid | 14 | 101 (82.2 to 119) | 17.1 (12.6 to 21.8) | 0.0174 (0.0136 to 0.0244) | -0.0101 (-0.018 to -0.00531) | 0.0169 (0.0126 to 0.0262) | -0.0124 (-0.024 to -0.00559) | 0 (0 to 0) | 0 (0 to 0) |
| amp3_KI_24w_mid | 14 | 107 (87.2 to 122) | 20 (14 to 28.6) | 0.015 (0.0121 to 0.0212) | -0.0125 (-0.0217 to -0.00626) | 0.0142 (0.0111 to 0.0218) | -0.0147 (-0.0283 to -0.00691) | 0 (0 to 0) | 0 (0 to 0) |
