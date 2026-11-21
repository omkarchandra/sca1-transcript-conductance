# Zang 2018 P4 metrics, E1 definitions on 10-25 s inclusive (25-s runs; PRIMARY used 30-45 s)

rate = window spikes / 15 s; CV ddof 0; CV2 adjacent-ISI; pause = share of ISIs > 3 x median; log-ISI mode = fullest 0.05-decade bin on [-1, 5]; dspikes = distal upward -30 mV crossings, 20 ms refractory; main = same rule at the smooth-dendrite site.

| condition | n | rate Hz | CV | CV2 | pause | log-ISI mode ms | dspikes | main crossings | wall h |
|---|---:|---|---|---|---|---|---|---|---|
| P8_atrophy_A0.6 | 1 | 38.4 (38.4 to 38.4) | 0.01787 (0.01787 to 0.01787) | 0.02071 (0.02071 to 0.02071) | 0 (0 to 0) | 26.61 (26.61 to 26.61) | 0 (0 to 0) | 0 (0 to 0) | 14.35 (14.35 to 14.35) |

Per run.

| folder | group | seed | spikes | rate Hz | CV | CV2 | pause | log-ISI mode ms | dspikes | main crossings | wall h |
|---|---|---|---:|---|---|---|---|---|---:|---:|---|
| P8_atrophy_A0.6_k01 | P8_atrophy_A0.6 | 20261014 | 576 | 38.4 | 0.01787 | 0.02071 | 0 | 26.61 | 0 | 0 | 14.35 |
