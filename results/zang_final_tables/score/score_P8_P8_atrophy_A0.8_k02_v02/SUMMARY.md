# Zang 2018 P4 metrics, E1 definitions on 10-25 s inclusive (25-s runs; PRIMARY used 30-45 s)

rate = window spikes / 15 s; CV ddof 0; CV2 adjacent-ISI; pause = share of ISIs > 3 x median; log-ISI mode = fullest 0.05-decade bin on [-1, 5]; dspikes = distal upward -30 mV crossings, 20 ms refractory; main = same rule at the smooth-dendrite site.

| condition | n | rate Hz | CV | CV2 | pause | log-ISI mode ms | dspikes | main crossings | wall h |
|---|---:|---|---|---|---|---|---|---|---|
| P8_atrophy_A0.8 | 1 | 55.47 (55.47 to 55.47) | 0.02342 (0.02342 to 0.02342) | 0.02783 (0.02783 to 0.02783) | 0 (0 to 0) | 18.84 (18.84 to 18.84) | 0 (0 to 0) | 0 (0 to 0) | 13.87 (13.87 to 13.87) |

Per run.

| folder | group | seed | spikes | rate Hz | CV | CV2 | pause | log-ISI mode ms | dspikes | main crossings | wall h |
|---|---|---|---:|---|---|---|---|---|---:|---:|---|
| P8_atrophy_A0.8_k02 | P8_atrophy_A0.8 | 20261015 | 832 | 55.47 | 0.02342 | 0.02783 | 0 | 18.84 | 0 | 0 | 13.87 |
