# Zang 2018 P4 metrics, E1 definitions on 10-25 s inclusive (25-s runs; PRIMARY used 30-45 s)

rate = window spikes / 15 s; CV ddof 0; CV2 adjacent-ISI; pause = share of ISIs > 3 x median; log-ISI mode = fullest 0.05-decade bin on [-1, 5]; dspikes = distal upward -30 mV crossings, 20 ms refractory; main = same rule at the smooth-dendrite site.

| condition | n | rate Hz | CV | CV2 | pause | log-ISI mode ms | dspikes | main crossings | wall h |
|---|---:|---|---|---|---|---|---|---|---|
| WT | 1 | 41.93 (41.93 to 41.93) | 0.01988 (0.01988 to 0.01988) | 0.02216 (0.02216 to 0.02216) | 0 (0 to 0) | 23.71 (23.71 to 23.71) | 0 (0 to 0) | 0 (0 to 0) | 15.84 (15.84 to 15.84) |

Per run.

| folder | group | seed | spikes | rate Hz | CV | CV2 | pause | log-ISI mode ms | dspikes | main crossings | wall h |
|---|---|---|---:|---|---|---|---|---|---:|---:|---|
| WT_k02 | WT | 20261015 | 629 | 41.93 | 0.01988 | 0.02216 | 0 | 23.71 | 0 | 0 | 15.84 |
