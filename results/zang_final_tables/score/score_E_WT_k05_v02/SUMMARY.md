# Zang 2018 P4 metrics, E1 definitions on 10-25 s inclusive (25-s runs; PRIMARY used 30-45 s)

rate = window spikes / 15 s; CV ddof 0; CV2 adjacent-ISI; pause = share of ISIs > 3 x median; log-ISI mode = fullest 0.05-decade bin on [-1, 5]; dspikes = distal upward -30 mV crossings, 20 ms refractory; main = same rule at the smooth-dendrite site.

| condition | n | rate Hz | CV | CV2 | pause | log-ISI mode ms | dspikes | main crossings | wall h |
|---|---:|---|---|---|---|---|---|---|---|
| WT | 1 | 41.8 (41.8 to 41.8) | 0.02038 (0.02038 to 0.02038) | 0.02371 (0.02371 to 0.02371) | 0 (0 to 0) | 23.71 (23.71 to 23.71) | 0 (0 to 0) | 0 (0 to 0) | 12.69 (12.69 to 12.69) |

Per run.

| folder | group | seed | spikes | rate Hz | CV | CV2 | pause | log-ISI mode ms | dspikes | main crossings | wall h |
|---|---|---|---:|---|---|---|---|---|---:|---:|---|
| WT_k05 | WT | 20261018 | 627 | 41.8 | 0.02038 | 0.02371 | 0 | 23.71 | 0 | 0 | 12.69 |
