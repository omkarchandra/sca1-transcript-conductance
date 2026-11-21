# Zang 2018 P4 metrics, E1 definitions on 10-25 s inclusive (25-s runs; PRIMARY used 30-45 s)

rate = window spikes / 15 s; CV ddof 0; CV2 adjacent-ISI; pause = share of ISIs > 3 x median; log-ISI mode = fullest 0.05-decade bin on [-1, 5]; dspikes = distal upward -30 mV crossings, 20 ms refractory; main = same rule at the smooth-dendrite site.

| condition | n | rate Hz | CV | CV2 | pause | log-ISI mode ms | dspikes | main crossings | wall h |
|---|---:|---|---|---|---|---|---|---|---|
| KI_18w_mid | 1 | 79.73 (79.73 to 79.73) | 0.0278 (0.0278 to 0.0278) | 0.03283 (0.03283 to 0.03283) | 0 (0 to 0) | 11.89 (11.89 to 11.89) | 0 (0 to 0) | 0 (0 to 0) | 13.36 (13.36 to 13.36) |

Per run.

| folder | group | seed | spikes | rate Hz | CV | CV2 | pause | log-ISI mode ms | dspikes | main crossings | wall h |
|---|---|---|---:|---|---|---|---|---|---:|---:|---|
| KI_18w_mid_k02 | KI_18w_mid | 20261015 | 1196 | 79.73 | 0.0278 | 0.03283 | 0 | 11.89 | 0 | 0 | 13.36 |
