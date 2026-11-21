# Zang 2018 P4 metrics, E1 definitions on 10-25 s inclusive (25-s runs; PRIMARY used 30-45 s)

rate = window spikes / 15 s; CV ddof 0; CV2 adjacent-ISI; pause = share of ISIs > 3 x median; log-ISI mode = fullest 0.05-decade bin on [-1, 5]; dspikes = distal upward -30 mV crossings, 20 ms refractory; main = same rule at the smooth-dendrite site.

| condition | n | rate Hz | CV | CV2 | pause | log-ISI mode ms | dspikes | main crossings | wall h |
|---|---:|---|---|---|---|---|---|---|---|
| KI_30w_mid | 1 | 63.13 (63.13 to 63.13) | 0.02651 (0.02651 to 0.02651) | 0.0314 (0.0314 to 0.0314) | 0 (0 to 0) | 14.96 (14.96 to 14.96) | 0 (0 to 0) | 0 (0 to 0) | 14.36 (14.36 to 14.36) |

Per run.

| folder | group | seed | spikes | rate Hz | CV | CV2 | pause | log-ISI mode ms | dspikes | main crossings | wall h |
|---|---|---|---:|---|---|---|---|---|---:|---:|---|
| KI_30w_mid_k04 | KI_30w_mid | 20261017 | 947 | 63.13 | 0.02651 | 0.0314 | 0 | 14.96 | 0 | 0 | 14.36 |
