# Zang 2018 P4 metrics, E1 definitions on 10-25 s inclusive (25-s runs; PRIMARY used 30-45 s)

rate = window spikes / 15 s; CV ddof 0; CV2 adjacent-ISI; pause = share of ISIs > 3 x median; log-ISI mode = fullest 0.05-decade bin on [-1, 5]; dspikes = distal upward -30 mV crossings, 20 ms refractory; main = same rule at the smooth-dendrite site.

| condition | n | rate Hz | CV | CV2 | pause | log-ISI mode ms | dspikes | main crossings | wall h |
|---|---:|---|---|---|---|---|---|---|---|
| KI_30w_mid | 1 | 63.07 (63.07 to 63.07) | 0.02714 (0.02714 to 0.02714) | 0.03346 (0.03346 to 0.03346) | 0 (0 to 0) | 16.79 (16.79 to 16.79) | 0 (0 to 0) | 0 (0 to 0) | 13.28 (13.28 to 13.28) |

Per run.

| folder | group | seed | spikes | rate Hz | CV | CV2 | pause | log-ISI mode ms | dspikes | main crossings | wall h |
|---|---|---|---:|---|---|---|---|---|---:|---:|---|
| KI_30w_mid_k05 | KI_30w_mid | 20261018 | 946 | 63.07 | 0.02714 | 0.03346 | 0 | 16.79 | 0 | 0 | 13.28 |
