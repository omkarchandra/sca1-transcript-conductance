# Zang 2018 P4 metrics, E1 definitions on 10-25 s inclusive (25-s runs; PRIMARY used 30-45 s)

rate = window spikes / 15 s; CV ddof 0; CV2 adjacent-ISI; pause = share of ISIs > 3 x median; log-ISI mode = fullest 0.05-decade bin on [-1, 5]; dspikes = distal upward -30 mV crossings, 20 ms refractory; main = same rule at the smooth-dendrite site.

| condition | n | rate Hz | CV | CV2 | pause | log-ISI mode ms | dspikes | main crossings | wall h |
|---|---:|---|---|---|---|---|---|---|---|
| KI_24w_mid | 1 | 64.07 (64.07 to 64.07) | 0.0256 (0.0256 to 0.0256) | 0.03069 (0.03069 to 0.03069) | 0 (0 to 0) | 14.96 (14.96 to 14.96) | 0 (0 to 0) | 0 (0 to 0) | 13.57 (13.57 to 13.57) |

Per run.

| folder | group | seed | spikes | rate Hz | CV | CV2 | pause | log-ISI mode ms | dspikes | main crossings | wall h |
|---|---|---|---:|---|---|---|---|---|---:|---:|---|
| KI_24w_mid_k02 | KI_24w_mid | 20261015 | 961 | 64.07 | 0.0256 | 0.03069 | 0 | 14.96 | 0 | 0 | 13.57 |
