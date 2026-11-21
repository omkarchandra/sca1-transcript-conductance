# Zang 2018 P4 metrics, E1 definitions on 10-25 s inclusive (25-s runs; PRIMARY used 30-45 s)

rate = window spikes / 15 s; CV ddof 0; CV2 adjacent-ISI; pause = share of ISIs > 3 x median; log-ISI mode = fullest 0.05-decade bin on [-1, 5]; dspikes = distal upward -30 mV crossings, 20 ms refractory; main = same rule at the smooth-dendrite site.

| condition | n | rate Hz | CV | CV2 | pause | log-ISI mode ms | dspikes | main crossings | wall h |
|---|---:|---|---|---|---|---|---|---|---|
| KI_18w_mid | 1 | 79.73 (79.73 to 79.73) | 0.02657 (0.02657 to 0.02657) | 0.03178 (0.03178 to 0.03178) | 0 (0 to 0) | 11.89 (11.89 to 11.89) | 0 (0 to 0) | 0 (0 to 0) | 13.1 (13.1 to 13.1) |

Per run.

| folder | group | seed | spikes | rate Hz | CV | CV2 | pause | log-ISI mode ms | dspikes | main crossings | wall h |
|---|---|---|---:|---|---|---|---|---|---:|---:|---|
| KI_18w_mid_k03 | KI_18w_mid | 20261016 | 1196 | 79.73 | 0.02657 | 0.03178 | 0 | 11.89 | 0 | 0 | 13.10 |
