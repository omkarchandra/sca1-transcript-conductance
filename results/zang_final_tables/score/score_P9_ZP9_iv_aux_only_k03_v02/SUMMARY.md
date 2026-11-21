# Zang 2018 P4 metrics, E1 definitions on 10-25 s inclusive (25-s runs; PRIMARY used 30-45 s)

rate = window spikes / 15 s; CV ddof 0; CV2 adjacent-ISI; pause = share of ISIs > 3 x median; log-ISI mode = fullest 0.05-decade bin on [-1, 5]; dspikes = distal upward -30 mV crossings, 20 ms refractory; main = same rule at the smooth-dendrite site.

| condition | n | rate Hz | CV | CV2 | pause | log-ISI mode ms | dspikes | main crossings | wall h |
|---|---:|---|---|---|---|---|---|---|---|
| ZP9_iv_aux_only | 1 | 59.33 (59.33 to 59.33) | 0.02695 (0.02695 to 0.02695) | 0.03294 (0.03294 to 0.03294) | 0 (0 to 0) | 16.79 (16.79 to 16.79) | 0 (0 to 0) | 0 (0 to 0) | 14.35 (14.35 to 14.35) |

Per run.

| folder | group | seed | spikes | rate Hz | CV | CV2 | pause | log-ISI mode ms | dspikes | main crossings | wall h |
|---|---|---|---:|---|---|---|---|---|---:|---:|---|
| ZP9_iv_aux_only_k03 | ZP9_iv_aux_only | 20261016 | 890 | 59.33 | 0.02695 | 0.03294 | 0 | 16.79 | 0 | 0 | 14.35 |
