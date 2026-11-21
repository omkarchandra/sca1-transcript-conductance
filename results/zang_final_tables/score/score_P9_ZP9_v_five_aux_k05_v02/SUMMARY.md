# Zang 2018 P4 metrics, E1 definitions on 10-25 s inclusive (25-s runs; PRIMARY used 30-45 s)

rate = window spikes / 15 s; CV ddof 0; CV2 adjacent-ISI; pause = share of ISIs > 3 x median; log-ISI mode = fullest 0.05-decade bin on [-1, 5]; dspikes = distal upward -30 mV crossings, 20 ms refractory; main = same rule at the smooth-dendrite site.

| condition | n | rate Hz | CV | CV2 | pause | log-ISI mode ms | dspikes | main crossings | wall h |
|---|---:|---|---|---|---|---|---|---|---|
| ZP9_v_five_aux | 1 | 120.6 (120.6 to 120.6) | 0.03611 (0.03611 to 0.03611) | 0.04406 (0.04406 to 0.04406) | 0 (0 to 0) | 8.414 (8.414 to 8.414) | 0 (0 to 0) | 0 (0 to 0) | 16.09 (16.09 to 16.09) |

Per run.

| folder | group | seed | spikes | rate Hz | CV | CV2 | pause | log-ISI mode ms | dspikes | main crossings | wall h |
|---|---|---|---:|---|---|---|---|---|---:|---:|---|
| ZP9_v_five_aux_k05 | ZP9_v_five_aux | 20261018 | 1809 | 120.6 | 0.03611 | 0.04406 | 0 | 8.414 | 0 | 0 | 16.09 |
