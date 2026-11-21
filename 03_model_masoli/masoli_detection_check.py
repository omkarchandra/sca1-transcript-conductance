"""Read-only detection-margin check for the Masoli P8 (28) and P9 (70) runs..
The runs saved SPIKES.csv only (no voltage trace; driver v06/v10 writes V_TRACE only with --save-v), so spikes cannot be
re-detected at -10 mV. What exists per run:
  SPIKES.csv: every upward 0 mV crossing (2 ms refractory) with its peak within 5 ms;
  RECORD.json: n_spikes_0mV and n_spikes_m20mV, upward crossings of 0 and -20 mV counted by the driver on the live trace
  over the whole 45 s.
Reported per run: spikes peaking below +5 mV (whole run and 30-45 s window), the minimum peak in the window, the 30-45 s
spike count at 0 mV, and n_spikes_m20mV - n_spikes_0mV (whole run). An event that crosses -10 mV but not 0 mV also crosses
-20 mV, so a zero difference at -20 mV bounds the -10 mV difference at zero for the whole run. Runs with a non-zero
difference or any peak below +5 mV are listed. Writes ONE new folder (CHECK.tsv, SHA256SUMS); no other file is touched.
python -I -B 03_model_masoli/masoli_detection_check.py <OUT_DIR>
"""

# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

PUB = Path(_release_paths.expand('${DATA_ROOT}/publication'))
ROOTS = {"P8": (PUB / 'masoli_P8_atrophy', 28), "P9": (PUB / 'masoli_P9_auxiliary', 70)}


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    out = Path(sys.argv[1])
    if out.exists():
        raise SystemExit(f"output exists: {out}")
    rows = []
    for tag, (root, n) in ROOTS.items():
        folders = sorted(p for p in root.iterdir() if p.is_dir())
        if len(folders) != n:
            raise SystemExit(f"{tag}: expected {n} folders, found {len(folders)}")
        for d in folders:
            for line in (d / "SHA256SUMS").read_text().splitlines():
                h, f = line.split("  ", 1)
                if sha256(d / f) != h:
                    raise SystemExit(f"checksum mismatch: {d / f}")
            sp = list(csv.DictReader(open(d / "SPIKES.csv")))
            t = np.array([float(x["time_ms"]) for x in sp])
            p = np.array([float(x["peak_mV"]) for x in sp])
            w = (t >= 30000) & (t <= 45000)
            rec = json.load(open(d / "RECORD.json"))
            n0, n20 = int(rec["n_spikes_0mV"]), int(rec["n_spikes_m20mV"])
            if n0 != len(t):
                raise SystemExit(f"SPIKES.csv rows differ from n_spikes_0mV: {d}")
            rows.append(dict(set=tag, run=d.name, n_spikes_0mV_run=n0, n_spikes_m20mV_run=n20, diff_m20_minus_0=n20 - n0,
                             n_window_0mV=int(w.sum()), peaks_below_5mV_run=int((p < 5.0).sum()),
                             peaks_below_5mV_window=int((p[w] < 5.0).sum()), min_peak_window_mV=repr(float(p[w].min())),
                             spikes_sha256=sha256(d / "SPIKES.csv")))
    out.mkdir()
    with open(out / "CHECK.tsv", "x", newline="") as f:
        w_ = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w_.writeheader()
        w_.writerows(rows)
    with open(out / "SHA256SUMS", "x") as f:
        f.write(f"{sha256(out / 'CHECK.tsv')}  CHECK.tsv\n")
    flag = [r for r in rows if r["diff_m20_minus_0"] != 0 or r["peaks_below_5mV_run"] > 0]
    print(f"runs: {len(rows)}; with -20 vs 0 mV difference: {sum(r['diff_m20_minus_0'] != 0 for r in rows)}; "
          f"with a peak below +5 mV: {sum(r['peaks_below_5mV_run'] > 0 for r in rows)}")
    print("lowest window peak:", min(rows, key=lambda r: float(r["min_peak_window_mV"]))["run"],
          min(float(r["min_peak_window_mV"]) for r in rows))
    for r in flag:
        print(r)


if __name__ == "__main__":
    main()

#%% debug
