"""Zang P9 auxiliary-subunit factors v02 (pre-specified) for Zang driver v05.
v02 = v01 with the beta4 labels set to the plan's definition: r = Kcnmb4 ratio = 0.38 = fraction of abBK that
remains; abBK gabkbar x r, mslo gbar + (1 - r) x abBK gabkbar (values in the files unchanged, r = 0.38).
(i) a2d2: newCaP pcabar x 0.85; (iii) KChIP1, partial mapping, inactivation-curve shift only: Kv4 and Kv4s
vshift +2 (V1/2 = vh - vshift, so V1/2 -2 mV); (b4) beta4 shift abBK -> mslo, r = 0.38 remains; (iv) aux-only =
(i) + (b4) with the five channels at WT; (v) five+aux = KI_18w_mid five-gene factors + (i) + (b4).
Kvbeta1 has no counterpart: Zang Kv1 inactivation (hinf, tauh = 1000/qt ms) is hard-coded, not RANGE.
Refuses to overwrite.
python -I -B 04_model_zang/zp9_make_factors.py <WT.json> <KI_18w_mid.json> <OUT_DIR>
"""
import hashlib, json, sys
from pathlib import Path

PINS = {"WT.json": "bdbb0f458069e46aefc13f6dfe0dd488de397fcf1546821dca2c3b804157c204",
        "KI_18w_mid.json": "660f8462aec5b6671a279c5041ee539909b8a8395de3c7eed30994ca283c937a"}
wt, ki, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
for p in (wt, ki):
    if hashlib.sha256(p.read_bytes()).hexdigest() != PINS[p.name]:
        raise SystemExit(f"hash mismatch: {p}")
W, K = json.loads(wt.read_text()), json.loads(ki.read_text())
A2D2 = {"newCaP.pcabar": 0.85}
KCHIP = {"Kv4.vshift": 2.0, "Kv4s.vshift": 2.0}
B4 = {"beta4_shift_r": 0.38}
conds = [("ZP9_i_a2d2", W, A2D2, "alpha2delta2 single: newCaP pcabar x 0.85"),
         ("ZP9_iii_kchip1", W, KCHIP, "KChIP1, partial mapping, inactivation-curve shift only: Kv4/Kv4s V1/2 -2 mV"),
         ("ZP9_b4_shift", W, B4, "beta4 shift single: abBK gabkbar x r, mslo gbar + (1 - r) x abBK gabkbar, r = 0.38"),
         ("ZP9_iv_aux_only", W, {**A2D2, **B4}, "aux-only: (i) + beta4 shift, five channels at WT"),
         ("ZP9_v_five_aux", K, {**A2D2, **B4}, "five+aux: KI_18w_mid five-gene factors + (i) + beta4 shift")]
names = []
for name, base, aux, label in conds:
    doc = {"condition": name, "label": label, "base_file": "WT.json" if base is W else "KI_18w_mid.json",
           "base_sha256": PINS["WT.json" if base is W else "KI_18w_mid.json"], "scales": base["scales"], "aux": aux}
    with open(out / f"{name}.json", "x", encoding="utf-8", newline="") as f:
        f.write(json.dumps(doc, indent=2) + "\n")
    names.append(f"{name}.json")
with open(out / "SHA256SUMS", "x", encoding="utf-8", newline="") as f:
    f.write("".join(f"{hashlib.sha256((out / n).read_bytes()).hexdigest()}  {n}\n" for n in names))
print("\n".join(names))

#%% debug
