"""Figure 5, panel A: one second of somatic voltage of the published parameter set, WT and with all five transcript
ratios, without and with synaptic input.

Usage: python3 figure5_panel_A.py TRACE_DIR OUT_PREFIX

TRACE_DIR holds the four trace files named below (time in ms, somatic voltage in mV; 10,001 rows each).
Writes OUT_PREFIX.png and OUT_PREFIX.pdf. Spikes are counted as upward crossings of -20 mV."""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
D, OUT = sys.argv[1:3]
files = [("Without synaptic input","WT","trace_noinput_WT_sig0p0065_n50261113_5000-6000ms.tsv"),
         ("Without synaptic input","All five transcript ratios","trace_noinput_KI_18w_mid_sig0p0065_n50261113_5000-6000ms.tsv"),
         ("With synaptic input","WT","trace_syn_WT_set000_5000-6000ms.tsv"),
         ("With synaptic input","All five transcript ratios","trace_syn_KI_18w_mid_set000_5000-6000ms.tsv")]
plt.rcParams.update({"font.size":11,"font.family":"DejaVu Sans"})
fig, axes = plt.subplots(2,2, figsize=(9.2,4.7), dpi=1950/9.2, sharex=True, sharey=True)
for ax,(cond,lab,fn) in zip(axes.ravel(), files):
    a = np.loadtxt(os.path.join(D,fn), skiprows=1)
    t, v = a[:,0]-a[0,0], a[:,1]
    assert a.shape==(10001,2)
    n = int(np.sum((v[:-1] < -20) & (v[1:] >= -20)))
    print(fn, "rows", a.shape[0], "spikes(-20 mV upward crossings)", n, "min %.1f max %.1f" % (v.min(), v.max()))
    col = "#333333" if lab=="WT" else "#c0392b"
    ax.plot(t, v, color=col, lw=0.7)
    ax.set_title(f"{lab}\n{cond.lower()}, {n} spikes in 1 s", fontsize=11, loc="left", linespacing=1.15)
    for s in ("top","right"): ax.spines[s].set_visible(False)
    ax.set_xlim(0,1000); ax.set_ylim(-75,10); ax.set_yticks([-75,-50,-25,0])
for ax in axes[:,0]: ax.set_ylabel("Somatic voltage (mV)")
for ax in axes[1,:]: ax.set_xlabel("Time (ms)")
fig.text(0.005,0.985,"A",fontsize=14,fontweight="bold",va="top")
fig.suptitle("Published parameter set: one second of somatic voltage", x=0.04, y=0.985, ha="left", va="top", fontsize=11)
fig.tight_layout(rect=(0.012,0,1,0.955))
fig.savefig(OUT+".png"); fig.savefig(OUT+".pdf")
