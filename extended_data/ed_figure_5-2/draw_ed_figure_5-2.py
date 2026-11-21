"""Extended Data Figure 5-2: one second of somatic voltage of the published parameter set with synaptic input,
for WT and five versions (v6 of 20 Nov 2026: spikes counted at -30 mV as in the Methods; titles give the 1-s count only). Input: data/TRACES_v4_with_synaptic_input.tsv and data/RECORD.json."""
import sys, os, json, hashlib
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
tsv, rec_p, out = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(out, exist_ok=False)
assert hashlib.sha256(open(tsv,"rb").read()).hexdigest()=="7c0f178a624a5c1b926c7a9ac6af43db065d92de605d332c0185614e02250ea0"
rec=json.load(open(rec_p))
D=np.genfromtxt(tsv,delimiter="\t",names=True)
t=D["time_ms"]-5000.0
ROWS=[("WT","WT"),
 ("all_five_transcript_ratios","All five transcript ratios"),
 ("Cav2.1_unchanged","Four transcript ratios, Cav2.1 unchanged"),
 ("staining_ratios","Staining ratios"),
 ("area_A0.8","Five transcript ratios divided by A = 0.8"),
 ("area_A0.75","Five transcript ratios divided by A = 0.75")]
names=D.dtype.names
plt.rcParams.update({"font.size":10,"font.family":"DejaVu Sans"})
fig,axes=plt.subplots(len(ROWS),1,figsize=(7.2,8.4),sharex=True,sharey=True)
counts={}
for ax,(key,lab) in zip(axes,ROWS):
    col=[n for n in names if n.replace("_","")==key.replace(".","").replace("_","")]
    assert len(col)==1,(key,col)
    v=D[col[0]]
    n=int(np.sum((v[1:]>=-30)&(v[:-1]<-30))); counts[key]=n
    assert n==rec["inputs"][key]["spikes_window"],(key,n)
    ax.plot(t,v,color="#333333" if key=="WT" else "#c0392b",lw=0.6)
    ax.set_title(f"{lab}: {n} spikes in 1 s",fontsize=9.5,loc="left",pad=3)
    ax.set_ylim(-75,12); ax.set_yticks([-50,0])
    for s in ("top","right"): ax.spines[s].set_visible(False)
axes[-1].set_xlim(0,1000); axes[-1].set_xlabel("Time (ms)")
fig.text(0.01,0.5,"Somatic voltage (mV)",rotation=90,va="center")
fig.suptitle("Published parameter set, with synaptic input: one second of somatic voltage",x=0.08,ha="left",fontsize=11)
fig.tight_layout(rect=(0.03,0,1,0.97))
fig.savefig(f"{out}/Extended_Data_Figure_5-2.png",dpi=300); fig.savefig(f"{out}/Extended_Data_Figure_5-2.pdf")
json.dump({"counts_window":counts,"png_sha256":hashlib.sha256(open(f"{out}/Extended_Data_Figure_5-2.png","rb").read()).hexdigest()},open(f"{out}/RECORD.json","x"),indent=1)
print(counts)
#%% debug
