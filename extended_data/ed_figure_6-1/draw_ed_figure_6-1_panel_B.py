"""Extended Data Figure 6-1, panel B redrawn for the v6 (bioRxiv) version, 20 Nov 2026.
Second rendering (b): smaller fonts and markers so that the x labels clear the key. Only change: the shaded band is the recorded rate range (0.86-1.33) instead of 0.75-1.25.
Panel A is kept as placed (top 897 rows of ED_Figure_6-1_2026-11-18.png); panel B is redrawn
from the per-run source table of the Masoli population (E_SUPPLEMENTARY_TABLE_S4.tsv)."""
import sys, os, csv, json, hashlib
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from PIL import Image
src_png, src_tsv, out = sys.argv[1], sys.argv[2], sys.argv[3]
assert hashlib.sha256(open(src_png,"rb").read()).hexdigest()=="15837b8a70065be57321e63a64faa9b485e1af40163a4c84fe116e8432cdd8c3"
assert hashlib.sha256(open(src_tsv,"rb").read()).hexdigest()=="d849ffb9fa0468788a7d3a5002d00dd36a1a13bf4e54e0ecd971e6732a413522"
os.makedirs(out, exist_ok=False)
COND=[("KI_18w_mid","unrestored KI",False),
 ("KI_18w_mid_restore_Kcnc3","restore Kcnc3/Kv3",False),
 ("KI_18w_mid_restore_Kcnma1","restore Kcnma1/BK",False),
 ("KI_18w_mid_restore_Kcnd3","restore Kcnd3/Kv4",False),
 ("KI_18w_mid_restore_Cacna1a","restore Cacna1a/CaP",True),
 ("KI_18w_mid_restore_Cacna1g","restore Cacna1g/CaT",False),
 ("KI_18w_mid_grp_BK_CaP","Kcnma1/BK + Cacna1a/CaP",True),
 ("KI_18w_mid_grp_Kv3_Kv4_CaT","Kcnc3/Kv3 + Kcnd3/Kv4 + Cacna1g/CaT",False),
 ("KI_18w_mid_grp_BK_CaT","Kcnma1/BK + Cacna1g/CaT",False),
 ("KI_18w_mid_grp_Kv3_Kv4_CaP","Kcnc3/Kv3 + Kcnd3/Kv4 + Cacna1a/CaP",True)]
rows=list(csv.DictReader(open(src_tsv),delimiter="\t"))
vals={}
for key,_,_ in COND:
    r=[x for x in rows if x["condition"]==key and x["stage"]=="primary"]
    assert len(r)==14,(key,len(r))
    vals[key]={int(x["member"]):float(x["rate_ratio"]) for x in r}
# spot checks against the workbook (Figure 6-4 sheet, 3 decimals)
for key,m,v in [("KI_18w_mid",2,1.147),("KI_18w_mid_restore_Cacna1a",2,1.045),("KI_18w_mid_grp_BK_CaP",9,0.997),("KI_18w_mid_grp_Kv3_Kv4_CaT",26,1.203),("KI_18w_mid_restore_Cacna1g",26,1.226)]:
    assert abs(round(vals[key][m],3)-v)<1e-9,(key,m,vals[key][m])
im=Image.open(src_png).convert("RGB"); W,H=im.size; YCUT=897
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8})
fig=plt.figure(figsize=(W/300,(H-YCUT)/300),dpi=300)
ax=fig.add_axes([218/W, 1-876/(H-YCUT), (2050-218)/W, (876-96)/(H-YCUT)])
ax.axhspan(0.86,1.33,color="#eceef2",zorder=0,lw=0)
ax.axhline(1.0,color="#555555",lw=1.0,zorder=1)
navy="#1d3557"; orange="#f4a261"
for i,(key,lab,cap) in enumerate(COND):
    v=np.array([vals[key][m] for m in sorted(vals[key])])
    order=np.argsort(v); off=np.zeros(14)
    pattern=[0,-0.08,0.08,-0.04,0.04,-0.12,0.12]
    for rank,j in enumerate(order): off[j]=pattern[rank%7]
    ax.scatter(i+off,v,s=20,color=orange if cap else navy,edgecolor="#222222",linewidth=0.5,zorder=3)
    ax.plot([i-0.22,i+0.22],[np.median(v)]*2,color="#222222",lw=2.0,zorder=4)
ax.set_xlim(-0.55,len(COND)-0.45); ax.set_ylim(0.80,1.38)
ax.set_yticks([0.8,0.9,1.0,1.1,1.2,1.3])
ax.set_xticks(range(len(COND))); ax.set_xticklabels([c[1] for c in COND],rotation=40,ha="right",rotation_mode="anchor")
ax.set_ylabel("KI-scaled / WT settings, firing rate",fontsize=8.5)
ax.set_title("KI 18 weeks",fontsize=9.5)
for s in ("top","right"): ax.spines[s].set_visible(False)
fig.text(60/W,1-38/(H-YCUT),"B",fontsize=13,fontweight="bold",va="center")
fig.legend(handles=[Patch(color=orange,label="orange = Cacna1a/CaP put back to WT in that run"),
                    Patch(color="#eceef2",label="recorded rate range (0.86–1.33)")],
           loc="lower center",ncol=2,frameon=False,bbox_to_anchor=(0.5,0.0),handlelength=2.2,fontsize=8.5)
fig.savefig(f"{out}/panelB.png",dpi=300)
plt.close(fig)
pb=Image.open(f"{out}/panelB.png").convert("RGB"); assert pb.size==(W,H-YCUT),pb.size
comp=Image.new("RGB",(W,H),"white"); comp.paste(im.crop((0,0,W,YCUT)),(0,0)); comp.paste(pb,(0,YCUT))
comp.save(f"{out}/Extended_Data_Figure_6-1_v6.png",dpi=(300,300))
rec={"source_png":src_png,"source_tsv":src_tsv,"ycut":YCUT,"medians":{k:float(np.median(list(vals[k].values()))) for k,_,_ in COND},
     "png_sha256":hashlib.sha256(open(f"{out}/Extended_Data_Figure_6-1_v6.png","rb").read()).hexdigest()}
json.dump(rec,open(f"{out}/RECORD.json","x"),indent=1)
print(json.dumps(rec["medians"],indent=0)); print(rec["png_sha256"])
#%% debug
