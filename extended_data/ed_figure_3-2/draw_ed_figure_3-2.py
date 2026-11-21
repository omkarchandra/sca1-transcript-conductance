import numpy as np, matplotlib, hashlib, os
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
# KI/WT staining ratio, molecular layer, 14 weeks (Bushart et al. 2021, read from their figure): anterior, nodular
stain = {'BK': (0.79, 0.77), 'Cav3.1': (0.77, 0.73), 'Cav2.1': (1.01, 0.93), 'IP3 receptor': (0.97, 0.98)}
# KI - WT log2 fold change, Purkinje nuclei, 18 weeks (Extended Data Figure 3-1)
l2 = {'BK': -0.476, 'Cav3.1': -0.623, 'Cav2.1': -0.573, 'IP3 receptor': -0.762}
fig, a = plt.subplots(figsize=(4.2, 4.0))
a.plot([0.5, 1.1], [0.5, 1.1], color='#999999', lw=0.8, ls='--')
a.text(0.905, 0.88, 'equal ratios', color='#777777', fontsize=7.5, ha='center', va='center', rotation=45)
a.axhline(1, color='#cccccc', lw=0.6); a.axvline(1, color='#cccccc', lw=0.6)
xs, ys = [], []
for n in stain:
    x = 2 ** l2[n]; y = float(np.mean(stain[n])); xs.append(x); ys.append(y)
    a.plot([x, x], list(stain[n]), color='#1f3b57', lw=1.2); a.plot(x, y, 'o', color='#1f3b57', ms=6)
    left = n == 'IP3 receptor'
    a.annotate(n, (x, y), xytext=((-4, 9) if left else (8, -3)), textcoords='offset points', fontsize=9, ha=('center' if left else 'left'))
a.set_xlim(0.5, 1.1); a.set_ylim(0.5, 1.1); a.set_aspect('equal')
a.set_xlabel('KI/WT transcript ratio (Purkinje nuclei, 18 weeks)'); a.set_ylabel('KI/WT staining ratio (molecular layer, 14 weeks)')
fig.tight_layout()
OUT = 'Extended_Data_Figure_3-2_2026-11-18.png'; assert not os.path.exists(OUT)
fig.savefig(OUT, dpi=300); fig.savefig('Extended_Data_Figure_3-2_2026-11-18_vector.pdf')
r = float(np.corrcoef(xs, ys)[0, 1]); print('points', [(round(x, 3), round(y, 3)) for x, y in zip(xs, ys)], 'Pearson r', round(r, 3))
from PIL import Image
im = Image.open(OUT); print(im.size, hashlib.sha256(open(OUT, 'rb').read()).hexdigest(), os.path.getsize(OUT))
