"""Figure 6 (version of 20 Nov 2026: the panel with Cav2.1 lowered in steps and the points of those versions removed).
A: median firing rate ratio (KI / WT) of the 16 combinations of the other four conductances, with Cav2.1 unchanged and reduced.
B: median CV ratio against median firing rate ratio, one point per version of the model (the 31 combinations of the five
   transcript ratios and the four staining-ratio versions without synaptic input; 14 versions with synaptic input).

Usage: python3 figure6_panels.py VERSIONS_WITHOUT_INPUT_JSON REFERENCE_VALUES_MD OUT_PREFIX VERSIONS_WITH_INPUT_TSV

Reads the scored table of the runs without synaptic input (one entry per version: medians and quartiles over the 130
parameter sets), the scored table of the runs with synaptic input, and the recorded values (firing frequency and CV of
the interspike interval, Bushart et al. 2021) from REFERENCE_VALUES_MD. Writes OUT_PREFIX.png (1950 px wide), OUT_PREFIX.pdf
and OUT_PREFIX.tsv (the plotted values) by exclusive create. The script stops if the tables do not hold the expected
versions, if a text leaves the page or if the two keys meet.
"""
import json, sys
from decimal import Decimal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

table, note, prefix, syn_tsv = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
R = json.load(open(table))
GENES = ["Kcnc3", "Kcnma1", "Kcnd3", "Cacna1a", "Cacna1g"]
OTHER = ["Kcnc3", "Kcnma1", "Kcnd3", "Cacna1g"]
FULL = {"Kcnc3": "0.667", "Kcnma1": "0.719", "Kcnd3": "0.588", "Cacna1a": "0.672", "Cacna1g": "0.649"}
CH = {"Kcnc3": "Kv3.3", "Kcnma1": "BK", "Kcnd3": "Kv4.3", "Cacna1a": "Cav2.1", "Cacna1g": "Cav3.1"}
# ---- recorded values from the reference note
def cells(l):
    return [c.strip() for c in l.strip().strip("|").split("|")]
def num(s):
    return float(s.split("±")[0].replace("−", "-").replace(",", "").strip())
ref = {}
for l in open(note, encoding="utf-8").read().split("\n"):
    x = cells(l)
    if len(x) == 7 and x[0] in ("Firing frequency", "CV of the interspike interval"):
        if x[0] in ref:
            sys.exit("STOP: reference row twice: " + x[0])
        ref[x[0]] = dict(wa=num(x[1]), ka=num(x[2]), wn=num(x[4]), kn=num(x[5]))
if len(ref) != 2:
    sys.exit("STOP: reference rows found: %r" % sorted(ref))
F, C = ref["Firing frequency"], ref["CV of the interspike interval"]
REC = [("anterior", F["ka"] / F["wa"], C["ka"] / C["wa"]), ("nodular", F["kn"] / F["wn"], C["kn"] / C["wn"])]
LO, HI = sorted(r[1] for r in REC)
if "%.2f-%.2f" % (LO, HI) != "0.86-1.33":
    sys.exit("STOP: recorded range")
# ---- versions that are drawn
cube = [v for v, r in R.items() if r["aux"] == "-" and all(r[g] in ("1", FULL[g]) for g in GENES)]
prot = ["KI_4_protein_Kcnma1_Cacna1g", "KI_4_protein_Kcnma1_Cacna1g_plus_Kcnd3", "KI_4_protein_Kcnma1_Cacna1g_plus_Kcnc3", "KI_4_protein_Kcnma1_Cacna1g_plus_Kcnc3_Kcnd3"]
dose = ["KI_5_protein_Kcnma1_Cacna1g_Cacna1a_0p95", "KI_5_protein_Kcnma1_Cacna1g_Cacna1a_0p90", "KI_5_protein_Kcnma1_Cacna1g_Cacna1a_0p85", "KI_5_protein_Kcnma1_Cacna1g_Cacna1a_0p80"]
alone = ["KI_6_aux_a2d2"]
if len(cube) != 31 or any(v not in R for v in prot):
    sys.exit("STOP: versions of the table are not as expected")
def ptype(r):
    f = Decimal(r["Cacna1a"]) * (Decimal("0.85") if "a2d2" in r["aux"] else Decimal(1))
    return float(f)
def val(v):
    r = R[v]
    return dict(v=v, p=ptype(r), q25=float(r["rate_q25"]), med=float(r["rate_median"]), q75=float(r["rate_q75"]), cv=float(r["cv_ratio_median"]), fire=int(r["both_fire"]))
V = dict((v, val(v)) for v in cube + prot)
V["WT"] = dict(v="WT", p=1.0, q25=1.0, med=1.0, q75=1.0, cv=1.0, fire=130)
if abs(V["KI_18w_mid"]["p"] - 0.672) > 1e-12:
    sys.exit("STOP: P-type factors")
# ---- values with synaptic input (130 parameter sets), gene versions only
SY = {}
_l = open(syn_tsv).read().rstrip("\n").split("\n"); _h = _l[0].split("\t")
for _r in _l[1:]:
    d = dict(zip(_h, _r.split("\t")))
    if d["variant"].startswith("P8_") or d["variant"] == "KI_5_protein_Kcnma1_Cacna1g_Cacna1a_0p85":
        continue
    if d["sets"] != "130":
        sys.exit("STOP: synaptic row not over 130 sets: " + d["variant"])
    SY[d["variant"]] = dict(v=d["variant"], p=float(d["Cacna1a"]), q25=float(d["rate_q25"]), med=float(d["rate_median"]), q75=float(d["rate_q75"]), cv=float(d["cv_ratio_median"]), fire=int(d["both_fire"]))
if len(SY) != 14 or any(v not in V for v in SY):
    sys.exit("STOP: synaptic versions not as expected: %d" % len(SY))
if any(abs(SY[v]["p"] - V[v]["p"]) > 1e-9 for v in SY):
    sys.exit("STOP: Cav2.1 factor differs between the two tables")
# ---- the 16 pairs of gene combinations that differ only in Cav2.1
key = lambda r: tuple(r[g] != "1" for g in OTHER)
at_wt = {key(R[v]): v for v in cube if R[v]["Cacna1a"] == "1"}
at_ki = {key(R[v]): v for v in cube if R[v]["Cacna1a"] != "1"}
at_wt[(False, False, False, False)] = "WT"
if len(at_wt) != 16 or len(at_ki) != 16 or set(at_wt) != set(at_ki):
    sys.exit("STOP: the 16 pairs were not found")
pairs = [(at_wt[k], at_ki[k]) for k in sorted(at_wt)]
if any(V[b]["med"] <= V[a]["med"] for a, b in pairs):
    sys.exit("STOP: a pair in which the Cav2.1 reduction does not raise the median firing rate ratio; the legend text would be wrong")
# ---- colours (one blue ramp for the three classes of the Cav2.1 conductance)
SURF, INK, INK2, MUTED, GRID, AXIS, BAND = "#ffffff", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#f0efec"
RAMP = {"at WT (1)": "#86b6ef", "at its transcript ratio (0.672)": "#104281"}
def klass(p):
    if p == 1.0:
        return "at WT (1)"
    if abs(p - 0.672) < 1e-12:
        return "at its transcript ratio (0.672)"
    sys.exit("STOP: a version with a Cav2.1 factor other than 1 or 0.672")
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 7.5, "axes.linewidth": 0.6, "axes.edgecolor": AXIS, "xtick.color": MUTED, "ytick.color": MUTED,
                     "xtick.labelcolor": INK2, "ytick.labelcolor": INK2, "axes.labelcolor": INK2, "text.color": INK, "pdf.fonttype": 42, "svg.fonttype": "none"})
fig = plt.figure(figsize=(5.0, 6.7), dpi=390, facecolor=SURF)
axc = fig.add_axes([0.118, 0.575, 0.285, 0.365])     # A: the 16 pairs
axd = fig.add_axes([0.118, 0.075, 0.467, 0.365])     # B: CV ratio against firing rate ratio
for ax in (axc, axd):
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(True, color=GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(length=2.5, width=0.5, labelsize=7)
# ---- panel A: every combination of the other four conductances, with Cav2.1 unchanged (left) and reduced (right)
axc.axhspan(LO, HI, color=BAND, zorder=0.5, linewidth=0)
KVI = OTHER.index("Kcnd3")
kred = dict(((at_wt[k], at_ki[k]), k[KVI]) for k in sorted(at_wt))
if sum(kred.values()) != 8:
    sys.exit("STOP: pairs with Kv4.3 reduced")
KRED, KUNCH = "#c2410c", "#0f8a6d"
hp = ("KI_1c_restore_Cacna1a", "KI_18w_mid")
if hp not in pairs:
    sys.exit("STOP: the named pair is not one of the 16 pairs")
def also(ab):
    g = [CH[x] for x in OTHER if R[ab[1]][x] != "1"]
    return "Cav2.1 alone" if not g else ("all four others" if len(g) == 4 else " + ".join(g))
for a_, b_ in pairs:
    if (a_, b_) == hp:
        continue
    _c = KRED if kred[(a_, b_)] else KUNCH
    axc.plot([0, 1], [V[a_]["med"], V[b_]["med"]], color=_c, linewidth=0.7, linestyle="-", zorder=2)
    axc.plot([0, 1], [V[a_]["med"], V[b_]["med"]], linestyle="none", marker="o", markersize=3.0, color=_c, markeredgecolor=SURF, markeredgewidth=0.4, zorder=3)
axc.plot([0, 1], [V[hp[0]]["med"], V[hp[1]]["med"]], color=INK2, linewidth=1.3, solid_capstyle="round", zorder=4)
for xx, v in zip((0, 1), hp):
    x = V[v]
    axc.plot([xx, xx], [x["q25"], x["q75"]], color=INK2, linewidth=1.0, solid_capstyle="round", zorder=4)
    axc.plot([xx], [x["med"]], marker="o", markersize=5.0, color=INK2, markeredgecolor=SURF, markeredgewidth=0.9, zorder=6)
    axc.plot([xx + 0.06], [SY[v]["med"]], marker="o", markersize=4.4, color=SURF, markeredgecolor=INK2, markeredgewidth=1.1, zorder=7)
axc.plot([0], [1.0], linestyle="none", marker="s", markersize=4.0, color=SURF, markeredgecolor=INK, markeredgewidth=0.9, zorder=6.5)
# every line is named at its right end by what else is reduced; the names are spread so that they do not touch
lab = sorted(((V[b_]["med"], also((a_, b_)), (a_, b_)) for a_, b_ in pairs))
if len(set(x[1] for x in lab)) != 16:
    sys.exit("STOP: the 16 names are not distinct")
STEP = 0.098
ypos = [x[0] for x in lab]
for _ in range(200):
    for k in range(1, len(ypos)):
        if ypos[k] - ypos[k - 1] < STEP:
            m_ = (ypos[k] + ypos[k - 1]) / 2; ypos[k - 1], ypos[k] = m_ - STEP / 2, m_ + STEP / 2
for (yv, name, ab), yl in zip(lab, ypos):
    axc.plot([1.03, 1.10, 1.16], [yv, yl, yl], color=(INK2 if ab == hp else (KRED if kred[ab] else KUNCH)), linewidth=0.5, clip_on=False, zorder=1)
    axc.text(1.19, yl, name, fontsize=6.3, color=INK, ha="left", va="center", clip_on=False, fontweight=("bold" if ab == hp else "normal"))
axc.text(1.19, ypos[-1] + 0.15, "also reduced:", fontsize=6.3, color=MUTED, ha="left", va="center", clip_on=False)
axc.text(1.1, LO + 0.02, "range recorded in\nKI Purkinje cells", fontsize=6, color=MUTED, ha="right", va="bottom", linespacing=1.05)
axc.plot([-0.045, -0.1], [1.0, 0.93], color=MUTED, linewidth=0.4, zorder=1)
axc.text(-0.11, 0.92, "WT", fontsize=6.3, color=INK, ha="right", va="center")
axc.set_xlim(-0.22, 1.12); axc.set_ylim(0.6, 2.6)
axc.set_xticks([0, 1]); axc.set_xticklabels(["Cav2.1\nunchanged", "Cav2.1 reduced\n(transcript\nratio 0.672)"])
axc.set_yticks([0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5]); axc.set_yticklabels(["0.75", "1.00", "1.25", "1.50", "1.75", "2.00", "2.25", "2.50"])
axc.set_ylabel("Firing rate ratio (KI / WT)")
axc.grid(axis="x", visible=False)
hc = [Line2D([0], [0], color=KRED, linewidth=0.7, marker="o", markersize=3.0, markeredgecolor=SURF, markeredgewidth=0.4, label="Kv4.3 reduced (8 lines)"),
      Line2D([0], [0], color=KUNCH, linewidth=0.7, marker="o", markersize=3.0, markeredgecolor=SURF, markeredgewidth=0.4, label="Kv4.3 unchanged (8)"),
      Line2D([0], [0], color=INK2, linewidth=1.3, marker="o", markersize=5.0, markeredgecolor=SURF, markeredgewidth=0.9, label="all four others reduced\n(right dot: all five)"),
      Line2D([0], [0], color=INK2, linestyle="none", marker="|", markersize=7.0, markeredgewidth=1.0, label="vertical bars: quartiles"),
      Line2D([0], [0], linestyle="none", marker="o", markersize=4.4, color=SURF, markeredgecolor=INK2, markeredgewidth=1.1, label="with synaptic input")]
fig.canvas.draw()
_xr = max(t.get_window_extent(fig.canvas.get_renderer()).x1 for t in axc.texts) / (fig.get_size_inches()[0] * fig.dpi)
lgc = fig.legend(handles=hc, loc="upper left", bbox_to_anchor=(_xr + 0.045, 0.94), bbox_transform=fig.transFigure, frameon=False, fontsize=6, handlelength=2.0, labelspacing=0.55, borderaxespad=0, labelcolor=INK2)
axc.set_title("One line per combination\nof the other four\nconductances", fontsize=7, color=INK2, loc="left", pad=6)
# ---- panel B (was C): CV ratio against firing rate ratio
axd.axvspan(LO, HI, color=BAND, zorder=0.5, linewidth=0)
order = ["at WT (1)", "at its transcript ratio (0.672)"]
SHOWN = {"at WT (1)": "unchanged (1)", "at its transcript ratio (0.672)": "at its transcript ratio (0.672)"}   # legend texts; the values file keeps the class names "at WT (1)" and so on
for k in order:
    pts = [x for x in V.values() if klass(x["p"]) == k and x["v"] != "WT"]
    axd.plot([x["med"] for x in pts], [x["cv"] for x in pts], linestyle="none", marker="o", markersize=4.8, color=RAMP[k], markeredgecolor=SURF, markeredgewidth=0.8, zorder=4)
for k in order:
    pts = [x for x in SY.values() if klass(x["p"]) == k]
    axd.plot([x["med"] for x in pts], [x["cv"] for x in pts], linestyle="none", marker="o", markersize=4.6, color="none", markeredgecolor=RAMP[k], markeredgewidth=1.4, zorder=5)
for name, rr, cr in REC:
    axd.plot([rr], [cr], linestyle="none", marker="D", markersize=4.6, color=INK, markeredgecolor=SURF, markeredgewidth=0.8, zorder=6)
    axd.annotate(name, (rr, cr), xytext=(rr + 0.06, cr), fontsize=7, color=INK, ha="left", va="center")
top = max((x for x in V.values() if x["v"] != "WT"), key=lambda x: x["cv"])
lab = " + ".join(CH[g] for g in GENES if R[top["v"]][g] != "1")
for x, txt, tx, ty, ha in ((top, "highest CV ratio\nwithout synaptic input:\n" + lab, 1.70, 2.08, "left"), (V["KI_18w_mid"], "all five transcript ratios", 1.86, 0.9, "center")):
    axd.annotate(txt, (x["med"], x["cv"]), xytext=(tx, ty), fontsize=7, color=INK, ha=ha, va="center", linespacing=1.1,
                 arrowprops=dict(arrowstyle="-", color=MUTED, linewidth=0.5, shrinkA=1, shrinkB=3), zorder=7)
axd.text((LO + HI) / 2, 2.32, "range recorded in\nKI Purkinje cells", fontsize=6.5, color=MUTED, ha="center", va="center", linespacing=1.1)
axd.set_xlim(0.6, 2.6); axd.set_ylim(0.7, 3.05)
axd.set_xticks([0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5]); axd.set_xticklabels(["0.75", "1.00", "1.25", "1.50", "1.75", "2.00", "2.25", "2.50"])
axd.set_yticks([1.0, 1.5, 2.0, 2.5, 3.0]); axd.set_yticklabels(["1.0", "1.5", "2.0", "2.5", "3.0"])
axd.set_xlabel("Firing rate ratio (KI / WT)"); axd.set_ylabel("CV ratio (KI / WT)")
hd = [Line2D([0], [0], linestyle="none", marker="o", markersize=4.8, color=RAMP[k], markeredgecolor=SURF, markeredgewidth=0.8, label=SHOWN[k]) for k in order]
hd.append(Line2D([0], [0], linestyle="none", marker="o", markersize=4.4, color="none", markeredgecolor=INK2, markeredgewidth=1.1, label="open: with synaptic input\n(14 versions)"))
hd.append(Line2D([0], [0], linestyle="none", marker="D", markersize=4.6, color=INK, markeredgecolor=SURF, markeredgewidth=0.8, label="recorded in KI Purkinje cells\n(Bushart et al. 2021)"))
lg = axd.legend(handles=hd, loc="upper left", bbox_to_anchor=(1.05, 1.0), frameon=False, fontsize=6.5, handletextpad=0.3, labelspacing=0.45, borderaxespad=0, labelcolor=INK2,
                title="One point per version;\nCav2.1 conductance:", title_fontsize=6.5, alignment="left")
lg.get_title().set_color(INK2)
for ax, letter, dx in ((axc, "A", -0.35), (axd, "B", -0.215)):
    ax.text(dx, 1.09, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", color=INK, va="bottom", ha="left")
fig.text(0.985, 0.012, "Filled: runs without synaptic input; open: with synaptic input.\nMedians over the 130 parameter sets (sets that fire in both runs).", fontsize=6, color=MUTED, ha="right", va="bottom")
import os
fig.canvas.draw()
_r = fig.canvas.get_renderer()
_c, _d = lgc.get_window_extent(_r), lg.get_window_extent(_r)
_gap = (_c.y0 - _d.y1) / fig.dpi * 2.54
_w, _h = fig.get_size_inches() * fig.dpi
_out = [t.get_text()[:40] for t in fig.findobj(matplotlib.text.Text) if t.get_text().strip() and t.get_visible() and ((e := t.get_window_extent(_r)).x0 < -0.5 or e.x1 > _w + 0.5 or e.y0 < -0.5 or e.y1 > _h + 0.5)]
print("gap between the key of A and the key of B: %.3f cm | right edge of the key of C: %.3f cm of %.3f | texts outside the page: %r" % (_gap, _d.x1 / fig.dpi * 2.54, _w / fig.dpi * 2.54, _out))
if _gap < 0.2 or _d.x1 > _w or _out:
    sys.exit("STOP: the legends meet or a text leaves the page")
for ext in (".png", ".pdf", ".tsv"):
    if os.path.exists(prefix + ext):
        sys.exit("STOP: output exists: " + prefix + ext)
fig.savefig(prefix + ".png", dpi=390, facecolor=SURF)
fig.savefig(prefix + ".pdf", facecolor=SURF)
rows = ["version\tptype_factor\tclass\trate_q25\trate_median\trate_q75\tcv_ratio_median\tsets_firing_in_both\tin_pair_of_panel_A"]
inpair = set(a for a, b in pairs) | set(b for a, b in pairs)
for v in ["WT"] + cube + prot:
    x = V[v]
    rows.append("%s\t%.3f\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%s" % (v, x["p"], klass(x["p"]), x["q25"], x["med"], x["q75"], x["cv"], x["fire"], "yes" if v in inpair else "no"))
rows.append("# with synaptic input")
for v in sorted(SY):
    x = SY[v]
    rows.append("%s\t%.3f\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%d\tsynaptic_input" % (v, x["p"], klass(x["p"]), x["q25"], x["med"], x["q75"], x["cv"], x["fire"]))
print("synaptic versions plotted in B:", len(SY), "| highest median CV ratio with synaptic input:", max(SY.values(), key=lambda x: x["cv"])["v"], max(x["cv"] for x in SY.values()))
open(prefix + ".tsv", "x").write("\n".join(rows) + "\n")
print("versions plotted:", len(V) - 1, "and WT | pairs:", len(pairs), "| recorded:", ["%s %.4f %.4f" % r for r in REC])
print("classes:", dict((k, sum(1 for x in V.values() if klass(x["p"]) == k and x["v"] != "WT")) for k in order))
print("highest median CV ratio among the plotted versions:", top["v"], top["cv"], "at firing rate ratio", top["med"])
print("pairs: smallest and largest quotient of the two medians: %.3f %.3f" % (min(V[b]["med"] / V[a]["med"] for a, b in pairs), max(V[b]["med"] / V[a]["med"] for a, b in pairs)))
