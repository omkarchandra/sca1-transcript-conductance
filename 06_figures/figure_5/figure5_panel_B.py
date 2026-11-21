"""Figure 5, panel B: the firing rate ratio (KI / WT) of every parameter set, version by version, without and with
synaptic input.

Usage: python3 figure5_panel_B.py PER_SET_VALUES_TXT REFERENCE_VALUES_MD OUT_PREFIX VERSION=label [VERSION=label ...]

PER_SET_VALUES_TXT holds the per-set values in blocks. Each block starts with a line "EXTRACT v1 | label <condition> |
versions <n>" and ends with "END EXTRACT | label <condition>"; inside, every version has the lines VERSION, SETS, RATE,
CV, AHP, WTHZ, KISTATE, WTSTATE and SUMS. Every block is checked against its own sums before anything is drawn; text
outside the blocks is ignored, and a line that was broken in two is joined again.
Each VERSION named on the command line is drawn, in that order, with one strip of points per condition (the label of
each block). Set 0, the published parameter set, is left out of the strips and of the medians and is drawn apart.
A block whose label ends in _area (the membrane-area versions) is joined to the block of the same label without that
ending; the two blocks must not hold the same version.
The recorded range is taken from the firing-frequency row of REFERENCE_VALUES_MD (Bushart et al. 2021).
Writes OUT_PREFIX.png, OUT_PREFIX.pdf and OUT_PREFIX.tsv (the medians and quartiles drawn) by exclusive create.
The count of sets inside the recorded range in the TSV is taken from the three-decimal values of the input and can
differ by a set at the bound from the count of the scored table, which uses the unrounded values.
"""
import os, re, sys
from decimal import Decimal

KEYS = ("EXTRACT", "VERSION", "SETS", "RATE", "CV", "AHP", "WTHZ", "KISTATE", "WTSTATE", "SUMS", "END")


def parse_extracts(text):
    """Returns {label: {version: dict}}; raises ValueError on any inconsistency. A line that was broken in two is joined
    again; blanks inside a line carry no information (no value, name or label holds one) and are ignored."""
    lines, inside = [], False
    for raw in text.split("\n"):
        s = raw.strip()
        if not s:
            continue
        m = re.match(r"(EXTRACT|VERSION|SETS|RATE|CV|AHP|WTHZ|KISTATE|WTSTATE|SUMS|END)\b", s)
        key = m.group(1) if m else None
        if key == "EXTRACT" and not re.match(r"EXTRACT\s*v1\b", s):
            key = None
        if not inside:
            if key == "EXTRACT":
                inside = True
                lines.append((key, s[len(key):]))
            continue
        if key:
            lines.append((key, s[len(key):]))
        else:
            lines[-1] = (lines[-1][0], lines[-1][1] + s)          # a line that was broken in two
        if key == "END":
            inside = False
    out, label, first, cur, want = {}, None, None, None, 0
    def close(b):
        n = len(b["sets"])
        for k in ("rate", "cv", "ahp", "wthz"):
            if len(b[k]) != n:
                raise ValueError("%s %s: %d values of %s for %d sets" % (label, b["v"], len(b[k]), k, n))
            num = [Decimal(x) for x in b[k] if x != "x"]
            if (str(sum(num, Decimal(0))), len(num)) != b["sums"][k]:
                raise ValueError("%s %s: the values of %s do not give the printed sum" % (label, b["v"], k))
        if len(b["kistate"]) != n or len(b["wtstate"]) != n or n != b["rows"] or len(set(b["sets"])) != n:
            raise ValueError("%s %s: states, rows or set numbers do not fit" % (label, b["v"]))
    for key, rest in lines:
        rest = re.sub(r"\s+", "", rest)
        if key == "EXTRACT":
            m = re.fullmatch(r"v1\|label(\w+)\|versions(\d+)", rest)
            if not m or m.group(1) in out or label is not None:
                raise ValueError("bad or repeated EXTRACT line: " + rest[:80])
            label, first, cur, want = m.group(1), None, None, int(m.group(2))
            out[label] = {}
            continue
        if label is None:
            raise ValueError("a %s line outside a block" % key)
        if key == "VERSION":
            m = re.fullmatch(r"(.+?)\|rows(\d+)\|tablesha256([0-9a-f]{64})", rest)
            if not m or m.group(1) in out[label]:
                raise ValueError("bad or repeated VERSION line: " + rest[:80])
            cur = dict(v=m.group(1), rows=int(m.group(2)), sha=m.group(3))
            out[label][cur["v"]] = cur
            if first is None:
                first = cur
        elif key in ("SETS", "RATE", "CV", "AHP", "WTHZ", "KISTATE", "WTSTATE"):
            k = key.lower()
            if cur is None or k in cur:
                raise ValueError("a %s line without its VERSION line, or twice" % key)
            if rest == "=first":
                if cur is first or k not in first:
                    raise ValueError("'= first' without a first version in " + cur["v"])
                cur[k] = first[k]
            elif key == "SETS":
                if re.fullmatch(r"\d+(,\d+)*", rest) is None:
                    raise ValueError("bad SETS line of " + cur["v"])
                cur[k] = [int(x) for x in rest.split(",")]
            elif key in ("KISTATE", "WTSTATE"):
                if re.fullmatch(r"[FPBR]+", rest) is None:
                    raise ValueError("bad state line of " + cur["v"])
                cur[k] = rest
            else:
                vals = rest.split(",")
                if any(re.fullmatch(r"x|-?\d+\.\d+", x) is None for x in vals):
                    raise ValueError("bad value in the %s line of %s" % (key, cur["v"]))
                cur[k] = vals
        elif key == "SUMS":
            m = re.fullmatch(r"rate(-?[\d.]+)n(\d+)\|cv(-?[\d.]+)n(\d+)\|ahp(-?[\d.]+)n(\d+)\|wthz(-?[\d.]+)n(\d+)", rest)
            if not m or cur is None or any(k not in cur for k in ("sets", "rate", "cv", "ahp", "wthz", "kistate", "wtstate")):
                raise ValueError("bad SUMS line, or lines missing before it")
            g = m.groups()
            cur["sums"] = dict((k, (g[2 * j], int(g[2 * j + 1]))) for j, k in enumerate(("rate", "cv", "ahp", "wthz")))
            close(cur)
        elif key == "END":
            if rest != "EXTRACT|label" + label or len(out[label]) != want or any("sums" not in b for b in out[label].values()):
                raise ValueError("block %s is not complete" % label)
            label = None
    if label is not None:
        raise ValueError("the last block has no END line")
    if not out:
        raise ValueError("no EXTRACT block found")
    return out


def q(xs, p):
    xs = sorted(xs)
    k = (len(xs) - 1) * p
    lo, hi = int(k // 1), int(-(-k // 1))
    return xs[lo] + (xs[hi] - xs[lo]) * (k - lo)


def main(extract, note, prefix, pairs):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    for ext in (".png", ".pdf", ".tsv"):
        if os.path.exists(prefix + ext):
            sys.exit("STOP: output exists: " + prefix + ext)
    try:
        E = parse_extracts(open(extract, encoding="utf-8").read())
    except ValueError as err:
        sys.exit("STOP: " + str(err))
    # recorded range from the reference note
    def cells(l):
        return [c.strip() for c in l.strip().strip("|").split("|")]
    num = lambda s: float(s.split("±")[0].replace("−", "-").replace(",", "").strip())
    rows = [cells(l) for l in open(note, encoding="utf-8").read().split("\n") if cells(l)[:1] == ["Firing frequency"] and len(cells(l)) == 7]
    if len(rows) != 1:
        sys.exit("STOP: firing-frequency row of the reference note")
    LO, HI = sorted([num(rows[0][5]) / num(rows[0][4]), num(rows[0][2]) / num(rows[0][1])])
    if "%.2f-%.2f" % (LO, HI) != "0.86-1.33":
        sys.exit("STOP: recorded range")
    for l in [x for x in E if x.endswith("_area")]:
        base = l[:-5]
        if base in E and set(E[base]) & set(E[l]):
            sys.exit("STOP: a version stands in both blocks: " + l)
        E.setdefault(base, {}).update(E.pop(l))
    labels = [x for x in ("no_input", "synaptic_input") if x in E] + [x for x in E if x not in ("no_input", "synaptic_input")]
    if len(labels) > 2:
        sys.exit("STOP: more than two conditions")
    for v, _ in pairs:
        if not any(v in E[l] for l in labels):
            sys.exit("STOP: version in no block: " + v)
    SURF, INK, INK2, MUTED, GRID, AXIS, BAND = "#ffffff", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#f0efec"
    COL = dict(zip(labels, ["#2a78d6", "#eb6834"]))      # blue and orange
    NAME = {"no_input": "without synaptic input", "synaptic_input": "with synaptic input"}
    plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 7.5, "axes.linewidth": 0.6, "axes.edgecolor": AXIS, "xtick.color": MUTED, "ytick.color": MUTED,
                         "xtick.labelcolor": INK2, "ytick.labelcolor": INK2, "axes.labelcolor": INK2, "text.color": INK, "pdf.fonttype": 42})
    fig, ax = plt.subplots(figsize=(5.0, 3.45), dpi=390, facecolor=SURF)
    fig.subplots_adjust(left=0.098, right=0.985, top=0.94, bottom=0.385)
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(True, axis="y", color=GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(length=2.5, width=0.5, labelsize=7)
    ax.axhspan(LO, HI, color=BAND, zorder=0.5, linewidth=0)
    rng = np.random.RandomState(20261114)
    out_rows = ["version\tcondition\tsets_with_a_ratio\tsets_without\tq25\tmedian\tq75\tin_recorded_range\tpublished_set"]
    top = 0.0
    off = {labels[0]: -0.19, labels[-1]: 0.19} if len(labels) == 2 else {labels[0]: 0.0}
    for i, (v, _) in enumerate(pairs):
        for l in labels:
            if v not in E[l]:
                continue
            b = E[l][v]
            vals = [(k, float(x)) for k, x in zip(b["sets"], b["rate"]) if x != "x" and k != 0]
            none = sum(1 for k, x in zip(b["sets"], b["rate"]) if x == "x" and k != 0)
            pub = [float(x) for k, x in zip(b["sets"], b["rate"]) if k == 0 and x != "x"]
            ys = [y for _, y in vals]
            xs = i + off[l] + rng.uniform(-0.13, 0.13, size=len(ys))
            ax.plot(xs, ys, linestyle="none", marker="o", markersize=2.3, color=COL[l], alpha=0.5, markeredgewidth=0, zorder=3)
            if ys:
                q1, md, q3 = q(ys, 0.25), q(ys, 0.5), q(ys, 0.75)
                ax.plot([i + off[l], i + off[l]], [q1, q3], color=INK, linewidth=0.7, solid_capstyle="butt", zorder=5)
                ax.plot([i + off[l] - 0.13, i + off[l] + 0.13], [md, md], color=INK, linewidth=1.2, solid_capstyle="butt", zorder=6)
                top = max(top, max(ys))
                out_rows.append("%s\t%s\t%d\t%d\t%.3f\t%.3f\t%.3f\t%d\t%s" % (v, l, len(ys), none, q1, md, q3, sum(1 for y in ys if LO <= y <= HI), ("%.3f" % pub[0]) if pub else ""))
            if pub:
                ax.plot([i + off[l]], [pub[0]], linestyle="none", marker="D", markersize=3.6, color=SURF, markeredgecolor=INK, markeredgewidth=0.8, zorder=7)
                top = max(top, pub[0])
            ax.text(i + off[l], 0.02, "%d" % len(ys), transform=ax.get_xaxis_transform(), fontsize=6, color=MUTED, ha="center", va="bottom")
    ax.set_xlim(-0.6, len(pairs) - 0.4)
    ymax = max(2.5, 0.25 * (int(top / 0.25) + 1))
    ax.set_ylim(0.35, ymax)
    ax.set_yticks([0.5 + 0.25 * k for k in range(int((ymax - 0.5) / 0.25) + 1) if (0.5 + 0.25 * k) <= ymax and k % 2 == 0])
    ax.set_xticks(range(len(pairs)))
    ax.set_xticklabels([n.replace("|", "\n") for _, n in pairs])
    ax.set_ylabel("Firing rate ratio (KI / WT)")
    h = [Line2D([0], [0], linestyle="none", marker="o", markersize=3.2, color=COL[l], alpha=0.8, markeredgewidth=0, label="one parameter set, " + NAME.get(l, l)) for l in labels]
    h.append(Line2D([0], [0], color=INK, linewidth=1.2, label="median; bar, quartiles"))
    h.append(Patch(facecolor=BAND, edgecolor="none", label="firing-rate range recorded in KI Purkinje cells"))
    if any(k == 0 for l in labels for b in E[l].values() for k in b["sets"]):
        h.append(Line2D([0], [0], linestyle="none", marker="D", markersize=3.6, color=SURF, markeredgecolor=INK, markeredgewidth=0.8, label="published parameter set"))
    lg = ax.legend(handles=h, loc="upper left", bbox_to_anchor=(-0.02, -0.27), frameon=False, fontsize=6.5, ncol=2, handletextpad=0.5, columnspacing=1.6, labelspacing=0.5, borderaxespad=0, labelcolor=INK2)
    ax.text(-0.02, -0.56, "Small numbers above the axis: parameter sets with a firing rate ratio (both runs fire).", transform=ax.transAxes, fontsize=6.5, color=INK2, ha="left", va="top")
    ax.text(-0.105, 1.03, "B", transform=ax.transAxes, fontsize=10, fontweight="bold", color=INK, va="bottom", ha="left")
    fig.savefig(prefix + ".png", dpi=390, facecolor=SURF)
    fig.savefig(prefix + ".pdf", facecolor=SURF)
    open(prefix + ".tsv", "x").write("\n".join(out_rows) + "\n")
    print("conditions:", labels, "| versions drawn:", [v for v, _ in pairs])
    for r in out_rows:
        print(r)


if __name__ == "__main__":
    a = sys.argv[1:]
    if len(a) < 4 or any("=" not in x for x in a[3:]):
        sys.exit(__doc__)
    main(a[0], a[1], a[2], [tuple(x.split("=", 1)) for x in a[3:]])
