"""Paper 1 Figure 2: characterisation of LDD MRI features.

A  baseline prevalence by age group (Wilson 95% CI), all.emergence.data.20201215.csv
B  feature-by-level correlation network at baseline (phi coefficients, same file)
C  cumulative risk of onset (SISE); curves digitised from Figures/Figure 2.tif
   because the SISE output object is not in DDD_Paper -> data/sise_cumrisk_digitized.csv
D  dynamic Bayesian networks (edges as reported in the text / 2021 figure)
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle, Patch
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import matplotlib.colors as mc
import style as S
import derive as Dv
from common import DATA, P1

S.setup()
ORDER = ["SS", "DB", "DN", "DH", "HIZ", "MC", "SN"]


def shade(c, f):
    """Blend colour c toward white by fraction f (0 = c)."""
    r = np.array(mc.to_rgb(c))
    return tuple(r + (1 - r) * f)


def panel_a(ax, e):
    pv = Dv.prevalence_by_age(e)
    pv["F"] = pv.feature.str.upper()
    ng = len(Dv.AGELAB)
    w = 0.8 / ng
    fades = np.linspace(0.72, 0.0, ng)
    for i, f in enumerate(ORDER):
        d = pv[pv.F == f].set_index("agegroup").reindex(Dv.AGELAB)
        xs = i + (np.arange(ng) - (ng - 1) / 2) * w
        cols = [shade(S.FEAT[f], fa) for fa in fades]
        ax.bar(xs, d.p, w * 0.92, color=cols, edgecolor="none")
        ax.errorbar(xs, d.p, yerr=[d.p - d.lo, d.hi - d.p], fmt="none", ecolor="#3A3A3A", elinewidth=0.5,
                    capsize=1.2, capthick=0.5)
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels(ORDER, fontweight="bold")
    for t, f in zip(ax.get_xticklabels(), ORDER):
        t.set_color(S.FEAT[f])
    ax.tick_params(axis="x", length=0)
    ax.set_ylim(0, 1.04)
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.set_yticklabels([f"{int(v*100)}" for v in np.arange(0, 1.01, 0.2)])
    ax.set_ylabel("Prevalence at first MRI (%)")
    ax.yaxis.grid(True); ax.set_axisbelow(True)
    n = pv[pv.feature == "ss"].set_index("agegroup").reindex(Dv.AGELAB).n
    h = [Patch(color=shade("#555555", fa), label=f"{a} (n = {int(k):,})") for a, fa, k in zip(Dv.AGELAB, fades, n)]
    ax.legend(handles=h, title="Age at first MRI (years)", ncol=1, loc="upper right",
              bbox_to_anchor=(1.0, 1.02), handlelength=1.0, handleheight=0.8, labelspacing=0.25)
    ax.set_xlim(-0.5, len(ORDER) - 0.5 + 1.15)


def panel_b(ax, e):
    F = ["ss", "db", "dn", "dh", "hiz", "mc", "sn"]
    cols = [f"{f}.L{l}.t1" for f in F for l in range(1, 6)]
    x = e[cols].dropna()
    C = x.corr().values
    prev = x.mean().values
    names = [f"{f.upper()}{l}" for f in F for l in range(1, 6)]
    n = len(names)
    gap = 0.35  # extra angular gap between features (in node units)
    pos_idx = np.array([i + gap * (i // 5) for i in range(n)])
    ang = np.pi / 2 - 2 * np.pi * pos_idx / (pos_idx.max() + 1 + gap)
    P = np.c_[np.cos(ang), np.sin(ang)]
    thr = 0.15
    pairs = [(i, j, C[i, j]) for i in range(n) for j in range(i + 1, n) if C[i, j] >= thr]
    pairs.sort(key=lambda t: t[2])
    for i, j, r in pairs:
        same = names[i][:-1] == names[j][:-1]
        ctrl = (P[i] + P[j]) / 2 * (0.75 if same else 0.15)
        path = Path([P[i], ctrl, P[j]], [Path.MOVETO, Path.CURVE3, Path.CURVE3])
        a = min(1, 0.25 + (r - thr) * 1.6)
        ax.add_patch(PathPatch(path, fc="none", ec=(0.18, 0.30, 0.45, a), lw=0.3 + 6.0 * (r - thr)))
    for k in range(n):
        f = names[k][:-1]
        c = S.FEAT[f]
        rad = 0.035 + 0.075 * np.sqrt(prev[k])
        ax.add_patch(Circle(P[k], rad, fc=c, ec="white", lw=0.6, zorder=3))
        lp = P[k] * (1.0 + rad + 0.07)
        ax.text(*lp, names[k], ha="center", va="center", fontsize=5.6, color=c, zorder=4,
                rotation=0)
    ax.set_xlim(-1.32, 1.32); ax.set_ylim(-1.32, 1.32); ax.set_aspect("equal"); ax.axis("off")
    # legend for edge width
    for k, r in enumerate([0.2, 0.4, 0.6]):
        y = -1.22 + 0.0 * k
        x0 = 0.55 + 0.27 * k
        ax.plot([x0, x0 + 0.17], [y, y], color=(0.18, 0.30, 0.45, min(1, 0.25 + (r - thr) * 1.6)),
                lw=0.3 + 6.0 * (r - thr), solid_capstyle="butt")
        ax.text(x0 + 0.085, y - 0.07, f"{r:.1f}", ha="center", va="top", fontsize=5.5)
    ax.text(0.53, -1.22, "r ", ha="right", va="center", fontsize=6, style="italic")
    ax.text(-1.3, -1.25, "Node size: prevalence\nEdges: r ≥ 0.15", fontsize=5.5, color=S.GREY, va="bottom")


def panel_c(ax):
    d = pd.read_csv(DATA + "sise_cumrisk_digitized.csv")
    ends = []
    for f in ORDER:
        g = d[d.feature == f]
        ax.plot(g.age, g.cumrisk, color=S.FEAT[f], lw=1.3)
        ends.append([f, g.age.iloc[-1], g.cumrisk.iloc[-1]])
    ends.sort(key=lambda t: -t[2])
    lab_y = []
    for f, x, y in ends:  # push labels apart vertically
        yy = min(y, lab_y[-1] - 0.06) if lab_y else y
        lab_y.append(yy)
        ax.plot([x + 0.5, 90.5], [y, yy], color=S.FEAT[f], lw=0.4)
        ax.text(91.5, yy, f, color=S.FEAT[f], fontsize=6.5, fontweight="bold", va="center")
    ax.axhline(0.5, color=S.GREY, lw=0.5, ls=(0, (2, 2)), zorder=0)
    ax.set_xlim(5, 97); ax.set_ylim(0, 1.02)
    ax.set_xticks(range(10, 91, 10))
    ax.set_xlabel("Age (years)"); ax.set_ylabel("Cumulative risk of onset")
    ax.yaxis.grid(True); ax.set_axisbelow(True)


def dbn(ax, nodes, colours, edges, title):
    """Two-slice DBN: grey = autoregressive (t -> t+1), black = cross edges."""
    n = len(nodes)
    ys = np.linspace(1, 0, n)
    xl, xr = 0.18, 0.82
    for k, (nd, c) in enumerate(zip(nodes, colours)):
        for x in (xl, xr):
            ax.add_patch(Circle((x, ys[k]), 0.052, fc=c, ec="white", lw=0.8, zorder=3, transform=ax.transData))
            ax.text(x, ys[k], nd, ha="center", va="center", fontsize=6, color="white", fontweight="bold", zorder=4)
        ax.add_patch(FancyArrowPatch((xl + 0.055, ys[k]), (xr - 0.06, ys[k]), arrowstyle="-|>", mutation_scale=6,
                                     color="#C4C8CE", lw=0.8, zorder=1))
    idx = {nd: k for k, nd in enumerate(nodes)}
    for a, b in edges:
        ax.add_patch(FancyArrowPatch((xl + 0.05, ys[idx[a]]), (xr - 0.058, ys[idx[b]]), arrowstyle="-|>",
                                     mutation_scale=7, color=S.INK, lw=0.9, zorder=2,
                                     connectionstyle="arc3,rad=0.0"))
    ax.text(xl, 1.13, "Scan t", ha="center", fontsize=6.5, color=S.GREY)
    ax.text(xr, 1.13, "Scan t + 1", ha="center", fontsize=6.5, color=S.GREY)
    ax.set_title(title, pad=14, fontweight="normal")
    ax.set_xlim(0, 1); ax.set_ylim(-0.1, 1.1); ax.set_aspect("equal"); ax.axis("off")


def main():
    e = Dv.emergence()
    fig = plt.figure(figsize=(S.COL2, 178 * S.MM))
    gs = fig.add_gridspec(3, 4, height_ratios=[0.8, 1.3, 1.05], hspace=0.30, wspace=0.35,
                          left=0.07, right=0.985, top=0.975, bottom=0.02)
    axA = fig.add_subplot(gs[0, :]); panel_a(axA, e); S.panel(axA, "A", x=-0.065, y=1.0)
    axB = fig.add_subplot(gs[1, :2]); panel_b(axB, e); S.panel(axB, "B", x=-0.02, y=0.96)
    axC = fig.add_subplot(gs[1, 2:]); panel_c(axC); S.panel(axC, "C", x=-0.16, y=0.98)
    axC.set_position([0.585, axC.get_position().y0 + 0.03, 0.37, axC.get_position().height - 0.05])
    f_nodes = ["SS", "DB", "DN", "DH", "HIZ", "MC", "SN"]
    f_edges = [("SS", "DB"), ("SS", "DH"), ("SS", "HIZ"), ("SS", "MC"), ("DB", "DN"), ("DB", "HIZ"),
               ("DN", "SS"), ("MC", "DN")]
    ax1 = fig.add_subplot(gs[2, :2])
    dbn(ax1, f_nodes, [S.FEAT[f] for f in f_nodes], f_edges, "Among MRI features")
    S.panel(ax1, "D", x=-0.02, y=1.02)
    lv = ["L1", "L2", "L3", "L4", "L5"]
    lc = [shade("#B5473A", f) for f in [0.55, 0.42, 0.28, 0.14, 0.0]]
    ax2 = fig.add_subplot(gs[2, 2:])
    dbn(ax2, lv, lc, [("L2", "L1"), ("L3", "L2"), ("L4", "L3")], "Among disc levels (any feature)")
    S.save(fig, P1, "Fig2_LDD_features")


if __name__ == "__main__":
    main()
