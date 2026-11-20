"""Paper 2 GWAS figures: Manhattan + regional plot (Fig 1), fine-mapping (Fig 2),
QQ plots (Fig S1) and Manhattan plots for the other three phenotypes (Fig S2).

Manhattan input keeps every SNP with p < 0.01 plus a random 5% of the rest
(prep_sumstats.py on DDD_Paper/GWAS_sumstat/summary.*.csv); QQ input keeps the
top 20,000 p-values and every 250th thereafter.
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import style as S
from common import SRC, WORK, P2

S.setup()
OUT = P2
PH = {"MC": "MC.avg", "Core": "Core.slope", "HIZ": "HIZ.avg", "DH": "DH.avg"}
LAMBDA = {"MC": 1.009, "Core": 1.012, "HIZ": 1.005, "DH": 1.005}
GWS, SUG = 5e-8, 1e-5


def manhattan(ax, ph, label_top=None, ymax=None):
    d = pd.read_csv(WORK + f"manh_{ph}.csv.gz")
    d = d[d.CHR.between(1, 22)].sort_values(["CHR", "BP"])
    # input keeps all p < 0.01 but only 5% of the rest; taper 2 < -log10 p < 3 so density stays continuous
    lp = -np.log10(d.P.values)
    keep = np.random.default_rng(2).random(len(d)) < np.clip(0.05 * 10 ** (1.3 * (lp - 2)), 0.05, 1)
    d = d[(lp < 2) | keep]
    off, mids, start = {}, [], 0
    for c in range(1, 23):
        s = d.loc[d.CHR == c, "BP"]
        off[c] = start
        mids.append(start + s.max() / 2)
        start += s.max() + 1.5e7
    x = d.BP.values + d.CHR.map(off).values
    y = -np.log10(d.P.values)
    c1, c2 = S.FEAT[ph], "#A9B4C2"
    col = np.where(d.CHR.values % 2 == 1, c1, c2)
    ax.scatter(x, y, c=col, s=1.2, lw=0, rasterized=True)
    hit = d.P.values < GWS
    if hit.any():
        ax.scatter(x[hit], y[hit], s=5, c="#D7263D", lw=0, zorder=3, rasterized=True)
    ax.axhline(-np.log10(GWS), color="#D7263D", lw=0.6, ls=(0, (4, 2)))
    ax.axhline(-np.log10(SUG), color=S.GREY, lw=0.5, ls=(0, (1.5, 1.5)))
    ax.set_xticks(mids)
    ax.set_xticklabels([str(c) if c <= 18 or c % 2 == 0 else "" for c in range(1, 23)])
    ax.set_xlim(-1e7, start)
    top = ymax or max(9, y.max() + 1)
    ax.set_ylim(0, top)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel(r"$-\log_{10}(P)$")
    ax.text(0.995, 0.97, PH[ph], transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
            fontweight="bold", color=S.FEAT[ph])
    if label_top:
        i = np.argmax(y)
        ax.annotate(label_top, (x[i], y[i]), xytext=(x[i] + 1.2e8, y[i] - 0.3), fontsize=6.5,
                    style="italic", arrowprops=dict(arrowstyle="-", lw=0.5, color=S.INK), va="center")


def regional(ax, gax):
    loc = pd.read_csv(SRC["locus"])
    ls = pd.read_csv(SRC["locus_snps"])
    r2 = dict(zip(ls.pos, ls.r2))
    loc["r2"] = loc.BP.map(r2)
    lead = 37495010
    loc["y"] = -np.log10(loc.P)
    bins = [(0.8, 1.01, "#D7263D", "0.8–1.0"), (0.6, 0.8, "#F28C28", "0.6–0.8")]
    ax.scatter(loc.BP / 1e6, loc.y, s=5, c="#B8C0CC", lw=0, rasterized=True)
    hs = [Line2D([], [], marker="o", ls="", ms=3.5, mfc="#B8C0CC", mec="none", label="< 0.6 / not in LD set")]
    for lo, hi, c, lab in bins[::-1]:
        m = (loc.r2 >= lo) & (loc.r2 < hi) & (loc.BP != lead)
        ax.scatter(loc.BP[m] / 1e6, loc.y[m], s=9, c=c, lw=0.2, edgecolors="white", zorder=3)
        hs.insert(0, Line2D([], [], marker="o", ls="", ms=3.5, mfc=c, mec="none", label=lab))
    L = loc[loc.BP == lead].iloc[0]
    ax.scatter(L.BP / 1e6, L.y, marker="D", s=26, c="#6A1B9A", edgecolors="white", lw=0.5, zorder=4)
    ax.annotate("rs1549596", (L.BP / 1e6, L.y), xytext=(L.BP / 1e6 - 0.36, L.y + 0.35), fontsize=6.5,
                arrowprops=dict(arrowstyle="-", lw=0.5))
    ax.axhline(-np.log10(GWS), color="#D7263D", lw=0.6, ls=(0, (4, 2)))
    ax.set_ylabel(r"$-\log_{10}(P)$")
    ax.set_ylim(0, 9.2)
    ax.legend(handles=hs, title="r² with rs1549596 (EAS)", loc="upper right", handletextpad=0.2, borderaxespad=0.1, fontsize=5.8, title_fontsize=6)
    ax.tick_params(labelbottom=False)
    # gene track
    g = pd.read_csv(SRC["fuma_genes"], sep="\t")
    g = g[(g.end > loc.BP.min()) & (g.start < loc.BP.max()) & (g.type == "protein_coding")].sort_values("start")
    lanes = []
    for _, r in g.iterrows():
        for k, e in enumerate(lanes):
            if r.start > e + 2.5e5:   # v3: wider gap so gene labels do not overlap
                lanes[k] = r.end; lane = k; break
        else:
            lanes.append(r.end); lane = len(lanes) - 1
        y = -lane
        hl = r.symbol == "ITGA9"
        col = "#6A1B9A" if hl else "#4D5B6B"
        gax.plot([r.start / 1e6, r.end / 1e6], [y, y], color=col, lw=2.2 if hl else 1.6, solid_capstyle="butt")
        gax.text((r.start + r.end) / 2e6, y - 0.45, r.symbol, ha="center", va="top", fontsize=5.6,
                 style="italic", color=col, fontweight="bold" if hl else "normal")
    gax.set_ylim(-len(lanes) - 0.3, 0.6)
    gax.set_yticks([])
    gax.spines["left"].set_visible(False)
    gax.set_xlim(loc.BP.min() / 1e6, loc.BP.max() / 1e6)
    ax.set_xlim(gax.get_xlim())
    gax.set_xlabel("Chromosome 3 position (Mb, GRCh37)")


def fig1():
    fig = plt.figure(figsize=(S.COL2, 140 * S.MM))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 0.22, 1.15], width_ratios=[1.1, 1], hspace=0.05, wspace=0.22,
                          left=0.07, right=0.99, top=0.97, bottom=0.07)
    ax = fig.add_subplot(gs[0, :])
    manhattan(ax, "MC", label_top="ITGA9")
    S.panel(ax, "A", x=-0.06)
    sub = gs[2, 0].subgridspec(2, 1, height_ratios=[3, 1], hspace=0.06)
    ra, ga = fig.add_subplot(sub[0]), fig.add_subplot(sub[1])
    regional(ra, ga)
    S.panel(ra, "B", x=-0.12)
    fm_axes = gs[2, 1].subgridspec(1, 2, width_ratios=[0.75, 1.6], wspace=0.75)
    finemap(fig.add_subplot(fm_axes[0]), fig.add_subplot(fm_axes[1]))
    S.save(fig, OUT, "Fig1_MC_GWAS_ITGA9_locus")


def finemap(a1, a2, letters=("C", "D")):
    nc = pd.read_csv(SRC["fm_ncausal"], sep=" ")
    post = nc[nc.type == "post"]
    prior = nc[nc.type == "prior"]
    k = post.ncausal_num.values
    a1.bar(k - 0.2, prior.set_index("ncausal_num").reindex(k).ncausal_prob, 0.4, color="#C9CCD1", label="Prior")
    a1.bar(k + 0.2, post.ncausal_prob, 0.4, color=S.FEAT["MC"], label="Posterior")
    a1.set_xticks(k)
    a1.set_xlabel("No. causal variants")
    a1.set_ylabel("Probability")
    a1.set_ylim(0, 1)
    a1.legend(loc="upper right", bbox_to_anchor=(1.05, 1.02), fontsize=5.8)
    S.panel(a1, letters[0], x=-0.55, y=1.03)
    sn = pd.read_csv(SRC["fm_snp"], sep=" ").sort_values("snp_prob", ascending=True).tail(12)
    y = np.arange(len(sn))
    c = [S.FEAT["MC"] if s == "rs1549596" else "#9AA5B4" for s in sn.snp]
    a2.hlines(y, 0, sn.snp_prob, color=c, lw=1.2)
    a2.scatter(sn.snp_prob, y, c=c, s=14, zorder=3)
    a2.set_yticks(y)
    a2.set_yticklabels(sn.snp, fontsize=6)
    a2.set_xlabel("Posterior inclusion probability")
    a2.set_xlim(0, 0.45)
    a2.tick_params(axis="y", length=0)
    a2.spines["left"].set_visible(False)
    for yy, p in zip(y, sn.snp_prob):
        a2.text(p + 0.012, yy, f"{p:.3f}", va="center", fontsize=5.6, color=S.INK)
    S.panel(a2, letters[1], x=-0.62, y=1.03)


def fig_qq():
    fig, axs = plt.subplots(1, 4, figsize=(S.COL2, 52 * S.MM), gridspec_kw=dict(wspace=0.35))
    for ax, ph, L in zip(axs, ["Core", "MC", "HIZ", "DH"], "ABCD"):
        q = pd.read_csv(WORK + f"qq_{ph}.csv")
        mx = max(q.exp.max(), q.obs.max()) + 0.3
        ax.plot([0, mx], [0, mx], color=S.GREY, lw=0.6)
        ax.scatter(q.exp, q.obs, s=2.5, c=S.FEAT[ph], lw=0, rasterized=True)
        ax.set_xlim(0, q.exp.max() + 0.3); ax.set_ylim(0, max(q.obs.max(), q.exp.max()) + 0.5)
        ax.set_title(PH[ph], color=S.FEAT[ph])
        ax.text(0.05, 0.93, f"λ = {LAMBDA[ph]:.3f}", transform=ax.transAxes, fontsize=6.5, va="top")
        ax.set_xlabel(r"Expected $-\log_{10}(P)$")
        if ph == "Core":
            ax.set_ylabel(r"Observed $-\log_{10}(P)$")
        ax.set_aspect("auto")
        S.panel(ax, L, x=-0.3, y=1.02)
    S.save(fig, OUT, "FigS1_QQ_plots")


def fig_s2():
    fig, axs = plt.subplots(3, 1, figsize=(S.COL2, 140 * S.MM), gridspec_kw=dict(hspace=0.45))
    for ax, ph, L in zip(axs, ["Core", "HIZ", "DH"], "ABC"):
        manhattan(ax, ph, ymax=9)
        S.panel(ax, L, x=-0.06)
    S.save(fig, OUT, "FigS2_Manhattan_Core_HIZ_DH")


if __name__ == "__main__":
    fig1(); fig_qq(); fig_s2()
