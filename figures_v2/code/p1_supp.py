"""Paper 1 supplementary figures S1 (sex-stratified distributions), S3 (correlation
of composite phenotypes) and S4 (level-wise change T1 -> T2)."""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import matplotlib.colors as mcol
import style as S
from common import SRC, P1

S.setup()
PH = ["Core", "DH", "HIZ", "MC"]
NAME = {"Core": "Core.slope", "DH": "DH.avg", "HIZ": "HIZ.avg", "MC": "MC.avg"}
SEXC = {"Female": "#C8443A", "Male": "#2C6FB7"}


def kde(x, grid, bw=None):
    x = x[~np.isnan(x)]
    bw = bw or 1.06 * x.std() * len(x) ** -0.2
    z = (grid[:, None] - x[None, :]) / bw
    return np.exp(-0.5 * z * z).sum(1) / (len(x) * bw * np.sqrt(2 * np.pi))


def s1(g):
    fig, axs = plt.subplots(1, 4, figsize=(S.COL2, 50 * S.MM), gridspec_kw=dict(wspace=0.35))
    for ax, p, L in zip(axs, PH, "ABCD"):
        lo, hi = np.nanpercentile(g[p], [0.2, 99.8])
        grid = np.linspace(lo, hi, 400)
        bw = (hi - lo) / 40
        for sx, code in (("Female", 1), ("Male", 0)):  # sex: 1 = female (59.8% of cohort)
            v = g.loc[g.sex == code, p].values.astype(float)
            y = kde(v, grid, bw)
            ax.fill_between(grid, y, color=SEXC[sx], alpha=0.12, lw=0)
            ax.plot(grid, y, color=SEXC[sx], lw=1, label=f"{sx} (n = {np.isfinite(v).sum():,})")
            ax.axvline(np.nanmean(v), color=SEXC[sx], lw=0.6, ls=(0, (2, 1.5)))
        ax.set_xlabel(NAME[p])
        ax.set_title(p, color=S.FEAT[p])
        if p == "Core":
            ax.set_ylabel("Density")
        S.panel(ax, L, x=-0.3, y=1.02)
    axs[0].legend(loc="upper right", fontsize=5.8)
    S.save(fig, P1, "FigS1_sex_distributions")


def s3(g):
    C = g[PH].corr(method="spearman").values  # pattern matches the text values; Pearson does not
    fig, ax = plt.subplots(figsize=(S.COL1 * 0.8, S.COL1 * 0.72))
    cmap = mcol.LinearSegmentedColormap.from_list("d", ["#B2443A", "#F7F7F7", "#2C6FB7"])
    m = np.triu(np.ones_like(C, bool), 1)
    for i in range(4):
        for j in range(4):
            if j <= i:
                continue
            ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, color=cmap((C[i, j] + 0.5) / 1.0), ec="white", lw=1.5))
            ax.text(j, i, f"{C[i, j]:.2f}".replace("-", "−"), ha="center", va="center", fontsize=7.5, fontweight="bold",
                    color=S.INK)
    ax.set_xlim(0.5, 3.5); ax.set_ylim(2.5, -0.5)
    ax.set_xticks([1, 2, 3]); ax.set_xticklabels([NAME[p] for p in PH[1:]]); ax.xaxis.tick_top()
    ax.set_yticks([0, 1, 2]); ax.set_yticklabels([NAME[p] for p in PH[:3]])
    for t in ax.get_xticklabels() + ax.get_yticklabels():
        t.set_color(S.FEAT[t.get_text().split(".")[0]]); t.set_fontweight("bold")
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=mcol.Normalize(-0.5, 0.5))
    cb = fig.colorbar(sm, ax=ax, fraction=0.05, pad=0.04, shrink=0.8)
    cb.set_label("Spearman ρ"); cb.outline.set_visible(False)
    S.save(fig, P1, "FigS3_composite_correlation")
    return C


def s4():
    d = pd.read_csv(SRC["level"])
    feats = ["ss", "db", "hiz", "mc", "sn"]
    lab = {"ss": "SS", "db": "DB", "hiz": "HIZ", "mc": "MC", "sn": "SN"}
    fig, axs = plt.subplots(1, 5, figsize=(S.COL2, 62 * S.MM), sharey=True, gridspec_kw=dict(wspace=0.12))
    cats = ["Regression", "", "No change", "", "Progression"]
    colors = ["#2C6FB7", "#2C6FB7", "#E4E6EA", "#C8443A", "#C8443A"]
    for ax, f, L in zip(axs, feats, "ABCDE"):
        for l in range(1, 6):
            x = (d[f"{f}{l}.2"] - d[f"{f}{l}.1"]).dropna()
            b = [(x <= -1).mean(), ((x > -1) & (x < 0)).mean(), (x == 0).mean(), ((x > 0) & (x < 1)).mean(), (x >= 1).mean()]
            # centre the stable share on zero
            left = -(b[0] + b[1] + b[2] / 2) * 100
            for v, c in zip(b, colors):
                ax.barh(l, v * 100, left=left, color=c, height=0.7, edgecolor="none")
                left += v * 100
            prog = (b[3] + b[4]) * 100; reg = (b[0] + b[1]) * 100
            ax.text(left + 1.5, l, f"{prog:.0f}%", va="center", fontsize=5.5, color="#C8443A")
        ax.axvline(0, color=S.GREY, lw=0.5)
        ax.set_title(lab[f], color=S.FEAT[lab[f]])
        ax.set_xlim(-75, 75)
        ax.set_xticks([-50, 0, 50]); ax.set_xticklabels(["50", "0", "50"])
        ax.set_xlabel("Participants (%)")
        S.panel(ax, L, x=-0.08 if f != "ss" else -0.3, y=1.02)
    axs[0].set_yticks(range(1, 6)); axs[0].set_yticklabels([f"L{l}/{l+1 if l < 5 else 'S1'}" for l in range(1, 6)])
    axs[0].invert_yaxis()
    from matplotlib.patches import Patch
    h = [Patch(color=c, label=k) for c, k in zip(colors, cats) if k]
    fig.legend(handles=h, title="Level-wise change in feature score between first and second MRI (% progressing shown at right)", ncol=3,
               loc="upper center", bbox_to_anchor=(0.5, 1.08))
    S.save(fig, P1, "FigS4_levelwise_change")


if __name__ == "__main__":
    g = pd.read_csv(SRC["gwas_pheno"], sep=r"\s+")
    s1(g); print(s3(g).round(2)); s4()
