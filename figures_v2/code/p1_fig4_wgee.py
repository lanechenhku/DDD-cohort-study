"""Paper 1 Figure 4: WGEE risk-factor estimates (Models 1-3).

Coefficients were digitised from the 2023 figures (Figures/wgee_m*_updated.tiff)
because the model output tables are not in DDD_Paper; replace data/wgee_*.csv
with the geepack output when available and rerun.
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import style as S

S.setup()
from common import DATA as D, P1
m1 = pd.read_csv(D + "wgee_m1_digitized.csv")
m23 = pd.read_csv(D + "wgee_m23_digitized.csv")
feats = ["Core", "HIZ", "MC", "DH"]
terms = ["Age", "Age²", "Sex (male)", "Core onset age", "BMI", "Height", "Smoking",
         "Lumbar injury", "Sports", "Age × sex (male)", "Age × BMI", "Age × height",
         "Age × smoking", "Age × injury", "Age × sports"]

fig = plt.figure(figsize=(S.COL2, 150 * S.MM))
gs = fig.add_gridspec(2, 4, height_ratios=[1.55, 1], hspace=0.42, wspace=0.10,
                      left=0.16, right=0.99, top=0.93, bottom=0.07)
y = np.arange(len(terms))[::-1]
axes = []
for j, f in enumerate(feats):
    ax = fig.add_subplot(gs[0, j])
    axes.append(ax)
    d = m1[m1.feature == f].set_index("term")
    ax.axvline(0, color=S.GREY, lw=0.6, ls=(0, (3, 2)), zorder=0)
    for t, yy in zip(terms, y):
        if yy % 2 == 0:
            ax.axhspan(yy - 0.5, yy + 0.5, color="#F4F5F7", zorder=-1, lw=0)
        if t not in d.index:
            ax.text(0, yy, "n/a", ha="center", va="center", fontsize=5.5, color=S.GREY)
            continue
        r = d.loc[t]
        sig = r.lo > 0 or r.hi < 0
        c = S.FEAT[f]
        ax.plot([r.lo, r.hi], [yy, yy], color=c, lw=1.1, solid_capstyle="round")
        ax.plot(r.est, yy, "o", ms=3.6, mfc=c if sig else "white", mec=c, mew=0.9)
    ax.set_title(f, color=S.FEAT[f], fontsize=8)
    ax.set_ylim(-0.6, len(terms) - 0.4)
    ax.set_yticks(y)
    ax.set_yticklabels(terms if j == 0 else [])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("β (95% CI)")
    ax.set_xlim(-0.58, 0.62)
    ax.set_xticks([-0.4, -0.2, 0, 0.2, 0.4])
S.panel(axes[0], "A", x=-0.75, y=1.02)
axes[0].text(-0.62, -0.13, "", transform=axes[0].transAxes)

def forest(ax, d, letter, title):
    d = d.reset_index(drop=True)
    yy = np.arange(len(d))[::-1]
    ax.axvline(0, color=S.GREY, lw=0.6, ls=(0, (3, 2)), zorder=0)
    for k, r in d.iterrows():
        c = S.FEAT[r.feature]
        if yy[k] % 2 == 0:
            ax.axhspan(yy[k] - 0.5, yy[k] + 0.5, color="#F4F5F7", zorder=-1, lw=0)
        ax.plot([r.lo, r.hi], [yy[k], yy[k]], color=c, lw=1.2, solid_capstyle="round")
        ax.plot(r.est, yy[k], "o", ms=4, mfc=c, mec=c)
        ax.text(1.01, yy[k], f"{r.est:.3f} ({r.lo:.3f}, {r.hi:.3f})", transform=ax.get_yaxis_transform(),
                va="center", fontsize=5.8, color=S.INK)
    ax.set_yticks(yy)
    ax.set_yticklabels([f"{t}  →  {f}" for t, f in zip(d.term, d.feature)])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_ylim(-0.6, len(d) - 0.4)
    ax.set_xlabel("β (95% CI)")
    ax.set_title(title, loc="left", fontweight="normal", x=-0.0)
    S.panel(ax, letter, x=-0.95 if letter == "B" else -1.05, y=1.04)

axB = fig.add_subplot(gs[1, 0:2])
forest(axB, m23[m23.model == "M2"], "B", "Model 2: each PGS modelled independently")
axC = fig.add_subplot(gs[1, 2:4])
forest(axC, m23[m23.model == "M3"], "C", "Model 3: PGS modelled jointly")
axB.set_xlim(-0.13, 0.15)
axC.set_xlim(-0.02, 0.15)
fig.subplots_adjust()
# second-row axes need room for long labels and the estimate column
axB.set_position([0.205, 0.07, 0.22, 0.25])
axC.set_position([0.695, 0.07, 0.17, 0.25])
h = [Line2D([], [], marker="o", ls="-", color=S.FEAT[f], mfc=S.FEAT[f], ms=4, label=f) for f in feats]
h += [Line2D([], [], marker="o", ls="", color=S.GREY, mfc="white", ms=4, label="95% CI includes 0")]
fig.legend(handles=h, loc="upper center", ncol=5, bbox_to_anchor=(0.57, 1.005), handlelength=1.6, columnspacing=1.6)
S.save(fig, P1, "Fig4_WGEE_risk_factors")
