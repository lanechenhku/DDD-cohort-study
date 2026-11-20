"""Supplementary Figure S2: distribution of imputed individual onset ages by feature (v3: SISE re-run, results/p1_imputed_onset_ages.csv)."""
import os, sys, pandas as pd, numpy as np, matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__)); import style as S
from common_v2 import RES
S.setup()
d = pd.read_csv(sys.argv[1] if len(sys.argv) > 1 else os.path.join(RES, "p1_imputed_onset_ages.csv"))
feats = ["SS", "DB", "DN", "DH", "HIZ", "MC", "SN"]
fig, axs = plt.subplots(2, 4, figsize=(S.COL2, 85 * S.MM), sharex=True)
for ax, f in zip(axs.flat, feats):
    v = d[f].dropna()
    ax.hist(v, bins=np.arange(0, 101, 2.5), color=S.FEAT[f], lw=0)
    ax.set_title(f"{f} (n = {len(v):,}; median {v.median():.1f})", fontweight="normal", fontsize=6.5)
    ax.yaxis.grid(True); ax.set_axisbelow(True)
for ax in axs[1]: ax.set_xlabel("Imputed onset age (years)")
for ax in axs[:, 0]: ax.set_ylabel("Participants")
axs.flat[-1].axis("off")
fig.tight_layout()
S.save(fig, os.path.join(RES, "..", "..", "figures_v2", "Paper1"), "FigS2_imputed_onset_ages")
print(d[feats].describe().round(1).to_string())
