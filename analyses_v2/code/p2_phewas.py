"""rs1549596 PheWAS: FinnGen R12 PheWeb (2,470 endpoints; API returns p only for p < ~0.05) and BioBank Japan PheWeb (324 traits).
Inputs downloaded 2026-10-23 from https://r12.finngen.fi/api/variant/3:37453519-C-T and https://pheweb.jp/api/variant/3-37495010-C-T"""
import json, sys, os, numpy as np, pandas as pd, matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__)); import style as S
from common_v2 import RES
S.setup()
SRC = sys.argv[1] if len(sys.argv) > 1 else "/tmp/lk"
fg = pd.DataFrame(json.load(open(f"{SRC}/fg_phewas.json"))["results"])
nfg = len(fg); fg = fg[fg.pval.notna()].copy()
bbj = pd.DataFrame(json.load(open(f"{SRC}/bbj_phewas.json"))["phenos"]); bbj = bbj[bbj.pval.notna()].copy()
bbj = bbj[bbj.category != "GxE"].copy()   # v3: drop the 47 gene-by-environment combination tests (not traits)
fg["src"], bbj["src"] = "FinnGen R12", "BioBank Japan"
fg["cat"] = fg.category.str.replace(r"\s*\(.*\)$", "", regex=True); bbj["cat"] = bbj.category
out = pd.concat([fg[["src","phenocode","phenostring","cat","beta","pval","n_case"]].rename(columns={"n_case":"cases"}),
                 bbj[["src","phenocode","phenostring","cat","beta","pval","num_cases"]].rename(columns={"num_cases":"cases"})])
out.sort_values("pval").to_csv(os.path.join(RES, "p2_rs1549596_phewas.csv"), index=False)
fig, axs = plt.subplots(1, 2, figsize=(S.COL2, 70 * S.MM), gridspec_kw=dict(wspace=0.25))
for ax, d, n, title in [(axs[0], fg, nfg, f"FinnGen R12 ({nfg:,} endpoints; p < 0.05 shown)"),
                        (axs[1], bbj, len(bbj), f"BioBank Japan ({len(bbj)} traits)")]:
    cats = list(dict.fromkeys(d.sort_values("cat").cat)); ci = {c: i for i, c in enumerate(cats)}
    d = d.sort_values("cat").reset_index(drop=True)
    cols = plt.cm.tab20(np.arange(len(cats)) % 20)
    ax.scatter(d.index, -np.log10(d.pval), s=7, c=[cols[ci[c]] for c in d.cat], lw=0)
    ax.axhline(-np.log10(0.05 / n), color=S.GREY, ls="--", lw=0.6)
    for k, (i, r) in enumerate(d.nsmallest(2, "pval").iterrows()):   # v3: arrows so each label points at its own dot
        right = i > len(d) / 2
        ax.annotate(r.phenostring if len(r.phenostring) <= 40 else r.phenostring[:38] + "…", (i, -np.log10(r.pval)),
                    xytext=(4 if right else -4, 10 + 12 * k), textcoords="offset points", fontsize=5,
                    ha="right" if right else "left", va="bottom",
                    arrowprops=dict(arrowstyle="-", lw=0.4, color=S.GREY, shrinkA=0, shrinkB=2))
    ax.set_xticks([]); ax.set_ylabel("−log$_{10}$ p"); ax.set_title(title, fontweight="normal")
    ax.set_ylim(0, max(-np.log10(0.05 / n) + 0.6, -np.log10(d.pval.min()) + 0.8))
S.save(fig, os.path.join(os.path.dirname(RES), "..", "figures_v2", "Paper2"), "FigS5_rs1549596_PheWAS")
print(out.sort_values("pval").head(8).to_string())
