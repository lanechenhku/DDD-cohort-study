"""ITGA9 expression across single-cell types (Human Protein Atlas v24 rna_single_cell_type.tsv, downloaded 2026-10-23)."""
import sys, os, pandas as pd, matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__)); import style as S
from common_v2 import RES
S.setup()
d = pd.read_csv(os.path.join(RES, "p2_ITGA9_HPA_singlecell.tsv"), sep="\t", header=None, names=["ens", "gene", "cell", "nCPM"])
d = d.sort_values("nCPM", ascending=False)
top = d.head(30).iloc[::-1]
hi = {"hematopoietic stem cells", "megakaryocyte-erythroid progenitors", "mast cells", "fibroblasts", "smooth muscle cells",
      "chondrocytes", "osteoblasts", "adipocytes", "fibro-adipogenic progenitors", "endothelial cells", "vascular endothelial cells", "lymphatic endothelial cells"}
fig, ax = plt.subplots(figsize=(S.COL1 * 1.15, 105 * S.MM))
ax.barh(range(len(top)), top.nCPM, color=["#7461AA" if c in hi else "#C4C8CE" for c in top.cell], height=0.7)
ax.set_yticks(range(len(top))); ax.set_yticklabels([c[0].upper() + c[1:] for c in top.cell], fontsize=5.5)
ax.set_xlabel("ITGA9 expression (nCPM)"); ax.set_title("ITGA9 across single-cell types (HPA; top 30 of %d)" % len(d), fontweight="normal")
S.save(fig, os.path.join(RES, "..", "..", "figures_v2", "Paper2"), "FigS4_ITGA9_singlecell")
print(d.head(12).to_string()); print("rank of chondrocytes etc:", d.reset_index(drop=True).query("cell.str.contains('chondro|osteo|fibro|adipo|smooth')", engine="python")[["cell","nCPM"]].to_string())
