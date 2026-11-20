"""Paper 2 Figure 2A (v3, 2026-11-06): ITGA9 median expression across the 54 GTEx v8 tissues.
Data: GTEx Portal API v2, expression/medianGeneExpression, gencodeId ENSG00000144668.11, datasetId gtex_v8
(saved as results/p2_ITGA9_GTEx_v8_median.json)."""
import os, sys, json, pandas as pd, matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__)); import style as S
from common_v2 import RES
S.setup()
d = pd.DataFrame(json.load(open(os.path.join(RES, "p2_ITGA9_GTEx_v8_median.json")))["data"])
d["label"] = d.tissueSiteDetailId.str.replace("_", " ").str.replace("Brain ", "Brain – ").str.replace("Cells ", "Cells – ")
d = d.sort_values("median", ascending=True)
def grp(t):
    if t.startswith("Artery") or t.startswith("Heart"): return "Cardiovascular", "#C8443A"
    if t.startswith("Brain") or t in ("Pituitary", "Nerve_Tibial"): return "Nervous system", "#6F8FA6"
    if t in ("Whole_Blood", "Spleen", "Cells_EBV-transformed_lymphocytes"): return "Blood and immune", "#7461AA"
    if t.startswith("Esophagus_Musc") or t.startswith("Esophagus_Gastro") or t.startswith("Colon_Sig") or t in ("Uterus", "Bladder", "Cervix_Endocervix", "Cervix_Ectocervix", "Fallopian_Tube", "Prostate", "Vagina"): return "Smooth-muscle-rich", "#E3912B"
    if t.startswith("Adipose") or t in ("Muscle_Skeletal", "Cells_Cultured_fibroblasts", "Breast_Mammary_Tissue"): return "Adipose, muscle, fibroblast", "#3A9A5B"
    return "Other", "#A6A6A6"
d[["grp", "col"]] = d.tissueSiteDetailId.apply(lambda t: pd.Series(grp(t)))
fig, ax = plt.subplots(figsize=(S.COL15, 150 * S.MM))
ax.barh(range(len(d)), d["median"], color=d.col, height=0.75, lw=0)
ax.set_yticks(range(len(d))); ax.set_yticklabels(d.label, fontsize=5.2); ax.set_ylim(-0.7, len(d) - 0.3)
ax.set_xlabel("ITGA9 expression, median TPM (GTEx v8)"); ax.xaxis.grid(True); ax.set_axisbelow(True)
ax.tick_params(axis="y", length=0)
from matplotlib.patches import Patch
seen = d.drop_duplicates("grp")[["grp", "col"]]
order = ["Cardiovascular", "Smooth-muscle-rich", "Blood and immune", "Adipose, muscle, fibroblast", "Nervous system", "Other"]
seen = seen.set_index("grp").loc[[g for g in order if g in set(seen.grp)]]
ax.legend([Patch(color=c) for c in seen.col], seen.index, loc="lower right", fontsize=5.8, frameon=False)
fig.subplots_adjust(left=0.36, right=0.98, top=0.99, bottom=0.06)
S.save(fig, os.path.join(RES, "..", "..", "figures_v2", "Paper2"), "Fig2A_ITGA9_GTEx")
