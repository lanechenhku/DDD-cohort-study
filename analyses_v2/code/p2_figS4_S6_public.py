"""Paper 2 Supplementary Figures S4 and S6 from public data only (v3, 2026-11-06; Lane: replace the unpublished
Cheah-group single-cell panel and the unsourced annotation screenshot).
S4  ITGA9 expression: (A) GTEx v8 tissues (p2_itga9_gtex.py); (B, C) human nucleus pulposus scRNA-seq, Tu et al. 2022
    (GSE165722; p2_np_scrna_load.py, p2_np_scrna_process.py); (D) human bone marrow, CZ CELLxGENE Discover gene-expression
    API (results/p2_ITGA9_CELLxGENE_bonemarrow.csv).
S6  Regional chromatin context of rs1549596 (GRCh37): MC.avg GWAS, GENCODE genes (UCSC knownGene), GeneHancer double-elite
    elements, ENCODE DNase clusters (V3), Roadmap 15-state ChromHMM (core marks) for 13 reference epigenomes, fetched via the
    UCSC API from the Roadmap Integrative Analysis Hub (results/p2_ITGA9_Roadmap_chromHMM15_region.csv, results/p2_s6_ucsc/)."""
import os, sys, json, numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.image import imread
sys.path.insert(0, os.path.dirname(__file__)); import style as S
from common_v2 import RES, UP
S.setup(); OUT = os.path.join(RES, "..", "..", "figures_v2", "Paper2")
SNP = 37495010
# ---------------- S4 ----------------
c = pd.read_csv(os.path.join(RES, "p2_ITGA9_NP_Tu2022_cells.csv.gz")); c = c[c.ct != "Mixed (doublets)"]
order = ["NP chondrocyte-like (DCN/CLU)", "NP chondrocyte-like (FMOD/PRG4)", "Fibroblast-like (COL1/COL3/TNC)", "Macrophages",
         "Neutrophils", "Myelocytes / immature neutrophils", "Promyelocytes (proliferating)", "Erythroid cells", "T / NK cells",
         "B cells, endothelial and other"]
pal = ["#2C6FB7", "#6F8FA6", "#3A9A5B", "#E3912B", "#D4508A", "#C8443A", "#7461AA", "#7A1F3D", "#A6A6A6", "#2A9D8F"]
col = dict(zip(order, pal))
fig = plt.figure(figsize=(S.COL2, 215 * S.MM))
gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1, 0.95], hspace=0.95, wspace=0.75, left=0.07, right=0.98, top=0.98, bottom=0.05)
ax = fig.add_subplot(gs[0, :])
import importlib.util as _iu
gt = pd.DataFrame(json.load(open(os.path.join(RES, "p2_ITGA9_GTEx_v8_median.json")))["data"]).sort_values("median", ascending=False)
_sp = _iu.spec_from_file_location("gx", os.path.join(os.path.dirname(__file__), "p2_itga9_gtex.py"))
src = open(os.path.join(os.path.dirname(__file__), "p2_itga9_gtex.py")).read(); ns = {}; exec(src[src.index("def grp"):src.index("d[[")], {"pd": pd}, ns)
gc = [ns["grp"](t) for t in gt.tissueSiteDetailId]
ax.bar(range(len(gt)), gt["median"], color=[c for _, c in gc], width=0.75, lw=0)
ax.set_xticks(range(len(gt))); ax.set_xticklabels(gt.tissueSiteDetailId.str.replace("_", " ").str.replace("Brain ", "Brain – "), rotation=60, ha="right", fontsize=4.8)
ax.set_xlim(-0.7, len(gt) - 0.3); ax.set_ylabel("ITGA9, median TPM"); ax.yaxis.grid(True); ax.set_axisbelow(True)
from matplotlib.patches import Patch as _P
seen = list(dict.fromkeys(gc)); ax.legend([_P(color=c) for _, c in seen], [g for g, _ in seen], fontsize=5.6, ncol=3, loc="upper right")
ax.set_title("Bulk tissue expression (GTEx v8, 54 tissues)", loc="left", fontsize=6.5); S.panel(ax, "A", x=-0.05, y=1.02)
# B UMAP
ax = fig.add_subplot(gs[1, 0])
for k in order:
    s = c[c.ct == k]; ax.scatter(s.umap1, s.umap2, s=0.4, c=col[k], lw=0, rasterized=True)
e = c[c.itga9_detected == 1]; ax.scatter(e.umap1, e.umap2, s=2.2, c="black", lw=0, rasterized=True, label="ITGA9 detected")
for k in order:
    s = c[c.ct == k]; ax.text(s.umap1.median(), s.umap2.median(), str(order.index(k) + 1), fontsize=6.5, fontweight="bold", ha="center", va="center",
                              color="white", bbox=dict(boxstyle="circle,pad=0.15", fc=col[k], ec="none"))
ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2"); ax.legend(loc="lower left", fontsize=5.6, markerscale=3, handletextpad=0.2)
ax.set_title(f"Human nucleus pulposus scRNA-seq\n(Tu et al. 2022; {len(c):,} cells, 8 donors)", loc="left", fontsize=6.5); S.panel(ax, "B", x=-0.08, y=1.02)
# C dot/bar of % detected by cell type
ax = fig.add_subplot(gs[1, 1])
t = c.groupby("ct").agg(n=("itga9_detected", "size"), pct=("itga9_detected", "mean")).reindex(order); t["pct"] *= 100
y = np.arange(len(order))[::-1]
ax.barh(y, t.pct, color=[col[k] for k in order], height=0.7, lw=0)
for yy, (k, r) in zip(y, t.iterrows()): ax.text(r.pct + 0.3, yy, f"{r.pct:.1f}%  (n={int(r.n):,})", va="center", fontsize=5.4)
ax.set_yticks(y); ax.set_yticklabels([f"{i + 1}  {k}" for i, k in enumerate(order)], fontsize=5.6); ax.tick_params(axis="y", length=0)
ax.set_xlim(0, 19); ax.set_xlabel("Cells with ITGA9 detected (%)"); ax.set_title("ITGA9 by cell type, nucleus pulposus", loc="left", fontsize=6.5)
S.panel(ax, "C", x=-0.9, y=1.02)
# D bone marrow (CELLxGENE)
ax = fig.add_subplot(gs[2, :])
b = pd.read_csv(os.path.join(RES, "p2_ITGA9_CELLxGENE_bonemarrow.csv"))
keep = ["endothelial cell", "mesenchymal cell", "neutrophil", "macrophage", "myelocyte", "preosteoblast", "hematopoietic stem cell",
        "common myeloid progenitor", "promyelocyte", "megakaryocyte-erythroid progenitor cell", "erythroblast", "monocyte",
        "plasmacytoid dendritic cell", "B cell", "T cell", "natural killer cell"]
b = b[b.label.isin(keep) & (b.order >= 0)].drop_duplicates("label").set_index("label").reindex([k for k in keep if k in set(b.label)]).dropna(subset=["total"])
bc = ["#C8443A" if "endothelial" in k else "#3A9A5B" if k in ("mesenchymal cell", "preosteoblast") else "#7461AA" for k in b.index]
x = np.arange(len(b)); ax.bar(x, b.pct, color=bc, width=0.7, lw=0)
for xx, (k, r) in zip(x, b.iterrows()): ax.text(xx, r.pct + 0.3, f"{r.pct:.1f}", ha="center", fontsize=5.2)
ax.set_xticks(x); ax.set_xticklabels([k[0].upper() + k[1:] + f"\n(n={int(r.total):,})" for k, r in b.iterrows()], rotation=40, ha="right", fontsize=5.3)
ax.set_ylabel("Cells with ITGA9\ndetected (%)"); ax.set_ylim(0, 21)
ax.set_title("ITGA9 by cell type, human bone marrow (CZ CELLxGENE Discover; red endothelial, green mesenchymal, purple haematopoietic)", loc="left", fontsize=6.5)
S.panel(ax, "D", x=-0.04, y=1.02)
S.save(fig, OUT, "FigS4_ITGA9_expression_public")

# ---------------- S6 ----------------
L0, L1 = 37470000, 37560000
g = pd.read_csv(os.path.join(UP, "Analyses_v2_2026-10-23/data/itga9_region_MC.csv")); g = g[(g.BP > L0) & (g.BP < L1)]
rm = pd.read_csv(os.path.join(RES, "p2_ITGA9_Roadmap_chromHMM15_region.csv"))
EPI = [("E035", "Primary HSC (CD34+)"), ("E050", "Mobilised HSC (CD34+)"), ("E029", "CD14+ monocytes"), ("E062", "PB mononuclear cells"),
       ("E032", "B cells"), ("E034", "T cells"), ("E026", "Bone-marrow MSC"), ("E049", "MSC-derived chondrocytes"), ("E129", "Osteoblasts"),
       ("E025", "Adipose-derived MSC"), ("E063", "Adipose nuclei"), ("E065", "Aorta"), ("E122", "HUVEC (endothelial)")]
U = os.path.join(RES, "p2_s6_ucsc")
fig = plt.figure(figsize=(S.COL2, 150 * S.MM))
gs = fig.add_gridspec(5, 1, height_ratios=[1.3, 0.35, 0.3, 0.3, 2.6], hspace=0.12, left=0.2, right=0.98, top=0.97, bottom=0.12)
ax = fig.add_subplot(gs[0]); ax.scatter(g.BP, -np.log10(g.P), s=5, c=S.GREY, lw=0)
l = g[g.BP == SNP]; ax.scatter(l.BP, -np.log10(l.P), s=26, marker="D", c=S.FEAT["MC"], zorder=3); ax.text(SNP + 900, -np.log10(l.P.iloc[0]), "rs1549596", fontsize=6, va="center")
ax.axhline(-np.log10(5e-8), color="#C8443A", ls="--", lw=0.6); ax.set_ylabel("MC.avg\n−log$_{10}$ P"); ax.set_xlim(L0, L1); ax.set_xticks([])
ax.axvspan(37495010, 37530978, color=S.FEAT["MC"], alpha=0.07, lw=0); S.panel(ax, "", x=0, y=1)
# genes
ax = fig.add_subplot(gs[1]); ax.set_xlim(L0, L1); ax.set_ylim(0, 1); ax.axis("off")
gj = json.load(open(os.path.join(U, "genes.json"))); gk = [k for k in gj if isinstance(gj[k], list)][0]
tx = [r for r in gj[gk] if r["geneName"] == "ITGA9" and r["chromEnd"] == 37865113][0]
ax.plot([L0, L1], [0.5, 0.5], color=S.INK, lw=0.6)
for s0, sz in zip(map(int, tx["chromStarts"].strip(",").split(",")), map(int, tx["blockSizes"].strip(",").split(","))):
    a0 = tx["chromStart"] + s0
    if a0 < L1: ax.add_patch(Rectangle((a0, 0.25), max(sz, 250), 0.5, color=S.INK))
ax.text(L0, 0.95, "ITGA9 (→ continues to 37.87 Mb; TSS 37,493,580)", fontsize=6, va="top", style="italic")
ax.text(-0.01, 0.5, "GENCODE", transform=ax.transAxes, ha="right", va="center", fontsize=6)
# GeneHancer
ax = fig.add_subplot(gs[2]); ax.set_xlim(L0, L1); ax.set_ylim(0, 1); ax.axis("off")
for r in json.load(open(os.path.join(U, "geneHancerRegElementsDoubleElite.json")))["geneHancerRegElementsDoubleElite"]:
    if r["chromEnd"] > L0 and r["chromStart"] < L1:
        ax.add_patch(Rectangle((r["chromStart"], 0.2), r["chromEnd"] - r["chromStart"] + 150, 0.6, color="#E3912B" if "Promoter" in r["elementType"] else "#F2C38A"))
ax.text(-0.01, 0.5, "GeneHancer", transform=ax.transAxes, ha="right", va="center", fontsize=6)
# DNase clusters
ax = fig.add_subplot(gs[3]); ax.set_xlim(L0, L1); ax.set_ylim(0, 1); ax.axis("off")
for r in json.load(open(os.path.join(U, "wgEncodeRegDnaseClusteredV3.json")))["wgEncodeRegDnaseClusteredV3"]:
    if r["chromEnd"] > L0 and r["chromStart"] < L1:
        ax.add_patch(Rectangle((r["chromStart"], 0.2), r["chromEnd"] - r["chromStart"] + 150, 0.6, color="black", alpha=min(1, r["sourceCount"] / 60)))
ax.text(-0.01, 0.5, "ENCODE DNase\n(125 cell types)", transform=ax.transAxes, ha="right", va="center", fontsize=6)
# Roadmap chromHMM
ax = fig.add_subplot(gs[4]); ax.set_xlim(L0, L1); ax.set_ylim(-0.5, len(EPI) - 0.5); ax.invert_yaxis()
for i, (e, lab) in enumerate(EPI):
    for _, r in rm[rm.eid == e].iterrows():
        if r.end > L0 and r.start < L1:
            ax.add_patch(Rectangle((r.start, i - 0.4), r.end - r.start, 0.8, color=tuple(int(v) / 255 for v in r.color.split(",")), lw=0))
ax.set_yticks(range(len(EPI))); ax.set_yticklabels([f"{lab} ({e})" for e, lab in EPI], fontsize=5.8); ax.tick_params(axis="y", length=0)
for yy in (5.5, 9.5): ax.axhline(yy, color="white", lw=1.2)
ax.set_xlabel("Chromosome 3 position (Mb, GRCh37)"); ax.set_xticks(np.arange(37.47e6, 37.561e6, 0.01e6)); ax.set_xticklabels([f"{v / 1e6:.2f}" for v in np.arange(37.47e6, 37.561e6, 0.01e6)])
for a in fig.axes: a.axvline(SNP, color=S.FEAT["MC"], lw=0.7, ls=":", zorder=5)
states = [("1_TssA", "255,0,0", "Active TSS"), ("2_TssAFlnk", "255,69,0", "Flanking active TSS"), ("4_Tx", "0,128,0", "Transcription"), ("5_TxWk", "0,100,0", "Weak transcription"),
          ("7_Enh", "255,255,0", "Enhancer"), ("6_EnhG", "194,225,5", "Genic enhancer"), ("10_TssBiv", "205,92,92", "Bivalent TSS"), ("11_BivFlnk", "233,150,122", "Flanking bivalent"),
          ("12_EnhBiv", "189,183,107", "Bivalent enhancer"), ("3_TxFlnk", "50,205,50", "Transcr. at 5'/3'"), ("8_ZNF/Rpts", "102,205,170", "ZNF genes & repeats"), ("9_Het", "138,145,208", "Heterochromatin"), ("13_ReprPC", "128,128,128", "Polycomb repressed"), ("14_ReprPCWk", "192,192,192", "Weak Polycomb"), ("15_Quies", "255,255,255", "Quiescent")]
from matplotlib.patches import Patch
ax.legend([Patch(fc=tuple(int(v) / 255 for v in cc.split(",")), ec="#999", lw=0.3) for _, cc, _ in states], [n for _, _, n in states], ncol=5, fontsize=5.4,
          loc="upper center", bbox_to_anchor=(0.42, -0.2), handlelength=1, columnspacing=0.8)
S.save(fig, OUT, "FigS6_rs1549596_chromatin_public")
print("ok")
