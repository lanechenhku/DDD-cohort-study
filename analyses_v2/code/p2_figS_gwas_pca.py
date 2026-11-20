"""Paper 2 Supplementary Figure S1b: genotype principal components of the 2,121 participants who passed genotype QC
(2026-11-07). Input: post-QC genotype files on Lane's Mac (DDD 2020-04/genotype_bfiles, chr1-22), thinned on the Mac to
bi-allelic, non-palindromic SNPs at least 50 kb apart (extract script recorded in the header of results/p2_gwas_pca_summary.txt);
the packed genotypes were copied to /tmp only (individual-level data are not kept in the project folder).
Long-range LD regions (MHC, chr8 inversion, chr17q21) are removed, SNPs with MAF < 0.05 or missingness > 2% are dropped,
genotypes are standardised and PCs computed by SVD. Pre-QC data (2,482 samples) are not on the Mac, so the PCA shows the
analysed sample only."""
import os, sys, numpy as np, pandas as pd, matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__)); import style as S
from common_v2 import RES
S.setup(); OUT = os.path.join(RES, "..", "..", "figures_v2", "Paper2")
SRC = sys.argv[1] if len(sys.argv) > 1 else "/tmp/claude-0/pca"
P = np.load(f"{SRC}/geno_packed.npy"); snps = pd.read_csv(f"{SRC}/snps.txt", sep="\t", header=None, names=["chr", "id", "cm", "bp", "a1", "a2"])
fam = pd.read_csv(f"{SRC}/samples.fam", sep=r"\s+", header=None); n = len(fam)
# decode PLINK .bed (SNP-major): 00 hom A1 -> 2 copies of A1, 01 missing, 10 het, 11 hom A2
lut = np.array([2, np.nan, 1, 0], dtype=float)
G = np.empty((P.shape[0], P.shape[1] * 4))
for k in range(4): G[:, k::4] = lut[(P >> (2 * k)) & 3]
G = G[:, :n]
lr = [(6, 25e6, 34e6), (8, 8e6, 12e6), (17, 40e6, 45e6), (5, 44e6, 51.5e6), (11, 45e6, 57e6)]
keep = np.ones(len(snps), bool)
for c, a, b in lr: keep &= ~((snps.chr == c) & (snps.bp >= a) & (snps.bp <= b))
miss = np.isnan(G).mean(1); p = np.nanmean(G, 1) / 2; maf = np.minimum(p, 1 - p)
keep &= (miss <= 0.02) & (maf >= 0.05)
G = G[keep]; p = p[keep]
Z = (G - 2 * p[:, None]) / np.sqrt(2 * p * (1 - p))[:, None]; Z[np.isnan(Z)] = 0
U, s, Vt = np.linalg.svd(Z.T, full_matrices=False)
pcs = U[:, :10] * s[:10]; ve = s ** 2 / (s ** 2).sum()
pcdf = pd.DataFrame(pcs, columns=[f"PC{i + 1}" for i in range(10)])
sd = (pcdf - pcdf.mean()) / pcdf.std()
out6 = (sd.iloc[:, :4].abs() > 6).any(axis=1)
summ = [f"SNPs used: {keep.sum():,} of {len(snps):,} thinned SNPs (>=50 kb apart; MAF>=0.05; missing<=2%; long-range LD regions removed)",
        f"Samples: {n:,}", "Variance explained PC1-10 (%): " + ", ".join(f"{v * 100:.2f}" for v in ve[:10]),
        f"Samples beyond 6 SD on PC1-4: {int(out6.sum())}"]
open(os.path.join(RES, "p2_gwas_pca_summary.txt"), "w").write("\n".join(summ) + "\n"); print("\n".join(summ))
fig, axs = plt.subplots(1, 3, figsize=(180 / 25.4, 62 / 25.4), gridspec_kw=dict(wspace=0.38))
for ax, (a, b), L in zip(axs[:2], [(0, 1), (2, 3)], "AB"):
    ax.scatter(pcs[~out6, a], pcs[~out6, b], s=2, lw=0, color=S.PAL[0] if hasattr(S, "PAL") else "#3b6ea5", alpha=0.6, rasterized=True)
    if out6.any(): ax.scatter(pcs[out6, a], pcs[out6, b], s=6, lw=0, color="#c0392b")
    ax.set_xlabel(f"PC{a + 1} ({ve[a] * 100:.2f}%)"); ax.set_ylabel(f"PC{b + 1} ({ve[b] * 100:.2f}%)")
    S.panel(ax, L, x=-0.2, y=1.04)
ax = axs[2]; ax.bar(range(1, 11), ve[:10] * 100, color="#7f8c8d"); ax.set_xticks(range(1, 11))
ax.set_xlabel("Principal component"); ax.set_ylabel("Variance explained (%)"); S.panel(ax, "C", x=-0.2, y=1.04)
for e in ("png", "pdf"): fig.savefig(os.path.join(OUT, f"FigS1b_genotype_PCA.{e}"), dpi=300, bbox_inches="tight")
