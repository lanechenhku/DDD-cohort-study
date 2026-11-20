"""Paper 2: permutation p for rs1549596 ~ MC.avg (Freedman–Lane, 1e6 perms) and robust count models."""
import numpy as np, pandas as pd, statsmodels.api as sm, statsmodels.formula.api as smf
from common_v2 import *
g = pd.read_csv(p(UP, "Analyses_v2_2026-10-23/data/itga9_locus_genotypes.csv")).rename(columns={"3:37495010_T_C": "g"})[["IID", "g"]]
ph = pd.read_csv(p(D8, "GWAS_phenotype.20210313.txt"), sep=" ").merge(g, on="IID").dropna(subset=["MC"])
C = sm.add_constant(ph[["sex", "mean_age_t1t2t3", "age_sq_t1t2t3", "age_cb_t1t2t3", "PC1", "PC2", "PC3"]]).values
H = C @ np.linalg.pinv(C)
ry = ph.MC.values - H @ ph.MC.values
rg = ph.g.values - H @ ph.g.values
obs = abs(rg @ ry)
rng = np.random.default_rng(20261023); B = 1_000_000; hit = 0
for _ in range(B // 50_000):
    P = np.argsort(rng.random((50_000, len(ry))), axis=1)
    hit += (np.abs(ry[P] @ rg) >= obs).sum()
print("perm p", (hit + 1) / (B + 1), "hits", hit)
open(p(RES, "p2_misc.txt"), "a").write(f"perm_p_1e6 {(hit+1)/(B+1)} hits {hit}\n")
