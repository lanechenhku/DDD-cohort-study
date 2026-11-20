"""Paper 2 new Figure 2 (v3, 2026-11-06): descriptive data for the rs1549596 x MC figure.
Any-MC prevalence (any scan, rater-harmonised emergence data, 69 undated T2 records excluded) by carrier status,
by age at last scan, by lumbar level and by sex; Wilson 95% CIs."""
import numpy as np, pandas as pd
from statsmodels.stats.proportion import proportion_confint
from common_v2 import *
g = pd.read_csv(p(UP, "Analyses_v2_2026-10-23/data/itga9_locus_genotypes.csv")).rename(columns={"3:37495010_T_C": "rs1549596"})[["IID", "rs1549596"]]
ph = pd.read_csv(p(D8, "GWAS_phenotype.20210313.txt"), sep=" ")
e = pd.read_csv(d8("all.emergence.data.20201215.csv"), index_col=0)
lv = pd.read_csv(d8("unadjusted.levelwise.data.seven_features20210202.csv"))
d = ph[["IID", "sex"]].merge(g, on="IID")
x = e[["IID", "mc.t1", "mc.t2", "mc.t3", "age1", "age2", "age3"]].merge(d, on="IID")
x["mc_ever"] = x[["mc.t1", "mc.t2", "mc.t3"]].max(axis=1); x["lastage"] = x[["age1", "age2", "age3"]].max(axis=1)
x = x.dropna(subset=["mc_ever"]); x["carrier"] = (x.rs1549596 > 0).astype(int)
L = lv[["IID"] + [f"mc{k}.{t}" for k in range(1, 6) for t in (1, 2, 3)]].copy()
for k in range(1, 6): L[f"L{k}"] = L[[f"mc{k}.{t}" for t in (1, 2, 3)]].max(axis=1)
x = x.merge(L[["IID"] + [f"L{k}" for k in range(1, 6)]], on="IID", how="left")
rows = []
def add(kind, grp, sub, col="mc_ever"):
    for c in (0, 1):
        s = sub[sub.carrier == c][col].dropna(); k, n = int(s.sum()), len(s)
        lo, hi = proportion_confint(k, n, method="wilson") if n else (np.nan, np.nan)
        rows.append(dict(kind=kind, group=grp, carrier=c, cases=k, n=n, prev=k / n if n else np.nan, lo=lo, hi=hi))
add("overall", "All", x)
bands = [(0, 40, "<40"), (40, 50, "40–49"), (50, 60, "50–59"), (60, 200, "≥60")]
for a, b, lab in bands: add("age", lab, x[(x.lastage >= a) & (x.lastage < b)])
for k, lab in zip(range(1, 6), ["L1/2", "L2/3", "L3/4", "L4/5", "L5/S1"]): add("level", lab, x, col=f"L{k}")
for s, lab in [(1, "Women"), (0, "Men")]: add("sex", lab, x[x.sex == s])
out = pd.DataFrame(rows); out.to_csv(p(RES, "p2_fig2_mc_prevalence_by_carrier.csv"), index=False)
print(out.to_string()); print("N", len(x), "carriers", x.carrier.sum())
