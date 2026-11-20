"""Paper 2: rs1549596 sensitivity analyses on hard-called genotypes (genotype_bfiles/chr3)."""
import numpy as np, pandas as pd, statsmodels.api as sm, statsmodels.formula.api as smf
from scipy import stats
from common_v2 import *

g = pd.read_csv(p(UP, "Analyses_v2_2026-10-23/data/itga9_locus_genotypes.csv"))
g = g.rename(columns={"3:37495010_T_C": "rs1549596"})[["IID", "rs1549596"]]
ph = pd.read_csv(p(D8, "GWAS_phenotype.20210313.txt"), sep=" ")
e = pd.read_csv(d8("all.emergence.data.20201215.csv"), index_col=0)
lv = pd.read_csv(d8("unadjusted.levelwise.data.seven_features20210202.csv"))
pr = pd.read_csv(p(D8, "all.progression.phenotype.20201215.csv"))
cov = pd.read_csv(p(D8, "average.cov.used20201021.csv"))
d = ph.merge(g, on="IID").merge(cov[["IID", "bmi", "height", "packyr"]], on="IID", how="left")
COV = "sex + mean_age_t1t2t3 + age_sq_t1t2t3 + age_cb_t1t2t3 + PC1 + PC2 + PC3"
out = []
def add(label, res, term="rs1549596", n=None, kind="linear"):
    b, se, pv = res.params[term], res.bse[term], res.pvalues[term]
    out.append(dict(analysis=label, model=kind, n=int(n if n is not None else res.nobs), beta=b, se=se, p=pv,
                    OR=np.exp(b) if kind == "logistic" else np.nan))

# phenotype descriptives
desc = []
for t in ["Core", "MC", "HIZ", "DH"]:
    x = d[t].dropna()
    desc.append(dict(phenotype=t, n=len(x), mean=x.mean(), sd=x.std(), median=x.median(), q1=x.quantile(.25), q3=x.quantile(.75),
                     pct_le_0_05=(x <= 0.05).mean() * 100, skew=stats.skew(x), kurt=stats.kurtosis(x)))
pd.DataFrame(desc).to_csv(p(RES, "p2_phenotype_descriptives.csv"), index=False)
d[["Core", "MC", "HIZ", "DH"]].corr("spearman").to_csv(p(RES, "p2_phenotype_spearman.csv"))
d[["Core", "MC", "HIZ", "DH"]].corr("pearson").to_csv(p(RES, "p2_phenotype_pearson.csv"))

m = d.dropna(subset=["MC"]).copy()
print("MC n", len(m), "T freq", m.rs1549596.mean() / 2, "genotype counts", m.rs1549596.value_counts().to_dict())
r0 = smf.ols(f"MC ~ rs1549596 + {COV}", m).fit(); add("Primary model, hard calls (replicates GWAS)", r0)
sdmc = m.MC.std()
r2 = r0.params["rs1549596"] ** 2 * 2 * (m.rs1549596.mean() / 2) * (1 - m.rs1549596.mean() / 2) / m.MC.var()
# robust SE
rr = smf.ols(f"MC ~ rs1549596 + {COV}", m).fit(cov_type="HC3"); add("Primary model, HC3 robust SE", rr)
# rank-based inverse normal transform
m["MC_int"] = stats.norm.ppf((stats.rankdata(m.MC) - 0.5) / len(m))
add("Rank-based inverse-normal transformed MC.avg", smf.ols(f"MC_int ~ rs1549596 + {COV}", m).fit())
# + BMI, height, Core.slope
add("+ BMI, height, pack-years", smf.ols(f"MC ~ rs1549596 + {COV} + bmi + height + packyr", m).fit())
add("+ Core.slope", smf.ols(f"MC ~ rs1549596 + {COV} + Core", m).fit())
# dominant model
m["dom"] = (m.rs1549596 > 0).astype(int)
add("Dominant coding (any T)", smf.ols(f"MC ~ dom + {COV}", m).fit(), term="dom")
# sex/age strata
for s, lab in [(1, "Women"), (0, "Men")]:
    mm = m[m.sex == s]; add(f"{lab} only", smf.ols(f"MC ~ rs1549596 + mean_age_t1t2t3 + age_sq_t1t2t3 + age_cb_t1t2t3 + PC1 + PC2 + PC3", mm).fit())
med = m.mean_age_t1t2t3.median()
for cond, lab in [(m.mean_age_t1t2t3 < med, "Younger half (mean age below median)"), (m.mean_age_t1t2t3 >= med, "Older half")]:
    add(lab, smf.ols(f"MC ~ rs1549596 + {COV}", m[cond]).fit())
ia = smf.ols(f"MC ~ rs1549596 * mean_age_t1t2t3 + sex + age_sq_t1t2t3 + age_cb_t1t2t3 + PC1 + PC2 + PC3", m).fit()
add("Interaction rs1549596 x mean age", ia, term="rs1549596:mean_age_t1t2t3")
ix = smf.ols(f"MC ~ rs1549596 * sex + mean_age_t1t2t3 + age_sq_t1t2t3 + age_cb_t1t2t3 + PC1 + PC2 + PC3", m).fit()
add("Interaction rs1549596 x sex", ix, term="rs1549596:sex")

# binary MC phenotypes from the rater-harmonised emergence data
e2 = e[["IID", "mc.t1", "mc.t2", "mc.t3", "age1", "age2", "age3"]].merge(d[["IID", "rs1549596", "sex", "PC1", "PC2", "PC3", "mean_age_t1t2t3", "age_sq_t1t2t3", "age_cb_t1t2t3"]], on="IID")
e2["mc_ever"] = e2[["mc.t1", "mc.t2", "mc.t3"]].max(axis=1)
e2["lastage"] = e2[["age1", "age2", "age3"]].max(axis=1)
x = e2.dropna(subset=["mc_ever"])
add("Any MC at any scan (logistic)", smf.logit(f"mc_ever ~ rs1549596 + sex + lastage + I(lastage**2) + PC1 + PC2 + PC3", x).fit(disp=0), kind="logistic")
x1 = e2.dropna(subset=["mc.t1"])
add("Any MC at first scan, cross-sectional (logistic)", smf.logit("Q('mc.t1') ~ rs1549596 + sex + age1 + I(age1**2) + PC1 + PC2 + PC3", x1).fit(disp=0), kind="logistic")
# incident MC among MC-free at T1 with T2
xi = e2[(e2["mc.t1"] == 0) & e2["mc.t2"].notna()].copy(); xi["gap"] = xi.age2 - xi.age1
add("Incident MC T1->T2 among MC-free at T1 (logistic)", smf.logit("Q('mc.t2') ~ rs1549596 + sex + age1 + gap + PC1 + PC2 + PC3", xi).fit(disp=0), kind="logistic")
# MC unconditioned on Core: mean adjusted composite over all scans
pr2 = pr[["IID", "mc.t1", "mc.t2", "mc.t3"]].copy(); pr2["mc_allscans"] = pr2[["mc.t1", "mc.t2", "mc.t3"]].mean(axis=1)
y = d.merge(pr2[["IID", "mc_allscans"]], on="IID").dropna(subset=["mc_allscans"])
add("MC averaged over all scans (not conditioned on Core)", smf.ols(f"mc_allscans ~ rs1549596 + {COV}", y).fit())
# number of levels with MC (any scan) — Poisson
L = lv[["IID"] + [f"mc{k}.{t}" for k in range(1, 6) for t in (1, 2, 3)]].copy()
for k in range(1, 6): L[f"mcL{k}"] = L[[f"mc{k}.{t}" for t in (1, 2, 3)]].max(axis=1)
L["nlev"] = L[[f"mcL{k}" for k in range(1, 6)]].sum(axis=1, min_count=5)
z = e2.merge(L[["IID", "nlev"] + [f"mcL{k}" for k in range(1, 6)]], on="IID").dropna(subset=["nlev"])
add("Number of levels with MC, any scan (Poisson)", smf.glm("nlev ~ rs1549596 + sex + lastage + I(lastage**2) + PC1 + PC2 + PC3", z, family=sm.families.Poisson()).fit(), kind="poisson")
add("Number of levels with MC, any scan (Poisson, robust SE)", smf.glm("nlev ~ rs1549596 + sex + lastage + I(lastage**2) + PC1 + PC2 + PC3", z, family=sm.families.Poisson()).fit(cov_type="HC0"), kind="poisson")
lev = []
for k in range(1, 6):
    zz = z.dropna(subset=[f"mcL{k}"])
    r = smf.logit(f"mcL{k} ~ rs1549596 + sex + lastage + I(lastage**2) + PC1 + PC2 + PC3", zz).fit(disp=0)
    lev.append(dict(level=f"L{k}", n=int(r.nobs), cases=int(zz[f"mcL{k}"].sum()), OR=np.exp(r.params.rs1549596),
                    lo=np.exp(r.conf_int().loc["rs1549596", 0]), hi=np.exp(r.conf_int().loc["rs1549596", 1]), p=r.pvalues.rs1549596))
pd.DataFrame(lev).to_csv(p(RES, "p2_rs1549596_by_level.csv"), index=False)
# genotype means
gm = m.groupby("rs1549596").MC.agg(["count", "mean", "std", "median"]); gm.to_csv(p(RES, "p2_rs1549596_genotype_means.csv"))
# leave-one-scanner? not available. Other phenotypes with hard calls
for t in ["HIZ", "DH"]:
    mm = d.dropna(subset=[t]); add(f"{t}.avg (specificity)", smf.ols(f"{t} ~ rs1549596 + {COV}", mm).fit())
mm = d.dropna(subset=["Core"]); add("Core.slope (specificity)", smf.ols("Core ~ rs1549596 + sex + PC1 + PC2 + PC3", mm).fit())
res = pd.DataFrame(out); res.to_csv(p(RES, "p2_rs1549596_sensitivity.csv"), index=False)
pd.set_option("display.width", 250); print(res.to_string()); print(gm); print(pd.DataFrame(lev))
print("SD MC", sdmc, "beta/SD (GWAS beta .144)", 0.144251 / sdmc, "R2 hardcall", r2)
# winner's curse: conditional MLE (Zhong & Prentice 2008) for z=5.7606, threshold 5e-8
zobs = 5.7606386377432; c = stats.norm.isf(2.5e-8)
def negll(mu): return -(stats.norm.logpdf(zobs - mu) - np.log(stats.norm.sf(c - mu) + stats.norm.cdf(-c - mu)))
from scipy.optimize import minimize_scalar
mu = minimize_scalar(negll, bounds=(-10, 10), method="bounded").x
print("winners-curse corrected z", mu, "beta", mu * 0.0250408)
with open(p(RES, "p2_misc.txt"), "w") as f:
    f.write(f"SD_MC {sdmc}\nbeta_SD {0.144251/sdmc}\nR2_hardcall {r2}\nR2_from_z {zobs**2/(zobs**2+1834-10)}\nWC_corrected_z {mu}\nWC_corrected_beta {mu*0.0250408}\n")
