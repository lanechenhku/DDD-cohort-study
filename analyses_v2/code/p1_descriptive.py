"""Paper 1: attrition, regression rates among affected, scanner effects, level-wise summary (Table S2), Turnbull onset."""
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy import stats
from common_v2 import *
e = pd.read_csv(d8("all.emergence.data.20201215.csv"), index_col=0)
L = pd.read_csv(p(RES, "p1_long.csv")); P = pd.read_csv(p(RES, "p1_person.csv"))
F = ["ss", "db", "dn", "dh", "hiz", "mc", "sn"]
def wilson(k, n, z=1.96):
    if n == 0: return np.nan, np.nan, np.nan
    ph = k / n; den = 1 + z*z/n; c = (ph + z*z/(2*n)) / den; h = z*np.sqrt(ph*(1-ph)/n + z*z/(4*n*n)) / den
    return ph, c - h, c + h

# ---- 1. attrition: attendance at T2 among T1 participants
b = L[L.t == 1].copy()
b["att2"] = b.IID.isin(L[L.t == 2].IID).astype(int)
b["anyHIZ"] = b.hiz_bin; b["anyMC"] = b.mc_bin; b["anyDH"] = b.dh_bin
b["agegrp"] = pd.cut(b.age, [0, 20, 30, 40, 50, 60, 70, 200], right=False, labels=["<20", "20-29", "30-39", "40-49", "50-59", "60-69", ">=70"])
tab = b.groupby("agegrp", observed=True).att2.agg(["sum", "count"]); tab["pct"] = 100 * tab["sum"] / tab["count"]
tab.to_csv(p(RES, "p1_attrition_by_agegroup.csv"))
b["lbp1n"] = b.lbp1.map({True: 1, False: 0, "True": 1, "False": 0})
b["zCore"] = (b.Core - b.Core.mean()) / b.Core.std()
rows = []
r0 = smf.logit("att2 ~ age + I(age**2)", b.dropna(subset=["zCore"])).fit(disp=0)
r1 = smf.logit("att2 ~ age + I(age**2) + female + zCore + anyHIZ + anyMC + anyDH", b).fit(disp=0)
r2 = smf.logit("att2 ~ age + I(age**2) + female + zCore + anyHIZ + anyMC + anyDH + lbp1n + oswestry1 + BMI_age40 + packyr", b).fit(disp=0)
for nm, r in [("age only", r0), ("+ sex + baseline MRI", r1), ("+ LBP, ODI, BMI, smoking", r2)]:
    for t in r.params.index:
        rows.append(dict(model=nm, term=t, OR=np.exp(r.params[t]), lo=np.exp(r.conf_int().loc[t, 0]), hi=np.exp(r.conf_int().loc[t, 1]), p=r.pvalues[t], n=int(r.nobs)))
# LR test of MRI block given age
r0b = smf.logit("att2 ~ age + I(age**2) + female", b.dropna(subset=["zCore", "anyHIZ", "anyMC", "anyDH"])).fit(disp=0)
lr = 2 * (r1.llf - r0b.llf); lrp = stats.chi2.sf(lr, 4)
att = pd.DataFrame(rows); att.to_csv(p(RES, "p1_attrition_logit.csv"), index=False)
open(p(RES, "p1_misc.txt"), "w").write(f"attrition_LR_MRIblock chi2={lr:.2f} df=4 p={lrp:.3g} n={int(r1.nobs)}\n")
print(tab); print(att[att.model != "age only"].round(3).to_string()); print("LR MRI block", lr, lrp)

# ---- 2. transitions with regression rate among affected
tr = []
for a_, b_ in [(1, 2), (2, 3)]:
    for f in F:
        x = e[[f"{f}.t{a_}", f"{f}.t{b_}"]].dropna()
        n10 = int(((x.iloc[:, 0] == 1) & (x.iloc[:, 1] == 0)).sum()); n11 = int(((x.iloc[:, 0] == 1) & (x.iloc[:, 1] == 1)).sum())
        n00 = int(((x.iloc[:, 0] == 0) & (x.iloc[:, 1] == 0)).sum()); n01 = int(((x.iloc[:, 0] == 0) & (x.iloc[:, 1] == 1)).sum())
        pa, la, ha = wilson(n10, n10 + n11); pi, li, hi = wilson(n01, n01 + n00); pt, lt, ht = wilson(n10, len(x))
        tr.append(dict(interval=f"T{a_}->T{b_}", feature=f.upper(), n=len(x), n10=n10, n00=n00, n11=n11, n01=n01,
                       regr_all=100*pt, regr_all_lo=100*lt, regr_all_hi=100*ht,
                       regr_affected=100*pa, regr_aff_lo=100*la, regr_aff_hi=100*ha, incidence=100*pi, inc_lo=100*li, inc_hi=100*hi))
tr = pd.DataFrame(tr); tr.to_csv(p(RES, "p1_transitions.csv"), index=False); print(tr.round(1).to_string())

# ---- 3. scanner distribution and effects (unadjusted composites at each scan)
sc = []
for t in (1, 2, 3):
    d = L[L.t == t].copy()
    for y in ["Core_unadj"]:
        r = smf.ols(f"{y} ~ C(scanner) + age + female + BMI_age40 + Height_age40", d).fit()
        for k in r.params.index:
            if "scanner" in k: sc.append(dict(t=t, outcome=y, term=k, est=r.params[k], se=r.bse[k], p=r.pvalues[k], n=int(r.nobs)))
    sc.append(dict(t=t, outcome="count", term=str(d.scanner.value_counts().sort_index().to_dict())))
pd.DataFrame(sc).to_csv(p(RES, "p1_scanner_effects.csv"), index=False); print(pd.DataFrame(sc).round(3).to_string())
# composite level after adjustment: residual association with scanner should be ~0
d1 = L[L.t == 1]; print("adj Core ~ scanner p", smf.ols("Core ~ C(scanner) + age + female", d1).fit().f_pvalue)

# ---- 4. Table S2: level-wise prevalence by scan
lv = pd.read_csv(d8("unadjusted.levelwise.data.seven_features20210202.csv"))
s2 = []
for t in (1, 2, 3):
    for f in ["ss", "db", "dn", "dh", "hiz", "mc", "sn"]:
        row = dict(scan=f"T{t}", feature=f.upper())
        for k in range(1, 6):
            col = f"{f}{k}.{t}" if f"{f}{k}.{t}" in lv else f"{f}.L{k}.t{t}"
            if col not in lv: row[f"L{k}"] = np.nan; continue
            v = lv[col].dropna(); row[f"L{k}"] = 100 * (v > 0).mean(); row["n"] = len(v)
            if f in ("ss", "db"): row[f"L{k}_mean"] = v.mean()
        s2.append(row)
s2 = pd.DataFrame(s2); s2.to_csv(p(RES, "p1_tableS2_levelwise.csv"), index=False); print(s2.round(1).to_string())
