"""Paper 1: MRI phenotypes, pain, disability and mental health (cross-sectional at baseline and prospective to Q2)."""
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_predict, StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from common_v2 import *
L = pd.read_csv(p(RES, "p1_long.csv")); P = pd.read_csv(p(RES, "p1_person.csv"))
yn = lambda s: s.astype(str).str.lower().map({"yes": 1, "no": 0, "true": 1, "false": 0})
for k in (1, 2):
    P[f"lbpyr{k}"] = yn(P[f"lbpyr{k}"]); P[f"sci{k}"] = yn(P[f"sciatica{k}"])
    for c in (f"mcs{k}", f"pcs{k}"): P.loc[(P[c] < 0) | (P[c] > 100), c] = np.nan     # one MCS of 437 is a data-entry error
t1 = L[L.t == 1][["IID", "age", "Core", "HIZ", "MC", "DH", "hiz_bin", "mc_bin", "dh_bin"]].rename(columns={"age": "mri_age1"})
t2 = L[L.t == 2][["IID", "age", "Core"]].rename(columns={"age": "mri_age2", "Core": "Core2"})
d = P.merge(t1, on="IID").merge(t2, on="IID", how="left")
z = lambda x: (x - x.mean()) / x.std()
for c in ["Core", "HIZ", "MC", "DH", "mcs1", "pcs1", "BMI_age40", "packyr"]: d["z" + c] = z(d[c])
d["dCore"] = z(d.Core2 - d.Core)
d["gap1"] = (d.quesage1 - d.mri_age1).abs(); d["gap2"] = (d.quesage2 - d.mri_age2).abs()
cs = d[d.gap1 <= 2].copy()          # questionnaire within 2 years of the first MRI
COV = "mri_age1 + I(mri_age1**2) + female + zBMI_age40 + zpackyr"
MRI = "zCore + zHIZ + zMC + zDH"
outs = [("lbpyr1", "logit", "LBP in past year"), ("sci1", "logit", "Sciatica ever"), ("oswestry1", "ols", "ODI"),
        ("vasToday1", "ols", "VAS today"), ("vasSeverest1", "ols", "VAS worst ever"), ("pcs1", "ols", "SF-36 PCS"), ("rmdq1", "ols", "RMDQ")]
rows, fitstat = [], []
def fit(y, kind, rhs, data):
    f = f"{y} ~ {rhs}"
    if kind == "logit": return smf.logit(f, data).fit(disp=0)
    return smf.ols(f, data).fit(cov_type="HC3")
for y, kind, lab in outs:
    dd = cs.dropna(subset=[y, "zCore", "zHIZ", "zMC", "zDH", "zmcs1", "zBMI_age40", "zpackyr", "female"])
    mods = {"A: covariates": COV, "B: + MRI": f"{COV} + {MRI}", "C: + MCS": f"{COV} + zmcs1", "D: + MRI + MCS": f"{COV} + {MRI} + zmcs1",
            "E: D + MCS x Core": f"{COV} + {MRI} + zmcs1 + zmcs1:zCore"}
    if y == "pcs1": mods = {k: v.replace(" + zmcs1:zCore", "") for k, v in mods.items()}
    res = {k: fit(y, kind, v, dd) for k, v in mods.items()}
    for k, r in res.items():
        for t in r.params.index:
            if t.startswith("z") and t not in ("zBMI_age40", "zpackyr") or ":" in t:
                b = r.params[t]; ci = r.conf_int().loc[t]
                rows.append(dict(outcome=lab, model=k, term=t, est=np.exp(b) if kind == "logit" else b,
                                 lo=np.exp(ci[0]) if kind == "logit" else ci[0], hi=np.exp(ci[1]) if kind == "logit" else ci[1],
                                 p=r.pvalues[t], scale="OR per SD" if kind == "logit" else "units per SD", n=int(r.nobs)))
        r2 = r.prsquared if kind == "logit" else r.rsquared
        fitstat.append(dict(outcome=lab, model=k, n=int(r.nobs), r2=r2, aic=r.aic))
    # block tests: MRI block given covariates + MCS; MCS given covariates + MRI
    from scipy import stats
    for blk, full, red, df in [("MRI block | cov+MCS", "D: + MRI + MCS", "C: + MCS", 4), ("MCS | cov+MRI", "D: + MRI + MCS", "B: + MRI", 1)]:
        if kind == "logit": stat = 2 * (res[full].llf - res[red].llf); pv = stats.chi2.sf(stat, df)
        else:
            rf = smf.ols(f"{y} ~ {mods[full]}", dd).fit(); rr = smf.ols(f"{y} ~ {mods[red]}", dd).fit()
            stat, pv, _ = rf.compare_f_test(rr)
        fitstat.append(dict(outcome=lab, model=blk, n=len(dd), stat=stat, p=pv))
pd.DataFrame(rows).to_csv(p(RES, "p1_pain_crosssectional.csv"), index=False)
pd.DataFrame(fitstat).to_csv(p(RES, "p1_pain_fitstats.csv"), index=False)

# prospective: Q2 outcomes on T1 MRI, change in Core, baseline outcome and MCS
pr_rows = []
lg = d[(d.gap1 <= 2) & d.quesage2.notna()].copy()
for y, base, kind, lab in [("lbpyr2", "lbpyr1", "logit", "LBP in past year at Q2"), ("oswestry2", "oswestry1", "ols", "ODI at Q2"),
                           ("vasToday2", "vasToday1", "ols", "VAS today at Q2"), ("vasSeverest2", "vasSeverest1", "ols", "VAS worst at Q2")]:
    lg["fu"] = lg.quesage2 - lg.quesage1
    for k, rhs in [("T1 MRI", f"{COV} + fu + {base} + {MRI} + zmcs1"), ("T1 MRI + change in Core", f"{COV} + fu + {base} + {MRI} + dCore + zmcs1")]:
        dd = lg.dropna(subset=[y, base, "zCore", "zHIZ", "zMC", "zDH", "zmcs1", "zBMI_age40", "zpackyr", "fu"] + (["dCore"] if "dCore" in rhs else []))
        r = fit(y, kind, rhs, dd)
        for t in ["zCore", "zHIZ", "zMC", "zDH", "dCore", "zmcs1"]:
            if t in r.params:
                b = r.params[t]; ci = r.conf_int().loc[t]
                pr_rows.append(dict(outcome=lab, model=k, term=t, est=np.exp(b) if kind == "logit" else b, lo=np.exp(ci[0]) if kind == "logit" else ci[0],
                                    hi=np.exp(ci[1]) if kind == "logit" else ci[1], p=r.pvalues[t], n=int(r.nobs)))
pd.DataFrame(pr_rows).to_csv(p(RES, "p1_pain_prospective.csv"), index=False)

# discrimination of LBP in past year at Q2: 10-fold CV AUC
dd = lg.dropna(subset=["lbpyr2", "lbpyr1", "zCore", "zHIZ", "zMC", "zDH", "zmcs1", "mri_age1", "female", "zBMI_age40", "zpackyr", "oswestry1"])
yv = dd.lbpyr2.astype(int).values; auc = []
sets = {"age+sex+BMI+smoking": ["mri_age1", "female", "zBMI_age40", "zpackyr"]}
sets["+ baseline LBP, ODI"] = sets["age+sex+BMI+smoking"] + ["lbpyr1", "oswestry1"]
sets["+ MCS"] = sets["+ baseline LBP, ODI"] + ["zmcs1"]
sets["+ MRI (Core, HIZ, MC, DH)"] = sets["+ MCS"] + ["zCore", "zHIZ", "zMC", "zDH"]
cv = StratifiedKFold(10, shuffle=True, random_state=1)
for k, cols in sets.items():
    pr = cross_val_predict(make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000)), dd[cols].values, yv, cv=cv, method="predict_proba")[:, 1]
    auc.append(dict(predictors=k, n=len(yv), cases=int(yv.sum()), cv_auc=roc_auc_score(yv, pr)))
pd.DataFrame(auc).to_csv(p(RES, "p1_pain_auc.csv"), index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 400)
R = pd.DataFrame(rows); print(R[R.model.isin(["D: + MRI + MCS", "E: D + MCS x Core"])].round(4).to_string())
print(pd.DataFrame(fitstat).round(4).to_string()); print(pd.DataFrame(pr_rows).round(4).to_string()); print(pd.DataFrame(auc).round(3))
print("cross-sectional n", len(cs), "prospective n", len(lg))
