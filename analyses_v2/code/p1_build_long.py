"""Paper 1: build the long (person x scan) analysis file used by the WGEE, attrition and pain models."""
import numpy as np, pandas as pd
from common_v2 import *

a = pd.read_csv(p(D8, "allcovariate.20200917.csv"), index_col=0); a["IID"] = a.index
cv = pd.read_csv(p(D8, "average.cov.used20201021.csv"))
core = pd.read_csv(p(D8, "core.machine_adjusted.20210202.csv"))
coreu = pd.read_csv(d8("core.unadjusted.20210202.csv"))
pr = pd.read_csv(p(D8, "all.progression.phenotype.20201215.csv"))
e = pd.read_csv(d8("all.emergence.data.20201215.csv"), index_col=0)
bmi = pd.read_csv(p(D8, "predictedBMI-20201129.txt"), sep="\t")[["DNA.No", "BMI_age40"]]
hgt = pd.read_csv(p(D8, "predictedHeight-20201129.txt"), sep="\t")[["DNA.No", "Height_age40"]]
prs = pd.read_csv(p(D4, "ALL.PRS.txt"), sep=" ").drop(columns=["BMI"])
import subprocess, io
on = pd.read_csv(io.StringIO(subprocess.run(["Rscript", "-e", f'x<-readRDS("{p(D8, "pred.onset8.RDS")}"); write.csv(x[,c("IID","SS","DB","SS_CensorType","DB_CensorType")], stdout(), row.names=FALSE)'], capture_output=True, text=True).stdout))
on["core_onset"] = on[["SS", "DB"]].min(axis=1)

base = a[["IID", "DNA.No", "sex", "proband", "cd1", "cd2", "cd3", "lbp1", "lbp2", "sciatica1", "sciatica2", "oswestry1", "oswestry2",
          "vasToday1", "vasToday2", "vasSeverest1", "vasSeverest2", "quesage1", "quesage2", "smoke1"]].copy()
base["female"] = (base.sex == "F").astype(float); base.loc[base.sex.isna(), "female"] = np.nan
base = base.merge(cv[["IID", "packyr", "sport", "injury", "workload"]], on="IID", how="left")
base = base.merge(bmi, on="DNA.No", how="left").merge(hgt, on="DNA.No", how="left").merge(prs, on="DNA.No", how="left")
base = base.merge(on[["IID", "core_onset", "SS", "DB"]], on="IID", how="left")
# SF-36 summary scores and LBP-in-past-year from the questionnaires (keyed by DNA No)
for k, f in [(1, "Questionnaire (1st).xlsx"), (2, "Questionnaire (2nd).xlsx")]:
    q = pd.read_excel(p(D4, f))
    q = q.rename(columns={"DNA No": "DNA.No", "Mental  Component (MCS)": f"mcs{k}", "Physical Component (PCS)": f"pcs{k}",
                          "LBP in this year": f"lbpyr{k}", "Roland Morris Disability total score": f"rmdq{k}"})
    base = base.merge(q[["DNA.No", f"mcs{k}", f"pcs{k}", f"lbpyr{k}", f"rmdq{k}"]].drop_duplicates("DNA.No"), on="DNA.No", how="left")
rows = []
for t in (1, 2, 3):
    d = pd.DataFrame({"IID": e.IID, "t": t, "age": e[f"age{t}"]})
    d = d.merge(core[["IID", f"core.t{t}.adj"]].rename(columns={f"core.t{t}.adj": "Core"}), on="IID", how="left")
    d = d.merge(coreu[["IID", f"core.t{t}"]].rename(columns={f"core.t{t}": "Core_unadj"}), on="IID", how="left")
    d = d.merge(pr[["IID", f"hiz.t{t}", f"mc.t{t}", f"dh.t{t}"]].rename(columns={f"hiz.t{t}": "HIZ", f"mc.t{t}": "MC", f"dh.t{t}": "DH"}), on="IID", how="left")
    d = d.merge(e[["IID", f"ss.t{t}", f"db.t{t}", f"hiz.t{t}", f"mc.t{t}", f"dh.t{t}"]].rename(columns=lambda c: c.replace(f".t{t}", "_bin")), on="IID", how="left")
    rows.append(d)
L = pd.concat(rows).merge(base, on="IID", how="left")
L["scanner"] = np.select([L.t == 1, L.t == 2, L.t == 3], [L.cd1, L.cd2, L.cd3])
L = L[L.age.notna()].sort_values(["IID", "t"])
# attendance weights: P(scan t | scan t-1) from logistic regression on age, age^2 at scan t-1
import statsmodels.formula.api as smf
w = pd.DataFrame({"IID": e.IID, "age1": e.age1, "age2": e.age2, "age3": e.age3})
w["att2"] = w.age2.notna().astype(int); w["att3"] = w.age3.notna().astype(int)
m2 = smf.logit("att2 ~ age1 + I(age1**2)", w[w.age1.notna()]).fit(disp=0)
m3 = smf.logit("att3 ~ age2 + I(age2**2)", w[w.age2.notna()]).fit(disp=0)
w.loc[w.age1.notna(), "p2"] = m2.predict(w[w.age1.notna()]); w.loc[w.age2.notna(), "p3"] = m3.predict(w[w.age2.notna()])
L = L.merge(w[["IID", "p2", "p3"]], on="IID", how="left")
L["ipw"] = np.select([L.t == 1, L.t == 2, L.t == 3], [1.0, 1 / L.p2, 1 / L.p3])
L["ipw_cum"] = np.select([L.t == 1, L.t == 2, L.t == 3], [1.0, 1 / L.p2, 1 / (L.p2 * L.p3)])
L.to_csv(p(RES, "p1_long.csv"), index=False)
base.to_csv(p(RES, "p1_person.csv"), index=False)
print(L.shape, L.groupby("t").size().to_dict())
print(L[["Core", "HIZ", "MC", "DH", "age", "packyr", "BMI_age40", "Height_age40", "core_onset", "mcs1"]].describe().T[["count", "mean", "std"]])
