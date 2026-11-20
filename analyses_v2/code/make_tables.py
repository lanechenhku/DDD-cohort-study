"""Render the v2 result CSVs as Markdown tables for the supplements (results/tables_md/)."""
import numpy as np, pandas as pd, os
from common_v2 import *
from figs_v2 import TM, PG
O = p(RES, "tables_md"); os.makedirs(O, exist_ok=True)
def fp(x):
    if pd.isna(x): return "–"
    if x < 1e-15: return "<1e−15"
    return f"{x:.2g}".replace("e-0", "e−").replace("e-", "e−") if x < 0.001 else f"{x:.3f}"
def ci(r, k=2): return f"{r.est:.{k}f} ({r.lo:.{k}f}, {r.hi:.{k}f})"
def w(name, s): open(p(O, name), "w").write(s); print("==", name); print(s[:1500])

m = pd.read_csv(p(RES, "p1_wgee_model1.csv"))
x = m[m.model == "M1"].copy(); x["t"] = x.term.map(TM); x = x[x.t.notna()]
rows = ["| Term | " + " | ".join(f"{o} β (95% CI) | p" for o in ["Core", "HIZ", "MC", "DH"]) + " |", "|---|" + "---|---|" * 4]
for t in TM.values():
    cells = []
    for o in ["Core", "HIZ", "MC", "DH"]:
        r = x[(x.t == t) & (x.outcome == o)]
        cells += [ci(r.iloc[0]), fp(r.iloc[0].p)] if len(r) else ["n/a", ""]
    rows.append(f"| {t} | " + " | ".join(cells) + " |")
n = x.groupby("outcome")[["n_obs", "n_id"]].first()
w("S4_model1.md", "\n".join(rows) + "\n\nObservations / participants: " + "; ".join(f"{o} {int(n.loc[o,'n_obs'])}/{int(n.loc[o,'n_id'])}" for o in ["Core", "HIZ", "MC", "DH"]) + "\n")

keys = ["male", "z_BMI_age40", "z_Height_age40", "z_injury", "z_sport", "z_core_onset", "zage:z_BMI_age40", "zage:z_packyr"]
models = ["M1", "M1_unweighted", "M1_exchangeable", "M1_cumulative_ipw", "M1_ipw_with_MRI", "M1_probands_only", "M1_baseline_age_18_55"]
lab = {"M1": "Primary", "M1_unweighted": "Unweighted", "M1_exchangeable": "Exchangeable corr.", "M1_cumulative_ipw": "Cumulative IPW",
       "M1_ipw_with_MRI": "IPW incl. sex + MRI", "M1_probands_only": "Probands only", "M1_baseline_age_18_55": "Baseline age 18–55"}
rows = ["| Outcome | Term | " + " | ".join(lab[k] for k in models) + " |", "|---|---|" + "---|" * len(models)]
for o in ["Core", "HIZ", "MC", "DH"]:
    for k in keys:
        r = m[(m.outcome == o) & (m.term == k)].set_index("model")
        if not len(r): continue
        rows.append(f"| {o} | {TM[k]} | " + " | ".join((f"{r.loc[q,'est']:.2f}{'*' if r.loc[q,'p']<0.05 else ''}" if q in r.index else "–") for q in models) + " |")
s = m[m.model == "M1_unadjCore_scanner_covariate"].set_index("term")
w("S4b_sensitivity.md", "\n".join(rows) + "\n\n\\* p < 0.05. Unadjusted Core with scanner as a covariate: " +
  ", ".join(f"{TM[k]} {s.loc[k,'est']:.2f}" for k in ["male", "z_BMI_age40", "z_Height_age40", "z_injury"]) + ".\n")

m2 = pd.read_csv(p(RES, "p1_wgee_model2.csv"))
rows = ["| PGS | Outcome | Main β (SE) | p | Age × PGS β (SE) | p |", "|---|---|---|---|---|---|"]
for v in PG:
    for o in ["Core", "HIZ", "MC", "DH"]:
        a = m2[(m2.pgs == v) & (m2.outcome == o) & (m2.term == f"z_{v}")]; b = m2[(m2.pgs == v) & (m2.outcome == o) & (m2.term == f"zage:z_{v}")]
        if not len(a): continue
        a, b = a.iloc[0], b.iloc[0]
        rows.append(f"| {PG[v]} | {o} | {a.est:.3f} ({a.se:.3f}) | {fp(a.p)} | {b.est:.3f} ({b.se:.3f}) | {fp(b.p)} |")
w("S5_model2.md", "\n".join(rows) + "\n")
m3 = pd.read_csv(p(RES, "p1_wgee_model3.csv")); m3 = m3[m3.term.str.contains("PRS")]
def pl(t):
    for k, v in PG.items():
        if k in t: return ("Age × " if t.startswith("zage:") else "") + v
rows = ["| Term | Outcome | Joint model β (95% CI) | p | + measured BMI/height β (95% CI) | p |", "|---|---|---|---|---|---|"]
for o in ["Core", "HIZ", "MC", "DH"]:
    for t in m3[(m3.model == "M3") & (m3.outcome == o)].term:
        a = m3[(m3.model == "M3") & (m3.outcome == o) & (m3.term == t)].iloc[0]; b = m3[(m3.model != "M3") & (m3.outcome == o) & (m3.term == t)].iloc[0]
        if min(a.p, b.p) < 0.05: rows.append(f"| {pl(t)} | {o} | {ci(a,3)} | {fp(a.p)} | {ci(b,3)} | {fp(b.p)} |")
w("S6_model3.md", "\n".join(rows) + "\n\nOnly terms with p < 0.05 in either model are shown; full output in results/p1_wgee_model3.csv.\n")

tr = pd.read_csv(p(RES, "p1_transitions.csv"))
rows = ["| Interval | Feature | n | 1→0 | 0→0 | 1→1 | 0→1 | Regression among affected, % (95% CI) | Regression among all, % | Incidence among unaffected, % (95% CI) |", "|---|---|---|---|---|---|---|---|---|---|"]
for _, r in tr.iterrows():
    rows.append(f"| {r.interval.replace('->','→')} | {r.feature} | {r.n:,} | {r.n10} | {r.n00:,} | {r.n11:,} | {r.n01} | {r.regr_affected:.1f} ({r.regr_aff_lo:.1f}–{r.regr_aff_hi:.1f}) | {r.regr_all:.1f} | {r.incidence:.1f} ({r.inc_lo:.1f}–{r.inc_hi:.1f}) |")
w("T2_transitions.md", "\n".join(rows) + "\n")

ir = pd.read_csv(p(RES, "p1_interrater_summary.csv"))
rows = ["| Feature | Scan | Rater pairs | Disc-level paired ratings | Exact agreement, % | κ median (range) |", "|---|---|---|---|---|---|"]
for _, r in ir.iterrows(): rows.append(f"| {r.feature} | {r.scan} | {int(r.pairs)} | {int(r.n):,} | {r.pct_agree:.1f} | {r.kappa_median:.2f} ({r.kappa_min:.2f}–{r.kappa_max:.2f}) |")
w("S1_interrater.md", "\n".join(rows) + "\n")

at = pd.read_csv(p(RES, "p1_attrition_logit.csv")); at = at[(at.model == "+ LBP, ODI, BMI, smoking") & (at.term != "Intercept")]
nm = {"age": "Age (per year)", "I(age ** 2)": "Age²", "female": "Female", "zCore": "Core composite (per SD)", "anyHIZ": "Any HIZ", "anyMC": "Any MC", "anyDH": "Any DH",
      "lbp1n": "LBP ever", "oswestry1": "ODI (per point)", "BMI_age40": "BMI at 40 (per kg/m²)", "packyr": "Pack-years"}
rows = ["| Baseline characteristic | OR (95% CI) for attending T2 | p |", "|---|---|---|"] + [f"| {nm[r.term]} | {r.OR:.3f} ({r.lo:.3f}–{r.hi:.3f}) | {fp(r.p)} |" for _, r in at.iterrows()]
ag = pd.read_csv(p(RES, "p1_attrition_by_agegroup.csv"))
w("S8_attrition.md", "\n".join(rows) + f"\n\nn = {int(at.n.iloc[0]):,}. Attendance by baseline age group: " + "; ".join(f"{r.agegrp} {r.pct:.1f}% ({int(r['sum'])}/{int(r['count'])})" for _, r in ag.iterrows()) + ".\n")

sc = pd.read_csv(p(RES, "p1_scanner_effects.csv")); sc = sc[sc.outcome != "count"]
rows = ["| Scan | Scanner (vs reference) | Core offset (SE) | p | n |", "|---|---|---|---|---|"] + [f"| T{int(r.t)} | {r.term.split('T.')[1].rstrip(']').replace('.0','')} | {r.est:.2f} ({r.se:.2f}) | {fp(r.p)} | {int(r.n)} |" for _, r in sc.iterrows()]
w("S2b_scanner.md", "\n".join(rows) + "\n")

s2 = pd.read_csv(p(RES, "p1_tableS2_levelwise.csv"))
rows = ["| Scan | Feature | n | L1/2 | L2/3 | L3/4 | L4/5 | L5/S1 |", "|---|---|---|---|---|---|---|---|"] + [f"| {r.scan} | {r.feature} | {int(r.n):,} | " + " | ".join(f"{r[f'L{k}']:.1f}" for k in range(1, 6)) + " |" for _, r in s2.iterrows()]
w("S2_levelwise.md", "\n".join(rows) + "\n")

c = pd.read_csv(p(RES, "p1_pain_crosssectional.csv")); c = c[c.model == "D: + MRI + MCS"]
tn = {"zCore": "Core", "zHIZ": "HIZ", "zMC": "MC", "zDH": "DH", "zmcs1": "SF-36 MCS"}
outs = ["LBP in past year", "Sciatica ever", "ODI", "RMDQ", "VAS today", "VAS worst ever", "SF-36 PCS"]
rows = ["| Predictor (per SD) | " + " | ".join(outs) + " |", "|---|" + "---|" * len(outs)]
for t, nm_ in tn.items():
    cells = []
    for o in outs:
        r = c[(c.outcome == o) & (c.term == t)].iloc[0]; k = 2 if r.scale.startswith("OR") else 1
        cells.append(f"{r.est:.{k}f} ({r.lo:.{k}f}, {r.hi:.{k}f}){'*' if r.p < 0.05 else ''}")
    rows.append(f"| {nm_} | " + " | ".join(cells) + " |")
fs = pd.read_csv(p(RES, "p1_pain_fitstats.csv"))
r2 = lambda o, mm: fs[(fs.outcome == o) & (fs.model == mm)].r2.iloc[0]
rows.append("| R² (covariates → + MRI → + MRI + MCS) | " + " | ".join(f"{r2(o,'A: covariates'):.3f} → {r2(o,'B: + MRI'):.3f} → {r2(o,'D: + MRI + MCS'):.3f}" for o in outs) + " |")
w("T3_pain.md", "\n".join(rows) + f"\n\nn = {int(c.n.max()):,}.\n")
pr = pd.read_csv(p(RES, "p1_pain_prospective.csv")); pr = pr[pr.model == "T1 MRI"]
rows = ["| Outcome at second questionnaire | " + " | ".join(tn.values()) + " | n |", "|---|" + "---|" * 6]
for o in pr.outcome.unique():
    cells = []
    for t in tn:
        r = pr[(pr.outcome == o) & (pr.term == t)].iloc[0]; k = 2 if "LBP" in o else 1
        cells.append(f"{r.est:.{k}f} ({r.lo:.{k}f}, {r.hi:.{k}f}){'*' if r.p < 0.05 else ''}")
    rows.append(f"| {o} | " + " | ".join(cells) + f" | {int(r.n)} |")
au = pd.read_csv(p(RES, "p1_pain_auc.csv"))
w("S10_pain_prospective.md", "\n".join(rows) + "\n\nCross-validated AUC for LBP in the past year at Q2: " + "; ".join(f"{r.predictors} {r.cv_auc:.3f}" for _, r in au.iterrows()) + f" (n = {au.n.iloc[0]}, {au.cases.iloc[0]} cases).\n")

d = pd.read_csv(p(RES, "p2_rs1549596_sensitivity.csv"))
rows = ["| Analysis | Model | n | β or log-OR (SE) | OR | p |", "|---|---|---|---|---|---|"] + [f"| {r.analysis} | {r.model} | {r.n:,} | {r.beta:.3f} ({r.se:.3f}) | {'' if pd.isna(r.OR) else f'{r.OR:.2f}'} | {fp(r.p)} |" for _, r in d.iterrows()]
w("P2_S4_sensitivity.md", "\n".join(rows) + "\n")
lv = pd.read_csv(p(RES, "p2_rs1549596_by_level.csv"))
w("P2_S5_level.md", "\n".join(["| Level | Cases | n | OR per T allele (95% CI) | p |", "|---|---|---|---|---|"] + [f"| {r.level} | {r.cases} | {r.n:,} | {r.OR:.2f} ({r.lo:.2f}–{r.hi:.2f}) | {fp(r.p)} |" for _, r in lv.iterrows()]) + "\n")
de = pd.read_csv(p(RES, "p2_phenotype_descriptives.csv")); sp = pd.read_csv(p(RES, "p2_phenotype_spearman.csv"), index_col=0)
rows = ["| Phenotype | n | Mean (SD) | Median (IQR) | % ≤ 0.05 | Skewness |", "|---|---|---|---|---|---|"] + [f"| {r.phenotype}{'.slope' if r.phenotype=='Core' else '.avg'} | {r.n:,} | {r['mean']:.2f} ({r.sd:.2f}) | {r['median']:.2f} ({r.q1:.2f}–{r.q3:.2f}) | {r.pct_le_0_05:.1f} | {r['skew']:.2f} |" for _, r in de.iterrows()]
w("P2_S2_descriptives.md", "\n".join(rows) + "\n\nSpearman ρ: " + "; ".join(f"{a}–{b} {sp.loc[a,b]:.2f}" for a, b in [("Core", "MC"), ("Core", "HIZ"), ("Core", "DH"), ("MC", "HIZ"), ("MC", "DH"), ("HIZ", "DH")]) + ".\n")
