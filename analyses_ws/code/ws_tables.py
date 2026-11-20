"""Markdown tables for the whole-spine section of Paper 1 (main Table 4; Supplementary Tables S12-S17). 2026-11-06."""
import pandas as pd, numpy as np
R = "/mnt/project-files/DDD_papers/analyses_ws/results"; T = R + "/tables_md"; import os; os.makedirs(T, exist_ok=True)
NAME = {'PS': 'Pfirrmann grade', 'SS': 'Schneiderman score*', 'PDB': 'Posterior bulge', 'ADB': 'Anterior bulge', 'HIZ': 'HIZ', 'OMC': 'Modic change (overall)', 'AMC': 'Modic change (anterior)', 'SN': "Schmorl's node"}
FE = ['PS', 'PDB', 'ADB', 'HIZ', 'OMC', 'AMC', 'SN']
lab = ['C2-3','C3-4','C4-5','C5-6','C6-7','C7-T1'] + [f'T{i}-{i+1}' for i in range(1, 12)] + ['T12-L1'] + [f'L{i}-{i+1}' for i in range(1, 5)] + ['L5-S1']
def md(df): return '| ' + ' | '.join(df.columns) + ' |\n|' + '|'.join(['---'] * len(df.columns)) + '|\n' + '\n'.join('| ' + ' | '.join(str(x) for x in r) + ' |' for r in df.values) + '\n'
f3 = lambda x: f'{x:.2f}'
pv = lambda p: ('<0.001' if p < 0.001 else f'{p:.3f}')
p = pd.read_csv(f"{R}/ws_prevalence_by_level.csv"); c = pd.read_csv(f"{R}/ws_corr_by_lag.csv"); rc = pd.read_csv(f"{R}/ws_factor_regioncor.csv")
lc = pd.read_csv(f"{R}/ws_lumbar_course_vs_regions.csv")
rc["pair"] = rc.pair.map({"C-T": "C–T", "L-T": "T–L", "C-L": "C–L"}); rc = rc.sort_values(["feature", "pair"])
# Table 4
rows = []
for f in FE:
    q = p[p.feature == f].copy(); q['v'] = q.prev_adv if f == 'PS' else q.prev_any
    top = q.sort_values('v', ascending=False).head(3); peak = ', '.join(f'{l} ({v*100:.0f}%)' for l, v in zip(top.level, top.v))
    reg = {r: q[q.region == r].v.mean() * 100 for r in 'CTL'}
    cr = lambda lg: c[(c.feature == f) & (c.lag == lg)].r_partial.values[0]
    rr = rc[rc.feature == f].set_index('pair').r
    rcs = ', '.join(f'{k} {v:.2f}' for k, v in rr.items())
    rows.append([NAME[f] + (' ≥ IV' if f == 'PS' else ''), f"{reg['C']:.0f} / {reg['T']:.0f} / {reg['L']:.0f}", peak, f3(cr(1)), f3(cr(5)), f3(cr(10)), rcs])
t4 = pd.DataFrame(rows, columns=['Feature', 'Mean prevalence, %: C / T / L', 'Three most affected levels', 'ρ, adjacent', 'ρ, 5 levels apart', 'ρ, 10 levels apart', 'Region-factor correlations'])
open(f"{T}/T4_whole_spine.md", 'w').write(md(t4))
# S12 agreement
a = pd.read_csv(f"{R}/ws_lumbar_agreement.csv")
a['Feature'] = a.feature.map({'SS': 'Schneiderman score', 'PDB': 'Posterior bulge (0–2)', 'HIZ': 'HIZ', 'OMC': 'Modic change', 'SN': "Schmorl's node"})
s12 = a[['Feature', 'level', 'n', 'exact', 'kappa', 'prev_ws', 'prev_lumbar']].copy()
s12['exact'] = (s12.exact * 100).round(0).astype(int); s12['prev_ws'] = (s12.prev_ws * 100).round(1); s12['prev_lumbar'] = (s12.prev_lumbar * 100).round(1); s12['kappa'] = s12.kappa.map(f3)
s12.columns = ['Feature', 'Level', 'n', 'Exact agreement, %', 'κ†', 'Prevalence, whole-spine reading, %', 'Prevalence, lumbar consensus, %']
open(f"{T}/S12_ws_agreement.md", 'w').write(md(s12))
# S13 prevalence by level
w = p.copy(); w['v'] = np.where(w.feature == 'PS', w.prev_adv, w.prev_any) * 100
s13 = w.pivot(index='level', columns='feature', values='v').reindex(lab)[FE].round(1).reset_index()
s13.columns = ['Level'] + [NAME[f] + (' ≥ IV' if f == 'PS' else '') for f in FE]
open(f"{T}/S13_ws_prevalence.md", 'w').write(md(s13))
# S14 factor models
fit = pd.read_csv(f"{R}/ws_factor_fit.csv"); ecv = pd.read_csv(f"{R}/ws_factor_ecv.csv").set_index('feature')
rows = []
for f in FE:
    q = fit[fit.feature == f].set_index('model')
    for m, mn in [('M1', 'One general factor'), ('M2', 'Three correlated region factors'), ('M3', 'Bifactor (general + region)')]:
        r = q.loc[m]
        if pd.isna(r["cfi.scaled"]): rows.append([NAME[f], mn, int(r.n_items), '–', '–', '–', f'{r.srmr:.3f}', 'not identified']); continue
        e = ecv.loc[f, 'ECV_G'] if (m == 'M3' and f in ['PS', 'PDB', 'ADB', 'OMC']) else ''
        rows.append([NAME[f], mn, int(r.n_items), f'{r["cfi.scaled"]:.3f}', f'{r["tli.scaled"]:.3f}', f'{r["rmsea.scaled"]:.3f}', f'{r.srmr:.3f}', f'ECV = {e:.2f}' if e != '' else ''])
s14 = pd.DataFrame(rows, columns=['Feature', 'Model', 'Levels', 'CFI', 'TLI', 'RMSEA', 'SRMR', 'Note'])
cov = pd.read_csv(f"{R}/ws_factor_covariates.csv"); cv = cov[(cov.model == 'M2')]
s14b = rc.assign(r=rc.r.round(2), Feature=rc.feature.map(NAME), CI=[f'{r - 1.96*s:.2f} to {r + 1.96*s:.2f}' for r, s in zip(rc.r, rc.se)])[['Feature', 'pair', 'r', 'CI']]
s14b.columns = ['Feature', 'Regions', 'Correlation', '95% CI']
s14c = cv[cv['cov'] == 'age10'].assign(Feature=lambda x: x.feature.map(NAME))[['Feature', 'factor', 'beta_std']]
s14c = s14c.pivot(index='Feature', columns='factor', values='beta_std').reindex([NAME[f] for f in FE]).reset_index().fillna('–')
s14c.columns = ['Feature'] + [f'Age → {x} factor (std. β per 10 y)' for x in s14c.columns[1:]]
open(f"{T}/S14_ws_factor.md", 'w').write(md(s14) + '\n' + md(s14b) + '\n' + md(s14c))
# S15 risk factors by region
rf = pd.read_csv(f"{R}/ws_riskfactors_by_region.csv"); tt = pd.read_csv(f"{R}/ws_riskfactors_region_tests.csv")
TN = {'age_z': 'Age (SD)', 'male': 'Male', 'BMI_age40_z': 'BMI at 40 (SD)', 'Height_age40_z': 'Height at 40 (SD)', 'smoke': 'Ever smoked', 'injury': 'Lumbar injury', 'workload_z': 'Workload (SD)'}
rf['cell'] = [f'{b:.2f} ({b-1.96*s:.2f}, {b+1.96*s:.2f})' + ('*' if pp < 0.05 else '') for b, s, pp in zip(rf.beta, rf.se, rf.p)]
s15 = rf.pivot_table(index=['feature', 'term'], columns='region', values='cell', aggfunc='first').reset_index()
s15 = s15.merge(tt[['feature', 'term', 'p']], on=['feature', 'term'], how='left')
s15['p'] = s15.p.map(lambda x: pv(x) if pd.notna(x) else '–'); s15['feature'] = pd.Categorical(s15.feature, FE); s15['ord'] = s15.term.map(list(TN).index)
s15 = s15.sort_values(['feature', 'ord']); s15['feature'] = s15.feature.map(NAME).astype(str); s15['term'] = s15.term.map(TN)
s15 = s15[['feature', 'term', 'C', 'T', 'L', 'p']].fillna('–'); s15.columns = ['Feature', 'Predictor', 'Cervical', 'Thoracic', 'Lumbar', 'p (equal across regions)']
open(f"{T}/S15_ws_riskfactors.md", 'w').write(md(s15))
# S16 lumbar course
tv = pd.read_csv(f"{R}/ws_sens_tv_lumbar_course.csv")
ON = {'PS_C': 'Pfirrmann, cervical', 'PS_T': 'Pfirrmann, thoracic', 'PS_L': 'Pfirrmann, lumbar', 'PDB_C': 'Posterior bulge, cervical', 'PDB_T': 'Posterior bulge, thoracic',
      'HIZ_C': 'HIZ, cervical', 'OMC_C': 'Modic change, cervical', 'SN_T': "Schmorl's node, thoracic"}
PN = {'core_onset': 'Lumbar Core onset age (later)', 'Core1': 'Lumbar Core at first MRI', 'core_rate': 'Lumbar Core progression rate (first to second MRI)'}
tv['Predictor'] = tv.predictor.map(PN); tv['Outcome'] = tv.outcome.map(ON)
tv['All'] = [f'{b:.2f} ({b-1.96*s:.2f}, {b+1.96*s:.2f}); p {pv(pp)}' for b, s, pp in zip(tv.beta_all, tv.se_all, tv.p_all)]
tv['Excluding transitional vertebrae'] = [f'{b:.2f}; p {pv(pp)}' for b, pp in zip(tv.beta_noTV, tv.p_noTV)]
s16 = tv[['Predictor', 'Outcome', 'n_all', 'All', 'n_noTV', 'Excluding transitional vertebrae']].rename(columns={'n_all': 'n', 'n_noTV': 'n '})
open(f"{T}/S16_ws_lumbar_course.md", 'w').write(md(s16))
# S17 symptoms
s = pd.read_csv(f"{R}/ws_symptoms_by_region.csv"); s = s[(s.spec == 'regions+MCS') & (s.outcome != 'lbpyr2')]
TT = {'PS_Cz': 'Pfirrmann, cervical', 'PS_Tz': 'Pfirrmann, thoracic', 'PS_Lz': 'Pfirrmann, lumbar', 'PDB_Cz': 'Posterior bulge, cervical', 'PDB_Lz': 'Posterior bulge, lumbar',
      'HIZ_Lz': 'HIZ, lumbar', 'OMC_Lz': 'Modic change, lumbar', 'OMC_Cz': 'Modic change, cervical', 'mcs2z': 'SF-36 MCS'}
OU = {'sciatica2': 'Sciatica ever (OR)', 'oswestry2': 'ODI (points)', 'rmdq2': 'RMDQ (points)', 'vasSeverest2': 'Worst-ever VAS (mm)'}
s['cell'] = [f'{e:.2f} ({lo:.2f}, {hi:.2f})' + ('*' if pp < 0.05 else '') for e, lo, hi, pp in zip(s.est, s.lo, s.hi, s.p)]
s17 = s.pivot_table(index='term', columns='outcome', values='cell', aggfunc='first').reindex(list(TT)).reset_index()
s17['term'] = s17.term.map(TT); s17 = s17[['term'] + list(OU)]; s17.columns = ['Per SD of'] + [f'{OU[o]}, n = {s[s.outcome == o].n.iloc[0]}' for o in OU]
open(f"{T}/S17_ws_symptoms.md", 'w').write(md(s17))
print(md(t4))
