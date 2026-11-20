"""Whole-spine regional scores: (a) risk factors by region (GEE with region interactions); (b) lumbar longitudinal course vs
cervical/thoracic degeneration; (c) symptoms vs regional degeneration. 2026-11-06."""
import pandas as pd, numpy as np, statsmodels.formula.api as smf, statsmodels.api as sm
R = "/mnt/project-files/DDD_papers/analyses_ws/results"; V2 = "/mnt/project-files/DDD_papers/analyses_v2/results"
w = pd.read_csv(f"{R}/ws_wide.csv")
lab = ['C2-3','C3-4','C4-5','C5-6','C6-7','C7-T1'] + [f'T{i}-{i+1}' for i in range(1, 12)] + ['T12-L1'] + [f'L{i}-{i+1}' for i in range(1, 5)] + ['L5-S1']
reg = dict(C=lab[:6], T=lab[6:18], L=lab[18:])
FE = ['PS', 'PDB', 'ADB', 'HIZ', 'OMC', 'AMC', 'SN']
sc = w[['DNA', 'sex', 'age']].copy()
for f in FE:
    for r, ls in reg.items():
        x = w[[f'{f}_{l}' for l in ls]]
        if f != 'PS': x = (x > 0).astype(float).where(x.notna())
        sc[f'{f}_{r}'] = x.mean(axis=1, skipna=True) * (len(ls) if f != 'PS' else 1)   # PS: mean grade; others: number of affected levels
pp = pd.read_csv(f"{V2}/p1_person.csv").rename(columns={'DNA.No': 'DNA'})
d = sc.merge(pp[['DNA', 'IID', 'smoke1', 'packyr', 'sport', 'injury', 'workload', 'BMI_age40', 'Height_age40', 'core_onset',
                 'lbpyr2', 'sciatica2', 'oswestry2', 'rmdq2', 'vasSeverest2', 'vasToday2', 'mcs2', 'quesage2']], on='DNA', how='left')
d['male'] = (d.sex == 'M').astype(int); d['smoke'] = d.smoke1.astype(float)
for c in ['age', 'BMI_age40', 'Height_age40', 'workload']: d[c + '_z'] = (d[c] - d[c].mean()) / d[c].std()
sc.to_csv(f"{R}/ws_regional_scores.csv", index=False)
# (a) risk factors by region: stack standardised regional scores
rows, tests = [], []
X = ['age_z', 'male', 'BMI_age40_z', 'Height_age40_z', 'smoke', 'injury', 'workload_z']
for f in FE:
    st = []
    for r in reg:
        if d[f'{f}_{r}'].std() == 0: continue
        z = (d[f'{f}_{r}'] - d[f'{f}_{r}'].mean()) / d[f'{f}_{r}'].std()
        st.append(d[['DNA'] + X].assign(y=z, region=r))
    s = pd.concat(st).dropna()
    s['region'] = pd.Categorical(s.region, ['L', 'C', 'T'])
    form = 'y ~ 0 + region + ' + ' + '.join(f'region:{x}' for x in X)
    m = smf.gee(form, groups='DNA', data=s, cov_struct=sm.cov_struct.Exchangeable()).fit()
    for nm, b in m.params.items():
        if ':' in nm:
            r, x = nm.split(':'); rows.append(dict(feature=f, region=r.replace('region[', '').replace(']', ''), term=x, beta=b, se=m.bse[nm], p=m.pvalues[nm], n=s.DNA.nunique()))
    for x in X:   # equality of the coefficient across regions
        nms = [n for n in m.params.index if n.endswith(':' + x)]
        if len(nms) < 2: continue
        L = np.zeros((len(nms) - 1, len(m.params)))
        for k in range(1, len(nms)): L[k - 1, list(m.params.index).index(nms[0])] = 1; L[k - 1, list(m.params.index).index(nms[k])] = -1
        t = m.wald_test(L, scalar=True); tests.append(dict(feature=f, term=x, chi2=float(t.statistic), df=len(nms) - 1, p=float(t.pvalue)))
pd.DataFrame(rows).to_csv(f"{R}/ws_riskfactors_by_region.csv", index=False)
pd.DataFrame(tests).to_csv(f"{R}/ws_riskfactors_region_tests.csv", index=False)
# (b) lumbar course vs cervical/thoracic degeneration at the 2nd MRI
pl = pd.read_csv(f"{V2}/p1_long.csv")
c1 = pl[pl.t == 1][['DNA.No', 'age', 'Core']].rename(columns={'DNA.No': 'DNA', 'age': 'age1', 'Core': 'Core1'})
c2 = pl[pl.t == 2][['DNA.No', 'age', 'Core']].rename(columns={'DNA.No': 'DNA', 'age': 'age2', 'Core': 'Core2'})
b = d.merge(c1, on='DNA', how='left').merge(c2, on='DNA', how='left')
b['core_rate'] = (b.Core2 - b.Core1) / (b.age2 - b.age1)
b.loc[(b.age2 - b.age1) < 1, 'core_rate'] = np.nan
zz = lambda v: (v - v.mean()) / v.std()
out = []
for pred in ['core_onset', 'Core1', 'core_rate']:
    for y in ['PS_C', 'PS_T', 'PDB_C', 'PDB_T', 'OMC_C', 'HIZ_C', 'SN_T', 'PS_L']:
        dd = b[['age', 'male', pred, y]].dropna().copy(); dd['x'] = zz(dd[pred]); dd['yz'] = zz(dd[y])
        adj = 'age + male' + (' + Core1' if pred == 'core_rate' else '')
        if pred == 'core_rate': dd = b[['age', 'male', pred, 'Core1', y]].dropna().copy(); dd['x'] = zz(dd[pred]); dd['yz'] = zz(dd[y])
        m = smf.ols(f'yz ~ x + {adj}', dd).fit(cov_type='HC1')
        out.append(dict(predictor=pred, outcome=y, n=len(dd), beta=m.params['x'], se=m.bse['x'], p=m.pvalues['x']))
pd.DataFrame(out).to_csv(f"{R}/ws_lumbar_course_vs_regions.csv", index=False)
# (c) symptoms at the 2nd questionnaire (within 2 years of the whole-spine MRI)
s = d[(d.quesage2 - d.age).abs() <= 2].copy()
for c in ['PS_C', 'PS_T', 'PS_L', 'PDB_C', 'PDB_L', 'HIZ_L', 'OMC_L', 'OMC_C', 'mcs2']: s[c + 'z'] = zz(s[c])
sym = []
for y, fam in [('lbpyr2', 'bin'), ('sciatica2', 'bin'), ('oswestry2', 'lin'), ('rmdq2', 'lin'), ('vasSeverest2', 'lin')]:
    for spec, terms in [('regions', 'PS_Cz + PS_Tz + PS_Lz + PDB_Cz + PDB_Lz + HIZ_Lz + OMC_Lz + OMC_Cz'),
                        ('regions+MCS', 'PS_Cz + PS_Tz + PS_Lz + PDB_Cz + PDB_Lz + HIZ_Lz + OMC_Lz + OMC_Cz + mcs2z')]:
        dd = s.copy()
        if fam == 'bin': dd[y] = dd[y].map({True: 1, False: 0, 'True': 1, 'False': 0, 1: 1, 0: 0, 'Yes': 1, 'No': 0})
        dd = dd[[y, 'age_z', 'male'] + [t.strip() for t in terms.split('+')]].dropna()
        m = (smf.logit if fam == 'bin' else smf.ols)(f'{y} ~ age_z + male + {terms}', dd).fit(disp=0, **({} if fam == 'bin' else {'cov_type': 'HC1'}))
        for t in [t.strip() for t in terms.split('+')]:
            sym.append(dict(outcome=y, spec=spec, term=t, n=len(dd), est=np.exp(m.params[t]) if fam == 'bin' else m.params[t],
                            lo=np.exp(m.conf_int().loc[t, 0]) if fam == 'bin' else m.conf_int().loc[t, 0],
                            hi=np.exp(m.conf_int().loc[t, 1]) if fam == 'bin' else m.conf_int().loc[t, 1], p=m.pvalues[t], scale='OR' if fam == 'bin' else 'beta'))
pd.DataFrame(sym).to_csv(f"{R}/ws_symptoms_by_region.csv", index=False)
print('symptom n window', len(s))
