"""Linkage to the lumbar cohort, agreement of the two lumbar gradings at the 2nd MRI, prevalence atlas. 2026-11-06."""
import pandas as pd, numpy as np, os
from sklearn.metrics import cohen_kappa_score
R = "/mnt/project-files/DDD_papers/analyses_ws/results"
V2 = "/mnt/project-files/DDD_papers/analyses_v2"
w = pd.read_csv(f"{R}/ws_wide.csv"); L = pd.read_csv(f"{R}/ws_long.csv")
pl = pd.read_csv(f"{V2}/results/p1_long.csv"); pp = pd.read_csv(f"{V2}/results/p1_person.csv")
t2 = pl[pl.t == 2][['DNA.No', 'IID', 'age']].rename(columns={'DNA.No': 'DNA', 'age': 'age_lumbarT2'})
m = w[['DNA', 'age']].merge(t2, on='DNA', how='left')
out = {}
out['linked_any'] = int(w.DNA.isin(pp['DNA.No']).sum())
out['with_lumbar_T2_record'] = int(m.age_lumbarT2.notna().sum())
dif = (m.age - m.age_lumbarT2).abs()
out['age_diff_median_days'] = round(float(dif.median() * 365.25), 1)
out['age_diff_le_0.1y'] = int((dif <= 0.1).sum())
# agreement: whole-spine lumbar grading vs consensus lumbar grading at T2 (levelwise, 69 undated T2 removed)
lw = pd.read_csv(f"{V2}/data_clean/DDD202008/unadjusted.levelwise.data.seven_features20210202.csv")
lw = lw.merge(pp[['IID', 'DNA.No']], on='IID', how='left').rename(columns={'DNA.No': 'DNA'})
mm = w.merge(lw, on='DNA', how='inner')
lab = ['L1-2', 'L2-3', 'L3-4', 'L4-5', 'L5-S1']
rows = []
for f, g, binar in [('SS', 'ss', False), ('PDB', 'db', False), ('HIZ', 'hiz', True), ('OMC', 'mc', True), ('SN', 'sn', True)]:
    for i, l in enumerate(lab, 1):
        a = mm[f'{f}_{l}']; b = mm[f'{g}{i}.2']; ok = a.notna() & b.notna()
        A, B = a[ok].astype(int), b[ok].astype(int)
        if f == 'PDB': A = A.clip(upper=2)
        rows.append(dict(feature=f, level=l, n=int(ok.sum()), exact=round((A == B).mean(), 3),
                         kappa=round(cohen_kappa_score(A, B, weights=None if binar else 'linear'), 3),
                         prev_ws=round((A > 0).mean(), 3), prev_lumbar=round((B > 0).mean(), 3)))
ag = pd.DataFrame(rows); ag.to_csv(f"{R}/ws_lumbar_agreement.csv", index=False)
# prevalence atlas: any grade (SS>0, DB>0, others =1) and advanced (PS>=4, SS>=2, DB>=2)
L['any'] = np.where(L.feature == 'PS', L.y >= 3, L.y > 0).astype(float); L.loc[L.y.isna(), 'any'] = np.nan
L['adv'] = np.select([L.feature == 'PS', L.feature == 'SS', L.feature.isin(['ADB', 'PDB'])], [L.y >= 4, L.y >= 2, L.y >= 2], np.nan).astype(float)
L.loc[L.y.isna(), 'adv'] = np.nan
prev = L.groupby(['feature', 'level', 'level_idx', 'region']).agg(n=('y', 'count'), mean=('y', 'mean'), prev_any=('any', 'mean'), prev_adv=('adv', 'mean')).reset_index().sort_values(['feature', 'level_idx'])
prev.to_csv(f"{R}/ws_prevalence_by_level.csv", index=False)
bins = [0, 35, 45, 55, 65, 100]; L['agegrp'] = pd.cut(L.age, bins, labels=['<35', '35-44', '45-54', '55-64', '65+'])
pa = L.groupby(['feature', 'level', 'level_idx', 'agegrp'], observed=True)['any'].mean().reset_index(); pa.to_csv(f"{R}/ws_prevalence_by_age_level.csv", index=False)
ps = L.groupby(['feature', 'level', 'level_idx', 'sex'])['any'].mean().reset_index(); ps.to_csv(f"{R}/ws_prevalence_by_sex_level.csv", index=False)
out['agegrp_n'] = w.assign(g=pd.cut(w.age, bins)).g.value_counts().sort_index().to_dict()
pd.Series(out).to_csv(f"{R}/ws_linkage_log.csv"); print(out)
print(ag.to_string())
print(prev[prev.feature.isin(['SS', 'PS', 'PDB'])].pivot(index='level', columns='feature', values='prev_any').reindex(L.drop_duplicates('level').sort_values('level_idx').level).round(3).T.to_string())
