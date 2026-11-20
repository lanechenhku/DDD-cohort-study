"""Whole-spine reader (J.P.Y.C.) vs his own earlier lumbar grading of the same 2nd-MRI scans (allpheno_v3 all_MRI, *.2.jason),
with the other lumbar T2 reader (jaro) and the lumbar consensus as benchmarks. Lumbar levels 1..5 = L1/2..L5/S1."""
import pandas as pd, numpy as np, os
from sklearn.metrics import cohen_kappa_score
R = "/mnt/project-files/DDD_papers/analyses_ws/results"; V2 = "/mnt/project-files/DDD_papers/analyses_v2"
x = pd.read_excel("/mnt/user-data/uploads/DDD_Paper/DDD 2020-04/allpheno_v3.xlsx", sheet_name="all_MRI")
x = x[x.date2.notna()].rename(columns={'DNA.No': 'DNA'})        # drop undated T2 records (v3 decision)
w = pd.read_csv(f"{R}/ws_wide.csv")
pp = pd.read_csv(f"{V2}/results/p1_person.csv")
lw = pd.read_csv(f"{V2}/data_clean/DDD202008/unadjusted.levelwise.data.seven_features20210202.csv")
lw = lw.merge(pp[['IID', 'DNA.No']], on='IID', how='left').rename(columns={'DNA.No': 'DNA'})
m = w.merge(x, on='DNA', how='inner', suffixes=('', '_x')).merge(lw, on='DNA', how='left', suffixes=('', '_c'))
lab = ['L1-2', 'L2-3', 'L3-4', 'L4-5', 'L5-S1']
F = [('SS', 'ss', False), ('PDB', 'db', False), ('HIZ', 'hiz', True), ('OMC', 'mc', True), ('SN', 'sn', True)]
def num(s): return pd.to_numeric(s, errors='coerce')
def src(name, f, g, i, l):
    if name == 'WS': a = num(m[f'{f}_{l}']); return a.clip(upper=2) if f == 'PDB' else a
    if name == 'cons': return num(m[f'{g}{i}.2'])
    c = f'{g}{i}.2.{name}'; return num(m[c]) if c in m else pd.Series(np.nan, index=m.index)
pairs = [('WS', 'jason', 'same reader, two sessions'), ('WS', 'jaro', 'different readers'),
         ('jason', 'jaro', 'different readers (lumbar)'), ('WS', 'cons', 'reader vs consensus')]
rows = []
for f, g, binar in F:
    for a_, b_, what in pairs:
        A_all, B_all = [], []
        for i, l in enumerate(lab, 1):
            a, b = src(a_, f, g, i, l), src(b_, f, g, i, l); ok = a.notna() & b.notna()
            if ok.sum() < 30: continue
            A, B = a[ok].astype(int), b[ok].astype(int); A_all += list(A); B_all += list(B)
            rows.append(dict(feature=f, pair=f'{a_}-{b_}', what=what, level=l, n=int(ok.sum()), exact=(A == B).mean(),
                             kappa=cohen_kappa_score(A, B, weights=None if binar else 'linear'),
                             mean_a=A.mean(), mean_b=B.mean(), prev_a=(A > 0).mean(), prev_b=(B > 0).mean()))
        if A_all:
            A, B = np.array(A_all), np.array(B_all)
            rows.append(dict(feature=f, pair=f'{a_}-{b_}', what=what, level='pooled', n=len(A), exact=(A == B).mean(),
                             kappa=cohen_kappa_score(A, B, weights=None if binar else 'linear'),
                             mean_a=A.mean(), mean_b=B.mean(), prev_a=(A > 0).mean(), prev_b=(B > 0).mean()))
d = pd.DataFrame(rows); d.to_csv(f"{R}/ws_reader_agreement.csv", index=False)
print('people with WS and Jason T2 lumbar:', int((m['ss3.2.jason'].notna() & m['SS_L3-4'].notna()).sum()))
print(d[d.level == 'pooled'].round(3).to_string(index=False))
print(d[(d.pair == 'WS-jason')].round(2).to_string(index=False))
