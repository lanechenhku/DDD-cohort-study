"""Between-level correlation of each feature, adjusted for age and sex (partial Spearman), by distance and region. 2026-11-06."""
import pandas as pd, numpy as np
from scipy.stats import rankdata
R = "/mnt/project-files/DDD_papers/analyses_ws/results"
w = pd.read_csv(f"{R}/ws_wide.csv")
lab = ['C2-3','C3-4','C4-5','C5-6','C6-7','C7-T1'] + [f'T{i}-{i+1}' for i in range(1, 12)] + ['T12-L1'] + [f'L{i}-{i+1}' for i in range(1, 5)] + ['L5-S1']
reg = ['C'] * 6 + ['T'] * 12 + ['L'] * 5
Z = np.column_stack([np.ones(len(w)), rankdata(w.age), (w.sex == 'M').astype(float)])
def resid(v):
    ok = ~np.isnan(v); r = np.full_like(v, np.nan); y = rankdata(v[ok]); X = Z[ok]
    r[ok] = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]; return r
rows = []
for f in ['PS', 'SS', 'PDB', 'ADB', 'HIZ', 'OMC', 'AMC', 'SN']:
    M = np.column_stack([w[f'{f}_{l}'].values.astype(float) for l in lab])
    raw = pd.DataFrame(M).corr(method='spearman').values
    Rz = np.column_stack([resid(M[:, j]) if np.nanstd(M[:, j]) > 0 else np.full(len(w), np.nan) for j in range(23)])
    par = pd.DataFrame(Rz).corr().values
    for i in range(23):
        for j in range(i + 1, 23):
            rows.append(dict(feature=f, a=lab[i], b=lab[j], lag=j - i, ra=reg[i], rb=reg[j], same_region=reg[i] == reg[j], r_raw=raw[i, j], r_partial=par[i, j]))
c = pd.DataFrame(rows); c.to_csv(f"{R}/ws_corr_pairs.csv", index=False)
s = c.groupby(['feature', 'lag'])[['r_raw', 'r_partial']].mean().reset_index(); s.to_csv(f"{R}/ws_corr_by_lag.csv", index=False)
c['blk'] = np.where(c.same_region, c.ra + c.ra, np.where(c.ra < c.rb, c.ra + c.rb, c.rb + c.ra))
b = c.groupby(['feature', 'blk'])['r_partial'].mean().unstack().round(3); b.to_csv(f"{R}/ws_corr_by_block.csv")
print(s[s.lag.isin([1, 2, 5, 10, 20])].pivot(index='feature', columns='lag', values='r_partial').round(3)); print(b)
