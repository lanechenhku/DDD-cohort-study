"""Whole-spine MRI (2nd lumbar MRI wave): clean and reshape. 2026-11-06.
Input: Mac DDD_Paper/DDD whole spine/whole_spine/whole_spine.xlsx (staged copy).
Changes vs the 2022 pipeline: drop the row flagged 'UNABLE TO COUNT LEVELS' and the 3 rows flagged 'duplicated'
(two are exact DNA-number duplicates); out-of-range codes set to missing (no median imputation);
columns parsed by name, so the lumbar HIZ/Schmorl order swap in the raw sheet is handled explicitly."""
import re, pandas as pd, numpy as np, os
SRC = "/mnt/user-data/uploads/DDD_Paper/DDD whole spine/whole_spine/whole_spine.xlsx"
OUT = "/mnt/project-files/DDD_papers/analyses_ws/results"
d = pd.read_excel(SRC)
cols = list(d.columns); ph = cols[7:cols.index('DNA No2')]
d = d[d['DNA no'].notna()].copy()
note = d['0'].astype(str).str.strip().str.lower()
log = {'rows_with_DNA': len(d)}
drop = note.isin(['duplicated', 'unable to count levels'])
log['dropped_flagged'] = int(drop.sum()); d = d[~drop]
assert not d['DNA no'].duplicated().any()
d['tv'] = np.select([note[~drop].str.contains('lumbar'), note[~drop].str.contains('sacral')], ['lumbarisation', 'sacralisation'], '')
FEAT = [('anterior bulge', 'ADB'), ('posterior bulge', 'PDB'), ('pfirrmann', 'PS'), ('schneiderman', 'SS'), ('hiz', 'HIZ'),
        ("schmorl", 'SN'), ('anterior modic', 'AMC'), ('ant marrow', 'AMC'), ('anterior marrow', 'AMC'), ('overall modic', 'OMC'), ('overall marrow', 'OMC')]
RANGE = {'ADB': (0, 3), 'PDB': (0, 3), 'PS': (1, 5), 'SS': (0, 3), 'HIZ': (0, 1), 'SN': (0, 1), 'AMC': (0, 1), 'OMC': (0, 1)}
ORDER = [f'C{i}' for i in range(2, 8)] + [f'T{i}' for i in range(1, 13)] + [f'L{i}' for i in range(1, 6)]
LAB = ['C2-3', 'C3-4', 'C4-5', 'C5-6', 'C6-7', 'C7-T1'] + [f'T{i}-{i+1}' for i in range(1, 12)] + ['T12-L1'] + [f'L{i}-{i+1}' for i in range(1, 5)] + ['L5-S1']
def parse(c):
    s = c.lower().replace('’', "'").replace('  ', ' ').strip()
    f = next((code for k, code in FEAT if k in s), None)
    if f is None: return None
    m = re.search(r'\b([ctl]) ?(\d+)', s) or re.search(r'grading j(\d)', s)
    if 'grading j' in s: lev = 'L' + re.search(r'j(\d)', s).group(1)
    else:
        num = m.group(2)
        if m.group(1) == 't' and int(num) > 12: num = num[0]  # 'T34 Schmorl's' = T3/4
        lev = m.group(1).upper() + num
    return f, lev
mp = {}
for c in ph:
    r = parse(c)
    if r: mp[c] = r
assert len(mp) == 184, len(mp)
assert len(set(mp.values())) == 184
wide = pd.DataFrame({'DNA': d['DNA no'].astype(int).values, 'sex': d['Sex'].values, 'age': d['Exact Age (2nd MRI)'].values,
                     'mri_date': d['Date of 2nd MRI'].values, 'tv': d['tv'].values})
oor = 0
for c, (f, lev) in mp.items():
    v = pd.to_numeric(d[c], errors='coerce').values.astype(float)
    lo, hi = RANGE[f]; bad = (~np.isnan(v)) & ((v < lo) | (v > hi) | (v != np.round(v)))
    oor += bad.sum(); v[bad] = np.nan
    wide[f'{f}_{LAB[ORDER.index(lev)]}'] = v
log['out_of_range_set_missing'] = int(oor)
log['n'] = len(wide); log['female'] = int((wide.sex == 'F').sum()); log['male'] = int((wide.sex == 'M').sum())
log['tv_lumbarisation'] = int((wide.tv == 'lumbarisation').sum()); log['tv_sacralisation'] = int((wide.tv == 'sacralisation').sum())
log['age_mean_sd_range'] = (round(wide.age.mean(), 1), round(wide.age.std(), 1), round(wide.age.min(), 1), round(wide.age.max(), 1))
# questionnaire block at the 2nd questionnaire
Q = {'q_date': 'DateOf2ndQuestionnaire', 'weight2': 'BodyWeightAt2ndQuestionnaire', 'height2': 'HeightAt2ndQuestionnaire', 'occupation': 'Occupation',
     'workload_hist': 'Historical Workload', 'workload_cur': 'Current Workload', 'smoke2': 'Smoke', 'sport2': 'Sport', 'injury2': 'Injury of lumbar',
     'lbp_ever2': 'Low  Back Pain', 'lbp_year2': 'LBP in this year', 'sciatica2': 'Sciatica', 'odi2': 'Oswestry Disbility total score',
     'rmdq2': 'Roland Morris Disability total score', 'vas_today2': 'Visual Analog Scale-Today (mm)', 'vas_worst2': 'Visual Analog Scale-severest (mm)',
     'pcs2': 'Physical Component (PCS)', 'mcs2': 'Mental  Component (MCS)'}
for k, c in Q.items(): wide[k] = d[c].values
wide = wide.copy()
wide.to_csv(os.path.join(OUT, 'ws_wide.csv'), index=False)
long = wide.melt(id_vars=['DNA', 'sex', 'age', 'tv'], value_vars=[f'{f}_{l}' for f in RANGE for l in LAB], var_name='var', value_name='y')
long[['feature', 'level']] = long['var'].str.split('_', expand=True)
long['level_idx'] = long['level'].map({l: i for i, l in enumerate(LAB)})
long['region'] = np.where(long.level_idx <= 5, 'C', np.where(long.level_idx <= 17, 'T', 'L'))
long.drop(columns='var').to_csv(os.path.join(OUT, 'ws_long.csv'), index=False)
pd.Series(log).to_csv(os.path.join(OUT, 'ws_clean_log.csv'))
print(log)
