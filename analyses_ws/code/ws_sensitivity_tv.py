"""Sensitivity: exclude participants with transitional lumbosacral vertebrae (lumbarisation/sacralisation). 2026-11-06."""
import pandas as pd, numpy as np, subprocess, os, shutil
RM = "/mnt/project-files/DDD_papers/analyses_ws/results"
w = pd.read_csv(f"{RM}/ws_wide.csv")
tmp = RM + "/_tv"; os.makedirs(tmp, exist_ok=True)
w[w.tv.isna() | (w.tv == '')].to_csv(f"{tmp}/ws_wide.csv", index=False)
src = open('ws_corr_decay.py').read().replace('R = "/mnt/project-files/DDD_papers/analyses_ws/results"', f'R = "{tmp}"')
exec(compile(src, 'corr', 'exec'))
src2 = open('ws_models.py').read().replace('R = "/mnt/project-files/DDD_papers/analyses_ws/results";', f'R = "{tmp}";')
exec(compile(src2, 'models', 'exec'))
a = pd.read_csv(f"{tmp}/ws_lumbar_course_vs_regions.csv"); b = pd.read_csv(f"{RM}/ws_lumbar_course_vs_regions.csv")
m = b.merge(a, on=['predictor', 'outcome'], suffixes=('_all', '_noTV')); m.to_csv(f"{RM}/ws_sens_tv_lumbar_course.csv", index=False)
print(m[['predictor', 'outcome', 'n_all', 'beta_all', 'n_noTV', 'beta_noTV', 'p_noTV']].round(3).to_string())
a = pd.read_csv(f"{tmp}/ws_corr_by_lag.csv"); b = pd.read_csv(f"{RM}/ws_corr_by_lag.csv")
m = b.merge(a, on=['feature', 'lag'], suffixes=('_all', '_noTV')); m.to_csv(f"{RM}/ws_sens_tv_corr_by_lag.csv", index=False)
print(m[m.lag.isin([1, 10])][['feature', 'lag', 'r_partial_all', 'r_partial_noTV']].round(3).to_string())
