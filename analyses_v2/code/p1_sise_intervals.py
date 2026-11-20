"""Paper 1 (v3): onset intervals per feature for SISE (same rules as p1_turnbull.py, ages to 0.01 y). Output: results/p1_sise_intervals.csv"""
import numpy as np, pandas as pd
from common_v2 import *
e = pd.read_csv(d8("all.emergence.data.20201215.csv"), index_col=0)
F = ["ss", "db", "dn", "dh", "hiz", "mc", "sn"]
rows = []
for f in F:
    for _, r in e.iterrows():
        obs = [(r[f"age{t}"], r[f"{f}.t{t}"]) for t in (1, 2, 3) if pd.notna(r[f"age{t}"]) and pd.notna(r[f"{f}.t{t}"])]
        if not obs: continue
        s = [v for _, v in obs]
        if any(s[i] == 1 and s[i + 1] == 0 for i in range(len(s) - 1)): continue   # implies regression: excluded
        if s[0] == 1: L, R, ty = np.nan, obs[0][0], "Left"
        else:
            k = next((i for i, v in enumerate(s) if v == 1), None)
            if k is None: L, R, ty = obs[-1][0], np.nan, "Right"
            else: L, R, ty = obs[k - 1][0], obs[k][0], "Interval"
        rows.append(dict(IID=r.IID, feature=f.upper(), left=L, right=R, type=ty))
d = pd.DataFrame(rows); d.to_csv(p(RES, "p1_sise_intervals.csv"), index=False)
print(d.groupby(["feature", "type"]).size().unstack())
