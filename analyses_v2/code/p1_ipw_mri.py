"""Paper 1: alternative attendance weights that also use sex and the previous scan's MRI findings (sensitivity)."""
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from common_v2 import *
L = pd.read_csv(p(RES, "p1_long.csv"))
w = []
for t in (2, 3):
    prev = L[L.t == t - 1].copy(); prev["att"] = prev.IID.isin(L[L.t == t].IID).astype(int)
    prev["zC"] = (prev.Core - prev.Core.mean()) / prev.Core.std()
    m = smf.logit("att ~ age + I(age**2) + female + zC + hiz_bin + mc_bin + dh_bin", prev).fit(disp=0)
    prev["pp"] = m.predict(prev); w.append(prev[["IID", "pp"]].assign(t=t))
w = pd.concat(w)
L = L.drop(columns=[c for c in ["pp", "ipw_mri"] if c in L]).merge(w, on=["IID", "t"], how="left")
L["ipw_mri"] = np.where(L.t == 1, 1.0, 1 / L.pp)
L.to_csv(p(RES, "p1_long.csv"), index=False); print(L.ipw_mri.describe())
