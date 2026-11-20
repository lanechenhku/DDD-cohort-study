"""Paper 1: Turnbull NPMLE of age at onset (interval/left/right censored) per feature — independent check of SISE."""
import numpy as np, pandas as pd
from common_v2 import *
e = pd.read_csv(d8("all.emergence.data.20201215.csv"), index_col=0)
a = pd.read_csv(p(D8, "allcovariate.20200917.csv"), index_col=0); a["IID"] = a.index
e = e.merge(a[["IID", "sex"]], on="IID", how="left")
F = ["ss", "db", "dn", "dh", "hiz", "mc", "sn"]

def intervals(f, sub):
    L, R = [], []
    for _, r in sub.iterrows():
        obs = [(r[f"age{t}"], r[f"{f}.t{t}"]) for t in (1, 2, 3) if pd.notna(r[f"age{t}"]) and pd.notna(r[f"{f}.t{t}"])]
        if not obs: continue
        s = [v for _, v in obs]
        if any(s[i] == 1 and s[i + 1] == 0 for i in range(len(s) - 1)): continue   # implies regression: excluded
        if s[0] == 1: L.append(0.0); R.append(obs[0][0]); continue
        k = next((i for i, v in enumerate(s) if v == 1), None)
        if k is None: L.append(obs[-1][0]); R.append(np.inf)
        else: L.append(obs[k - 1][0]); R.append(obs[k][0])
    return np.round(np.array(L) * 4) / 4, np.round(np.array(R) * 4) / 4   # ages to 0.25 y keeps the EM tractable

def turnbull(L, R, tol=1e-6, it=3000):
    lefts = np.unique(L); rights = np.unique(R[np.isfinite(R)])
    pts = sorted(set(lefts) | set(rights))
    # Turnbull innermost intervals [q, p]
    ev = sorted([(x, 0) for x in L] + [(x, 1) for x in R if np.isfinite(x)])
    inner = []; last_left = None
    for x, typ in ev:
        if typ == 0: last_left = x
        elif last_left is not None: inner.append((last_left, x)); last_left = None
    inner = sorted(set(inner)) + [(np.inf, np.inf)]
    q = np.array([i[0] for i in inner]); pp = np.array([i[1] for i in inner])
    A = ((q[None, :] >= L[:, None]) & (pp[None, :] <= R[:, None])).astype(float)
    A[:, -1] = np.isinf(R)
    s = np.full(len(inner), 1 / len(inner))
    for _ in range(it):
        d = A @ s; s_new = s * (A / d[:, None]).mean(0)
        if np.abs(s_new - s).max() < tol: s = s_new; break
        s = s_new
    return pp, s

def summarise(pp, s, ages=(30, 40, 50, 60, 70, 80)):
    F_ = np.cumsum(s); fin = np.isfinite(pp)
    def at(x): return F_[fin][pp[fin] <= x].max() if (pp[fin] <= x).any() else 0.0
    med = pp[fin][F_[fin] >= 0.5].min() if (F_[fin] >= 0.5).any() else np.nan
    return med, {f"F{x}": at(x) for x in ages}, F_[fin].max()

rows, curves = [], []
rng = np.random.default_rng(1)
for grp, sub in [("All", e), ("Women", e[e.sex == "F"]), ("Men", e[e.sex == "M"])]:
    for f in F:
        Lf, Rf = intervals(f, sub); pp, s = turnbull(Lf, Rf); med, cr, fmax = summarise(pp, s)
        bmeds = []
        if grp == "All":
            for b in range(30):
                i = rng.integers(0, len(Lf), len(Lf)); pb, sb = turnbull(Lf[i], Rf[i], tol=1e-5, it=800); bmeds.append(summarise(pb, sb)[0])
            fin = np.isfinite(pp); curves += [dict(feature=f.upper(), age=a_, F=c) for a_, c in zip(pp[fin], np.cumsum(s)[fin])]
        lo, hi = (np.nanpercentile(bmeds, 2.5), np.nanpercentile(bmeds, 97.5)) if bmeds else (np.nan, np.nan)
        rows.append(dict(group=grp, feature=f.upper(), n=len(Lf), n_left=int((Lf == 0).sum()), n_interval=int(((Lf > 0) & np.isfinite(Rf)).sum()),
                         n_right=int(np.isinf(Rf).sum()), median=med, med_lo=lo, med_hi=hi, F_max_observed=fmax, **cr))
        print(grp, f, len(Lf), med, cr, flush=True)
pd.DataFrame(rows).to_csv(p(RES, "p1_turnbull_onset.csv"), index=False)
pd.DataFrame(curves).to_csv(p(RES, "p1_turnbull_curves.csv"), index=False)
