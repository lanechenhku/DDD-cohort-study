"""Derived tables from the cleaned cohort files (all.emergence.data / core.unadjusted)."""
import numpy as np, pandas as pd
from common import SRC

FEATS = ["ss", "db", "dn", "dh", "hiz", "mc", "sn"]
AGEBINS = [0, 20, 30, 40, 50, 60, 200]
AGELAB = ["<20", "20–29", "30–39", "40–49", "50–59", "≥60"]
STATES = ["Normal", "Core only", "Core + HIZ", "Core + MC", "Core + DH", "Core + two", "All four"]


def emergence():
    return pd.read_csv(SRC["emerg"], index_col=0)


def wilson(k, n, z=1.96):
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return p, c - h, c + h


def prevalence_by_age(e):
    e = e.copy()
    e["any.t1"] = e[[f + ".t1" for f in FEATS]].max(axis=1, skipna=False)
    e["ag"] = pd.cut(e.age1, AGEBINS, right=False, labels=AGELAB)
    rows = []
    for f in FEATS + ["any"]:
        g = e.dropna(subset=[f + ".t1", "ag"]).groupby("ag", observed=True)[f + ".t1"].agg(["sum", "count"])
        for a, r in g.iterrows():
            p, lo, hi = wilson(r["sum"], r["count"])
            rows.append(dict(feature=f, agegroup=a, k=int(r["sum"]), n=int(r["count"]), p=p, lo=lo, hi=hi))
    return pd.DataFrame(rows)


def states(e):
    """Composite state at each scan: Core = SS or DB present at any level."""
    out = {}
    for t in (1, 2, 3):
        core = ((e[f"ss.t{t}"] > 0) | (e[f"db.t{t}"] > 0)).astype(float)
        core[e[f"ss.t{t}"].isna()] = np.nan
        h, m, d = e[f"hiz.t{t}"], e[f"mc.t{t}"], e[f"dh.t{t}"]
        k = h + m + d
        s = np.select([core.isna() | k.isna(), (core == 0) & (k == 0), core == 0, k == 0,
                       (k == 1) & (h == 1), (k == 1) & (m == 1), (k == 1) & (d == 1), k == 3],
                      [None, "Normal", "Other", "Core only", "Core + HIZ", "Core + MC", "Core + DH", "All four"],
                      "Core + two")
        out[t] = pd.Series(s, index=e.index)
    return pd.DataFrame(out)


def state_freq(st):
    rows = []
    for t in (1, 2, 3):
        s = st[t].dropna()
        n = len(s)
        for k in STATES:
            c = int((s == k).sum())
            p, lo, hi = wilson(c, n)
            rows.append(dict(time=f"T{t}", state=k, k=c, n=n, p=p, lo=lo, hi=hi))
    return pd.DataFrame(rows)


def transitions(st):
    pairs = [st[[1, 2]].set_axis(["a", "b"], axis=1), st[[2, 3]].set_axis(["a", "b"], axis=1)]
    x = pd.concat(pairs, ignore_index=True).dropna()
    x = x[x.a.isin(STATES) & x.b.isin(STATES)]
    return pd.crosstab(x.a, x.b).reindex(index=STATES, columns=STATES, fill_value=0)
