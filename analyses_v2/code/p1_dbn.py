"""Paper 1: re-estimation of the dynamic Bayesian networks (feature level and disc level) with bootstrap arc strengths.
Two-slice model learned on stacked consecutive scan pairs (T1->T2, T2->T3); arcs only from the earlier to the later slice
(blacklist: arcs into the earlier slice and within a slice); hill-climbing with BIC (pgmpy). Reconstruction, not the original code."""
import numpy as np, pandas as pd, itertools, warnings, logging
warnings.filterwarnings("ignore"); logging.disable(logging.WARNING)
from pgmpy.estimators import HillClimbSearch
from pgmpy.causal_discovery import ExpertKnowledge
from common_v2 import *
e = pd.read_csv(d8("all.emergence.data.20201215.csv"), index_col=0)
lv = pd.read_csv(d8("unadjusted.levelwise.data.seven_features20210202.csv"))
FE = ["ss", "db", "dn", "dh", "hiz", "mc", "sn"]

def pairs_feature():
    out = []
    for a, b in [(1, 2), (2, 3)]:
        d = pd.DataFrame({**{f"{f.upper()}_0": e[f"{f}.t{a}"] for f in FE}, **{f"{f.upper()}_1": e[f"{f}.t{b}"] for f in FE}})
        out.append(d.dropna())
    return pd.concat(out).astype(int).astype(str).astype("category")

def pairs_level():
    out = []
    for a, b in [(1, 2), (2, 3)]:
        cols = {}
        for k in range(1, 6):
            for s, t in [("0", a), ("1", b)]:
                c = [f"{f}{k}.{t}" for f in ["ss", "db", "hiz", "mc", "sn"] if f"{f}{k}.{t}" in lv]
                cols[f"L{k}_{s}"] = (lv[c] > 0).any(axis=1).where(lv[c].notna().all(axis=1))
        out.append(pd.DataFrame(cols).dropna())
    return pd.concat(out).astype(int).astype(str).astype("category")

def learn(d):
    nodes = list(d.columns)
    bl = [(u, v) for u, v in itertools.permutations(nodes, 2) if v.endswith("_0") or (u.endswith("_1") and v.endswith("_1"))]
    hc = HillClimbSearch(d)
    m = hc.estimate(scoring_method="bic-d", expert_knowledge=ExpertKnowledge(forbidden_edges=bl), show_progress=False)
    return set(m.edges())

for name, d in [("feature", pairs_feature()), ("level", pairs_level())]:
    full = learn(d)
    rng = np.random.default_rng(2026); B = 200; cnt = {}
    for _ in range(B):
        for ed in learn(d.iloc[rng.integers(0, len(d), len(d))]): cnt[ed] = cnt.get(ed, 0) + 1
    rows = [dict(model=name, frm=u.rsplit("_", 1)[0], to=v.rsplit("_", 1)[0], autoregressive=u.rsplit("_", 1)[0] == v.rsplit("_", 1)[0],
                 in_full_data_network=(u, v) in full, bootstrap_strength=c / B) for (u, v), c in cnt.items()]
    r = pd.DataFrame(rows).sort_values("bootstrap_strength", ascending=False)
    r.to_csv(p(RES, f"p1_dbn_{name}.csv"), index=False)
    print(name, "n pairs", len(d)); print(r[r.bootstrap_strength >= 0.3].to_string())
