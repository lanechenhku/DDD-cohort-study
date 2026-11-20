"""Paper 1: inter-rater agreement (Table S1) from per-rater grades in allpheno_v3.xlsx (sheet all_MRI)."""
import numpy as np, pandas as pd, itertools
from common_v2 import *
x = pd.read_excel(p(D4, "allpheno_v3.xlsx"), sheet_name="all_MRI")
R = ["flo", "jaro", "jason", "kc"]
def kappa(a, b, weighted):
    cats = np.union1d(np.unique(a), np.unique(b)); k = len(cats); idx = {c: i for i, c in enumerate(cats)}
    O = np.zeros((k, k))
    for u, v in zip(a, b): O[idx[u], idx[v]] += 1
    O /= O.sum(); E = np.outer(O.sum(1), O.sum(0))
    W = (np.subtract.outer(np.arange(k), np.arange(k)) ** 2 / max((k - 1) ** 2, 1)) if weighted else (1 - np.eye(k))
    return 1 - (W * O).sum() / (W * E).sum() if (W * E).sum() > 0 else np.nan
rows = []
for f in ["ss", "db", "hiz", "mc", "sn"]:
    for t in (1, 2):
        for r1, r2 in itertools.combinations(R, 2):
            a, b = [], []
            for lv in range(1, 6):
                c1, c2 = f"{f}{lv}.{t}.{r1}", f"{f}{lv}.{t}.{r2}"
                if c1 in x and c2 in x:
                    d = x[[c1, c2]].apply(pd.to_numeric, errors="coerce").dropna(); a += list(d[c1]); b += list(d[c2])
            if len(a) < 50: continue
            a, b = np.array(a), np.array(b)
            rows.append(dict(feature=f.upper(), scan=f"T{t}", pair=f"{r1}-{r2}", n_disc_ratings=len(a),
                             pct_agree=100 * (a == b).mean(), kappa=kappa(a, b, f in ("ss", "db"))))
k = pd.DataFrame(rows)
# anonymise rater names
lab = {r: f"R{i+1}" for i, r in enumerate(R)}
k["pair"] = k.pair.apply(lambda s: "–".join(lab[q] for q in s.split("-")))
k.to_csv(p(RES, "p1_interrater_pairs.csv"), index=False)
s = k.groupby(["feature", "scan"]).apply(lambda g: pd.Series(dict(pairs=len(g), n=g.n_disc_ratings.sum(), pct_agree=np.average(g.pct_agree, weights=g.n_disc_ratings),
      kappa_median=g.kappa.median(), kappa_min=g.kappa.min(), kappa_max=g.kappa.max()))).reset_index()
s.to_csv(p(RES, "p1_interrater_summary.csv"), index=False); print(s.round(2).to_string()); print(k.round(2).to_string())
