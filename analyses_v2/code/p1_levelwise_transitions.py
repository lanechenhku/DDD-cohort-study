"""Level-wise transitions (binary feature present at disc level) between consecutive scans; regression among affected and among all."""
import os, pandas as pd, numpy as np
from common_v2 import D8, RES, d8
from statsmodels.stats.proportion import proportion_confint
e = pd.read_csv(d8("all.emergence.data.20201215.csv"), index_col=0)
rows = []
for f in ["ss", "db", "dn", "dh", "hiz", "mc", "sn"]:
    for a, b in [("t1", "t2"), ("t2", "t3")]:
        for L in range(1, 6):
            x, y = e[f"{f}.L{L}.{a}"], e[f"{f}.L{L}.{b}"]
            ok = x.notna() & y.notna(); x, y = x[ok], y[ok]
            n10 = int(((x == 1) & (y == 0)).sum()); n11 = int(((x == 1) & (y == 1)).sum())
            n01 = int(((x == 0) & (y == 1)).sum()); n00 = int(((x == 0) & (y == 0)).sum())
            aff = n10 + n11
            lo, hi = proportion_confint(n10, aff, method="wilson") if aff else (np.nan, np.nan)
            rows.append(dict(feature=f.upper(), interval=f"{a.upper()}→{b.upper()}", level=f"L{L}", n=len(x), n00=n00, n01=n01, n11=n11, n10=n10,
                             regress_all=n10 / len(x), regress_affected=n10 / aff if aff else np.nan, ci_lo=lo, ci_hi=hi))
d = pd.DataFrame(rows); d.to_csv(os.path.join(RES, "p1_levelwise_transitions.csv"), index=False)
lab = {1: "L1/2", 2: "L2/3", 3: "L3/4", 4: "L4/5", 5: "L5/S1"}
out = ["| Feature | Interval | " + " | ".join(lab.values()) + " |", "|---|---|" + "---|" * 5]
for (f, iv), g in d.groupby(["feature", "interval"], sort=False):
    cells = [f"{r.regress_affected*100:.1f} ({r.n10}/{r.n10+r.n11})" if r.n10 + r.n11 else "–" for r in g.itertuples()]
    out.append(f"| {f} | {iv} | " + " | ".join(cells) + " |")
open(os.path.join(RES, "tables_md", "S7_levelwise_transitions.md"), "w").write("\n".join(out) + "\n")
print("\n".join(out)); print(d[d.feature == "SS"][["interval", "level", "regress_all"]].round(3).to_string())
