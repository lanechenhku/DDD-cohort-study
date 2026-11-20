"""v3 (2026-11-06, Lane decision): drop second-scan (T2) gradings of participants with no T2 scan date.

69 participants have T2 consensus grades in the 2020 rater/emergence files but no T2 date, age or scanner
in allpheno_v3 (and no T3). Their T2 grades change about twice as often as those of dated participants
(e.g. SS regression 8/43 vs 14/1,592), so they are treated as mis-linked records and set to missing.
Usage: python3 p1_exclude_undated_t2.py SRC_DIR OUT_DIR   (SRC_DIR = .../DDD 2020-04/DDD202008)
Writes cleaned copies of the files below into OUT_DIR (scripts read them via common_v2.d8()) and
writes OUT_DIR/excluded_T2_IIDs.csv. Adjusted files (core.machine_adjusted, complication.adjusted,
t2_phenotype.machine.adjusted, GWAS_phenotype) already have no T2 values for these people.
"""
import os, re, sys, pandas as pd
src, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
F = ["ss", "db", "dn", "dh", "hiz", "mc", "sn"]
e = pd.read_csv(os.path.join(src, "all.emergence.data.20201215.csv"), index_col=0)
ids = e.loc[e[[f + ".t2" for f in F]].notna().any(axis=1) & e.age2.isna(), "IID"]
pd.DataFrame({"IID": sorted(ids)}).to_csv(os.path.join(out, "excluded_T2_IIDs.csv"), index=False)
T2 = re.compile(r"(\.t2$|^[a-z]+\d\.2$|^core\.t2$|t1t2|t2t3|^age\.t2$)")
CLEAN = ["all.emergence.data.20201215.csv", "unadjusted.levelwise.data.seven_features20210202.csv",
         "data.for.feature_onset.20201210.csv", "core.unadjusted.20210202.csv", "rater_consensus_phenotype_20201006.csv"]
for f in CLEAN:
    d = pd.read_csv(os.path.join(src, f), index_col=0 if f.startswith(("all.emergence", "data.for")) else None)
    cols = [c for c in d.columns if T2.search(c)]
    m = d.IID.isin(ids); n = int(d.loc[m, cols].notna().any(axis=1).sum())
    d.loc[m, cols] = pd.NA
    d.to_csv(os.path.join(out, f), index=not isinstance(d.index, pd.RangeIndex) or f.startswith(("all.emergence", "data.for")))
    print(f"{f}: blanked {len(cols)} T2 columns for {n} participants")
print("excluded", len(ids))
