import pandas as pd, numpy as np, sys
ph=sys.argv[1]
d=pd.read_csv(f"GWAS_sumstat/summary.{ph}.csv",usecols=["CHR","BP","P","Z"],dtype={"CHR":"int8","BP":"int32","P":"float64","Z":"float32"})
n=len(d); lam=np.nanmedian(d.Z.values.astype(float)**2)/0.454936
rng=np.random.default_rng(1)
keep=(d.P<0.01)|(rng.random(n)<0.05)
d[keep][["CHR","BP","P"]].to_csv(f"Figures_redrawn_2026-10/_work/data/manh_{ph}.csv.gz",index=False)
p=np.sort(d.P.dropna().values); m=len(p)
e=(np.arange(1,m+1)-0.5)/m
idx=np.unique(np.concatenate([np.arange(min(20000,m)),np.arange(0,m,250),[m-1]]))
pd.DataFrame({"exp":-np.log10(e[idx]),"obs":-np.log10(p[idx]),"rank":idx+1}).to_csv(f"Figures_redrawn_2026-10/_work/data/qq_{ph}.csv",index=False)
print(ph,n,m,round(lam,4),keep.sum(),p[0])
