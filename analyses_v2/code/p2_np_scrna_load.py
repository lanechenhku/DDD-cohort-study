# Paper 2 (v3, 2026-11-06): public human NP scRNA-seq, Tu et al. Adv Sci 2022 (GEO GSE165722, 8 NP samples, author grades I-IV, BD Rhapsody).
# Step 1: read per-sample counts into one AnnData (run in /tmp/claude-0/geo/work; input = GSE165722_RAW.tar extracted).
import pandas as pd, scanpy as sc, anndata as ad, numpy as np, scipy.sparse as sp, glob, sys
D="/tmp/claude-0/geo/GSE165722"
grade={1:"I",2:"I",3:"II",4:"II",5:"III",6:"III",7:"IV",8:"IV"}
L=[]
for i in range(1,9):
    f=glob.glob(f"{D}/*_Sample{i}.counts.tsv.gz")[0]
    df=pd.read_csv(f,sep="\t",index_col=0)
    a=ad.AnnData(sp.csr_matrix(df.T.values.astype(np.float32)),obs=pd.DataFrame(index=[f"S{i}_{c}" for c in df.columns]),var=pd.DataFrame(index=df.index))
    a.obs["sample"]=f"S{i}"; a.obs["grade"]=grade[i]; L.append(a); print(i,a.shape,flush=True)
A=ad.concat(L,join="outer",fill_value=0); A.var_names_make_unique()
A.write_h5ad("raw.h5ad"); print(A)
