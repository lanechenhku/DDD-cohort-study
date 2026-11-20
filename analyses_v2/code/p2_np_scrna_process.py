# Paper 2 (v3, 2026-11-06): QC (>=300 genes, <20% mito), log-normalisation, 2,000 HVGs, PCA, Harmony by sample, Leiden 0.5, UMAP.
# Cluster labels were assigned by canonical markers (see p2_figS4_S6_public.py, LABELS); output results/p2_ITGA9_NP_Tu2022_cells.csv.gz.
# Labels: 0,4 neutrophils; 6 myelocytes/immature neutrophils (LTF, CAMP); 11 promyelocytes (MS4A3, DEFA, MKI67); 5 macrophages (CD163, C1QA);
# 1,12 T/NK (CD3E, NKG7); 10 B cells, endothelial and other (CD74, MS4A1, CDH5); 9 erythroid (HBB, AHSP); 3 NP chondrocyte-like FMOD/PRG4;
# 8 NP chondrocyte-like DCN/CLU; 7 fibroblast-like COL1/COL3/TNC; 2 mixed (doublets).
import scanpy as sc, numpy as np, pandas as pd
A=sc.read_h5ad("raw.h5ad")
A.var["mt"]=A.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(A,qc_vars=["mt"],inplace=True,percent_top=None,log1p=False)
print(A.obs[["n_genes_by_counts","total_counts","pct_counts_mt"]].describe().round(1))
A=A[(A.obs.n_genes_by_counts>=300)&(A.obs.pct_counts_mt<20)].copy(); sc.pp.filter_genes(A,min_cells=10)
A.layers["counts"]=A.X.copy(); sc.pp.normalize_total(A,1e4); sc.pp.log1p(A)
sc.pp.highly_variable_genes(A,n_top_genes=2000,batch_key="sample",flavor="seurat")
sc.tl.pca(A,n_comps=30,mask_var="highly_variable")
import harmonypy as hm; Z=np.asarray(hm.run_harmony(A.obsm["X_pca"],A.obs,"sample",max_iter_harmony=20).Z_corr); A.obsm["X_pca_harmony"]=Z if Z.shape[0]==A.n_obs else Z.T; sc.pp.neighbors(A,use_rep="X_pca_harmony",n_neighbors=20)
sc.tl.leiden(A,resolution=0.5,flavor="igraph",n_iterations=2); sc.tl.umap(A)
A.write_h5ad("proc.h5ad"); print(A.shape, A.obs.leiden.value_counts().sort_index().to_dict())
