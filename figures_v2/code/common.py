"""Paths. Scripts run on Lane's Mac VM from DDD_Paper/Figures_redrawn_2026-10/code."""
import os
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get("DDD_ROOT", os.path.abspath(os.path.join(HERE, "..", "..")))   # DDD_Paper
OUTROOT = os.path.abspath(os.path.join(HERE, ".."))              # Figures_redrawn_2026-10
D8 = os.path.join(ROOT, "DDD 2020-04", "DDD202008")
G21 = os.path.join(ROOT, "DDD 2020-04", "GWAS20210313")
# v3 (2026-11-06): grading files with the 69 undated second-scan records removed (analyses_v2/code/p1_exclude_undated_t2.py)
CLEAN = os.environ.get("DDD_CLEAN", os.path.join(OUTROOT, "data", "clean_DDD202008"))
def _c(f):
    c = os.path.join(CLEAN, f)
    return c if os.path.exists(c) else os.path.join(D8, f)
SRC = {
    "emerg": _c("all.emergence.data.20201215.csv"),
    "core": _c("core.unadjusted.20210202.csv"),
    "level": _c("unadjusted.levelwise.data.seven_features20210202.csv"),
    "fig2d": os.path.join(D8, "fig2d.csv"),
    "fm_ncausal": os.path.join(D8, "fm.ncausal.txt"),
    "fm_snp": os.path.join(D8, "fm.snp.txt"),
    "gwas_pheno": os.path.join(D8, "GWAS_phenotype.20210313.txt"),
    "locus": os.path.join(G21, "locus.3.MC.csv"),
    "locus_snps": os.path.join(G21, "locus_snps.csv"),
    "fuma_genes": os.path.join(G21, "FUMA_job122508", "genes.txt"),
}
WORK = os.path.join(OUTROOT, "_work", "data") + "/"     # thinned sumstats
DATA = os.path.join(OUTROOT, "data") + "/"              # derived / digitised tables
P1 = os.path.join(OUTROOT, "Paper1")
P2 = os.path.join(OUTROOT, "Paper2")
