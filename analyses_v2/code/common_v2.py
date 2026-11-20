"""Shared paths for the v2 (2026-10-23) re-analyses. Inputs are the staged copies of Lane's Mac files."""
import os
UP = "/mnt/user-data/uploads/DDD_Paper"
D8 = os.path.join(UP, "DDD 2020-04", "DDD202008")
D4 = os.path.join(UP, "DDD 2020-04")
RES = "/mnt/project-files/DDD_papers/analyses_v2/results"
def p(*a): return os.path.join(*a)
# v3 (2026-11-06): cleaned copies of the unadjusted grading files with the 69 undated T2 records removed
D8C = "/mnt/project-files/DDD_papers/analyses_v2/data_clean/DDD202008"
def d8(f):
    c = os.path.join(D8C, f)
    return c if os.path.exists(c) else os.path.join(D8, f)
