# Paper 1 (v3, 2026-11-06): SISE onset curves re-run with the published SISE package (github.com/tubbsjd/SISE)
# on the cleaned data (intervals from p1_sise_intervals.py). Main penalty logNe (n.obs = 2); logNm and logN as robustness.
# Usage: Rscript p1_sise.R [FEATURE ...]   (default: all seven). Writes results/p1_sise_<F>.rds and p1_sise_curves.csv
.libPaths(c("/tmp/claude-0/Rlib", .libPaths()))
suppressMessages(library(SISE))
RES <- "/mnt/project-files/DDD_papers/analyses_v2/results"
iv <- read.csv(file.path(RES, "p1_sise_intervals.csv"))
args <- commandArgs(TRUE); FE <- if (length(args)) args else c("SS", "DB", "DN", "DH", "HIZ", "MC", "SN")
for (f in FE) {
  d <- round(iv[iv$feature == f, c("left", "right")], 1)
  t0 <- Sys.time()
  r <- smoothTB(d, n.obs = 2, left.bound = 0, penalty = c("logNe", "logNm", "logN"), n.dec = 1)
  r$dat <- NULL
  saveRDS(r, file.path(RES, paste0("p1_sise_", f, ".rds")))
  # as in the original analysis (pred.onset8.RDS): SS, DB and DN imputed for all censoring types, other features for left- and interval-censored only
  imp <- impute.time(d, r$TB.o.smooth.logNe, r$time.o, type = if (f %in% c("SS", "DB", "DN")) "LIR" else "LI", n.dec = 1)
  write.csv(data.frame(IID = iv$IID[iv$feature == f], onset = imp, type = iv$type[iv$feature == f]),
            file.path(RES, paste0("p1_sise_imputed_", f, ".csv")), row.names = FALSE)
  cat(f, "N", r$N, "Ne", round(r$Ne), "bw logNe/logNm/logN", round(r$logNe_bw, 2), round(r$logNm_bw, 2), round(r$logN_bw, 2),
      "time", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
}
