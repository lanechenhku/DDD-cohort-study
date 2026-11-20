# Whole-spine latent structure per feature (2026-11-06).
# Ordinal/binary CFA (lavaan, WLSMV) with age and sex as covariates of every factor (MIMIC).
# M1 one general factor; M2 three correlated region factors (C, T, L); M3 bifactor (G + orthogonal C/T/L).
# Levels with < 1.5% prevalence (fewer than ~30 cases) are dropped for that feature; a region with < 3 usable levels is dropped.
.libPaths(c("/tmp/claude-0/Rlib", .libPaths()))
suppressMessages(library(lavaan))
R <- "/mnt/project-files/DDD_papers/analyses_ws/results"
w <- read.csv(file.path(R, "ws_wide.csv"), check.names = FALSE)
w$male <- as.numeric(w$sex == "M"); w$age10 <- (w$age - 50) / 10
lab <- c("C2-3","C3-4","C4-5","C5-6","C6-7","C7-T1", paste0("T",1:11,"-",2:12), "T12-L1", paste0("L",1:4,"-",2:5), "L5-S1")
reg <- c(rep("C",6), rep("T",12), rep("L",5))
recode <- function(f, y) {
  if (f == "PS") return(cut(y, c(-Inf, 2, 3, Inf), labels = FALSE) - 1)   # <=II, III, >=IV
  as.numeric(y > 0)
}
feats <- c("PS","PDB","ADB","HIZ","OMC","AMC","SN")
fitrows <- list(); corrows <- list(); ecvrows <- list(); agerows <- list()
for (f in feats) {
  d <- data.frame(age10 = w$age10, male = w$male)
  keep <- c(); kreg <- c()
  for (i in seq_along(lab)) {
    y <- recode(f, w[[paste0(f, "_", lab[i])]])
    p <- mean(y > 0, na.rm = TRUE); if (f == "PS") p <- min(mean(y == 0, na.rm=TRUE), mean(y == 2, na.rm=TRUE))
    if (p >= 0.015) { v <- paste0("x", i); d[[v]] <- y; keep <- c(keep, v); kreg <- c(kreg, reg[i]) }
  }
  regs <- names(which(table(kreg) >= 3)); sel <- kreg %in% regs; keep <- keep[sel]; kreg <- kreg[sel]
  ind <- function(r) paste(keep[kreg == r], collapse = " + ")
  m1 <- paste0("G =~ ", paste(keep, collapse = " + "), "\nG ~ age10 + male")
  m2 <- paste(c(sapply(regs, function(r) paste0(r, " =~ ", ind(r))), sapply(regs, function(r) paste0(r, " ~ age10 + male"))), collapse = "\n")
  m3 <- paste(c(paste0("G =~ ", paste(keep, collapse = " + ")), sapply(regs, function(r) paste0(r, "s =~ ", ind(r))),
                "G ~ age10 + male", sapply(regs, function(r) paste0(r, "s ~ age10 + male"))), collapse = "\n")
  for (mn in c("M1","M2","M3")) {
    mod <- switch(mn, M1 = m1, M2 = m2, M3 = m3)
    fit <- tryCatch(sem(mod, data = d, ordered = keep, estimator = "WLSMV", std.lv = TRUE, orthogonal = (mn == "M3"), missing = "pairwise"),
                    error = function(e) NULL, warning = function(w) NULL)
    if (is.null(fit)) fit <- tryCatch(sem(mod, data = d, ordered = keep, estimator = "WLSMV", std.lv = TRUE, orthogonal = (mn == "M3"), missing = "pairwise"), error = function(e) NULL)
    if (is.null(fit) || !lavInspect(fit, "converged")) { fitrows[[length(fitrows)+1]] <- data.frame(feature=f, model=mn, converged=FALSE); next }
    fm <- fitMeasures(fit, c("chisq.scaled","df.scaled","cfi.scaled","tli.scaled","rmsea.scaled","srmr"))
    fitrows[[length(fitrows)+1]] <- data.frame(feature=f, model=mn, converged=TRUE, n_items=length(keep), regions=paste(regs, collapse=""), as.list(round(fm, 4)))
    ss <- standardizedSolution(fit)
    if (mn == "M2") { cc <- ss[ss$op == "~~" & ss$lhs %in% regs & ss$rhs %in% regs & ss$lhs != ss$rhs, ]
      corrows[[length(corrows)+1]] <- data.frame(feature=f, pair=paste(cc$lhs, cc$rhs, sep="-"), r=round(cc$est.std,3), se=round(cc$se,3)) }
    if (mn == "M3") { ld <- ss[ss$op == "=~", ]; g <- sum(ld$est.std[ld$lhs == "G"]^2); s <- sum(ld$est.std[ld$lhs != "G"]^2)
      byreg <- sapply(regs, function(r) { it <- keep[kreg == r]; gg <- sum(ld$est.std[ld$lhs=="G" & ld$rhs %in% it]^2); gg / (gg + sum(ld$est.std[ld$lhs==paste0(r,"s") & ld$rhs %in% it]^2)) })
      ecvrows[[length(ecvrows)+1]] <- data.frame(feature=f, ECV_G=round(g/(g+s),3), as.list(round(byreg,3)), mean_loading_G=round(mean(ld$est.std[ld$lhs=="G"]),3)) }
    if (mn %in% c("M2","M3")) { rr <- ss[ss$op == "~" & ss$rhs %in% c("age10","male"), ]
      agerows[[length(agerows)+1]] <- data.frame(feature=f, model=mn, factor=rr$lhs, cov=rr$rhs, beta_std=round(rr$est.std,3), p=signif(rr$pvalue,3)) }
  }
  cat(f, "done\n")
}
bind <- function(x) do.call(rbind, lapply(x, function(z) { z[setdiff(unique(unlist(lapply(x, names))), names(z))] <- NA; z }))
write.csv(bind(fitrows), file.path(R, "ws_factor_fit.csv"), row.names = FALSE)
write.csv(bind(corrows), file.path(R, "ws_factor_regioncor.csv"), row.names = FALSE)
write.csv(bind(ecvrows), file.path(R, "ws_factor_ecv.csv"), row.names = FALSE)
write.csv(bind(agerows), file.path(R, "ws_factor_covariates.csv"), row.names = FALSE)
