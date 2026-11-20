## Paper 1: re-estimation of the WGEE risk-factor models (Models 1-3 + sensitivity) with geepack.
suppressMessages(library(geepack))
RES <- "/mnt/project-files/DDD_papers/analyses_v2/results"
L <- read.csv(file.path(RES, "p1_long.csv"))
z <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
L$male <- 1 - L$female
L$zage <- z(L$age); L$zage2 <- L$zage^2
for (v in c("packyr", "injury", "sport", "BMI_age40", "Height_age40", "core_onset")) L[[paste0("z_", v)]] <- z(L[[v]])
for (y in c("Core", "HIZ", "MC", "DH", "Core_unadj")) L[[paste0("y_", y)]] <- z(L[[y]])
L$base_age <- ave(L$age, L$IID, FUN = min)
prs <- grep("PRS$", names(L), value = TRUE)
for (v in prs) L[[paste0("z_", v)]] <- z(L[[v]])
L <- L[order(L$IID, L$t), ]

fit <- function(y, rhs, d, w = "ipw", cs = "ar1") {
  f <- as.formula(paste0("y_", y, " ~ ", rhs))
  d <- d[complete.cases(d[, all.vars(f)]) & !is.na(d[[w]]), ]
  d <- d[order(d$IID, d$t), ]
  m <- geeglm(f, id = IID, waves = t, data = d, weights = d[[w]], corstr = cs, family = gaussian)
  s <- summary(m)$coefficients
  data.frame(outcome = y, term = rownames(s), est = s[, 1], se = s[, 2], p = s[, 4],
             lo = s[, 1] - 1.96 * s[, 2], hi = s[, 1] + 1.96 * s[, 2], n_obs = nrow(d), n_id = length(unique(d$IID)))
}
cov_full <- "zage + zage2 + male + z_packyr + z_injury + z_sport + z_BMI_age40 + z_Height_age40 + zage:male + zage:z_packyr + zage:z_injury + zage:z_sport + zage:z_BMI_age40 + zage:z_Height_age40"
rhs1 <- function(y) if (y == "Core") cov_full else paste(cov_full, "+ z_core_onset")
out <- list()
for (y in c("Core", "HIZ", "MC", "DH")) {
  out[[length(out) + 1]] <- cbind(model = "M1", fit(y, rhs1(y), L))
  out[[length(out) + 1]] <- cbind(model = "M1_unweighted", fit(y, rhs1(y), transform(L, one = 1), w = "one"))
  out[[length(out) + 1]] <- cbind(model = "M1_exchangeable", fit(y, rhs1(y), L, cs = "exchangeable"))
  out[[length(out) + 1]] <- cbind(model = "M1_cumulative_ipw", fit(y, rhs1(y), L, w = "ipw_cum"))
  out[[length(out) + 1]] <- cbind(model = "M1_probands_only", fit(y, rhs1(y), L[L$proband == "True", ]))
  out[[length(out) + 1]] <- cbind(model = "M1_baseline_age_18_55", fit(y, rhs1(y), L[L$base_age >= 18 & L$base_age <= 55, ]))
}
L$scanner <- factor(L$scanner)
out[[length(out) + 1]] <- cbind(model = "M1_unadjCore_scanner_covariate", fit("Core_unadj", paste(cov_full, "+ scanner"), L))
M1 <- do.call(rbind, out); write.csv(M1, file.path(RES, "p1_wgee_model1.csv"), row.names = FALSE)

## Model 2: one PGS at a time (+ PGS x age); reduced covariates; BMI/height omitted for the BMI/height PGS
cov_red <- "zage + zage2 + male + z_packyr + z_injury + z_BMI_age40 + z_Height_age40 + zage:z_packyr + zage:z_injury + zage:z_BMI_age40"
o2 <- list()
for (v in prs) for (y in c("Core", "HIZ", "MC", "DH")) {
  cr <- cov_red
  if (grepl("BMI|height", v)) cr <- "zage + zage2 + male + z_packyr + z_injury + zage:z_packyr + zage:z_injury"
  if (y != "Core") cr <- paste(cr, "+ z_core_onset")
  r <- fit(y, paste0(cr, " + z_", v, " + zage:z_", v), L)
  o2[[length(o2) + 1]] <- cbind(model = "M2", pgs = v, r[grepl(v, r$term), ])
}
M2 <- do.call(rbind, o2); write.csv(M2, file.path(RES, "p1_wgee_model2.csv"), row.names = FALSE)

## Model 3: joint (East Asian score where available)
j <- c("EAS.BMI.PRS", "EUR.height.PRS", "EAS.LDL_C.PRS", "EAS.HDL_C.PRS", "EAS.T2DM.PRS", "EUR.lsBMD.PRS", "CHI.MD.PRS", "EUR.backpain.PRS", "EUR.multisitechrpain.PRS", "EUR.SpineOA.PRS")
jt <- paste(c(paste0("z_", j), paste0("zage:z_", j)), collapse = " + ")
o3 <- list()
for (y in c("Core", "HIZ", "MC", "DH")) {
  base <- "zage + zage2 + male + z_packyr + z_injury + zage:z_packyr + zage:z_injury"
  if (y != "Core") base <- paste(base, "+ z_core_onset")
  o3[[length(o3) + 1]] <- cbind(model = "M3", fit(y, paste(base, "+", jt), L))
  o3[[length(o3) + 1]] <- cbind(model = "M3_plus_measured_BMI_height", fit(y, paste(base, "+ z_BMI_age40 + z_Height_age40 +", jt), L))
}
M3 <- do.call(rbind, o3); write.csv(M3, file.path(RES, "p1_wgee_model3.csv"), row.names = FALSE)
cat("done\n")
