## Paper 1: WGEE Model 1 with MRI-informed attendance weights (sensitivity); appended to p1_wgee_model1.csv
suppressMessages(library(geepack))
RES <- "/mnt/project-files/DDD_papers/analyses_v2/results"
L <- read.csv(file.path(RES, "p1_long.csv")); z <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
L$male <- 1 - L$female; L$zage <- z(L$age); L$zage2 <- L$zage^2
for (v in c("packyr", "injury", "sport", "BMI_age40", "Height_age40", "core_onset")) L[[paste0("z_", v)]] <- z(L[[v]])
for (y in c("Core", "HIZ", "MC", "DH")) L[[paste0("y_", y)]] <- z(L[[y]])
cov_full <- "zage + zage2 + male + z_packyr + z_injury + z_sport + z_BMI_age40 + z_Height_age40 + zage:male + zage:z_packyr + zage:z_injury + zage:z_sport + zage:z_BMI_age40 + zage:z_Height_age40"
o <- list()
for (y in c("Core", "HIZ", "MC", "DH")) {
  f <- as.formula(paste0("y_", y, " ~ ", if (y == "Core") cov_full else paste(cov_full, "+ z_core_onset")))
  d <- L[complete.cases(L[, all.vars(f)]) & !is.na(L$ipw_mri), ]; d <- d[order(d$IID, d$t), ]
  s <- summary(geeglm(f, id = IID, waves = t, data = d, weights = d$ipw_mri, corstr = "ar1"))$coefficients
  o[[y]] <- data.frame(model = "M1_ipw_with_MRI", outcome = y, term = rownames(s), est = s[, 1], se = s[, 2], p = s[, 4],
                       lo = s[, 1] - 1.96 * s[, 2], hi = s[, 1] + 1.96 * s[, 2], n_obs = nrow(d), n_id = length(unique(d$IID)))
}
M <- read.csv(file.path(RES, "p1_wgee_model1.csv")); M <- M[M$model != "M1_ipw_with_MRI", ]
write.csv(rbind(M, do.call(rbind, o)), file.path(RES, "p1_wgee_model1.csv"), row.names = FALSE)
print(do.call(rbind, o)[, c("outcome", "term", "est", "p")][grepl("male$|BMI_age40$|Height|injury$|onset", do.call(rbind, o)$term), ])
