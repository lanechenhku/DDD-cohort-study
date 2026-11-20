# Check of the 2020 BMI/height-at-age-40 prediction (predicted*-20201129.txt):
# refit LMMs on the two questionnaire visits and compare with the 2020 predictions.
suppressMessages(library(lme4))
d0 <- "/mnt/user-data/uploads/DDD_Paper/DDD 2020-04/DDD202008"
out <- "/mnt/project-files/DDD_papers/analyses_v2/results"
cv <- read.csv(file.path(d0, "allcovariate.20200917.csv"))
pb <- read.delim(file.path(d0, "predictedBMI-20201129.txt"))
ph <- read.delim(file.path(d0, "predictedHeight-20201129.txt"))
cv$yob <- as.numeric(substr(cv$dob, 1, 4))
long <- do.call(rbind, lapply(1:2, function(v) data.frame(
  id = cv$DNA.No, sex = cv$sex, yob = cv$yob,
  age = cv[[paste0("quesage", v)]], bmi = cv[[paste0("bmi", v)]],
  ht = cv[[paste0("height", v)]], visit = v)))
long <- long[!is.na(long$age) & !is.na(long$sex), ]
i <- which(long$ht > 3); long$ht[i] <- long$ht[i] / 100   # cm -> m if needed
cat("rows with BMI:", sum(!is.na(long$bmi)), " people:", length(unique(long$id[!is.na(long$bmi)])),
    " with 2 visits:", sum(table(long$id[!is.na(long$bmi)]) == 2), "\n")
print(summary(long$age)); print(summary(long$bmi)); print(summary(long$ht))
long$a <- (long$age - 40) / 10; long$y <- (long$yob - 1960) / 10
fit_one <- function(resp, with_yob) {
  f <- if (with_yob) as.formula(paste(resp, "~ a + I(a^2) + sex + sex:a + y + (1 + a | id)"))
       else as.formula(paste(resp, "~ a + I(a^2) + sex + sex:a + (1 + a | id)"))
  lmer(f, data = long[!is.na(long[[resp]]), ], REML = TRUE,
       control = lmerControl(check.nobs.vs.nRE = "ignore"))
}
res <- list()
for (resp in c("bmi", "ht")) for (wy in c(FALSE, TRUE)) {
  m <- fit_one(resp, wy); cat("\n====", resp, "yob =", wy, "\n"); print(summary(m)$coefficients); print(VarCorr(m))
  nd <- unique(long[, c("id", "sex", "yob")]); nd$a <- 0; nd$y <- (nd$yob - 1960) / 10
  nd$pred <- predict(m, newdata = nd, allow.new.levels = TRUE)
  old <- if (resp == "bmi") pb[, c("DNA.No", "BMI_age40")] else ph[, c("DNA.No", "Height_age40")]
  names(old) <- c("id", "old")
  mm <- merge(nd, old, by = "id"); mm <- mm[!is.na(mm$old), ]
  # observed at visit nearest age 40
  ob <- long[!is.na(long[[resp]]), ]; ob <- ob[order(ob$id, abs(ob$age - 40)), ]; ob <- ob[!duplicated(ob$id), c("id", resp)]
  names(ob)[2] <- "obs"; mm <- merge(mm, ob, by = "id", all.x = TRUE)
  vc <- as.data.frame(VarCorr(m))
  res[[paste(resp, wy)]] <- data.frame(trait = resp, yob_term = wy, n_obs = nobs(m), n_id = ngrps(m)[["id"]],
    sd_int = vc$sdcor[1], sd_slope = vc$sdcor[2], sd_resid = vc$sdcor[nrow(vc)],
    r_with_2020 = cor(mm$pred, mm$old), mean_abs_diff_2020 = mean(abs(mm$pred - mm$old)),
    r_with_observed = cor(mm$pred, mm$obs, use = "complete"), AIC = AIC(refitML(m)))
}
# does the 2020 prediction depend on yob beyond age/sex? (would indicate a cohort term)
pb2 <- merge(pb, unique(long[, c("id", "yob")]), by.x = "DNA.No", by.y = "id", suffixes = c("", ".cv"))
cat("\n2020 BMI_age40 ~ yob + sex:\n"); print(summary(lm(BMI_age40 ~ yob + sex, data = pb))$coefficients)
cat("2020 Height_age40 ~ yob + sex:\n"); print(summary(lm(Height_age40 ~ yob + sex, data = ph))$coefficients)
r <- do.call(rbind, res); print(r, digits = 3)
write.csv(r, file.path(out, "p1_bmi_height_model_check.csv"), row.names = FALSE)
