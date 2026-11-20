## Paper 1 (v3, 2026-11-06): sensitivity of Models 2-3 to keeping sport and age x height,
## which were dropped from Models 2-3 under the 2021 plan but have p < 0.05 for Core and MC in the re-estimated Model 1.
src <- readLines("/mnt/project-files/DDD_papers/analyses_v2/code/p1_wgee.R")
eval(parse(text = src[1:23]))   # data preparation and fit()
prs <- grep("PRS$", names(L), value = TRUE); prs <- prs[!grepl("^z_", prs)]
o <- list()
for (v in prs) for (y in c("Core", "HIZ", "MC", "DH")) {
  cr <- "zage + zage2 + male + z_packyr + z_injury + z_sport + z_BMI_age40 + z_Height_age40 + zage:z_packyr + zage:z_injury + zage:z_BMI_age40 + zage:z_Height_age40"
  if (grepl("BMI|height", v)) cr <- "zage + zage2 + male + z_packyr + z_injury + z_sport + zage:z_packyr + zage:z_injury"
  if (y != "Core") cr <- paste(cr, "+ z_core_onset")
  r <- fit(y, paste0(cr, " + z_", v, " + zage:z_", v), L)
  o[[length(o) + 1]] <- cbind(model = "M2_keep_sport_ageheight", pgs = v, r[grepl(v, r$term), ])
}
j <- c("EAS.BMI.PRS", "EUR.height.PRS", "EAS.LDL_C.PRS", "EAS.HDL_C.PRS", "EAS.T2DM.PRS", "EUR.lsBMD.PRS", "CHI.MD.PRS", "EUR.backpain.PRS", "EUR.multisitechrpain.PRS", "EUR.SpineOA.PRS")
jt <- paste(c(paste0("z_", j), paste0("zage:z_", j)), collapse = " + ")
for (y in c("Core", "HIZ", "MC", "DH")) {
  base <- "zage + zage2 + male + z_packyr + z_injury + z_sport + zage:z_packyr + zage:z_injury"
  if (y != "Core") base <- paste(base, "+ z_core_onset")
  r <- fit(y, paste(base, "+", jt), L)
  o[[length(o) + 1]] <- cbind(model = "M3_keep_sport", pgs = NA, r[grepl("PRS", r$term), ])
}
write.csv(do.call(rbind, o), file.path(RES, "p1_wgee_m23_keep_sport_ageheight.csv"), row.names = FALSE)
