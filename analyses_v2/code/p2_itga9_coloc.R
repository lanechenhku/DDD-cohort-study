# Paper 2 (v3, 2026-11-06): colocalisation of the MC.avg GWAS signal at ITGA9 with ITGA9 eQTLs (eQTL Catalogue).
# GWAS: hg19 region 3:37.2-37.8 Mb (data/itga9_region_<trait>.csv); eQTL: results/p2_ITGA9_region_eQTLCatalogue.csv (hg38).
# hg38 = hg19 - 41,491 throughout this region (checked for rs1549596, rs56902407, rs17227306, rs73055216).
# Caveat: GWAS in East Asians, eQTL in Europeans; coloc.abf assumes one causal variant per trait and shared LD.
.libPaths(c("/tmp/claude-0/Rlib", .libPaths())); suppressMessages(library(coloc))
B <- "/mnt/project-files/DDD_papers/analyses_v2"
ds <- read.csv(file.path(B, "results/p2_phenotype_descriptives.csv")); sdY <- setNames(ds$sd, ds$phenotype)
nG <- setNames(ds$n, ds$phenotype)
eq <- read.csv(file.path(B, "results/p2_ITGA9_region_eQTLCatalogue.csv")); eq$pos19 <- eq$position + 41491
lab <- c(QTD000021 = "BLUEPRINT monocyte", QTD000356 = "GTEx whole blood", QTD000373 = "Lepik 2017 blood",
         QTD000090 = "FUSION adipose", QTD000603 = "Perez 2022 classical monocyte")
nE <- c(QTD000021 = 191, QTD000356 = 670, QTD000373 = 471, QTD000090 = 271, QTD000603 = 193)
out <- NULL
for (tr in c("MC", "HIZ", "DH", "Core")) {
  g <- read.csv(file.path(B, sprintf("data/itga9_region_%s.csv", tr)))
  for (d in names(lab)) {
    e <- eq[eq$dataset == d, ]; e <- e[!duplicated(e$pos19), ]
    m <- merge(g, e, by.x = "BP", by.y = "pos19")
    m <- m[is.finite(m$Beta) & is.finite(m$SE) & m$SE > 0 & is.finite(m$se) & m$se > 0, ]
    m$snp <- as.character(m$BP)
    r <- coloc.abf(list(beta = m$Beta, varbeta = m$SE^2, type = "quant", sdY = sdY[[tr]], N = nG[[tr]], snp = m$snp),
                   list(beta = m$beta, varbeta = m$se^2, type = "quant", sdY = 1, N = nE[[d]], snp = m$snp))
    s <- r$summary
    out <- rbind(out, data.frame(trait = tr, dataset = d, label = lab[[d]], nsnps = s[["nsnps"]],
                                 PP.H0 = s[["PP.H0.abf"]], PP.H1 = s[["PP.H1.abf"]], PP.H2 = s[["PP.H2.abf"]],
                                 PP.H3 = s[["PP.H3.abf"]], PP.H4 = s[["PP.H4.abf"]],
                                 top_shared = r$results$snp[which.max(r$results$SNP.PP.H4)],
                                 top_shared_pp = max(r$results$SNP.PP.H4)))
  }
}
write.csv(out, file.path(B, "results/p2_itga9_coloc.csv"), row.names = FALSE)
print(out[, c("trait", "label", "nsnps", "PP.H3", "PP.H4", "top_shared", "top_shared_pp")], digits = 3)
