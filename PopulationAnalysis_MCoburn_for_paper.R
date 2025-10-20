library(Seurat)
library(ggplot2)
library(multcomp)
library(emmeans)
library(scDC)
library(dplyr)

# Setup

set.seed(20250929)

scRNA <- readRDS("~/RNA_Seq/MCoburn_TBIComparison/PS195X_AllMerge/050120_PS195X_PS195X_AllMerge_Seurat_Object.rds")
scRNA <- DietSeurat(scRNA, counts = FALSE, data = TRUE, scale.data = TRUE, 
                    assays=c(DefaultAssay(scRNA)),
                    dimreducs = NULL)

metadata <- data.frame(Sample = scRNA$orig.ident, 
                       Cluster = scRNA$Clusts)
table(metadata$Sample)
metadata$Genotype <- relevel(factor(sub("_.*", "", metadata$Sample)), ref="WT")
metadata$Cluster <- relevel(metadata$Cluster, ref="Homeostatic")

exprsMat <- GetAssayData(scRNA, slot="scale.data")
rm(scRNA)

# Run scDC without clustering

res.nc = scDC_noClustering(cellTypes = metadata$Cluster,
                           subject = metadata$Sample,
                           calCI = TRUE,
                           calCI_method = c("multinom", "percentile"),
                           nboot = 10000,
                           conf_level = 0.95,
                           ncores = 4,
                           verbose = TRUE)

res.nc$results$genotype <- relevel(factor(sub("_.*", "", res.nc$results$subject)), ref="WT")
res.nc$info$genotype <- relevel(factor(sub("_.*", "", res.nc$info$subject)), ref="WT")
res.nc$results$cellTypes <- relevel(res.nc$results$cellTypes, ref="Homeostatic")
res.nc$info$cellTypes <- relevel(res.nc$info$cellTypes, ref="Homeostatic")

saveRDS(res.nc, '~/RNA_Seq/MCoburn_TBIComparison/noclust_results_2025.rds')

# Load pre-bootstrapped no-clustered data

res.nc <- readRDS('~/RNA_Seq/MCoburn_TBIComparison/noclust_results_2025.rds')

# Normalize cell count across samples

res2.nc <- res.nc
res2.nc$nstar <- round(res2.nc$thetastar * min(table(metadata$Sample))) + 1

fit.nc <- fitGLM(res2.nc, res2.nc$info$genotype, subject_effect = FALSE, pairwise = FALSE, 
                 fixed_only = TRUE, verbose = TRUE)

saveRDS(fit.nc, '~/RNA_Seq/MCoburn_TBIComparison/glm_noclust_2025.rds')

fit.nc <- readRDS('~/RNA_Seq/MCoburn_TBIComparison/glm_noclust_2025.rds')

summary(fit.nc$pool_res_fixed)
subset(summary(fit.nc$pool_res_fixed), p.value < 0.05)


# Pairwise comparisons

all_conts <- lapply(fit.nc$fit_fixed, function(i) {
  glht(i, linfct = emm(pairwise~cond|cellTypes, adjust="tukey"))
})

all_tests <- lapply(1:length(all_conts), function(i) {
  ls <- lapply(all_conts[[i]], function(x) {
    data <- summary(x)$test
  })
  names(ls) <- sub("cellTypes = ", "", names(all_conts[[i]]))
  ls
})

cols <- names(all_tests[[1]]$Homeostatic$coefficients)
clust_ests <- lapply(names(all_tests[[1]]), function(i) {
  data <- do.call(rbind, lapply(all_tests, function(a) {
    a[[i]]$coefficients[cols]
  }))
})
names(clust_ests) <- names(all_tests[[1]])

cols <- names(all_tests[[1]]$Homeostatic$sigma)
clust_ses <- lapply(names(all_tests[[1]]), function(i) {
  data <- do.call(rbind, lapply(all_tests, function(a) {
    a[[i]]$sigma[cols]
  }))
})
names(clust_ses) <- names(all_tests[[1]])

tmp <- lapply(names(all_tests[[1]]), function(i) {
  data <- do.call(rbind, lapply(all_tests, function(a) {
    a[[i]]$pvalues
  }))
})
names(tmp) <- names(all_tests[[1]])

sig <- lapply(names(clust_ests), function(C) {
  data <- do.call(rbind, lapply(1:ncol(clust_ests[[C]]), function(P) {
    pooled.s <- mice::pool.scalar(clust_ests[[C]][,P], clust_ses[[C]][,P]**2, n=48)
    se <- sqrt(pooled.s$ubar + pooled.s$b + pooled.s$b / pooled.s$m)
    pooled.s <- t(cbind(pooled.s[4:length(pooled.s)]))
    rownames(pooled.s) <- c(colnames(clust_ests[[C]])[P])
    pooled.s
  }))
  
  data <- as.data.frame(data)
  data$p.val <- (1 - ptukey(sqrt(2)*abs(unlist(data$qbar) / sqrt(unlist(data$t))), nmeans=32, df=47))
  data$sig <- ' '
  data$sig[data$p.val < 0.05] <- "*"
  data
})

names(sig) <- names(clust_ests)
sig

sig.df <- do.call(rbind, lapply(sig, function(S) { t(S$p.val) }))
rownames(sig.df) <- names(sig)
colnames(sig.df) <- rownames(sig[[1]])

write.csv(t(sig.df), file='~/RNA_Seq/MCoburn_TBIComparison/sig_results_noclust.csv')

saveRDS(list(sig = sig, all_conts = all_conts, clust_ests = clust_ests, 
             clust_ses = clust_ses), '~/RNA_Seq/MCoburn_TBIComparison/sig_results_noclust.rds')
