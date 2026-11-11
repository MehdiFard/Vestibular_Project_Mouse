# 06 - Why the DE models keep batch in the design instead of testing
#      limma::removeBatchEffect output.
#
# removeBatchEffect subtracts the fitted batch effect but does not tell the
# downstream test that batch parameters were estimated. Two things follow, and
# the tests below show they matter very differently here:
#
#   - lost degrees of freedom. Real, but negligible at this sample size
#     (123 vs 126 residual df on 128 libraries).
#
#   - lost collinearity. This is the one that bites. When batch overlaps the
#     factor of interest, the correct model inflates the standard error of that
#     factor to reflect how little independent information there is. Testing
#     "~ group" on a corrected matrix has no batch term left to be collinear
#     with, so the standard error comes out far too small and the effect looks
#     much more certain than the data support.
#
# Three checks, on the real saccule data:
#   1. the real w22-vs-w8 contrast, correct model vs corrected-matrix model
#   2. a null label balanced within batch (orthogonal - the harmless case)
#   3. a null label confounded with batch (like age really is here)
# In 2 and 3 the truth is "no difference", so every rejection is a false positive.

source("../../analyses/05_dge/scripts/00_config.R")

library(data.table)
library(DESeq2)
library(limma)
library(edgeR)

set.seed(1)

samples <- fread(file.path(DATA, "sample_table.tsv"))
samples[, organ := factor(organ, levels = c("Saccule", "Utricle"))]
samples[, age   := factor(age,   levels = c("w8", "w22"))]
samples[, batch := factor(batch)]
txi  <- readRDS(file.path(DATA, "txi_gene.rds"))
anno <- fread(file.path(DATA, "gene_annotation.tsv"))

# Same saccule subset and same filters as contrast C.
s <- samples[organ == "Saccule"]
s[, batch := droplevels(batch)]
idx <- match(s$library_id, colnames(txi$counts))

counts <- txi$counts[, idx]
counts <- counts[rownames(counts) %in% anno[biotype %in% KEEP_BIOTYPES, gene_id], ]
counts <- counts[rowSums(counts >= 10) >= min(table(s$age)), ]
cat("Saccule libraries:", ncol(counts), " genes:", nrow(counts), "\n\n")

dge <- calcNormFactors(DGEList(round(counts)))

# ---------------------------------------------------------------- the two ways
# CORRECT: batch stays in the design, so lmFit knows what it cost.
fit_in_model <- function(group) {
  design <- model.matrix(~ s$batch + group)
  v   <- voom(dge, design)
  fit <- eBayes(lmFit(v, design))
  topTable(fit, coef = ncol(design), number = Inf, sort.by = "none")
}

# WRONG: subtract the batch effect first, then test as if batch never existed.
fit_on_corrected <- function(group) {
  design_full <- model.matrix(~ s$batch + group)
  v <- voom(dge, design_full)
  v$E <- removeBatchEffect(v$E, batch = s$batch, design = model.matrix(~ group))
  design_naive <- model.matrix(~ group)
  fit <- eBayes(lmFit(v, design_naive))
  topTable(fit, coef = 2, number = Inf, sort.by = "none")
}

# ---------------------------------------------------------------- 1. real contrast
cat("=== Real contrast: w22 vs w8 in Saccule ===\n")
real_ok  <- fit_in_model(s$age)
real_bad <- fit_on_corrected(s$age)

cmp <- data.table(
  approach = c("batch in the design (~ batch + age)",
               "removeBatchEffect, then test ~ age"),
  residual_df = c(ncol(counts) - length(unique(s$batch)) - 1,
                  ncol(counts) - 2),
  sig_fdr05  = c(sum(real_ok$adj.P.Val  < 0.05), sum(real_bad$adj.P.Val  < 0.05)),
  sig_p001   = c(sum(real_ok$P.Value    < 0.001), sum(real_bad$P.Value   < 0.001)))
print(cmp)
cat("\nInflation at FDR < 0.05:",
    sprintf("%.1f%% more genes called\n",
            100 * (cmp$sig_fdr05[2] / cmp$sig_fdr05[1] - 1)))

# ---------------------------------------------------------------- 2. null test
# A random label, balanced inside every batch, so the truth is "no difference".
# Anything called significant here is a false positive.
cat("\n=== Null test: random label balanced within batch, 20 draws ===\n")
null_run <- function(i) {
  fake <- factor(unlist(lapply(split(seq_len(nrow(s)), s$batch), function(ii) {
    x <- rep(c("a", "b"), length.out = length(ii)); sample(x)
  }))[order(unlist(split(seq_len(nrow(s)), s$batch)))])

  ok  <- fit_in_model(fake)
  bad <- fit_on_corrected(fake)
  data.table(draw = i,
             ok_fdr05  = sum(ok$adj.P.Val  < 0.05),
             bad_fdr05 = sum(bad$adj.P.Val < 0.05),
             ok_p05    = mean(ok$P.Value  < 0.05),
             bad_p05   = mean(bad$P.Value < 0.05))
}

null_res <- rbindlist(lapply(1:20, null_run))

cat("\nFalse positives at FDR < 0.05 (should be ~0):\n")
cat("  batch in the design    : mean", sprintf("%6.1f", mean(null_res$ok_fdr05)),
    " max", max(null_res$ok_fdr05), "\n")
cat("  removeBatchEffect first: mean", sprintf("%6.1f", mean(null_res$bad_fdr05)),
    " max", max(null_res$bad_fdr05), "\n")

cat("\nFraction of genes with raw p < 0.05 (should be ~0.05):\n")
cat("  batch in the design    :", sprintf("%.3f", mean(null_res$ok_p05)), "\n")
cat("  removeBatchEffect first:", sprintf("%.3f", mean(null_res$bad_p05)), "\n")

# ---------------------------------------------------------------- 3. confounded null
# The test above used a label balanced inside every batch, i.e. orthogonal to
# batch - the easy case. Age in this study is NOT like that: batch 1 is w8-only
# and batch 4 is w22-only. So repeat the null with a random label that is
# deliberately correlated with batch, matching the real confounding.
cat("\n=== Null test 2: random label CONFOUNDED with batch, 20 draws ===\n")
cat("(mimics age here: batch 1 is w8-only, batch 4 is w22-only)\n")

null_conf <- function(i) {
  # 85% of one batch gets label "a", 15% of another - a strong batch/label overlap
  p_a  <- c(0.85, 0.6, 0.4, 0.15)[as.integer(s$batch)]
  fake <- factor(ifelse(runif(nrow(s)) < p_a, "a", "b"))
  if (length(unique(fake)) < 2) return(NULL)

  ok  <- fit_in_model(fake)
  bad <- fit_on_corrected(fake)
  data.table(draw = i,
             ok_fdr05  = sum(ok$adj.P.Val  < 0.05),
             bad_fdr05 = sum(bad$adj.P.Val < 0.05))
}
conf_res <- rbindlist(lapply(1:20, null_conf))

cat("\nFalse positives at FDR < 0.05 (should be ~0):\n")
cat("  batch in the design    : mean", sprintf("%7.1f", mean(conf_res$ok_fdr05)),
    " max", max(conf_res$ok_fdr05), "\n")
cat("  removeBatchEffect first: mean", sprintf("%7.1f", mean(conf_res$bad_fdr05)),
    " max", max(conf_res$bad_fdr05), "\n")

fwrite(conf_res, file.path(TABLES, "batch_correction_null_confounded.tsv"), sep = "\t")
fwrite(null_res, file.path(TABLES, "batch_correction_null_test.tsv"), sep = "\t")
fwrite(cmp, file.path(TABLES, "batch_correction_comparison.tsv"), sep = "\t")

cat("\nWritten to", TABLES, "\n")
