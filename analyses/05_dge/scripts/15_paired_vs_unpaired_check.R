# 15 - Why does the UNPAIRED organ contrast report more genes than the paired one?
#
# The headline count in DE_summary.tsv is "FDR < 0.05 AND |log2FC| > 1". For the
# organ contrasts the FDR half of that is not binding at all - every gene that
# clears |log2FC| > 1 also clears FDR, in both models. So the headline number is
# counting effect sizes, not significance, and a noisier estimate produces MORE
# of them: noise widens the spread of the fold-change estimates, which pushes
# extra genes past the +/-1 line from both sides.
#
# Two things differ between the paired and unpaired runs in 03: the model
# (strain in or out) and the sample set (the unpaired one also keeps libraries
# whose partner organ is missing). This script separates them by refitting both
# models on exactly the same libraries.

source("../../analyses/05_dge/scripts/00_config.R")

library(data.table)
library(DESeq2)
library(BiocParallel)

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
BPP <- MulticoreParam(workers = 16)
register(BPP)

samples <- fread(file.path(DATA, "sample_table.tsv"))
txi     <- readRDS(file.path(DATA, "txi_gene.rds"))
anno    <- fread(file.path(DATA, "gene_annotation.tsv"))

samples[, organ := factor(organ, levels = c("Saccule", "Utricle"))]
samples[, age   := factor(age,   levels = c("w8", "w22"))]
samples[, batch := factor(batch)]

KEEP_IDS <- anno[biotype %in% KEEP_BIOTYPES, gene_id]

subset_txi <- function(txi, ids) {
  idx <- match(ids, colnames(txi$counts))
  list(counts = txi$counts[, idx, drop = FALSE],
       abundance = txi$abundance[, idx, drop = FALSE],
       length = txi$length[, idx, drop = FALSE],
       countsFromAbundance = txi$countsFromAbundance)
}

# Fit one organ contrast and return the shrunken result table.
fit_organ <- function(s_use, design, label) {
  cat("\n---", label, "- libraries:", nrow(s_use), "\n")
  dds <- DESeqDataSetFromTximport(
    subset_txi(txi, s_use$library_id),
    colData = data.frame(s_use, row.names = s_use$library_id),
    design  = design)
  dds  <- dds[rownames(dds) %in% KEEP_IDS, ]
  keep <- rowSums(counts(dds) >= 10) >= min(table(s_use$organ))
  dds  <- dds[keep, ]
  dds  <- DESeq(dds, quiet = TRUE, parallel = TRUE, BPPARAM = BPP)
  res  <- results(dds, contrast = c("organ", "Utricle", "Saccule"), alpha = FDR_CUT)
  shr  <- lfcShrink(dds, contrast = c("organ", "Utricle", "Saccule"), res = res,
                    type = "ashr", quiet = TRUE, parallel = TRUE, BPPARAM = BPP)
  data.table(gene_id = rownames(res), baseMean = res$baseMean,
             lfc = shr$log2FoldChange, lfc_raw = res$log2FoldChange,
             se = res$lfcSE, padj = res$padj,
             resid_df = ncol(dds) - qr(model.matrix(design, data = colData(dds)))$rank)
}

summarise <- function(d, label, n_lib, model) {
  data.table(comparison = label, model = model, n_libraries = n_lib,
             residual_df = d$resid_df[1],
             genes_tested = nrow(d),
             sig_fdr = sum(d$padj < FDR_CUT, na.rm = TRUE),
             lfc_over_1 = sum(abs(d$lfc) > 1, na.rm = TRUE),
             both = sum(d$padj < FDR_CUT & abs(d$lfc) > 1, na.rm = TRUE),
             median_se = round(median(d$se, na.rm = TRUE), 4),
             sd_lfc = round(sd(d$lfc_raw, na.rm = TRUE), 3))
}

rows <- list()
for (a in c("w8", "w22")) {
  s_age  <- samples[age == a]
  paired <- s_age[, .(n = .N, k = uniqueN(organ)), by = pair_id][n == 2 & k == 2, pair_id]
  s_pair <- s_age[pair_id %in% paired]
  s_pair[, strain_f := factor(pair_id)]

  # 1. the paired model as run in 03
  d1 <- fit_organ(s_pair, ~ strain_f + organ, paste(a, "paired, ~ strain + organ"))
  # 2. SAME libraries, strain dropped - isolates the model
  d2 <- fit_organ(s_pair, ~ batch + organ, paste(a, "same libraries, ~ batch + organ"))
  # 3. all libraries, strain dropped - the sensitivity run from 03
  s_all <- copy(s_age)[, batch := droplevels(batch)]
  d3 <- fit_organ(s_all, ~ batch + organ, paste(a, "all libraries, ~ batch + organ"))

  rows <- c(rows, list(
    summarise(d1, a, nrow(s_pair), "~ strain + organ (paired)"),
    summarise(d2, a, nrow(s_pair), "~ batch + organ (same libraries)"),
    summarise(d3, a, nrow(s_all),  "~ batch + organ (all libraries)")))

  # Where do the extra genes come from? Compare the same-library fits gene by gene.
  m <- merge(d1[, .(gene_id, lfc_paired = lfc, padj_paired = padj, se_paired = se)],
             d2[, .(gene_id, lfc_unpair = lfc, padj_unpair = padj, se_unpair = se)],
             by = "gene_id")
  cat("\n[", a, "] same libraries, paired vs strain-free:\n")
  cat("  genes over |lfc|>1 in both:", m[abs(lfc_paired) > 1 & abs(lfc_unpair) > 1, .N],
      "| only strain-free:", m[abs(lfc_paired) <= 1 & abs(lfc_unpair) > 1, .N],
      "| only paired:", m[abs(lfc_paired) > 1 & abs(lfc_unpair) <= 1, .N], "\n")
  cat("  SE ratio (strain-free / paired), median:",
      round(median(m$se_unpair / m$se_paired, na.rm = TRUE), 3), "\n")
  cat("  among genes near the threshold (0.7 < |lfc_paired| < 1.3), median SE ratio:",
      round(median(m[abs(lfc_paired) > 0.7 & abs(lfc_paired) < 1.3, se_unpair / se_paired],
                   na.rm = TRUE), 3), "\n")
  fwrite(m, file.path(TABLES, paste0("paired_vs_unpaired_genes_", a, ".tsv.gz")), sep = "\t")
}

out <- rbindlist(rows)
fwrite(out, file.path(TABLES, "paired_vs_unpaired_summary.tsv"), sep = "\t")
cat("\n\n================ SUMMARY ================\n")
print(out)
cat("\nDone.\n")
