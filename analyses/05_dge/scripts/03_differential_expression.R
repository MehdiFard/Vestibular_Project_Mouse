# 03 - Differential expression with DESeq2.
#
# Five analyses:
#   A  Utricle vs Saccule at w8      (paired within strain)
#   B  Utricle vs Saccule at w22     (paired within strain)
#   C  w22 vs w8 in Saccule          (batch-adjusted)
#   D  w22 vs w8 in Utricle          (batch-adjusted)
#   E  organ x age interaction       (does the organ difference change with age?)
#
# Sign convention: positive log2FC means higher in Utricle (A, B) or higher at
# w22 (C, D).
#
# Why paired for the organ contrasts: saccule and utricle were taken from the
# same animal, so strain identifies the pair. Putting strain in the model
# removes between-animal variation and the genetic background of this BXD-style
# panel, which is by far the largest nuisance source here.

source("../../analyses/05_dge/scripts/00_config.R")

library(data.table)
library(DESeq2)
library(BiocParallel)

# DESeq2 fits one small GLM per gene. Letting the BLAS spread each of those
# across all 128 cores is far slower than fitting many genes at once, so pin
# the BLAS to one thread and parallelise over genes instead.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
BPP <- MulticoreParam(workers = 16)
register(BPP)

samples <- fread(file.path(DATA, "sample_table.tsv"))
txi     <- readRDS(file.path(DATA, "txi_gene.rds"))
anno    <- fread(file.path(DATA, "gene_annotation.tsv"))

samples[, organ := factor(organ, levels = c("Saccule", "Utricle"))]
samples[, age   := factor(age,   levels = c("w8", "w22"))]
samples[, batch := factor(batch)]

stopifnot(identical(colnames(txi$counts), samples$library_id))

# ---------------------------------------------------------------- helpers

# Subset the tximport object to a set of libraries.
subset_txi <- function(txi, ids) {
  idx <- match(ids, colnames(txi$counts))
  list(counts             = txi$counts[, idx, drop = FALSE],
       abundance          = txi$abundance[, idx, drop = FALSE],
       length             = txi$length[, idx, drop = FALSE],
       countsFromAbundance = txi$countsFromAbundance)
}

# Two filters, applied in this order:
#   1. biotype - only protein-coding and lincRNA genes are analysed
#   2. expression - >= 10 counts in >= (size of the smaller group) libraries
# Dropping the other biotypes first keeps them out of the multiple-testing
# burden, which makes the FDR correction less conservative for the genes we
# actually care about.
KEEP_IDS <- anno[biotype %in% KEEP_BIOTYPES, gene_id]

filter_genes <- function(dds, group_sizes) {
  dds  <- dds[rownames(dds) %in% KEEP_IDS, ]
  keep <- rowSums(counts(dds) >= 10) >= min(group_sizes)
  dds[keep, ]
}

# Run one contrast and return an annotated, ordered results table.
run_contrast <- function(dds, contrast, label) {
  cat("\n=====", label, "=====\n")
  dds <- DESeq(dds, quiet = TRUE, parallel = TRUE, BPPARAM = BPP)
  res <- results(dds, contrast = contrast, alpha = FDR_CUT)
  # Shrink the fold changes so the ranking is not dominated by noisy low-count genes.
  res_s <- lfcShrink(dds, contrast = contrast, res = res, type = "ashr", quiet = TRUE,
                     parallel = TRUE, BPPARAM = BPP)

  out <- data.table(
    gene_id        = rownames(res),
    baseMean       = res$baseMean,
    log2FoldChange = res_s$log2FoldChange,   # shrunken
    lfcSE          = res_s$lfcSE,
    log2FC_raw     = res$log2FoldChange,     # unshrunken, for reference
    stat           = res$stat,
    pvalue         = res$pvalue,
    padj           = res$padj)

  out <- merge(out, anno[, .(gene_id, gene_name, biotype, chr)],
               by = "gene_id", all.x = TRUE, sort = FALSE)
  # setorder cannot take an expression, so sort on an explicit helper column.
  out[, abs_lfc := abs(log2FoldChange)]
  setorder(out, padj, -abs_lfc, na.last = TRUE)
  out[, abs_lfc := NULL]

  sig <- out[!is.na(padj) & padj < FDR_CUT & abs(log2FoldChange) > LFC_CUT]
  cat("Genes tested:", nrow(out), "\n")
  cat("FDR <", FDR_CUT, ":", sum(out$padj < FDR_CUT, na.rm = TRUE), "\n")
  cat("FDR <", FDR_CUT, "and |log2FC| >", LFC_CUT, ":", nrow(sig),
      " (up:", sum(sig$log2FoldChange > 0), " down:", sum(sig$log2FoldChange < 0), ")\n")

  list(dds = dds, res = out)
}

results_list <- list()
summary_rows <- list()

record <- function(key, label, model, res, n_samples, note = "") {
  fwrite(res, file.path(TABLES, paste0("DE_", key, ".tsv.gz")), sep = "\t", na = "NA")
  sig <- res[!is.na(padj) & padj < FDR_CUT & abs(log2FoldChange) > LFC_CUT]
  fwrite(sig, file.path(TABLES, paste0("DE_", key, "_significant.tsv")), sep = "\t", na = "NA")
  results_list[[key]] <<- res
  summary_rows[[key]] <<- data.table(
    analysis = key, description = label, model = model, n_libraries = n_samples,
    genes_tested = nrow(res),
    sig_fdr      = sum(res$padj < FDR_CUT, na.rm = TRUE),
    sig_fdr_lfc  = nrow(sig),
    up           = sum(sig$log2FoldChange > 0),
    down         = sum(sig$log2FoldChange < 0),
    note         = note)
}

# ================================================================ A & B
# Utricle vs Saccule, within each age, paired by strain.
for (a in c("w8", "w22")) {
  # A pair is one animal: same strain, same age, same sex, both organs.
  s_age  <- samples[age == a]
  paired <- s_age[, .(n = .N, k = uniqueN(organ)), by = pair_id][n == 2 & k == 2, pair_id]
  s_use  <- s_age[pair_id %in% paired]
  s_use[, strain := factor(pair_id)]

  cat("\n[", a, "] animal pairs:", length(paired),
      "-> libraries:", nrow(s_use), "\n")

  dds <- DESeqDataSetFromTximport(
    subset_txi(txi, s_use$library_id),
    colData = data.frame(s_use, row.names = s_use$library_id),
    design  = ~ strain + organ)
  dds <- filter_genes(dds, table(s_use$organ))

  key <- paste0("A_utricle_vs_saccule_", a)
  if (a == "w22") key <- paste0("B_utricle_vs_saccule_", a)
  fit <- run_contrast(dds, c("organ", "Utricle", "Saccule"),
                      paste("Utricle vs Saccule at", a, "(paired by strain)"))
  record(key, paste("Utricle vs Saccule at", a), "~ strain + organ",
         fit$res, nrow(s_use),
         paste0(length(paired), " strains with both organs"))
}

# ================================================================ C & D
# w22 vs w8, within each organ. These cannot be paired across most strains
# (only ~25 strains were sampled at both ages), so batch is used as a covariate
# on the full sample set. Batch and age are partly confounded: batch 1 is
# w8-only and batch 4 is w22-only, so the age effect is identified mainly from
# batches 2 and 3.
for (o in c("Saccule", "Utricle")) {
  s_use <- samples[organ == o]
  s_use[, batch := droplevels(factor(batch))]

  cat("\n[", o, "] libraries:", nrow(s_use), "\n")
  print(table(s_use$batch, s_use$age))

  dds <- DESeqDataSetFromTximport(
    subset_txi(txi, s_use$library_id),
    colData = data.frame(s_use, row.names = s_use$library_id),
    design  = ~ batch + age)
  dds <- filter_genes(dds, table(s_use$age))

  key <- paste0(ifelse(o == "Saccule", "C", "D"), "_w22_vs_w8_", tolower(o))
  fit <- run_contrast(dds, c("age", "w22", "w8"),
                      paste("w22 vs w8 in", o, "(batch-adjusted)"))
  record(key, paste("w22 vs w8 in", o), "~ batch + age", fit$res, nrow(s_use),
         "batch and age partly confounded")
}

# ================================================================ E
# Organ x age interaction on all libraries: is the Utricle-vs-Saccule
# difference itself different between w8 and w22?
cat("\n===== organ x age interaction =====\n")
s_all <- copy(samples)
dds_int <- DESeqDataSetFromTximport(
  subset_txi(txi, s_all$library_id),
  colData = data.frame(s_all, row.names = s_all$library_id),
  design  = ~ batch + organ + age + organ:age)
dds_int <- filter_genes(dds_int, table(s_all$group))
dds_int <- DESeq(dds_int, quiet = TRUE, parallel = TRUE, BPPARAM = BPP)
cat("Coefficients:", paste(resultsNames(dds_int), collapse = ", "), "\n")

res_int <- results(dds_int, name = "organUtricle.agew22", alpha = FDR_CUT)
out_int <- data.table(
  gene_id        = rownames(res_int),
  baseMean       = res_int$baseMean,
  log2FoldChange = res_int$log2FoldChange,
  lfcSE          = res_int$lfcSE,
  log2FC_raw     = res_int$log2FoldChange,
  stat           = res_int$stat,
  pvalue         = res_int$pvalue,
  padj           = res_int$padj)
out_int <- merge(out_int, anno[, .(gene_id, gene_name, biotype, chr)],
                 by = "gene_id", all.x = TRUE, sort = FALSE)
setorder(out_int, padj, na.last = TRUE)

cat("Interaction FDR <", FDR_CUT, ":", sum(out_int$padj < FDR_CUT, na.rm = TRUE), "\n")
record("E_organ_by_age_interaction",
       "Organ x age interaction (does the organ difference change with age)",
       "~ batch + organ + age + organ:age", out_int, nrow(s_all),
       "log2FC = change in the Utricle-Saccule difference from w8 to w22")

saveRDS(dds_int, file.path(DATA, "dds_interaction.rds"))

# ================================================================ sensitivity
# Unpaired organ contrasts, i.e. the model without strain, so the effect of
# pairing can be seen directly.
for (a in c("w8", "w22")) {
  s_use <- samples[age == a]
  s_use[, batch := droplevels(factor(batch))]
  dds <- DESeqDataSetFromTximport(
    subset_txi(txi, s_use$library_id),
    colData = data.frame(s_use, row.names = s_use$library_id),
    design  = ~ batch + organ)
  dds <- filter_genes(dds, table(s_use$organ))
  fit <- run_contrast(dds, c("organ", "Utricle", "Saccule"),
                      paste("Utricle vs Saccule at", a, "(UNPAIRED, sensitivity)"))
  record(paste0("S_unpaired_utricle_vs_saccule_", a),
         paste("Utricle vs Saccule at", a, "- unpaired sensitivity check"),
         "~ batch + organ", fit$res, nrow(s_use), "sensitivity analysis")
}

# ================================================================ summary
summary_dt <- rbindlist(summary_rows)
fwrite(summary_dt, file.path(TABLES, "DE_summary.tsv"), sep = "\t")
saveRDS(results_list, file.path(DATA, "de_results.rds"))

cat("\n\n================ SUMMARY ================\n")
print(summary_dt[, .(analysis, n_libraries, genes_tested, sig_fdr, sig_fdr_lfc, up, down)])
cat("\nDone.\n")
