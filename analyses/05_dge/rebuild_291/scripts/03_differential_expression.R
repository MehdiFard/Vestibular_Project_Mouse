# 03 - Differential expression on the 291-library set.
# Same five contrasts as v1. The only change is that the blocking factor for the
# organ contrasts is the ANIMAL from the submission record, not strain+age+sex.
source("../../analyses/05_dge/rebuild_291/scripts/00_config.R")
library(data.table); library(DESeq2); library(BiocParallel)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
BPP <- MulticoreParam(workers = 16); register(BPP)

samples <- fread(file.path(DATA, "sample_table.tsv"))
txi <- readRDS(file.path(DATA, "txi_gene.rds"))
anno <- fread(file.path(DATA, "gene_annotation.tsv"))
samples[, organ := factor(organ, levels = c("Saccule", "Utricle"))]
samples[, age := factor(age, levels = c("w8", "w22"))]
samples[, batch := factor(batch)]
stopifnot(identical(colnames(txi$counts), samples$library_id))

subset_txi <- function(txi, ids) {
  i <- match(ids, colnames(txi$counts))
  list(counts = txi$counts[, i, drop = FALSE], abundance = txi$abundance[, i, drop = FALSE],
       length = txi$length[, i, drop = FALSE], countsFromAbundance = txi$countsFromAbundance)
}
KEEP_IDS <- anno[biotype %in% KEEP_BIOTYPES, gene_id]
filter_genes <- function(dds, sizes) {
  dds <- dds[rownames(dds) %in% KEEP_IDS, ]
  dds[rowSums(counts(dds) >= 10) >= min(sizes), ]
}
run <- function(dds, contrast, label) {
  cat("\n=====", label, "=====\n")
  dds <- DESeq(dds, quiet = TRUE, parallel = TRUE, BPPARAM = BPP)
  res <- results(dds, contrast = contrast, alpha = FDR_CUT)
  shr <- lfcShrink(dds, contrast = contrast, res = res, type = "ashr", quiet = TRUE,
                   parallel = TRUE, BPPARAM = BPP)
  out <- data.table(gene_id = rownames(res), baseMean = res$baseMean,
                    log2FoldChange = shr$log2FoldChange, lfcSE = shr$lfcSE,
                    log2FC_raw = res$log2FoldChange, stat = res$stat,
                    pvalue = res$pvalue, padj = res$padj)
  out <- merge(out, anno[, .(gene_id, gene_name, biotype, chr)], by = "gene_id", all.x = TRUE, sort = FALSE)
  out[, abs_lfc := abs(log2FoldChange)]; setorder(out, padj, -abs_lfc, na.last = TRUE); out[, abs_lfc := NULL]
  sig <- out[!is.na(padj) & padj < FDR_CUT & abs(log2FoldChange) > LFC_CUT]
  cat("tested", nrow(out), "| FDR<0.05:", sum(out$padj < FDR_CUT, na.rm = TRUE),
      "| + |lfc|>1:", nrow(sig), "(up", sum(sig$log2FoldChange > 0), "down", sum(sig$log2FoldChange < 0), ")\n")
  list(dds = dds, res = out)
}
results_list <- list(); summary_rows <- list()
record <- function(key, label, model, res, n, note = "") {
  fwrite(res, file.path(TABLES, paste0("DE_", key, ".tsv.gz")), sep = "\t", na = "NA")
  sig <- res[!is.na(padj) & padj < FDR_CUT & abs(log2FoldChange) > LFC_CUT]
  fwrite(sig, file.path(TABLES, paste0("DE_", key, "_significant.tsv")), sep = "\t", na = "NA")
  results_list[[key]] <<- res
  summary_rows[[key]] <<- data.table(analysis = key, description = label, model = model,
    n_libraries = n, genes_tested = nrow(res), sig_fdr = sum(res$padj < FDR_CUT, na.rm = TRUE),
    sig_fdr_lfc = nrow(sig), up = sum(sig$log2FoldChange > 0), down = sum(sig$log2FoldChange < 0), note = note)
}

# ---- A & B: organ, blocked on animal
for (a in c("w8", "w22")) {
  s <- samples[age == a]
  paired <- s[, .(n = .N, k = uniqueN(organ)), by = pair_id][n == 2 & k == 2, pair_id]
  s <- s[pair_id %in% paired]; s[, animal := factor(pair_id)]
  cat("\n[", a, "] animals:", length(paired), "-> libraries:", nrow(s), "\n")
  dds <- DESeqDataSetFromTximport(subset_txi(txi, s$library_id),
           colData = data.frame(s, row.names = s$library_id), design = ~ animal + organ)
  dds <- filter_genes(dds, table(s$organ))
  key <- paste0(ifelse(a == "w8", "A", "B"), "_utricle_vs_saccule_", a)
  fit <- run(dds, c("organ", "Utricle", "Saccule"), paste("Utricle vs Saccule at", a))
  record(key, paste("Utricle vs Saccule at", a), "~ animal + organ", fit$res, nrow(s),
         paste(length(paired), "animals with both organs"))
}
# ---- C & D: age within organ
for (o in c("Saccule", "Utricle")) {
  s <- samples[organ == o]; s[, batch := droplevels(factor(batch))]
  cat("\n[", o, "] libraries:", nrow(s), "\n"); print(table(s$batch, s$age))
  dds <- DESeqDataSetFromTximport(subset_txi(txi, s$library_id),
           colData = data.frame(s, row.names = s$library_id), design = ~ batch + age)
  dds <- filter_genes(dds, table(s$age))
  key <- paste0(ifelse(o == "Saccule", "C", "D"), "_w22_vs_w8_", tolower(o))
  fit <- run(dds, c("age", "w22", "w8"), paste("w22 vs w8 in", o))
  record(key, paste("w22 vs w8 in", o), "~ batch + age", fit$res, nrow(s), "batch and age partly confounded")
}
# ---- E: interaction
cat("\n===== interaction =====\n")
s <- copy(samples)
dds <- DESeqDataSetFromTximport(subset_txi(txi, s$library_id),
         colData = data.frame(s, row.names = s$library_id), design = ~ batch + organ + age + organ:age)
dds <- filter_genes(dds, table(paste(s$organ, s$age)))
dds <- DESeq(dds, quiet = TRUE, parallel = TRUE, BPPARAM = BPP)
res <- results(dds, name = "organUtricle.agew22", alpha = FDR_CUT)
out <- data.table(gene_id = rownames(res), baseMean = res$baseMean,
                  log2FoldChange = res$log2FoldChange, lfcSE = res$lfcSE,
                  log2FC_raw = res$log2FoldChange, stat = res$stat,
                  pvalue = res$pvalue, padj = res$padj)
out <- merge(out, anno[, .(gene_id, gene_name, biotype, chr)], by = "gene_id", all.x = TRUE, sort = FALSE)
setorder(out, padj, na.last = TRUE)
cat("interaction FDR<0.05:", sum(out$padj < FDR_CUT, na.rm = TRUE), "\n")
record("E_organ_by_age_interaction", "Organ x age interaction",
       "~ batch + organ + age + organ:age", out, nrow(s), "change in the organ gap from w8 to w22")
saveRDS(dds, file.path(DATA, "dds_interaction.rds"))

summary_dt <- rbindlist(summary_rows)
fwrite(summary_dt, file.path(TABLES, "DE_summary.tsv"), sep = "\t")
saveRDS(results_list, file.path(DATA, "de_results.rds"))
cat("\n\n============ SUMMARY ============\n")
print(summary_dt[, .(analysis, n_libraries, genes_tested, sig_fdr, sig_fdr_lfc, up, down)])
