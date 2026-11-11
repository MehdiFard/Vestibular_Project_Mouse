# 14 - Now that sex is known for all 260 libraries, ask what not having it cost.
#
# Three things:
#   1. how many genes respond to sex at all, tested inside w22 where both
#      sexes are present
#   2. a female-only age contrast - w8 is entirely female, so restricting w22
#      to females removes the sex confound from the age comparison completely
#   3. how much the original age hits overlap the sex hits

source("../../analyses/05_dge/scripts/00_config.R")

library(data.table)
library(DESeq2)
library(BiocParallel)
library(ggplot2)

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
BPP <- MulticoreParam(workers = 8); register(BPP)

samples <- fread(file.path(DATA, "sample_table.tsv"))
samples[, organ := factor(organ, levels = c("Saccule", "Utricle"))]
samples[, age   := factor(age,   levels = c("w8", "w22"))]
samples[, batch := factor(batch)]
samples[, sex   := factor(sex_final, levels = c("F", "M"))]

txi  <- readRDS(file.path(DATA, "txi_gene.rds"))
anno <- fread(file.path(DATA, "gene_annotation.tsv"))
de   <- readRDS(file.path(DATA, "de_results.rds"))
KEEP_IDS <- anno[biotype %in% KEEP_BIOTYPES, gene_id]

cat("=== sex by organ and age, after recovery ===\n")
print(samples[, .N, by = .(organ, age, sex)][order(organ, age, sex)])

subset_txi <- function(ids) {
  i <- match(ids, colnames(txi$counts))
  list(counts = txi$counts[, i, drop = FALSE], abundance = txi$abundance[, i, drop = FALSE],
       length = txi$length[, i, drop = FALSE], countsFromAbundance = txi$countsFromAbundance)
}

run <- function(s, design, contrast, label, group_sizes) {
  cat("\n=====", label, "=====\n")
  cat("libraries:", nrow(s), "\n")
  dds <- DESeqDataSetFromTximport(
    subset_txi(s$library_id),
    colData = data.frame(s, row.names = s$library_id), design = design)
  dds <- dds[rownames(dds) %in% KEEP_IDS, ]
  dds <- dds[rowSums(counts(dds) >= 10) >= min(group_sizes), ]
  dds <- DESeq(dds, quiet = TRUE, parallel = TRUE, BPPARAM = BPP)
  res <- results(dds, contrast = contrast, alpha = FDR_CUT)
  out <- data.table(gene_id = rownames(res), baseMean = res$baseMean,
                    log2FoldChange = res$log2FoldChange, lfcSE = res$lfcSE,
                    pvalue = res$pvalue, padj = res$padj)
  out <- merge(out, anno[, .(gene_id, gene_name, biotype, chr)], by = "gene_id",
               all.x = TRUE, sort = FALSE)
  setorder(out, padj, na.last = TRUE)
  sig <- out[!is.na(padj) & padj < FDR_CUT & abs(log2FoldChange) > LFC_CUT]
  cat("tested:", nrow(out), " FDR<0.05:", sum(out$padj < FDR_CUT, na.rm = TRUE),
      " +|LFC|>1:", nrow(sig), "\n")
  out
}

res <- list()

# ------------------------------------------------------------------ 1. sex effect
# Inside w22 only, where both sexes are present.
for (o in c("Saccule", "Utricle")) {
  s <- samples[organ == o & age == "w22"][, batch := droplevels(batch)][]
  res[[paste0("sex_w22_", tolower(o))]] <- run(
    s, ~ batch + sex, c("sex", "M", "F"),
    paste("Male vs Female, w22", o), table(s$sex))
}

# ------------------------------------------------------------------ 2. female-only age
for (o in c("Saccule", "Utricle")) {
  s <- samples[organ == o & sex == "F"][, batch := droplevels(batch)][]
  res[[paste0("age_females_", tolower(o))]] <- run(
    s, ~ batch + age, c("age", "w22", "w8"),
    paste("w22 vs w8 in", o, "- FEMALES ONLY"), table(s$age))
}

for (k in names(res))
  fwrite(res[[k]], file.path(TABLES, paste0("DE_", k, ".tsv.gz")), sep = "\t", na = "NA")

# ------------------------------------------------------------------ 3. overlap
sig_of <- function(d) d[!is.na(padj) & padj < FDR_CUT & abs(log2FoldChange) > LFC_CUT, gene_id]

cmp <- rbindlist(lapply(c("Saccule", "Utricle"), function(o) {
  orig <- de[[paste0(ifelse(o == "Saccule", "C", "D"), "_w22_vs_w8_", tolower(o))]]
  fem  <- res[[paste0("age_females_", tolower(o))]]
  sx   <- res[[paste0("sex_w22_", tolower(o))]]
  a <- sig_of(orig); f <- sig_of(fem); x <- sig_of(sx)
  m <- merge(orig[, .(gene_id, o_lfc = log2FoldChange)],
             fem[, .(gene_id, f_lfc = log2FoldChange)], by = "gene_id")
  data.table(organ = o,
             age_hits_all_samples = length(a),
             age_hits_females_only = length(f),
             shared = length(intersect(a, f)),
             pct_recovered = round(100 * length(intersect(a, f)) / max(length(a), 1), 1),
             r_lfc = round(cor(m$o_lfc, m$f_lfc, use = "complete.obs"), 3),
             sex_hits_w22 = length(x),
             age_hits_that_are_also_sex_hits = length(intersect(a, x)))
}))
cat("\n\n================ WHAT SEX COST US ================\n")
print(cmp)
fwrite(cmp, file.path(TABLES, "sex_effect_summary.tsv"), sep = "\t")

# ------------------------------------------------------------------ figure
pd <- rbindlist(lapply(c("Saccule", "Utricle"), function(o) {
  orig <- de[[paste0(ifelse(o == "Saccule", "C", "D"), "_w22_vs_w8_", tolower(o))]]
  fem  <- res[[paste0("age_females_", tolower(o))]]
  m <- merge(orig[, .(gene_id, o_lfc = log2FoldChange, o_padj = padj)],
             fem[, .(gene_id, f_lfc = log2FoldChange)], by = "gene_id")
  m <- m[!is.na(o_lfc) & !is.na(f_lfc)]
  m[, hit := !is.na(o_padj) & o_padj < FDR_CUT & abs(o_lfc) > LFC_CUT]
  m[, lab := sprintf("%s  (r = %.2f)", o, cor(o_lfc, f_lfc))]
  m
}))

p <- ggplot(pd[order(hit)], aes(o_lfc, f_lfc, colour = hit)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = 0, colour = "grey88") + geom_vline(xintercept = 0, colour = "grey88") +
  geom_point(size = 0.8, alpha = 0.6) +
  facet_wrap(~ lab) +
  scale_colour_manual(values = c(`FALSE` = "grey80", `TRUE` = "#B2182B"),
                      labels = c("not significant", "age hit in the main analysis")) +
  labs(title = "Does the age effect survive when sex is held constant?",
       subtitle = paste("x: w22 vs w8 using all libraries (mixed sexes at w22).",
                        "y: the same contrast using females only."),
       x = "log2FC - all libraries", y = "log2FC - females only", colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

ggsave(file.path(FIGS, "31_age_effect_females_only.png"), p,
       width = 10, height = 5.5, dpi = 300, bg = "white")

cat("\nDone.\n")
