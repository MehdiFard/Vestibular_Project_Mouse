# 07 - Does keeping all biotypes really separate the groups better?
#
# Compares the PCA built from all biotypes against the protein-coding + lincRNA
# one, using measures that do not depend on how much total variance is in the
# matrix. Also asks WHICH biotypes end up in the top-2000 variable genes that
# the PCA actually uses, and what those genes track.

source("../../analyses/05_dge/scripts/00_config.R")

library(data.table)
library(DESeq2)
library(ggplot2)

samples <- fread(file.path(DATA, "sample_table.tsv"))
samples[, organ := factor(organ, levels = c("Saccule", "Utricle"))]
samples[, age   := factor(age,   levels = c("w8", "w22"))]
samples[, batch := factor(batch)]
samples[, group := factor(paste(organ, age), levels = names(GROUP_COLS))]

txi  <- readRDS(file.path(DATA, "txi_gene.rds"))
anno <- fread(file.path(DATA, "gene_annotation.tsv"))
stopifnot(identical(colnames(txi$counts), samples$library_id))

coldata <- data.frame(samples, row.names = samples$library_id)
dds_all <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ organ)
min_n   <- min(table(samples$group))

# Two gene sets, each with the same expression filter applied.
make_vst <- function(dds, ids = NULL) {
  if (!is.null(ids)) dds <- dds[rownames(dds) %in% ids, ]
  dds <- dds[rowSums(counts(dds) >= 10) >= min_n, ]
  list(mat = assay(vst(dds, blind = TRUE)), n = nrow(dds))
}

set_all  <- make_vst(dds_all)
set_pcl  <- make_vst(dds_all, anno[biotype %in% KEEP_BIOTYPES, gene_id])

cat("Genes after filtering - all biotypes:", set_all$n,
    " protein_coding + lincRNA:", set_pcl$n, "\n\n")

# ---------------------------------------------------------------- separation
# Standardised centroid distance: how far apart two groups are in units of the
# pooled within-group spread. Unlike "% variance explained by PC1" this does
# not change just because the matrix has more noise genes in it.
sep <- function(x, f) {
  f <- droplevels(factor(f)); lv <- levels(f)
  m1 <- mean(x[f == lv[1]]); m2 <- mean(x[f == lv[2]])
  s1 <- var(x[f == lv[1]]);  s2 <- var(x[f == lv[2]])
  n1 <- sum(f == lv[1]);     n2 <- sum(f == lv[2])
  abs(m1 - m2) / sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
}

pca_of <- function(mat, ntop = 2000) {
  v   <- matrixStats::rowVars(mat)
  top <- order(v, decreasing = TRUE)[seq_len(min(ntop, length(v)))]
  p   <- prcomp(t(mat[top, ]), center = TRUE, scale. = FALSE)
  list(x = p$x, pct = 100 * p$sdev^2 / sum(p$sdev^2), genes = rownames(mat)[top])
}

p_all <- pca_of(set_all$mat)
p_pcl <- pca_of(set_pcl$mat)

r2 <- function(x, f) { d <- !is.na(f); summary(lm(x[d] ~ factor(f[d])))$r.squared }

res <- data.table(
  gene_set = c("all biotypes", "protein_coding + lincRNA"),
  n_genes  = c(set_all$n, set_pcl$n),
  PC1_pct  = c(p_all$pct[1], p_pcl$pct[1]),
  PC2_pct  = c(p_all$pct[2], p_pcl$pct[2]),
  # separation measured in the standardised way
  organ_sep_PC1 = c(sep(p_all$x[, 1], samples$organ), sep(p_pcl$x[, 1], samples$organ)),
  age_sep_PC2   = c(sep(p_all$x[, 2], samples$age),   sep(p_pcl$x[, 2], samples$age)),
  # what each PC tracks
  PC1_r2_organ  = c(r2(p_all$x[, 1], samples$organ), r2(p_pcl$x[, 1], samples$organ)),
  PC2_r2_age    = c(r2(p_all$x[, 2], samples$age),   r2(p_pcl$x[, 2], samples$age)),
  PC2_r2_batch  = c(r2(p_all$x[, 2], samples$batch), r2(p_pcl$x[, 2], samples$batch)))

cat("=== PCA comparison ===\n")
print(res[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)])

# ---------------------------------------------------------------- what's in the top 2000
cat("\n=== Biotype makeup of the top-2000 variable genes (all-biotype run) ===\n")
bt <- anno[match(p_all$genes, gene_id), biotype]
tab <- sort(table(bt), decreasing = TRUE)
print(head(tab, 10))
cat("\nNon (protein_coding/lincRNA) share of the top 2000:",
    sprintf("%.1f%%\n", 100 * mean(!bt %in% KEEP_BIOTYPES)))
cat("Their share of the filtered matrix overall:",
    sprintf("%.1f%%\n", 100 * mean(!anno[match(rownames(set_all$mat), gene_id), biotype]
                                   %in% KEEP_BIOTYPES)))

# ---------------------------------------------------------------- what do they track?
# For every gene in the filtered all-biotype matrix, how much of its variance is
# explained by organ, by age and by batch. Then compare across biotype classes.
cat("\n=== What each biotype class tracks (mean R-squared per gene) ===\n")
mat <- set_all$mat
keep_bt <- anno[match(rownames(mat), gene_id), biotype]
cls <- ifelse(keep_bt %in% KEEP_BIOTYPES, "protein_coding + lincRNA", "other biotypes")

mm_organ <- model.matrix(~ samples$organ)
mm_age   <- model.matrix(~ samples$age)
mm_batch <- model.matrix(~ samples$batch)

# R-squared of every gene against a design, all genes at once.
# Y is genes x samples, X is samples x p; work in samples x genes internally.
fast_r2 <- function(Y, X) {
  Yt  <- t(Y)                                     # samples x genes
  fit <- X %*% qr.solve(crossprod(X), crossprod(X, Yt))
  rss <- colSums((Yt - fit)^2)
  tss <- colSums(sweep(Yt, 2, colMeans(Yt))^2)
  1 - rss / tss
}
r2_organ <- fast_r2(mat, mm_organ)
r2_age   <- fast_r2(mat, mm_age)
r2_batch <- fast_r2(mat, mm_batch)

byclass <- data.table(class = cls, organ = r2_organ, age = r2_age, batch = r2_batch)
print(byclass[, .(n = .N,
                  mean_r2_organ = round(mean(organ), 3),
                  mean_r2_age   = round(mean(age),   3),
                  mean_r2_batch = round(mean(batch), 3)), by = class])

cat("\nSame, restricted to the top-2000 variable genes of the all-biotype run:\n")
top_idx <- match(p_all$genes, rownames(mat))
print(byclass[top_idx][, .(n = .N,
                           mean_r2_organ = round(mean(organ), 3),
                           mean_r2_age   = round(mean(age),   3),
                           mean_r2_batch = round(mean(batch), 3)), by = class])

# Where did the age signal go? Age does not disappear when the noise biotypes
# are dropped - it just stops dominating an early PC. Report the best PC for
# each factor in each gene set.
cat("\n=== Which PC best tracks each factor (R-squared, first 10 PCs) ===\n")
best_pc <- rbindlist(lapply(list(`all biotypes` = p_all,
                                 `protein_coding + lincRNA` = p_pcl), function(p) {
  rbindlist(lapply(c("organ", "age", "batch"), function(f) {
    v <- sapply(1:10, function(k) r2(p$x[, k], samples[[f]]))
    data.table(factor = f, best_PC = paste0("PC", which.max(v)),
               r2 = round(max(v), 3), pct_var_of_that_PC = round(p$pct[which.max(v)], 1))
  }))
}), idcol = "gene_set")
print(best_pc)
fwrite(best_pc, file.path(TABLES, "biotype_best_pc_per_factor.tsv"), sep = "\t")

fwrite(res, file.path(TABLES, "biotype_separation_comparison.tsv"), sep = "\t")
fwrite(byclass[, .(mean_r2_organ = mean(organ), mean_r2_age = mean(age),
                   mean_r2_batch = mean(batch), n = .N), by = class],
       file.path(TABLES, "biotype_variance_tracking.tsv"), sep = "\t")

# ---------------------------------------------------------------- figure
# Three panels that answer the question directly: the organ split is identical,
# the age band on PC2 of the all-biotype run is really batch, and age itself is
# still there in the filtered data - on PC3, where it is cleaner.
panel <- function(p, xi, yi, lab) {
  data.table(library_id = rownames(p$x), x = p$x[, xi], y = p$x[, yi],
             xlab = sprintf("PC%d (%.1f%%)", xi, p$pct[xi]),
             ylab = sprintf("PC%d (%.1f%%)", yi, p$pct[yi]), lab = lab)
}

plot_dt <- rbindlist(list(
  panel(p_all, 1, 2, sprintf("all biotypes (%d genes)\nPC2: batch R2=%.2f, age R2=%.2f",
                             set_all$n, r2(p_all$x[, 2], samples$batch),
                             r2(p_all$x[, 2], samples$age))),
  panel(p_pcl, 1, 2, sprintf("protein_coding + lincRNA (%d genes)\nPC2: batch R2=%.2f, age R2=%.2f",
                             set_pcl$n, r2(p_pcl$x[, 2], samples$batch),
                             r2(p_pcl$x[, 2], samples$age))),
  panel(p_pcl, 1, 3, sprintf("protein_coding + lincRNA, PC3\nPC3: batch R2=%.2f, age R2=%.2f",
                             r2(p_pcl$x[, 3], samples$batch),
                             r2(p_pcl$x[, 3], samples$age)))))
plot_dt[, lab := factor(lab, levels = unique(lab))]
plot_dt <- merge(plot_dt, samples[, .(library_id, organ, age)], by = "library_id")

p <- ggplot(plot_dt, aes(x, y, colour = organ, shape = age)) +
  geom_point(size = 2.1, alpha = 0.85) +
  facet_wrap(~ lab, scales = "free", nrow = 1) +
  scale_colour_manual(values = ORGAN_COLS) +
  scale_shape_manual(values = c(w8 = 16, w22 = 17)) +
  labs(title = "The organ split is identical; only the age axis moves",
       subtitle = paste("Left: with all biotypes, the vertical age band on PC2 tracks batch more than age.",
                        "Right: after filtering, age reappears on PC3, stronger and cleaner."),
       x = "PC1", y = "PC2 (left, middle) / PC3 (right)",
       colour = "Organ", shape = "Age") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

ggsave(file.path(FIGS, "22_biotype_pca_comparison.png"), p,
       width = 14, height = 5.0, dpi = 300, bg = "white")

cat("\nDone.\n")
