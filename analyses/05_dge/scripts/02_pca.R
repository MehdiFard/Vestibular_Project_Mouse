# 02 - PCA on expression across all four organ-age groups.
#
# Colour = organ, shape = age, as requested.
# Values are VST (variance-stabilising transform) of the gene counts, which is
# DESeq2's recommended input for distance/PCA work - it behaves better at low
# counts than log2(TPM + 1).

source("../../analyses/05_dge/scripts/00_config.R")

library(data.table)
library(DESeq2)
library(ggplot2)

samples <- fread(file.path(DATA, "sample_table.tsv"))
samples[, organ := factor(organ, levels = c("Saccule", "Utricle"))]
samples[, age   := factor(age,   levels = c("w8", "w22"))]
samples[, batch := factor(batch)]
samples[, group := factor(paste(organ, age), levels = names(GROUP_COLS))]

txi      <- readRDS(file.path(DATA, "txi_gene.rds"))
gene_anno <- fread(file.path(DATA, "gene_annotation.tsv"))

stopifnot(identical(colnames(txi$counts), samples$library_id))

# ---------------------------------------------------------------- build + filter
coldata <- data.frame(samples, row.names = samples$library_id)
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ organ)

cat("Genes in the annotation:", nrow(dds), "\n")

# Biotype filter first: only protein-coding and lincRNA genes are analysed.
keep_bt <- rownames(dds) %in% gene_anno[biotype %in% KEEP_BIOTYPES, gene_id]
dds <- dds[keep_bt, ]
cat("After biotype filter (", paste(KEEP_BIOTYPES, collapse = ", "), "):",
    nrow(dds), "\n")

# Expression filter: keep a gene if it reaches 10 counts in at least as many
# samples as the smallest group (n = 55). This is the usual edgeR/DESeq2
# heuristic and is stricter than "non-zero in >= 20 samples".
min_n <- min(table(samples$group))
keep  <- rowSums(counts(dds) >= 10) >= min_n

cat("Genes kept (>=10 counts in >=", min_n, "samples):", sum(keep), "\n")

dds <- dds[keep, ]
saveRDS(dds, file.path(DATA, "dds_all_samples.rds"))

# ---------------------------------------------------------------- VST + PCA
vsd <- vst(dds, blind = TRUE)          # blind: don't let the design inform the transform
mat <- assay(vsd)

# PCA on the 2000 most variable genes, the usual convention
vars     <- matrixStats::rowVars(mat)
top      <- order(vars, decreasing = TRUE)[seq_len(min(2000, length(vars)))]
pca      <- prcomp(t(mat[top, ]), center = TRUE, scale. = FALSE)
pct      <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

pcs <- data.table(library_id = rownames(pca$x), as.data.table(pca$x[, 1:10]))
pcs <- merge(samples, pcs, by = "library_id", sort = FALSE)
fwrite(pcs, file.path(TABLES, "pca_scores.tsv"), sep = "\t", na = "NA")

fwrite(data.table(PC = paste0("PC", seq_along(pct)), percent_variance = pct),
       file.path(TABLES, "pca_variance_explained.tsv"), sep = "\t")

# ---------------------------------------------------------------- main PCA plot
axis_lab <- function(i) paste0("PC", i, " (", pct[i], "%)")

p_main <- ggplot(pcs, aes(PC1, PC2, colour = organ, shape = age)) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_point(size = 2.6, alpha = 0.85, stroke = 0.9) +
  scale_colour_manual(values = ORGAN_COLS) +
  scale_shape_manual(values = c(w8 = 16, w22 = 17)) +
  labs(title    = "Expression PCA - all 260 libraries",
       subtitle = "Variance-stabilised expression (DESeq2 vst), top 2000 most variable genes",
       x = axis_lab(1), y = axis_lab(2),
       colour = "Organ", shape = "Age") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(FIGS, "01_pca_pc1_pc2_organ_age.png"), p_main,
       width = 8, height = 6, dpi = 300, bg = "white")

# Same plot with 95% ellipses per organ-age group, easier to read the overlap
p_ell <- p_main +
  stat_ellipse(aes(group = group, colour = organ, linetype = age),
               type = "norm", level = 0.95, linewidth = 0.5, show.legend = FALSE) +
  labs(title = "Expression PCA with 95% group ellipses")

ggsave(file.path(FIGS, "02_pca_pc1_pc2_ellipses.png"), p_ell,
       width = 8, height = 6, dpi = 300, bg = "white")

# ---------------------------------------------------------------- same PCA, by sex
# Colour = sex, shape = organ. Sex is known for all 260 libraries: it was
# recovered from Xist and the Y genes where the metadata sheet had none
# (scripts/13_infer_sex.R). Two panels, because sex means different things in
# the two halves of the design: every 8-week animal is female, so only the
# 22-week libraries can say anything about sex on their own.
w22_lab <- sprintf("22 weeks only (%d female, %d male)",
                   pcs[age == "w22" & sex_final == "F", .N],
                   pcs[age == "w22" & sex_final == "M", .N])
all_lab <- "All 260 libraries (every 8-week animal is female)"

pcs_sex <- rbind(copy(pcs)[, panel := all_lab],
                 copy(pcs)[age == "w22"][, panel := w22_lab])
pcs_sex[, panel := factor(panel, levels = c(all_lab, w22_lab))]
pcs_sex[, sex_final := factor(sex_final, levels = c("F", "M"))]

# R-squared of a PC against a covariate, used in the subtitles below.
r2_of <- function(d, pc, cv) summary(lm(d[[pc]] ~ factor(d[[cv]])))$r.squared
w22   <- pcs[age == "w22"]

sex_panels <- function(dat, mapping, title, subtitle, xlab, ylab) {
  ggplot(dat, mapping) +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
    geom_point(size = 2.4, alpha = 0.85, stroke = 0.9) +
    scale_colour_manual(values = SEX_COLS) +
    scale_shape_manual(values = c(Saccule = 16, Utricle = 17)) +
    facet_wrap(~ panel) +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab,
         colour = "Sex", shape = "Organ") +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA))
}

p_sex12 <- sex_panels(
  pcs_sex[panel == all_lab],
  aes(PC1, PC2, colour = sex_final, shape = organ),
  "Expression PCA coloured by sex",
  sprintf(paste0("The organ axis is free of sex: PC1 is explained by organ (R2 = %.2f)\n",
                 "and not by sex (R2 = %.3f across all libraries, %.3f within 22 weeks).\n",
                 "The two sexes sit on top of each other within each organ."),
          r2_of(pcs, "PC1", "organ"), r2_of(pcs, "PC1", "sex_final"),
          r2_of(w22, "PC1", "sex_final")),
  axis_lab(1), axis_lab(2))

ggsave(file.path(FIGS, "08_pca_by_sex.png"), p_sex12,
       width = 8, height = 6, dpi = 300, bg = "white")

# PC3 is the component sex actually loads on, so plot it directly.
p_sex13 <- sex_panels(
  pcs_sex,
  aes(PC1, PC3, colour = sex_final, shape = organ),
  "Where sex does show up: PC3",
  sprintf(paste0("Across all libraries PC3 cannot separate sex (R2 = %.2f) from age (R2 = %.2f) - ",
                 "they are the same split.\n",
                 "Within 22 weeks sex explains R2 = %.2f of PC3, but batch explains %.2f of the ",
                 "same component:\n",
                 "batch 4 is entirely female and batch 3 is %d male to %d female, so sex and batch stay tangled."),
          r2_of(pcs, "PC3", "sex_final"), r2_of(pcs, "PC3", "age"),
          r2_of(w22, "PC3", "sex_final"), r2_of(w22, "PC3", "batch"),
          w22[batch == 3 & sex_final == "M", .N], w22[batch == 3 & sex_final == "F", .N]),
  axis_lab(1), axis_lab(3))

ggsave(file.path(FIGS, "09_pca_pc1_pc3_by_sex.png"), p_sex13,
       width = 11, height = 5.2, dpi = 300, bg = "white")

# The sex/batch tangle inside the 22-week libraries, in one small table.
sex_batch <- dcast(w22[, .N, by = .(batch, sex_final)], batch ~ sex_final,
                   value.var = "N", fill = 0)
fwrite(sex_batch, file.path(TABLES, "w22_sex_by_batch.tsv"), sep = "\t")
cat("\n=== sex by batch within the 22-week libraries ===\n"); print(sex_batch)

# ---------------------------------------------------------------- scree plot
scree <- data.table(PC = 1:10, var = pct[1:10])
p_scree <- ggplot(scree, aes(PC, var)) +
  geom_col(fill = "#4C72B0") +
  geom_text(aes(label = paste0(var, "%")), vjust = -0.4, size = 3.2) +
  scale_x_continuous(breaks = 1:10) +
  labs(title = "Variance explained by the first 10 PCs",
       x = "Principal component", y = "% variance") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(FIGS, "03_pca_scree.png"), p_scree,
       width = 7, height = 5, dpi = 300, bg = "white")

# ---------------------------------------------------------------- batch check
# Same coordinates, coloured by sequencing batch, to show what PC1/PC2 are
# really tracking. Age and batch are partly confounded in this design.
p_batch <- ggplot(pcs, aes(PC1, PC2, colour = batch, shape = age)) +
  geom_point(size = 2.4, alpha = 0.85) +
  scale_shape_manual(values = c(w8 = 16, w22 = 17)) +
  labs(title = "Same PCA coloured by sequencing batch",
       subtitle = "Batch 1 is w8-only and batch 4 is w22-only, so batch and age overlap",
       x = axis_lab(1), y = axis_lab(2), colour = "Batch", shape = "Age") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(FIGS, "04_pca_by_batch.png"), p_batch,
       width = 8, height = 6, dpi = 300, bg = "white")

# ---------------------------------------------------------------- batch-corrected
# Remove the batch effect while protecting organ and age, then redo the PCA.
mat_adj <- limma::removeBatchEffect(
  mat, batch = samples$batch,
  design = model.matrix(~ organ * age, data = samples))

vars_adj <- matrixStats::rowVars(mat_adj)
top_adj  <- order(vars_adj, decreasing = TRUE)[seq_len(min(2000, length(vars_adj)))]
pca_adj  <- prcomp(t(mat_adj[top_adj, ]), center = TRUE, scale. = FALSE)
pct_adj  <- round(100 * pca_adj$sdev^2 / sum(pca_adj$sdev^2), 1)

pcs_adj <- data.table(library_id = rownames(pca_adj$x), as.data.table(pca_adj$x[, 1:5]))
pcs_adj <- merge(samples, pcs_adj, by = "library_id", sort = FALSE)
fwrite(pcs_adj, file.path(TABLES, "pca_scores_batch_corrected.tsv"), sep = "\t", na = "NA")

p_adj <- ggplot(pcs_adj, aes(PC1, PC2, colour = organ, shape = age)) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_point(size = 2.6, alpha = 0.85) +
  scale_colour_manual(values = ORGAN_COLS) +
  scale_shape_manual(values = c(w8 = 16, w22 = 17)) +
  labs(title    = "Expression PCA after removing the batch effect",
       subtitle = "limma::removeBatchEffect, organ and age protected",
       x = paste0("PC1 (", pct_adj[1], "%)"), y = paste0("PC2 (", pct_adj[2], "%)"),
       colour = "Organ", shape = "Age") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(FIGS, "05_pca_batch_corrected.png"), p_adj,
       width = 8, height = 6, dpi = 300, bg = "white")

# ---------------------------------------------------------------- PC ~ covariate
# How much of each of the first 6 PCs is explained by each known covariate.
# sex_final is the recovered sex, complete for all 260; the raw "sex" column
# from the metadata sheet is missing for 52 libraries and would bias this.
covars <- c("organ", "age", "batch", "sex_final", "strain")
r2 <- rbindlist(lapply(paste0("PC", 1:6), function(pc) {
  rbindlist(lapply(covars, function(cv) {
    d <- pcs[!is.na(get(cv)), ]
    if (uniqueN(d[[cv]]) < 2) return(NULL)
    fit <- lm(d[[pc]] ~ factor(d[[cv]]))
    data.table(PC = pc, covariate = cv, r2 = summary(fit)$r.squared)
  }))
}))
r2[covariate == "sex_final", covariate := "sex"]
r2[, PC := factor(PC, levels = paste0("PC", 1:6))]
fwrite(r2, file.path(TABLES, "pca_covariate_r2.tsv"), sep = "\t")

p_r2 <- ggplot(r2, aes(PC, covariate, fill = r2)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r2)), size = 3.4,
            colour = ifelse(r2$r2 > 0.5, "white", "black")) +
  scale_fill_gradient(low = "#F7F7F7", high = "#B2182B", limits = c(0, 1)) +
  labs(title = "Variance in each PC explained by each covariate",
       subtitle = "R-squared from PC ~ covariate",
       x = NULL, y = NULL, fill = expression(R^2)) +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank())

ggsave(file.path(FIGS, "06_pca_covariate_r2.png"), p_r2,
       width = 8, height = 4.5, dpi = 300, bg = "white")

# ---------------------------------------------------------------- sample distances
# Correlation heatmap of the group means, a compact view of how the four
# organ-age groups relate to each other.
grp_mean <- sapply(levels(samples$group), function(g)
  rowMeans(mat[, samples$group == g, drop = FALSE]))
cor_mat <- cor(grp_mean, method = "spearman")

cor_dt <- as.data.table(as.table(cor_mat))
setnames(cor_dt, c("g1", "g2", "rho"))
p_cor <- ggplot(cor_dt, aes(g1, g2, fill = rho)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.3f", rho)), size = 3.4) +
  scale_fill_gradient(low = "#FDDBC7", high = "#2166AC") +
  labs(title = "Spearman correlation between organ-age group means",
       x = NULL, y = NULL, fill = "rho") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(FIGS, "07_group_correlation.png"), p_cor,
       width = 7, height = 5.5, dpi = 300, bg = "white")
fwrite(data.table(group = rownames(cor_mat), cor_mat),
       file.path(TABLES, "group_mean_correlation.tsv"), sep = "\t")

cat("\nPC1:", pct[1], "%  PC2:", pct[2], "%  PC3:", pct[3], "%\n")
cat("Top covariate associations:\n")
print(r2[order(-r2)][1:8])
cat("\nPCA done.\n")
