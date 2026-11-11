# 09 - Two questions about the PCA:
#   (a) what VST actually does to the numbers, versus raw counts and log2(x+1)
#   (b) whether restricting the PCA to the top-N variable genes changes anything
#
# Note this only concerns the PCA. The differential expression in script 03 uses
# every gene that survives the biotype and expression filters - there is no
# top-N selection anywhere in the DE path.

source("../../analyses/05_dge/scripts/00_config.R")

library(data.table)
library(DESeq2)
library(ggplot2)

samples <- fread(file.path(DATA, "sample_table.tsv"))
samples[, organ := factor(organ, levels = c("Saccule", "Utricle"))]
samples[, age   := factor(age,   levels = c("w8", "w22"))]
samples[, batch := factor(batch)]

dds <- readRDS(file.path(DATA, "dds_all_samples.rds"))   # already filtered, 18273 genes
# The cached object was saved before DESeq() ran, so the normalisation factors
# (from the tximport length offsets) still need to be computed.
dds <- estimateSizeFactors(dds)
cat("Genes in the filtered matrix:", nrow(dds), "\n")

vsd     <- vst(dds, blind = TRUE)
mat_vst <- assay(vsd)
mat_raw <- counts(dds, normalized = TRUE)
mat_log <- log2(mat_raw + 1)

# ---------------------------------------------------------------- (a) what VST does
# The point of VST is that a gene's variance should not depend on how highly it
# is expressed. Bin genes by mean expression and look at the spread in each bin.
mean_expr <- rowMeans(mat_raw)
bin <- cut(log10(mean_expr + 1), breaks = quantile(log10(mean_expr + 1), 0:8 / 8),
           include.lowest = TRUE, labels = paste0("Q", 1:8))

sd_by_bin <- rbindlist(list(
  data.table(transform = "normalised counts", bin = bin, sd = matrixStats::rowSds(mat_raw)),
  data.table(transform = "log2(count + 1)",   bin = bin, sd = matrixStats::rowSds(mat_log)),
  data.table(transform = "VST",               bin = bin, sd = matrixStats::rowSds(mat_vst))
))
summ_bin <- sd_by_bin[, .(median_sd = median(sd)), by = .(transform, bin)]
# scale each transform to its own max so three very different scales are comparable
summ_bin[, rel := median_sd / max(median_sd), by = transform]

cat("\n=== Median SD by expression octile (scaled to each transform's own max) ===\n")
cat("A flat line means the transform did its job.\n")
print(dcast(summ_bin, bin ~ transform, value.var = "rel")[
  , lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)])

summ_bin[, transform := factor(transform,
         levels = c("normalised counts", "log2(count + 1)", "VST"))]
p_var <- ggplot(summ_bin, aes(bin, rel, colour = transform, group = transform)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
  scale_colour_manual(values = c("normalised counts" = "#B2182B",
                                 "log2(count + 1)" = "#F4A582", "VST" = "#2166AC")) +
  labs(title = "What the variance-stabilising transform is for",
       subtitle = paste("Median per-gene SD across expression octiles, each curve scaled to its own maximum.",
                        "\nFlat is the goal: a gene's spread should not depend on how highly it is expressed."),
       x = "Mean expression octile (low to high)", y = "Relative median SD", colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

ggsave(file.path(FIGS, "25_vst_variance_stabilisation.png"), p_var,
       width = 8.5, height = 5.5, dpi = 300, bg = "white")

# A few concrete values, so the scale of VST is not mysterious.
cat("\n=== The same three genes on each scale (first 4 libraries) ===\n")
ex <- c(which.min(abs(mean_expr - quantile(mean_expr, .10))),
        which.min(abs(mean_expr - quantile(mean_expr, .50))),
        which.min(abs(mean_expr - quantile(mean_expr, .99))))
for (i in ex) {
  cat("\n", rownames(mat_raw)[i], " (mean normalised count ",
      round(mean_expr[i], 1), ")\n", sep = "")
  cat("  normalised counts:", paste(round(mat_raw[i, 1:4], 1), collapse = "  "), "\n")
  cat("  log2(count + 1)  :", paste(round(mat_log[i, 1:4], 2), collapse = "  "), "\n")
  cat("  VST              :", paste(round(mat_vst[i, 1:4], 2), collapse = "  "), "\n")
}

# ---------------------------------------------------------------- (b) does ntop matter
cat("\n\n=== Does the top-N choice change the PCA? ===\n")
r2 <- function(x, f) { d <- !is.na(f); summary(lm(x[d] ~ factor(f[d])))$r.squared }
sep <- function(x, f) {
  f <- droplevels(factor(f)); lv <- levels(f)
  m1 <- mean(x[f == lv[1]]); m2 <- mean(x[f == lv[2]])
  n1 <- sum(f == lv[1]); n2 <- sum(f == lv[2])
  s <- ((n1 - 1) * var(x[f == lv[1]]) + (n2 - 1) * var(x[f == lv[2]])) / (n1 + n2 - 2)
  abs(m1 - m2) / sqrt(s)
}

ntops <- c(250, 500, 1000, 2000, 5000, 10000, nrow(mat_vst))
sens <- rbindlist(lapply(ntops, function(n) {
  v   <- matrixStats::rowVars(mat_vst)
  top <- order(v, decreasing = TRUE)[seq_len(min(n, length(v)))]
  p   <- prcomp(t(mat_vst[top, ]), center = TRUE, scale. = FALSE)
  pct <- 100 * p$sdev^2 / sum(p$sdev^2)
  age_pc <- which.max(sapply(1:6, function(k) r2(p$x[, k], samples$age)))
  data.table(ntop = n,
             PC1_pct = round(pct[1], 1),
             PC1_r2_organ = round(r2(p$x[, 1], samples$organ), 3),
             organ_sep_PC1 = round(sep(p$x[, 1], samples$organ), 2),
             best_age_PC = paste0("PC", age_pc),
             age_r2 = round(r2(p$x[, age_pc], samples$age), 3))
}))
sens[ntop == nrow(mat_vst), ntop := NA]        # mark the "all genes" row
print(sens)
cat("\n(last row, ntop = NA, is all", nrow(mat_vst), "filtered genes)\n")
fwrite(sens, file.path(TABLES, "pca_ntop_sensitivity.tsv"), sep = "\t", na = "all")

p_sens <- ggplot(melt(sens[!is.na(ntop), .(ntop, `PC1 R2 organ` = PC1_r2_organ,
                                           `PC1 % variance` = PC1_pct / 100)],
                      id.vars = "ntop"),
                 aes(ntop, value, colour = variable)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
  scale_x_log10(breaks = c(250, 500, 1000, 2000, 5000, 10000)) +
  scale_colour_manual(values = c("PC1 R2 organ" = "#1B7837", "PC1 % variance" = "#762A83")) +
  labs(title = "The top-N choice barely moves the PCA",
       subtitle = "PC1 keeps tracking organ whether 250 genes or all 18,273 are used",
       x = "Number of most-variable genes used (log scale)", y = NULL, colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

ggsave(file.path(FIGS, "26_pca_ntop_sensitivity.png"), p_sens,
       width = 8.5, height = 5.5, dpi = 300, bg = "white")

cat("\nDone.\n")
