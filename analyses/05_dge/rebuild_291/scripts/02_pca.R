# 02 - PCA on the 291-library set.
source("../../analyses/05_dge/rebuild_291/scripts/00_config.R")
library(data.table); library(DESeq2); library(ggplot2)

samples <- fread(file.path(DATA, "sample_table.tsv"))
samples[, organ := factor(organ, levels = c("Saccule", "Utricle"))]
samples[, age := factor(age, levels = c("w8", "w22"))]
samples[, batch := factor(batch)]
samples[, group := factor(paste(organ, age), levels = names(GROUP_COLS))]
txi <- readRDS(file.path(DATA, "txi_gene.rds"))
anno <- fread(file.path(DATA, "gene_annotation.tsv"))
stopifnot(identical(colnames(txi$counts), samples$library_id))

dds <- DESeqDataSetFromTximport(txi, colData = data.frame(samples, row.names = samples$library_id),
                                design = ~ organ)
dds <- dds[rownames(dds) %in% anno[biotype %in% KEEP_BIOTYPES, gene_id], ]
min_n <- min(table(samples$group))
keep  <- rowSums(counts(dds) >= 10) >= min_n
cat("smallest group:", min_n, "| genes kept:", sum(keep), "\n")
dds <- dds[keep, ]
saveRDS(dds, file.path(DATA, "dds_all_samples.rds"))

vsd <- vst(dds, blind = TRUE); mat <- assay(vsd)
top <- order(matrixStats::rowVars(mat), decreasing = TRUE)[seq_len(min(2000, nrow(mat)))]
pca <- prcomp(t(mat[top, ]), center = TRUE, scale. = FALSE)
pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pcs <- merge(samples, data.table(library_id = rownames(pca$x), as.data.table(pca$x[, 1:10])),
             by = "library_id", sort = FALSE)
fwrite(pcs, file.path(TABLES, "pca_scores.tsv"), sep = "\t", na = "NA")
fwrite(data.table(PC = paste0("PC", seq_along(pct))[1:20], percent_variance = pct[1:20]),
       file.path(TABLES, "pca_variance_explained.tsv"), sep = "\t")

axis_lab <- function(i) paste0("PC", i, " (", pct[i], "%)")
p <- ggplot(pcs, aes(PC1, PC2, colour = organ, shape = age)) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_point(size = 2.6, alpha = 0.85) +
  scale_colour_manual(values = ORGAN_COLS) +
  scale_shape_manual(values = c(w8 = 16, w22 = 17)) +
  labs(title = paste0("Expression PCA - all ", nrow(pcs), " libraries"),
       subtitle = "Variance-stabilised expression, top 2000 most variable genes",
       x = axis_lab(1), y = axis_lab(2), colour = "Organ", shape = "Age") +
  theme_bw(base_size = 13) + theme(panel.grid.minor = element_blank())
ggsave(file.path(FIGS, "01_pca_pc1_pc2.png"), p, width = 8, height = 6, dpi = 300, bg = "white")

covars <- c("organ", "age", "batch", "sex_final", "strain")
r2 <- rbindlist(lapply(paste0("PC", 1:6), function(pc) rbindlist(lapply(covars, function(cv) {
  d <- pcs[!is.na(get(cv))]
  if (uniqueN(d[[cv]]) < 2) return(NULL)
  data.table(PC = pc, covariate = cv, r2 = summary(lm(d[[pc]] ~ factor(d[[cv]])))$r.squared)
}))))
r2[covariate == "sex_final", covariate := "sex"]
r2[, PC := factor(PC, levels = paste0("PC", 1:6))]
fwrite(r2, file.path(TABLES, "pca_covariate_r2.tsv"), sep = "\t")

p2 <- ggplot(r2, aes(PC, covariate, fill = r2)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r2)), size = 3.4,
            colour = ifelse(r2$r2 > 0.5, "white", "black")) +
  scale_fill_gradient(low = "#F7F7F7", high = "#B2182B", limits = c(0, 1)) +
  labs(title = "Variance in each PC explained by each covariate",
       subtitle = paste0("R-squared, ", nrow(pcs), " libraries"), x = NULL, y = NULL,
       fill = expression(R^2)) +
  theme_minimal(base_size = 13) + theme(panel.grid = element_blank())
ggsave(file.path(FIGS, "02_pca_covariate_r2.png"), p2, width = 8, height = 4.5, dpi = 300, bg = "white")

cat("\nPC1:", pct[1], "% PC2:", pct[2], "% PC3:", pct[3], "%\n")
print(r2[order(-r2)][1:8])
cat("\nPCA done.\n")
