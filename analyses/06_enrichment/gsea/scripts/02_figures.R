# 02 - Figures for the GSEA results.
#
# Square panels wherever the content allows (bar charts, dot plots, the
# contrast-similarity matrix, and the 2x2 enrichment grids). The NES heatmap
# is the one genuinely rectangular figure: it has 5 columns and ~40 rows.

BASE <- "../../../analyses/05_dge"
ROOT <- file.path(dirname(BASE), "06_enrichment")   # logs are shared
PA   <- file.path(ROOT, "gsea")                # this analysis's own outputs
Sys.setenv(TMPDIR = file.path(ROOT, "logs"))

library(data.table)
library(fgsea)
library(ggplot2)
library(patchwork)

set.seed(42)
FDR_CUT <- 0.05
FIGS <- file.path(PA, "figures")

obj  <- readRDS(file.path(PA, "tables", "gsea_objects.rds"))
comb <- obj$combined
SETS <- obj$sets
de   <- readRDS(file.path(BASE, "data", "de_results.rds"))

CONTRAST_KEY <- c(A_organ_w8 = "A_utricle_vs_saccule_w8",
                  B_organ_w22 = "B_utricle_vs_saccule_w22",
                  C_age_sac = "C_w22_vs_w8_saccule",
                  D_age_utr = "D_w22_vs_w8_utricle",
                  E_interaction = "E_organ_by_age_interaction")
CONTRAST_ORDER <- names(CONTRAST_KEY)
COLL_NAME <- c(hallmark = "Hallmark", go_bp = "GO Biological Process", reactome = "Reactome")

comb[, contrast := factor(contrast, levels = CONTRAST_ORDER)]
lab_of <- function(cid) comb[contrast == cid, contrast_label[1]]

# Keep long GO names readable on an axis.
trim <- function(x, n = 52) ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "…"), x)

theme_sq <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey30", size = 10))

NES_COLS <- c(up = "#B2182B", down = "#2166AC")

# ================================================================ 1. NES bars
# Top 10 up and top 10 down by adjusted p, per contrast per collection.
for (cid in CONTRAST_ORDER) {
  for (cn in names(COLL_NAME)) {
    d <- comb[contrast == cid & collection == cn & padj < FDR_CUT]
    if (!nrow(d)) next
    top <- rbind(head(d[NES > 0][order(padj)], 10), head(d[NES < 0][order(padj)], 10))
    if (!nrow(top)) next
    top[, lab := trim(pathway_label)]
    top[, lab := factor(lab, levels = top[order(NES), lab])]
    top[, dir := fifelse(NES > 0, "up", "down")]

    p <- ggplot(top, aes(NES, lab, fill = dir)) +
      geom_col(width = 0.72) +
      geom_vline(xintercept = 0, colour = "grey30") +
      scale_fill_manual(values = NES_COLS, guide = "none") +
      labs(title = lab_of(cid),
           subtitle = paste0(COLL_NAME[cn], " - top 10 each way of ",
                             nrow(d), " significant at FDR < ", FDR_CUT,
                             "\npositive NES = ", d[NES > 0, direction[1]],
                             "   |   negative = ", d[NES < 0, direction[1]]),
           x = "Normalised enrichment score (NES)", y = NULL) +
      theme_sq + theme(axis.text.y = element_text(size = 9))

    ggsave(file.path(FIGS, cn, paste0(cid, "_nes_bar.png")), p,
           width = 9, height = 9, dpi = 300, bg = "white")
  }
}
cat("NES bar charts done\n")

# ================================================================ 2. dot plots
# Same information with set size and significance, Hallmark only.
for (cid in CONTRAST_ORDER) {
  d <- comb[contrast == cid & collection == "hallmark" & padj < FDR_CUT]
  if (!nrow(d)) next
  d[, lab := factor(trim(pathway_label), levels = trim(pathway_label)[order(d$NES)])]
  p <- ggplot(d, aes(NES, lab, size = size, colour = -log10(padj))) +
    geom_vline(xintercept = 0, colour = "grey70", linetype = "dashed") +
    geom_point() +
    scale_colour_gradient(low = "#F4A582", high = "#67001F", name = "-log10 FDR") +
    scale_size_continuous(range = c(2, 7), name = "Genes in set") +
    labs(title = lab_of(cid),
         subtitle = paste0("Hallmark - all ", nrow(d), " sets at FDR < ", FDR_CUT),
         x = "Normalised enrichment score (NES)", y = NULL) +
    theme_sq + theme(axis.text.y = element_text(size = 9))
  ggsave(file.path(FIGS, "hallmark", paste0(cid, "_dotplot.png")), p,
         width = 9, height = 9, dpi = 300, bg = "white")
}
cat("dot plots done\n")

# ================================================================ 3. enrichment
# The classic running-score curve for the four strongest Hallmark sets,
# arranged 2x2 so the panel is square.
for (cid in CONTRAST_ORDER) {
  res <- as.data.table(de[[CONTRAST_KEY[cid]]])
  r <- res[!is.na(stat), .(gene_id, stat)]
  ranks <- sort(setNames(r$stat, r$gene_id), decreasing = TRUE)

  top4 <- head(comb[contrast == cid & collection == "hallmark" & padj < FDR_CUT][order(padj)], 4)
  if (!nrow(top4)) next
  plots <- lapply(seq_len(nrow(top4)), function(i) {
    pw <- top4$pathway[i]
    plotEnrichment(SETS$hallmark[[pw]], ranks) +
      labs(title = trim(top4$pathway_label[i], 40),
           subtitle = sprintf("NES %+.2f   FDR %.1e", top4$NES[i], top4$padj[i]),
           x = "Rank in the ordered gene list", y = "Enrichment score") +
      theme_bw(base_size = 10) +
      theme(plot.title = element_text(face = "bold", size = 10),
            plot.subtitle = element_text(size = 8.5, colour = "grey30"),
            panel.grid.minor = element_blank())
  })
  dd <- comb[contrast == cid & collection == "hallmark"]
  up_dir <- dd[NES > 0, direction[1]]
  dn_dir <- dd[NES < 0, direction[1]]
  p <- wrap_plots(plots, ncol = 2) +
    plot_annotation(title = paste0(lab_of(cid), " - four strongest Hallmark sets"),
                    subtitle = paste0("Genes ranked by the DESeq2 Wald statistic: ",
                                      "top of the list = ", up_dir,
                                      ", bottom = ", dn_dir, ".\n",
                                      "A curve that climbs above zero means the set is ",
                                      "concentrated at the top; one that dips below zero, ",
                                      "at the bottom."),
                    theme = theme(plot.title = element_text(face = "bold", size = 13),
                                  plot.subtitle = element_text(size = 10, colour = "grey30")))
  ggsave(file.path(FIGS, "hallmark", paste0(cid, "_enrichment_curves.png")), p,
         width = 9, height = 9, dpi = 300, bg = "white")
}
cat("enrichment curves done\n")

# ================================================================ 4. summary
# 4a. how many pathways each contrast turns up
cnt <- comb[padj < FDR_CUT, .(n = .N), by = .(contrast, collection,
                                              dir = fifelse(NES > 0, "up", "down"))]
cnt[, collection := factor(COLL_NAME[collection], levels = COLL_NAME)]
p <- ggplot(cnt, aes(contrast, n, fill = dir)) +
  geom_col(position = "stack", width = 0.7) +
  facet_wrap(~ collection, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = NES_COLS, name = NULL,
                    labels = c(up = "positive NES", down = "negative NES")) +
  labs(title = "Significantly enriched gene sets per contrast",
       subtitle = paste0("FDR < ", FDR_CUT, ". Note the y scale differs by collection."),
       x = NULL, y = "Gene sets") +
  theme_sq + theme(axis.text.x = element_text(angle = 20, hjust = 1),
                   legend.position = "top")
ggsave(file.path(FIGS, "summary", "01_significant_per_contrast.png"), p,
       width = 9, height = 9, dpi = 300, bg = "white")

# 4b. Hallmark NES across all five contrasts - the one rectangular figure
h <- comb[collection == "hallmark"]
keep <- h[padj < FDR_CUT, unique(pathway)]
h <- h[pathway %in% keep]
ord <- h[, .(m = mean(NES)), by = pathway_label][order(m), pathway_label]
h[, lab := factor(trim(pathway_label), levels = trim(ord))]
h[, star := fifelse(padj < 0.001, "***", fifelse(padj < 0.01, "**",
             fifelse(padj < FDR_CUT, "*", "")))]
p <- ggplot(h, aes(contrast, lab, fill = NES)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = star), size = 3, vjust = 0.75) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(title = "Hallmark enrichment across all five contrasts",
       subtitle = paste0("NES; * FDR<0.05, ** <0.01, *** <0.001. ",
                         length(keep), " sets significant in at least one contrast."),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 1),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(face = "bold", size = 13))
ggsave(file.path(FIGS, "summary", "02_hallmark_nes_heatmap.png"), p,
       width = 9, height = 11, dpi = 300, bg = "white")

# 4c. how similar are the contrasts at the pathway level - naturally square
w <- dcast(comb[collection == "go_bp"], pathway ~ contrast, value.var = "NES")
mat <- cor(as.matrix(w[, -1]), use = "pairwise.complete.obs", method = "spearman")
cm <- as.data.table(as.table(mat)); setnames(cm, c("x", "y", "rho"))
cm[, x := factor(x, levels = CONTRAST_ORDER)][, y := factor(y, levels = rev(CONTRAST_ORDER))]
p <- ggplot(cm, aes(x, y, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 4) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Do the contrasts agree at the pathway level?",
       subtitle = "Spearman correlation of GO:BP NES between contrasts",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 1),
        plot.title = element_text(face = "bold", size = 13))
ggsave(file.path(FIGS, "summary", "03_contrast_similarity.png"), p,
       width = 8, height = 8, dpi = 300, bg = "white")
fwrite(cm, file.path(PA, "tables", "contrast_nes_correlation.tsv"), sep = "\t")

cat("summary figures done\n")
cat("\nFigures written under", FIGS, "\n")
print(list.files(FIGS, recursive = TRUE))
