# 04 - Figures for the rebuild, and the direct comparison against the 260-library run.
source("../../analyses/05_dge/rebuild_291/scripts/00_config.R")
library(data.table); library(ggplot2)

v2 <- readRDS(file.path(DATA, "de_results.rds"))
v1 <- readRDS(file.path(V1, "data", "de_results.rds"))
s2 <- fread(file.path(TABLES, "DE_summary.tsv"))
s1 <- fread(file.path(V1, "results", "tables", "DE_summary.tsv"))

KEYS <- c(A_utricle_vs_saccule_w8 = "Utricle vs Saccule, w8",
          B_utricle_vs_saccule_w22 = "Utricle vs Saccule, w22",
          C_w22_vs_w8_saccule = "w22 vs w8, Saccule",
          D_w22_vs_w8_utricle = "w22 vs w8, Utricle",
          E_organ_by_age_interaction = "Organ x age interaction")

# ---------------------------------------------------------------- volcanoes
for (k in names(KEYS)) {
  d <- copy(v2[[k]])
  d[, sig := !is.na(padj) & padj < FDR_CUT & abs(log2FoldChange) > LFC_CUT]
  d[, y := -log10(pmax(pvalue, 1e-300))]
  cap <- quantile(d$y, 0.9995, na.rm = TRUE); d[, y := pmin(y, cap)]
  lab <- rbind(head(d[sig == TRUE][order(padj)], 8),
               head(d[sig == TRUE][order(-abs(log2FoldChange))], 8))
  lab <- unique(lab, by = "gene_id")
  p <- ggplot(d, aes(log2FoldChange, y)) +
    geom_point(aes(colour = sig), size = 0.7, alpha = 0.5) +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", colour = "grey50") +
    scale_colour_manual(values = c(`FALSE` = "grey78", `TRUE` = "#B2182B"), guide = "none") +
    geom_text(data = lab, aes(label = gene_name), size = 2.7, vjust = -0.6, check_overlap = TRUE) +
    labs(title = KEYS[[k]],
         subtitle = paste0(s2[analysis == k, sig_fdr_lfc], " genes at FDR < ", FDR_CUT,
                           " and |log2FC| > ", LFC_CUT, "   (", s2[analysis == k, n_libraries], " libraries)"),
         x = "log2 fold change", y = "-log10 p") +
    theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())
  ggsave(file.path(FIGS, paste0("10_volcano_", k, ".png")), p, width = 7, height = 7, dpi = 300, bg = "white")
}

# ---------------------------------------------------------------- v1 vs v2 counts
cmp <- merge(s1[analysis %in% names(KEYS), .(analysis, n1 = n_libraries, fdr1 = sig_fdr, sig1 = sig_fdr_lfc)],
             s2[, .(analysis, n2 = n_libraries, fdr2 = sig_fdr, sig2 = sig_fdr_lfc)], by = "analysis")
cmp[, label := KEYS[analysis]]
cmp[, `:=`(d_lib = n2 - n1, d_fdr = fdr2 - fdr1, d_sig = sig2 - sig1)]
fwrite(cmp, file.path(TABLES, "v1_vs_v2_summary.tsv"), sep = "\t")
print(cmp[, .(label, n1, n2, fdr1, fdr2, sig1, sig2)])

m <- melt(cmp[, .(label, `260 libraries` = sig1, `291 libraries` = sig2)],
          id.vars = "label", variable.name = "run", value.name = "n")
p <- ggplot(m, aes(label, n, fill = run)) +
  geom_col(position = position_dodge(0.8), width = 0.72) +
  geom_text(aes(label = n), position = position_dodge(0.8), vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = c(`260 libraries` = "#9E9E9E", `291 libraries` = "#1B7837"), name = NULL) +
  labs(title = "Differentially expressed genes: 260 vs 291 libraries",
       subtitle = paste0("FDR < ", FDR_CUT, " and |log2FC| > ", LFC_CUT),
       x = NULL, y = "Genes") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 18, hjust = 1), legend.position = "top",
        panel.grid.minor = element_blank())
ggsave(file.path(FIGS, "20_v1_vs_v2_counts.png"), p, width = 9, height = 7, dpi = 300, bg = "white")

m2 <- melt(cmp[, .(label, `260 libraries` = fdr1, `291 libraries` = fdr2)],
           id.vars = "label", variable.name = "run", value.name = "n")
p <- ggplot(m2, aes(label, n, fill = run)) +
  geom_col(position = position_dodge(0.8), width = 0.72) +
  geom_text(aes(label = n), position = position_dodge(0.8), vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = c(`260 libraries` = "#9E9E9E", `291 libraries` = "#1B7837"), name = NULL) +
  labs(title = "Genes at FDR < 0.05 (no fold-change cut): 260 vs 291",
       subtitle = "This is where extra libraries buy power",
       x = NULL, y = "Genes") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 18, hjust = 1), legend.position = "top",
        panel.grid.minor = element_blank())
ggsave(file.path(FIGS, "21_v1_vs_v2_fdr_counts.png"), p, width = 9, height = 7, dpi = 300, bg = "white")

# ---------------------------------------------------------------- concordance
conc <- rbindlist(lapply(names(KEYS), function(k) {
  a <- v1[[k]][, .(gene_id, lfc1 = log2FoldChange, p1 = padj)]
  b <- v2[[k]][, .(gene_id, lfc2 = log2FoldChange, p2 = padj)]
  m <- merge(a, b, by = "gene_id")
  s1g <- m[!is.na(p1) & p1 < FDR_CUT & abs(lfc1) > LFC_CUT, gene_id]
  s2g <- m[!is.na(p2) & p2 < FDR_CUT & abs(lfc2) > LFC_CUT, gene_id]
  data.table(analysis = k, label = KEYS[[k]], genes = nrow(m),
             r = cor(m$lfc1, m$lfc2, use = "complete"),
             sig_v1 = length(s1g), sig_v2 = length(s2g),
             shared = length(intersect(s1g, s2g)),
             pct_v1_kept = round(100 * length(intersect(s1g, s2g)) / max(1, length(s1g)), 1))
}))
fwrite(conc, file.path(TABLES, "v1_vs_v2_concordance.tsv"), sep = "\t")
print(conc)

pl <- rbindlist(lapply(names(KEYS), function(k) {
  m <- merge(v1[[k]][, .(gene_id, lfc1 = log2FoldChange)], v2[[k]][, .(gene_id, lfc2 = log2FoldChange)], by = "gene_id")
  m[, label := KEYS[[k]]]; m
}))
pl[, label := factor(label, levels = unname(KEYS))]
p <- ggplot(pl, aes(lfc1, lfc2)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(size = 0.5, alpha = 0.25, colour = "#1B7837") +
  facet_wrap(~ label, scales = "free", ncol = 3) +
  labs(title = "Fold changes agree between the two runs",
       subtitle = "Each point is a gene. x = 260-library run, y = 291-library run.",
       x = "log2FC, 260 libraries", y = "log2FC, 291 libraries") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92", colour = NA))
ggsave(file.path(FIGS, "22_v1_vs_v2_lfc.png"), p, width = 10, height = 7, dpi = 300, bg = "white")
cat("\nFigures done.\n")
