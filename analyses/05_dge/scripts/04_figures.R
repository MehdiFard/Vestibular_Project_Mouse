# 04 - Figures for the differential expression results.
#
# Volcano plots per contrast, then the cross-contrast comparisons:
#   - is the organ difference the same at w8 and w22?
#   - is the ageing response the same in the two organs?

source("../../analyses/05_dge/scripts/00_config.R")

library(data.table)
library(ggplot2)

de      <- readRDS(file.path(DATA, "de_results.rds"))
summ    <- fread(file.path(TABLES, "DE_summary.tsv"))
samples <- fread(file.path(DATA, "sample_table.tsv"))

# Short display names for each contrast
NICE <- c(
  A_utricle_vs_saccule_w8   = "Utricle vs Saccule, w8",
  B_utricle_vs_saccule_w22  = "Utricle vs Saccule, w22",
  C_w22_vs_w8_saccule       = "w22 vs w8, Saccule",
  D_w22_vs_w8_utricle       = "w22 vs w8, Utricle",
  E_organ_by_age_interaction = "Organ x age interaction")

label_of <- function(res, n = 12) {
  # Genes to annotate on a volcano. Half are the most significant and half have
  # the largest fold change, split evenly between the two directions, so the
  # labels spread out instead of stacking on the p-value cap.
  s <- res[!is.na(padj) & padj < FDR_CUT & abs(log2FoldChange) > LFC_CUT &
             !is.na(gene_name) & gene_name != ""]
  if (nrow(s) == 0) return(s)
  s[, abs_lfc := abs(log2FoldChange)]

  pick <- function(d, col, k) {
    if (nrow(d) == 0) return(character())
    setorderv(d, col, if (col == "padj") 1L else -1L)
    head(d$gene_id, k)
  }
  k <- max(1L, n %/% 4L)
  ids <- unique(c(
    pick(s[log2FoldChange > 0], "padj", k), pick(s[log2FoldChange < 0], "padj", k),
    pick(s[log2FoldChange > 0], "abs_lfc", k), pick(s[log2FoldChange < 0], "abs_lfc", k)))
  out <- s[gene_id %in% ids]
  out[, abs_lfc := NULL]
  out[]
}

# ---------------------------------------------------------------- volcano
volcano <- function(res, title, up_lab, down_lab, file, n_label = 12) {
  d <- copy(res)[!is.na(padj)]
  d[, status := "not significant"]
  d[padj < FDR_CUT & log2FoldChange >  LFC_CUT, status := up_lab]
  d[padj < FDR_CUT & log2FoldChange < -LFC_CUT, status := down_lab]
  d[, status := factor(status, levels = c(down_lab, "not significant", up_lab))]

  cols <- setNames(c("#2166AC", "grey75", "#B2182B"),
                   c(down_lab, "not significant", up_lab))

  # Cap -log10(p) so a few astronomically significant genes don't flatten the rest
  cap <- quantile(-log10(d$padj), 0.999, na.rm = TRUE)
  d[, y := pmin(-log10(padj), cap)]
  lab <- label_of(res, n_label)
  lab[, y := pmin(-log10(padj), cap)]

  n_up   <- sum(d$status == up_lab)
  n_down <- sum(d$status == down_lab)

  p <- ggplot(d, aes(log2FoldChange, y, colour = status)) +
    geom_hline(yintercept = -log10(FDR_CUT), linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", colour = "grey40") +
    geom_point(size = 0.9, alpha = 0.6) +
    geom_text(data = lab, aes(label = gene_name), size = 2.9, colour = "black",
              vjust = -0.7, check_overlap = TRUE, show.legend = FALSE) +
    scale_colour_manual(values = cols) +
    labs(title = title,
         subtitle = sprintf("FDR < %.2f and |log2FC| > %g   |   %s: %d    %s: %d",
                            FDR_CUT, LFC_CUT, up_lab, n_up, down_lab, n_down),
         x = "log2 fold change (shrunken)", y = "-log10 adjusted p", colour = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top", panel.grid.minor = element_blank())

  ggsave(file.path(FIGS, file), p, width = 8, height = 6.5, dpi = 300, bg = "white")
  p
}

volcano(de$A_utricle_vs_saccule_w8,  "Utricle vs Saccule at 8 weeks (paired by strain)",
        "up in Utricle", "up in Saccule", "10_volcano_utricle_vs_saccule_w8.png")
volcano(de$B_utricle_vs_saccule_w22, "Utricle vs Saccule at 22 weeks (paired by strain)",
        "up in Utricle", "up in Saccule", "11_volcano_utricle_vs_saccule_w22.png")
volcano(de$C_w22_vs_w8_saccule, "22 weeks vs 8 weeks in Saccule (batch-adjusted)",
        "up at w22", "up at w8", "12_volcano_w22_vs_w8_saccule.png")
volcano(de$D_w22_vs_w8_utricle, "22 weeks vs 8 weeks in Utricle (batch-adjusted)",
        "up at w22", "up at w8", "13_volcano_w22_vs_w8_utricle.png")
volcano(de$E_organ_by_age_interaction,
        "Organ x age interaction",
        "organ gap widens with age", "organ gap narrows with age",
        "14_volcano_interaction.png")

# ---------------------------------------------------------------- summary bars
bar_dt <- summ[analysis %in% names(NICE)]
bar_dt[, label := NICE[analysis]]
bar_long <- melt(bar_dt[, .(label, up, down)], id.vars = "label",
                 variable.name = "direction", value.name = "n")
bar_long[, label := factor(label, levels = rev(NICE))]

p_bar <- ggplot(bar_long, aes(n, label, fill = direction)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = n), position = position_dodge(width = 0.75),
            hjust = -0.15, size = 3.3) +
  scale_fill_manual(values = c(up = "#B2182B", down = "#2166AC"),
                    labels = c(up = "positive log2FC", down = "negative log2FC")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Significant genes per contrast",
       subtitle = sprintf("FDR < %.2f and |log2FC| > %g", FDR_CUT, LFC_CUT),
       x = "Number of genes", y = NULL, fill = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

ggsave(file.path(FIGS, "15_significant_gene_counts.png"), p_bar,
       width = 9, height = 5, dpi = 300, bg = "white")

# ---------------------------------------------------------------- concordance
# Scatter of log2FC from two contrasts over the genes tested in both.
concordance <- function(res1, res2, lab1, lab2, title, file) {
  m <- merge(res1[, .(gene_id, gene_name, lfc1 = log2FoldChange, p1 = padj)],
             res2[, .(gene_id,            lfc2 = log2FoldChange, p2 = padj)],
             by = "gene_id")
  m <- m[!is.na(lfc1) & !is.na(lfc2)]
  m[, sig := "neither"]
  m[!is.na(p1) & p1 < FDR_CUT & abs(lfc1) > LFC_CUT, sig := paste("only", lab1)]
  m[!is.na(p2) & p2 < FDR_CUT & abs(lfc2) > LFC_CUT, sig := paste("only", lab2)]
  m[!is.na(p1) & p1 < FDR_CUT & abs(lfc1) > LFC_CUT &
    !is.na(p2) & p2 < FDR_CUT & abs(lfc2) > LFC_CUT, sig := "both"]

  lv <- c("neither", paste("only", lab1), paste("only", lab2), "both")
  m[, sig := factor(sig, levels = lv)]
  cols <- setNames(c("grey80", "#4393C3", "#D6604D", "#1A9850"), lv)

  r_all  <- cor(m$lfc1, m$lfc2, method = "pearson")
  both   <- m[sig == "both"]
  r_both <- if (nrow(both) > 2) cor(both$lfc1, both$lfc2) else NA_real_

  p <- ggplot(m[order(sig)], aes(lfc1, lfc2, colour = sig)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = 0, colour = "grey85") +
    geom_vline(xintercept = 0, colour = "grey85") +
    geom_point(size = 0.9, alpha = 0.6) +
    scale_colour_manual(values = cols) +
    labs(title = title,
         subtitle = sprintf("r = %.3f over %d genes tested in both\nr = %s over the %d significant in both",
                            r_all, nrow(m),
                            ifelse(is.na(r_both), "NA", sprintf("%.3f", r_both)), nrow(both)),
         x = paste("log2FC -", lab1), y = paste("log2FC -", lab2), colour = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top", panel.grid.minor = element_blank())

  ggsave(file.path(FIGS, file), p, width = 7.5, height = 7, dpi = 300, bg = "white")
  list(plot = p, merged = m, r_all = r_all, r_both = r_both,
       n_both = nrow(both),
       n_only1 = sum(m$sig == paste("only", lab1)),
       n_only2 = sum(m$sig == paste("only", lab2)))
}

c_organ <- concordance(de$A_utricle_vs_saccule_w8, de$B_utricle_vs_saccule_w22,
                       "w8", "w22",
                       "Organ difference (Utricle vs Saccule) at the two ages",
                       "16_concordance_organ_w8_vs_w22.png")

c_age <- concordance(de$C_w22_vs_w8_saccule, de$D_w22_vs_w8_utricle,
                     "Saccule", "Utricle",
                     "Ageing response (w22 vs w8) in the two organs",
                     "17_concordance_age_saccule_vs_utricle.png")

fwrite(c_organ$merged, file.path(TABLES, "organ_effect_w8_vs_w22_merged.tsv.gz"),
       sep = "\t", na = "NA")
fwrite(c_age$merged, file.path(TABLES, "age_effect_saccule_vs_utricle_merged.tsv.gz"),
       sep = "\t", na = "NA")

overlap <- data.table(
  comparison = c("Organ effect: w8 vs w22", "Age effect: Saccule vs Utricle"),
  r_all_genes      = c(c_organ$r_all,  c_age$r_all),
  r_shared_sig     = c(c_organ$r_both, c_age$r_both),
  sig_in_both      = c(c_organ$n_both, c_age$n_both),
  sig_only_first   = c(c_organ$n_only1, c_age$n_only1),
  sig_only_second  = c(c_organ$n_only2, c_age$n_only2))
fwrite(overlap, file.path(TABLES, "contrast_overlap_summary.tsv"), sep = "\t")

# ---------------------------------------------------------------- shared vs specific
ov_long <- rbindlist(list(
  data.table(comparison = "Organ effect (Utricle vs Saccule)",
             set = c("w8 only", "shared", "w22 only"),
             n   = c(c_organ$n_only1, c_organ$n_both, c_organ$n_only2)),
  data.table(comparison = "Age effect (w22 vs w8)",
             set = c("Saccule only", "shared", "Utricle only"),
             n   = c(c_age$n_only1, c_age$n_both, c_age$n_only2))))

p_ov <- ggplot(ov_long, aes(set, n, fill = set)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = n), vjust = -0.35, size = 3.4) +
  facet_wrap(~ comparison, scales = "free") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  scale_fill_manual(values = c("w8 only" = "#4393C3", "w22 only" = "#D6604D",
                               "Saccule only" = "#1B7837", "Utricle only" = "#762A83",
                               "shared" = "#F4A582")) +
  labs(title = "Shared and condition-specific significant genes",
       x = NULL, y = "Number of genes") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 15, hjust = 1))

ggsave(file.path(FIGS, "18_shared_vs_specific.png"), p_ov,
       width = 9, height = 5, dpi = 300, bg = "white")

# ---------------------------------------------------------------- top gene panels
# Normalised expression of the strongest organ-marker genes, per group.
dds_int <- readRDS(file.path(DATA, "dds_interaction.rds"))
norm    <- DESeq2::counts(dds_int, normalized = TRUE)

plot_top_genes <- function(res, n, title, file) {
  s <- res[!is.na(padj) & padj < FDR_CUT & abs(log2FoldChange) > LFC_CUT]
  setorder(s, padj)
  s <- head(s[gene_id %in% rownames(norm) & !is.na(gene_name) & gene_name != ""], n)
  if (nrow(s) == 0) return(NULL)

  d <- rbindlist(lapply(seq_len(nrow(s)), function(i) {
    data.table(gene = s$gene_name[i],
               library_id = colnames(norm),
               value = log2(norm[s$gene_id[i], ] + 1))
  }))
  d <- merge(d, samples[, .(library_id, organ, age)], by = "library_id")
  d[, gene := factor(gene, levels = s$gene_name)]
  d[, age  := factor(age, levels = c("w8", "w22"))]

  p <- ggplot(d, aes(age, value, fill = organ)) +
    geom_boxplot(outlier.size = 0.4, linewidth = 0.35, position = position_dodge(0.8)) +
    facet_wrap(~ gene, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = ORGAN_COLS) +
    labs(title = title, x = NULL, y = "log2(normalised count + 1)", fill = "Organ") +
    theme_bw(base_size = 11) +
    theme(legend.position = "top", panel.grid.minor = element_blank(),
          strip.text = element_text(face = "italic"))

  ggsave(file.path(FIGS, file), p, width = 10, height = 6.5, dpi = 300, bg = "white")
  p
}

plot_top_genes(de$A_utricle_vs_saccule_w8, 8,
               "Top organ-discriminating genes (from the w8 contrast)",
               "19_top_organ_genes.png")
plot_top_genes(de$E_organ_by_age_interaction, 8,
               "Top organ x age interaction genes",
               "20_top_interaction_genes.png")

# ---------------------------------------------------------------- paired vs unpaired
# What the strain pairing bought us, at w8.
pu <- merge(de$A_utricle_vs_saccule_w8[, .(gene_id, paired = -log10(padj))],
            de$S_unpaired_utricle_vs_saccule_w8[, .(gene_id, unpaired = -log10(padj))],
            by = "gene_id")
pu <- pu[is.finite(paired) & is.finite(unpaired)]

p_pu <- ggplot(pu, aes(unpaired, paired)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(size = 0.8, alpha = 0.4, colour = "#4C72B0") +
  labs(title = "Pairing by strain sharpens the organ contrast (w8)",
       subtitle = sprintf(paste("points above the line gained power from pairing;",
                                "%d vs %d genes at FDR < 0.05 (%d vs %d once |log2FC| > 1 is also required)"),
                          summ[analysis == "A_utricle_vs_saccule_w8", sig_fdr],
                          summ[analysis == "S_unpaired_utricle_vs_saccule_w8", sig_fdr],
                          summ[analysis == "A_utricle_vs_saccule_w8", sig_fdr_lfc],
                          summ[analysis == "S_unpaired_utricle_vs_saccule_w8", sig_fdr_lfc]),
       x = "-log10 adjusted p, ~ batch + organ (unpaired)",
       y = "-log10 adjusted p, ~ strain + organ (paired)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(FIGS, "21_paired_vs_unpaired.png"), p_pu,
       width = 7, height = 6.5, dpi = 300, bg = "white")

cat("\nFigures written to", FIGS, "\n")
print(overlap)
