# 16 - What the organ x age interaction tests, drawn out.
#
# Contrast E asks a different question from A/B/C/D. A and B ask "is this gene
# different between the organs?", C and D ask "does this gene change with age?".
# The interaction asks "is the ORGAN DIFFERENCE itself different at the two ages?"
#
#   gap(w8)  = log2( Utricle / Saccule ) at 8 weeks
#   gap(w22) = log2( Utricle / Saccule ) at 22 weeks
#   interaction = gap(w22) - gap(w8)       <- this is what is tested
#
# A gene can be hugely organ-specific and still have zero interaction: that just
# means the organ difference is the same at both ages (two parallel lines).
#
# Two figures:
#   31  six hand-picked genes, one per shape the interaction can take, including
#       one NON-significant gene as the control
#   32  the nine most significant interaction genes
#
# Numbers come from the fitted model rather than from raw group averages, because
# the model works on log2 of the mean count while a raw average of log2(count+1)
# behaves very differently for genes that are near zero in one group.

source("../../analyses/05_dge/scripts/00_config.R")

library(data.table)
library(DESeq2)
library(ggplot2)

dds  <- readRDS(file.path(DATA, "dds_interaction.rds"))
anno <- fread(file.path(DATA, "gene_annotation.tsv"))
e    <- readRDS(file.path(DATA, "de_results.rds"))[["E_organ_by_age_interaction"]]

# With design ~ batch + organ + age + organ:age the four group means are simple
# sums of coefficients, and the interaction is exactly gap(w22) - gap(w8).
cf <- coef(dds)
fit <- data.table(
  gene_id = rownames(cf),
  sac_w8  = cf[, "Intercept"],
  sac_w22 = cf[, "Intercept"] + cf[, "age_w22_vs_w8"],
  utr_w8  = cf[, "Intercept"] + cf[, "organ_Utricle_vs_Saccule"],
  utr_w22 = cf[, "Intercept"] + cf[, "organ_Utricle_vs_Saccule"] +
            cf[, "age_w22_vs_w8"] + cf[, "organUtricle.agew22"],
  gap_w8  = cf[, "organ_Utricle_vs_Saccule"],
  interaction = cf[, "organUtricle.agew22"])
fit[, gap_w22 := gap_w8 + interaction]
fit <- merge(fit, e[, .(gene_id, gene_name, baseMean, padj)], by = "gene_id")

fdr_lab <- function(p) ifelse(p < 0.001, "FDR < 0.001", sprintf("FDR = %s", signif(p, 2)))

# ---------------------------------------------------------------- the drawing
# One panel per gene: the four fitted group means, the two organ lines, and a
# grey bar marking the gap that is being compared at each age.
draw_interaction <- function(d, title, subtitle, file, ncol, width, height) {
  d[, panel := factor(panel, levels = panel)]

  long <- melt(d[, .(panel, Saccule_w8 = sac_w8, Saccule_w22 = sac_w22,
                     Utricle_w8 = utr_w8, Utricle_w22 = utr_w22)],
               id.vars = "panel", variable.name = "key", value.name = "expr")
  long[, organ := sub("_.*$", "", key)]
  long[, age   := factor(sub("^.*_", "", key), levels = c("w8", "w22"))]

  seg <- rbind(
    d[, .(panel, age = factor("w8",  levels = c("w8", "w22")),
          lo = pmin(sac_w8, utr_w8),   hi = pmax(sac_w8, utr_w8),   lab = gap_w8)],
    d[, .(panel, age = factor("w22", levels = c("w8", "w22")),
          lo = pmin(sac_w22, utr_w22), hi = pmax(sac_w22, utr_w22), lab = gap_w22)])

  p <- ggplot(long, aes(age, expr, colour = organ, group = organ)) +
    geom_segment(data = seg, inherit.aes = FALSE,
                 aes(x = age, xend = age, y = lo, yend = hi),
                 colour = "grey45", linewidth = 2.6, alpha = 0.28) +
    geom_text(data = seg, inherit.aes = FALSE,
              aes(x = age, y = (lo + hi) / 2, label = sprintf("gap\n%+.2f", lab)),
              colour = "grey25", size = 2.8, fontface = "bold") +
    geom_line(linewidth = 0.9) +
    geom_point(size = 3.2) +
    scale_colour_manual(values = ORGAN_COLS) +
    facet_wrap(~ panel, scales = "free_y", ncol = ncol) +
    labs(title = title, subtitle = subtitle,
         x = NULL, y = "Model-fitted expression, log2 scale", colour = "Organ") +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text = element_text(size = 9.5, lineheight = 1.05),
          legend.position = "top")

  ggsave(file.path(FIGS, file), p, width = width, height = height, dpi = 300, bg = "white")
  p
}

# ================================================================ figure 31
# One gene for each shape the interaction can take. Gpnmb is deliberately NOT
# significant - it is the control that shows a large, constant organ difference
# scoring zero on this test. The other five all pass FDR < 0.05 and |log2FC| > 1.
show <- data.table(
  gene_name = c("Gpnmb",  "Kcne1", "Jchain", "S100a8", "Slc14a1", "Cxcl10"),
  shape     = c("gap unchanged", "gap widens (utricle)", "gap widens (saccule)",
                "gap narrows", "gap narrows", "gap reverses sign"),
  meaning   = c("No interaction: the organ difference is the same at both ages",
                "Utricle pulls further ahead of saccule with age",
                "Saccule pulls further ahead of utricle with age",
                "The two organs converge with age",
                "The two organs converge with age",
                "Utricle-high at w8 becomes saccule-high at w22"))
d31 <- merge(show, fit, by = "gene_name", sort = FALSE)
d31[, panel := paste0(gene_name, " - ", shape, "\ninteraction = ",
                      sprintf("%+.2f", interaction), "   (", fdr_lab(padj), ")")]

draw_interaction(
  d31,
  "What the organ x age interaction tests",
  paste("The grey bar is the organ gap, log2(Utricle / Saccule), at each age.",
        "The interaction is the CHANGE in that gap:\ninteraction = gap(w22) - gap(w8).",
        "It is zero whenever the two lines are parallel, however far apart they are.",
        "\nGpnmb is the control: a large organ difference that does NOT change with age.",
        "The other five are significant."),
  "31_interaction_illustration.png", ncol = 3, width = 11, height = 7.3)

# ================================================================ figure 32
# The nine most significant interaction genes, same drawing.
sig <- e[!is.na(padj) & padj < FDR_CUT & abs(log2FoldChange) > LFC_CUT]
top9 <- head(sig[order(padj), gene_id], 9)
d32  <- fit[gene_id %in% top9]
setorder(d32, padj)
d32[, panel := paste0(gene_name, "   interaction = ", sprintf("%+.2f", interaction),
                      "\n", fdr_lab(padj), "   |   mean count ", round(baseMean))]

draw_interaction(
  d32,
  "The nine most significant organ x age interaction genes",
  paste("Ranked by adjusted p among the", nrow(sig),
        "genes passing FDR < 0.05 and |log2FC| > 1. Same drawing as the previous figure:",
        "\ngrey bar = log2(Utricle / Saccule) at that age; the interaction is how much that gap moved.",
        "\nMean count is each gene's baseMean - several of the top-ranked genes are lowly expressed."),
  "32_interaction_top9.png", ncol = 3, width = 11, height = 8.4)

# ---------------------------------------------------------------- tables
tab <- d31[, .(Gene = gene_name,
               `Utricle/Saccule at w8`  = round(gap_w8, 2),
               `Utricle/Saccule at w22` = round(gap_w22, 2),
               `Interaction (change in the gap)` = round(interaction, 2),
               `FDR` = signif(padj, 2),
               `What it means` = meaning)]
fwrite(tab, file.path(TABLES, "interaction_examples.tsv"), sep = "\t")
print(tab)

tab9 <- d32[, .(Gene = gene_name,
                `Mean count` = round(baseMean),
                `gap w8` = round(gap_w8, 2), `gap w22` = round(gap_w22, 2),
                Interaction = round(interaction, 2), FDR = signif(padj, 2))]
fwrite(tab9, file.path(TABLES, "interaction_top9.tsv"), sep = "\t")
cat("\n=== top 9 significant interaction genes ===\n"); print(tab9)

cat("\nInteraction genes at FDR < 0.05 and |log2FC| > 1:", nrow(sig), "of", nrow(e), "tested\n")
