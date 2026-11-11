# 02 - Figures for the ORA results, plus the direct comparison against GSEA.
# Square panels throughout.

BASE <- "../../../analyses/05_dge"
ROOT <- file.path(dirname(BASE), "06_enrichment")
PA   <- file.path(ROOT, "over_rep")
Sys.setenv(TMPDIR = file.path(ROOT, "logs"))

library(data.table)
library(ggplot2)
library(msigdbr)

FDR_CUT <- 0.05
FIGS <- file.path(PA, "figures")

ora  <- readRDS(file.path(PA, "tables", "ora_objects.rds"))
oc   <- ora$combined
cnts <- ora$counts
gsea <- fread(file.path(ROOT, "gsea", "tables", "gsea_all_contrasts.tsv.gz"))

CONTRAST_ORDER <- c("A_organ_w8", "B_organ_w22", "C_age_sac", "D_age_utr", "E_interaction")
COLL_NAME <- c(hallmark = "Hallmark", go_bp = "GO Biological Process", reactome = "Reactome")
oc[, contrast := factor(contrast, levels = CONTRAST_ORDER)]
gsea[, contrast := factor(contrast, levels = CONTRAST_ORDER)]

trim <- function(x, n = 52) ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "…"), x)
theme_sq <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey30", size = 10))
DIR_COLS <- c(up = "#B2182B", down = "#2166AC", `up or down` = "grey45")

# ================================================================ 1. ORA bars
# Top terms from the up and the down list, per contrast per collection.
for (cid in CONTRAST_ORDER) {
  for (cn in names(COLL_NAME)) {
    d <- oc[contrast == cid & collection == cn & gene_list != "all" & p.adjust < FDR_CUT]
    if (!nrow(d)) next
    top <- rbind(head(d[gene_list == "up"][order(p.adjust)], 10),
                 head(d[gene_list == "down"][order(p.adjust)], 10))
    top[, score := -log10(p.adjust) * fifelse(gene_list == "up", 1, -1)]
    top[, lab := trim(pathway_label)]
    top[, lab := factor(lab, levels = unique(top[order(score), lab]))]
    p <- ggplot(top, aes(score, lab, fill = gene_list)) +
      geom_col(width = 0.72) +
      geom_vline(xintercept = 0, colour = "grey30") +
      scale_fill_manual(values = DIR_COLS, guide = "none") +
      labs(title = oc[contrast == cid, contrast_label[1]],
           subtitle = paste0(COLL_NAME[cn], " - ORA on the DE genes, top 10 per direction\n",
                             "right = enriched among ", d[gene_list == "up", direction[1]],
                             "   |   left = among ", d[gene_list == "down", direction[1]]),
           x = "-log10 FDR   (signed by which gene list)", y = NULL) +
      theme_sq + theme(axis.text.y = element_text(size = 9))
    ggsave(file.path(FIGS, cn, paste0(cid, "_ora_bar.png")), p,
           width = 9, height = 9, dpi = 300, bg = "white")
  }
}
cat("ORA bar charts done\n")

# ================================================================ 2. input sizes
m <- melt(cnts, id.vars = c("contrast", "contrast_label", "universe"),
          measure.vars = c("n_up", "n_down"), variable.name = "dir", value.name = "n")
m[, dir := factor(fifelse(dir == "n_up", "up", "down"), levels = c("down", "up"))]
m[, contrast := factor(contrast, levels = CONTRAST_ORDER)]
p <- ggplot(m, aes(contrast, n, fill = dir)) +
  geom_col(width = 0.7) +
  geom_text(data = cnts[, .(contrast = factor(contrast, CONTRAST_ORDER), n_all)],
            aes(contrast, n_all, label = n_all), inherit.aes = FALSE,
            vjust = -0.5, size = 3.6, fontface = "bold") +
  scale_fill_manual(values = DIR_COLS, name = NULL) +
  labs(title = "How many genes ORA actually gets to work with",
       subtitle = paste0("DE genes at FDR < ", FDR_CUT, " and |log2FC| > 1, out of ~17,000-18,000 tested.\n",
                         "GSEA uses the whole ranked list; ORA sees only these."),
       x = NULL, y = "DE genes") +
  theme_sq + theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "top")
ggsave(file.path(FIGS, "summary", "01_de_gene_counts.png"), p,
       width = 9, height = 9, dpi = 300, bg = "white")

# ================================================================ 3. ORA vs GSEA counts
g_sig <- gsea[padj < FDR_CUT, .(n = .N), by = .(contrast, collection)][, method := "GSEA"]
o_sig <- oc[gene_list == "all" & p.adjust < FDR_CUT, .(n = .N), by = .(contrast, collection)][, method := "ORA"]
cmp <- rbind(g_sig, o_sig)
# fill in the combinations where a method found nothing, so the bar reads 0 not absent
cmp <- merge(CJ(contrast = factor(CONTRAST_ORDER, levels = CONTRAST_ORDER),
                collection = names(COLL_NAME), method = c("GSEA", "ORA"), unique = TRUE),
             cmp, by = c("contrast", "collection", "method"), all.x = TRUE)
cmp[is.na(n), n := 0]
cmp[, collection := factor(COLL_NAME[collection], levels = COLL_NAME)]
p <- ggplot(cmp, aes(contrast, n, fill = method)) +
  geom_col(position = position_dodge(0.78), width = 0.72) +
  facet_wrap(~ collection, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(GSEA = "#1B7837", ORA = "#762A83"), name = NULL) +
  labs(title = "GSEA finds more, except where the DE list is large",
       subtitle = paste0("Gene sets at FDR < ", FDR_CUT, ". ORA run on the combined DE list.\n",
                         "GSEA wins in Hallmark and Reactome everywhere. GO:BP is mixed: ORA wins ",
                         "in contrast A,\nwhose 756-gene DE list is big enough to hit many small, ",
                         "overlapping GO terms."),
       x = NULL, y = "Significant gene sets") +
  theme_sq + theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "top")
ggsave(file.path(FIGS, "summary", "02_ora_vs_gsea_counts.png"), p,
       width = 9, height = 9, dpi = 300, bg = "white")
fwrite(dcast(cmp, contrast + collection ~ method, value.var = "n"),
       file.path(PA, "tables", "ora_vs_gsea_counts.tsv"), sep = "\t")

# ================================================================ 4. worked example
# Why ORA misses the strongest GSEA result in the whole analysis.
de <- readRDS(file.path(BASE, "data", "de_results.rds"))
d  <- as.data.table(de[["D_w22_vs_w8_utricle"]])[!is.na(padj)]
h  <- as.data.table(msigdbr(db_species = "MM", species = "Mus musculus", collection = "MH"))
ox <- unique(h[gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION", ensembl_gene])
d[, grp := fifelse(gene_id %in% ox, "Oxidative phosphorylation", "All other tested genes")]
d[, grp := factor(grp, levels = c("All other tested genes", "Oxidative phosphorylation"))]
n_set <- d[grp != "All other tested genes", .N]
n_fdr <- d[grp != "All other tested genes" & padj < FDR_CUT, .N]

p <- ggplot(d, aes(log2FoldChange, fill = grp)) +
  geom_density(alpha = 0.55, colour = NA, adjust = 1.2) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "#B2182B", linewidth = 0.7) +
  geom_vline(xintercept = 0, colour = "grey55") +
  annotate("text", x = 1.05, y = Inf, label = "the |log2FC| > 1 cut", hjust = 0, vjust = 1.6,
           size = 3.4, colour = "#B2182B") +
  coord_cartesian(xlim = c(-2.2, 2.2)) +
  scale_fill_manual(values = c("grey70", "#2166AC"), name = NULL) +
  labs(title = "Why ORA misses what GSEA finds",
       subtitle = paste0("Contrast D (w22 vs w8, Utricle). Of the ", n_set,
                         " oxidative phosphorylation genes tested, ", n_fdr,
                         " pass FDR < 0.05\nbut NOT ONE passes |log2FC| > 1. ",
                         "88% move down, median log2FC -0.13.\n",
                         "GSEA: NES -3.28, FDR 2e-37 - the strongest result in the analysis. ",
                         "ORA: nothing, because\nnone of these genes ever enter the DE list."),
       x = "log2 fold change", y = "Density") +
  theme_sq + theme(legend.position = "top",
                   plot.subtitle = element_text(size = 9.5, lineheight = 1.2))
ggsave(file.path(FIGS, "summary", "03_worked_example_oxphos.png"), p,
       width = 9, height = 9, dpi = 300, bg = "white")

# ================================================================ 5. agreement
# For Hallmark, do the two methods agree on which sets matter?
g <- gsea[collection == "hallmark", .(contrast, pathway, pathway_label, NES, gsea_padj = padj)]
o <- oc[collection == "hallmark" & gene_list == "all",
        .(contrast, pathway, ora_padj = p.adjust, fold_enrichment)]
mg <- merge(g, o, by = c("contrast", "pathway"), all.x = TRUE)
mg[is.na(ora_padj), ora_padj := 1]
mg[, call := fifelse(gsea_padj < FDR_CUT & ora_padj < FDR_CUT, "both",
              fifelse(gsea_padj < FDR_CUT, "GSEA only",
               fifelse(ora_padj < FDR_CUT, "ORA only", "neither")))]
mg[, call := factor(call, levels = c("neither", "ORA only", "GSEA only", "both"))]
fwrite(mg, file.path(PA, "tables", "hallmark_ora_vs_gsea.tsv"), sep = "\t")

p <- ggplot(mg[order(call)], aes(NES, -log10(ora_padj), colour = call)) +
  geom_hline(yintercept = -log10(FDR_CUT), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(size = 2.2, alpha = 0.85) +
  scale_colour_manual(values = c(neither = "grey80", `ORA only` = "#762A83",
                                 `GSEA only` = "#1B7837", both = "#D95F02"), name = NULL) +
  labs(title = "Where the two methods agree, per Hallmark set",
       subtitle = paste0("Each point is one Hallmark set in one contrast (", nrow(mg), " in total).\n",
                         "x = GSEA NES, y = ORA significance. ",
                         mg[call == "both", .N], " called by both, ",
                         mg[call == "GSEA only", .N], " by GSEA only, ",
                         mg[call == "ORA only", .N], " by ORA only."),
       x = "GSEA normalised enrichment score", y = "-log10 ORA FDR") +
  theme_sq + theme(legend.position = "top")
ggsave(file.path(FIGS, "summary", "04_hallmark_agreement.png"), p,
       width = 9, height = 9, dpi = 300, bg = "white")

cat("summary figures done\n")
print(mg[, .N, by = call])
print(list.files(FIGS, recursive = TRUE))
