# annot6_tissue_DE.R — utricle-vs-saccule differential expression WITHIN each
#   shared cell type, using the combined all-samples integrated object (both
#   tissues co-embedded, one consistent cell_type label).  Produces per-cell-type
#   DE tables, volcano plots, a regional-marker panel (Otx1/Otx2/Emx2/...), and a
#   pseudobulk replicate-consistency check.
# setwd("analyses/07_single_cell_rna")   # from the project root

rm(list = ls())
suppressMessages({
    library(Seurat); library(tidyverse); library(patchwork)
    library(pheatmap); library(viridis); library(ggrepel)
})
source("scripts/vestibular_markers.R")

# FindMarkers on the large HC_type_I population (~24k cells) exceeds the default
# 500 MiB future globals cap; raise it so every cell type can be tested.
options(future.globals.maxSize = 8 * 1024^3)   # 8 GiB

outdir <- "results/tissue_DE"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

MIN_CELLS <- 40      # per tissue, within a cell type, to attempt DE
LFC       <- 0.5     # |log2FC| threshold for calling / labelling
PADJ      <- 0.05



### 0. Load combined object; DE must use the RNA (not integrated) assay -------
intg <- readRDS("results/cell_annotation/integrated_all_samples_annotated.rds")
DefaultAssay(intg) <- "RNA"
intg <- tryCatch(JoinLayers(intg), error = function(e) intg)
intg <- NormalizeData(intg, verbose = FALSE)
stopifnot(all(c("tissue", "cell_type", "sample") %in% colnames(intg@meta.data)))
cat("tissues:", paste(unique(intg$tissue), collapse = ", "), "\n")
print(table(intg$cell_type, intg$tissue))



### 1. Per-cell-type DE: utricle vs saccule ---------------------------------
tab <- table(intg$cell_type, intg$tissue)
testable <- rownames(tab)[apply(tab, 1, function(x) all(x >= MIN_CELLS))]
cat("cell types testable (>=", MIN_CELLS, "cells/tissue):",
    paste(testable, collapse = ", "), "\n")

de_all <- list()
for (ct in testable) {
    sub <- subset(intg, subset = cell_type == ct)
    Idents(sub) <- "tissue"
    de <- tryCatch(
        FindMarkers(sub, ident.1 = "utricle", ident.2 = "saccule",
                    min.pct = 0.1, logfc.threshold = 0.1, verbose = FALSE),
        error = function(e) { message("  DE failed for ", ct, ": ", e$message); NULL })
    if (is.null(de) || !nrow(de)) next
    de$gene <- rownames(de); de$cell_type <- ct
    de$direction <- ifelse(de$p_val_adj < PADJ & de$avg_log2FC >  LFC, "up_utricle",
                    ifelse(de$p_val_adj < PADJ & de$avg_log2FC < -LFC, "up_saccule", "ns"))
    de_all[[ct]] <- de
    write.csv(de, file.path(outdir, paste0("DE_utricle_vs_saccule_", ct, ".csv")),
              row.names = FALSE)

    # volcano
    d <- de; d$neglog10 <- -log10(pmax(d$p_val_adj, 1e-300))
    lab <- d %>% filter(direction != "ns") %>% arrange(p_val_adj) %>% head(20)
    reg_here <- d %>% filter(gene %in% regional_markers)
    p <- ggplot(d, aes(avg_log2FC, neglog10)) +
        geom_point(aes(color = direction), size = 0.7, alpha = 0.6) +
        scale_color_manual(values = c(up_utricle = "#c0392b", up_saccule = "#2e86ab", ns = "grey80")) +
        geom_vline(xintercept = c(-LFC, LFC), linetype = "dashed", color = "grey60") +
        geom_hline(yintercept = -log10(PADJ), linetype = "dashed", color = "grey60") +
        geom_text_repel(data = lab, aes(label = gene), size = 2.6, max.overlaps = 25) +
        geom_point(data = reg_here, color = "black", size = 1.6, shape = 21, fill = "gold") +
        geom_text_repel(data = reg_here, aes(label = gene), size = 3, fontface = "bold",
                        color = "black", max.overlaps = 30) +
        labs(title = paste0(ct, " — utricle vs saccule"),
             subtitle = "red = up in utricle · blue = up in saccule · gold = regional candidate",
             x = "avg log2FC (utricle / saccule)", y = "-log10 adj.p") +
        theme_classic()
    ggsave(file.path(outdir, paste0("volcano_", ct, ".png")),
           p, width = 7, height = 6, dpi = 200, bg = "white")
}

de_tbl <- bind_rows(de_all)
write.csv(de_tbl, file.path(outdir, "DE_utricle_vs_saccule_ALL.csv"), row.names = FALSE)

# summary: number of significant genes per cell type / direction
de_summary <- de_tbl %>% filter(direction != "ns") %>%
    count(cell_type, direction) %>% pivot_wider(names_from = direction, values_from = n, values_fill = 0)
write.csv(de_summary, file.path(outdir, "DE_significant_counts_by_celltype.csv"), row.names = FALSE)
cat("---- significant genes per cell type ----\n"); print(as.data.frame(de_summary))



### 2. Regional-marker panel: are Otx1/Otx2/Emx2/... different by tissue? ----
reg <- intersect(regional_markers, rownames(intg))
cat("regional markers present:", paste(reg, collapse = ", "), "\n")
epith <- c("HC_type_I", "Hair_cell", "HC_type_II", "Supporting_cell", "Transitional_epith")
epith <- intersect(epith, unique(intg$cell_type))
sub_e <- subset(intg, subset = cell_type %in% epith)
sub_e$ct_tissue <- paste(sub_e$cell_type, sub_e$tissue, sep = " | ")

if (length(reg)) {
    png(file.path(outdir, "regional_markers_DotPlot_by_tissue.png"),
        width = 2000, height = 1500, res = 200, bg = "white")
    print(DotPlot(sub_e, features = reg, group.by = "ct_tissue") + RotatedAxis() +
          ggtitle("Regional candidates in epithelial types (cell type | tissue)"))
    dev.off()

    # violins of the two headline regional genes across tissue, within epithelium
    for (g in intersect(c("Otx1", "Otx2", "Emx2"), reg)) {
        png(file.path(outdir, paste0("violin_", g, "_by_tissue.png")),
            width = 1800, height = 1100, res = 200, bg = "white")
        print(VlnPlot(sub_e, features = g, group.by = "cell_type", split.by = "tissue",
                      pt.size = 0) + ggtitle(paste0(g, " — utricle vs saccule")))
        dev.off()
    }
}



### 3. Pseudobulk replicate-consistency check (no DE pkgs needed) ------------
# average expression per (cell_type x sample); shows whether tissue differences
# are consistent across the 2 replicates rather than driven by one sample.
avg <- AverageExpression(intg, assays = "RNA", group.by = c("cell_type", "sample"),
                         layer = "data")$RNA
saveRDS(avg, file.path(outdir, "pseudobulk_avg_expr_celltype_sample.rds"))

# for the two most abundant epithelial types, scatter mean(utricle) vs mean(saccule)
smp_tissue <- c(FTL_1146 = "saccule", FTL_1147 = "utricle",
                FTL_1148 = "saccule", FTL_1149 = "utricle")
for (ct in intersect(c("HC_type_I", "Supporting_cell"), testable)) {
    # AverageExpression sanitizes "_" -> "-" in group names (HC_type_I -> HC-type-I,
    # FTL_1146 -> FTL-1146); match on the sanitized form and map samples back.
    ct_safe <- gsub("_", "-", ct)
    cols <- colnames(avg)[startsWith(colnames(avg), paste0(ct_safe, "_"))]
    if (length(cols) < 4) next
    smp <- gsub("-", "_", sub(paste0("^", ct_safe, "_"), "", cols))   # FTL-1146 -> FTL_1146
    ut  <- cols[smp_tissue[smp] == "utricle"]; sa <- cols[smp_tissue[smp] == "saccule"]
    df <- data.frame(gene = rownames(avg),
                     utricle = rowMeans(avg[, ut, drop = FALSE]),
                     saccule = rowMeans(avg[, sa, drop = FALSE]))
    df$l2fc <- log2((df$utricle + 1) / (df$saccule + 1))
    hi <- df %>% filter(abs(l2fc) > 1, pmax(utricle, saccule) > 1) %>%
        arrange(desc(abs(l2fc))) %>% head(20)
    rr <- df %>% filter(gene %in% reg)
    p <- ggplot(df, aes(log2(saccule + 1), log2(utricle + 1))) +
        geom_point(size = 0.5, alpha = 0.4, color = "grey70") +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        geom_point(data = hi, color = "#c0392b", size = 1.2) +
        geom_text_repel(data = hi, aes(label = gene), size = 2.5, max.overlaps = 25) +
        geom_point(data = rr, color = "gold3", size = 2) +
        geom_text_repel(data = rr, aes(label = gene), size = 3, fontface = "bold", color = "gold4") +
        labs(title = paste0(ct, " — pseudobulk mean expression"),
             x = "log2 saccule (mean of reps)", y = "log2 utricle (mean of reps)") +
        theme_classic()
    ggsave(file.path(outdir, paste0("pseudobulk_scatter_", ct, ".png")),
           p, width = 6.5, height = 6, dpi = 200, bg = "white")
}



### 4. Heatmap of top tissue-DE genes across cell types ---------------------
top_genes <- de_tbl %>% filter(direction != "ns") %>%
    group_by(cell_type) %>% slice_max(abs(avg_log2FC), n = 8) %>% ungroup() %>%
    pull(gene) %>% unique()
if (length(top_genes) > 2) {
    ae <- AverageExpression(subset(intg, subset = cell_type %in% testable),
                            assays = "RNA", features = top_genes,
                            group.by = c("cell_type", "tissue"), layer = "data")$RNA
    png(file.path(outdir, "heatmap_top_tissueDE_genes.png"),
        width = 2000, height = min(4000, 200 + 22 * length(top_genes)), res = 200, bg = "white")
    pheatmap(log1p(ae), scale = "row", color = viridis(100),
             cluster_cols = TRUE, fontsize_row = 7,
             main = "Top utricle-vs-saccule DE genes (row z-score of log mean expr)")
    dev.off()
}

message("==== annot6 tissue DE done ====")
### END
