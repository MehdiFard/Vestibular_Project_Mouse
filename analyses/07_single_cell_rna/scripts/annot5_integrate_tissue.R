# annot5_integrate_tissue.R — integrate the 2 replicates of ONE tissue, annotate,
#   and find per-cell-type markers.  Run once per tissue:
#     Rscript annot5_integrate_tissue.R saccule
#     Rscript annot5_integrate_tissue.R utricle
# setwd("analyses/07_single_cell_rna")   # from the project root

rm(list = ls())
suppressMessages({
    library(Seurat); library(tidyverse); library(patchwork)
    library(pheatmap); library(viridis)
})
source("scripts/vestibular_markers.R")

args   <- commandArgs(trailingOnly = TRUE)
tissue <- if (length(args)) args[1] else stop("give tissue: saccule | utricle")
stopifnot(tissue %in% c("saccule", "utricle"))

# sample -> replicate map per tissue (see scRNA_integration.R)
sample_map <- if (tissue == "saccule") {
    c(FTL_1146 = "rep1", FTL_1148 = "rep2")
} else {
    c(FTL_1147 = "rep1", FTL_1149 = "rep2")
}

outdir    <- "results/per_tissue"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
clust_res <- "integrated_snn_res.0.3"
message("==== integrating ", tissue, ": ", paste(names(sample_map), collapse = ", "), " ====")



### 1. Read + per-sample QC / normalize ------------------------------------
obj.list <- list()
for (s in names(sample_map)) {
    m <- Read10X_h5(file.path("intermediate/cell_ranger_output", s,
                              "outs/filtered_feature_bc_matrix.h5"),
                    use.names = TRUE, unique.features = TRUE)
    o <- CreateSeuratObject(m, project = s, min.cells = 3, min.features = 200)
    o$mitoPercent <- PercentageFeatureSet(o, pattern = "^mt-")     # mouse: mt-
    o <- subset(o, subset = nCount_RNA > 800 & nFeature_RNA > 500 & mitoPercent < 5)
    o$sample <- s; o$tissue <- tissue; o$replicate <- sample_map[[s]]
    o <- NormalizeData(o, verbose = FALSE)
    o <- FindVariableFeatures(o, nfeatures = 2000, verbose = FALSE)
    obj.list[[s]] <- o
}



### 2. CCA integration (batch-correct the 2 replicates) --------------------
features <- SelectIntegrationFeatures(object.list = obj.list)
anchors  <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)
intg     <- IntegrateData(anchorset = anchors)

DefaultAssay(intg) <- "integrated"
intg <- ScaleData(intg, verbose = FALSE)
intg <- RunPCA(intg, npcs = 30, verbose = FALSE)
intg <- FindNeighbors(intg, dims = 1:30)
intg <- FindClusters(intg, resolution = c(0.1, 0.3, 0.5))
intg <- RunUMAP(intg, dims = 1:30)



### 3. Marker-based annotation (scored on RNA) -----------------------------
DefaultAssay(intg) <- "RNA"
intg <- tryCatch(JoinLayers(intg), error = function(e) intg)
intg <- NormalizeData(intg, verbose = FALSE)

ann <- annotate_by_markers(intg, clust_res, assay = "RNA")
intg <- ann$seu
cat("---- ", tissue, " cluster -> cell type ----\n"); print(ann$labels)
cat("---- cell-type counts ----\n"); print(table(intg$cell_type))



### 4. Per-cell-type marker genes (DE that defines each type) --------------
Idents(intg) <- "cell_type"
ct_markers <- FindAllMarkers(intg, only.pos = TRUE, min.pct = 0.25,
                             logfc.threshold = 0.25, verbose = FALSE)
readr::write_tsv(ct_markers, file.path(outdir, paste0(tissue, "_celltype_markers.tsv")))
top_ct <- ct_markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10) %>% ungroup()
readr::write_tsv(top_ct, file.path(outdir, paste0(tissue, "_celltype_markers_top10.tsv")))

# composition per replicate
comp <- intg@meta.data %>% count(tissue, sample, replicate, cell_type) %>%
    group_by(sample) %>% mutate(prop = n / sum(n)) %>% ungroup()
readr::write_csv(comp, file.path(outdir, paste0(tissue, "_composition_by_replicate.csv")))



### 4b. Save the annotated object NOW (before plotting) --------------------
saveRDS(intg, file.path(outdir, paste0(tissue, "_integrated_annotated.rds")))



### 5. Plots ---------------------------------------------------------------
p_ct <- DimPlot(intg, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE) +
            ggtitle(paste0(tissue, " — marker cell types"))
p_rep <- DimPlot(intg, reduction = "umap", group.by = "replicate") + ggtitle("replicate (batch)")
png(file.path(outdir, paste0(tissue, "_UMAP_celltype.png")),
    width = 1800, height = 1500, res = 200, bg = "white"); print(p_ct); dev.off()
png(file.path(outdir, paste0(tissue, "_UMAP_replicate.png")),
    width = 1600, height = 1400, res = 200, bg = "white"); print(p_rep); dev.off()

png(file.path(outdir, paste0(tissue, "_DotPlot_markers.png")),
    width = 2000, height = 1400, res = 200, bg = "white")
print(DotPlot(intg, features = intersect(rep_markers, rownames(intg)),
              group.by = "cell_type") + RotatedAxis() +
      ggtitle(paste0(tissue, " — representative markers")))
dev.off()

p_bar <- ggplot(comp, aes(sample, prop, fill = cell_type)) + geom_col() +
    labs(y = "proportion of cells", x = NULL, title = paste0(tissue, " — composition")) +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(outdir, paste0(tissue, "_composition_bar.png")),
       p_bar, width = 6, height = 6, dpi = 200, bg = "white")



message("==== done ", tissue, " ====")

### END
