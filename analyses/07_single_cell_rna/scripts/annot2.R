# script to annotate cell types in the mouse utricle by canonical markers
# sample: FTL_1147 (utricle_1)
# approach: cluster -> module-score canonical vestibular markers -> label clusters
# setwd("analyses/07_single_cell_rna")   # from the project root

rm(list = ls())

library(Seurat)
library(tidyverse)
library(pheatmap)
library(viridis)

s      <- "FTL_1147"
path   <- file.path("intermediate/cell_ranger_output", s,
                    "outs/filtered_feature_bc_matrix.h5")
outdir <- "results/cell_annotation"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# resolution used to define clusters that then get a cell-type label
clust_res <- "RNA_snn_res.0.3"



### 1. Load, QC-filter, and pre-process -------------------------------------
mtx <- Read10X_h5(path, use.names = TRUE, unique.features = TRUE)
seu <- CreateSeuratObject(counts = mtx, project = s,
                          min.cells = 3, min.features = 200)

# mouse mitochondrial genes are lowercase-prefixed "mt-" (NOT "MT-")
seu$mitoPercent <- PercentageFeatureSet(seu, pattern = "^mt-")
seu <- subset(seu, subset = nCount_RNA > 800 &
                            nFeature_RNA > 500 &
                            mitoPercent < 5)
seu$sample <- s
seu$tissue <- "utricle"

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu, features = VariableFeatures(seu), verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:15)
seu <- FindClusters(seu, resolution = c(0.05, 0.1, 0.3, 0.5, 0.7, 1))
seu <- RunUMAP(seu, dims = 1:10)

Idents(seu) <- clust_res



### 2. Canonical marker panels for mouse vestibular macula (utricle) ---------
# sensory-epithelial + non-sensory / stromal cell types expected in a dissected
# macula. all genes below were confirmed present in this dataset.
marker_list <- list(
    Hair_cell            = c("Myo7a", "Pou4f3", "Gfi1", "Myo6", "Ccer2", "Tmc1"),
    HC_type_I            = c("Spp1", "Ocm", "Adgrv1", "Cabp2"),
    HC_type_II           = c("Mapt", "Calb2", "Anxa4"),
    Supporting_cell      = c("Sox2", "Gjb2", "Gjb6", "Hey1", "Hey2", "Otog", "Otogl", "S100b"),
    Transitional_epith   = c("Krt8", "Krt18", "Cldn3", "Cldn6"),
    Roof_dark_cell       = c("Kcnq1", "Kcne1", "Kcnj10", "Atp1b2", "Slc12a2"),
    Glia_Schwann         = c("Plp1", "Mpz", "Sox10", "Mbp"),
    Neuron               = c("Tubb3", "Snap25", "Nefl", "Nefm", "Pvalb", "Stmn2"),
    Immune_macrophage    = c("Ptprc", "C1qa", "C1qb", "Aif1", "Cd68", "Lyz2"),
    Endothelial          = c("Pecam1", "Cldn5", "Flt1", "Cdh5"),
    Melanocyte           = c("Mlana", "Dct", "Pmel", "Tyrp1"),
    Fibrocyte_mesench    = c("Col1a1", "Col1a2", "Dcn", "Lum", "Pdgfra"),
    Pericyte_SMC         = c("Acta2", "Pdgfrb", "Rgs5", "Myh11"),
    Erythrocyte          = c("Hba-a1", "Hbb-bs")
)

# keep only genes present in the object (defensive; drops empty panels)
marker_list <- lapply(marker_list, intersect, y = rownames(seu))
marker_list <- marker_list[lengths(marker_list) > 0]



### 3. Module-score each cell type, then label clusters by argmax ------------
seu <- AddModuleScore(seu, features = marker_list, name = "modscore", seed = 42)
score_cols <- paste0("modscore", seq_along(marker_list))
# rename the auto-generated columns to the cell-type names
seu@meta.data <- seu@meta.data %>%
    rename_with(~ names(marker_list), all_of(score_cols))

# mean module score per cluster
mean_scores <- seu@meta.data %>%
    group_by(cluster = .data[[clust_res]]) %>%
    summarise(across(all_of(names(marker_list)), mean), .groups = "drop")

score_mat <- as.matrix(mean_scores[, names(marker_list)])
rownames(score_mat) <- mean_scores$cluster

# assign each cluster to its highest-scoring cell type
cluster_label <- setNames(colnames(score_mat)[max.col(score_mat, ties.method = "first")],
                          rownames(score_mat))
print(cluster_label)

# unname(): drop the cluster-ID names so Seurat assigns positionally
# (a named vector would be matched against cell barcodes -> "no cell overlap")
seu$cell_type <- unname(cluster_label[as.character(seu@meta.data[[clust_res]])])
print(table(seu$cell_type, useNA = "ifany"))



### 4. Plots ----------------------------------------------------------------
# 4a. Annotated UMAP
png(file.path(outdir, paste0("UMAP_markers_annotated_", s, ".png")),
    width = 1700, height = 1400, res = 200, bg = "white")
print(DimPlot(seu, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE) +
          ggtitle(paste0(s, " — marker-based cell types")))
dev.off()

# 4b. UMAP colored by cluster (for cross-check)
png(file.path(outdir, paste0("UMAP_clusters_", s, ".png")),
    width = 1600, height = 1400, res = 200, bg = "white")
print(DimPlot(seu, reduction = "umap", group.by = clust_res, label = TRUE) +
          ggtitle(paste0(s, " — clusters (", clust_res, ")")))
dev.off()

# 4c. DotPlot: one representative marker per cell type across clusters
rep_markers <- unique(unlist(lapply(marker_list, `[`, 1)))
png(file.path(outdir, paste0("DotPlot_markers_", s, ".png")),
    width = 2000, height = 1400, res = 200, bg = "white")
print(DotPlot(seu, features = rep_markers, group.by = clust_res) +
          RotatedAxis() + ggtitle(paste0(s, " — marker expression by cluster")))
dev.off()

# 4d. Heatmap of mean module score (cluster x cell type) -> shows the assignment
png(file.path(outdir, paste0("heatmap_moduleScores_", s, ".png")),
    width = 1800, height = 1400, res = 200, bg = "white")
pheatmap(t(score_mat), color = viridis(100),
         cluster_rows = FALSE, cluster_cols = TRUE,
         main = paste0(s, " — mean module score (cell type x cluster)"))
dev.off()

# 4e. FeaturePlots for a few canonical markers (visual sanity check)
key_feats <- intersect(c("Myo7a", "Sox2", "Spp1", "Anxa4", "Ptprc", "Pecam1"),
                       rownames(seu))
png(file.path(outdir, paste0("FeaturePlot_key_markers_", s, ".png")),
    width = 2000, height = 1400, res = 200, bg = "white")
print(FeaturePlot(seu, features = key_feats, reduction = "umap", order = TRUE))
dev.off()



### 5. Cluster DE genes (to manually verify / refine the labels) -------------
Idents(seu) <- clust_res
cluster_markers <- FindAllMarkers(seu, only.pos = TRUE,
                                    min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- cluster_markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 15) %>%
    ungroup()
readr::write_tsv(top_markers, file.path(outdir, paste0(s, "_cluster_top_markers.tsv")))



### 6. Save annotated object ------------------------------------------------
saveRDS(seu, file.path(outdir, paste0(s, "_seurat_annotated_markers.rds")))

### END
