# subcluster the hair-cell lineage of FTL_1147 to force a Type I / Type II split
# input : results/cell_annotation/FTL_1147_seurat_annotated_markers.rds  (from annot2.R)
# setwd("analyses/07_single_cell_rna")   # from the project root

rm(list = ls())

library(Seurat)
library(tidyverse)
library(pheatmap)
library(viridis)

s         <- "FTL_1147"
outdir    <- "results/cell_annotation"
hc_res    <- "hc_snn_res.0.4"   # subclustering resolution used for the I/II call
hc_margin <- 0.10               # min gap (winner - loser module score) to call a subtype;
                                # below this the subcluster stays "unresolved I/II"

seu <- readRDS(file.path(outdir, paste0(s, "_seurat_annotated_markers.rds")))



### 1. Subset the hair-cell lineage -----------------------------------------
# includes the generic "Hair_cell" clusters + the already-split "HC_type_II"
hc <- subset(seu, subset = cell_type %in% c("Hair_cell", "HC_type_II"))
cat("hair-cell lineage cells:", ncol(hc), "\n")

# re-embed on this subset only (fresh HVGs/PCA/graph so subtle I-vs-II
# structure is not dominated by the whole-tissue variation)
DefaultAssay(hc) <- "RNA"
hc <- FindVariableFeatures(hc, nfeatures = 2000)
hc <- ScaleData(hc)
hc <- RunPCA(hc, features = VariableFeatures(hc), verbose = FALSE)
hc <- FindNeighbors(hc, dims = 1:15, graph.name = c("hc_nn", "hc_snn"))
hc <- FindClusters(hc, resolution = 0.4, graph.name = "hc_snn")   # -> hc_snn_res.0.4
hc <- RunUMAP(hc, dims = 1:15)
Idents(hc) <- hc_res



### 2. Type I / Type II marker panels ---------------------------------------
# NOTE: Adgrv1 dropped from Type I - it is a pan-hair-cell stereocilia (Usher)
# gene expressed in BOTH subtypes, so it does not discriminate. Kept the
# subtype-specific markers only.
hc_markers <- list(
    HC_type_I  = c("Spp1", "Ocm", "Cabp2"),             # Type I: Spp1/oncomodulin
    HC_type_II = c("Anxa4", "Mapt", "Calb2")            # Type II: Anxa4/Mapt/calretinin
)
hc_markers <- lapply(hc_markers, intersect, y = rownames(hc))

hc <- AddModuleScore(hc, features = hc_markers, name = "hcmod", seed = 42)
score_cols  <- paste0("hcmod", seq_along(hc_markers))
# distinct names so we don't clash with the inherited module-score columns
# (annot2.R already left columns literally named "HC_type_I"/"HC_type_II")
score_names <- paste0("hcscore_", names(hc_markers))
hc@meta.data <- hc@meta.data %>%
    rename_with(~ score_names, all_of(score_cols))



### 3. Assign each subcluster to Type I or Type II (argmax) ------------------
mean_scores <- hc@meta.data %>%
    group_by(cluster = .data[[hc_res]]) %>%
    summarise(across(all_of(score_names), mean), .groups = "drop")

score_mat <- as.matrix(mean_scores[, score_names])
rownames(score_mat) <- mean_scores$cluster
colnames(score_mat) <- names(hc_markers)   # map score cols back to subtype labels

# confident-only: label a subcluster only if its winning module score is
# positive AND beats the other subtype by >= hc_margin; otherwise "unresolved".
# 3' scRNA captures Type I markers (Ocm/Spp1) poorly, so most subclusters do
# not separate cleanly and are honestly left unresolved.
best       <- colnames(score_mat)[max.col(score_mat, ties.method = "first")]
winner_val <- apply(score_mat, 1, max)
gap        <- apply(score_mat, 1, function(x) max(x) - min(x))  # 2 cols: winner - loser
confident  <- winner_val > 0 & gap >= hc_margin

cluster_label <- setNames(ifelse(confident, best, "Hair_cell (unresolved I/II)"),
                          rownames(score_mat))

# transparency: decision table (scores, gap, call) per subcluster
decision <- data.frame(subcluster = rownames(score_mat),
                       round(score_mat, 3),
                       gap = round(gap, 3),
                       confident = confident,
                       label = cluster_label,
                       row.names = NULL)
print(decision)

hc$hc_subtype <- unname(cluster_label[as.character(hc@meta.data[[hc_res]])])
print(table(hc$hc_subtype))



### 4. Write the refined labels back into the full object -------------------
seu$cell_type_refined <- as.character(seu$cell_type)
hc_sub  <- setNames(hc$hc_subtype, colnames(hc))            # barcode -> subtype
is_hc   <- colnames(seu) %in% names(hc_sub)                 # hair-cell barcodes in seu
seu$cell_type_refined[is_hc] <- hc_sub[colnames(seu)[is_hc]]
print(table(seu$cell_type_refined))



### 5. Plots ----------------------------------------------------------------
# 5a. subset UMAP colored by Type I/II
png(file.path(outdir, paste0("UMAP_HC_subtype_", s, ".png")),
    width = 1600, height = 1400, res = 200, bg = "white")
print(DimPlot(hc, reduction = "umap", group.by = "hc_subtype", label = TRUE) +
          ggtitle(paste0(s, " — hair-cell Type I / II")))
dev.off()

# 5b. subcluster UMAP (raw hair-cell subclusters)
png(file.path(outdir, paste0("UMAP_HC_subclusters_", s, ".png")),
    width = 1600, height = 1400, res = 200, bg = "white")
print(DimPlot(hc, reduction = "umap", group.by = hc_res, label = TRUE) +
          ggtitle(paste0(s, " — hair-cell subclusters")))
dev.off()

# 5c. I vs II marker DotPlot across subclusters
png(file.path(outdir, paste0("DotPlot_HC_IvsII_", s, ".png")),
    width = 1500, height = 1200, res = 200, bg = "white")
print(DotPlot(hc, features = unlist(hc_markers, use.names = FALSE),
              group.by = hc_res) + RotatedAxis() +
          ggtitle(paste0(s, " — Type I vs II markers by subcluster")))
dev.off()

# 5d. FeaturePlots of the canonical I (Spp1/Ocm) vs II (Anxa4/Mapt) markers
key <- intersect(c("Spp1", "Ocm", "Anxa4", "Mapt"), rownames(hc))
png(file.path(outdir, paste0("FeaturePlot_HC_IvsII_", s, ".png")),
    width = 1600, height = 1400, res = 200, bg = "white")
print(FeaturePlot(hc, features = key, reduction = "umap", order = TRUE))
dev.off()

# 5e. mean module score heatmap (subcluster x subtype)
png(file.path(outdir, paste0("heatmap_HC_moduleScores_", s, ".png")),
    width = 1400, height = 1000, res = 200, bg = "white")
pheatmap(t(score_mat), color = viridis(100),
         cluster_rows = FALSE, cluster_cols = TRUE,
         main = paste0(s, " — HC module score (subtype x subcluster)"))
dev.off()

# 5f. full-tissue UMAP with the refined labels (Type I/II now explicit)
png(file.path(outdir, paste0("UMAP_markers_refined_", s, ".png")),
    width = 1700, height = 1400, res = 200, bg = "white")
print(DimPlot(seu, reduction = "umap", group.by = "cell_type_refined",
              label = TRUE, repel = TRUE) +
          ggtitle(paste0(s, " — cell types (hair cells split I/II)")))
dev.off()



### 6. Save -----------------------------------------------------------------
saveRDS(hc,  file.path(outdir, paste0(s, "_hair_cells_subclustered.rds")))
saveRDS(seu, file.path(outdir, paste0(s, "_seurat_annotated_markers_refined.rds")))

### END
