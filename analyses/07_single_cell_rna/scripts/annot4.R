# annot4.R — integrate all 4 samples, marker-annotate, check cross-sample consistency
# samples : FTL_1146 saccule_1, FTL_1147 utricle_1, FTL_1148 saccule_2, FTL_1149 utricle_2
# goal    : one batch-corrected cell-type map + validation that the marker calls are
#           reproducible across the saccule/utricle replicates (not batch artifacts).
# setwd("analyses/07_single_cell_rna")   # from the project root

rm(list = ls())

library(Seurat)
library(tidyverse)
library(patchwork)
library(pheatmap)
library(viridis)

outdir <- "results/cell_annotation"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

samples <- c(FTL_1146 = "saccule", FTL_1147 = "utricle",
             FTL_1148 = "saccule", FTL_1149 = "utricle")
clust_res <- "integrated_snn_res.0.3"



### 1. Read + per-sample QC -------------------------------------------------
obj.list <- list()
for (s in names(samples)) {
    m <- Read10X_h5(file.path("intermediate/cell_ranger_output", s,
                              "outs/filtered_feature_bc_matrix.h5"),
                    use.names = TRUE, unique.features = TRUE)
    o <- CreateSeuratObject(m, project = s, min.cells = 3, min.features = 200)
    o$mitoPercent <- PercentageFeatureSet(o, pattern = "^mt-")   # mouse: mt-
    o <- subset(o, subset = nCount_RNA > 800 & nFeature_RNA > 500 & mitoPercent < 5)
    o$sample <- s
    o$tissue <- samples[[s]]
    o$replicate <- ifelse(s %in% c("FTL_1146", "FTL_1147"), "rep1", "rep2")
    o <- NormalizeData(o, verbose = FALSE)
    o <- FindVariableFeatures(o, nfeatures = 2000, verbose = FALSE)
    obj.list[[s]] <- o
}



### 2. CCA integration (batch correction) -----------------------------------
features <- SelectIntegrationFeatures(object.list = obj.list)
anchors  <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)
intg     <- IntegrateData(anchorset = anchors)

DefaultAssay(intg) <- "integrated"
intg <- ScaleData(intg, verbose = FALSE)
intg <- RunPCA(intg, npcs = 30, verbose = FALSE)
intg <- FindNeighbors(intg, dims = 1:30)
intg <- FindClusters(intg, resolution = c(0.1, 0.3, 0.5))
intg <- RunUMAP(intg, dims = 1:30)
Idents(intg) <- clust_res



### 3. Marker-based annotation (same panels as annot2), scored on RNA --------
DefaultAssay(intg) <- "RNA"
intg <- tryCatch(JoinLayers(intg), error = function(e) intg)  # Seurat v5 layers
intg <- NormalizeData(intg, verbose = FALSE)

marker_list <- list(
    Hair_cell          = c("Myo7a", "Pou4f3", "Gfi1", "Myo6", "Ccer2", "Tmc1"),
    HC_type_I          = c("Spp1", "Ocm", "Adgrv1", "Cabp2"),
    HC_type_II         = c("Mapt", "Calb2", "Anxa4"),
    Supporting_cell    = c("Sox2", "Gjb2", "Gjb6", "Hey1", "Hey2", "Otog", "Otogl", "S100b"),
    Transitional_epith = c("Krt8", "Krt18", "Cldn3", "Cldn6"),
    Roof_dark_cell     = c("Kcnq1", "Kcne1", "Kcnj10", "Atp1b2", "Slc12a2"),
    Glia_Schwann       = c("Plp1", "Mpz", "Sox10", "Mbp"),
    Neuron             = c("Tubb3", "Snap25", "Nefl", "Nefm", "Pvalb", "Stmn2"),
    Immune_macrophage  = c("Ptprc", "C1qa", "C1qb", "Aif1", "Cd68", "Lyz2"),
    Endothelial        = c("Pecam1", "Cldn5", "Flt1", "Cdh5"),
    Melanocyte         = c("Mlana", "Dct", "Pmel", "Tyrp1"),
    Fibrocyte_mesench  = c("Col1a1", "Col1a2", "Dcn", "Lum", "Pdgfra"),
    Pericyte_SMC       = c("Acta2", "Pdgfrb", "Rgs5", "Myh11"),
    Erythrocyte        = c("Hba-a1", "Hbb-bs")
)
marker_list <- lapply(marker_list, intersect, y = rownames(intg))
marker_list <- marker_list[lengths(marker_list) > 0]

intg <- AddModuleScore(intg, features = marker_list, name = "modscore", seed = 42)
score_cols <- paste0("modscore", seq_along(marker_list))
intg@meta.data <- intg@meta.data %>% rename_with(~ names(marker_list), all_of(score_cols))

mean_scores <- intg@meta.data %>%
    group_by(cluster = .data[[clust_res]]) %>%
    summarise(across(all_of(names(marker_list)), mean), .groups = "drop")
score_mat <- as.matrix(mean_scores[, names(marker_list)]); rownames(score_mat) <- mean_scores$cluster

cluster_label <- setNames(colnames(score_mat)[max.col(score_mat, ties.method = "first")],
                          rownames(score_mat))
print(cluster_label)
intg$cell_type <- unname(cluster_label[as.character(intg@meta.data[[clust_res]])])
print(table(intg$cell_type))



### 4. Cross-sample / cross-tissue consistency ------------------------------
comp <- intg@meta.data %>%
    count(tissue, sample, replicate, cell_type) %>%
    group_by(sample) %>% mutate(prop = n / sum(n)) %>% ungroup()
write_csv(comp, file.path(outdir, "integrated_celltype_composition_by_sample.csv"))

# per-cluster sample mixing: a well-integrated cluster draws from multiple samples
mix <- table(cluster = intg@meta.data[[clust_res]], sample = intg$sample)
write.csv(as.data.frame.matrix(mix), file.path(outdir, "integrated_cluster_sample_mixing.csv"))

# is a cell type seen in BOTH replicates of each tissue? (reproducibility)
repro <- intg@meta.data %>%
    distinct(tissue, replicate, cell_type) %>%
    count(tissue, cell_type, name = "n_reps_present")   # 2 = in both replicates
print(as.data.frame(repro %>% arrange(tissue, desc(n_reps_present))))
write_csv(repro, file.path(outdir, "integrated_celltype_reproducibility.csv"))



### 5. Plots ----------------------------------------------------------------
p_ct <- DimPlot(intg, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE) +
            ggtitle("integrated — marker cell types")
p_ti <- DimPlot(intg, reduction = "umap", group.by = "tissue")   + ggtitle("tissue")
p_sa <- DimPlot(intg, reduction = "umap", group.by = "sample")   + ggtitle("sample (batch)")

png(file.path(outdir, "UMAP_integrated_celltype_FTL_all.png"),
    width = 1800, height = 1500, res = 200, bg = "white"); print(p_ct); dev.off()
png(file.path(outdir, "UMAP_integrated_tissue_sample_FTL_all.png"),
    width = 2400, height = 1100, res = 200, bg = "white"); print(p_ti + p_sa); dev.off()

# stacked bar: cell-type composition per sample (visual consistency check)
p_bar <- ggplot(comp, aes(x = sample, y = prop, fill = cell_type)) +
    geom_col() + facet_grid(~ tissue, scales = "free_x", space = "free_x") +
    labs(y = "proportion of cells", x = NULL,
         title = "cell-type composition across samples (utricle vs saccule)") +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(outdir, "barplot_celltype_composition_FTL_all.png"),
       p_bar, width = 9, height = 6, dpi = 200, bg = "white")



### 6. Save -----------------------------------------------------------------
saveRDS(intg, file.path(outdir, "integrated_all_samples_annotated.rds"))

### END
