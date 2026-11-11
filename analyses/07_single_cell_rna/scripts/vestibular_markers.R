# Shared marker panels + cluster-annotation helper for the vestibular maculae.
# Sourced by annot5_integrate_tissue.R and annot6_tissue_DE.R.
suppressMessages({ library(Seurat); library(dplyr) })

# canonical macular marker panels (same as annot2 / annot4)
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

# regional / positional-identity candidates that may distinguish saccule vs utricle
# (developmental patterning genes; see literature — Otx1/Otx2, Emx2, Lfng, Tbx1, Bmp4)
regional_markers <- c("Otx1", "Otx2", "Emx2", "Lfng", "Tbx1", "Bmp4", "Gbx2", "Gata3")

# representative single marker per cell type (for compact DotPlots)
rep_markers <- vapply(marker_list, `[`, character(1), 1)

# assign each cluster to its highest mean-module-score cell type; scored on RNA assay
annotate_by_markers <- function(seu, clust_col, assay = "RNA") {
    DefaultAssay(seu) <- assay
    ml <- lapply(marker_list, intersect, y = rownames(seu))
    ml <- ml[lengths(ml) > 0]
    seu <- AddModuleScore(seu, features = ml, name = "modscore", seed = 42)
    score_cols <- paste0("modscore", seq_along(ml))
    seu@meta.data <- seu@meta.data %>% rename_with(~ names(ml), all_of(score_cols))
    ms <- seu@meta.data %>%
        group_by(cluster = .data[[clust_col]]) %>%
        summarise(across(all_of(names(ml)), mean), .groups = "drop")
    mat <- as.matrix(ms[, names(ml)]); rownames(mat) <- ms$cluster
    lab <- setNames(colnames(mat)[max.col(mat, ties.method = "first")], rownames(mat))
    seu$cell_type <- unname(lab[as.character(seu@meta.data[[clust_col]])])
    list(seu = seu, labels = lab, score_mat = mat)
}
