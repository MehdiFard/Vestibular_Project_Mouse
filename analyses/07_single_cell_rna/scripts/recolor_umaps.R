# recolor_umaps.R  (seurat_env)
# Redraw the cleaned cell-type UMAPs (saccule / utricle / all) with ONE shared
# palette so colours are consistent across panels (and HC_type_I != Immune).
suppressMessages({ library(Seurat); library(ggplot2) })
source("scripts/celltype_palette.R")
pal <- celltype_pal

jobs <- list(
    saccule = c("saccule_integrated_clean.rds", "saccule_clean_UMAP_celltype.png", "saccule (cleaned) — cell types"),
    utricle = c("utricle_integrated_clean.rds", "utricle_clean_UMAP_celltype.png", "utricle (cleaned) — cell types"),
    all     = c("all_integrated_clean.rds",     "all_clean_UMAP_celltype.png",     "all (cleaned) — cell types"))

dir <- "results/cell_annotation_clean"
for (nm in names(jobs)) {
    o <- readRDS(file.path(dir, jobs[[nm]][1]))
    stopifnot(all(unique(o$cell_type) %in% names(pal)))
    png(file.path(dir, jobs[[nm]][2]), width = 1800, height = 1500, res = 200, bg = "white")
    print(DimPlot(o, reduction = "umap", group.by = "cell_type", cols = pal,
                  label = TRUE, repel = TRUE) + ggtitle(jobs[[nm]][3]))
    dev.off()
    cat("recolored:", nm, "\n"); rm(o); gc(verbose = FALSE)
}
### END
