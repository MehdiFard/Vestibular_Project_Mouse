# annot10_integrate_clean.R  (seurat_env)
# Integrate the CLEANED per-sample objects (annot9) and annotate:
#   - saccule (2 reps)   - utricle (2 reps)   - all 4 samples
# Produces UMAPs + cell-type composition for each, into results/cell_annotation_clean/.
# setwd("analyses/07_single_cell_rna")   # from the project root

rm(list = ls())
suppressMessages({ library(Seurat); library(tidyverse); library(patchwork) })
source("scripts/vestibular_markers.R")

outdir <- "results/cell_annotation_clean"
samples <- c("FTL_1146", "FTL_1147", "FTL_1148", "FTL_1149")
clean <- lapply(samples, function(s) readRDS(file.path(outdir, paste0(s, "_clean.rds"))))
names(clean) <- samples

clust_res <- "integrated_snn_res.0.3"

integrate_annotate <- function(objs, name) {
    message("==== integrating ", name, " (", length(objs), " samples) ====")
    for (i in seq_along(objs)) {
        objs[[i]] <- NormalizeData(objs[[i]], verbose = FALSE)
        objs[[i]] <- FindVariableFeatures(objs[[i]], nfeatures = 2000, verbose = FALSE)
    }
    feats <- SelectIntegrationFeatures(objs)
    anch  <- FindIntegrationAnchors(objs, anchor.features = feats)
    intg  <- IntegrateData(anch)
    DefaultAssay(intg) <- "integrated"
    intg <- ScaleData(intg, verbose = FALSE)
    intg <- RunPCA(intg, npcs = 30, verbose = FALSE)
    intg <- FindNeighbors(intg, dims = 1:30)
    intg <- FindClusters(intg, resolution = c(0.1, 0.3, 0.5))
    intg <- RunUMAP(intg, dims = 1:30)

    DefaultAssay(intg) <- "RNA"
    intg <- tryCatch(JoinLayers(intg), error = function(e) intg)
    intg <- NormalizeData(intg, verbose = FALSE)
    ann  <- annotate_by_markers(intg, clust_res); intg <- ann$seu
    cat("----", name, "cell types ----\n"); print(sort(table(intg$cell_type), decreasing = TRUE))

    # UMAP by cell type
    png(file.path(outdir, paste0(name, "_clean_UMAP_celltype.png")),
        width = 1800, height = 1500, res = 200, bg = "white")
    print(DimPlot(intg, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE) +
              ggtitle(paste0(name, " (cleaned) — cell types")))
    dev.off()
    # UMAP by tissue+sample (all) or replicate (per organ)
    png(file.path(outdir, paste0(name, "_clean_UMAP_batch.png")),
        width = 2400, height = 1100, res = 200, bg = "white")
    if (length(unique(intg$tissue)) > 1)
        print((DimPlot(intg, reduction="umap", group.by="tissue") + ggtitle("tissue")) +
              (DimPlot(intg, reduction="umap", group.by="sample") + ggtitle("sample")))
    else
        print(DimPlot(intg, reduction="umap", group.by="replicate") + ggtitle("replicate (batch)"))
    dev.off()

    # composition
    comp <- intg@meta.data %>% count(tissue, sample, replicate, cell_type) %>%
        group_by(sample) %>% mutate(prop = n / sum(n)) %>% ungroup()
    write_csv(comp, file.path(outdir, paste0(name, "_clean_composition.csv")))
    p_bar <- ggplot(comp, aes(sample, prop, fill = cell_type)) + geom_col() +
        facet_grid(~ tissue, scales = "free_x", space = "free_x") +
        labs(y = "proportion of cells", x = NULL, title = paste0(name, " (cleaned) — composition")) +
        theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(outdir, paste0(name, "_clean_composition_bar.png")),
           p_bar, width = 8, height = 6, dpi = 200, bg = "white")

    saveRDS(intg, file.path(outdir, paste0(name, "_integrated_clean.rds")))
    intg
}

invisible(integrate_annotate(clean[c("FTL_1146", "FTL_1148")], "saccule"))
invisible(integrate_annotate(clean[c("FTL_1147", "FTL_1149")], "utricle"))
invisible(integrate_annotate(clean, "all"))

message("==== annot10 integration/annotation done ====")
### END
