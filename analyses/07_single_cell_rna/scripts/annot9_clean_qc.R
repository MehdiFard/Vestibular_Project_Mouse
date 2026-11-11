# annot9_clean_qc.R  (seurat_env + SoupX + scDblFinder)
# Per-sample cleaning:  ambient-RNA correction (SoupX) -> QC filter -> doublet
# removal (scDblFinder).  Saves one cleaned Seurat object per sample plus a
# filtering-stats table.  All outputs go to results/cell_annotation_clean/.
# setwd("analyses/07_single_cell_rna")   # from the project root

rm(list = ls())
suppressMessages({
    library(Seurat); library(SoupX); library(scDblFinder)
    library(SingleCellExperiment); library(Matrix)
})
set.seed(42)

outdir <- "results/cell_annotation_clean"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

samples <- c(FTL_1146 = "saccule", FTL_1147 = "utricle",
             FTL_1148 = "saccule", FTL_1149 = "utricle")
reps    <- c(FTL_1146 = "rep1", FTL_1147 = "rep1", FTL_1148 = "rep2", FTL_1149 = "rep2")

# same QC thresholds as the earlier (uncleaned) runs, for a fair comparison
MIN_COUNT <- 800; MIN_FEAT <- 500; MAX_MITO <- 5

stats <- list()
for (s in names(samples)) {
    message("==== cleaning ", s, " ====")
    base <- file.path("intermediate/cell_ranger_output", s, "outs")
    toc  <- Read10X_h5(file.path(base, "filtered_feature_bc_matrix.h5"))  # cells
    tod  <- Read10X_h5(file.path(base, "raw_feature_bc_matrix.h5"))       # all droplets
    tod  <- tod[rownames(toc), ]                                          # align genes

    ## --- quick clustering (needed by SoupX autoEstCont) ---
    q <- CreateSeuratObject(toc)
    q <- NormalizeData(q, verbose = FALSE)
    q <- FindVariableFeatures(q, verbose = FALSE)
    q <- ScaleData(q, verbose = FALSE)
    q <- RunPCA(q, npcs = 30, verbose = FALSE)
    q <- FindNeighbors(q, dims = 1:20, verbose = FALSE)
    q <- FindClusters(q, resolution = 0.5, verbose = FALSE)
    clusters <- setNames(as.character(Idents(q)), colnames(q))

    ## --- SoupX ambient-RNA correction ---
    sc  <- SoupChannel(tod, toc)
    sc  <- setClusters(sc, clusters)
    sc  <- tryCatch(autoEstCont(sc, doPlot = FALSE),
                    error = function(e) {
                        message("  autoEstCont failed (", e$message, ") -> rho=0.1")
                        setContaminationFraction(sc, 0.1)
                    })
    rho <- mean(sc$metaData$rho)
    adj <- adjustCounts(sc, roundToInt = TRUE)                            # corrected counts

    ## --- build object + QC filter ---
    o <- CreateSeuratObject(adj, project = s, min.cells = 3, min.features = 200)
    o$mitoPercent <- PercentageFeatureSet(o, pattern = "^mt-")
    n_ambient <- ncol(o)
    o <- subset(o, subset = nCount_RNA > MIN_COUNT & nFeature_RNA > MIN_FEAT &
                            mitoPercent < MAX_MITO)
    n_qc <- ncol(o)

    ## --- doublet removal (scDblFinder) ---
    sce <- SingleCellExperiment(list(counts = GetAssayData(o, layer = "counts")))
    sce <- scDblFinder(sce)
    o$scDblFinder.class <- sce$scDblFinder.class
    o$scDblFinder.score <- sce$scDblFinder.score
    n_doublet <- sum(o$scDblFinder.class == "doublet")
    o <- subset(o, subset = scDblFinder.class == "singlet")

    o$sample <- s; o$tissue <- samples[[s]]; o$replicate <- reps[[s]]
    saveRDS(o, file.path(outdir, paste0(s, "_clean.rds")))

    stats[[s]] <- data.frame(
        sample = s, tissue = samples[[s]], replicate = reps[[s]],
        n_cellranger    = ncol(toc),
        soupx_rho       = round(rho, 4),
        n_after_ambient = n_ambient,
        n_after_qc      = n_qc,
        n_doublets      = n_doublet,
        n_final         = ncol(o),
        pct_removed     = round(100 * (ncol(toc) - ncol(o)) / ncol(toc), 1))
    message("  ", s, ": rho=", round(rho,3), " cellranger=", ncol(toc),
            " -> final=", ncol(o), " (", n_doublet, " doublets)")
}

st <- do.call(rbind, stats)
write.csv(st, file.path(outdir, "filtering_stats.csv"), row.names = FALSE)
cat("\n==== FILTERING SUMMARY ====\n"); print(st)
message("==== annot9 cleaning done ====")
### END
