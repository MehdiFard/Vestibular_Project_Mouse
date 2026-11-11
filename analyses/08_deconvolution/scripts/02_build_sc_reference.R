#!/usr/bin/env Rscript
# 02_build_sc_reference.R  (run in seurat_env)
# ---------------------------------------------------------------------------
# Build the single-cell reference for deconvolution from the cleaned, jointly
# integrated snRNA object. We keep the 6 abundant cell types and collapse the
# hair-cell annotation to a single "Hair_cell" (the HC_type_I/II split is an
# annotation-threshold artifact, not a real organ-specific type).
#
# Outputs (results/inputs/):
#   sc_counts.rds       raw RNA counts, genes(symbol) x cells  (dgCMatrix)
#   sc_meta.rds         data.frame: cell_type, SubjectName (sample), tissue
#   sc_markers.rds/.csv per-cell-type marker genes (for ssGSEA enrichment)
#   sc_composition.csv  observed proportions per sample (validation reference)
#
# Heavy (loads a ~1.4 GB object) -> run on a compute node, not the login node.
# ---------------------------------------------------------------------------
suppressMessages({ library(Seurat); library(Matrix) })

IN   <- "../07_single_cell_rna/results/cell_annotation_clean/all_integrated_clean.rds"
OUT  <- "results/inputs"
KEEP <- c("Hair_cell", "Supporting_cell", "Transitional_epith",
          "Glia_Schwann", "Fibrocyte_mesench", "Roof_dark_cell")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cat("[1] loading", IN, "\n")
o <- readRDS(IN)

# subset to the 6 abundant types (Hair_cell already collapsed in the joint object)
o <- subset(o, subset = cell_type %in% KEEP)
o$cell_type <- factor(as.character(o$cell_type), levels = KEEP)
cat("    cells kept:", ncol(o), "\n"); print(table(o$cell_type))

# raw counts from the RNA assay (join split v5 layers first)
DefaultAssay(o) <- "RNA"
o <- JoinLayers(o)
counts <- GetAssayData(o, assay = "RNA", layer = "counts")
cat("[2] counts:", nrow(counts), "genes x", ncol(counts), "cells\n")

meta <- data.frame(
    cell_type   = o$cell_type,
    SubjectName = as.character(o$sample),   # FTL_1146/1147/1148/1149 = 4 donors
    tissue      = as.character(o$tissue),
    row.names   = colnames(o), stringsAsFactors = FALSE)
cat("    donors (SubjectName):\n"); print(table(meta$SubjectName, meta$cell_type))

saveRDS(counts, file.path(OUT, "sc_counts.rds"))
saveRDS(meta,   file.path(OUT, "sc_meta.rds"))

# observed sc composition per sample (validation reference for the deconvolution)
comp <- as.data.frame(prop.table(table(meta$SubjectName, meta$cell_type), margin = 1))
colnames(comp) <- c("SubjectName", "cell_type", "prop")
write.csv(comp, file.path(OUT, "sc_composition.csv"), row.names = FALSE)

# ---- marker genes per cell type (for the ssGSEA enrichment approach) --------
# downsample to <=1500 cells/type for a fast, balanced marker scan
cat("[3] marker genes (FindAllMarkers, downsampled)\n")
set.seed(1)
Idents(o) <- "cell_type"
o.sub <- subset(o, downsample = 1500)
o.sub <- NormalizeData(o.sub, verbose = FALSE)
mk <- FindAllMarkers(o.sub, assay = "RNA", only.pos = TRUE,
                     min.pct = 0.25, logfc.threshold = 0.5, verbose = FALSE)
# top 50 specific markers per type by log2FC
markers <- lapply(KEEP, function(ct) {
    m <- mk[mk$cluster == ct, ]
    head(m$gene[order(-m$avg_log2FC)], 50)
})
names(markers) <- KEEP
saveRDS(markers, file.path(OUT, "sc_markers.rds"))
write.csv(mk, file.path(OUT, "sc_markers.csv"), row.names = FALSE)
cat("    markers/type:", sapply(markers, length), "\n")

cat("DONE\n")
