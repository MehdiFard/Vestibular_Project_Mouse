# annot3.R — label transfer from a published reference + reconcile with markers
# reference : Burns et al. 2015, GSE71982 (P1 mouse utricle sensory epithelium)
# query     : FTL_1147, already marker-annotated by annot2.R
# approach  : SingleR (primary) + Seurat TransferData (secondary), then compare
#             the transferred labels against the annot2.R marker labels.
# setwd("analyses/07_single_cell_rna")   # from the project root

rm(list = ls())

library(SingleR)
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(pheatmap)
library(viridis)

outdir <- "results/cell_annotation"
refdir <- "src/annot_ref/burns2015"
dir.create(refdir, showWarnings = FALSE, recursive = TRUE)

# IMPORTANT CAVEAT baked into the interpretation below:
# Burns 2015 is P1 (neonatal) SMART-seq of the SENSORY EPITHELIUM only -> its
# labels are Hair cell (HC), Supporting cell (SC), Transitional epithelial (TEC).
# It contains NO stromal / immune / endothelial / glial cells, and its HC (i/ii/
# iii-iv) groups are DEVELOPMENTAL stages, not vestibular Type I/II. So transfer
# is only trustworthy for the epithelial compartment; non-epithelial query cells
# will be forced into HC/SC/TEC and should show LOW confidence (we flag those).



### 1. Download Burns 2015 reference (skip if present) -----------------------
base <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE71nnn/GSE71982/suppl"
for (f in c("GSE71982_P1_Utricle_PhenoData.txt", "GSE71982_RSEM_Counts_Matrix.txt")) {
    if (!file.exists(file.path(refdir, f))) {
        download.file(file.path(base, paste0(f, ".gz")),
                      destfile = file.path(refdir, paste0(f, ".gz")), mode = "wb")
        system(paste("gunzip -f", file.path(refdir, paste0(f, ".gz"))))
    }
}



### 2. Build the reference SingleCellExperiment -----------------------------
ref_counts <- read.delim(file.path(refdir, "GSE71982_RSEM_Counts_Matrix.txt"),
                         row.names = 1, check.names = FALSE)
rownames(ref_counts) <- gsub('"', "", rownames(ref_counts))    # strip stray quotes
colnames(ref_counts) <- gsub('"', "", colnames(ref_counts))

pheno <- read.delim(file.path(refdir, "GSE71982_P1_Utricle_PhenoData.txt"),
                    check.names = FALSE)

# keep the utricle cells that have a label AND are in the matrix
keep <- intersect(pheno$Short_Name, colnames(ref_counts))
pheno <- pheno[match(keep, pheno$Short_Name), ]
ref_counts <- as.matrix(ref_counts[, keep])
cat("reference cells:", length(keep), "\n")

# collapse the fine GroupID (HC (i)/(ii)/(iii-iv), SC (i)/(ii), TEC) to major types
coarse <- dplyr::case_when(
    grepl("^HC",  pheno$GroupID) ~ "HC",
    grepl("^SC",  pheno$GroupID) ~ "SC",
    grepl("^TEC", pheno$GroupID) ~ "TEC",
    TRUE ~ pheno$GroupID
)
print(table(coarse))

ref <- SingleCellExperiment(assays = list(counts = ref_counts),
                            colData = DataFrame(celltype = coarse,
                                                celltype_fine = pheno$GroupID,
                                                row.names = keep))
# log-normalize (RSEM expected counts -> CPM/10k -> log1p); SingleR wants logcounts
sf <- colSums(assay(ref, "counts")); sf[sf == 0] <- 1
assay(ref, "logcounts") <- log1p(t(t(assay(ref, "counts")) / sf) * 1e4)



### 3. Load the marker-annotated query (from annot2.R) ----------------------
seu <- readRDS(file.path(outdir, "FTL_1147_seurat_annotated_markers_refined.rds"))

# map the marker labels onto the reference's epithelial vocabulary for comparison
marker_to_coarse <- c(Hair_cell = "HC", HC_type_II = "HC",
                      "Hair_cell (unresolved I/II)" = "HC",
                      Supporting_cell = "SC", Transitional_epith = "TEC")
seu$marker_coarse <- unname(marker_to_coarse[seu$cell_type_refined])
seu$marker_coarse[is.na(seu$marker_coarse)] <- "other (non-epithelial)"
# is this query cell in a compartment the reference can actually explain?
seu$in_ref_compartment <- seu$marker_coarse %in% c("HC", "SC", "TEC")



### 4. SingleR label transfer (primary; robust across platforms) ------------
pred <- SingleR(test      = GetAssayData(seu, layer = "data"),
                ref       = ref,
                labels    = ref$celltype,
                de.method = "classic")               # "classic" -> no scran needed
seu$burns_singleR        <- pred$labels[match(colnames(seu), rownames(pred))]
seu$burns_singleR_pruned <- pred$pruned.labels[match(colnames(seu), rownames(pred))]
cat("\n== SingleR labels (pruned; NA = low confidence) ==\n")
print(table(seu$burns_singleR_pruned, useNA = "ifany"))



### 5. Seurat TransferData (secondary; wrapped - small cross-platform ref) ---
seu$burns_seurat    <- NA_character_
seu$burns_seurat_sc <- NA_real_
try({
    ref.seu <- CreateSeuratObject(counts = assay(ref, "counts"))
    ref.seu$celltype <- ref$celltype
    ref.seu <- NormalizeData(ref.seu, verbose = FALSE)
    ref.seu <- FindVariableFeatures(ref.seu, nfeatures = 2000, verbose = FALSE)
    ref.seu <- ScaleData(ref.seu, verbose = FALSE)
    ref.seu <- RunPCA(ref.seu, npcs = 30, verbose = FALSE)

    anchors <- FindTransferAnchors(reference = ref.seu, query = seu,
                                   dims = 1:30, reference.reduction = "pca",
                                   k.filter = NA)      # ref has only ~158 cells
    tp <- TransferData(anchorset = anchors, refdata = ref.seu$celltype,
                       dims = 1:30, k.weight = 30)
    seu$burns_seurat    <- tp$predicted.id
    seu$burns_seurat_sc <- tp$prediction.score.max
    cat("\n== Seurat TransferData labels ==\n"); print(table(seu$burns_seurat))
}, silent = FALSE)



### 6. Reconcile: markers vs transferred labels -----------------------------
# focus on the epithelial compartment, where the reference is valid
ep <- seu@meta.data %>% filter(in_ref_compartment)

cm_singleR <- table(marker = ep$marker_coarse, SingleR = ep$burns_singleR)
cat("\n== marker vs SingleR (epithelial cells) ==\n"); print(cm_singleR)
agree_singleR <- mean(ep$marker_coarse == ep$burns_singleR, na.rm = TRUE)
cat(sprintf("epithelial agreement (marker vs SingleR): %.1f%%\n", 100 * agree_singleR))

write.csv(as.data.frame.matrix(cm_singleR),
          file.path(outdir, "FTL_1147_confusion_marker_vs_SingleR.csv"))
if (any(!is.na(seu$burns_seurat))) {
    cm_seurat <- table(marker = ep$marker_coarse, Seurat = ep$burns_seurat)
    write.csv(as.data.frame.matrix(cm_seurat),
              file.path(outdir, "FTL_1147_confusion_marker_vs_Seurat.csv"))
}



### 7. Plots ----------------------------------------------------------------
# 7a. UMAP colored by transferred (SingleR) label
png(file.path(outdir, "UMAP_burns_SingleR_FTL_1147.png"),
    width = 1700, height = 1400, res = 200, bg = "white")
print(DimPlot(seu, reduction = "umap", group.by = "burns_singleR", label = TRUE, repel = TRUE) +
          ggtitle("FTL_1147 — Burns 2015 labels via SingleR"))
dev.off()

# 7b. side-by-side: marker labels vs transferred labels
png(file.path(outdir, "UMAP_marker_vs_burns_FTL_1147.png"),
    width = 2400, height = 1200, res = 200, bg = "white")
print((DimPlot(seu, reduction = "umap", group.by = "cell_type_refined", label = TRUE, repel = TRUE) +
           ggtitle("marker-based (annot2)")) +
      (DimPlot(seu, reduction = "umap", group.by = "burns_singleR", label = TRUE, repel = TRUE) +
           ggtitle("Burns transfer (SingleR)")))
dev.off()

# 7c. confusion heatmap (row-normalized) on epithelial cells
cm_norm <- sweep(cm_singleR, 1, rowSums(cm_singleR), "/")
png(file.path(outdir, "heatmap_confusion_marker_vs_SingleR_FTL_1147.png"),
    width = 1400, height = 1200, res = 200, bg = "white")
pheatmap(cm_norm, color = viridis(100), cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.2f",
         main = "FTL_1147 epithelial: marker (rows) vs SingleR (cols)")
dev.off()

# 7d. SingleR delta / score heatmap (assignment confidence)
png(file.path(outdir, "heatmap_burns_SingleR_scores_FTL_1147.png"),
    width = 1800, height = 1200, res = 200, bg = "white")
print(plotScoreHeatmap(pred))
dev.off()



### 8. Save -----------------------------------------------------------------
saveRDS(seu, file.path(outdir, "FTL_1147_seurat_labeltransfer.rds"))

### END
