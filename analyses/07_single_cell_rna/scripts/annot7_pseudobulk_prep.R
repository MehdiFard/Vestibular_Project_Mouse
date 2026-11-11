# annot7_pseudobulk_prep.R  (run with seurat_env)
# Build pseudobulk count matrices from the combined integrated object, restricted
# to protein_coding + lncRNA genes:
#   - organ level   : summed counts per SAMPLE (all cells pooled)
#   - cell level    : summed counts per (cell_type x SAMPLE)
# Saves a plain-R list (base matrices) for annot8 to run PCA + edgeR DE in Env1.
# setwd("analyses/07_single_cell_rna")   # from the project root

rm(list = ls())
suppressMessages({ library(Seurat); library(Matrix) })

outdir <- "results/pseudobulk_DE"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# --- protein_coding + lncRNA gene set (from the 10x GRCm39 GTF) ---
bt <- read.delim("src/annot_ref/gene_biotypes_GRCm39.tsv", header = FALSE,
                 col.names = c("gene", "type"))
keep_genes <- unique(bt$gene[bt$type %in% c("protein_coding", "lncRNA")])
cat("biotype keep-list (pc + lncRNA):", length(keep_genes), "genes\n")

# --- combined object (both tissues, matched cell_type labels) ---
o <- readRDS("results/cell_annotation/integrated_all_samples_annotated.rds")
DefaultAssay(o) <- "RNA"
o <- tryCatch(JoinLayers(o), error = function(e) o)

keep <- intersect(keep_genes, rownames(o))
cat("genes kept in data:", length(keep), "of", nrow(o), "\n")
counts <- GetAssayData(o, assay = "RNA", layer = "counts")[keep, ]

# --- organ-level pseudobulk: sum counts per sample (sparse indicator multiply) ---
fs <- factor(o$sample)
mm <- Matrix::sparse.model.matrix(~ 0 + fs); colnames(mm) <- levels(fs)
pb_sample <- as.matrix(counts %*% mm)                       # genes x samples

# --- cell-level pseudobulk: sum per cell_type x sample ---
cts <- paste(o$cell_type, o$sample, sep = "__")
fc  <- factor(cts)
mmc <- Matrix::sparse.model.matrix(~ 0 + fc); colnames(mmc) <- levels(fc)
pb_ct <- as.matrix(counts %*% mmc)                          # genes x (celltype__sample)

# --- metadata ---
smeta <- unique(o@meta.data[, c("sample", "tissue", "replicate")])
rownames(smeta) <- smeta$sample
ncell <- as.data.frame(table(cell_type = o$cell_type, sample = o$sample))  # cells per group

saveRDS(list(pb_sample = pb_sample, pb_ct = pb_ct, smeta = smeta,
             ncell = ncell, genes = keep),
        file.path(outdir, "pseudobulk_data.rds"))

cat("samples:", paste(colnames(pb_sample), collapse = ", "), "\n")
cat("cell_type x sample groups:", ncol(pb_ct), "\n")
message("==== annot7 pseudobulk prep done ====")
### END
