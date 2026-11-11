# annot12_clean_pseudobulk.R  (seurat_env)
# Pseudobulk prep from the CLEANED combined object (all_integrated_clean.rds),
# protein_coding + lncRNA only. Mirrors annot7 but on cleaned data.
rm(list = ls())
suppressMessages({ library(Seurat); library(Matrix) })

indir  <- "results/cell_annotation_clean"
outdir <- file.path(indir, "pseudobulk_DE")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

bt <- read.delim("src/annot_ref/gene_biotypes_GRCm39.tsv", header = FALSE,
                 col.names = c("gene", "type"))
keep_genes <- unique(bt$gene[bt$type %in% c("protein_coding", "lncRNA")])

o <- readRDS(file.path(indir, "all_integrated_clean.rds"))
DefaultAssay(o) <- "RNA"
o <- tryCatch(JoinLayers(o), error = function(e) o)
keep <- intersect(keep_genes, rownames(o))
cat("genes kept (pc+lncRNA):", length(keep), "of", nrow(o), "\n")
counts <- GetAssayData(o, assay = "RNA", layer = "counts")[keep, ]

fs <- factor(o$sample)
mm <- Matrix::sparse.model.matrix(~ 0 + fs); colnames(mm) <- levels(fs)
pb_sample <- as.matrix(counts %*% mm)

cts <- paste(o$cell_type, o$sample, sep = "__"); fc <- factor(cts)
mmc <- Matrix::sparse.model.matrix(~ 0 + fc); colnames(mmc) <- levels(fc)
pb_ct <- as.matrix(counts %*% mmc)

smeta <- unique(o@meta.data[, c("sample", "tissue", "replicate")])
rownames(smeta) <- smeta$sample
ncell <- as.data.frame(table(cell_type = o$cell_type, sample = o$sample))

saveRDS(list(pb_sample = pb_sample, pb_ct = pb_ct, smeta = smeta,
             ncell = ncell, genes = keep),
        file.path(outdir, "pseudobulk_data_clean.rds"))
cat("samples:", paste(colnames(pb_sample), collapse = ", "), "\n")
cat("cell_type x sample groups:", ncol(pb_ct), "\n")
message("==== annot12 clean pseudobulk prep done ====")
### END
