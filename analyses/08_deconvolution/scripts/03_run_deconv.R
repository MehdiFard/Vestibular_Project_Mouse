#!/usr/bin/env Rscript
# 03_run_deconv.R  (run in deconv_env)
# ---------------------------------------------------------------------------
# Estimate cell-type composition of each bulk sample, per tissue, with 3 methods:
#   * Bisque  (PRIMARY)     -- learns a reference->bulk transform (ReferenceBasedDecomposition)
#   * MuSiC   (SENSITIVITY) -- weights genes by cross-donor consistency
#   * ssGSEA  (enrichment)  -- GTEx-style marker enrichment score (GSVA)
#
# Bulk and single-cell are DIFFERENT animals, so Bisque runs with use.overlap=FALSE.
# The sc reference is downsampled to <=2000 cells/type (keeps all 4 donors) so the
# dense reference stays small; results are written per method x tissue.
#
# Inputs : results/inputs/{sc_counts.rds, sc_meta.rds, sc_markers.rds,
#          bulk_tpm_<tissue>.csv}
# Outputs: results/tables/{fractions_bisque,fractions_music,enrichment_ssgsea}_<tissue>.csv
#
# Run on a compute node (not the login node).
# ---------------------------------------------------------------------------
suppressMessages({
    library(Biobase); library(BisqueRNA); library(MuSiC)
    library(SingleCellExperiment); library(GSVA); library(Matrix)
})

IN  <- "results/inputs"
OUT <- "results/tables"
BIOTYPE <- "../07_single_cell_rna/src/annot_ref/gene_biotypes_GRCm39.tsv"   # symbol <tab> biotype
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
KEEP <- c("Hair_cell", "Supporting_cell", "Transitional_epith",
          "Glia_Schwann", "Fibrocyte_mesench", "Roof_dark_cell")
MAX_PER_TYPE <- 2000
set.seed(1)

# Gene whitelist for deconvolution. A nuclear (snRNA) reference deconvolving
# whole-cell bulk is badly distorted by technical genes that differ between the
# two: mitochondrial genes (huge in bulk, ~0 in nuclei), ribosomal genes, and
# miRNA/snoRNA/snRNA/pseudogenes (e.g. Mir6236, Gm* small RNAs that dominate the
# bulk TPM). We restrict to protein-coding genes and drop mt-/ribosomal.
bt <- read.table(BIOTYPE, sep = "\t", header = FALSE,
                 col.names = c("symbol", "biotype"), quote = "", stringsAsFactors = FALSE)
PC <- bt$symbol[bt$biotype == "protein_coding"]
gene_ok <- function(g) {
    g <- intersect(g, PC)                                   # protein-coding only
    g[!grepl("^mt-", g, ignore.case = TRUE) &               # drop mitochondrial
      !grepl("^Rp[sl]", g)]                                 # drop ribosomal
}

cat("[load] single-cell reference\n")
sc_counts <- readRDS(file.path(IN, "sc_counts.rds"))          # genes x cells
sc_meta   <- readRDS(file.path(IN, "sc_meta.rds"))            # cell_type, SubjectName, tissue
markers   <- readRDS(file.path(IN, "sc_markers.rds"))         # list per cell type

# downsample reference to <= MAX_PER_TYPE cells / type (balanced, all donors kept)
keep_cells <- unlist(lapply(split(rownames(sc_meta), sc_meta$cell_type), function(ids)
    if (length(ids) > MAX_PER_TYPE) sample(ids, MAX_PER_TYPE) else ids))
sc_counts <- sc_counts[, keep_cells]
sc_meta   <- sc_meta[keep_cells, ]
cat("    reference cells:", ncol(sc_counts), "\n"); print(table(sc_meta$cell_type))

run_tissue <- function(tissue) {
    cat("\n==================  ", tissue, "  ==================\n")
    bulk <- as.matrix(read.csv(file.path(IN, paste0("bulk_tpm_", tissue, ".csv")),
                               row.names = 1, check.names = FALSE))
    common <- gene_ok(intersect(rownames(bulk), rownames(sc_counts)))
    cat("[genes] bulk", nrow(bulk), "| sc", nrow(sc_counts),
        "| common protein-coding (no mt/ribo)", length(common), "\n")
    bulk_c <- bulk[common, , drop = FALSE]
    sc_c   <- sc_counts[common, , drop = FALSE]

    # ------- Bisque (PRIMARY) -------
    cat("[Bisque] ReferenceBasedDecomposition (use.overlap=FALSE)\n")
    bulk.eset <- ExpressionSet(assayData = bulk_c)
    pheno <- AnnotatedDataFrame(data.frame(cellType = sc_meta$cell_type,
                                           SubjectName = sc_meta$SubjectName,
                                           row.names = colnames(sc_c)))
    sc.eset <- ExpressionSet(assayData = as.matrix(sc_c), phenoData = pheno)
    bq <- tryCatch({
        r <- ReferenceBasedDecomposition(bulk.eset, sc.eset, markers = NULL,
                                         use.overlap = FALSE)
        t(r$bulk.props)                                   # samples x cell types
    }, error = function(e) { cat("  Bisque FAILED:", conditionMessage(e), "\n"); NULL })
    if (!is.null(bq)) write.csv(bq, file.path(OUT, paste0("fractions_bisque_", tissue, ".csv")))

    # ------- MuSiC (SENSITIVITY) -------
    cat("[MuSiC] music_prop\n")
    mu <- tryCatch({
        sce <- SingleCellExperiment(assays = list(counts = sc_c),
                  colData = DataFrame(cellType = sc_meta$cell_type,
                                      SubjectName = sc_meta$SubjectName))
        r <- music_prop(bulk.mtx = bulk_c, sc.sce = sce,
                        clusters = "cellType", samples = "SubjectName",
                        select.ct = KEEP, verbose = FALSE)
        r$Est.prop.weighted                               # samples x cell types
    }, error = function(e) { cat("  MuSiC FAILED:", conditionMessage(e), "\n"); NULL })
    if (!is.null(mu)) write.csv(mu, file.path(OUT, paste0("fractions_music_", tissue, ".csv")))

    # ------- ssGSEA enrichment (GTEx-style) -------
    cat("[ssGSEA] GSVA on log2(TPM+1) with sc marker sets\n")
    gs <- lapply(markers, function(g) intersect(g, common))
    logtpm <- log2(bulk_c + 1)
    en <- tryCatch({
        # version-robust GSVA call (new param API vs legacy signature)
        r <- tryCatch(
            GSVA::gsva(GSVA::ssgseaParam(logtpm, gs), verbose = FALSE),
            error = function(e) GSVA::gsva(logtpm, gs, method = "ssgsea", verbose = FALSE))
        t(r)                                              # samples x cell types
    }, error = function(e) { cat("  ssGSEA FAILED:", conditionMessage(e), "\n"); NULL })
    if (!is.null(en)) write.csv(en, file.path(OUT, paste0("enrichment_ssgsea_", tissue, ".csv")))

    cat("[done]", tissue, "\n")
}

for (t in c("saccule", "utricle")) run_tissue(t)
cat("\nALL DONE\n")
