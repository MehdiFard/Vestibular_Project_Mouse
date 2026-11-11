#!/usr/bin/env Rscript
# 03b_pseudobulk_check.R  (deconv_env)  -- method sanity check
# Build pseudobulk from the single cells themselves (sum counts per donor), whose
# TRUE composition is known, then deconvolve with Bisque & MuSiC and compare to
# truth. This isolates method mechanics from the nuclei-vs-whole-cell domain shift:
#   * recovers truth  -> mechanics fine; real-bulk failure = domain shift
#   * fails on truth  -> method setup is the problem
suppressMessages({ library(Biobase); library(BisqueRNA); library(MuSiC)
                   library(SingleCellExperiment); library(Matrix); library(ggplot2) })
IN  <- "results/inputs"
TAB <- "results/tables"
PLT <- "results/plots"
dir.create(PLT, showWarnings = FALSE, recursive = TRUE)
KEEP <- c("Hair_cell","Supporting_cell","Transitional_epith",
          "Glia_Schwann","Fibrocyte_mesench","Roof_dark_cell")
sc <- readRDS(file.path(IN,"sc_counts.rds")); meta <- readRDS(file.path(IN,"sc_meta.rds"))

# TRUE composition per donor
truth <- prop.table(table(meta$SubjectName, meta$cell_type), 1)[, KEEP]

# pseudobulk = sum counts per donor (genes x 4), from ALL cells (true composition)
subs <- sort(unique(meta$SubjectName))
pb <- sapply(subs, function(s) Matrix::rowSums(sc[, rownames(meta)[meta$SubjectName==s], drop=FALSE]))
rownames(pb) <- rownames(sc)
pb <- as.matrix(pb)
cat("=== TRUTH (donor x cell type) ===\n"); print(round(truth,3))

# reference: downsample to <=2000 cells/type (as in the real run), keeps memory sane
set.seed(1)
ref_cells <- unlist(lapply(split(rownames(meta), meta$cell_type), function(ids)
    if (length(ids) > 2000) sample(ids, 2000) else ids))
sc_c <- sc[, ref_cells]; meta_r <- meta[ref_cells, ]

# Bisque
pheno <- AnnotatedDataFrame(data.frame(cellType=meta_r$cell_type, SubjectName=meta_r$SubjectName,
                                       row.names=colnames(sc_c)))
bq <- ReferenceBasedDecomposition(ExpressionSet(assayData=pb),
        ExpressionSet(assayData=as.matrix(sc_c), phenoData=pheno),
        markers=NULL, use.overlap=FALSE)$bulk.props
cat("\n=== BISQUE estimate ===\n"); print(round(t(bq)[,KEEP],3))
cat("Bisque Pearson vs truth:", round(cor(as.vector(truth), as.vector(t(bq)[,KEEP])),3), "\n")

# MuSiC
sce <- SingleCellExperiment(assays=list(counts=sc_c),
        colData=DataFrame(cellType=meta_r$cell_type, SubjectName=meta_r$SubjectName))
mu <- music_prop(bulk.mtx=pb, sc.sce=sce, clusters="cellType", samples="SubjectName",
                 select.ct=KEEP, verbose=FALSE)$Est.prop.weighted
cat("\n=== MuSiC estimate ===\n"); print(round(mu[,KEEP],3))
cat("MuSiC Pearson vs truth:", round(cor(as.vector(truth), as.vector(mu[,KEEP])),3), "\n")

# ---- save validation table + estimated-vs-truth plot ------------------------
mk_long <- function(mat, method) data.frame(
    donor = rep(rownames(mat), length(KEEP)),
    cell_type = rep(KEEP, each = nrow(mat)),
    estimated = as.vector(as.matrix(mat[, KEEP])),
    truth = as.vector(as.matrix(truth[rownames(mat), KEEP])),
    method = method)
val <- rbind(mk_long(t(bq), "Bisque"), mk_long(mu, "MuSiC"))
write.csv(val, file.path(TAB, "pseudobulk_validation.csv"), row.names = FALSE)

rlab <- sapply(c("Bisque","MuSiC"), function(m){ d<-val[val$method==m,]
    sprintf("%s: r=%.2f", m, cor(d$truth,d$estimated)) })
val$facet <- rlab[val$method]
p <- ggplot(val, aes(truth, estimated, color = cell_type)) +
    geom_abline(lty = 2, color = "grey60") +
    geom_point(size = 2.5) + facet_wrap(~ facet) +
    labs(title = "Pseudobulk self-recovery (known truth from the single cells)",
         subtitle = "Deconvolving pseudobulk validates method mechanics",
         x = "true proportion", y = "estimated proportion") +
    coord_equal() + theme_bw()
ggsave(file.path(PLT, "00_pseudobulk_validation.png"), p, width = 9, height = 4.8, dpi = 200, bg = "white")
cat("\nsaved pseudobulk_validation.csv + 00_pseudobulk_validation.png\nDONE\n")
