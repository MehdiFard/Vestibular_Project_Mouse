# annot8_pseudobulk_DE.R  (run with Env1 — has edgeR/DESeq2/limma/ggplot2)
# Pseudobulk PCA + proper (edgeR) differential expression, utricle vs saccule:
#   1. organ level  : PCA of the 4 samples + DE (all cells pooled per sample)
#   2. cell level   : DE within each cell type (pseudobulk 2 vs 2)
# Genes already restricted to protein_coding + lncRNA in annot7.
# setwd("analyses/07_single_cell_rna")   # from the project root

rm(list = ls())
suppressMessages({ library(edgeR); library(ggplot2); library(matrixStats) })

outdir <- "results/pseudobulk_DE"
d <- readRDS(file.path(outdir, "pseudobulk_data.rds"))
pb_sample <- d$pb_sample; pb_ct <- d$pb_ct; smeta <- d$smeta; ncell <- d$ncell

FDR <- 0.05; LFC <- 1                       # significance thresholds
regional <- c("Otx1", "Otx2", "Emx2", "Tbx1", "Lfng", "Bmp4", "Gbx2", "Gata3")
COL <- c(up_utricle = "#c0392b", up_saccule = "#2e86ab", ns = "grey80")

# volcano helper (base ggplot2, no ggrepel; regional genes highlighted gold)
volcano <- function(res, title, file, highlight = regional) {
    res$dir <- ifelse(res$FDR < FDR & res$logFC >  LFC, "up_utricle",
               ifelse(res$FDR < FDR & res$logFC < -LFC, "up_saccule", "ns"))
    res$nl  <- -log10(pmax(res$FDR, 1e-300))
    ro  <- res[order(res$FDR), ]
    lab <- head(ro[ro$dir != "ns", ], 20)
    hl  <- res[res$gene %in% highlight, ]
    p <- ggplot(res, aes(logFC, nl, color = dir)) +
        geom_point(size = 0.7, alpha = 0.6) +
        scale_color_manual(values = COL) +
        geom_vline(xintercept = c(-LFC, LFC), linetype = "dashed", color = "grey60") +
        geom_hline(yintercept = -log10(FDR), linetype = "dashed", color = "grey60") +
        geom_text(data = lab, aes(label = gene), size = 2.5, vjust = -0.6,
                  check_overlap = TRUE, show.legend = FALSE) +
        labs(title = title, subtitle = "red = up in utricle · blue = up in saccule · gold = regional",
             x = "log2FC (utricle / saccule)", y = "-log10 FDR") +
        theme_classic()
    if (nrow(hl)) p <- p +
        geom_point(data = hl, color = "gold3", size = 2) +
        geom_text(data = hl, aes(label = gene), color = "gold4", fontface = "bold",
                  size = 3, vjust = 1.4, check_overlap = TRUE, show.legend = FALSE)
    ggsave(file, p, width = 7.5, height = 6, dpi = 200, bg = "white")
}

# PCA showed the dominant axis is the batch (rep1 vs rep2), not tissue — so we
# BLOCK on replicate: design ~ batch + tissue, and test the tissue coefficient.
edgeR_DE <- function(counts, tissue, batch) {    # both factors, aligned to columns
    y <- DGEList(counts, group = tissue)
    y <- y[filterByExpr(y, group = tissue), , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y)
    des <- if (nlevels(droplevels(batch)) > 1) model.matrix(~ batch + tissue)
           else                                  model.matrix(~ tissue)
    y <- estimateDisp(y, des)
    qlf <- glmQLFTest(glmQLFit(y, des), coef = ncol(des))   # last col = tissue effect
    list(y = y, res = topTags(qlf, n = Inf)$table)
}


### 1. ORGAN LEVEL ---------------------------------------------------------
sc     <- colnames(pb_sample)
tissue <- factor(smeta[sc, "tissue"], levels = c("saccule", "utricle"))
repl   <- smeta[sc, "replicate"]

batch <- factor(repl)
fit <- edgeR_DE(pb_sample, tissue, batch)
lcpm <- cpm(fit$y, log = TRUE, prior.count = 1)

# helper: draw a PCA plot from a logCPM matrix
pca_plot <- function(mat, title, file) {
    v  <- rowVars(mat); tg <- head(order(v, decreasing = TRUE), 2000)
    pc <- prcomp(t(mat[tg, ]), scale. = TRUE)
    ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
    df <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                     sample = sc, tissue = tissue, replicate = repl)
    p <- ggplot(df, aes(PC1, PC2, color = tissue, shape = replicate)) +
        geom_point(size = 6) +
        geom_text(aes(label = sample), vjust = -1.4, size = 3.2, show.legend = FALSE) +
        scale_color_manual(values = c(saccule = "#c0392b", utricle = "#2e86ab")) +
        labs(title = title, x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)")) +
        coord_cartesian(clip = "off") + theme_classic()
    ggsave(file, p, width = 7.5, height = 6, dpi = 200, bg = "white")
}
# raw PCA (batch dominates) and batch-corrected PCA (does tissue separate?)
pca_plot(lcpm, "Pseudobulk PCA — organ level (protein-coding + lncRNA)",
         file.path(outdir, "PCA_organ_pseudobulk.png"))
lcpm_bc <- limma::removeBatchEffect(lcpm, batch = batch)
pca_plot(lcpm_bc, "Pseudobulk PCA — after removing replicate/batch effect",
         file.path(outdir, "PCA_organ_batch_corrected.png"))

# --- sample-sample correlation heatmap ---
cc <- cor(lcpm); ccdf <- expand.grid(s1 = rownames(cc), s2 = colnames(cc)); ccdf$r <- as.vector(cc)
p <- ggplot(ccdf, aes(s1, s2, fill = r)) + geom_tile() +
    geom_text(aes(label = round(r, 3)), size = 3) +
    scale_fill_gradient(low = "white", high = "#2e86ab") +
    labs(title = "Sample correlation (organ pseudobulk logCPM)", x = NULL, y = NULL) +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(outdir, "sample_correlation.png"), p, width = 6, height = 5, dpi = 200, bg = "white")

# --- organ DE ---
res <- fit$res; res$gene <- rownames(res)
write.csv(res, file.path(outdir, "organ_DE_utricle_vs_saccule_edgeR.csv"), row.names = FALSE)
n_up_u <- sum(res$FDR < FDR & res$logFC >  LFC)
n_up_s <- sum(res$FDR < FDR & res$logFC < -LFC)
cat(sprintf("ORGAN DE (edgeR, ~batch+tissue): up_utricle=%d  up_saccule=%d  (FDR<%.2f,|log2FC|>%d)\n",
            n_up_u, n_up_s, FDR, LFC))
cat("top nominal organ genes (batch-blocked):\n")
print(head(res[order(res$PValue), c("gene", "logFC", "PValue", "FDR")], 15))
volcano(res, "Organ level — utricle vs saccule (edgeR pseudobulk 2v2)",
        file.path(outdir, "volcano_organ.png"))


### 2. CELL LEVEL (per cell type, pseudobulk 2 vs 2) -----------------------
parts <- strsplit(colnames(pb_ct), "__", fixed = TRUE)
ctv <- vapply(parts, `[`, "", 1); smv <- vapply(parts, `[`, "", 2)
nc  <- setNames(ncell$Freq, paste(ncell$cell_type, ncell$sample, sep = "__"))
MINC <- 20                                    # min cells per pseudobulk sample

summ <- data.frame()
for (type in unique(ctv)) {
    idx <- which(ctv == type)
    idx <- idx[nc[colnames(pb_ct)[idx]] >= MINC]
    if (length(idx) < 4) next
    tt <- factor(smeta[smv[idx], "tissue"], levels = c("saccule", "utricle"))
    bb <- factor(smeta[smv[idx], "replicate"])
    if (sum(tt == "saccule") < 2 || sum(tt == "utricle") < 2) next
    fitc <- tryCatch(edgeR_DE(pb_ct[, idx, drop = FALSE], tt, bb),
                     error = function(e) { message("  ", type, ": ", e$message); NULL })
    if (is.null(fitc)) next
    rc <- fitc$res; rc$gene <- rownames(rc); rc$cell_type <- type
    write.csv(rc, file.path(outdir, paste0("cellDE_", type, "_edgeR.csv")), row.names = FALSE)
    volcano(rc, paste0(type, " — utricle vs saccule (edgeR)"),
            file.path(outdir, paste0("volcano_cell_", type, ".png")))
    summ <- rbind(summ, data.frame(cell_type = type,
        up_saccule = sum(rc$FDR < FDR & rc$logFC < -LFC),
        up_utricle = sum(rc$FDR < FDR & rc$logFC >  LFC)))
}
write.csv(summ, file.path(outdir, "cellDE_significant_counts.csv"), row.names = FALSE)
cat("---- cell-level significant genes (edgeR) ----\n"); print(summ)

if (nrow(summ)) {
    long <- rbind(data.frame(cell_type = summ$cell_type, direction = "up_saccule", n = summ$up_saccule),
                  data.frame(cell_type = summ$cell_type, direction = "up_utricle", n = summ$up_utricle))
    p <- ggplot(long, aes(reorder(cell_type, n), n, fill = direction)) +
        geom_col(position = "dodge") + coord_flip() +
        scale_fill_manual(values = c(up_saccule = "#2e86ab", up_utricle = "#c0392b")) +
        labs(title = "Cell-level DE (edgeR pseudobulk 2v2): significant genes/cell type",
             x = NULL, y = paste0("n significant (FDR<", FDR, ", |log2FC|>", LFC, ")")) +
        theme_classic()
    ggsave(file.path(outdir, "cellDE_summary_bar.png"), p, width = 8, height = 5, dpi = 200, bg = "white")
}

message("==== annot8 pseudobulk DE done ====")
### END
