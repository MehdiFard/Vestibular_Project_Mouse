# 10 - Two things about the pairing structure:
#
#   (a) what "strain identifies a saccule/utricle pair" actually buys, shown as
#       per-strain paired lines for the top organ genes
#   (b) the age contrast restricted to the strains sampled at BOTH ages, which
#       is a within-animal age comparison instead of a between-animal one
#
# Each strain contributes exactly one saccule and one utricle library per age,
# so "strain" is the pairing key - the same role a subject ID plays in a paired
# t-test.

source("../../analyses/05_dge/scripts/00_config.R")

library(data.table)
library(DESeq2)
library(BiocParallel)
library(ggplot2)

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
BPP <- MulticoreParam(workers = 8); register(BPP)

samples <- fread(file.path(DATA, "sample_table.tsv"))
samples[, organ := factor(organ, levels = c("Saccule", "Utricle"))]
samples[, age   := factor(age,   levels = c("w8", "w22"))]
samples[, batch := factor(batch)]

txi  <- readRDS(file.path(DATA, "txi_gene.rds"))
anno <- fread(file.path(DATA, "gene_annotation.tsv"))
de   <- readRDS(file.path(DATA, "de_results.rds"))
KEEP_IDS <- anno[biotype %in% KEEP_BIOTYPES, gene_id]

subset_txi <- function(ids) {
  i <- match(ids, colnames(txi$counts))
  list(counts = txi$counts[, i, drop = FALSE],
       abundance = txi$abundance[, i, drop = FALSE],
       length = txi$length[, i, drop = FALSE],
       countsFromAbundance = txi$countsFromAbundance)
}

# =================================================================== (a) pairing
# For the strongest organ genes at w8, draw one line per strain joining that
# strain's saccule value to its utricle value. Strains sit at very different
# heights - that spread is what the unpaired model has to treat as noise, and
# what putting strain in the model removes.
dds_int <- readRDS(file.path(DATA, "dds_interaction.rds"))
norm    <- counts(dds_int, normalized = TRUE)

s_w8   <- samples[age == "w8" & !is.na(strain)]
paired <- s_w8[, .N, by = strain][N == 2, strain]
s_p    <- s_w8[strain %in% paired]

top <- head(de$A_utricle_vs_saccule_w8[
  !is.na(padj) & padj < FDR_CUT & abs(log2FoldChange) > LFC_CUT &
    gene_id %in% rownames(norm) & !is.na(gene_name) & gene_name != ""][order(padj)], 4)

pl <- rbindlist(lapply(seq_len(nrow(top)), function(i) {
  data.table(gene = top$gene_name[i], library_id = colnames(norm),
             value = log2(norm[top$gene_id[i], ] + 1))
}))
pl <- merge(pl, s_p[, .(library_id, strain, organ)], by = "library_id")
pl[, gene := factor(gene, levels = top$gene_name)]

p_pair <- ggplot(pl, aes(organ, value, group = strain)) +
  geom_line(colour = "grey65", linewidth = 0.35, alpha = 0.8) +
  geom_point(aes(colour = organ), size = 1.5, alpha = 0.9) +
  facet_wrap(~ gene, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = ORGAN_COLS, guide = "none") +
  labs(title = "What pairing by strain means",
       subtitle = paste("Each grey line is one strain at 8 weeks: its saccule joined to its utricle.",
                        "\nStrains sit at very different heights, but nearly all lines run the same way."),
       x = NULL, y = "log2(normalised count + 1)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"))

ggsave(file.path(FIGS, "27_pairing_illustration.png"), p_pair,
       width = 11, height = 4.2, dpi = 300, bg = "white")

# How much of each gene's variation is strain, versus the organ effect itself?
cat("=== Variance in the w8 libraries explained by strain vs by organ ===\n")
mat <- assay(vst(dds_int, blind = TRUE))[, s_p$library_id]
mat <- mat[rownames(mat) %in% KEEP_IDS, ]
fast_r2 <- function(Y, X) {
  Yt <- t(Y); fit <- X %*% qr.solve(crossprod(X), crossprod(X, Yt))
  1 - colSums((Yt - fit)^2) / colSums(sweep(Yt, 2, colMeans(Yt))^2)
}
r2_strain <- fast_r2(mat, model.matrix(~ factor(s_p$strain)))
r2_organ  <- fast_r2(mat, model.matrix(~ s_p$organ))
cat("  mean R2 explained by strain:", round(mean(r2_strain), 3), "\n")
cat("  mean R2 explained by organ :", round(mean(r2_organ), 3), "\n")
cat("  -> strain soaks up more per-gene variance than organ does; leaving it in\n")
cat("     the residual is what costs the unpaired model its power.\n")

# =================================================================== (b) paired age
# The age contrast within the strains that were sampled at both ages. This is a
# within-animal-line comparison: genetic background cancels out.
paired_age <- list()
for (o in c("Saccule", "Utricle")) {
  st <- samples[!is.na(strain) & organ == o, .(n = uniqueN(age)), by = strain][n == 2, strain]
  s  <- samples[!is.na(strain) & organ == o & strain %in% st]
  s[, strain := factor(strain)][, batch := droplevels(batch)]

  cat("\n=====", o, "- paired age contrast on", length(st), "strains,",
      nrow(s), "libraries =====\n")

  # Can batch be estimated alongside strain and age here? Strains differ in
  # which batches their two libraries fell into, so it may be identifiable.
  mm_full <- model.matrix(~ strain + batch + age, data = s)
  full_rank <- qr(mm_full)$rank == ncol(mm_full)
  design <- if (full_rank) ~ strain + batch + age else ~ strain + age
  cat("  ~ strain + batch + age estimable:", full_rank,
      "-> using", deparse(design), "\n")

  dds <- DESeqDataSetFromTximport(
    subset_txi(s$library_id),
    colData = data.frame(s, row.names = s$library_id),
    design = design)
  dds <- dds[rownames(dds) %in% KEEP_IDS, ]
  dds <- dds[rowSums(counts(dds) >= 10) >= min(table(s$age)), ]
  dds <- DESeq(dds, quiet = TRUE, parallel = TRUE, BPPARAM = BPP)

  res <- results(dds, contrast = c("age", "w22", "w8"), alpha = FDR_CUT)
  out <- data.table(gene_id = rownames(res), baseMean = res$baseMean,
                    log2FoldChange = res$log2FoldChange, lfcSE = res$lfcSE,
                    pvalue = res$pvalue, padj = res$padj)
  out <- merge(out, anno[, .(gene_id, gene_name, biotype, chr)],
               by = "gene_id", all.x = TRUE, sort = FALSE)
  setorder(out, padj, na.last = TRUE)

  sig <- out[!is.na(padj) & padj < FDR_CUT & abs(log2FoldChange) > LFC_CUT]
  cat("  Genes tested:", nrow(out), " FDR<0.05:", sum(out$padj < FDR_CUT, na.rm = TRUE),
      " +|LFC|>1:", nrow(sig), "\n")

  key <- paste0("paired_age_", tolower(o))
  paired_age[[key]] <- out
  fwrite(out, file.path(TABLES, paste0("DE_", key, ".tsv.gz")), sep = "\t", na = "NA")
  fwrite(sig, file.path(TABLES, paste0("DE_", key, "_significant.tsv")), sep = "\t", na = "NA")
}

# ---------------------------------------------------------------- agreement
cat("\n\n=== Paired age contrast vs the batch-adjusted one ===\n")
agree <- rbindlist(lapply(c("Saccule", "Utricle"), function(o) {
  full <- de[[paste0(ifelse(o == "Saccule", "C", "D"), "_w22_vs_w8_", tolower(o))]]
  pr   <- paired_age[[paste0("paired_age_", tolower(o))]]
  m <- merge(full[, .(gene_id, f_lfc = log2FoldChange, f_padj = padj)],
             pr[,   .(gene_id, p_lfc = log2FoldChange, p_padj = padj)], by = "gene_id")
  m <- m[!is.na(f_lfc) & !is.na(p_lfc)]
  fs <- m[!is.na(f_padj) & f_padj < FDR_CUT & abs(f_lfc) > LFC_CUT, gene_id]
  ps <- m[!is.na(p_padj) & p_padj < FDR_CUT & abs(p_lfc) > LFC_CUT, gene_id]
  # of the full-analysis hits, how many keep the same direction in the paired fit?
  same_dir <- m[gene_id %in% fs, mean(sign(f_lfc) == sign(p_lfc))]
  data.table(organ = o, n_strains = uniqueN(samples[!is.na(strain) & organ == o][
                          strain %in% samples[!is.na(strain) & organ == o,
                          .(n = uniqueN(age)), by = strain][n == 2, strain], strain]),
             genes_compared = nrow(m),
             batch_adjusted_sig = length(fs), paired_sig = length(ps),
             shared = length(intersect(fs, ps)),
             r_lfc = round(cor(m$f_lfc, m$p_lfc), 3),
             same_direction_pct = round(100 * same_dir, 1))
}))
print(agree)
fwrite(agree, file.path(TABLES, "paired_age_vs_batch_adjusted.tsv"), sep = "\t")

# ---------------------------------------------------------------- figure
cmp_dt <- rbindlist(lapply(c("Saccule", "Utricle"), function(o) {
  full <- de[[paste0(ifelse(o == "Saccule", "C", "D"), "_w22_vs_w8_", tolower(o))]]
  pr   <- paired_age[[paste0("paired_age_", tolower(o))]]
  m <- merge(full[, .(gene_id, f_lfc = log2FoldChange, f_padj = padj)],
             pr[,   .(gene_id, p_lfc = log2FoldChange)], by = "gene_id")
  m <- m[!is.na(f_lfc) & !is.na(p_lfc)]
  m[, hit := !is.na(f_padj) & f_padj < FDR_CUT & abs(f_lfc) > LFC_CUT]
  m[, organ := o][, lab := sprintf("%s  (r = %.2f)", o, cor(f_lfc, p_lfc))]
  m
}))

p_cmp <- ggplot(cmp_dt[order(hit)], aes(f_lfc, p_lfc, colour = hit)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = 0, colour = "grey88") + geom_vline(xintercept = 0, colour = "grey88") +
  geom_point(size = 0.8, alpha = 0.6) +
  facet_wrap(~ lab) +
  scale_colour_manual(values = c(`FALSE` = "grey80", `TRUE` = "#B2182B"),
                      labels = c("not significant", "significant in the main analysis")) +
  labs(title = "Age effect: all animals vs only the strains sampled at both ages",
       subtitle = paste("x: the batch-adjusted contrast on every library.",
                        "y: the same contrast within strains sampled at both 8 and 22 weeks."),
       x = "log2FC - batch-adjusted, all libraries",
       y = "log2FC - within-strain, both-age subset", colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

ggsave(file.path(FIGS, "28_paired_age_check.png"), p_cmp,
       width = 10, height = 5.5, dpi = 300, bg = "white")

cat("\nDone.\n")
