# 08 - limma-voom cross-check of the DESeq2 results.
#
# Same libraries, same filters, same designs - a different statistical engine.
# DESeq2 fits a negative-binomial GLM to counts; limma-voom transforms to
# log-CPM, estimates a mean-variance weight per observation and fits a normal
# linear model. If the two agree, the results are not an artefact of either
# package's assumptions.
#
# Batch stays in the design here too, for the reasons measured in script 06.

source("../../analyses/05_dge/scripts/00_config.R")

library(data.table)
library(limma)
library(edgeR)
library(ggplot2)

samples <- fread(file.path(DATA, "sample_table.tsv"))
samples[, organ := factor(organ, levels = c("Saccule", "Utricle"))]
samples[, age   := factor(age,   levels = c("w8", "w22"))]
samples[, batch := factor(batch)]
samples[, group := factor(paste(organ, age), levels = names(GROUP_COLS))]

txi  <- readRDS(file.path(DATA, "txi_gene.rds"))
anno <- fread(file.path(DATA, "gene_annotation.tsv"))
de   <- readRDS(file.path(DATA, "de_results.rds"))

KEEP_IDS <- anno[biotype %in% KEEP_BIOTYPES, gene_id]

run_limma <- function(s, design, coef, label) {
  cat("\n=====", label, "=====\n")
  m <- txi$counts[, match(s$library_id, colnames(txi$counts))]
  m <- m[rownames(m) %in% KEEP_IDS, ]
  gs <- if (grepl("organ", label)) table(s$organ) else table(s$age)
  m <- m[rowSums(m >= 10) >= min(gs), ]

  dge <- calcNormFactors(DGEList(round(m)))
  v   <- voom(dge, design)
  fit <- eBayes(lmFit(v, design))
  tt  <- topTable(fit, coef = coef, number = Inf, sort.by = "none")

  out <- data.table(gene_id = rownames(tt),
                    log2FoldChange = tt$logFC,
                    AveExpr = tt$AveExpr,
                    pvalue = tt$P.Value,
                    padj = tt$adj.P.Val)
  out <- merge(out, anno[, .(gene_id, gene_name, biotype, chr)],
               by = "gene_id", all.x = TRUE, sort = FALSE)
  setorder(out, padj, na.last = TRUE)

  sig <- out[padj < FDR_CUT & abs(log2FoldChange) > LFC_CUT]
  cat("Genes tested:", nrow(out), " FDR<0.05:", sum(out$padj < FDR_CUT),
      " +|LFC|>1:", nrow(sig), "\n")
  out
}

limma_res <- list()

# ---------------------------------------------------------------- A & B  organ, paired
for (a in c("w8", "w22")) {
  s_age  <- samples[age == a & !is.na(strain)]
  # pairs = animals with both organs, as in 03_differential_expression.R (pair_id from 01_build_inputs.R)
  paired <- s_age[, .(n = .N, k = uniqueN(organ)), by = pair_id][n == 2 & k == 2, pair_id]
  s <- s_age[pair_id %in% paired][, strain := factor(pair_id)][]
  design <- model.matrix(~ strain + organ, data = s)
  key <- paste0(ifelse(a == "w8", "A", "B"), "_utricle_vs_saccule_", a)
  limma_res[[key]] <- run_limma(s, design, "organUtricle",
                                paste("Utricle vs Saccule at", a, "(organ, paired)"))
}

# ---------------------------------------------------------------- C & D  age
for (o in c("Saccule", "Utricle")) {
  s <- samples[organ == o][, batch := droplevels(batch)][]
  design <- model.matrix(~ batch + age, data = s)
  key <- paste0(ifelse(o == "Saccule", "C", "D"), "_w22_vs_w8_", tolower(o))
  limma_res[[key]] <- run_limma(s, design, "agew22",
                                paste("w22 vs w8 in", o, "(age)"))
}

# ---------------------------------------------------------------- E  interaction
s <- copy(samples)
design <- model.matrix(~ batch + organ * age, data = s)
limma_res[["E_organ_by_age_interaction"]] <-
  run_limma(s, design, "organUtricle:agew22", "Organ x age interaction (age)")

for (k in names(limma_res))
  fwrite(limma_res[[k]], file.path(TABLES, paste0("limma_", k, ".tsv.gz")),
         sep = "\t", na = "NA")

# ---------------------------------------------------------------- agreement
NICE <- c(A_utricle_vs_saccule_w8   = "Utricle vs Saccule, w8",
          B_utricle_vs_saccule_w22  = "Utricle vs Saccule, w22",
          C_w22_vs_w8_saccule       = "w22 vs w8, Saccule",
          D_w22_vs_w8_utricle       = "w22 vs w8, Utricle",
          E_organ_by_age_interaction = "Organ x age interaction")

cmp <- rbindlist(lapply(names(NICE), function(k) {
  d <- de[[k]]; l <- limma_res[[k]]
  m <- merge(d[, .(gene_id, d_lfc = log2FoldChange, d_padj = padj)],
             l[, .(gene_id, l_lfc = log2FoldChange, l_padj = padj)], by = "gene_id")
  m <- m[!is.na(d_lfc) & !is.na(l_lfc)]
  ds <- m[!is.na(d_padj) & d_padj < FDR_CUT & abs(d_lfc) > LFC_CUT, gene_id]
  ls <- m[!is.na(l_padj) & l_padj < FDR_CUT & abs(l_lfc) > LFC_CUT, gene_id]
  data.table(analysis = k, contrast = NICE[k],
             genes_compared = nrow(m),
             deseq2_sig = length(ds), limma_sig = length(ls),
             shared = length(intersect(ds, ls)),
             jaccard = round(length(intersect(ds, ls)) / length(union(ds, ls)), 3),
             r_lfc = round(cor(m$d_lfc, m$l_lfc), 3),
             rho_rank = round(cor(m$d_lfc, m$l_lfc, method = "spearman"), 3))
}))
cmp[, pct_of_deseq2_recovered := round(100 * shared / pmax(deseq2_sig, 1), 1)]

cat("\n\n================ DESeq2 vs limma-voom ================\n")
print(cmp[, .(contrast, deseq2_sig, limma_sig, shared, pct_of_deseq2_recovered,
              jaccard, r_lfc)])
fwrite(cmp, file.path(TABLES, "limma_vs_deseq2_summary.tsv"), sep = "\t")

# ---------------------------------------------------------------- figures
lfc_dt <- rbindlist(lapply(names(NICE), function(k) {
  m <- merge(de[[k]][, .(gene_id, d_lfc = log2FoldChange, d_padj = padj)],
             limma_res[[k]][, .(gene_id, l_lfc = log2FoldChange, l_padj = padj)],
             by = "gene_id")
  m <- m[!is.na(d_lfc) & !is.na(l_lfc)]
  m[, sig := fifelse(!is.na(d_padj) & d_padj < FDR_CUT & abs(d_lfc) > LFC_CUT &
                     !is.na(l_padj) & l_padj < FDR_CUT & abs(l_lfc) > LFC_CUT, "both",
              fifelse(!is.na(d_padj) & d_padj < FDR_CUT & abs(d_lfc) > LFC_CUT, "DESeq2 only",
              fifelse(!is.na(l_padj) & l_padj < FDR_CUT & abs(l_lfc) > LFC_CUT, "limma only",
                      "neither")))]
  m[, contrast := NICE[k]]
  m[, lab := sprintf("%s\nr = %.3f", NICE[k], cor(d_lfc, l_lfc))]
  m
}))
lfc_dt[, sig := factor(sig, levels = c("neither", "DESeq2 only", "limma only", "both"))]

p_lfc <- ggplot(lfc_dt[order(sig)], aes(d_lfc, l_lfc, colour = sig)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(size = 0.7, alpha = 0.55) +
  facet_wrap(~ lab, scales = "free", nrow = 2) +
  scale_colour_manual(values = c(neither = "grey82", `DESeq2 only` = "#4393C3",
                                 `limma only` = "#D6604D", both = "#1A9850")) +
  labs(title = "DESeq2 vs limma-voom: log2 fold changes agree across every contrast",
       subtitle = "Same libraries, same filters, same designs; dashed line is y = x",
       x = "log2FC - DESeq2 (ashr-shrunken)", y = "log2FC - limma-voom", colour = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

ggsave(file.path(FIGS, "23_limma_vs_deseq2_lfc.png"), p_lfc,
       width = 12, height = 7, dpi = 300, bg = "white")

ov <- melt(cmp[, .(contrast, `DESeq2 only` = deseq2_sig - shared,
                   shared, `limma only` = limma_sig - shared)],
           id.vars = "contrast", variable.name = "set", value.name = "n")
ov[, contrast := factor(contrast, levels = rev(NICE))]

p_ov <- ggplot(ov, aes(n, contrast, fill = set)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = ifelse(n > 0, n, "")), position = position_stack(vjust = 0.5),
            size = 3, colour = "white") +
  scale_fill_manual(values = c(`DESeq2 only` = "#4393C3", shared = "#1A9850",
                               `limma only` = "#D6604D")) +
  labs(title = "Significant genes: recovered by both engines, or only one",
       subtitle = sprintf("FDR < %.2f and |log2FC| > %g", FDR_CUT, LFC_CUT),
       x = "Number of genes", y = NULL, fill = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

ggsave(file.path(FIGS, "24_limma_vs_deseq2_overlap.png"), p_ov,
       width = 9.5, height = 5, dpi = 300, bg = "white")

cat("\nDone.\n")
