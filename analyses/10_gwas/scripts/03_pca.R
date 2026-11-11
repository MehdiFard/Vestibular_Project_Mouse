#!/usr/bin/env Rscript

### 03 Genotype PCA of the phenotyped strains: one row per strain (mice of a strain share its genotype, so a PCA on
# mice would weight strains by their number of mice), LD-pruned SNPs, standardised genotypes.
# Each mouse gets its strain's scores as GWAS covariates.
# In:  intermediate/genotype_pca/<trait>_strains.raw (plink2 --export A, one mouse per strain),
#      intermediate/genotype/<trait>_clean_RefAltOk.fam (all mice)
# Out: results/genotype_pca/<trait>.eigenvec.tsv (strain scores, PC1-20), <trait>.variance.tsv (% variance per PC),
#      <trait>.covar.pcs (FID IID PC1..PC<n>, all mice, no header; FaST-LMM covariates),
#      figures/pca/<trait>_scree.jpeg, figures/pca/<trait>_PC1_PC2.jpeg
# Usage (from analyses/10_gwas): Rscript scripts/03_pca.R <trait> <number of PCs as covariates>
suppressPackageStartupMessages({ library(data.table); library(dplyr); library(ggplot2) })

args <- commandArgs(trailingOnly = TRUE)
trait <- args[1]; n_pcs <- as.integer(args[2])
out <- file.path("results/genotype_pca", trait)

raw <- fread(paste0("intermediate/genotype_pca/", trait, "_strains.raw"))
g <- as.matrix(raw[, -(1:6)])
rownames(g) <- raw$FID
g <- g[, apply(g, 2, function(x) var(x, na.rm = TRUE) > 0)]          # SNPs that vary among these strains
g <- apply(g, 2, function(x) { x[is.na(x)] <- mean(x, na.rm = TRUE); x })   # mean-impute missing calls
pc <- prcomp(g, center = TRUE, scale. = TRUE)

var_pct <- 100 * pc$sdev^2 / sum(pc$sdev^2)
variance <- data.table(PC = seq_along(var_pct), pct_variance = round(var_pct, 3), cumulative_pct = round(cumsum(var_pct), 3))
fwrite(variance, paste0(out, ".variance.tsv"), sep = "\t")
k <- min(20, ncol(pc$x))
scores <- data.table(FID = rownames(pc$x), round(pc$x[, 1:k], 6))
fwrite(scores, paste0(out, ".eigenvec.tsv"), sep = "\t")

# covariates for every mouse: its strain's PC1..PC<n>
fam <- fread(paste0("intermediate/genotype/", trait, "_clean_RefAltOk.fam"), header = FALSE, select = 1:2, col.names = c("FID", "IID"))
covar <- merge(fam, scores[, c("FID", paste0("PC", 1:n_pcs)), with = FALSE], by = "FID", sort = FALSE)
stopifnot(nrow(covar) == nrow(fam), !anyNA(covar))
fwrite(covar, paste0(out, ".covar.pcs"), sep = "\t", col.names = FALSE)

message(sprintf("%s: %d strains, %d SNPs; PC1-%d explain %.1f%% (PC1 %.1f%%)", trait, nrow(g), ncol(g), n_pcs,
                sum(var_pct[1:n_pcs]), var_pct[1]))

# Scree plot: % variance of PC1-20, PCs used as covariates highlighted
dir.create("figures/pca", recursive = TRUE, showWarnings = FALSE)
used_lab <- paste0("PC1-", n_pcs, " (covariates)")
sc <- variance[1:k] %>% mutate(used = factor(ifelse(PC <= n_pcs, used_lab, "not used"), levels = c(used_lab, "not used")))
pl <- ggplot(sc, aes(PC, pct_variance)) +
  geom_line(colour = "grey60") + geom_point(aes(colour = used), size = 2.5) +
  scale_colour_manual(values = setNames(c("firebrick", "grey40"), c(used_lab, "not used")), name = NULL) +
  scale_x_continuous(breaks = 1:k) +
  labs(title = paste0("Scree plot - ", trait, " (", nrow(g), " strains)"), x = "PC", y = "% variance explained") +
  theme_minimal() + theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(paste0("figures/pca/", trait, "_scree.jpeg"), pl, width = 6, height = 6, bg = "white")

# PC1 vs PC2, BXD strains vs the other inbred strains
scores[, panel := fifelse(grepl("^BXD", FID), "BXD", fifelse(FID %in% c("C57BL.6J", "DBA.2J"), "BXD parent", "other inbred"))]
pl <- ggplot(scores, aes(PC1, PC2, colour = panel)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dotted") + geom_vline(xintercept = 0, linetype = "dotted") +
  labs(title = paste0("PCA - ", trait), colour = NULL,
       x = sprintf("PC 1 (%.1f %%)", var_pct[1]), y = sprintf("PC 2 (%.1f %%)", var_pct[2])) +
  theme_minimal() + theme(legend.position = "bottom")
ggsave(paste0("figures/pca/", trait, "_PC1_PC2.jpeg"), pl, width = 6, height = 6, bg = "white")

### END
