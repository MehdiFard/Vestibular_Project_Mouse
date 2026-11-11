#!/usr/bin/env Rscript

# Set wd
# runs in analyses/01_rna_genotype_matching (the Snakefile's working directory)
rm(list = ls())

library(dplyr)
library(stringr)
library(DESeq2)
library(data.table)
library(ggplot2)

# Making expression matrix using raw counts
files_list <- list.files("results/count_rsem/", pattern = "*.genes.results.gz")
lib_id <- str_extract(files_list, ".*(?=\\.genes\\.results\\.gz)")

# Run the loop to assembel the expresion data from all samples
for (i in 1:length(files_list)) {
    cat("Reading sample:", lib_id[i], "- loop No. =", i, "\n")

    sample_expr <- fread(paste0("results/count_rsem/", files_list[i]))[, c("gene_id", "expected_count")]
    sample_expr$expected_count <- as.integer(sample_expr$expected_count)
    names(sample_expr)[2] <- lib_id[i]

# Merging samples
    if(i == 1) {
    expr_count = sample_expr
    } else {
    expr_count <- merge(expr_count, sample_expr, by = "gene_id", sort = F)
    }
}
# Putting the count in expr_table to keep it safe
expr_table <- expr_count

# Reading GTF file
library(rtracklayer)
gtf <- readGFF("../../resources/genome/Mus_musculus.GRCm38.102.chr.gtf")

# Getting a list of genes that are protein coding or lncRNA 
genes_protein_lncRNA <- gtf %>%
    dplyr::filter(type == "gene") %>% 
    dplyr::filter(gene_biotype %in% c("protein_coding", "lncRNA")) %>%
    dplyr::pull("gene_id")

# Filter for protein genes and lncRNAs
expr_table <- expr_table[expr_table$gene_id %in% genes_protein_lncRNA, ]
dim(expr_table)

# Metadata
rna_metadata <- read.table("src/metadata_mouse_all_batches.tsv", header = T, sep = "\t") # Reading metadata
sample_data <- data.frame(row.names = rna_metadata$library_id, rna_metadata[, c("tissue", "age", "sex", "batch")])
sample_data <- sample_data %>% mutate(sex = ifelse(is.na(sex), "no_info", sex))

# Filtering the genes for normalization
expr_table_df <-  data.frame(expr_table[, -1], row.names = expr_table$gene_id)
expr_table_filtered_4normalization <- expr_table_df %>%
    select(-FTL_775, -FTL_752) %>% # remove these two libraries (FTL_775, FTL_752) from expr table cause we have other libs for these samples and I earlier removed them from metadata
    select(all_of(row.names(sample_data))) %>% # Order by library order in metadata
    arrange(row.names(.)) %>% # Sort by gene id
    filter(rowSums(.) != 0) %>% # Remove genes with zero expression in all libraries
    filter(rowSums(. != 0) >= 50) %>% # Remove genes with fewer than 50 non-zero values
    {
        gene_sds <- apply(., 1, sd) # Calculate standard deviations for all genes
        first_qu <- quantile(gene_sds, 0.25) # Calculate the first quartile
        filter(., gene_sds > first_qu) # Keep genes with SD greater than the first quartile
    }
    
dim(expr_table_filtered_4normalization)

# Checking the IDs in sample_data ad expression file
table(row.names(sample_data) == colnames(expr_table_filtered_4normalization))

# Making DEseq2 data  
deseq2Data <- DESeqDataSetFromMatrix(countData=expr_table_filtered_4normalization, colData=sample_data,
                                    design = ~ tissue + age + sex + batch)

# Normalizing the data
deseq2Data <- estimateSizeFactors(deseq2Data)
# sizeFactors(deseq2Data)
normalized_counts <- counts(deseq2Data, normalized=TRUE)

###----------- Running PCA on Expression data
normalized_counts_log <- log2(normalized_counts + 1)
normalized_counts_final <- as.data.frame(t(normalized_counts_log))

# Running PCA
pca <- prcomp(normalized_counts_final, scale = F)
vars <- apply(pca$x, 2, var)
var.prop <- vars / sum(vars)

# scree plot
var.prop_df <- as.data.frame(var.prop)
var.prop_df$pc <- seq(1:dim(var.prop_df)[1])
p1_scree = ggplot(data =  var.prop_df[1:10, ], aes(x = pc, y=var.prop)) +
    geom_point(color="red", fill="#69b3a2", shape=21, size = 3) +
    geom_line(color="#69b3a2") +
    labs(title = "Expression PCA - Scree plot",
        x = "PC",
        y = "Proportion of variation") +
    theme_classic() +
    theme(
        axis.line = element_blank(),   # Removes the axis lines
        axis.ticks = element_blank()) +
    scale_x_continuous(breaks = 1:10); p1_scree

ggsave(filename = "01_expr_pc1_vs_pc2_DEseq2_normalized_scree.jpeg", plot =  p1_scree,
        units = "in", path = "figures/", dpi = 300)
# PC scatter plot
pcs <- as.data.frame(pca$x)
pc_levels <- merge(sample_data, pcs, by = 'row.names', sort = F)
colnames(pc_levels)[1] <- "library_id"
var.percent <- round(vars * 100 / sum(vars), 2) # Getting the percentage of variance for each PC level

# Based on the previous results these samples are mislabeled, highlighting them
mislabeled_samples <- c("FTL_600", "FTL_601", "FTL_900", "FTL_901")
# mislabeled_samples <- c("FTL_600", "FTL_601", "FTL_900", "FTL_901", "FTL_880", "FTL_881")

# PCA plot
p2_pca = ggplot(data = pc_levels, aes(x = PC1, y = PC2, color = tissue, shape = age), show.legend = T) +
    geom_point(size = 3) + # Adjusts dot size in the plot
    geom_hline(yintercept = 0, linetype = "dotted", col = "black") +
    geom_vline(xintercept = 0, linetype = "dotted", col = "black") +
    labs(
        title = "Expression PCA - PC1 vs PC2, DEseq2_normalized",
        x = paste0("PC 1 (", var.percent[1], " %)"),
        y = paste0("PC 2 (", var.percent[2], " %)"),
        color = "Tissue", shape = "Age") +
    theme_classic() +
    theme(
        axis.line = element_blank(), # Removes the axis lines
        axis.ticks = element_blank()) +
    geom_text(data = pc_levels[pc_levels$library_id %in% mislabeled_samples, ],
            aes(label = library_id), 
            nudge_y = 1, # Adjust position of labels
            nudge_x = 5,
            color = "black") +
    guides(color = guide_legend(override.aes = list(size = 2))) # Adjusts legend dot size
p2_pca

ggsave(filename = "02.3_expr_pc1_vs_pc2_DEseq2_normalized_scatter_mislabeled_tagged.jpeg", plot =  p2_pca,
        units = "cm", width = 30, height = 30, path = "figures/", dpi = 300)


###---------------- Using TPM
# Making expression matrix using raw counts
files_list <- list.files("results/count_rsem/", pattern = "*.genes.results.gz")
lib_id <- str_extract(files_list, ".*(?=\\.genes\\.results\\.gz)")

# Run the loop to assembel the expresion data from all samples
for (i in 1:length(files_list)) {
    cat("Reading sample:", lib_id[i], "- loop No. =", i, "\n")

    sample_expr <- fread(paste0("results/count_rsem/", files_list[i]))[, c("gene_id", "TPM")]
    sample_expr$TPM <- as.integer(sample_expr$TPM)
    names(sample_expr)[2] <- lib_id[i]

# Merging samples
    if(i == 1) {
    expr_count = sample_expr
    } else {
    expr_count <- merge(expr_count, sample_expr, by = "gene_id", sort = F)
    }
}
expr_table_tpm <- expr_count

# Filter for protein genes and lncRNAs
expr_table_tpm <- expr_table_tpm[expr_table_tpm$gene_id %in% genes_protein_lncRNA, ]
dim(expr_table_tpm)

# Filtering the genes for normalization
expr_table__tpm_df <-  data.frame(expr_table_tpm[, -1], row.names = expr_table_tpm$gene_id)
expr_table__tpm_filtered <- expr_table__tpm_df %>%
    select(-FTL_775, -FTL_752) %>% # remove these two libraries (FTL_775, FTL_752) from expr table cause we have other libs for these samples and I earlier removed them from metadata
    select(all_of(row.names(sample_data))) %>% # Order by library order in metadata
    arrange(row.names(.)) %>% # Sort by gene id
    filter(rowSums(.) != 0) %>% # Remove genes with zero expression in all libraries
    filter(rowSums(. != 0) >= 50) %>% # Remove genes with fewer than 50 non-zero values
    {
        gene_sds <- apply(., 1, sd) # Calculate standard deviations for all genes
        first_qu <- quantile(gene_sds, 0.25) # Calculate the first quartile
        filter(., gene_sds > first_qu) # Keep genes with SD greater than the first quartile
    }

dim(expr_table__tpm_filtered)

# Checking the IDs in sample_data ad expression file
table(row.names(sample_data) == colnames(expr_table__tpm_filtered))

# Running PCA
counts_tpm_log <- log2(expr_table__tpm_filtered + 1)
counts_tpm_final <- as.data.frame(t(counts_tpm_log))

# Running PCA
pca_tpm <- prcomp(counts_tpm_final, scale = F)
vars_tpm <- apply(pca_tpm$x, 2, var)
var.prop_tpm <- vars_tpm / sum(vars_tpm)

# PC scatter plot
pcs_tpm <- as.data.frame(pca_tpm$x)
pc_levels_tpm <- merge(sample_data, pcs_tpm, by = 'row.names', sort = F)
colnames(pc_levels_tpm)[1] <- "library_id"
var.percent_tpm <- round(vars_tpm * 100 / sum(vars_tpm), 2) # Getting the percentage of variance for each PC level

# PCA plot
p2_pca_tpm = ggplot(data = pc_levels_tpm, aes(x = PC1, y = PC2, color = tissue, shape = age), show.legend = T) +
    geom_point(size = 3) + # Adjusts dot size in the plot
    geom_hline(yintercept = 0, linetype = "dotted", col = "black") +
    geom_vline(xintercept = 0, linetype = "dotted", col = "black") +
    labs(
        title = "Expression PCA - PC1 vs PC2, TPM",
        x = paste0("PC 1 (", var.percent_tpm[1], " %)"),
        y = paste0("PC 2 (", var.percent_tpm[2], " %)"),
        color = "Tissue", shape = "Age") +
    theme_classic() +
    theme(
        axis.line = element_blank(),   # Removes the axis lines
        axis.ticks = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 2))) # Adjusts legend dot size
p2_pca_tpm

ggsave(filename = "02.2_expr_pc1_vs_pc2_tpm_scatter.jpeg", plot =  p2_pca_tpm,
        units = "in", path = "figures/", dpi = 300)

### END