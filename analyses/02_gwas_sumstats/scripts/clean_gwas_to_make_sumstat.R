#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
rm(list = ls())

# Reading variant info
vars <- fread("results/genotype/geno_all_refalt_ok.pvar")

# Reading GWAS result file name
gwas_fn <- list.files(path = "input/gwas_results/", pattern = "*.gwas$")
gwas_fn_cln <- gsub("\\.rsID\\.gwas$|\\.gwas$", "", gwas_fn)


for(i in 1:length(gwas_fn)){

    cat("Cleaning:", gwas_fn[i], "\n")

    # Reading gwas results
    gwas <- fread(paste0("input/gwas_results/", gwas_fn[i]))

    # Correct column name
    if (gwas_fn[i] %in% c("aged_vsep_f.20250228.rsID.gwas", "aged_vsep_m.20250228.rsID.gwas")) {
        colnames(gwas)[5] <- "Beta"
    }

    # Add chr_pos column
    gwas <- gwas %>%
        mutate(chr_pos = paste0(CHR, "_", BP))
    # zscore = round(sign(Beta) * abs(qnorm(P / 2)), 6)

    # Merge and making 
    gwas_cln <- vars %>%
        mutate(CHR = gsub("X", 20, `#CHROM`),
                chr_pos = paste0(CHR, "_", POS)) %>%
        inner_join(gwas, by = "chr_pos") %>%
        select(SNP = ID, A1 = ALT, A2 = REF, BETA = Beta, P)

    # Write output
    write.table(gwas_cln, file = paste0("results/gwas_results/", gwas_fn_cln[i], ".gwas.cln"), quote = F, col.names = T, row.names = F, sep = "\t")
}

# Making metadata file for twas
gwas_fn_cln_no_number <- gsub("\\.[0-9]+$", "", gwas_fn_cln)
gwas_meta <- data.frame(trait = gwas_fn_cln_no_number, sample_size = rep(100, length(gwas_fn_cln_no_number)))
write.table(gwas_meta, "results/gwas_results/gwas_metadata.tsv", quote = F, col.names = T, row.names = F, sep = "\t")

# Copy to  twas dir
file.copy("results/gwas_results/gwas_metadata.tsv", "../../analyses/03_pantry/pheast/input/gwas_sumstat/")

### END