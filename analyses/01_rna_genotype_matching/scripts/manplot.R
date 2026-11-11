library(data.table)
library(qqman)
library(gridExtra)
library(dplyr)

# Set wd
setwd("../../data/gwas_results/")

# List of the files
fn <- list.files()

for(i in 1:length(fn)){
# read gwas result
gwas <- fread(fn[i])

# set thr
n_snp = nrow(gwas)
thr_bonf = 0.05/n_snp

# numeric chr column
gwas <- gwas %>%
    mutate(
        CHR = ifelse(CHR == "X", 20, CHR),
        CHR = as.numeric(CHR)
    )

# Set output file
pdf(paste0(fn[i], ".plot.pdf"), width = 10, height = 8)

# Create Manhattan plot
manhattan(gwas, chr = "CHR", bp = "BP", p = "P", snp = "SNP", 
            main = paste0(fn[i]), 
            col = c("dodgerblue4", "darkorange3"))

# Create Q-Q plot
qq(gwas$P, main = paste0(fn[i]))

# Close PDF device
dev.off()
}

# Threshold FDR
gwas$FDR <- p.adjust(gwas$P, method = "BH")
signif_snps <- gwas[gwas$FDR < 0.05, ]
nrow(signif_snps)  # number of significant SNPs
p_cutoff <- max(gwas$P[gwas$FDR < 0.05])
-log10(p_cutoff)  # for plotting
sum(gwas$FDR < 0.05)


# scp <user>@<login_node>:/path/to/data/gwas_results/*.pdf .