#!/usr/bin/env Rscript

### 05 Manhattan and QQ plots of one trait and model
# Usage (from analyses/10_gwas): Rscript scripts/05_manhattan_qq.R <trait> <noPCs|withPCs>
suppressPackageStartupMessages({ library(dplyr); library(qqman) })

args <- commandArgs(trailingOnly = TRUE)
trait <- args[1]; model <- args[2]
tag <- paste0(trait, ", FaST-LMM LOCO, ", model)
dir.create(file.path("figures", model), recursive = TRUE, showWarnings = FALSE)

gwas.res <- read.table(paste0("results/", model, "/", trait, ".gwas.loco.tsv"), header = T, sep = "\t") %>%
  select(SNP, Chr, ChrPos, PValue)
names(gwas.res) <- c("SNP", "CHR", "BP", "P")

SNP_numbers <- dim(gwas.res)[1]
chr_numbers <- length(unique(gwas.res$CHR))

Genome_Wide_tr <- 0.05/SNP_numbers
Chr_Wide_tr <- 0.05/(SNP_numbers/chr_numbers)

# Manhattan plot (chr 20 = X)
jpeg(filename = paste0("figures/", model, "/", trait, ".manhattan.jpeg"), width = 1200, height = 600, units = "px", quality = 100)
manhattan(gwas.res, chr="CHR",bp="BP",p="P",snp="SNP",col=c("dodgerblue","darkgreen"), logp = T, xlab="Chromosome",ylab="-log10 pvalue",main=tag,
          chrlabs = c(1:19, "X"), suggestiveline = -log10(Chr_Wide_tr), genomewideline = -log10(Genome_Wide_tr))
dev.off()

# QQ-Plot
p <- gwas.res$P
chisq <- qchisq(p ,1, lower.tail = F) # or qchisq(1-p ,1)
lambda <- round(median(chisq)/qchisq(0.5,1), 2)
jpeg(filename = paste0("figures/", model, "/", trait, ".QQ.jpeg"), width = 700, height = 700, units = "px", quality = 100)
qq(p,  main = paste0(tag, " - ", "GIF=", lambda))
dev.off()

### END
