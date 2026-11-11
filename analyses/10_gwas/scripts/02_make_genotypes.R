#!/usr/bin/env Rscript

### 02 Genotypes of the phenotyped mice: each mouse gets its strain's genotype (inbred strains)
# Out: intermediate/genotype/<trait>.{tped,tfam}, intermediate/phenotype/<trait>.pheno.txt (FaST-LMM: FID IID value, no header)
# Usage (from analyses/10_gwas): Rscript scripts/02_make_genotypes.R intermediate/phenotype/<trait>.phen input/all_strains_mm10.tped
suppressPackageStartupMessages({ library(dplyr); library(stringr); library(tibble) })

args <- commandArgs(trailingOnly = TRUE)
strain_fn <- args[1]
geno <- data.table::fread(args[2], sep = '\t', header = T)   # snp_chr snp_id centimorgans snp_bp_mm10 <strain> ...
src_strains <- colnames(geno)[-c(1:4)]

trait <- str_remove(basename(strain_fn), '\\.phen$')
dst_stem <- file.path("intermediate/genotype", trait)
dir.create("intermediate/genotype", recursive = TRUE, showWarnings = FALSE)

pheno <- data.table::fread(strain_fn, sep = '\t', header = T)
targets <- pheno %>%
  select(1:2) %>%
  setNames(c('FID', 'IID'))

absent <- setdiff(targets$FID, src_strains)
if (length(absent)) {
  message(paste(c('Strains not in source file:', absent), collapse = '\n'))
  targets <- targets %>% filter(FID %in% src_strains)
  pheno <- pheno %>% filter(.[[1]] %in% src_strains)
}

# write .tped file (one genotype column per mouse)
geno[, c(colnames(geno)[1:4], targets$FID), with = F] %>%
  data.table::fwrite(str_glue('{dst_stem}.tped'), sep = '\t', col.names = F)

# write .tfam file
tibble(FID = targets$FID %>% str_replace_all('[ /]', '.'),
       IID = targets$IID %>% str_replace_all('[ /]', '.'),
       fa = 0,
       ma = 0,
       sex = 0,
       pheno = -9) %>%
  data.table::fwrite(str_glue('{dst_stem}.tfam'), sep = '\t', col.names = F)

# write .pheno.txt file
pheno %>% mutate(across(1:2, ~str_replace_all(.x, '[ /]', '.'))) %>%
  data.table::fwrite(str_glue('intermediate/phenotype/{trait}.pheno.txt'), col.names = FALSE, sep = '\t')

message(trait, ": ", nrow(targets), " mice, ", nrow(geno), " SNPs")

### END
