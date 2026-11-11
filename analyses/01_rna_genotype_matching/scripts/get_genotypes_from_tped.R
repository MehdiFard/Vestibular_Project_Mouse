library(tidyverse)
library(magrittr)
library(data.table)

make_geno <- function(strain_fn, out_dir, geno_file) {
  geno <- fread(geno_file, sep = "\t")
  src_strains <- colnames(geno)[-c(1:4)]

  # Create the output filenames
  dst_stem <- file.path(out_dir, str_remove(basename(strain_fn), '\\.list'))
  dst_tped <- str_glue('{dst_stem}.tped')
  dst_tfam <- str_glue('{dst_stem}.tfam')

  # Read phenotype data
  pheno <- fread(strain_fn, sep = '\t', header = TRUE)
  targets <- pheno %>%
    select(1:2) %>%
    set_names(c('FID', 'IID'))

  # Check for absent strains
  absent <- setdiff(targets$FID, src_strains)
  if (length(absent)) {
    message(paste(c('Strains not in source file:', absent), collapse = '\n'))
    targets %<>% filter(FID %in% src_strains)
  }

  # Write .tped file
  geno[, c(colnames(geno)[1:4], targets$FID), with = FALSE] %>%
    fwrite(dst_tped, sep = '\t', col.names = FALSE)

  # Write .tfam file
  tibble(FID = targets$FID %>% str_replace_all('[ /]', '.'),
         IID = targets$IID %>% str_replace_all('[ /]', '.'),
         fa = 0,
         ma = 0,
         sex = 0,
         pheno = -9) %>%
    fwrite(dst_tfam, sep = '\t', col.names = FALSE)

  # Write .pheno.txt file
  pheno %>% mutate_at(1:2, ~str_replace_all(.x, '[ /]', '.')) %>%
    fwrite(str_glue('{dst_stem}.pheno.txt'), sep = '\t')
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Please provide at least three arguments: the input list file, the output directory, and the genotype file path.")
}

input_list <- args[1]
output_dir <- args[2]
geno_file <- args[3]

make_geno(input_list, output_dir, geno_file)

### END