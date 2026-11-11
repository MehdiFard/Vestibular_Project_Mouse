#!/usr/bin/env Rscript

### 01 Phenotype files: one per trait from the young-mice VsEP + beam test workbook (sheets 1-4, one row per mouse)
# Out: intermediate/phenotype/<trait>.phen (header; FID = strain as in the genotype file, IID = <strain>_<n>, trait value)
# Usage (from analyses/10_gwas): Rscript scripts/01_make_phenotypes.R input/<workbook.xlsx>
suppressPackageStartupMessages({ library(readxl); library(dplyr); library(stringr) })

xlsx <- commandArgs(trailingOnly = TRUE)[1]
trait <- c("threshold", "latency", "amplitude", "beamtest")
sheets <- c("VsEP_Threshold_young mice", "P2_latency", "P1_amplitude", "Beam_test")   # the order the traits rely on
stopifnot(identical(excel_sheets(xlsx)[1:4], sheets))

dir.create("intermediate/phenotype", recursive = TRUE, showWarnings = FALSE)
for (i in 1:4) {
  phen <- read_excel(xlsx, sheet = i)
  names(phen) <- c("FID", "jax", "n", trait[i])
  rep_str <- c("/" = ".", " " = ".") # Order of replacement character is important and shouldn't be changed.
  phen <- phen %>%
    filter(!is.na(.data[[trait[i]]])) %>%
    mutate(IID = paste0(str_replace_all(FID, rep_str), "_", n)) %>%
    select(FID, IID, all_of(trait[i]))
  stopifnot(!anyDuplicated(phen$IID), !anyNA(phen$IID))
  write.table(phen, file = paste0("intermediate/phenotype/", trait[i], ".phen"), col.names = T, row.names = F, quote = F, sep = "\t")
  message(trait[i], ": ", nrow(phen), " mice, ", n_distinct(phen$FID), " strains")
}

### END
