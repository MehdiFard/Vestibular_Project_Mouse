### Makeing matadata file for single-cell RNA-seq data

rm(list = ls())

library(data.table)
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)


# List FASTQ files
fqs <- list.files("data/fastq_scRNA", full.names = TRUE)

# Extract metadata from filenames
fq_df <- tibble(path = fqs) %>%
    mutate(
        file = basename(path),
        lib_id = str_extract(file, "FTL_\\d+"),
        read = str_extract(file, "R[12]")
    ) %>%
    select(lib_id, read, path)

# Pivot wider to get R1 and R2 in columns
sample_df <- fq_df %>%
    pivot_wider(names_from = read, values_from = path) %>%
    rename(R1 = R1, R2 = R2) %>%
    arrange(lib_id)


# Read metadata from Excel file
metadata <- as.data.frame(readxl::read_excel("src/single_cell_rnaseq_metadata_mouse_b1.xlsx", sheet = 1))
metadata <- janitor::clean_names(metadata)

# Merge metadata with sample_df
metadata_fq <- metadata %>%
    select(library_id, tissue, read_length = seq_length,
            barcode_i7, barcode_i5) %>%
    mutate(
        # Standardize tissue names
        tissue = case_when(
        str_detect(tissue, regex("saccule", ignore_case = TRUE)) ~ "saccule",
        str_detect(tissue, regex("utricle", ignore_case = TRUE)) ~ "utricle",
        TRUE ~ tissue
        ),
        # Clean up library_id and read_length
        library_id = str_remove(library_id, "_snRNAseq"),
        read_length = str_remove(read_length, "PE")) %>%
        inner_join(sample_df, by = c("library_id" = "lib_id"))


# Optionally write to TSV
write.table(metadata_fq, file = "src/rna_metadata.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
