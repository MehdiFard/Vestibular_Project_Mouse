#!/usr/bin/env Rscript

# Set wd
# runs in analyses/01_rna_genotype_matching (the Snakefile's working directory)
rm(list = ls())

# Load libraries
library(readxl)
library(janitor)
library(dplyr)
library(data.table)
library(stringr)

# Making a list of fastq files of all batches
fqs <- list.files("../../data/fastq_bulk", pattern = "^FTL.*fastq\\.gz$", recursive = TRUE, full.names = F)  # batch1-4/
fqs <- basename(fqs)
length(fqs)

# Removing fastq files FTL_946 to FTL_968 (from batch 3), because they're not mouse RNA data (they are some human samples for another project !!!)
fqs_to_exclude <- "946|947|948|949|950|951|952|953|954|955|956|957|958|959|960|961|962|963|964|965|966|967|968"
fqs_clean <- fqs[!grepl(fqs_to_exclude, fqs)]
fqs_clean <- gsub("_R[12]_001.fastq.gz", "", fqs_clean)
fqs_clean <- unique(fqs_clean)
length(fqs_clean)
# # Removing "FTL_1046_RNAseq_oldmouse_S30_L005". It doesn't have date
# fqs_clean <- fqs_clean[!grepl("FTL_1046_RNAseq_oldmouse_S30_L005", fqs_clean)]


# Write a list containing address to all fastqs
write.table(fqs_clean, "src/fastqs_all.list", col.names = F, row.names = F, quote = F)

###----------- Reading the metadata for all plates
# First defining some functions to rename the age and sex
# Rename the age 
rename_age <- function(x){
    x <- switch(as.character(x), "8 weeks"="w8", "22 weeks"="w22")
    return(x)
}

# Rename the sex (considering NAs)
rename_sex <- function(x){
    if (is.na(x)) {
    return(NA)
    }
    x <- switch(as.character(x), 
    "female" = "F", "male" = "M", x) # Return the original value if no match
    return(x)
}

# Batch 1
batch01 <- read_excel("src/mouse_samples_info.xlsx", sheet = 1)
batch01 <- clean_names(batch01)
batch01$age <- unlist(lapply(batch01$age, rename_age))
batch01$sex <- unlist(lapply(batch01$sex, rename_sex))
batch01 <- batch01 %>%
    dplyr::select(library_id, mouse_strain, tissue, sex, age) %>% 
    mutate(batch = rep("b1", n())) %>% 
    mutate(insert_size = rep(101, n()))

# Batch 2
batch02 <- read_excel("src/mouse_samples_info.xlsx", sheet = 2)
batch02 <- clean_names(batch02)
batch02$age <- unlist(lapply(batch02$age, rename_age))
batch02$sex <- unlist(lapply(batch02$sex, rename_sex))
batch02 <- batch02 %>%
    dplyr::select(library_id, mouse_strain, tissue, sex, age) %>% 
    mutate(batch = rep("b2", n())) %>% 
    mutate(insert_size = rep(101, n()))
head(batch02)
# View(batch02)

## Some info
# BXD78 >>> It has been sequenced 2 times for saccule in w22 and 2 times for utricle in w22
# Library IDs for saccule
# FTL_751 (Sac-w22) *** Has bigger size
# FTL_775 (Sac-w22)

# Library IDs for utricle
# FTL_752 (Utr-w22)
# FTL_776 (Utr-w22) *** Has bigger size

# I keep those fastq files with bigger size, so removing FTL_775 and FTL_752
indx3 <- which(batch02$library_id %in% c("FTL_775", "FTL_752"))
batch02 <- batch02[-indx3, ]

# Batch 3
batch03 <- read_excel("src/mouse_samples_info.xlsx", sheet = 3)
batch03 <- clean_names(batch03)
batch03$age <- unlist(lapply(batch03$age, rename_age))
batch03 <- batch03 %>%
    dplyr::select(library_id, mouse_strain, tissue, sex, age) %>% 
    mutate(batch = rep("b3", n())) %>% 
    mutate(insert_size = rep(151, n()))

# Batch 4
batch04 <- read_excel("src/mouse_samples_info.xlsx", sheet = 4)
batch04 <- clean_names(batch04)
batch04$age <- unlist(lapply(batch04$age, rename_age))
batch04 <- batch04 %>%
    dplyr::select(library_id, mouse_strain, tissue, sex, age) %>% 
    mutate(batch = rep("b4", n())) %>% 
    mutate(insert_size = rep(101, n())) %>%
    mutate(library_id = str_replace(library_id, "_RNAseq_oldmouse", ""))
print(batch04, n = dim(batch04)[1])
# View(batch04)
length(unique(batch04$mouse_strain))


# Have to remove "/RwwJ" and "/TyJ" from some strain ids to make them the same with other batches format
batch04 <- batch04 %>%
    mutate(mouse_strain = str_replace_all(mouse_strain, c("/RwwJ" = "", "/TyJ" = "")))
print(batch04, n = dim(batch04)[1])


# Binding all batches metadata to make a full metadata
meta_data_all_batches <- rbind(batch01, batch02, batch03, batch04)
head(meta_data_all_batches)
dim(meta_data_all_batches)
# View(meta_data_all_batches)

# Writing metadata file
write.table(meta_data_all_batches, "src/metadata_mouse_all_batches.tsv", row.names = F, col.names = T, quote = F, sep = "\t")

# Writing a file for getting Read Group information
write.table(meta_data_all_batches[, c("library_id", "mouse_strain", "insert_size")],
            file = "src/lib_sam_ins_4RG.tsv", col.names = F, row.names = F, quote = F, sep = "\t")


###-------------------------- Matching RNAseq data IDs with genotype IDs
# Reading raw genotype file and getting strain ids
geno <- fread("results/genotype/all_strains.tped")
strain_ids_geno <- names(geno)[-1:-4]
strain_ids_geno_df <- as.data.frame(strain_ids_geno)
dim(strain_ids_geno_df)

# Matching IDs
strain_ids_rna <- unique(meta_data_all_batches$mouse_strain); length(strain_ids_rna)
which(strain_ids_rna == "BXD174")
which(strain_ids_rna == "BXD174a")

# Matching ids in RNA-seq data with ids in genotype file
table(strain_ids_rna %in% strain_ids_geno)
strain_ids_rna[strain_ids_rna %in% strain_ids_geno]
strain_ids_rna[!strain_ids_rna %in% strain_ids_geno]

# Seems most of the strain ids in RNA-seq data don't have letters after "/". Changing strain ids from genotype data based ids in mouse RNA-seq data
pattern <- "BXD"
modified_ids_geno <- c()

for (i in 1:length(strain_ids_geno)) {
    if (grepl(pattern, strain_ids_geno[i], ignore.case = T)) { # Check if the pattern is present in the string
        modified_ids_geno[i] <- sub("/.*", "", strain_ids_geno[i]) # Remove the pattern from the string
    } else {
        modified_ids_geno[i] <- strain_ids_geno[i]
        }
}

# Changing BTBR "T<+> tf/J" to "BTBRT<+>tf/J" in genotype ids
modified_ids_geno <- stringr::str_replace_all(string =  modified_ids_geno, pattern =   c(" " = ""))

# Problem with "BXD24". Changing to original one in modified_ids_geno
indx5 <- which(modified_ids_geno == "BXD24")
modified_ids_geno[indx5[2]] <- "BXD24/TyJ-Cep290<rd16>/J"

# Checking match/unmatched ids
table(strain_ids_rna %in% modified_ids_geno)
strain_ids_rna[!strain_ids_rna %in% modified_ids_geno] # IDs from mouse RNA-seq data that are not in the genotype data

# Correcting some IDs in "strain_ids_rna" because they're wrong there
### Note: There is no close match to "BXD139a", so I assign a random ID (BXD16/TyJ) based on previous mixup analysis result.  
replacements <- c("BL6KS/J" = "C57BLKS/J",  "NOD/ShiLTJ" = "NOD/ShiLtJ", "NZB/BINJ" = "NZB/BlNJ", "KK/HlJ+A12" = "KK/HlJ",
                "C3H/H2J" = "C3H/HeJ", "BXF213" = "BXD213", "BXD131a" = "BXD131", "BXD174a" = "BXD174",
                "C57BI6J" = "C57BL/6J", "BXD162a" = "BXD162", "BXD147a" = "BXD147",
                "BXD128a" = "BXD128", "BXD183a" = "BXD183", "BXD139a" = "BXD16")


strain_ids_rna_modified <- stringr::str_replace_all(string =  strain_ids_rna, pattern =  replacements)
# The "str_replace_all" command couldn't replace "KK/HlJ+A12"! So, doing that using "which" command
strain_ids_rna_modified[which(strain_ids_rna_modified == "KK/HlJ+A12")] <- "KK/HlJ"

# Checking match/unmatched ids
table(strain_ids_rna_modified %in% modified_ids_geno)
length(unique(strain_ids_rna_modified))
strain_ids_rna_modified[duplicated(strain_ids_rna_modified)]
# "BXD174"    "C57BLKS/J" "BXD128"    "BXD213"
strain_ids_rna_modified_clean = janitor::make_clean_names(strain_ids_rna_modified)

# Need to match original and modified geno IDs
geno_ids_mod_ori = data.frame(modified_ids_geno, strain_ids_geno)
geno_id_mod <- c()
geno_id_ori <- c()
for (i in 1:length(strain_ids_rna_modified)) {
  for (j in 1:dim(geno_ids_mod_ori)[1]) {
    if(strain_ids_rna_modified[i] == geno_ids_mod_ori$modified_ids_geno[j]){
      geno_id_mod[i] <- geno_ids_mod_ori$modified_ids_geno[j]
      geno_id_ori[i] <- geno_ids_mod_ori$strain_ids_geno[j]
    }
  }
}
geno_id_ori <- gsub(" ", "", geno_id_ori) # Removing spaces between characters
geno_id_ori <- gsub("-", "_", geno_id_ori) # "-" to "_"

# Making a df with both original and modified RNA IDs
rna_ids_original_corrected <- data.frame(strain_ids_rna,
                                        strain_ids_ran_corrected_matched_with_geno = strain_ids_rna_modified,
                                        strain_ids_ran_corrected_matched_with_geno_clean = strain_ids_rna_modified_clean,
                                        geno_id_original = geno_id_ori,
                                        geno_id_modified = geno_id_mod)

# Matching and adding modified strain IDs to the Metadata
ids_matched <- c()
ids_matched_clean <- c()
geno_id_original <- c()
geno_id_modified <- c()

for (i in 1:dim(meta_data_all_batches)[1]) {
  for (j in 1:dim(rna_ids_original_corrected)[1]) {
    if(meta_data_all_batches$mouse_strain[i] == rna_ids_original_corrected$strain_ids_rna[j]) {
      ids_matched[i] <- rna_ids_original_corrected$strain_ids_ran_corrected_matched_with_geno[j]
      ids_matched_clean[i] <- rna_ids_original_corrected$strain_ids_ran_corrected_matched_with_geno_clean[j]
      geno_id_original[i] <- rna_ids_original_corrected$geno_id_original[j]
      geno_id_modified[i] <- rna_ids_original_corrected$geno_id_modified[j]
    }
  }
}

# Making final metadata
meta_data_all_plates_with_geno_ids <- data.frame(meta_data_all_batches[, 1:2], strids_matche_with_geno = ids_matched,
                                                strids_matche_with_geno_clean = ids_matched_clean,
                                                geno_id_original, geno_id_modified,
                                                meta_data_all_batches[, 3:7])
colnames(meta_data_all_plates_with_geno_ids)[2] <- "mouse_strain_provided_by_RNA_data"
write.table(meta_data_all_plates_with_geno_ids, "src/metadata_mouse_all_plates_with_matched_geno_ids.tsv", row.names = F, col.names = T, quote = F, sep = "\t")
# View(meta_data_all_plates_with_geno_ids)

# Write a look up table for getting sample id for ASE analysis
write.table(meta_data_all_plates_with_geno_ids[, c("library_id", "geno_id_original")],
            "src/lookup_table.tsv", col.names = F, row.names = F,
            quote = F, sep = "\t")

###-------------------------- Preparing genotype file
geno <- fread("results/genotype/all_strains.tped")
strain_ids_geno <- names(geno)[-1:-4]
strain_ids_geno <- gsub(" ", "", strain_ids_geno) # Removing spaces between character
strain_ids_geno <- gsub("-", "_", strain_ids_geno) # "-" to "_"

# Making .tped and .tfam files
l <- length(strain_ids_geno)
tfam <- data.frame(FID = rep(0, l), IID = strain_ids_geno, SIRE = rep(0, l),
                  DAM = rep(0, l), SEX = rep(0, l), PHEN = rep(-9, l))
# Writing .tfam file
write.table(tfam, file = "results/genotype/mouse_geno_all.tfam",
            col.names = F, row.names = F, quote = F, sep = "\t")

# Writing .tped file
write.table(geno, file = "results/genotype/mouse_geno_all.tped",
            col.names = F, row.names = F, quote = F, sep = "\t")

### END

