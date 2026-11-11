# runs in analyses/01_rna_genotype_matching (the Snakefile's working directory)
rm(list = ls())
library(data.table)
library(dplyr)
library(stringr)


# Reading ASE analysis results (mixup analysis)
mixup_result <- read.csv("results/ase_counter_result/mixup_table_forcheck.csv")
colnames(mixup_result)[1] <- "library_id"

# Reading metadata for all plates
meta_all <- read.table("src/metadata_mouse_all_batches.tsv", header = T, sep = "\t")

# Based on the PCA analysis, tissue labels for FTL_900 and FTL_901 are not right. Here I change them. Original Labels = Saccule, Clustered with: Utricle
indx1.1 <- which(meta_all$library_id %in% c("FTL_900", "FTL_901"))
meta_all[indx1.1, "tissue"] <- "Utricle"

# Based on the PCA analysis, tissue label for FTL_600 and FTL_601 are swapped. Here I'm correcting them.
indx1.2 <- which(meta_all$library_id == "FTL_600")
indx1.3 <- which(meta_all$library_id == "FTL_601")
meta_all[indx1.2, "tissue"] <- "Saccule"
meta_all[indx1.3, "tissue"] <- "Utricle"

# Merging metadata with mixup analysis result
meta_mixup <- merge(meta_all, mixup_result[, c(1:3,6,8)], by = "library_id", sort = F)

# I have to replace "BTBRT<+>tf/J" with "BTBR T<+> tf/J" to be matched with IDs in genotype file
indx2 <- which(meta_mixup$Top_Simil_ID_ByMixAna == "BTBRT<+>tf/J")
meta_mixup[indx2, "Top_Simil_ID_ByMixAna"] <- "BTBR T<+> tf/J"
View(meta_mixup)


###---------------- Subsetting Utricle - W8 - only unique ones with high similarity
utricle_w8_HighSimil <- meta_mixup %>%
    filter(tissue == "Utricle" & age == "w8") %>%
    filter(Top_Simil_ByMixAna >= 0.9) %>%
    mutate(FID = rep(0, n())) %>%
    mutate(IID_tmp = seq(1, n())) %>% 
    dplyr::select(library_id, Top_Simil_ID_ByMixAna, FID, IID_tmp, batch, insert_size)
dim(utricle_w8_HighSimil)
View(utricle_w8_HighSimil)

# Check
length(unique(utricle_w8_HighSimil$Top_Simil_ID_ByMixAna))
utricle_w8_HighSimil$Top_Simil_ID_ByMixAna[duplicated(utricle_w8_HighSimil$Top_Simil_ID_ByMixAna)]
# View(utricle_w8_HighSimil)
# View(meta_mixup)

#- I need to unique the df for each tissue-age, so just don't want to remove one of the duplicates in random,
# want to remove that one that duplicated because of replacing with mix-up results

# Filter duplicated strains in utricle_w8_HighSimil
dups_utr_w8_HighSimil <- utricle_w8_HighSimil %>%
    filter(Top_Simil_ID_ByMixAna %in% Top_Simil_ID_ByMixAna[duplicated(Top_Simil_ID_ByMixAna)])
tissue <- "Utricle"
age <- "w8"

# Create the dups_check_utr_w8_HighSimil dataframe with "YES" or "NO"
dups_check_utr_w8_HighSimil <- dups_utr_w8_HighSimil %>%
    mutate(Original_ID = ifelse(
    library_id %in% meta_mixup$library_id &
    meta_mixup$tissue[match(library_id, meta_mixup$library_id)] == tissue &
    meta_mixup$age[match(library_id, meta_mixup$library_id)] == age &
    meta_mixup$Original_IDs_In_Geno[match(library_id, meta_mixup$library_id)] == Top_Simil_ID_ByMixAna,
    "YES", "NO"
    )) %>%
    select(Strain = Top_Simil_ID_ByMixAna, Library_ID = library_id, Original_ID)

# Filter utricle_w8_HighSimil to keep only rows with "YES" in Original_ID
utricle_w8_unik_HighSimil <- utricle_w8_HighSimil %>%
    filter(!(library_id %in% dups_check_utr_w8_HighSimil$Library_ID[dups_check_utr_w8_HighSimil$Original_ID == "NO"]))

# Check
dim(utricle_w8_unik_HighSimil)
length(unique(utricle_w8_unik_HighSimil$Top_Simil_ID_ByMixAna))
View(utricle_w8_unik_HighSimil)
utricle_w8_unik_HighSimil$Top_Simil_ID_ByMixAna[duplicated(utricle_w8_unik_HighSimil$Top_Simil_ID_ByMixAna)]

# We still have some duplicated samples because of being sequenced two times one in batch 3 and on in batch 4. I remove one of them in random.
# "BXD51/RwwJ"  "BXD71/RwwJ"  "BXD79/RwwJ"  "BXD213/RwwJ" "BXD128/RwwJ" "BXD183/RwwJ"
utricle_w8_unik_HighSimil <- utricle_w8_unik_HighSimil %>%
    distinct(Top_Simil_ID_ByMixAna, .keep_all = T) %>%
    mutate(IID_tmp = row_number())
dim(utricle_w8_unik_HighSimil)
View(utricle_w8_unik_HighSimil)

# utricle_w8_HighSimil sample IDs to subset in plink
write.table(utricle_w8_unik_HighSimil[, c("Top_Simil_ID_ByMixAna", "IID_tmp")], "src/utricle_w8_similarity0.9_unik.list", quote = F, col.names = T, row.names = F, sep = "\t")  

# Due to repeted samples, I have to use the following script to subset genotype for saccule_w8.
# In bash - order of inputs and outputs: list of samples, output dir, input genotype.
system("Rscript scripts/get_genotypes_from_tped.R src/utricle_w8_similarity0.9_unik.list results/genotype/ results/genotype/all_strains.tped")

# Reading tfam for adding library IDs
tfam.utr_w8_unik_HighSimil <- fread("results/genotype/utricle_w8_similarity0.9_unik.tfam")
tfam.utr_w8_unik_HighSimil$V1 <- stringr::str_replace_all(tfam.utr_w8_unik_HighSimil$V1, pattern =  c("/" = ".", " " = "."))
utricle_w8_unik_HighSimil$Top_Simil_ID_ByMixAna <- stringr::str_replace_all(utricle_w8_unik_HighSimil$Top_Simil_ID_ByMixAna, pattern =  c("/" = ".", " " = "."))

# Matching Ids
matching_indices_utr_w8_unik_HighSimil <- match(utricle_w8_unik_HighSimil$Top_Simil_ID_ByMixAna, tfam.utr_w8_unik_HighSimil$V1)
tfam.libID.utr_w8_unik_HighSimil  <- data.frame(V1 = utricle_w8_unik_HighSimil$Top_Simil_ID_ByMixAna,
                                libID = utricle_w8_unik_HighSimil$library_id,
                                V2_V6 = tfam.utr_w8_unik_HighSimil[, 2:6],
                                geno_ID = tfam.utr_w8_unik_HighSimil$V1[matching_indices_utr_w8_unik_HighSimil])
# Check
table(tfam.utr_w8_unik_HighSimil$V1 == tfam.libID.utr_w8_unik_HighSimil$geno_ID)
table(tfam.utr_w8_unik_HighSimil$V1 == utricle_w8_unik_HighSimil$Top_Simil_ID_ByMixAna)
table(tfam.libID.utr_w8_unik_HighSimil$geno_ID == utricle_w8_unik_HighSimil$Top_Simil_ID_ByMixAna)
table(paste0(tfam.libID.utr_w8_unik_HighSimil$libID, "_", tfam.libID.utr_w8_unik_HighSimil$geno_ID) == paste0(utricle_w8_unik_HighSimil$library_id, "_", utricle_w8_unik_HighSimil$Top_Simil_ID_ByMixAna))

# Making tfam with library IDs, i'm over writting the tfam file.
tfam.libID.utr_w8_unik_HighSimil <- tfam.libID.utr_w8_unik_HighSimil[, c(4,2,4:7)]
write.table(tfam.libID.utr_w8_unik_HighSimil, "results/genotype/utricle_w8_similarity0.9_unik.tfam", quote = F, col.names = F, row.names = F, sep = "\t") 



# Making the genotype
system("plink2 --tfile results/genotype/utricle_w8_similarity0.9_unik --maf 0.05 --geno 0.05 --make-pgen --not-chr X Y --snps-only 'just-acgt' --sort-vars --nonfounders --out results/pantry_inputs/utricle_w8_sorted")
system("mkdir -p ../../analyses/03_pantry/pheast/input/utricle_w8")
system("plink2 --pfile results/pantry_inputs/utricle_w8_sorted --make-bed --nonfounders --out ../../analyses/03_pantry/pheast/input/utricle_w8/utricle_w8")



# ### Sub-setting expression data
# expr_bed <- fread("intermediate/count_table/expr_moues.iqn.filtered.bed")
# dim(expr_bed)

# ### Subset BED file - Utricle - w8
# expr_bed_utricle_w8 <- expr_bed %>% 
#     dplyr::select("#chr", start, end, gene_id, tfam.libID.utr_w8_unik_HighSimil$libID)
# dim(expr_bed_utricle_w8)

# # Write out the bed file - Utricle - w8
# write.table(expr_bed_utricle_w8, "intermediate/count_table/utricle_w8.iqn.bed",
#             quote = F, col.names = T, row.names = F, sep = "\t")
# system("bgzip intermediate/count_table/utricle_w8.iqn.bed")
# system("tabix -p bed intermediate/count_table/utricle_w8.iqn.bed.gz")



### Making fastq_map

library(dplyr)
library(stringr)
library(tidyr)

lib_ids_to_sub <- tfam.libID.utr_w8_unik_HighSimil$libID
pattern_lin_ides <- paste0("^(", paste(lib_ids_to_sub, collapse = "|"), ")")

# List all matching fastq files
fastas <- list.files(
  path = "../../data/fastq_bulk",   # batch1-4/
  recursive = TRUE,
  pattern = pattern_lin_ides,
  full.names = TRUE
)
length(fastas)

# Create output in 3 columns: R1, R2, libID
df_out <- tibble(path = fastas) %>%
  mutate(
    libID = str_extract(basename(path), "^FTL_\\d+"),
    read  = if_else(str_detect(path, "_R1_"), "R1", "R2")
  ) %>%
  pivot_wider(names_from = read, values_from = path) %>%
  select(R1, R2, libID)  # final column order

# # Function to get read length from gzipped FASTQ (first read only)
# get_read_length <- function(fq) {
#   con <- gzfile(fq, "rt")
#   readLines(con, n = 2) -> lines
#   close(con)
#   nchar(lines[2])
# }

# df_out <- df_out %>%
#   rowwise() %>%
#   mutate(read_length = get_read_length(R1)) %>%
#   ungroup()


# Or, can get the isert size from "utricle_w8_unik_HighSimil"
df_out <- utricle_w8_unik_HighSimil %>%
  select(libID = library_id, insert_size) %>%
  right_join(df_out, by = "libID") %>%
  select(R1, R2, libID, insert_size)
View(df_out)

# Write to file
# Pantry reads the FASTQs through links in its input/utricle_w8/fastq/ (to the files in data/fastq_bulk/batch<n>/)
pantry_fq <- "../../analyses/03_pantry/phenotyping/input/utricle_w8/fastq"
dir.create(pantry_fq, showWarnings = FALSE, recursive = TRUE)
for (f in c(df_out$R1, df_out$R2)) {
    link <- file.path(pantry_fq, basename(f))
    if (!file.exists(link)) file.symlink(normalizePath(f), link)
}
df_out <- df_out %>% mutate(across(c(R1, R2), ~ file.path("input/utricle_w8/fastq", basename(.))))

write.table(df_out, "results/pantry_inputs/utricle_w8.fastq_map.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(df_out$libID, "results/pantry_inputs/utricle_w8.samples.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

system("cp results/pantry_inputs/utricle_w8.fastq_map.txt ../../analyses/03_pantry/phenotyping/input/utricle_w8/")
system("cp results/pantry_inputs/utricle_w8.samples.txt ../../analyses/03_pantry/phenotyping/input/utricle_w8/")


# extract extra info to include in covariate file
batch_for_covar <- utricle_w8_unik_HighSimil %>%
  select(library_id, batch) %>%
  mutate(batch_numeric = case_when(
    batch == "b1" ~ 1,
    batch == "b2" ~ 2,
    batch == "b3" ~ 3,
    batch == "b4" ~ 4,
    TRUE ~ NA_real_
  ))
head(batch_for_covar)

# Write
write.table(batch_for_covar[, c("library_id", "batch_numeric")], "results/pantry_inputs/utricle_w8.batches.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
system("cp results/pantry_inputs/utricle_w8.batches.txt ../../analyses/03_pantry/pheast/input/utricle_w8/")
