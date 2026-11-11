# runs in analyses/01_rna_genotype_matching (the Snakefile's working directory)
rm(list = ls())
library(data.table)
library(dplyr)
library(stringr)
library(tidyr)



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


###---------------- Subsetting Utricle - w22 - only unique ones with high similarity
saccule_w22_HighSimil <- meta_mixup %>%
    filter(tissue == "Saccule" & age == "w22") %>%
    filter(Top_Simil_ByMixAna >= 0.9) %>%
    mutate(FID = rep(0, n())) %>%
    mutate(IID_tmp = seq(1, n())) %>% 
    dplyr::select(library_id, Top_Simil_ID_ByMixAna, FID, IID_tmp, batch, insert_size)
dim(saccule_w22_HighSimil)
View(saccule_w22_HighSimil)

# Check
length(unique(saccule_w22_HighSimil$Top_Simil_ID_ByMixAna))
saccule_w22_HighSimil$Top_Simil_ID_ByMixAna[duplicated(saccule_w22_HighSimil$Top_Simil_ID_ByMixAna)]
# View(saccule_w22_HighSimil)
# View(meta_mixup)

#- I need to unique the df for each tissue-age, so just don't want to remove one of the duplicates in random,
# want to remove that one that duplicated because of replacing with mix-up results

# Filter duplicated strains in saccule_w22_HighSimil
dups_sac_w22_HighSimil <- saccule_w22_HighSimil %>%
    filter(Top_Simil_ID_ByMixAna %in% Top_Simil_ID_ByMixAna[duplicated(Top_Simil_ID_ByMixAna)])
tissue <- "Saccule"
age <- "w22"

# Create the dups_check_sac_w22_HighSimil dataframe with "YES" or "NO"
dups_check_sac_w22_HighSimil <- dups_sac_w22_HighSimil %>%
    mutate(Original_ID = ifelse(
    library_id %in% meta_mixup$library_id &
    meta_mixup$tissue[match(library_id, meta_mixup$library_id)] == tissue &
    meta_mixup$age[match(library_id, meta_mixup$library_id)] == age &
    meta_mixup$Original_IDs_In_Geno[match(library_id, meta_mixup$library_id)] == Top_Simil_ID_ByMixAna,
    "YES", "NO"
    )) %>%
    select(Strain = Top_Simil_ID_ByMixAna, Library_ID = library_id, Original_ID)

# Filter saccule_w22_HighSimil to keep only rows with "YES" in Original_ID
saccule_w22_unik_HighSimil <- saccule_w22_HighSimil %>%
    filter(!(library_id %in% dups_check_sac_w22_HighSimil$Library_ID[dups_check_sac_w22_HighSimil$Original_ID == "NO"]))

# Check
dim(saccule_w22_unik_HighSimil)
length(unique(saccule_w22_unik_HighSimil$Top_Simil_ID_ByMixAna))
View(saccule_w22_unik_HighSimil)
saccule_w22_unik_HighSimil$Top_Simil_ID_ByMixAna[duplicated(saccule_w22_unik_HighSimil$Top_Simil_ID_ByMixAna)]

# We still have some duplicated samples because of being sequenced two times one in batch 3 and on in batch 4. I remove one of them in random.
# "BXD79/RwwJ"  "BXD213/RwwJ"
saccule_w22_unik_HighSimil <- saccule_w22_unik_HighSimil %>%
    distinct(Top_Simil_ID_ByMixAna, .keep_all = T) %>%
    mutate(IID_tmp = row_number())
dim(saccule_w22_unik_HighSimil)
View(saccule_w22_unik_HighSimil)

# saccule_w22_HighSimil sample IDs to subset in plink
write.table(saccule_w22_unik_HighSimil[, c("Top_Simil_ID_ByMixAna", "IID_tmp")], "src/saccule_w22_similarity0.9_unik.list", quote = F, col.names = T, row.names = F, sep = "\t")  

# Due to repeted samples, I have to use the following script to subset genotype for saccule_w22.
# In bash - order of inputs and outputs: list of samples, output dir, input genotype.
system("Rscript scripts/get_genotypes_from_tped.R src/saccule_w22_similarity0.9_unik.list results/genotype/ results/genotype/all_strains.tped")

# Reading tfam for adding library IDs
tfam.sac_w22_unik_HighSimil <- fread("results/genotype/saccule_w22_similarity0.9_unik.tfam")
tfam.sac_w22_unik_HighSimil$V1 <- stringr::str_replace_all(tfam.sac_w22_unik_HighSimil$V1, pattern =  c("/" = ".", " " = "."))
saccule_w22_unik_HighSimil$Top_Simil_ID_ByMixAna <- stringr::str_replace_all(saccule_w22_unik_HighSimil$Top_Simil_ID_ByMixAna, pattern =  c("/" = ".", " " = "."))

# Matching Ids
matching_indices_sac_w22_unik_HighSimil <- match(saccule_w22_unik_HighSimil$Top_Simil_ID_ByMixAna, tfam.sac_w22_unik_HighSimil$V1)
tfam.libID.sac_w22_unik_HighSimil  <- data.frame(V1 = saccule_w22_unik_HighSimil$Top_Simil_ID_ByMixAna,
                                libID = saccule_w22_unik_HighSimil$library_id,
                                V2_V6 = tfam.sac_w22_unik_HighSimil[, 2:6],
                                geno_ID = tfam.sac_w22_unik_HighSimil$V1[matching_indices_sac_w22_unik_HighSimil])
# Check
table(tfam.sac_w22_unik_HighSimil$V1 == tfam.libID.sac_w22_unik_HighSimil$geno_ID)
table(tfam.sac_w22_unik_HighSimil$V1 == saccule_w22_unik_HighSimil$Top_Simil_ID_ByMixAna)
table(tfam.libID.sac_w22_unik_HighSimil$geno_ID == saccule_w22_unik_HighSimil$Top_Simil_ID_ByMixAna)
table(paste0(tfam.libID.sac_w22_unik_HighSimil$libID, "_", tfam.libID.sac_w22_unik_HighSimil$geno_ID) == paste0(saccule_w22_unik_HighSimil$library_id, "_", saccule_w22_unik_HighSimil$Top_Simil_ID_ByMixAna))

# Making tfam with library IDs, i'm over writting the tfam file.
tfam.libID.sac_w22_unik_HighSimil <- tfam.libID.sac_w22_unik_HighSimil[, c(4,2,4:7)]
write.table(tfam.libID.sac_w22_unik_HighSimil, "results/genotype/saccule_w22_similarity0.9_unik.tfam", quote = F, col.names = F, row.names = F, sep = "\t") 

# Making the genotype
system("plink2 --tfile results/genotype/saccule_w22_similarity0.9_unik --maf 0.05 --geno 0.05 --make-pgen --not-chr X Y --snps-only 'just-acgt' --sort-vars --nonfounders --out results/pantry_inputs/saccule_w22_sorted")
system("mkdir -p ../../analyses/03_pantry/pheast/input/saccule_w22")
system("plink2 --pfile results/pantry_inputs/saccule_w22_sorted --make-bed --nonfounders --out ../../analyses/03_pantry/pheast/input/saccule_w22/saccule_w22")


# ### Sub-setting expression data
# expr_bed <- fread("intermediate/count_table/expr_moues.iqn.filtered.bed")
# dim(expr_bed)

# ### Subset BED file - saccule - w22
# expr_bed_saccule_w22 <- expr_bed %>% 
#     dplyr::select("#chr", start, end, gene_id, tfam.libID.sac_w22_unik_HighSimil$libID)
# dim(expr_bed_saccule_w22)

# # Write out the bed file - Utricle - w22
# write.table(expr_bed_saccule_w22, "intermediate/count_table/saccule_w22.iqn.bed",
#             quote = F, col.names = T, row.names = F, sep = "\t")
# system("bgzip intermediate/count_table/saccule_w22.iqn.bed")
# system("tabix -p bed intermediate/count_table/saccule_w22.iqn.bed.gz")



### Making fastq_map
lib_ids_to_sub <- tfam.libID.sac_w22_unik_HighSimil$libID
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
#     con <- gzfile(fq, "rt")
#     readLines(con, n = 2) -> lines
#     close(con)
#     nchar(lines[2])
# }

# df_out <- df_out %>%
#     rowwise() %>%
#     mutate(read_length = get_read_length(R1)) %>%
#     ungroup() %>%
#     mutate(across(everything(), ~ gsub("../../data/fastq_bulk/batch[0-9]+", "input/saccule_w22/fastq", .)))


# Or, can get the isert size from "saccule_w22_unik_HighSimil"
df_out <- saccule_w22_unik_HighSimil %>%
    select(libID = library_id, insert_size) %>%
    right_join(df_out, by = "libID") %>%
    select(R1, R2, libID, insert_size)
View(df_out)

# Write to file
# Pantry reads the FASTQs through links in its input/saccule_w22/fastq/ (to the files in data/fastq_bulk/batch<n>/)
pantry_fq <- "../../analyses/03_pantry/phenotyping/input/saccule_w22/fastq"
dir.create(pantry_fq, showWarnings = FALSE, recursive = TRUE)
for (f in c(df_out$R1, df_out$R2)) {
    link <- file.path(pantry_fq, basename(f))
    if (!file.exists(link)) file.symlink(normalizePath(f), link)
}
df_out <- df_out %>% mutate(across(c(R1, R2), ~ file.path("input/saccule_w22/fastq", basename(.))))

write.table(df_out, "results/pantry_inputs/saccule_w22.fastq_map.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(df_out$libID, "results/pantry_inputs/saccule_w22.samples.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# copy to pantry
system("cp results/pantry_inputs/saccule_w22.fastq_map.txt ../../analyses/03_pantry/phenotyping/input/saccule_w22/")
system("cp results/pantry_inputs/saccule_w22.samples.txt ../../analyses/03_pantry/phenotyping/input/saccule_w22/")


# extract extra info to include in covariate file
batch_for_covar <- saccule_w22_unik_HighSimil %>%
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
write.table(batch_for_covar[, c("library_id", "batch_numeric")], "results/pantry_inputs/saccule_w22.batches.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
system("cp results/pantry_inputs/saccule_w22.batches.tsv ../../analyses/03_pantry/pheast/input/saccule_w22/")
