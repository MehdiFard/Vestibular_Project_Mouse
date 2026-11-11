#!/usr/bin/env Rscript
# Preparing genotype data for eQTL analysis - sub-setting samples for Saccule and Utricle

# runs in analyses/01_rna_genotype_matching (the Snakefile's working directory)
rm(list = ls())

library(data.table)
library(dplyr)

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


###---------------- Subsetting Saccule
saccule <- meta_mixup %>%
    filter(tissue == 'Saccule') %>%
    mutate(FID = rep(0, n())) %>%
    mutate(IID_tmp = seq(1, n())) %>% 
    dplyr::select(library_id, Top_Simil_ID_ByMixAna, FID, IID_tmp)

dim(saccule)
length(unique(saccule$Top_Simil_ID_ByMixAna))

# A table of saccule sample IDs to subset in plink
write.table(saccule[, c("Top_Simil_ID_ByMixAna", "IID_tmp")], "src/saccule.list", quote = F, col.names = T, row.names = F, sep = "\t")  

# Subsetting samples for Saccule tissue (it contains both w8 and w22 data for Saccule tissue)
# Due to repeated samples, I have to use the following script to subset genotype for saccule.
# In bash - order of inputs and outputs: list of samples, output dir, input genotype.
system("Rscript scripts/get_genotypes_from_tped.R src/saccule.list results/genotype/ results/genotype/all_strains.tped")

# Reading tfam for adding library IDs
tfam.sub.sac <- fread("results/genotype/saccule.tfam")
tfam.sub.sac$V1 <- stringr::str_replace_all(tfam.sub.sac$V1, pattern =  c("/" = ".", " " = "."))
saccule$Top_Simil_ID_ByMixAna <- stringr::str_replace_all(saccule$Top_Simil_ID_ByMixAna, pattern =  c("/" = ".", " " = "."))

# Matching Ids
matching_indices_sac <- match(saccule$Top_Simil_ID_ByMixAna, tfam.sub.sac$V1)
tfam.libID.sac  <- data.frame(V1 = saccule$Top_Simil_ID_ByMixAna,
                              libID = saccule$library_id,
                              V2_V6 = tfam.sub.sac[, 2:6],
                              geno_ID = tfam.sub.sac$V1[matching_indices_sac])
# Check
table(tfam.sub.sac$V1 == tfam.libID.sac$geno_ID)
table(tfam.sub.sac$V1 == saccule$Top_Simil_ID_ByMixAna)
table(tfam.libID.sac$geno_ID == saccule$Top_Simil_ID_ByMixAna)
table(paste0(tfam.libID.sac$libID, "_", tfam.libID.sac$geno_ID) == paste0(saccule$library_id, "_", saccule$Top_Simil_ID_ByMixAna))

# Making tfam with library IDs.
tfam.libID.sac <- tfam.libID.sac[, c(4,2,4:7)]
write.table(tfam.libID.sac, "results/genotype/saccule_geno_libid.tfam", quote = F, col.names = F, row.names = F, sep = "\t") 

# Copy .tped file with the base name of "saccule_geno_libid"
file.copy("results/genotype/saccule.tped", "results/genotype/saccule_geno_libid.tped")



###---------------- Sub-setting Utricle
utricle <- meta_mixup %>%
  filter(tissue == 'Utricle') %>%
  mutate(FID = rep(0, n())) %>%
  mutate(IID_tmp = seq(1, n())) %>% 
  dplyr::select(library_id, Top_Simil_ID_ByMixAna, FID, IID_tmp)

# Some check
dim(utricle)
length(unique(utricle$Top_Simil_ID_ByMixAna))
table(unique(utricle$Top_Simil_ID_ByMixAna) %in% unique(saccule$Top_Simil_ID_ByMixAna))
table(unique(saccule$Top_Simil_ID_ByMixAna) %in% unique(utricle$Top_Simil_ID_ByMixAna))

# A table of utricle sample IDs to subset in plink
write.table(utricle[, c("Top_Simil_ID_ByMixAna", "IID_tmp")], "src/utricle.list", quote = F, col.names = T, row.names = F, sep = "\t")  

# Due to repeated samples, I have to use the following script to subset genotype for utricle.
# In bash - order of inputs and outputs: list of samples, output dir, input genotype.
system("Rscript scripts/get_genotypes_from_tped.R src/utricle.list results/genotype/ results/genotype/all_strains.tped")

# Reading tfam for adding library IDs
tfam.sub.utr <- fread("results/genotype/utricle.tfam")
tfam.sub.utr$V1 <- stringr::str_replace_all(tfam.sub.utr$V1, pattern =  c("/" = ".", " " = "."))
utricle$Top_Simil_ID_ByMixAna <- stringr::str_replace_all(utricle$Top_Simil_ID_ByMixAna, pattern =  c("/" = ".", " " = "."))

# Matching IDs
matching_indices_utr <- match(utricle$Top_Simil_ID_ByMixAna, tfam.sub.utr$V1)
tfam.libID.utr  <- data.frame(V1 = utricle$Top_Simil_ID_ByMixAna,
                              libID = utricle$library_id,
                              V2_V6 = tfam.sub.utr[, 2:6],
                              geno_ID = tfam.sub.utr$V1[matching_indices_utr])
# Check
table(tfam.sub.utr$V1 == tfam.libID.utr$geno_ID)
table(tfam.sub.utr$V1 == utricle$Top_Simil_ID_ByMixAna)
table(tfam.libID.utr$geno_ID == utricle$Top_Simil_ID_ByMixAna)
table(paste0(tfam.libID.utr$libID, "_", tfam.libID.utr$geno_ID) == paste0(utricle$library_id, "_", utricle$Top_Simil_ID_ByMixAna))

# Making tfam with library IDs.
tfam.libID.utr <- tfam.libID.utr[, c(4,2,4:7)]
write.table(tfam.libID.utr, "results/genotype/utricle_geno_libid.tfam", quote = F, col.names = F, row.names = F, sep = "\t") 

# Copy .tped file with the base name of "utricle_geno_libid"
if(file.copy("results/genotype/utricle.tped", "results/genotype/utricle_geno_libid.tped")) {
  cat("The '.tped' file has been renamed and copied properly!")
  } else {
      cat("The copy wasn't successful!", "\n")
}


###---------------- Subsetting Saccule - w8 - all (not unique)
saccule_w8 <- meta_mixup %>%
  filter(tissue == 'Saccule' & age == "w8") %>%
  # filter(Top_Simil_ByMixAna >= 0.9) %>%
  mutate(FID = rep(0, n())) %>%
  mutate(IID_tmp = seq(1, n())) %>% 
  dplyr::select(library_id, Top_Simil_ID_ByMixAna, FID, IID_tmp) # %>%
  # distinct(Top_Simil_ID_ByMixAna, .keep_all = TRUE)

dim(saccule_w8)
length(unique(saccule_w8$Top_Simil_ID_ByMixAna))
# table(saccule_w8$Top_Simil_ID_ByMixAna)
saccule_w8$Top_Simil_ID_ByMixAna[duplicated(saccule_w8$Top_Simil_ID_ByMixAna)]
View(saccule_w8)
View(meta_mixup)

dups <- saccule_w8$Top_Simil_ID_ByMixAna[duplicated(saccule_w8$Top_Simil_ID_ByMixAna)]

which(dups[1] %in% saccule_w8$Top_Simil_ID_ByMixAna)


# saccule_w8 sample IDs to subset in plink
write.table(saccule_w8[, c("Top_Simil_ID_ByMixAna", "IID_tmp")], "src/saccule_w8.list", quote = F, col.names = T, row.names = F, sep = "\t")  

# Due to repeted samples, I have to use the following script to subset genotype for saccule_w8.
# In bash - order of inputs and outputs: list of samples, output dir, input genotype.
system("Rscript scripts/get_genotypes_from_tped.R src/saccule_w8.list results/genotype/ results/genotype/all_strains.tped")

# Reading tfam for adding library IDs
tfam.sub.sac_w8 <- fread("results/genotype/saccule_w8.tfam")
tfam.sub.sac_w8$V1 <- stringr::str_replace_all(tfam.sub.sac_w8$V1, pattern =  c("/" = ".", " " = "."))
saccule_w8$Top_Simil_ID_ByMixAna <- stringr::str_replace_all(saccule_w8$Top_Simil_ID_ByMixAna, pattern =  c("/" = ".", " " = "."))

# Matching Ids
matching_indices_sac <- match(saccule_w8$Top_Simil_ID_ByMixAna, tfam.sub.sac_w8$V1)
tfam.libID.sac_w8  <- data.frame(V1 = saccule_w8$Top_Simil_ID_ByMixAna,
                              libID = saccule_w8$library_id,
                              V2_V6 = tfam.sub.sac_w8[, 2:6],
                              geno_ID = tfam.sub.sac_w8$V1[matching_indices_sac])
# Check
table(tfam.sub.sac_w8$V1 == tfam.libID.sac_w8$geno_ID)
table(tfam.sub.sac_w8$V1 == saccule_w8$Top_Simil_ID_ByMixAna)
table(tfam.libID.sac_w8$geno_ID == saccule_w8$Top_Simil_ID_ByMixAna)
table(paste0(tfam.libID.sac_w8$libID, "_", tfam.libID.sac_w8$geno_ID) == paste0(saccule_w8$library_id, "_", saccule_w8$Top_Simil_ID_ByMixAna))

# Making tfam with library IDs.
tfam.libID.sac_w8 <- tfam.libID.sac_w8[, c(4,2,4:7)]
write.table(tfam.libID.sac_w8, "results/genotype/saccule_w8_withdups.tfam", quote = F, col.names = F, row.names = F, sep = "\t") 


# Copy .tped file with the base name of "saccule_w8_geno_libid"
if(file.copy("results/genotype/saccule_w8.tped", "results/genotype/saccule_w8_withdups.tped")) {
  cat("The '.tped' file has been renamed and copied properly!", "\n")
  } else {
      cat("The copy wasn't successful!", "\n")
}


###---------------- Subsetting Saccule - w8 - only unique ones
saccule_w8 <- meta_mixup %>%
  filter(tissue == 'Saccule' & age == "w8") %>%
  # filter(Top_Simil_ByMixAna >= 0.9) %>%
  mutate(FID = rep(0, n())) %>%
  mutate(IID_tmp = seq(1, n())) %>% 
  dplyr::select(library_id, Top_Simil_ID_ByMixAna, FID, IID_tmp) # %>%
  # distinct(Top_Simil_ID_ByMixAna, .keep_all = TRUE)

dim(saccule_w8)
length(unique(saccule_w8$Top_Simil_ID_ByMixAna))
# table(saccule_w8$Top_Simil_ID_ByMixAna)
saccule_w8$Top_Simil_ID_ByMixAna[duplicated(saccule_w8$Top_Simil_ID_ByMixAna)]
View(saccule_w8)
View(meta_mixup)

#- I need to unique the df for each tissue-age, so juat dont want to remove one of the duplicates in random,
# want to remove that one that duplicated because of replacing with mix-up results

# Filter duplicated strains in saccule_w8
dups <- saccule_w8 %>%
  filter(Top_Simil_ID_ByMixAna %in% Top_Simil_ID_ByMixAna[duplicated(Top_Simil_ID_ByMixAna)])
tissue <- "Saccule"
age <- "w8"

# Create the result dataframe with "YES" or "NO"
result <- dups %>%
  mutate(Original_ID = ifelse(
    library_id %in% meta_mixup$library_id &
    meta_mixup$tissue[match(library_id, meta_mixup$library_id)] == tissue &
    meta_mixup$age[match(library_id, meta_mixup$library_id)] == age &
    meta_mixup$Original_IDs_In_Geno[match(library_id, meta_mixup$library_id)] == Top_Simil_ID_ByMixAna,
    "YES", "NO"
  )) %>%
  select(Strain = Top_Simil_ID_ByMixAna, Library_ID = library_id, Original_ID)

# Filter saccule_w8 to keep only rows with "YES" in Original_ID
saccule_w8_unik <- saccule_w8 %>%
  filter(!(library_id %in% result$Library_ID[result$Original_ID == "NO"]))

dim(saccule_w8_unik)
length(unique(saccule_w8_unik$Top_Simil_ID_ByMixAna))
View(saccule_w8_unik)


# saccule_w8 sample IDs to subset in plink
write.table(saccule_w8_unik[, c("Top_Simil_ID_ByMixAna", "IID_tmp")], "src/saccule_w8_unik.list", quote = F, col.names = T, row.names = F, sep = "\t")  

# Due to repeted samples, I have to use the following script to subset genotype for saccule_w8.
# In bash - order of inputs and outputs: list of samples, output dir, input genotype.
system("Rscript scripts/get_genotypes_from_tped.R src/saccule_w8_unik.list results/genotype/ results/genotype/all_strains.tped")

# Reading tfam for adding library IDs
tfam.sac_w8_unik <- fread("results/genotype/saccule_w8_unik.tfam")
tfam.sac_w8_unik$V1 <- stringr::str_replace_all(tfam.sac_w8_unik$V1, pattern =  c("/" = ".", " " = "."))
saccule_w8_unik$Top_Simil_ID_ByMixAna <- stringr::str_replace_all(saccule_w8_unik$Top_Simil_ID_ByMixAna, pattern =  c("/" = ".", " " = "."))

# Matching Ids
matching_indices_sacw8_unik <- match(saccule_w8_unik$Top_Simil_ID_ByMixAna, tfam.sac_w8_unik$V1)
tfam.libID.sac_w8_unik  <- data.frame(V1 = saccule_w8_unik$Top_Simil_ID_ByMixAna,
                              libID = saccule_w8_unik$library_id,
                              V2_V6 = tfam.sac_w8_unik[, 2:6],
                              geno_ID = tfam.sac_w8_unik$V1[matching_indices_sacw8_unik])
# Check
table(tfam.sac_w8_unik$V1 == tfam.libID.sac_w8_unik$geno_ID)
table(tfam.sac_w8_unik$V1 == saccule_w8_unik$Top_Simil_ID_ByMixAna)
table(tfam.libID.sac_w8_unik$geno_ID == saccule_w8_unik$Top_Simil_ID_ByMixAna)
table(paste0(tfam.libID.sac_w8_unik$libID, "_", tfam.libID.sac_w8_unik$geno_ID) == paste0(saccule_w8_unik$library_id, "_", saccule_w8_unik$Top_Simil_ID_ByMixAna))

# Making tfam with library IDs, i'm over writting the tfam file.
tfam.libID.sac_w8_unik <- tfam.libID.sac_w8_unik[, c(4,2,4:7)]
write.table(tfam.libID.sac_w8_unik, "results/genotype/saccule_w8_unik.tfam", quote = F, col.names = F, row.names = F, sep = "\t") 






###---------------- Subsetting Saccule - w22
saccule_w22 <- meta_mixup %>%
  filter(tissue == 'Saccule' & age == "w22") %>%
  filter(Top_Simil_ByMixAna >= 0.9) %>%
  mutate(FID = rep(0, n())) %>%
  mutate(IID_tmp = seq(1, n())) %>% 
  dplyr::select(library_id, Top_Simil_ID_ByMixAna, FID, IID_tmp)

# Check
dim(saccule_w22)
length(unique(saccule_w22$Top_Simil_ID_ByMixAna))
# table(saccule_w22$Top_Simil_ID_ByMixAna)
saccule_w22$Top_Simil_ID_ByMixAna[duplicated(saccule_w22$Top_Simil_ID_ByMixAna)]
View(saccule_w22)
View(meta_mixup)


# saccule_w22 sample IDs to subset in plink
write.table(saccule_w22[, c("Top_Simil_ID_ByMixAna", "IID_tmp")], "src/saccule_w22.list", quote = F, col.names = T, row.names = F, sep = "\t")  

# Due to repeted samples, I have to use the following script to subset genotype for saccule_w22.
# In bash - order of inputs and outputs: list of samples, output dir, input genotype.
system("Rscript scripts/get_genotypes_from_tped.R src/saccule_w22.list results/genotype/ results/genotype/all_strains.tped")

# Reading tfam for adding library IDs
tfam.sub.sac_w22 <- fread("results/genotype/saccule_w22.tfam")
tfam.sub.sac_w22$V1 <- stringr::str_replace_all(tfam.sub.sac_w22$V1, pattern =  c("/" = ".", " " = "."))
saccule_w22$Top_Simil_ID_ByMixAna <- stringr::str_replace_all(saccule_w22$Top_Simil_ID_ByMixAna, pattern =  c("/" = ".", " " = "."))

# Matching Ids
matching_indices_sac <- match(saccule_w22$Top_Simil_ID_ByMixAna, tfam.sub.sac_w22$V1)
tfam.libID.sac_w22  <- data.frame(V1 = saccule_w22$Top_Simil_ID_ByMixAna,
                              libID = saccule_w22$library_id,
                              V2_V6 = tfam.sub.sac_w22[, 2:6],
                              geno_ID = tfam.sub.sac_w22$V1[matching_indices_sac])
# Check
table(tfam.sub.sac_w22$V1 == tfam.libID.sac_w22$geno_ID)
table(tfam.sub.sac_w22$V1 == saccule_w22$Top_Simil_ID_ByMixAna)
table(tfam.libID.sac_w22$geno_ID == saccule_w22$Top_Simil_ID_ByMixAna)
table(paste0(tfam.libID.sac_w22$libID, "_", tfam.libID.sac_w22$geno_ID) == paste0(saccule_w22$library_id, "_", saccule_w22$Top_Simil_ID_ByMixAna))

# Making tfam with library IDs.
tfam.libID.sac_w22 <- tfam.libID.sac_w22[, c(4,2,4:7)]
write.table(tfam.libID.sac_w22, "results/genotype/saccule_w22_geno_libid.tfam", quote = F, col.names = F, row.names = F, sep = "\t") 

# Copy .tped file with the base name of "saccule_w22_geno_libid"
if(file.copy("results/genotype/saccule_w22.tped", "results/genotype/saccule_w22_geno_libid.tped")) {
  cat("The '.tped' file has been renamed and copied properly!")
  } else {
      cat("The copy wasn't successful!")
}



###---------------- Subsetting Utricle - W8
utricle_w8 <- meta_mixup %>%
  filter(tissue == "Utricle" & age == "w8") %>%
  filter(Top_Simil_ByMixAna >= 0.9) %>%
  mutate(FID = rep(0, n())) %>%
  mutate(IID_tmp = seq(1, n())) %>% 
  dplyr::select(library_id, Top_Simil_ID_ByMixAna, FID, IID_tmp)

# Check
dim(utricle_w8)
length(unique(utricle_w8$Top_Simil_ID_ByMixAna))
# table(utricle_w8$Top_Simil_ID_ByMixAna)
utricle_w8$Top_Simil_ID_ByMixAna[duplicated(utricle_w8$Top_Simil_ID_ByMixAna)]
View(utricle_w8)
View(meta_mixup)

# utricle_w8 sample IDs to subset in plink
write.table(utricle_w8[, c("Top_Simil_ID_ByMixAna", "IID_tmp")], "src/utricle_w8.list", quote = F, col.names = T, row.names = F, sep = "\t")  

# Due to repeted samples, I have to use the following script to subset genotype for utricle_w8.
# In bash - order of inputs and outputs: list of samples, output dir, input genotype.
system("Rscript scripts/get_genotypes_from_tped.R src/utricle_w8.list results/genotype/ results/genotype/all_strains.tped")

# Reading tfam for adding library IDs
tfam.sub.utr_w8 <- fread("results/genotype/utricle_w8.tfam")
tfam.sub.utr_w8$V1 <- stringr::str_replace_all(tfam.sub.utr_w8$V1, pattern =  c("/" = ".", " " = "."))
utricle_w8$Top_Simil_ID_ByMixAna <- stringr::str_replace_all(utricle_w8$Top_Simil_ID_ByMixAna, pattern =  c("/" = ".", " " = "."))

# Matching IDs
matching_indices_utr_w8 <- match(utricle_w8$Top_Simil_ID_ByMixAna, tfam.sub.utr_w8$V1)
tfam.libID.utr_w8  <- data.frame(V1 = utricle_w8$Top_Simil_ID_ByMixAna,
                              libID = utricle_w8$library_id,
                              V2_V6 = tfam.sub.utr_w8[, 2:6],
                              geno_ID = tfam.sub.utr_w8$V1[matching_indices_utr_w8])
# Check
table(tfam.sub.utr_w8$V1 == tfam.libID.utr_w8$geno_ID)
table(tfam.sub.utr_w8$V1 == utricle_w8$Top_Simil_ID_ByMixAna)
table(tfam.libID.utr_w8$geno_ID == utricle_w8$Top_Simil_ID_ByMixAna)
table(paste0(tfam.libID.utr_w8$libID, "_", tfam.libID.utr_w8$geno_ID) == paste0(utricle_w8$library_id, "_", utricle_w8$Top_Simil_ID_ByMixAna))

# Making tfam with library IDs.
tfam.libID.utr_w8 <- tfam.libID.utr_w8[, c(4,2,4:7)]
write.table(tfam.libID.utr_w8, "results/genotype/utricle_w8_geno_libid.tfam", quote = F, col.names = F, row.names = F, sep = "\t") 


# Copy .tped file with the base name of "utricle_w8_geno_libid"
if(file.copy("results/genotype/utricle_w8.tped", "results/genotype/utricle_w8_geno_libid.tped")) {
  cat("The '.tped' file has been renamed and copied properly!")
  } else {
      cat("The copy wasn't successful!")
}


###---------------- Subsetting Utricle - W22
utricle_w22 <- meta_mixup %>%
  filter(tissue == "Utricle" & age == "w22") %>%
  filter(Top_Simil_ByMixAna >= 0.9) %>%
  mutate(FID = rep(0, n())) %>%
  mutate(IID_tmp = seq(1, n())) %>% 
  dplyr::select(library_id, Top_Simil_ID_ByMixAna, FID, IID_tmp)

# Check
dim(utricle_w22)
length(unique(utricle_w22$Top_Simil_ID_ByMixAna))
utricle_w22$Top_Simil_ID_ByMixAna[duplicated(utricle_w22$Top_Simil_ID_ByMixAna)]
View(utricle_w22)
View(meta_mixup)


# A table of utricle_w22 sample IDs to subset in plink
write.table(utricle_w22[, c("Top_Simil_ID_ByMixAna", "IID_tmp")], "src/utricle_w22.list", quote = F, col.names = T, row.names = F, sep = "\t")  

# Due to repeted samples, I have to use the following script to subset genotype for utricle_w22.
# In bash - order of inputs and outputs: list of samples, output dir, input genotype.
system("Rscript scripts/get_genotypes_from_tped.R src/utricle_w22.list results/genotype/ results/genotype/all_strains.tped")

# Reading tfam for adding library IDs
tfam.sub.utr_w22 <- fread("results/genotype/utricle_w22.tfam")
tfam.sub.utr_w22$V1 <- stringr::str_replace_all(tfam.sub.utr_w22$V1, pattern =  c("/" = ".", " " = "."))
utricle_w22$Top_Simil_ID_ByMixAna <- stringr::str_replace_all(utricle_w22$Top_Simil_ID_ByMixAna, pattern =  c("/" = ".", " " = "."))

# Matching IDs
matching_indices_utr_w22 <- match(utricle_w22$Top_Simil_ID_ByMixAna, tfam.sub.utr_w22$V1)
tfam.libID.utr_w22  <- data.frame(V1 = utricle_w22$Top_Simil_ID_ByMixAna,
                                  libID = utricle_w22$library_id,
                                  V2_V6 = tfam.sub.utr_w22[, 2:6],
                                  geno_ID = tfam.sub.utr_w22$V1[matching_indices_utr_w22])
# Check
table(tfam.sub.utr_w22$V1 == tfam.libID.utr_w22$geno_ID)
table(tfam.sub.utr_w22$V1 == utricle_w22$Top_Simil_ID_ByMixAna)
table(tfam.libID.utr_w22$geno_ID == utricle_w22$Top_Simil_ID_ByMixAna)
table(paste0(tfam.libID.utr_w22$libID, "_", tfam.libID.utr_w22$geno_ID) == paste0(utricle_w22$library_id, "_", utricle_w22$Top_Simil_ID_ByMixAna))

# Making tfam with library IDs.
tfam.libID.utr_w22 <- tfam.libID.utr_w22[, c(4,2,4:7)]
write.table(tfam.libID.utr_w22, "results/genotype/utricle_w22_geno_libid.tfam", quote = F, col.names = F, row.names = F, sep = "\t") 

# Copy .tped file with the base name of "utricle_w22_geno_libid"
if(file.copy("results/genotype/utricle_w22.tped", "results/genotype/utricle_w22_geno_libid.tped")) {
  cat("The '.tped' file has been renamed and copied properly!")
  } else {
      cat("The copy wasn't successful!")
}








###-------------------------- Sub-setting expression data for Saccule and Utricle
expr_bed <- fread("count_table/expr.qn.filtered.bed")

### Saccule
expr_bed_saccule <- expr_bed %>% 
  dplyr::select("#chr", start, end, gene_id, tfam.libID.sac$libID)

# Write out the bed file - Saccule
write.table(expr_bed_saccule, "count_table/expr.iqn.filtered.saccule.bed",
            quote = F, col.names = T, row.names = F, sep = "\t")

### Utricle
expr_bed_utricle <- expr_bed %>% 
  dplyr::select("#chr", start, end, gene_id, tfam.libID.utr$libID)

# Write out the bed file - Utricle
write.table(expr_bed_utricle, "count_table/expr.iqn.filtered.utricle.bed",
            quote = F, col.names = T, row.names = F, sep = "\t")




### Subset BED file - Saccule - w8
expr_bed_saccule_w8 <- expr_bed %>% 
  dplyr::select("#chr", start, end, gene_id, tfam.libID.sac_w8$libID)

# Write out the bed file - Saccule -w8
write.table(expr_bed_saccule_w8, "count_table/expr.iqn.filtered.saccule_w8.bed",
            quote = F, col.names = T, row.names = F, sep = "\t")



### Subset BED file - Saccule - w22
expr_bed_saccule_w22 <- expr_bed %>% 
  dplyr::select("#chr", start, end, gene_id, tfam.libID.sac_w22$libID)

# Write out the bed file - Saccule - w22
write.table(expr_bed_saccule_w22, "count_table/expr.iqn.filtered.saccule_w22.bed",
            quote = F, col.names = T, row.names = F, sep = "\t")




### Subset BED file - Utricle - w8
expr_bed_utricle_w8 <- expr_bed %>% 
  dplyr::select("#chr", start, end, gene_id, tfam.libID.utr_w8$libID)

# Write out the bed file - Utricle - w8
write.table(expr_bed_utricle_w8, "count_table/expr.iqn.filtered.utricle_w8.bed",
            quote = F, col.names = T, row.names = F, sep = "\t")


### Subset BED file - Utricle - w22
expr_bed_utricle_w22 <- expr_bed %>% 
  dplyr::select("#chr", start, end, gene_id, tfam.libID.utr_w22$libID)

# Write out the bed file - Utricle - w22
write.table(expr_bed_utricle_w22, "count_table/expr.iqn.filtered.utricle_w22.bed",
            quote = F, col.names = T, row.names = F, sep = "\t")



### END