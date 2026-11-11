# Making a summary table to look at matched and unmatched strain IDs

### Mix-up table investigation.
# runs in analyses/01_rna_genotype_matching (the Snakefile's working directory)
rm(list = ls())

# Read mix-up data
mix_up <- read.table("results/ase_counter_result/rna_to_geno_summary.tsv", header = T)

# Reading meta data
metadata <- read.table("src/metadata_mouse_all_plates_with_matched_geno_ids.tsv", header = T, sep = "\t")
colnames(metadata)[1] <- "RNA_ID"

# Merging metadata and mix-up table
mix_meta <- merge(mix_up, metadata[, c(1:3,5,6)], by = "RNA_ID", sort = F)
mix_meta <- mix_meta[, -c(4,10)]
names(mix_meta) <- c("Library_ID", "Top_Simil_ID_ByMixAna", "Top_Simil_ByMixAna", "Next_Simil_ID_ByMixAna", "Next_Simil_ByMixAna", "Strain_ID_ByUCSD","Closest_Match_In_Geno", "Original_IDs_In_Geno")
write.csv(mix_meta, file = "results/ase_counter_result/mixup_table_forcheck.csv", row.names = F)

### END