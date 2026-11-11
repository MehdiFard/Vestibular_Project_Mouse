# 12 - Assemble one master metadata file for all 260 libraries.
#
# Pulls together everything scattered across the project into a single table:
# identity, genetics (with the provenance of every strain call), technical
# details, sequencing QC, expression QC, design flags and PCA coordinates.
#
# Writes results/tables/master_metadata.tsv (plus a .csv copy) and a small
# data-dictionary file describing every column.

source("../../analyses/05_dge/scripts/00_config.R")

library(data.table)

# ------------------------------------------------------------------ base
md <- fread(file.path(DATA, "sample_table.tsv"))
setnames(md, "strain_meta", "strain_metadata_sheet", skip_absent = TRUE)

# ------------------------------------------------------------------ genetics
# Carry the runner-up strain as well, so every call can be second-guessed.
gm <- fread(GENO_MATCH)
gm <- unique(gm[, .(library_id = RNA_ID,
                    strain_runnerup = next_simil_ID,
                    strain_runnerup_simil = next_simil)])
md <- merge(md, gm, by = "library_id", all.x = TRUE, sort = FALSE)

# The original RNA-submission label and the raw genotype ID, for traceability.
ref <- fread(REF_META)
ref <- unique(ref[, .(library_id,
                      strain_label_from_rna_submission = mouse_strain_provided_by_RNA_data,
                      geno_id_original,
                      insert_size)])
md <- merge(md, ref, by = "library_id", all.x = TRUE, sort = FALSE)

# ------------------------------------------------------------------ fastq
fq <- rbindlist(lapply(DATASETS, function(ds) {
  f <- fread(file.path(PANTRY, "phenotyping", "input", ds,
                       paste0(ds, ".fastq_map.txt")), header = FALSE)
  data.table(library_id = f$V3,
             fastq_r1 = basename(f$V1), fastq_r2 = basename(f$V2),
             read_length = f$V4)
}))
fq <- fq[, .(fastq_r1 = paste(unique(fastq_r1), collapse = ";"),
             fastq_r2 = paste(unique(fastq_r2), collapse = ";"),
             read_length = first(read_length),
             n_fastq_pairs = .N), by = library_id]
md <- merge(md, fq, by = "library_id", all.x = TRUE, sort = FALSE)
md[, paired_end := TRUE]

# ------------------------------------------------------------------ kallisto QC
ri <- rbindlist(lapply(seq_len(nrow(md)), function(i) {
  f <- file.path(PANTRY, "phenotyping", "intermediate", md$dataset[i],
                 "expression", md$library_id[i], "run_info.json")
  if (!file.exists(f)) return(NULL)
  j <- readLines(f, warn = FALSE)
  num <- function(key) {
    ln <- grep(paste0('"', key, '"'), j, value = TRUE)
    if (!length(ln)) return(NA_real_)
    as.numeric(sub('.*:\\s*([0-9.]+).*', '\\1', ln[1]))
  }
  data.table(library_id = md$library_id[i],
             n_reads_processed = num("n_processed"),
             n_pseudoaligned   = num("n_pseudoaligned"),
             pct_pseudoaligned = num("p_pseudoaligned"),
             pct_unique        = num("p_unique"))
}))
md <- merge(md, ri, by = "library_id", all.x = TRUE, sort = FALSE)

# ------------------------------------------------------------------ expression QC
txi  <- readRDS(file.path(DATA, "txi_gene.rds"))
anno <- fread(file.path(DATA, "gene_annotation.tsv"))
cnt  <- txi$counts[, md$library_id]
keep <- rownames(cnt) %in% anno[biotype %in% KEEP_BIOTYPES, gene_id]

md[, total_gene_counts     := round(colSums(cnt))]
md[, genes_detected_min10  := colSums(cnt >= 10)]
md[, genes_detected_min10_pc_linc := colSums(cnt[keep, ] >= 10)]

# ------------------------------------------------------------------ design flags
# pair_id = animal from 01_build_inputs.R (submission workbook; strain + age + sex where the workbook has no pair)
md[, n_organs_for_pair := .N, by = pair_id]
md[, n_ages_for_strain_organ := uniqueN(age), by = .(strain, organ)]
md[, sampled_at_both_ages := n_ages_for_strain_organ == 2]

md[, in_paired_organ_contrast := has_organ_pair]
md[, in_paired_age_contrast   := sampled_at_both_ages]

# ------------------------------------------------------------------ PCA
pcs <- fread(file.path(TABLES, "pca_scores.tsv"))
md <- merge(md, pcs[, .(library_id, PC1, PC2, PC3, PC4, PC5)],
            by = "library_id", all.x = TRUE, sort = FALSE)

# ------------------------------------------------------------------ tidy + write
setcolorder(md, c(
  # identity
  "library_id", "dataset", "organ", "age", "group",
  # genetics
  "strain", "strain_source", "strain_simil", "strain_gap",
  "strain_runnerup", "strain_runnerup_simil",
  "strain_metadata_sheet", "strain_label_from_rna_submission", "geno_id_original",
  # technical
  "batch", "sex_final", "sex_source", "sex", "insert_size", "read_length", "paired_end",
  "n_fastq_pairs", "fastq_r1", "fastq_r2",
  # sequencing QC
  "n_reads_processed", "n_pseudoaligned", "pct_pseudoaligned", "pct_unique",
  # expression QC
  "total_gene_counts", "genes_detected_min10", "genes_detected_min10_pc_linc",
  # design
  "pair_id", "pair_source", "has_organ_pair", "sampled_at_both_ages",
  "in_paired_organ_contrast", "in_paired_age_contrast",
  # coordinates
  "PC1", "PC2", "PC3", "PC4", "PC5"))
md[, c("n_organs_for_pair", "n_ages_for_strain_organ") := NULL]
setorder(md, organ, age, library_id)

stopifnot(nrow(md) == 260, uniqueN(md$library_id) == 260, !any(is.na(md$strain)))

fwrite(md, file.path(TABLES, "master_metadata.tsv"), sep = "\t", na = "NA")
fwrite(md, file.path(TABLES, "master_metadata.csv"), na = "NA")

# ------------------------------------------------------------------ dictionary
dict <- data.table(
  column = names(md),
  description = c(
    "Sequencing library ID, the key used everywhere in this project",
    "Which of the four Pantry datasets this library belongs to",
    "Saccule or Utricle - taken from the dataset name",
    "w8 or w22 - taken from the dataset name",
    "Organ and age combined, for convenience",
    "Mouse strain, from RNA-to-genotype concordance (the authoritative label)",
    "How the strain was assigned, and whether it agrees with the metadata sheet",
    "Concordance with the best-matching strain (0-1); build threshold was 0.90",
    "Concordance minus the runner-up's; a large gap means an unambiguous call",
    "Second-best matching strain",
    "Concordance with that runner-up",
    "Strain as recorded in metadata_mouse_all_plates_with_matched_geno_ids.csv (NA for the 24 later libraries)",
    "Strain as written on the original RNA submission sheet",
    "Raw genotype sample ID before name cleaning",
    "Sequencing plate, 1-4 (batch 4 is the later 'oldmouse' run)",
    "Sex, complete for all 260 - from expression (Xist vs Y genes) where the metadata had none",
    "How sex was determined and whether it agrees with the metadata",
    "Sex as recorded in the metadata sheet; NA for 52 libraries",
    "Library insert size",
    "Read length in bp",
    "Whether the library is paired-end",
    "Number of FASTQ pairs (lanes) for this library",
    "R1 FASTQ filename(s), semicolon-separated if more than one lane",
    "R2 FASTQ filename(s)",
    "Reads processed by kallisto",
    "Reads pseudoaligned by kallisto",
    "Percent pseudoaligned - the main mapping-rate QC metric",
    "Percent uniquely pseudoaligned",
    "Total gene-level counts after tximport",
    "Genes with at least 10 counts in this library, all biotypes",
    "Genes with at least 10 counts, protein-coding and lincRNA only",
    "Animal ID: from the submission workbook (consecutive library IDs), else strain_age_sex; its two organs form a pair",
    "Where pair_id comes from: submission workbook, or strain + age + sex",
    "TRUE if both organs of this animal are present",
    "TRUE if this strain/organ was sampled at both 8 and 22 weeks",
    "TRUE if the library is used in the paired organ contrast (analyses A and B)",
    "TRUE if the library is used in the within-strain age contrast",
    "Expression PC1 (VST, top 2000 variable genes) - tracks organ",
    "Expression PC2", "Expression PC3 - tracks age", "Expression PC4 - tracks strain",
    "Expression PC5"))
stopifnot(nrow(dict) == ncol(md))
fwrite(dict, file.path(TABLES, "master_metadata_dictionary.tsv"), sep = "\t")

cat("Wrote master_metadata.tsv:", nrow(md), "libraries x", ncol(md), "columns\n\n")
cat("Strain provenance:\n"); print(table(md$strain_source))
cat("\nMapping rate (percent pseudoaligned):\n"); print(summary(md$pct_pseudoaligned))
cat("\nLibraries in each contrast:\n")
cat("  paired organ contrast:", sum(md$in_paired_organ_contrast), "\n")
cat("  within-strain age contrast:", sum(md$in_paired_age_contrast), "\n")
