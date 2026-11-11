# 01 - Build the gene-level count matrix and the sample table.
#
# Counts come from the kallisto runs that Pantry already did
# (phenotyping/intermediate/<dataset>/expression/<sample>/abundance.tsv).
# tximport sums transcript est_counts to genes and carries the average
# transcript length along as an offset, which is what DESeq2 wants.
#
# Organ and age come from which Pantry dataset a sample belongs to.
# Strain, sex and plate come from the old project's metadata file.

source("../../analyses/05_dge/scripts/00_config.R")

library(data.table)
library(tximport)

# ---------------------------------------------------------------- gene annotation
# Pull gene_id -> gene_name / biotype / chr out of the GTF once and cache it.
anno_file <- file.path(DATA, "gene_annotation.tsv")

if (!file.exists(anno_file)) {
  cat("Parsing gene records from the GTF (one-off, takes a minute)...\n")
  awk <- sprintf(
    "awk -F'\\t' '$3==\"gene\"' %s | sed -E 's/.*gene_id \"([^\"]+)\".*gene_name \"([^\"]+)\".*gene_biotype \"([^\"]+)\".*/\\1\\t\\2\\t\\3/'",
    shQuote(GTF))
  gene_lines <- fread(cmd = awk, header = FALSE, sep = "\t",
                      col.names = c("gene_id", "gene_name", "biotype"))
  chr <- fread(cmd = sprintf("awk -F'\\t' '$3==\"gene\"{print $1}' %s", shQuote(GTF)), header = FALSE)
  gene_lines$chr <- chr$V1
  fwrite(unique(gene_lines), anno_file, sep = "\t")
}
gene_anno <- fread(anno_file)
cat("Gene annotation:", nrow(gene_anno), "genes\n")

# ---------------------------------------------------------------- tx2gene
# Built from the GTF transcript records so it matches the kallisto index.
t2g_file <- file.path(DATA, "tx2gene.tsv")

if (!file.exists(t2g_file)) {
  cat("Building tx2gene from the GTF...\n")
  awk <- sprintf(
    "awk -F'\\t' '$3==\"transcript\"' %s | sed -E 's/.*gene_id \"([^\"]+)\".*transcript_id \"([^\"]+)\".*/\\2\\t\\1/'",
    shQuote(GTF))
  t2g <- fread(cmd = awk, header = FALSE, sep = "\t", col.names = c("tx_id", "gene_id"))
  fwrite(unique(t2g), t2g_file, sep = "\t")
}
tx2gene <- fread(t2g_file)
cat("tx2gene:", nrow(tx2gene), "transcripts\n")

# ---------------------------------------------------------------- sample table
# Organ + age are taken from the dataset name, which is the definitive source.
samples <- rbindlist(lapply(DATASETS, function(ds) {
  ids <- readLines(file.path(PANTRY, "phenotyping", "input", ds, paste0(ds, ".samples.txt")))
  ids <- trimws(ids[nzchar(trimws(ids))])
  data.table(library_id = ids,
             dataset    = ds,
             organ      = ifelse(grepl("^saccule", ds), "Saccule", "Utricle"),
             age        = ifelse(grepl("w8$", ds), "w8", "w22"))
}))

# Batch comes from the batches.tsv files already prepared for Pheast (plate 1-4).
batches <- rbindlist(lapply(DATASETS, function(ds) {
  f <- file.path(PANTRY, "pheast", "input", ds, paste0(ds, ".batches.tsv"))
  b <- fread(f, header = FALSE, col.names = c("library_id", "batch"))
  b
}))
samples <- merge(samples, batches, by = "library_id", all.x = TRUE, sort = FALSE)

# Sex from the old metadata. That file stops at FTL_945, so the 24 later
# "oldmouse" w22 libraries are absent from it entirely.
meta <- fread(REF_META)
meta <- unique(meta[, .(library_id, strain_meta = geno_id_modified, sex)])
samples <- merge(samples, meta, by = "library_id", all.x = TRUE, sort = FALSE)
samples[is.na(strain_meta) | strain_meta == "", strain_meta := NA_character_]
samples[sex == "" | sex == "NA", sex := NA_character_]

# Strain comes from the RNA-to-genotype concordance table rather than the
# metadata spreadsheet, for three reasons:
#   1. it covers all 260 libraries, including the 24 missing from the metadata
#   2. it is what the Pantry pipeline itself used to attach genotypes, so the
#      strain label and the genotype behind each library agree by construction
#   3. it is 99.1% self-consistent - two organs from one animal, matched
#      independently, land on the same strain 112 times out of 113
# top_simil is the concordance with the best-matching strain and next_simil the
# runner-up; both are kept so any call can be audited.
gm <- fread(GENO_MATCH)
gm <- unique(gm[, .(library_id = RNA_ID, strain = top_simil_ID,
                    strain_simil = top_simil, strain_runnerup = next_simil)])
gm[, strain_gap := strain_simil - strain_runnerup]
samples <- merge(samples, gm[, .(library_id, strain, strain_simil, strain_gap)],
                 by = "library_id", all.x = TRUE, sort = FALSE)

# Flag where the genotype match and the old spreadsheet disagree. These are
# label swaps the concordance analysis was built to catch, not match failures.
norm_strain <- function(x) gsub("[^A-Z0-9]", "", sub("/.*$", "", toupper(trimws(x))))
samples[, strain_source := fifelse(
  is.na(strain_meta), "genotype match (absent from metadata)",
  fifelse(norm_strain(strain_meta) == norm_strain(strain),
          "genotype match (agrees with metadata)",
          "genotype match (DISAGREES with metadata)"))]

cat("\nStrain source:\n"); print(table(samples$strain_source))
cat("\nMatch confidence: min similarity", round(min(samples$strain_simil), 3),
    " min gap to runner-up", round(min(samples$strain_gap), 3), "\n")

samples[, group := factor(paste(organ, age), levels = names(GROUP_COLS))]
samples[, organ := factor(organ, levels = c("Saccule", "Utricle"))]
samples[, age   := factor(age,   levels = c("w8", "w22"))]
samples[, batch := factor(batch)]

setkey(samples, library_id)

cat("\nSample table:", nrow(samples), "libraries\n")
print(table(samples$organ, samples$age))
cat("\nBatch by age:\n"); print(table(samples$batch, samples$age))
cat("\nSex by age (NA = not in the old metadata file):\n")
print(table(samples$sex, samples$age, useNA = "ifany"))
cat("\nLibraries with a strain:", sum(!is.na(samples$strain)), "of", nrow(samples), "\n")
for (a in c("w8", "w22"))
  cat("  saccule/utricle pairs at", a, ":",
      nrow(samples[age == a & !is.na(strain), .N, by = strain][N == 2]), "\n")

fwrite(samples, file.path(DATA, "sample_table.tsv"), sep = "\t", na = "NA")

# ---------------------------------------------------------------- import counts
files <- file.path(PANTRY, "phenotyping", "intermediate", samples$dataset,
                   "expression", samples$library_id, "abundance.tsv")
names(files) <- samples$library_id
stopifnot(all(file.exists(files)))

txi_file <- file.path(DATA, "txi_gene.rds")

if (file.exists(txi_file)) {
  cat("\nReusing the cached tximport result. Delete", basename(txi_file), "to redo it.\n")
  txi <- readRDS(txi_file)
} else {
  cat("\nImporting", length(files), "kallisto runs with tximport...\n")
  txi <- tximport(files, type = "kallisto", tx2gene = tx2gene,
                  ignoreTxVersion = TRUE, ignoreAfterBar = TRUE)
  saveRDS(txi, txi_file)
}

stopifnot(identical(colnames(txi$counts), samples$library_id))
cat("Gene-level matrix:", nrow(txi$counts), "genes x", ncol(txi$counts), "samples\n")

# ---------------------------------------------------------------- sex
# The metadata has no sex for 52 libraries and no other file on disk carries it,
# but the expression data settles it directly: Xist is female-specific and the Y
# genes Ddx3y/Uty/Eif2s3y/Kdm5d are male-specific. Validated at 100% against the
# 208 libraries whose sex IS recorded - see scripts/13_infer_sex.R.
gene_anno_sex <- fread(anno_file)
xist  <- gene_anno_sex[gene_name == "Xist", gene_id][1]
ygene <- gene_anno_sex[gene_name %in% c("Ddx3y", "Uty", "Eif2s3y", "Kdm5d"), gene_id]
cpm   <- t(t(txi$counts) / colSums(txi$counts)) * 1e6

samples[, sex_score := log2(colSums(cpm[ygene, , drop = FALSE]) + 1) - log2(cpm[xist, ] + 1)]
samples[, sex_final := ifelse(sex_score > 0, "M", "F")]
samples[, sex_source := fifelse(is.na(sex), "expression (absent from metadata)",
                         fifelse(sex == sex_final, "metadata (confirmed by expression)",
                                 "expression (disagrees with metadata)"))]
cat("\nSex after recovery:\n"); print(table(samples$sex_final, samples$age))

# A saccule/utricle pair must come from one animal. The submission workbook lists the two organs of an animal as
# consecutive library IDs (same recorded strain and age, one saccule and one utricle), so that is the animal ID.
# Strain + age + sex alone is not enough: it merged two different BXD128 w22 mice (FTL_1031, FTL_1036; their
# partners are FTL_1030 and FTL_1037) into one pair (fixed October 2026, the same rule as rebuild_291/).
# Libraries without a consecutive partner in the workbook fall back to strain + age + sex.
animal_ids <- function(sheet_file) {
  sh <- rbindlist(lapply(readxl::excel_sheets(sheet_file), function(s) {
    d <- as.data.table(readxl::read_excel(sheet_file, sheet = s))
    setnames(d, trimws(names(d)))
    d[, .(library_id = regmatches(`Library ID`, regexpr("FTL_[0-9]+", `Library ID`)),
          sheet_strain = `mouse strain`, organ = tools::toTitleCase(trimws(tolower(tissue))),
          age = paste0("w", regmatches(Age, regexpr("[0-9]+", Age))))]
  }))
  # tissue labels corrected by expression earlier in the project (PCA)
  sh[library_id %in% c("FTL_900", "FTL_901", "FTL_601"), organ := "Utricle"]
  sh[library_id == "FTL_600", organ := "Saccule"]
  sh[, n := as.integer(sub("FTL_", "", library_id))]
  setorder(sh, n)
  aid <- rep(NA_character_, nrow(sh)); a <- 0L; i <- 1L
  while (i < nrow(sh)) {
    r <- sh[i]; s2 <- sh[i + 1]
    if (s2$n == r$n + 1 && identical(r$sheet_strain, s2$sheet_strain) && r$age == s2$age &&
        setequal(c(r$organ, s2$organ), c("Saccule", "Utricle"))) {
      a <- a + 1L; aid[c(i, i + 1)] <- sprintf("animal_%03d", a); i <- i + 2L
    } else i <- i + 1L
  }
  data.table(library_id = sh$library_id, animal_sheet = aid)
}
samples <- merge(samples, animal_ids(SUBMISSION_SHEET), by = "library_id", all.x = TRUE, sort = FALSE)
samples[, pair_id := fifelse(!is.na(animal_sheet), animal_sheet, paste(strain, age, sex_final, sep = "_"))]
samples[, pair_source := fifelse(!is.na(animal_sheet), "submission workbook", "strain + age + sex")]
samples[, has_organ_pair := .N == 2 && uniqueN(organ) == 2, by = pair_id]
cat("\nOrgan pairs (one animal, both organs):\n")
for (a in c("w8", "w22"))
  cat("  ", a, ":", uniqueN(samples[age == a & has_organ_pair, pair_id]), "\n")

fwrite(samples, file.path(DATA, "sample_table.tsv"), sep = "\t", na = "NA")

# A plain rounded count matrix is handy for quick checks outside DESeq2.
fwrite(data.table(gene_id = rownames(txi$counts), round(txi$counts)),
       file.path(DATA, "gene_counts_raw.tsv.gz"), sep = "\t")

cat("\nDone. Written to", DATA, "\n")
