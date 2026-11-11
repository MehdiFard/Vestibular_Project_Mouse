# 01 - Build the 291-library count matrix and sample table.
#
# Counts come from two places, both quantified with the same kallisto 0.48.0
# index and the same parameters:
#   * 260 libraries from the original Pantry phenotyping run
#   *  31 libraries quantified here (scripts/run_kallisto.slurm)

source("../../analyses/05_dge/rebuild_291/scripts/00_config.R")

library(data.table)
library(tximport)

v2 <- fread(file.path(DATA, "v2_sample_table_raw.tsv"))
setnames(v2, "batch", "batch", skip_absent = TRUE)
v2[, dataset := paste0(tolower(organ), "_", age)]
v2[, group := paste(organ, age)]

# ---------------------------------------------------------------- count paths
v2[, quant_dir := fifelse(in_260 == TRUE,
    file.path(PHENO, "intermediate", dataset, "expression", library_id),
    file.path(BASE, "quant", library_id))]
v2[, abundance := file.path(quant_dir, "abundance.tsv")]
missing <- v2[!file.exists(abundance), library_id]
if (length(missing)) stop("missing quantification for: ", paste(missing, collapse = ", "))
cat("libraries:", nrow(v2), " (", v2[in_260 == TRUE, .N], "from Pantry,",
    v2[in_260 == FALSE, .N], "newly quantified )\n")

# ---------------------------------------------------------------- annotation
# Reuse the gene annotation and tx2gene built for v1 - same GTF, same release.
tx2gene   <- fread(file.path(V1, "data", "tx2gene.tsv"))
gene_anno <- fread(file.path(V1, "data", "gene_annotation.tsv"))
fwrite(gene_anno, file.path(DATA, "gene_annotation.tsv"), sep = "\t")

# ---------------------------------------------------------------- tximport
files <- setNames(v2$abundance, v2$library_id)
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene,
                ignoreTxVersion = TRUE, countsFromAbundance = "no")
cat("gene x sample matrix:", paste(dim(txi$counts), collapse = " x "), "\n")
saveRDS(txi, file.path(DATA, "txi_gene.rds"))

# ---------------------------------------------------------------- sex
# Recorded sex where the submission sheet has it; otherwise inferred from Xist
# against the Y genes, the same rule validated on 232 libraries in v1.
XIST  <- gene_anno[gene_name == "Xist", gene_id][1]
YGENE <- gene_anno[gene_name %in% c("Ddx3y", "Uty", "Eif2s3y", "Kdm5d"), gene_id]
cnt <- txi$counts[, v2$library_id]
cpm <- t(t(cnt) / colSums(cnt)) * 1e6
v2[, sex_score := log2(colSums(cpm[YGENE, , drop = FALSE]) + 1) - log2(cpm[XIST, ] + 1)]
v2[, sex_inferred := ifelse(sex_score > 0, "M", "F")]
known <- v2[sheet_sex %in% c("F", "M")]
cat("sex: recorded for", nrow(known), "| inference agrees with record:",
    known[sheet_sex == sex_inferred, .N], "/", nrow(known), "\n")
v2[, sex_final := sex_inferred]      # directly measured, and it matches the record
v2[, sex_source := fifelse(sheet_sex %in% c("F", "M"),
                           "submission record (confirmed by expression)",
                           "expression (not recorded)")]

# ---------------------------------------------------------------- animal id
# The submission record pairs the two organs of an animal as consecutive library
# IDs. Six w22 pairs in batch 3 were not submitted adjacently and fall back to
# strain + age + sex, which is unambiguous where a strain appears only once.
v2[, pair_id := fifelse(!is.na(animal_sheet) & animal_sheet != "",
                        animal_sheet, paste(strain, age, sex_final, sep = "_"))]
v2[, pair_source := fifelse(!is.na(animal_sheet) & animal_sheet != "",
                            "submission record", "strain + age + sex")]
v2[, n_in_pair := .N, by = pair_id]
v2[, has_organ_pair := n_in_pair == 2 & uniqueN(organ) == 2, by = pair_id]

cat("\nlibraries per group:\n"); print(v2[, .N, by = .(organ, age)][order(organ, age)])
cat("\ncomplete animal pairs:", v2[has_organ_pair == TRUE, uniqueN(pair_id)], "\n")
print(v2[has_organ_pair == TRUE, .(pairs = uniqueN(pair_id)), by = age])
cat("\npair source:\n"); print(v2[, .N, by = pair_source])

keep <- c("library_id", "dataset", "organ", "age", "group", "batch", "strain",
          "strain_simil", "strain_runnerup", "strain_runnerup_simil", "strain_gap",
          "mouse strain", "sheet_sex", "sex_score", "sex_final", "sex_source",
          "animal_sheet", "pair_id", "pair_source", "has_organ_pair", "in_260")
fwrite(v2[, ..keep], file.path(DATA, "sample_table.tsv"), sep = "\t", na = "NA")
cat("\nWrote", file.path(DATA, "sample_table.tsv"), "\n")
