# 00 - Shared settings for the 291-library rebuild.
#
# This is a SEPARATE analysis from the 260-library one in ../. Nothing here
# writes outside rebuild_291/. The v1 results are left untouched for comparison.
#
# What changes relative to v1:
#   * 291 libraries instead of 260. Everything with RNA-to-genotype concordance
#     >= 0.90 is kept. v1 additionally kept only one library per strain per
#     organ-age group, a rule that exists because QTL mapping needs one genotype
#     per sample; differential expression has no such requirement and a second
#     animal of the same strain is a valid biological replicate.
#   * Animal identity comes from the submission record (consecutive library IDs,
#     same strain, the two organs) rather than being inferred from strain, age
#     and sex. That fixes one w22 pair in v1 that merged two different mice.

V1     <- "../../analyses/05_dge"
BASE   <- file.path(V1, "rebuild_291")
PANTRY <- file.path(dirname(V1), "03_pantry")
PHENO  <- file.path(PANTRY, "phenotyping")   # first run read the plain Pantry run here (see README)

DATA   <- file.path(BASE, "data")
TABLES <- file.path(BASE, "results", "tables")
FIGS   <- file.path(BASE, "results", "figures")

DATASETS <- c("saccule_w8", "saccule_w22", "utricle_w8", "utricle_w22")
KEEP_BIOTYPES <- c("protein_coding", "lincRNA")
FDR_CUT <- 0.05
LFC_CUT <- 1

ORGAN_COLS <- c(Saccule = "#1B7837", Utricle = "#762A83")
GROUP_COLS <- c("Saccule w8"  = "#7FBF7B", "Saccule w22" = "#1B7837",
                "Utricle w8"  = "#C2A5CF", "Utricle w22" = "#762A83")
SEX_COLS   <- c(F = "#D6604D", M = "#4393C3")

for (d in c(DATA, TABLES, FIGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
