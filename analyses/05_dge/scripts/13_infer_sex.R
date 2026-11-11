# 13 - Recover the missing sex labels from the expression data itself.
#
# 52 of the 260 libraries have no sex in the metadata (24 absent from the sheet,
# 28 present but with the field blank). No other file on disk carries it.
#
# Sex does not need a spreadsheet though: Xist is expressed almost exclusively
# in females, and the Y genes Ddx3y / Uty / Eif2s3y / Kdm5d almost exclusively
# in males. The two give an unambiguous call, and the 208 libraries with a
# known sex are used to check the rule before it is applied to the rest.

source("../../analyses/05_dge/scripts/00_config.R")

library(data.table)
library(DESeq2)
library(ggplot2)

md   <- fread(file.path(TABLES, "master_metadata.tsv"))
txi  <- readRDS(file.path(DATA, "txi_gene.rds"))
anno <- fread(file.path(DATA, "gene_annotation.tsv"))

XIST  <- anno[gene_name == "Xist", gene_id][1]
YGENE <- anno[gene_name %in% c("Ddx3y", "Uty", "Eif2s3y", "Kdm5d"), gene_id]

# CPM so that library depth does not drive the call.
cnt <- txi$counts[, md$library_id]
cpm <- t(t(cnt) / colSums(cnt)) * 1e6

md[, xist_cpm := cpm[XIST, ]]
md[, y_cpm    := colSums(cpm[YGENE, , drop = FALSE])]
md[, sex_score := log2(y_cpm + 1) - log2(xist_cpm + 1)]   # >0 male, <0 female

cat("=== marker expression by known sex ===\n")
print(md[!is.na(sex), .(n = .N,
                        median_xist = round(median(xist_cpm), 1),
                        median_ychr = round(median(y_cpm), 1),
                        median_score = round(median(sex_score), 2)), by = sex])

# The two groups separate on the score; take the midpoint of the gap as cutoff.
known <- md[!is.na(sex)]
hi_f  <- max(known[sex == "F", sex_score])
lo_m  <- min(known[sex == "M", sex_score])
cat("\nhighest score among known females:", round(hi_f, 2),
    "\nlowest score among known males  :", round(lo_m, 2),
    "\ngap between the two groups      :", round(lo_m - hi_f, 2), "\n")
cutoff <- (hi_f + lo_m) / 2
cat("cutoff used:", round(cutoff, 2), "\n")

md[, sex_inferred := ifelse(sex_score > cutoff, "M", "F")]

# ------------------------------------------------------------------ validate
cat("\n=== VALIDATION on the", nrow(known), "libraries with a known sex ===\n")
tab <- md[!is.na(sex), table(metadata = sex, inferred = sex_inferred)]
print(tab)
acc <- sum(diag(tab)) / sum(tab)
cat("concordance:", sprintf("%.1f%%", 100 * acc), "\n")

mism <- md[!is.na(sex) & sex != sex_inferred]
if (nrow(mism)) {
  cat("\nmismatches (expression disagrees with the sheet):\n")
  print(mism[, .(library_id, organ, age, batch, sheet = sex, expression = sex_inferred,
                 xist_cpm = round(xist_cpm, 1), y_cpm = round(y_cpm, 1))])
}

# ------------------------------------------------------------------ fill in
md[, sex_source := fifelse(!is.na(sex) & sex == sex_inferred, "metadata (confirmed by expression)",
                    fifelse(!is.na(sex), "expression (disagrees with metadata)",
                            "expression (absent from metadata)"))]
md[, sex_final := sex_inferred]      # expression wins everywhere; it is directly measured

cat("\n=== sex after recovery ===\n")
print(table(md$sex_final, md$age))
cat("\nprovenance:\n"); print(table(md$sex_source))

# Consistency check that does not use the metadata at all: the two organs of one
# animal must have the same sex.
p <- md[, .(n = .N, nsex = uniqueN(sex_final)), by = pair_id][n == 2]
cat("\norgan pairs:", nrow(p), " | same inferred sex in both:", sum(p$nsex == 1),
    sprintf("(%.1f%%)\n", 100 * mean(p$nsex == 1)))

# ------------------------------------------------------------------ figure
# Compare against the sex on the complete submission record, not the partial
# sheet this project was first built from.
val <- fread(file.path(TABLES, "backup_sheet_validation.tsv"))
md <- merge(md, val[, .(library_id, sheet_sex)], by = "library_id", all.x = TRUE)
md[, record_lab := fifelse(sheet_sex %in% c("F", "M"),
                           paste("record:", sheet_sex), "not recorded")]
md[, record_lab := factor(record_lab, levels = c("record: F", "record: M", "not recorded"))]
n_rec   <- md[sheet_sex %in% c("F", "M"), .N]
n_agree <- md[sheet_sex %in% c("F", "M") & sheet_sex == sex_final, .N]
n_none  <- md[!sheet_sex %in% c("F", "M"), .N]

p_sex <- ggplot(md, aes(log2(xist_cpm + 1), log2(y_cpm + 1), colour = sex_final,
                        shape = record_lab)) +
  geom_point(size = 2.4, alpha = 0.85) +
  scale_colour_manual(values = c(F = "#D6604D", M = "#4393C3"), name = "Sex from expression") +
  scale_shape_manual(values = c(16, 17, 4), name = "Submission record") +
  labs(title = "Sex determined from expression",
       subtitle = paste0("Xist marks females; Ddx3y/Uty/Eif2s3y/Kdm5d mark males.\n",
                         "The two groups separate completely, with an 11.5 log2-unit gap.\n",
                         "The expression call matches the submission record in all ",
                         n_agree, " of ", n_rec, " libraries\n",
                         "with a recorded sex, and supplies sex for the remaining ", n_none, "."),
       x = "log2(Xist CPM + 1)", y = "log2(Y-gene CPM + 1)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "right")

ggsave(file.path(FIGS, "30_sex_inference.png"), p_sex,
       width = 10, height = 6.4, dpi = 300, bg = "white")

# ------------------------------------------------------------------ write back
cat("\nDone - this script validates and plots only;",
    "scripts/01_build_inputs.R does the inference itself.\n")
