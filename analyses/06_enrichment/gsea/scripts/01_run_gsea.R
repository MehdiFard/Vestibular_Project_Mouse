# 01 - GSEA for every differential-expression contrast.
#
# Input : the DESeq2 results already computed in ../data/de_results.rds
# Ranking: the Wald statistic. It carries the direction of the effect and its
#          precision in one number, which is what GSEA wants - ranking on raw
#          log2FC alone lets noisy low-count genes sit at the extremes.
# Method : fgsea (multilevel), run on the whole ranked list. No FDR/fold-change
#          cut is applied first - thresholding before GSEA throws away exactly
#          the graded signal the method is built to detect.
# Sets   : mouse-native MSigDB 2026.1.Mm, keyed by Ensembl gene ID so no
#          ID conversion is needed anywhere.
#
#   MH            Hallmark            50 well-curated summary sets
#   M5 GO:BP      GO Biological Process
#   M2 CP:REACTOME Reactome pathways

BASE <- "../../../analyses/05_dge"
ROOT <- file.path(dirname(BASE), "06_enrichment")   # logs are shared
PA   <- file.path(ROOT, "gsea")                # this analysis's own outputs
Sys.setenv(TMPDIR = file.path(ROOT, "logs"))

library(data.table)
library(fgsea)
library(msigdbr)

set.seed(42)
FDR_CUT <- 0.05

# ------------------------------------------------------------------ contrasts
de <- readRDS(file.path(BASE, "data", "de_results.rds"))

CONTRASTS <- list(
  A_organ_w8  = list(key = "A_utricle_vs_saccule_w8",
                     label = "Utricle vs Saccule, w8",
                     up = "higher in Utricle", down = "higher in Saccule"),
  B_organ_w22 = list(key = "B_utricle_vs_saccule_w22",
                     label = "Utricle vs Saccule, w22",
                     up = "higher in Utricle", down = "higher in Saccule"),
  C_age_sac   = list(key = "C_w22_vs_w8_saccule",
                     label = "w22 vs w8, Saccule",
                     up = "higher at w22", down = "higher at w8"),
  D_age_utr   = list(key = "D_w22_vs_w8_utricle",
                     label = "w22 vs w8, Utricle",
                     up = "higher at w22", down = "higher at w8"),
  E_interaction = list(key = "E_organ_by_age_interaction",
                     label = "Organ x age interaction",
                     up = "organ gap widens toward Utricle",
                     down = "organ gap widens toward Saccule"))

# ------------------------------------------------------------------ gene sets
COLLECTIONS <- list(
  hallmark = list(coll = "MH", sub = NULL,           name = "Hallmark"),
  go_bp    = list(coll = "M5", sub = "GO:BP",        name = "GO Biological Process"),
  reactome = list(coll = "M2", sub = "CP:REACTOME",  name = "Reactome"))

get_sets <- function(cn) {
  s <- COLLECTIONS[[cn]]
  df <- if (is.null(s$sub)) msigdbr(db_species = "MM", species = "Mus musculus",
                                    collection = s$coll)
        else msigdbr(db_species = "MM", species = "Mus musculus",
                     collection = s$coll, subcollection = s$sub)
  df <- as.data.table(df)[!is.na(ensembl_gene)]
  split(df$ensembl_gene, df$gs_name)
}

cat("Loading gene sets...\n")
SETS <- lapply(names(COLLECTIONS), get_sets)
names(SETS) <- names(COLLECTIONS)
for (cn in names(SETS)) cat(sprintf("  %-9s %5d sets\n", cn, length(SETS[[cn]])))

# Tidy the MSigDB names for printing: HALLMARK_TNFA_SIGNALING -> "Tnfa signaling"
pretty_name <- function(x) {
  x <- sub("^(HALLMARK|GOBP|REACTOME)_", "", x)
  x <- gsub("_", " ", x)
  paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2)))
}

# ------------------------------------------------------------------ run
all_rows <- list()

for (cid in names(CONTRASTS)) {
  ct  <- CONTRASTS[[cid]]
  res <- as.data.table(de[[ct$key]])

  # One rank per gene, highest stat first.
  r <- res[!is.na(stat), .(gene_id, stat)]
  r <- r[, .(stat = mean(stat)), by = gene_id]          # no duplicate ids expected
  ranks <- setNames(r$stat, r$gene_id)
  ranks <- sort(ranks, decreasing = TRUE)

  cat("\n=====", ct$label, "- ranked genes:", length(ranks), "=====\n")

  for (cn in names(SETS)) {
    fg <- fgsea(pathways = SETS[[cn]], stats = ranks,
                minSize = 15, maxSize = 500, eps = 0)
    fg <- as.data.table(fg)
    fg[, leadingEdge := vapply(leadingEdge, function(g) paste(head(g, 50), collapse = ","), "")]
    fg[, `:=`(contrast = cid, contrast_label = ct$label,
              collection = cn, pathway_label = pretty_name(pathway),
              direction = fifelse(NES > 0, ct$up, ct$down))]
    # setorder cannot take an expression, so sort on an explicit helper column.
    fg[, abs_nes := abs(NES)]
    setorder(fg, padj, -abs_nes, na.last = TRUE)
    fg[, abs_nes := NULL]
    setcolorder(fg, c("contrast", "contrast_label", "collection",
                      "pathway", "pathway_label", "NES", "pval", "padj",
                      "size", "direction"))
    fwrite(fg, file.path(PA, "tables", cn, paste0("gsea_", cid, ".tsv")), sep = "\t")
    all_rows[[paste(cid, cn)]] <- fg

    cat(sprintf("  %-9s tested %5d | FDR<%.2f: %4d  (up %3d / down %3d)\n",
                cn, nrow(fg), FDR_CUT, fg[padj < FDR_CUT, .N],
                fg[padj < FDR_CUT & NES > 0, .N], fg[padj < FDR_CUT & NES < 0, .N]))
  }
}

# ------------------------------------------------------------------ combined
comb <- rbindlist(all_rows)
fwrite(comb, file.path(PA, "tables", "gsea_all_contrasts.tsv.gz"), sep = "\t")

summ <- comb[, .(pathways_tested = .N,
                 sig = sum(padj < FDR_CUT),
                 sig_up = sum(padj < FDR_CUT & NES > 0),
                 sig_down = sum(padj < FDR_CUT & NES < 0),
                 top = pathway_label[which.min(padj)],
                 top_NES = round(NES[which.min(padj)], 2)),
             by = .(contrast, contrast_label, collection)]
fwrite(summ, file.path(PA, "tables", "gsea_summary.tsv"), sep = "\t")

cat("\n\n================ SUMMARY ================\n")
print(summ[, .(contrast, collection, pathways_tested, sig, sig_up, sig_down)])
saveRDS(list(combined = comb, sets = SETS), file.path(PA, "tables", "gsea_objects.rds"))
cat("\nDone.\n")
