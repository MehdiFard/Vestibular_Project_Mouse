# 01 - Over-representation analysis (ORA) for every contrast.
#
# The classical approach, and the counterpart to ../gsea:
#   GSEA  ranks ALL genes and asks whether a set clusters at one end
#   ORA   takes only the DE genes and asks whether a set is over-represented
#         in that list, by a hypergeometric (Fisher) test
#
# Two choices matter and are easy to get wrong:
#
#   1. The universe. The background is the genes TESTED in that contrast, not
#      every gene in the genome and not every gene in MSigDB. Using a bigger
#      universe than was actually tested inflates every p-value's significance.
#   2. Direction. Up and down genes are tested separately as well as together.
#      A set that is half up and half down looks enriched in the combined list
#      while being biologically incoherent.
#
# Same gene sets and same significance rule as the GSEA run, so the two are
# directly comparable.

BASE <- "../../../analyses/05_dge"
ROOT <- file.path(dirname(BASE), "06_enrichment")
PA   <- file.path(ROOT, "over_rep")
Sys.setenv(TMPDIR = file.path(ROOT, "logs"))

library(data.table)
library(clusterProfiler)
library(msigdbr)

set.seed(42)
FDR_CUT <- 0.05      # for calling a gene DE, and for calling a set enriched
LFC_CUT <- 1

de <- readRDS(file.path(BASE, "data", "de_results.rds"))

CONTRASTS <- list(
  A_organ_w8    = list(key = "A_utricle_vs_saccule_w8",   label = "Utricle vs Saccule, w8",
                       up = "higher in Utricle", down = "higher in Saccule"),
  B_organ_w22   = list(key = "B_utricle_vs_saccule_w22",  label = "Utricle vs Saccule, w22",
                       up = "higher in Utricle", down = "higher in Saccule"),
  C_age_sac     = list(key = "C_w22_vs_w8_saccule",       label = "w22 vs w8, Saccule",
                       up = "higher at w22", down = "higher at w8"),
  D_age_utr     = list(key = "D_w22_vs_w8_utricle",       label = "w22 vs w8, Utricle",
                       up = "higher at w22", down = "higher at w8"),
  E_interaction = list(key = "E_organ_by_age_interaction", label = "Organ x age interaction",
                       up = "gap widens toward Utricle", down = "gap widens toward Saccule"))

COLLECTIONS <- list(
  hallmark = list(coll = "MH", sub = NULL,          name = "Hallmark"),
  go_bp    = list(coll = "M5", sub = "GO:BP",       name = "GO Biological Process"),
  reactome = list(coll = "M2", sub = "CP:REACTOME", name = "Reactome"))

# TERM2GENE: two columns, gene set name and Ensembl gene id.
get_t2g <- function(cn) {
  s <- COLLECTIONS[[cn]]
  df <- if (is.null(s$sub)) msigdbr(db_species = "MM", species = "Mus musculus",
                                    collection = s$coll)
        else msigdbr(db_species = "MM", species = "Mus musculus",
                     collection = s$coll, subcollection = s$sub)
  df <- as.data.table(df)[!is.na(ensembl_gene)]
  unique(df[, .(gs_name, ensembl_gene)])
}

cat("Loading gene sets...\n")
T2G <- lapply(names(COLLECTIONS), get_t2g)
names(T2G) <- names(COLLECTIONS)
for (cn in names(T2G)) cat(sprintf("  %-9s %5d sets\n", cn, uniqueN(T2G[[cn]]$gs_name)))

pretty_name <- function(x) {
  x <- sub("^(HALLMARK|GOBP|REACTOME)_", "", x)
  x <- gsub("_", " ", x)
  paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2)))
}

all_rows <- list()
gene_counts <- list()

for (cid in names(CONTRASTS)) {
  ct  <- CONTRASTS[[cid]]
  res <- as.data.table(de[[ct$key]])

  # Universe = everything that was actually tested in this contrast.
  universe <- res[!is.na(padj), gene_id]
  sig      <- res[!is.na(padj) & padj < FDR_CUT & abs(log2FoldChange) > LFC_CUT]
  lists <- list(all  = sig$gene_id,
                up   = sig[log2FoldChange > 0, gene_id],
                down = sig[log2FoldChange < 0, gene_id])

  gene_counts[[cid]] <- data.table(contrast = cid, contrast_label = ct$label,
                                   universe = length(universe),
                                   n_all = length(lists$all),
                                   n_up = length(lists$up), n_down = length(lists$down))

  cat(sprintf("\n===== %s =====\n  universe %d | DE %d (up %d / down %d)\n",
              ct$label, length(universe), length(lists$all),
              length(lists$up), length(lists$down)))

  for (cn in names(T2G)) {
    for (dir in names(lists)) {
      genes <- lists[[dir]]
      if (length(genes) < 5) {
        cat(sprintf("  %-9s %-4s  skipped, only %d genes\n", cn, dir, length(genes)))
        next
      }
      e <- enricher(gene = genes, universe = universe, TERM2GENE = T2G[[cn]],
                    pvalueCutoff = 1, qvalueCutoff = 1,
                    minGSSize = 15, maxGSSize = 500)
      if (is.null(e) || !nrow(as.data.frame(e))) next
      d <- as.data.table(as.data.frame(e))
      setnames(d, "ID", "pathway")
      d[, `:=`(contrast = cid, contrast_label = ct$label, collection = cn,
               gene_list = dir, pathway_label = pretty_name(pathway),
               direction = switch(dir, up = ct$up, down = ct$down, all = "up or down"))]
      # GeneRatio / BgRatio are "k/n" strings; turn them into a fold enrichment.
      d[, k := as.numeric(sub("/.*", "", GeneRatio))]
      d[, n := as.numeric(sub(".*/", "", GeneRatio))]
      d[, K := as.numeric(sub("/.*", "", BgRatio))]
      d[, N := as.numeric(sub(".*/", "", BgRatio))]
      d[, fold_enrichment := (k / n) / (K / N)]
      setorder(d, p.adjust, -fold_enrichment, na.last = TRUE)
      setcolorder(d, c("contrast", "contrast_label", "collection", "gene_list",
                       "pathway", "pathway_label", "fold_enrichment",
                       "pvalue", "p.adjust", "k", "n", "K", "N", "direction"))
      fwrite(d, file.path(PA, "tables", cn, paste0("ora_", cid, "_", dir, ".tsv")), sep = "\t")
      all_rows[[paste(cid, cn, dir)]] <- d
      cat(sprintf("  %-9s %-4s  tested %4d | FDR<%.2f: %4d\n",
                  cn, dir, nrow(d), FDR_CUT, d[p.adjust < FDR_CUT, .N]))
    }
  }
}

comb <- rbindlist(all_rows, fill = TRUE)
fwrite(comb, file.path(PA, "tables", "ora_all_contrasts.tsv.gz"), sep = "\t")
fwrite(rbindlist(gene_counts), file.path(PA, "tables", "ora_input_gene_counts.tsv"), sep = "\t")

summ <- comb[, .(sets_tested = .N,
                 sig = sum(p.adjust < FDR_CUT),
                 top = pathway_label[which.min(p.adjust)],
                 top_fold = round(fold_enrichment[which.min(p.adjust)], 2)),
             by = .(contrast, contrast_label, collection, gene_list)]
fwrite(summ, file.path(PA, "tables", "ora_summary.tsv"), sep = "\t")

saveRDS(list(combined = comb, counts = rbindlist(gene_counts)),
        file.path(PA, "tables", "ora_objects.rds"))

cat("\n\n================ SUMMARY (significant sets) ================\n")
print(dcast(summ, contrast + collection ~ gene_list, value.var = "sig"))
cat("\nDone.\n")
