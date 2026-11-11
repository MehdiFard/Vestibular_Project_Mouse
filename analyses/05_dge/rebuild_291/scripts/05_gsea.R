# 05 - GSEA on the 291-library results, and the comparison against the 260-library GSEA.
# Runs in the enrich_r environment with the shared project-local library.
V1   <- "../../analyses/05_dge"
BASE <- file.path(V1, "rebuild_291")
PA   <- file.path(dirname(V1), "06_enrichment")
.libPaths(c(file.path(PA, "rlib"), .libPaths()))
Sys.setenv(TMPDIR = file.path(BASE, "logs"))
DATA <- file.path(BASE, "data"); TABLES <- file.path(BASE, "results", "tables"); FIGS <- file.path(BASE, "results", "figures")

library(data.table); library(fgsea); library(msigdbr); library(ggplot2)
set.seed(42); FDR_CUT <- 0.05

de <- readRDS(file.path(DATA, "de_results.rds"))
CONTRASTS <- list(
  A_organ_w8 = list(key = "A_utricle_vs_saccule_w8", label = "Utricle vs Saccule, w8",
                    up = "higher in Utricle", down = "higher in Saccule"),
  B_organ_w22 = list(key = "B_utricle_vs_saccule_w22", label = "Utricle vs Saccule, w22",
                     up = "higher in Utricle", down = "higher in Saccule"),
  C_age_sac = list(key = "C_w22_vs_w8_saccule", label = "w22 vs w8, Saccule",
                   up = "higher at w22", down = "higher at w8"),
  D_age_utr = list(key = "D_w22_vs_w8_utricle", label = "w22 vs w8, Utricle",
                   up = "higher at w22", down = "higher at w8"),
  E_interaction = list(key = "E_organ_by_age_interaction", label = "Organ x age interaction",
                       up = "gap widens toward Utricle", down = "gap widens toward Saccule"))
COLLECTIONS <- list(hallmark = list(coll = "MH", sub = NULL),
                    go_bp = list(coll = "M5", sub = "GO:BP"),
                    reactome = list(coll = "M2", sub = "CP:REACTOME"))
get_sets <- function(cn) {
  s <- COLLECTIONS[[cn]]
  df <- if (is.null(s$sub)) msigdbr(db_species = "MM", species = "Mus musculus", collection = s$coll)
        else msigdbr(db_species = "MM", species = "Mus musculus", collection = s$coll, subcollection = s$sub)
  df <- as.data.table(df)[!is.na(ensembl_gene)]
  split(df$ensembl_gene, df$gs_name)
}
SETS <- lapply(names(COLLECTIONS), get_sets); names(SETS) <- names(COLLECTIONS)
pretty <- function(x) { x <- sub("^(HALLMARK|GOBP|REACTOME)_", "", x); x <- gsub("_", " ", x)
  paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2))) }

all_rows <- list()
for (cid in names(CONTRASTS)) {
  ct <- CONTRASTS[[cid]]; res <- as.data.table(de[[ct$key]])
  r <- res[!is.na(stat), .(gene_id, stat)]
  ranks <- sort(setNames(r$stat, r$gene_id), decreasing = TRUE)
  cat("\n=====", ct$label, "- ranked genes:", length(ranks), "=====\n")
  for (cn in names(SETS)) {
    fg <- as.data.table(fgsea(SETS[[cn]], ranks, minSize = 15, maxSize = 500, eps = 0))
    fg[, leadingEdge := vapply(leadingEdge, function(g) paste(head(g, 50), collapse = ","), "")]
    fg[, `:=`(contrast = cid, contrast_label = ct$label, collection = cn,
              pathway_label = pretty(pathway), direction = fifelse(NES > 0, ct$up, ct$down))]
    fg[, abs_nes := abs(NES)]; setorder(fg, padj, -abs_nes, na.last = TRUE); fg[, abs_nes := NULL]
    all_rows[[paste(cid, cn)]] <- fg
    cat(sprintf("  %-9s tested %5d | FDR<0.05: %4d\n", cn, nrow(fg), fg[padj < FDR_CUT, .N]))
  }
}
comb <- rbindlist(all_rows)
fwrite(comb, file.path(TABLES, "gsea_all_contrasts.tsv.gz"), sep = "\t")
summ <- comb[, .(pathways_tested = .N, sig = sum(padj < FDR_CUT),
                 sig_up = sum(padj < FDR_CUT & NES > 0), sig_down = sum(padj < FDR_CUT & NES < 0),
                 top = pathway_label[which.min(padj)], top_NES = round(NES[which.min(padj)], 2)),
             by = .(contrast, contrast_label, collection)]
fwrite(summ, file.path(TABLES, "gsea_summary.tsv"), sep = "\t")
print(summ[, .(contrast, collection, sig, top)])

# ---- comparison with the 260-library GSEA
g1 <- fread(file.path(PA, "gsea", "tables", "gsea_all_contrasts.tsv.gz"))
h1 <- g1[collection == "hallmark", .(contrast, pathway, nes1 = NES, p1 = padj)]
h2 <- comb[collection == "hallmark", .(contrast, pathway, pathway_label, nes2 = NES, p2 = padj)]
hm <- merge(h1, h2, by = c("contrast", "pathway"))
hm[, call := fifelse(p1 < FDR_CUT & p2 < FDR_CUT, "both",
              fifelse(p1 < FDR_CUT, "260 only", fifelse(p2 < FDR_CUT, "291 only", "neither")))]
fwrite(hm, file.path(TABLES, "gsea_v1_vs_v2_hallmark.tsv"), sep = "\t")
cat("\nHallmark sets by call:\n"); print(hm[, .N, by = call])
cat("NES correlation between runs:", round(cor(hm$nes1, hm$nes2), 3), "\n")

cnt <- rbind(g1[padj < FDR_CUT, .(n = .N), by = .(contrast, collection)][, run := "260 libraries"],
             comb[padj < FDR_CUT, .(n = .N), by = .(contrast, collection)][, run := "291 libraries"])
CN <- c(hallmark = "Hallmark", go_bp = "GO Biological Process", reactome = "Reactome")
cnt[, collection := factor(CN[collection], levels = CN)]
p <- ggplot(cnt, aes(contrast, n, fill = run)) +
  geom_col(position = position_dodge(0.8), width = 0.72) +
  facet_wrap(~ collection, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(`260 libraries` = "#9E9E9E", `291 libraries` = "#1B7837"), name = NULL) +
  labs(title = "Enriched gene sets: 260 vs 291 libraries", subtitle = "FDR < 0.05",
       x = NULL, y = "Gene sets") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 18, hjust = 1), legend.position = "top",
        panel.grid.minor = element_blank())
ggsave(file.path(FIGS, "30_gsea_v1_vs_v2.png"), p, width = 9, height = 9, dpi = 300, bg = "white")

p <- ggplot(hm, aes(nes1, nes2, colour = call)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(size = 1.9, alpha = 0.85) +
  scale_colour_manual(values = c(neither = "grey80", `260 only` = "#9E9E9E",
                                 `291 only` = "#1B7837", both = "#D95F02"), name = NULL) +
  labs(title = "Hallmark enrichment scores agree between runs",
       subtitle = paste0("r = ", round(cor(hm$nes1, hm$nes2), 3), " over ", nrow(hm),
                         " set-by-contrast comparisons"),
       x = "NES, 260 libraries", y = "NES, 291 libraries") +
  theme_bw(base_size = 12) + theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave(file.path(FIGS, "31_gsea_nes_agreement.png"), p, width = 8, height = 8, dpi = 300, bg = "white")
cat("\nGSEA done.\n")
