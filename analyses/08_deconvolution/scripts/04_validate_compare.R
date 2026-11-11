#!/usr/bin/env Rscript
# 04_validate_compare.R  (run in deconv_env; needs only ggplot2 + base)
# ---------------------------------------------------------------------------
# Validate and compare the three deconvolution methods (Bisque, MuSiC, ssGSEA)
# for both tissues. Produces comparison plots + a correlation summary table.
#
# Checks:
#   1. mean composition per method vs OBSERVED single-cell composition
#   2. Bisque vs MuSiC agreement (per-sample, per cell type; Spearman r)
#   3. Bisque(fraction) vs ssGSEA(enrichment) agreement
#   4. per-sample fraction distributions (boxplots)
#   5. Hair_cell vs Supporting_cell anti-correlation (biological sanity)
#
# Inputs : results/tables/*, results/inputs/sc_composition.csv
# Outputs: results/plots/*.png, results/tables/method_correlation_summary.csv
# ---------------------------------------------------------------------------
suppressMessages({ library(ggplot2) })
source("../07_single_cell_rna/scripts/celltype_palette.R")   # shared celltype_pal (kept in 07)

TAB <- "results/tables"
PLT <- "results/plots"
IN  <- "results/inputs"
dir.create(PLT, showWarnings = FALSE, recursive = TRUE)
TISSUES <- c("saccule", "utricle")
KEEP <- c("Hair_cell", "Supporting_cell", "Transitional_epith",
          "Glia_Schwann", "Fibrocyte_mesench", "Roof_dark_cell")

# read a method table -> long df (sample, cell_type, value); NULL if missing
load_long <- function(path, method, tissue) {
    if (!file.exists(path)) { cat("  missing:", path, "\n"); return(NULL) }
    m <- read.csv(path, row.names = 1, check.names = FALSE)
    df <- data.frame(sample = rep(rownames(m), ncol(m)),
                     cell_type = rep(colnames(m), each = nrow(m)),
                     value = as.vector(as.matrix(m)),
                     method = method, tissue = tissue, stringsAsFactors = FALSE)
    df
}

long <- do.call(rbind, unlist(lapply(TISSUES, function(t) list(
    load_long(file.path(TAB, paste0("fractions_bisque_", t, ".csv")),  "Bisque", t),
    load_long(file.path(TAB, paste0("fractions_music_",  t, ".csv")),  "MuSiC",  t),
    load_long(file.path(TAB, paste0("enrichment_ssgsea_",t, ".csv")),  "ssGSEA", t))),
    recursive = FALSE))
long$cell_type <- factor(long$cell_type, levels = KEEP)

# ---- 1. mean composition per method vs observed single-cell -----------------
sc <- read.csv(file.path(IN, "sc_composition.csv"))                 # SubjectName,cell_type,prop
sc <- sc[sc$cell_type %in% KEEP, ]
sc$tissue <- ifelse(sc$SubjectName %in% c("FTL_1146", "FTL_1148"), "saccule", "utricle")
sc_mean <- aggregate(prop ~ cell_type + tissue, sc, mean)
sc_mean$method <- "scRNA (observed)"; names(sc_mean)[names(sc_mean) == "prop"] <- "mean"

frac_mean <- aggregate(value ~ cell_type + tissue + method,
                       long[long$method %in% c("Bisque", "MuSiC"), ], mean)
names(frac_mean)[names(frac_mean) == "value"] <- "mean"
comp <- rbind(frac_mean, sc_mean[, names(frac_mean)])
comp$cell_type <- factor(comp$cell_type, levels = KEEP)
comp$method <- factor(comp$method, levels = c("scRNA (observed)", "Bisque", "MuSiC"))

p1 <- ggplot(comp, aes(method, mean, fill = cell_type)) +
    geom_col() + facet_wrap(~ tissue) +
    scale_fill_manual(values = celltype_pal, drop = FALSE) +
    labs(title = "Mean cell-type composition: deconvolution vs observed scRNA",
         y = "mean proportion", x = NULL) +
    theme_classic() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(PLT, "01_mean_composition_by_method.png"), p1, width = 9, height = 5, dpi = 200, bg = "white")

# ---- helper: wide matrix (sample x cell_type) for one method/tissue ---------
wide <- function(method, tissue) {
    d <- long[long$method == method & long$tissue == tissue, ]
    if (nrow(d) == 0) return(NULL)
    reshape(d[, c("sample", "cell_type", "value")], idvar = "sample",
            timevar = "cell_type", direction = "wide")
}

# ---- 2/3. pairwise method agreement (per-sample) + correlation summary ------
pairs <- list(c("Bisque", "MuSiC"), c("Bisque", "ssGSEA"), c("MuSiC", "ssGSEA"))
scatter_df <- NULL; cor_rows <- list()
for (t in TISSUES) for (pr in pairs) {
    a <- long[long$method == pr[1] & long$tissue == t, ]
    b <- long[long$method == pr[2] & long$tissue == t, ]
    if (nrow(a) == 0 || nrow(b) == 0) next
    m <- merge(a, b, by = c("sample", "cell_type", "tissue"))    # .x = pr1, .y = pr2
    m$pair <- paste(pr[1], "vs", pr[2])
    scatter_df <- rbind(scatter_df, m)
    for (ct in KEEP) {
        sub <- m[m$cell_type == ct, ]
        if (nrow(sub) > 3 && sd(sub$value.x) > 0 && sd(sub$value.y) > 0)
            cor_rows[[length(cor_rows) + 1]] <- data.frame(
                tissue = t, pair = m$pair[1], cell_type = ct,
                spearman = cor(sub$value.x, sub$value.y, method = "spearman"),
                pearson  = cor(sub$value.x, sub$value.y, method = "pearson"))
    }
}
cor_summary <- do.call(rbind, cor_rows)
write.csv(cor_summary, file.path(TAB, "method_correlation_summary.csv"), row.names = FALSE)
cat("[corr] median Spearman by pair:\n")
print(aggregate(spearman ~ pair, cor_summary, median))

# Bisque vs MuSiC scatter (both are fractions -> comparable on same axis)
bm <- scatter_df[scatter_df$pair == "Bisque vs MuSiC", ]
p2 <- ggplot(bm, aes(value.x, value.y, color = tissue)) +
    geom_abline(lty = 2, color = "grey60") + geom_point(alpha = 0.5, size = 1) +
    facet_wrap(~ cell_type, scales = "free") +
    labs(title = "Bisque vs MuSiC estimated fractions (per bulk sample)",
         x = "Bisque fraction", y = "MuSiC fraction") +
    theme_bw()
ggsave(file.path(PLT, "02_bisque_vs_music_scatter.png"), p2, width = 10, height = 6.5, dpi = 200, bg = "white")

# Bisque (fraction) vs ssGSEA (enrichment) -> z-score each within cell type
bs <- scatter_df[scatter_df$pair == "Bisque vs ssGSEA", ]
if (!is.null(bs) && nrow(bs) > 0) {
    bs <- do.call(rbind, by(bs, list(bs$tissue, bs$cell_type), function(d) {
        d$zx <- scale(d$value.x); d$zy <- scale(d$value.y); d }))
    p3 <- ggplot(bs, aes(zx, zy, color = tissue)) +
        geom_point(alpha = 0.5, size = 1) +
        facet_wrap(~ cell_type, scales = "free") +
        labs(title = "Bisque fraction vs ssGSEA enrichment (z-scored per cell type)",
             x = "Bisque fraction (z)", y = "ssGSEA enrichment (z)") +
        theme_bw()
    ggsave(file.path(PLT, "03_bisque_vs_ssgsea_scatter.png"), p3, width = 10, height = 6.5, dpi = 200, bg = "white")
}

# ---- 4. per-sample fraction distributions (Bisque + MuSiC) ------------------
fr <- long[long$method %in% c("Bisque", "MuSiC"), ]
p4 <- ggplot(fr, aes(cell_type, value, fill = cell_type)) +
    geom_boxplot(outlier.size = 0.4) +
    facet_grid(method ~ tissue) +
    scale_fill_manual(values = celltype_pal, drop = FALSE) +
    labs(title = "Estimated cell-type fractions across bulk samples",
         y = "estimated fraction", x = NULL) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "none")
ggsave(file.path(PLT, "04_fraction_boxplots.png"), p4, width = 9, height = 6.5, dpi = 200, bg = "white")

# ---- 5. Hair vs Supporting anti-correlation (Bisque) ------------------------
hc <- do.call(rbind, lapply(TISSUES, function(t) {
    w <- wide("Bisque", t); if (is.null(w)) return(NULL)
    data.frame(sample = w$sample, tissue = t,
               Hair = w$`value.Hair_cell`, Support = w$`value.Supporting_cell`)
}))
if (!is.null(hc)) {
    rr <- sapply(TISSUES, function(t) { s <- hc[hc$tissue == t, ]
        if (nrow(s) > 3) cor(s$Hair, s$Support, method = "spearman") else NA })
    lab <- paste0(TISSUES, ": rho=", round(rr, 2), collapse = "   ")
    p5 <- ggplot(hc, aes(Hair, Support, color = tissue)) +
        geom_point(alpha = 0.6) + geom_smooth(method = "lm", se = FALSE) +
        labs(title = "Hair-cell vs supporting-cell fraction (Bisque)",
             subtitle = lab, x = "Hair_cell fraction", y = "Supporting_cell fraction") +
        theme_bw()
    ggsave(file.path(PLT, "05_hair_vs_supporting.png"), p5, width = 7, height = 5, dpi = 200, bg = "white")
}

# ---- 6. ssGSEA enrichment scores (the robust, GTEx-style output) -----------
en <- long[long$method == "ssGSEA", ]
if (nrow(en) > 0) {
    en$cell_type <- factor(en$cell_type, levels = KEEP)
    p6 <- ggplot(en, aes(cell_type, value, fill = cell_type)) +
        geom_boxplot(outlier.size = 0.4) + facet_wrap(~ tissue) +
        scale_fill_manual(values = celltype_pal, drop = FALSE) +
        labs(title = "ssGSEA marker-enrichment scores across bulk samples",
             subtitle = "Per-sample cell-type abundance proxy (robust to the nuclei-vs-bulk shift; feeds interaction-eQTL)",
             y = "ssGSEA enrichment score", x = NULL) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "none")
    ggsave(file.path(PLT, "06_ssgsea_enrichment.png"), p6, width = 9, height = 5, dpi = 200, bg = "white")
}

cat("DONE — plots in", PLT, "\n")
