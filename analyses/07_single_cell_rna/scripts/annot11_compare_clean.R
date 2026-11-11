# annot11_compare_clean.R  (seurat_env)
# Compare cleaned vs previous (uncleaned) annotations: did cleaning recover any
# cell type?  Builds a per-cell-type before/after table for each integration.
# setwd("analyses/07_single_cell_rna")   # from the project root

rm(list = ls())
suppressMessages({ library(Seurat); library(tidyverse) })
outdir <- "results/cell_annotation_clean"

compare <- function(old_path, new_path, label) {
    old <- readRDS(old_path); new <- readRDS(new_path)
    to <- table(old$cell_type); tn <- table(new$cell_type)
    types <- sort(union(names(to), names(tn)))
    df <- data.frame(comparison = label, cell_type = types,
                     n_before = as.integer(to[types]), n_after = as.integer(tn[types]))
    df$n_before[is.na(df$n_before)] <- 0; df$n_after[is.na(df$n_after)] <- 0
    df$status <- ifelse(df$n_before == 0 & df$n_after > 0, "recovered (new)",
                 ifelse(df$n_before > 0 & df$n_after == 0, "lost", "present both"))
    attr(df, "totals") <- c(before = ncol(old), after = ncol(new))
    df
}

pairs <- list(
    saccule     = c("results/per_tissue/saccule_integrated_annotated.rds",
                    file.path(outdir, "saccule_integrated_clean.rds")),
    utricle     = c("results/per_tissue/utricle_integrated_annotated.rds",
                    file.path(outdir, "utricle_integrated_clean.rds")),
    all_samples = c("results/cell_annotation/integrated_all_samples_annotated.rds",
                    file.path(outdir, "all_integrated_clean.rds")))

cmp <- list(); totals <- list()
for (nm in names(pairs)) {
    d <- compare(pairs[[nm]][1], pairs[[nm]][2], nm)
    cmp[[nm]] <- d; totals[[nm]] <- attr(d, "totals")
}
cmp_all <- do.call(rbind, cmp)
write.csv(cmp_all, file.path(outdir, "celltype_recovery_comparison.csv"), row.names = FALSE)

cat("\n==== total cells before vs after cleaning ====\n")
for (nm in names(totals)) cat(sprintf("  %-12s before=%d  after=%d  (%.1f%% removed)\n",
    nm, totals[[nm]]["before"], totals[[nm]]["after"],
    100 * (totals[[nm]]["before"] - totals[[nm]]["after"]) / totals[[nm]]["before"]))
cat("\n==== recovered / lost cell types ====\n")
print(cmp_all[cmp_all$status != "present both", ])
cat("\n==== full comparison ====\n"); print(cmp_all)

# before/after bar for the all-samples comparison
a <- cmp[["all_samples"]]
long <- rbind(data.frame(cell_type = a$cell_type, when = "before (raw)",     n = a$n_before),
              data.frame(cell_type = a$cell_type, when = "after (cleaned)",  n = a$n_after))
p <- ggplot(long, aes(reorder(cell_type, n), n, fill = when)) +
    geom_col(position = "dodge") + coord_flip() +
    scale_fill_manual(values = c("before (raw)" = "grey60", "after (cleaned)" = "#2e86ab")) +
    labs(title = "All samples: cell-type counts before vs after cleaning",
         x = NULL, y = "cells") + theme_classic()
ggsave(file.path(outdir, "celltype_recovery_all_samples.png"), p, width = 9, height = 6, dpi = 200, bg = "white")

message("==== annot11 comparison done ====")
### END
