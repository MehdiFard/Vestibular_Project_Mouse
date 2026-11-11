# recolor_composition.R  (seurat_env or Env1 — only needs ggplot2)
# Redraw the cleaned composition bar charts with the shared cell-type palette,
# reading the composition CSVs (no need to reload the Seurat objects).
suppressMessages(library(ggplot2))
source("scripts/celltype_palette.R")

dir <- "results/cell_annotation_clean"
for (nm in c("saccule", "utricle", "all")) {
    comp <- read.csv(file.path(dir, paste0(nm, "_clean_composition.csv")))
    ## keep only cell types actually present, in shared-palette order, so unused
    ## levels (e.g. HC_type_I/II in the combined run) don't show blank legend keys.
    present <- intersect(names(celltype_pal), unique(comp$cell_type))
    comp$cell_type <- factor(comp$cell_type, levels = present)
    p <- ggplot(comp, aes(sample, prop, fill = cell_type)) + geom_col() +
        facet_grid(~ tissue, scales = "free_x", space = "free_x") +
        scale_fill_manual(values = celltype_pal[present], drop = FALSE) +
        labs(y = "proportion of cells", x = NULL,
             title = paste0(nm, " (cleaned) — composition")) +
        theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(dir, paste0(nm, "_clean_composition_bar.png")),
           p, width = if (nm == "all") 8 else 6, height = 6, dpi = 200, bg = "white")
    cat("recolored bar:", nm, "\n")
}
### END
