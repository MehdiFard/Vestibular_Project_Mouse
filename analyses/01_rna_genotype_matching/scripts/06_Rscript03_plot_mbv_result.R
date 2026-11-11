#!/usr/bin/env Rscript
# Plotting the MBV output

# Set wd
# runs in analyses/01_rna_genotype_matching (the Snakefile's working directory)
rm(list = ls())
library(optparse)

# Define options
option_list <- list(
    make_option("--input", type = "character", help = "MBV output file"),
    make_option("--sample_id", type = "character", help = "Sample ID"),
    make_option("--library_id", type = "character", help = "Library ID")
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Print parsed options
cat("Input file:", opt$input, "\n")
cat("Sample ID:", opt$sample_id, "\n")
cat("Library ID:", opt$library_id, "\n")

# Reading input data
mbv_res <- read.table(opt$input, header = T)

# Define the true sample to highlight
true_sample <- opt$sample_id
true_sample_index <- which(mbv_res$SampleID == true_sample)

# Plot the data with the specified y-axis range and suppressed axes
sample_id_edited <- gsub("/", "_", opt$sample_id) # Replacing "/" with "_"
jpeg(paste0("figures/mbv/Mm_",opt$library_id, "_", sample_id_edited, ".jpeg"), width = 25, height = 25, units = "cm", quality = 100, res = 300)
plot(
    mbv_res$n_het_consistent, 
    mbv_res$n_hom_consistent, 
    xlab = "Concordance at hets",
    ylab = "Concordance at homs", 
    main = paste0("Mm_", opt$library_id, "_", opt$sample_id),
    axes = FALSE, 
    # ylim = c(3400, 4500),
    # xlim = c()
    pch = 16,
    cex = 1
)
# Add custom axes
axis(side = 1)
axis(side = 2)

# Highlight the specific sample
points(
    mbv_res$n_het_consistent[true_sample_index], 
    mbv_res$n_hom_consistent[true_sample_index], 
    col = "red", 
    pch = 16, # Filled circle
    cex = 1.5 # Larger size for emphasis
)

# # Add the sample name as a label above the highlighted point
# text(
#     mbv_res$n_het_consistent[true_sample_index], 
#     mbv_res$n_hom_consistent[true_sample_index] + 100, # Offset above the point
#     labels = mbv_res$SampleID[true_sample_index], 
#     col = "red"
# )
dev.off()

### END