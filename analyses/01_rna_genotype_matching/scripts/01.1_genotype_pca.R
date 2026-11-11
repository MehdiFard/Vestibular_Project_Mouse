#!/usr/bin/env Rscript

# Set wd
# runs in analyses/01_rna_genotype_matching (the Snakefile's working directory)
rm(list = ls())

### Running pca on all samples
system("plink2 --tfile results/genotype/mouse_geno_all --maf 0.05 --geno 0.05 --pca 10 --out intermediate/pca_geno/pca_all_strains")

# PC plot
pca_eigenvec <- data.table::fread("intermediate/pca_geno/pca_all_strains.eigenvec")
pca_eigenval <- data.table::fread("intermediate/pca_geno/pca_all_strains.eigenval")$V1
pca_eigenval_percent <- pca_eigenval / sum(pca_eigenval) * 100

metadata <- data.table::fread("src/metadata_mouse_all_plates_with_matched_geno_ids.tsv")
dim(metadata)
length(unique(metadata$geno_id_original))
table(unique(metadata$geno_id_original) %in% pca_eigenvec$IID)

library(dplyr)
# Add a new column indicating if IID is a study sample
pca_eigenvec <- pca_eigenvec %>%
    mutate(is_study_sample = ifelse(IID %in% unique(metadata$geno_id_original), "YES", "NO"))
table(pca_eigenvec$is_study_sample)

plot1_pca = ggplot(data = pca_eigenvec, aes(x = PC1, y = PC2, color = is_study_sample), show.legend = T) +
    geom_point(size = 2) + # Adjusts dot size in the plot
    geom_hline(yintercept = 0, linetype = "dotted", col = "black") +
    geom_vline(xintercept = 0, linetype = "dotted", col = "black") +
    labs(
        title = "Genotype PCA - PC1 vs PC2",
        x = paste0("PC 1 (", round(pca_eigenval_percent[1], 2), " %)"),
        y = paste0("PC 2 (", round(pca_eigenval_percent[2], 2), " %)"),
        color = "Study sample") + 
    theme_classic() +
    theme(
        axis.line = element_blank(), # Removes the axis lines
        axis.ticks = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 2))) # Adjusts legend dot size
plot1_pca

ggsave(filename = "00.1_geno_pc1_vs_pc2_scatter.jpeg", plot =  plot1_pca,
        units = "in", path = "figures/", dpi = 600)

### 3D plot using plotly (fully interactive)
library(plotly)
library(reticulate)
x <- pca_eigenvec$PC1
y <- pca_eigenvec$PC2
z <- pca_eigenvec$PC3

plot2_pca_3d <- plot_ly(
    x = ~x,
    y = ~y,
    z = ~z,
    type = 'scatter3d',
    mode = 'markers',
    marker = list(size = 3)) %>%
    layout(title = "Genotype PCA",
        scene = list(xaxis = list(title = paste0("PC1 (", round(pca_eigenval_percent[1], 2), " %)")),
                    yaxis = list(title = paste0("PC2 (", round(pca_eigenval_percent[2], 2), " %)")),
                    zaxis = list(title = paste0("PC3 (", round(pca_eigenval_percent[3], 2), " %)"))))
plot2_pca_3d

# Save plot - jpeg
plotly::save_image(plot2_pca_3d, "figures/00.2_geno_pc1_vs_pc2_vs_pc3_3D.jpeg", size = 0.5)

# Save plot - HTML
htmlwidgets::saveWidget(as_widget(plot2_pca_3d), "figures/00.2_geno_pc1_vs_pc2_vs_pc3_3D.html")



###------ Other Packages for 3D plots
## 3D plot using lattice (not interactive)
library(lattice)
# Create 3D scatter plot
cloud(z ~ x * y, pch = 16, col = "blue", main = "3D Scatter Plot",
        xlab = "X", ylab = "Y", zlab = "Z")


## 3D plot using ggplot2 and plotly (interactive, can't move the plot but can select dots and ...)
library(ggplot2)
library(plotly)
# Sample data
df <- data.frame(x = pca_eigenvec$PC1, y = pca_eigenvec$PC2, z = pca_eigenvec$PC3)
# ggplot2 3D scatter plot with plotly
p <- ggplot(df, aes(x = x, y = y, z = z)) +
    geom_point() +
    labs(title = "3D Scatter Plot")
ggplotly(p)
