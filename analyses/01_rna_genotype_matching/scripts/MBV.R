#!/usr/bin/env Rscript

# Set wd
# runs in analyses/01_rna_genotype_matching (the Snakefile's working directory)
rm(list = ls())

mbv_res_mm_602 <- read.table("results/mbv_result/test/Mm_FTL_602_BXD45_RwwJ.bamstat.txt", header = T)
View(mbv_res_mm_602)


# Define the true sample to highlight
true_sample_mm_602 <- "BXD45/RwwJ"
true_sample_index_mm_602 <- which(mbv_res_mm_602$SampleID == true_sample_mm_602)

# Plot the data with the specified y-axis range and suppressed axes
jpeg("results/mbv_result/Mm_FTL_602_BXD45_test.jpeg", width = 25, height = 25, units = "cm", quality = 100, res = 300)
plot(
    mbv_res_mm_602$n_het_consistent, 
    mbv_res_mm_602$n_hom_consistent, 
    xlab = "Concordance at hets",
    ylab = "Concordance at homs", 
    main = "Mm_FTL_602_BXD45",
    axes = FALSE, 
    ylim = c(3400, 4500),
    pch = 16,
    cex = 1
)
# Add custom axes
axis(side = 1)
axis(side = 2)

# Highlight the specific sample
points(
    mbv_res_mm_602$n_het_consistent[true_sample_index_mm_602], 
    mbv_res_mm_602$n_hom_consistent[true_sample_index_mm_602], 
    col = "red", 
    pch = 16, # Filled circle
    cex = 1.5 # Larger size for emphasis
)

# Add the sample name as a label above the highlighted point
text(
    mbv_res_mm_602$n_het_consistent[true_sample_index_mm_602] + 1, 
    mbv_res_mm_602$n_hom_consistent[true_sample_index_mm_602] + 50, # Offset above the point
    labels = mbv_res_mm_602$SampleID[true_sample_index_mm_602], 
    col = "red"
)
dev.off()


### FTL_872
mbv_res_mm_872 <- read.table("results/mbv_result/Mm_FTL_872_BXD77_RwwJ.bamstat.txt", header = T)
View(mbv_res_mm_872)


# Define the true sample to highlight
true_sample_mm_872 <- "BXD77/RwwJ"
true_sample_index_mm_872 <- which(mbv_res_mm_872$SampleID == true_sample_mm_872)

# Plot the data with the specified y-axis range and suppressed axes
jpeg("results/mbv_result/Mm_FTL_872_BXD77_test.jpeg", width = 25, height = 25, units = "cm", quality = 100, res = 300)
plot(
    mbv_res_mm_872$n_het_consistent, 
    mbv_res_mm_872$n_hom_consistent, 
    xlab = "Concordance at hets",
    ylab = "Concordance at homs", 
    main = "Mm_FTL_872_BXD77",
    axes = FALSE, 
    ylim = c(5500, 7500),
    pch = 16,
    cex = 1
)
# Add custom axes
axis(side = 1)
axis(side = 2)

# Highlight the specific sample
points(
    mbv_res_mm_872$n_het_consistent[true_sample_index_mm_872], 
    mbv_res_mm_872$n_hom_consistent[true_sample_index_mm_872], 
    col = "red", 
    pch = 16, # Filled circle
    cex = 1.5 # Larger size for emphasis
)

# Add the sample name as a label above the highlighted point
text(
    mbv_res_mm_872$n_het_consistent[true_sample_index_mm_872], 
    mbv_res_mm_872$n_hom_consistent[true_sample_index_mm_872] + 100, # Offset above the point
    labels = mbv_res_mm_872$SampleID[true_sample_index_mm_872], 
    col = "red"
)

# Highlighting BXD75
bxd75_index = which(mbv_res_mm_872$SampleID == "BXD75/RwwJ")
points(
    mbv_res_mm_872$n_het_consistent[bxd75_index], 
    mbv_res_mm_872$n_hom_consistent[bxd75_index], 
    col = "chartreuse4", 
    pch = 16, # Filled circle
    cex = 1.5 # Larger size for emphasis
)

# Add the sample name as a label above the highlighted point
text(
    mbv_res_mm_872$n_het_consistent[bxd75_index], 
    mbv_res_mm_872$n_hom_consistent[bxd75_index] + 100, # Offset above the point
    labels = mbv_res_mm_872$SampleID[bxd75_index], 
    col = "chartreuse4"
)
dev.off()


###----------- Human data
# FTL_817
mbv_res_hs_817 <- read.table("results/mbv_result/Hs_FTL_817_FTL_32287621.bamstat.txt", header = T)
mbv_res_hs_817$SampleID <- gsub("FTL", "FTL_", mbv_res_hs_817$SampleID)
View(mbv_res_hs_817)

# Define the true sample to highlight
true_sample_hs_817 <- "FTL_32287621"
true_sample_index_hs_817 <- which(mbv_res_hs_817$SampleID == true_sample_hs_817)

# Plot the data with the specified y-axis range and suppressed axes
jpeg("results/mbv_result/Hs_FTL_817_FTL_32287621_test.jpeg", width = 25, height = 25, units = "cm", quality = 100, res = 300)
plot(
    mbv_res_hs_817$n_het_consistent, 
    mbv_res_hs_817$n_hom_consistent, 
    xlab = "Concordance at hets",
    ylab = "Concordance at homs", 
    main = "Hs_FTL_817_FTL_32287621",
    axes = F, 
    xlim = c(20000, 65000),
    pch = 16,
    cex = 1
)
# Add custom axes
axis(side = 1)
axis(side = 2)

# Highlight the specific sample
points(
    mbv_res_hs_817$n_het_consistent[true_sample_index_hs_817], 
    mbv_res_hs_817$n_hom_consistent[true_sample_index_hs_817], 
    col = "red", 
    pch = 16, # Filled circle
    cex = 1.5 # Larger size for emphasis
)

# Add the sample name as a label above the highlighted point
text(
    mbv_res_hs_817$n_het_consistent[true_sample_index_hs_817] - 1000, 
    mbv_res_hs_817$n_hom_consistent[true_sample_index_hs_817] + 1000, # Offset above the point
    labels = mbv_res_hs_817$SampleID[true_sample_index_hs_817], 
    col = "red"
)
dev.off()


###-------------------------- FTL_844 - Plotting using number of concordance
mbv_res_hs_844 <- read.table("results/mbv_result/test/Hs_FTL_844_FTL_23884570.bamstat.txt", header = T)
mbv_res_hs_844$SampleID <- gsub("FTL", "FTL_", mbv_res_hs_844$SampleID)
View(mbv_res_hs_844)

# Define the true sample to highlight
true_sample_hs_844 <- "FTL_23884570"
true_sample_index_hs_844 <- which(mbv_res_hs_844$SampleID == true_sample_hs_844)

# Plot the data with the specified y-axis range and suppressed axes
jpeg("results/mbv_result/test/Hs_FTL_844_FTL_23884570_test.jpeg", width = 30, height = 30, units = "cm", quality = 100, res = 300)
plot(
    mbv_res_hs_844$n_het_consistent, 
    mbv_res_hs_844$n_hom_consistent, 
    xlab = "Concordance at hets",
    ylab = "Concordance at homs", 
    main = "Hs_FTL_844_FTL_23884570",
    axes = F, 
    # xlim = c(20000, 40000),
    # ylim = c(120000, 140000),
    pch = 16,
    cex = 1
)
# Add custom axes
axis(side = 1)
axis(side = 2)

# Highlight the specific sample
points(
    mbv_res_hs_844$n_het_consistent[true_sample_index_hs_844], 
    mbv_res_hs_844$n_hom_consistent[true_sample_index_hs_844], 
    col = "red", 
    pch = 16, # Filled circle
    cex = 1.5 # Larger size for emphasis
)

# Add the sample name as a label above the highlighted point
text(
    mbv_res_hs_844$n_het_consistent[true_sample_index_hs_844] + 50, 
    mbv_res_hs_844$n_hom_consistent[true_sample_index_hs_844] + 150, # Offset above the point
    labels = mbv_res_hs_844$SampleID[true_sample_index_hs_844], 
    col = "red"
)

# Highlight the specific sample
ase_matched_index <- which(mbv_res_hs_844$SampleID == "FTL_31819005")


points(
    mbv_res_hs_844$n_het_consistent[ase_matched_index], 
    mbv_res_hs_844$n_hom_consistent[ase_matched_index], 
    col = "chartreuse4", 
    pch = 16, # Filled circle
    cex = 1.5 # Larger size for emphasis
)

# Add the sample name as a label above the highlighted point
text(
    mbv_res_hs_844$n_het_consistent[ase_matched_index], 
    mbv_res_hs_844$n_hom_consistent[ase_matched_index] + 150, # Offset above the point
    labels = mbv_res_hs_844$SampleID[ase_matched_index], 
    col = "chartreuse4"
)

# Highlighting the sample with highest concordance 
# max_val_het <- max(mbv_res_hs_844$n_het_consistent)
# max_val_hom <- max(mbv_res_hs_844$n_hom_consistent)
max_val_het = mbv_res_hs_844[which(mbv_res_hs_844$perc_het_consistent == max(mbv_res_hs_844$perc_het_consistent)), "n_het_consistent"]
max_val_hom = mbv_res_hs_844[which(mbv_res_hs_844$perc_hom_consistent == max(mbv_res_hs_844$perc_hom_consistent)), "n_hom_consistent"]


points(
    max_val_het, 
    max_val_hom, 
    col = "blue",
    pch = 16, # Filled circle
    cex = 1.5 # Larger size for emphasis
)

text(
    max_val_het - 150, 
    max_val_hom + 150, # Offset above the point
    labels = mbv_res_hs_844$SampleID[which(mbv_res_hs_844$n_het_consistent == max_val_het)], 
    col = "blue"
)
dev.off()


###-------------------------- FTL_844 - Plotting using percentage of concordance to covered
mbv_res_hs_844 <- read.table("results/mbv_result/test/Hs_FTL_844_FTL_23884570.bamstat.txt", header = T)
mbv_res_hs_844$SampleID <- gsub("FTL", "FTL_", mbv_res_hs_844$SampleID)
View(mbv_res_hs_844)

# Define the true sample to highlight
true_sample_hs_844 <- "FTL_23884570"
true_sample_index_hs_844 <- which(mbv_res_hs_844$SampleID == true_sample_hs_844)

# Plotting
jpeg("results/mbv_result/test/Hs_FTL_844_FTL_23884570_perc_extended_axis_test.jpeg", width = 30, height = 30, units = "cm", quality = 100, res = 300)
plot(
    mbv_res_hs_844$perc_het_consistent, 
    mbv_res_hs_844$perc_hom_consistent, 
    xlab = "Ratio of concordance/covered at hets",
    ylab = "Ratio of concordance/covered at homs", 
    main = "Hs_FTL_844_FTL_23884570",
    axes = F, 
    xlim = c(0.3, 0.5),
    ylim = c(0.7, 0.8),
    pch = 16,
    cex = 1
)
# Add custom axes
axis(side = 1)
axis(side = 2)

# Highlight the specific sample
points(
    mbv_res_hs_844$perc_het_consistent[true_sample_index_hs_844], 
    mbv_res_hs_844$perc_hom_consistent[true_sample_index_hs_844], 
    col = "red", 
    pch = 16, # Filled circle
    cex = 1.5 # Larger size for emphasis
)

# Add the sample name as a label above the highlighted point
text(
    mbv_res_hs_844$perc_het_consistent[true_sample_index_hs_844] - 0.001, 
    mbv_res_hs_844$perc_hom_consistent[true_sample_index_hs_844] + 0.0003, # Offset above the point
    labels = mbv_res_hs_844$SampleID[true_sample_index_hs_844], 
    col = "red"
)

# Highlight the specific sample
ase_matched_index <- which(mbv_res_hs_844$SampleID == "FTL_31819005")


points(
    mbv_res_hs_844$perc_het_consistent[ase_matched_index], 
    mbv_res_hs_844$perc_hom_consistent[ase_matched_index], 
    col = "chartreuse4", 
    pch = 16, # Filled circle
    cex = 1.5 # Larger size for emphasis
)

# Add the sample name as a label above the highlighted point
text(
    mbv_res_hs_844$perc_het_consistent[ase_matched_index], 
    mbv_res_hs_844$perc_hom_consistent[ase_matched_index] + 0.0005, # Offset above the point
    labels = mbv_res_hs_844$SampleID[ase_matched_index], 
    col = "chartreuse4"
)

# Highlighting the sample with highest concordance 
max_val_het <- max(mbv_res_hs_844$perc_het_consistent)
max_val_hom <- max(mbv_res_hs_844$perc_hom_consistent)

points(
    max_val_het, 
    max_val_hom, 
    col = "blue",
    pch = 16, # Filled circle
    cex = 1.5 # Larger size for emphasis
)

text(
    max_val_het - 0.002, 
    max_val_hom + 0.0005, # Offset above the point
    labels = mbv_res_hs_844$SampleID[which(mbv_res_hs_844$perc_het_consistent == max_val_het)], 
    col = "blue"
)
dev.off()



###-------------------------- FTL_834
mbv_res_hs_834 <- read.table("results/mbv_result/Hs_FTL_834_FTL_17465303.bamstat.txt", header = T)
mbv_res_hs_834$SampleID <- gsub("FTL", "FTL_", mbv_res_hs_834$SampleID)
View(mbv_res_hs_834)

# Define the true sample to highlight
true_sample_hs_834 <- "FTL_17465303"
true_sample_index_hs_834 <- which(mbv_res_hs_834$SampleID == true_sample_hs_834)

# Plot the data with the specified y-axis range and suppressed axes
jpeg("results/mbv_result/Hs_FTL_834_FTL_17465303_test.jpeg", width = 25, height = 25, units = "cm", quality = 100, res = 300)
plot(
    mbv_res_hs_834$n_het_consistent, 
    mbv_res_hs_834$n_hom_consistent, 
    xlab = "Concordance at hets",
    ylab = "Concordance at homs", 
    main = "Hs_FTL_834_FTL_17465303",
    axes = F, 
    ylim = c(120000, 180000),
    xlim = c(60000, 115000),
    pch = 16,
    cex = 1
)
# Add custom axes
axis(side = 1)
axis(side = 2)

# Highlight the specific sample
points(
    mbv_res_hs_834$n_het_consistent[true_sample_index_hs_834], 
    mbv_res_hs_834$n_hom_consistent[true_sample_index_hs_834], 
    col = "red", 
    pch = 16, # Filled circle
    cex = 1.5 # Larger size for emphasis
)

# Add the sample name as a label above the highlighted point
text(
    mbv_res_hs_834$n_het_consistent[true_sample_index_hs_834], 
    mbv_res_hs_834$n_hom_consistent[true_sample_index_hs_834] + 1000, # Offset above the point
    labels = mbv_res_hs_834$SampleID[true_sample_index_hs_834], 
    col = "red"
)

# Highlighting the sample with highest concordance 
max_val_het <- max(mbv_res_hs_834$n_het_consistent)
max_val_hom <- max(mbv_res_hs_834$n_hom_consistent)

points(
    max_val_het, 
    max_val_hom, 
    col = "chartreuse4",
    pch = 16, # Filled circle
    cex = 1.5 # Larger size for emphasis
)

text(
    max_val_het, 
    max_val_hom + 1000, # Offset above the point
    labels = mbv_res_hs_834$SampleID[which(mbv_res_hs_834$n_het_consistent == max_val_het)], 
    col = "chartreuse4"
)

dev.off()

###-------------------------- FTL_551 - Plotting using number of concordance
mbv_res_hs_551 <- read.table("results/mbv_result/Hs_FTL_551_FTL_31907773.bamstat.txt", header = T)
mbv_res_hs_551$SampleID <- gsub("FTL", "FTL_", mbv_res_hs_551$SampleID)
View(mbv_res_hs_551)

# Define the true sample to highlight
true_sample_hs_551 <- "FTL_31907773"
true_sample_index_hs_551 <- which(mbv_res_hs_551$SampleID == true_sample_hs_551)

# Plot the data with the specified y-axis range and suppressed axes
jpeg("results/mbv_result/Hs_FTL_551_FTL_31907773_test.jpeg", width = 25, height = 25, units = "cm", quality = 100, res = 300)
plot(
    mbv_res_hs_551$n_het_consistent, 
    mbv_res_hs_551$n_hom_consistent, 
    xlab = "Concordance at hets",
    ylab = "Concordance at homs", 
    main = "Hs_FTL_551_FTL_31907773",
    axes = F, 
    ylim = c(35000, 55000),
    xlim = c(10000, 30000),
    pch = 16,
    cex = 1
)
# Add custom axes
axis(side = 1)
axis(side = 2)

# Highlight the specific sample
points(
    mbv_res_hs_551$n_het_consistent[true_sample_index_hs_551], 
    mbv_res_hs_551$n_hom_consistent[true_sample_index_hs_551], 
    col = "red", 
    pch = 16, # Filled circle
    cex = 1.5 # Larger size for emphasis
)

# Add the sample name as a label above the highlighted point
text(
    mbv_res_hs_551$n_het_consistent[true_sample_index_hs_551] + 2200, 
    mbv_res_hs_551$n_hom_consistent[true_sample_index_hs_551] + 500, # Offset above the point
    labels = mbv_res_hs_551$SampleID[true_sample_index_hs_551], 
    col = "red"
)

# Highlighting the sample with highest concordance 
max_val_het <- max(mbv_res_hs_551$n_het_consistent)
max_val_hom <- max(mbv_res_hs_551$n_hom_consistent)

points(
    max_val_het, 
    max_val_hom, 
    col = "chartreuse4",
    pch = 16, # Filled circle
    cex = 1.5 # Larger size for emphasis
)

text(
    max_val_het, 
    max_val_hom + 500, # Offset above the point
    labels = mbv_res_hs_551$SampleID[which(mbv_res_hs_551$n_het_consistent == max_val_het)], 
    col = "chartreuse4"
)

dev.off()


###-------------------------- FTL_551 - Plotting using percentage of concordance to covered
mbv_res_hs_551 <- read.table("results/mbv_result/test/Hs_FTL_551_FTL_31907773.bamstat.txt", header = T)
mbv_res_hs_551$SampleID <- gsub("FTL", "FTL_", mbv_res_hs_551$SampleID)
View(mbv_res_hs_551)

# Define the true sample to highlight
true_sample_hs_551 <- "FTL_31907773"
true_sample_index_hs_551 <- which(mbv_res_hs_551$SampleID == true_sample_hs_551)

# Plot the data with the specified y-axis range and suppressed axes
jpeg("results/mbv_result/test/Hs_FTL_551_FTL_31907773_perc_test.jpeg", width = 25, height = 25, units = "cm", quality = 100, res = 300)
plot(
    mbv_res_hs_551$perc_het_consistent, 
    mbv_res_hs_551$perc_hom_consistent, 
    xlab = "Concordance at hets",
    ylab = "Concordance at homs", 
    main = "Hs_FTL_551_FTL_31907773",
    axes = F, 
    # ylim = c(35000, 55000),
    # xlim = c(10000, 30000),
    pch = 16,
    cex = 1
)
# Add custom axes
axis(side = 1)
axis(side = 2)

# Highlight the specific sample
points(
    mbv_res_hs_551$perc_het_consistent[true_sample_index_hs_551], 
    mbv_res_hs_551$perc_hom_consistent[true_sample_index_hs_551], 
    col = "red", 
    pch = 16, # Filled circle
    cex = 1.5 # Larger size for emphasis
)

# Add the sample name as a label above the highlighted point
text(
    mbv_res_hs_551$perc_het_consistent[true_sample_index_hs_551], 
    mbv_res_hs_551$perc_hom_consistent[true_sample_index_hs_551] + 0.01, # Offset above the point
    labels = mbv_res_hs_551$SampleID[true_sample_index_hs_551], 
    col = "red"
)

# Highlighting the sample with highest concordance 
max_val_het <- max(mbv_res_hs_551$perc_het_consistent)
max_val_hom <- max(mbv_res_hs_551$perc_hom_consistent)

points(
    max_val_het, 
    max_val_hom, 
    col = "chartreuse4",
    pch = 16, # Filled circle
    cex = 1.5 # Larger size for emphasis
)

text(
    max_val_het - 0.02, 
    max_val_hom + 0.005, # Offset above the point
    labels = mbv_res_hs_551$SampleID[which(mbv_res_hs_551$perc_het_consistent == max_val_het)], 
    col = "chartreuse4"
)
dev.off()


###------ FTL_940
mbv_res_mm_940 <- read.table("results/mbv_result/Mm_FTL_940_BXD224_RwwJ.mbv.bamstat.txt", header = T)
View(mbv_res_mm_940)


# Define the true sample to highlight
true_sample_mm_940 <- "BXD224/RwwJ"
true_sample_index_mm_940 <- which(mbv_res_mm_940$SampleID == true_sample_mm_940)

# Plot the data with the specified y-axis range and suppressed axes
jpeg("results/mbv_result/test/Mm_FTL_940_BXD224_test.jpeg", width = 25, height = 25, units = "cm", quality = 100, res = 300)
plot(
    mbv_res_mm_940$n_het_consistent, 
    mbv_res_mm_940$n_hom_consistent, 
    xlab = "Concordance at hets",
    ylab = "Concordance at homs", 
    main = "Mm_FTL_940_BXD224",
    axes = FALSE, 
    # ylim = c(5500, 7500),
    pch = 16,
    cex = 1
)
# Add custom axes
axis(side = 1)
axis(side = 2)

# Highlight the specific sample
points(
    mbv_res_mm_940$n_het_consistent[true_sample_index_mm_940], 
    mbv_res_mm_940$n_hom_consistent[true_sample_index_mm_940], 
    col = "red", 
    pch = 16, # Filled circle
    cex = 1.5 # Larger size for emphasis
)

# Add the sample name as a label above the highlighted point
text(
    mbv_res_mm_940$n_het_consistent[true_sample_index_mm_940], 
    mbv_res_mm_940$n_hom_consistent[true_sample_index_mm_940] + 25, # Offset above the point
    labels = mbv_res_mm_940$SampleID[true_sample_index_mm_940], 
    col = "red"
)

# Highlighting BXD16
bxd16_index = which(mbv_res_mm_940$SampleID == "BXD16/TyJ")
points(
    mbv_res_mm_940$n_het_consistent[bxd16_index], 
    mbv_res_mm_940$n_hom_consistent[bxd16_index], 
    col = "chartreuse4", 
    pch = 16, # Filled circle
    cex = 1.5 # Larger size for emphasis
)

# Add the sample name as a label above the highlighted point
text(
    mbv_res_mm_940$n_het_consistent[bxd16_index] + 2, 
    mbv_res_mm_940$n_hom_consistent[bxd16_index] + 20, # Offset above the point
    labels = mbv_res_mm_940$SampleID[bxd16_index], 
    col = "chartreuse4"
)

# Highlighting BXD12
bxd12_index = which(mbv_res_mm_940$SampleID == "BXD12/TyJ")
points(
    mbv_res_mm_940$n_het_consistent[bxd12_index], 
    mbv_res_mm_940$n_hom_consistent[bxd12_index], 
    col = "blue", 
    pch = 16, # Filled circle
    cex = 1.5 # Larger size for emphasis
)

# Add the sample name as a label above the highlighted point
text(
    mbv_res_mm_940$n_het_consistent[bxd12_index] + 1, 
    mbv_res_mm_940$n_hom_consistent[bxd12_index] + 40, # Offset above the point
    labels = mbv_res_mm_940$SampleID[bxd12_index], 
    col = "blue"
)

dev.off()



