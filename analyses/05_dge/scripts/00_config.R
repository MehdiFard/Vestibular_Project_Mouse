# Shared paths and settings for the saccule-vs-utricle analysis.
# Every other script starts by sourcing this file.

PROJ    <- "../.."
PANTRY  <- file.path(PROJ, "analyses", "03_pantry")
BASE    <- file.path(PROJ, "analyses", "05_dge")

DATA    <- file.path(BASE, "data")
TABLES  <- file.path(BASE, "results", "tables")
FIGS    <- file.path(BASE, "results", "figures")

# Sample metadata (strain, sex, plate), copied from the 2024 project
# (cybertron_data/project01_vestibular/mouse/vestibular_mouse/src/)
REF_META <- file.path(BASE, "src", "metadata_mouse_all_plates_with_matched_geno_ids.csv")

# RNA-to-genotype concordance results (ASEReadCounter vs the strain VCF).
# This is the same table Pantry itself used to attach genotypes to libraries,
# so it is the authoritative strain label for this project - see README.
GENO_MATCH <- file.path(PROJ, "analyses", "01_rna_genotype_matching", "results", "ase_counter_result",
                        "rna_to_geno_summary.tsv")

GTF <- file.path(PROJ, "resources", "genome", "Mus_musculus.GRCm38.102.chr.gtf")

# Submission workbook (batches 1-4, one sheet each): library, strain, tissue, sex, age. The two organs of one animal
# were submitted as consecutive library IDs, which is how organ pairs are identified (01_build_inputs.R).
SUBMISSION_SHEET <- file.path(PROJ, "data", "fastq_bulk", "bulk_rnaseq_metadata_mouse_b1_to_b4.xlsx")

# The four Pantry datasets. Organ and age are read straight from these names.
DATASETS <- c("saccule_w8", "saccule_w22", "utricle_w8", "utricle_w22")

# Gene biotypes to analyse. Everything else (pseudogenes, TEC, miRNA, snoRNA,
# antisense and the other small-RNA classes) is dropped before the expression
# filter, so it never enters the multiple-testing burden.
# Note: this GTF is Ensembl 102, where lncRNAs are still split into several
# biotypes. "lincRNA" here is the intergenic class only; add "antisense",
# "bidirectional_promoter_lncRNA", "sense_intronic" and "processed_transcript"
# if you want lncRNA in the broad sense.
KEEP_BIOTYPES <- c("protein_coding", "lincRNA")

# Significance thresholds used throughout
FDR_CUT  <- 0.05
LFC_CUT  <- 1      # |log2FC| > 1, i.e. a two-fold change

# Consistent colours for the four organ-age groups
ORGAN_COLS <- c(Saccule = "#1B7837", Utricle = "#762A83")
GROUP_COLS <- c("Saccule w8"  = "#7FBF7B", "Saccule w22" = "#1B7837",
                "Utricle w8"  = "#C2A5CF", "Utricle w22" = "#762A83")
SEX_COLS   <- c(F = "#D6604D", M = "#4393C3")

dir.create(DATA,   showWarnings = FALSE, recursive = TRUE)
dir.create(TABLES, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGS,   showWarnings = FALSE, recursive = TRUE)
