#!/bin/bash
# Running MBV to ensure good match between the sequence and genotype data

#SBATCH --job-name=mbv
#SBATCH --output=logs/%x-%j.log
#SBATCH --error=logs/%x-%j.err 
#SBATCH --time=12:00:00  
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=cpu-core-sponsored
#SBATCH --account=cpu-core-sponsored
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<email>
#SBATCH --array=1-321


# Go to WD
cd analyses/01_rna_genotype_matching

# Activate package environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate Env1


# Getting line from source file (each line is the path to a bam file) and defining the path to the input/output files
BAM=$(awk 'NR=='${SLURM_ARRAY_TASK_ID}' && NF' src/bam_aligned_to_genome.list) # Ignores empty lines
LIB_ID=$(echo $BAM | cut -f3 -d"/" | cut -f1-2 -d"_")
SAM_ID=$(awk -v id="$LIB_ID" '$1 == id {print $2}' src/lookup_table.tsv)
SAM_ID_FIXED=${SAM_ID//\//_}
QTLtools="QTLtools" # Path to QTLtools
VCF="results/genotype/mouse_geno_all_cln_sorted_maf0.05_geno0.05_RefAltOk_fixed.vcf.gz"
OUT="results/mbv_result/Mm_${LIB_ID}_${SAM_ID_FIXED}.mbv.bamstat.txt"

# Index bam
samtools index $BAM

# Run the mbv
$QTLtools   mbv \
            --bam $BAM \
            --vcf $VCF \
            --filter-mapping-quality 150 \
            --filter-base-quality 20 \
            --filter-minimal-coverage 10 \
            --out ${OUT}


# Plotting
Rscript scripts/06_Rscript03_plot_mbv_result.R --input $OUT --sample_id $SAM_ID --library_id $LIB_ID


### END


