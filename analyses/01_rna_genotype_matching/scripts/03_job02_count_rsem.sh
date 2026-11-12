#!/bin/bash

#SBATCH --job-name=count_rsem
#SBATCH --output=logs/%x-%j.log
#SBATCH --error=logs/%x-%j.err 
#SBATCH --time=12:00:00  
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=64G
#SBATCH --partition=cpu-core-sponsored
#SBATCH --account=cpu-core-sponsored
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<email>
#SBATCH --array=1-323

###---------------- Count - RSEM
# WD
cd analyses/01_rna_genotype_matching

# Activate package environment
source "$(conda info --base)/etc/profile.d/conda.sh"
# source $(conda info --base)/etc/profile.d/mamba.sh
conda activate Env1

# Getting file path and library ids for each job
BAM_FILE=$(awk 'NR=='${SLURM_ARRAY_TASK_ID}' && NF' src/bam_aligned_to_transcriptome.list) # Ignores empty lines
LIBRARY_ID=$(echo "$BAM_FILE" | cut -f3 -d"/" | cut -f1-2 -d"_")

# Path to Ref and fai index
REF="../../resources/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa"
REF_FAI_INDX="../../resources/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.fai"
INDX_RSEM="../../resources/index_rsem/rsem_indx_m38"
OUT="results/count_rsem/"
mkdir -p $OUT

# echo " "
# echo "########## Analyzing ${LIBRARY_ID} #########"
# echo " "

# RSEM
rsem-calculate-expression \
        --num-threads 14 \
        --quiet \
        --no-bam-output \
        --estimate-rspd \
        --alignments ${BAM_FILE} \
        --bam \
        --paired-end \
        --fai ${REF_FAI_INDX} \
        ${INDX_RSEM} \
        ${OUT}${LIBRARY_ID}

gzip $OUT/${LIBRARY_ID}.genes.results
gzip $OUT/${LIBRARY_ID}.isoforms.results

### END