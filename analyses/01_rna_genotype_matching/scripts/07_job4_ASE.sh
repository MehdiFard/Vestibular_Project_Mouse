#!/bin/bash
# Running ASEReadCounter to ensure good match between the sequence and genotype data

#SBATCH --job-name=ase
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


# Defining inputs and outputs
REF="../../resources/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa"
VCF_IN="results/genotype/snps_for_mixup_test.vcf.gz"
OUT="results/ase_counter_result/"

# Getting path to BAM file
BAM=$(awk 'NR=='${SLURM_ARRAY_TASK_ID}' && NF' src/bam_aligned_to_genome.list) # Ignores empty lines
LIB_ID=$(echo $BAM | cut -f3 -d"/" | cut -f1-2 -d"_")


# Running ASE
gatk ASEReadCounter \
            -R $REF \
            -I $BAM \
            -V $VCF_IN \
            -O ${OUT}${LIB_ID}.readcounts.txt \
            --verbosity ERROR

### END