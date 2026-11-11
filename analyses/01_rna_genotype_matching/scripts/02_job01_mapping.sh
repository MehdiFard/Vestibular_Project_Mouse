#!/bin/bash

#SBATCH --job-name=mapping
#SBATCH --output=logs/%x-%j.log
#SBATCH --error=logs/%x-%j.err 
#SBATCH --time=12:00:00  
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=128G
#SBATCH --partition=cpu-core-sponsored
#SBATCH --account=cpu-core-sponsored
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<email>
#SBATCH --array=1-324

###---------------- Mapping - STAR
# Go to WD
cd analyses/01_rna_genotype_matching

# Activate package environment
source "$(conda info --base)/etc/profile.d/conda.sh"
# source $(conda info --base)/etc/profile.d/mamba.sh
conda activate Env1

# Defining input/output
REF="../../resources/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa"
INDEX="/path/to/refs/index_star/"
STAR="STAR"
OUT_DIR="intermediate/bams/"
mkdir -p ${OUT_DIR}

# Getting line from source file - ech line is a library id
FASTQ=$(awk 'NR=='${SLURM_ARRAY_TASK_ID}' && NF' src/fastqs_all.list) # Ignores empty lines
PASS="../../data/fastq_bulk/"   # batch1-4/


# Making path to triimed reads
R1=$(ls ${PASS}batch*/${FASTQ}_R1_001.fastq.gz)
R2=$(ls ${PASS}batch*/${FASTQ}_R2_001.fastq.gz)

# getting required info for RG
LIBRARY_ID=$(echo $FASTQ | cut -f1-2 -d"_")
SAMPLE_ID=$(awk -v search="$LIBRARY_ID" -F'\t' '$1 == search { print $2 }' src/lib_sam_ins_4RG.tsv)
INSERT_SIZE=$(awk -v search="$LIBRARY_ID" -F'\t' '$1 == search { print $3 }' src/lib_sam_ins_4RG.tsv)

# getting identifier reads from fastq file and extracting the info for RG
IDENTIFIRE=$(head -n 1 <(zcat ${R1}))
INSTRO_NAME=$(echo $IDENTIFIRE | awk -F ':' '{print $1}' | cut -f 2 -d " " | sed 's/@//')
RUN_ID=$(echo $IDENTIFIRE | awk -F ':' '{print $2}')
FLOWCELL=$(echo $IDENTIFIRE | awk -F ':' '{print $3}')
LANE=$(echo $IDENTIFIRE | awk -F ':' '{print $4}')
TILE=$(echo $IDENTIFIRE | awk -F ':' '{print $5}')

RGID=$(echo $FLOWCELL.$LANE)
RGSM=$SAMPLE_ID
RGLB=$LIBRARY_ID
RGPU=$(echo $FLOWCELL.$LANE.$TILE)
RGPL="ILLUMINA"

# echo " "
# echo "########## Analyzing ${LIBRARY_ID} ##########"
# echo " "

$STAR   --runMode alignReads \
        --runThreadN 14 \
        --genomeDir ${INDEX} \
        --readFilesIn ${R1} ${R2} \
        --readFilesCommand zcat \
        --twopassMode Basic \
        --quantMode TranscriptomeSAM GeneCounts \
        --outSAMstrandField intronMotif \
        --outSAMattrRGline ID:$RGID SM:$RGSM LB:$RGLB PL:$RGPL PU:$RGPU \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outFileNamePrefix ${OUT_DIR}${LIBRARY_ID}_

# --outSAMmapqUnique 60 \ # If we need to change the MAPQ max defulte from 255 to 60

### END