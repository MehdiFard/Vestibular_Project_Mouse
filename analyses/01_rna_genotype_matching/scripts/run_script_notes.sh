#!/bin/bash

ssh user@login_node
ssh user@login_node
sshare -o "Account,Partition%20"
srun -p <cpu_partition> -N 1 -n 1 -c 4 --mem=16g --time=24:00:00 -A <account> --pty bash -i
srun -p <gpu_partition> -N 1 -n 1 -c 1 --mem=256g --time=24:00:00 -A <account> --pty bash -i
mamba activate Env1
cd analyses/01_rna_genotype_matching


cp /path/to/cybertron_data/project01_vestibular/mouse/vestibular_mouse/scripts/Rscript05_subsetting_sac_and_utr_samples_for_eqtl_mapping.R ./scripts/

# Creating a soft link of data
ln -s /path/to/cybertron_data/project01_vestibular/vestibular_project_data/mouse/plate1/fastq_p1/FTL_* data/batch_01

###---------------- Preparing a list of fastqs, making metedata, and preparing genotype file
Rscript scripts/01_Rscript01_making_fsqlist_metadata_genotype.R



####--------------- Mapping - STAR

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


# Go to WD
cd analyses/01_rna_genotype_matching

# Activate package environment
source "$(conda info --base)/etc/profile.d/conda.sh"
# source $(conda info --base)/etc/profile.d/mamba.sh
conda activate Env1

# Defining input/output
REF="../../resources/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa"
INDEX="../../resources/index_star/"
STAR="STAR"
OUT_DIR="intermediate/bams/"
mkdir -p ${OUT_DIR}

# Getting line from source file - ech line is a library id
FASTQ=$(awk 'NR=='1' && NF' src/fastqs_all.list)
FASTQ=$(awk 'NR=='${SLURM_ARRAY_TASK_ID}' && NF' src/fastqs_all.list) # Ignores empty lines
PASS="../../data/fastq_bulk/"   # batch1-4/


# Making path to triimed reads
R1=${PASS}${FASTQ}_R1_001.fastq.gz
R2=${PASS}${FASTQ}_R2_001.fastq.gz

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
sbatch scripts/02_job01_mapping.sh
589311

# Cheking read groups
BAM="${OUT_DIR}${LIBRARY_ID}_Aligned.sortedByCoord.out.bam"
samtools view -H $BAM | grep '^@RG'

# Cheking if the BAM files are in correct format
gatk ValidateSamFile \
        -I $BAM \
        -R $REF \
        -M SUMMARY

###---------------- Quantify expression - RSEM

# Getting a list of bam files
ls -1 intermediate/bams/*_Aligned.toTranscriptome.out.bam > src/bam_aligned_to_transcriptome.list

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


# Go to WD
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
sbatch scripts/03_job02_count_rsem.sh
589809

### Assemble_expression:
COUNT_DIR="results/count_rsem/"
LIBRARIES="src/lib_ids.list"
ANNOT="../../resources/genome/Mus_musculus.GRCm38.102.chr.collapsed.genes.gtf"
OUT_PREFIX="intermediate/count_table/expr_moues"

# Note 1: Order of inputs is important.
# Note 2: Scrip works with pandas version 1.5.3
python3 scripts/assemble_expression_modified.py ${COUNT_DIR} ${LIBRARIES} ${ANNOT} ${OUT_PREFIX}



###---------------- RUN ASE
# Running bcftools to filter and extract exonic SNPs.
# Getting the exon regions from GTF
GTF_GRCm38="../../resources/genome/Mus_musculus.GRCm38.102.chr.collapsed.genes.gtf"
grep -v '^#' $GTF_GRCm38 | awk '$3=="exon"' | cut -f1,4,5 | gzip -c > src/GRCm38_exon_regions.tsv.gz

bcftools view results/genotype/mouse_geno_all_cln_sorted_maf0.05_geno0.05_RefAltOk_fixed.vcf.gz \
                --regions-file src/GRCm38_exon_regions.tsv.gz \
                -Ou | bcftools view \
                --min-af 0.05:minor \
                -i 'F_MISSING<0.1 && N_ALT=1 && TYPE="snp"' \
                -Oz -o results/genotype/snps_for_mixup_test.vcf.gz

tabix -p vcf results/genotype/snps_for_mixup_test.vcf.gz
# Checking samples
bcftools query -l results/genotype/snps_for_mixup_test.vcf.gz

# Making the batch script

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

# Submitting ASE analysis job
sbatch scripts/07_job4_ASE.sh
647840

## Now, checking the similarity using ASE results
# Making a list of library ids
ls -1 results/ase_counter_result/ | cut -f1 -d"." > src/lib_ids.list

# Define vars
VCF_IN="results/genotype/snps_for_mixup_test.vcf.gz"
ASE_COUNT_DIR="results/ase_counter_result/"
LIB_ID="src/lib_ids.list"
MATRIX="results/ase_counter_result/rna_to_geno_similarity.tsv"
SUMMARY="results/ase_counter_result/rna_to_geno_summary.tsv"

# Runing the similarity check
python3 scripts/qc_rna_to_geno_similarity.py ${VCF_IN} ${ASE_COUNT_DIR} ${LIB_ID} ${MATRIX} ${SUMMARY}

# Making a summary table to look at matched and unmatched strain IDs
Rscript scripts/08_Rscript04_summary_mixup.R



rsync -av user@login_node:/path/to/project/analyses/01_rna_genotype_matching/results/mbv_result .
rsync -av user@login_node:/path/to/project/analyses/01_rna_genotype_matching/figures/02.3_expr_pc1_vs_pc2_DEseq2_normalized_scatter_mislabeled_tagged.jpeg .


# # Extract duplicated strains
# dups <- saccule_w8$Top_Simil_ID_ByMixAna[duplicated(saccule_w8$Top_Simil_ID_ByMixAna)]
# tissue <- "Saccule"
# age <- "w8"

# # Initialize result vectors
# original_id <- c()
# strain <- c()
# lib <- c()

# # Index for storing results
# index <- 1

# # Loop through duplicated strains
# for (i in seq_along(dups)) {
#   # Extract library IDs for the current strain
#   lib_ids <- saccule_w8$library_id[which(saccule_w8$Top_Simil_ID_ByMixAna == dups[i])]
  
#   # Loop through each library ID
#   for (j in seq_along(lib_ids)) {
#     # Check if the library ID meets the condition
#     if (any(meta_mixup$library_id == lib_ids[j] &
#             meta_mixup$tissue == tissue &
#             meta_mixup$age == age &
#             meta_mixup$Original_IDs_In_Geno == dups[i])) {
      
#       # Store data for matching case
#       original_id[index] <- "YES"
#     } else {
#       # Store data for non-matching case
#       original_id[index] <- "NO"
#     }
    
#     # Add strain and library ID to result
#     strain[index] <- dups[i]
#     lib[index] <- lib_ids[j]
    
#     # Increment index
#     index <- index + 1
#   }
# }

# # Combine results into a data frame for easier viewing
# result <- data.frame(
#   Strain = strain,
#   Library_ID = lib,
#   Original_ID = original_id
# )

# print(result)


###---------------- Test QTL mapping
# gzip and index bed files utricle w22
bgzip intermediate/count_table/expr.iqn.filtered.utricle_w22_withdups.bed
bgzip intermediate/count_table/expr.iqn.filtered.utricle_w22_unik.bed
bgzip intermediate/count_table/expr.iqn.filtered.utricle_w22_high_simil_withdups.bed
bgzip intermediate/count_table/expr.iqn.filtered.utricle_w22_high_simil_unik.bed

tabix -p bed intermediate/count_table/expr.iqn.filtered.utricle_w22_withdups.bed.gz
tabix -p bed intermediate/count_table/expr.iqn.filtered.utricle_w22_unik.bed.gz
tabix -p bed intermediate/count_table/expr.iqn.filtered.utricle_w22_high_simil_withdups.bed.gz
tabix -p bed intermediate/count_table/expr.iqn.filtered.utricle_w22_high_simil_unik.bed.gz


# Cleaning genotypes - utricle w22
REF="../../resources/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa"
plink2 --tfile results/genotype/utricle_w22_high_simil_unik \
        --mouse \
        --fa $REF \
        --maf 0.05 \
        --geno 0.05 \
        --mind 0.05 \
        --make-pgen \
        --nonfounders \
        --not-chr X Y \
        --chr 1-19 \
        --snps-only 'just-acgt' \
        --sort-vars \
        --out results/genotype/utricle_w22_high_simil_unik_maf0.05
# 103529 variants remaining after main filters >>> utricle - w22 - with dups
# 103391 variants remaining after main filters >>> utricle - w22 - unik
# 103887 variants remaining after main filters >>> utricle - w22 - with dups - high simil
# 103627 variants remaining after main filters >>> utricle - w22 - unik - high simil

# Correcting REF and ALT alleles - Saccule
plink2 --pfile results/genotype/utricle_w22_high_simil_unik_maf0.05 \
        --fa $REF \
        --maf 0.05 \
        --ref-from-fa force \
        --make-bed \
        --out results/genotype/utricle_w22_high_simil_unik_maf0.05_RefAltOk

# --ref-from-fa force: 53038 variants changed, 50471 validated >>> utricle - w22 - with dups
# --ref-from-fa force: 54728 variants changed, 48643 validated >>> utricle - w22 - unik
# --ref-from-fa force: 54790 variants changed, 49077 validated >>> utricle - w22 - with dups - high simil
# --ref-from-fa force: 56496 variants changed, 47111 validated >>> utricle - w22 - unik - high simil



# For Utricle w22
plink2 --bfile results/genotype/utricle_w22_high_simil_unik_maf0.05_RefAltOk \
        --indep-pairwise 200 100 0.5 \
        --out intermediate/covar/utricle_w22_high_simil_unik_prunlist

plink2 --bfile results/genotype/utricle_w22_high_simil_unik_maf0.05_RefAltOk \
        --extract intermediate/covar/utricle_w22_high_simil_unik_prunlist.prune.in \
        --export vcf bgz id-paste=iid \
        --fa $REF \
        --out intermediate/covar/utricle_w22_high_simil_unik_prunned_for_pca

# Covariate file - Utricle
VCF="intermediate/covar/utricle_w22_high_simil_unik_prunned_for_pca.vcf.gz"
BED="intermediate/count_table/expr.iqn.filtered.utricle_w22_high_simil_unik.bed.gz"
N_GENO_PCS="3" 
N_EXPR_PCS="8"
OUTPUT="intermediate/covar/utricle_w22_high_simil_unik_covars.txt"
Rscript scripts/covariates.R $VCF $BED $N_GENO_PCS $N_EXPR_PCS $OUTPUT




# VCFs
results/genotype/utricle_w22_withdups_maf0.05_RefAltOk
results/genotype/utricle_w22_unik_maf0.05_RefAltOk
results/genotype/utricle_w22_high_simil_withdups_maf0.05_RefAltOk
results/genotype/utricle_w22_high_simil_unik_maf0.05_RefAltOk


# BEDs
intermediate/count_table/expr.iqn.filtered.utricle_w22_withdups.bed.gz
intermediate/count_table/expr.iqn.filtered.utricle_w22_unik.bed.gz
intermediate/count_table/expr.iqn.filtered.utricle_w22_high_simil_withdups.bed.gz
intermediate/count_table/expr.iqn.filtered.utricle_w22_high_simil_unik.bed.gz

# Covars
intermediate/covar/utricle_w22_withdups_covars.txt
intermediate/covar/utricle_w22_unik_covars.txt
intermediate/covar/utricle_w22_high_simil_withdups_covars.txt
intermediate/covar/utricle_w22_high_simil_unik_covars.txt



## Map cis-eQTLs, determining significance using permutations, outputs the top association per gene.
# mkdir -p eQTL/cis_eQTL/
GENO="results/genotype/utricle_w22_withdups_maf0.05_RefAltOk"
BED="intermediate/count_table/expr.iqn.filtered.utricle_w22_withdups.bed.gz"
COVAR="intermediate/covar/utricle_w22_withdups_covars.txt"
OUTPUT="results/qtl_mapping/cis_eqtl_utricle_w22_withdups.txt.gz"

python3 scripts/run_tensorqtl.py \
        ${GENO} \
        ${BED} \
        ${OUTPUT} \
        --covariates ${COVAR} \
        --mode cis

sbatch scripts/10_job5_eqtl_mapping_test.sh



### Downloading mouse single cell data
HOST="ftp://igm-storage.ucsd.edu/250203_LH00444_0284_A22TGLTLT3/"
USER="telese"
PASS="JsdYkQBKiW"

# List files using curl
curl -u $USER:$PASS $HOST --list-only

# Or
# Print the directory listing to STDOUT
wget --ftp-user="$USER" \
        --ftp-password="$PASS" \
        $HOST -O -

HOST="ftp://igm-storage.ucsd.edu/250203_LH00444_0284_A22TGLTLT3"
USER="telese"
PASS="JsdYkQBKiW"
DIR="/path/to/cybertron_data/project01_vestibular/scRNA_data_Feb10_25/"
wget -P $DIR --ftp-user="$USER" --ftp-password="$PASS" $HOST/*

