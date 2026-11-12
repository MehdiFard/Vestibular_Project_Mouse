#!/bin/bash

cd analyses/07_single_cell_rna
# srun -p <cpu_partition> -N 1 -n 1 -c 16 --mem=64g --time=24:00:00 -A <account> --pty bash -i
srun -p <cpu_partition> -N 1 -n 1 -c 32 --mem=128g --time=8:00:00 -A <account> --pty bash -i
srun -p cpu-core -N 1 -n 1 -c 16 --mem=64g --time=8:00:00 -A core --pty bash -i

mamba activate Env
mamba activate seurat_env

# Making soft links to single-cell data
ln -s /path/to/cybertron_data/project01_vestibular/new_data_Nov12_24/241029_LH00444_0227_B22NW5CLT3/*_csRNAseq_rat_* data/fastqs_raw/

####--------------- Quality control - fastqc
# Run FastQC
OUT="intermediate/fastqc_report_raw_fqs/"
fastqc -t 8 data/fastqs_raw/* -o $OUT

# Running MultiQC
FQ_REPORT_DIR="intermediate/fastqc_report_raw_fqs/"
multiqc $FQ_REPORT_DIR -o $FQ_REPORT_DIR

# Copying to local system (Should be ran on local system)
PATH_TO_DEST="/path/to/project/analyses/07_single_cell_rna/intermediate/fastqc_report_raw_fqs"
rsync -avz user@login_node:${PATH_TO_DEST}/*.html .


### Runing qc by Snakemake
snakemake -j 16 --printshellcmds --scheduler greedy -n
sbatch run_workflow.slurm
2830104


snakemake -j 8 --printshellcmds --scheduler greedy --default-resources mem_mb=128000 cpus=16 -R cellranger_count -n



# install Cell Ranger
cd src/tools
tar -xzvf cellranger-9.0.1.tar.gz

# Unpack the reference data for mouse
cd ..
tar -xzvf refdata-gex-GRCm39-2024-A.tar.gz

export PATH=$PWD/src/tools/cellranger-9.0.1:$PATH

# or better, putting it in bashrc
echo 'export PATH=/path/to/project/analyses/07_single_cell_rna/src/tools/cellranger-9.0.1:$PATH' >> ~/.bashrc
source ~/.bashrc


mkdir -p intermediate/cell_ranger_output/
cellranger count \
    --id FTL_1146 \
    --output-dir intermediate/cell_ranger_output/FTL_1146 \
    --create-bam true \
    --transcriptome src/refdata-gex-GRCm39-2024-A \
    --fastqs data/fastq_scRNA \
    --sample FTL_1146_snRNAseq \
    --include-introns=true \
    --localcores 16 \
    --localmem 64

# --chemistry=SC3Pv4
# --id FTL_1146 \ # If using "-output-dir", the "--id" will be ignored, but mus be provided.

# Dry run
snakemake -j 8 --printshellcmds --scheduler greedy --default-resources mem_mb=64000 cpus=16  -n

# Submit
sbatch run_workflow.slurm
2840162

cellranger_count


###---------------- Filtering reads with overrepresented sequences
# Getting read names with overrepresented sequences from MultiQC output
MULTIQC_OVERREP="intermediate/fastqc_report_raw_fqs/multiqc_data/fastqc_overrepresented_sequences_plot.txt"
awk -F'\t' 'NR > 1 && ($2 != 0.0 || $3 != 0.0)' ${MULTIQC_OVERREP} | cut -f1 > src/ReadsWithDups.txt

touch src/overrepresented_seq.txt
PATH_to_FASTQC="intermediate/fastqc_report_raw_fqs/"

while IFS= read -r READ; do
    unzip -p ${PATH_to_FASTQC}${READ}_fastqc.zip ${READ}_fastqc/fastqc_data.txt > src/fastqc_data.txt
    awk '/>>Overrepresented sequences/,/>>END_MODULE/' src/fastqc_data.txt | \
    grep -v '>>Overrepresented sequences\|>>END_MODULE' | \
    tail -n +2 | cut -f1,4 >> src/overrepresented_seq.txt
done < "src/ReadsWithDups.txt"
rm src/fastqc_data.txt

# Run Rscript to make a fasta file contains overrepresented seqs to be used in trimmomatic
Rscript scripts/01_Rscript01_MakingFastaWithOverrepresentedSeqs.R "src/overrepresented_seq.txt"



###---------------- Trimming - Trimmomatic
# making a list of unique IDs (with path) to the fastqs for trimming
ls -1 data/fastqs_raw/ | cut -f1-5 -d"_" | sort -u > src/library_ids_unik.txt


#!/bin/bash

#SBATCH --job-name=trimming
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
#SBATCH --array=1-8

# WD
cd analyses/07_single_cell_rna

# Activate package environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate Env1

# Trimming the reads using trimmomatic
PATH_FQ="data/fastqs_raw/"
OUT="data/fastqs_trimmed/"

# Getting run ID for each job
LIB_ID=$(awk 'NR=='${SLURM_ARRAY_TASK_ID}' && NF' src/library_ids_unik.txt) # Ignores empty lines

R1="${LIB_ID}_L005_R1_001.fastq.gz"
R2="${LIB_ID}_L005_R2_001.fastq.gz"

trimmomatic PE -threads 14 -phred33 ${PATH_FQ}${R1} ${PATH_FQ}${R2} \
    ${OUT}${LIB_ID}_R1_trimmed.fq.gz ${OUT}${LIB_ID}_R1_single.fq.gz \
    ${OUT}${LIB_ID}_R2_trimmed.fq.gz ${OUT}${LIB_ID}_R2_single.fq.gz \
    ILLUMINACLIP:${SLURM_SUBMIT_DIR}/src/adaptor_overrepseq.fasta:2:30:10 MAXINFO:90:0.8  MINLEN:85

# Submitting this job
sbatch scripts/02_job01_trimming.sh
639174



####--------------- Quality control, round 2, on trimmed fastqs - fastqc
# Run FastQC
OUT="intermediate/fastqc_report_trimmed_fqs/"
fastqc -t 8 data/fastqs_trimmed/*_trimmed.fq.gz -o $OUT

# Running MultiQC
FQ_REPORT_DIR="intermediate/fastqc_report_trimmed_fqs/"
multiqc $FQ_REPORT_DIR -o $FQ_REPORT_DIR



####--------------- QC and trim using fastp
PATH_TO_FQS="data/fastqs_raw"
OUT_DIR="data/fastqs_fastp_out"

while IFS= read -r LIB_ID_BASE; do

    R1="${PATH_TO_FQS}/${LIB_ID_BASE}_L005_R1_001.fastq.gz"
    R2="${PATH_TO_FQS}/${LIB_ID_BASE}_L005_R2_001.fastq.gz"

    fastp -i "$R1" -I "$R2" \
    -o "${OUT_DIR}/${LIB_ID_BASE}.fastp_R1.fq.gz" \
    -O "${OUT_DIR}/${LIB_ID_BASE}.fastp_R2.fq.gz" \
    -h ${OUT_DIR}/${LIB_ID_BASE}.fastp.html \
    -j ${OUT_DIR}/${LIB_ID_BASE}.fastp.json

done < src/library_ids_unik.txt


####--------------- Quality control, round 3, on trimmed fastqs by fastp
# Run FastQC
OUT="intermediate/fastqc_report_trimmed_fqs_by_fastp"
fastqc -t 8 data/fastqs_fastp_out/*.fastp_R[12].fq.gz -o $OUT


# Running MultiQC
FQ_REPORT_DIR="intermediate/fastqc_report_trimmed_fqs_by_fastp/"
multiqc $FQ_REPORT_DIR -o $FQ_REPORT_DIR

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
cd /path/to/vestibulart_project_mouse

# Activate package environment
source "$(conda info --base)/etc/profile.d/conda.sh"
# source $(conda info --base)/etc/profile.d/mamba.sh
conda activate Env1

# Defining input/output
REF="../../resources/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa"
INDEX="../../resources/index_star/"
STAR="STAR"
OUT_DIR="data/bams/"
mkdir -p ${OUT_DIR}

# Getting line from source file - ech line is a library id
FASTQ=$(awk 'NR=='1' && NF' src/fastqs_all.list)
FASTQ=$(awk 'NR=='${SLURM_ARRAY_TASK_ID}' && NF' src/fastqs_all.list) # Ignores empty lines
PATH_to_FASTQC="data/rna_fastqs/all_batches_pooled/"


# Making path to triimed reads
R1=${PATH_to_FASTQC}${FASTQ}_R1_001.fastq.gz
R2=${PATH_to_FASTQC}${FASTQ}_R2_001.fastq.gz

# getting required info for RG
LIBRARY_ID=$(echo $FASTQ | cut -f1-2 -d"_")
SAMPLE_ID=$(awk -v search="$LIBRARY_ID" -F'\t' '$1 == search { print $2 }' src/lib_sam_ins_4RG.tsv)
INSERT_SIZE=$(awk -v search="$LIBRARY_ID" -F'\t' '$1 == search { print $3 }' src/lib_sam_ins_4RG.tsv)

# getting identifier reads from fastq file and extracting the info for RG
IDENTIFIRE=$(head -n 1 <(zcat ${R1}))
INSTRO_NAME=$(echo $IDENTIFIRE | awk -F ':' '{print $1}' | cut -f 2 -d " " | sed 's/@//')
LIB_ID=$(echo $IDENTIFIRE | awk -F ':' '{print $2}')
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
        --twoPATH_to_FASTQCMode Basic \
        --quantMode TranscriptomeSAM GeneCounts \
        --outSAMstrandField intronMotif \
        --outSAMattrRGline ID:$RGID SM:$RGSM LB:$RGLB PL:$RGPL PU:$RGPU \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outFileNamePrefix ${OUT_DIR}${LIBRARY_ID}_
# --outSAMmapqUnique 60 \ # If we need to change the MAPQ max defulte from 255 to 60

### END
sbatch scripts/02_job01_mapping.sh


mamba create -n seurat_env -c conda-forge \
    r-base r-seurat r-seuratobject r-rspectra r-sp r-tidyverse r-hdf5r


### Transfer data to local machine
# Generate a key (on local machine)
ssh-keygen -t ed25519 -C "<email>"
ssh-copy-id user@login_node
# test
ssh user@login_node
exit

# Copy Result dir to local
rsync -avz user@login_node:/path/to/vestibulart_project_mouse_single_cell/result/ ./figures
# Copy scripts to local
rsync -avz user@login_node:/path/to/vestibulart_project_mouse_single_cell/scripts .

# Copy cellranger output to local
for sample in FTL_1146 FTL_1147 FTL_1148 FTL_1149; do
    echo "Copying $sample"
    mkdir -p $sample
    rsync -avz user@login_node:/path/to/project/analyses/07_single_cell_rna/intermediate/cell_ranger_output/${sample}/outs/filtered_feature_bc_matrix ./${sample}
    rsync -avz user@login_node:/path/to/project/analyses/07_single_cell_rna/intermediate/cell_ranger_output/${sample}/outs/raw_feature_bc_matrix ./${sample}
    rsync -avz user@login_node:/path/to/project/analyses/07_single_cell_rna/intermediate/cell_ranger_output/${sample}/outs/filtered_feature_bc_matrix.h5 ./${sample}
    rsync -avz user@login_node:/path/to/project/analyses/07_single_cell_rna/intermediate/cell_ranger_output/${sample}/outs/web_summary.html./${sample}
done




