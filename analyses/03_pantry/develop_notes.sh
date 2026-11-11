#!/bin/bash
srun -p <cpu_partition> -N 1 -n 1 -c 16 --mem=64g --time=8:00:00 -A <account> --pty bash -i
srun -p cpu-core -N 1 -n 1 -c 16 --mem=64g --time=8:00:00 -A core --pty bash -i



### Some useful codes
# If cancelled a job and need to unlock snakemake:
snakemake --configfile config.yml --config dataset=utricle_w22 --unlock

snakemake   --configfile config.yml \
            --config dataset=saccule_w22 \
            --jobs 34 \
            --default-resources mem_mb=185000 cpus=34 \
            --resources gpu=1 \
            --printshellcmds \
            --latency-wait 60 

# -R tensorqtl_trans > sac_w8_qtl_trans.log 2>&1
# --keep-going
# --rerun-incomplete
# #SBATCH --job-name=pantry_pheast_utr_w22

scontrol show job 2058155 | grep JobName
scontrol show job 2058154 | grep JobName


### Run phenotyping
# Making the pantry env
cd phenotyping/
mamba env create -n pantry --file environment.yml
mamba activate pantry

# Downloading the input and dry run, still in "phenotyping/"
wget -O test_input_phenotyping.tar.gz "https://www.dropbox.com/scl/fi/5c4yjrqji9catsazc8nxd/test_input_phenotyping.tar.gz?rlkey=jsppnzu4cdls9fp5z3xxm1i37&dl=0"
tar -xzvf test_input_phenotyping.tar.gz
rm test_input_phenotyping.tar.gz
snakemake -j 24 -n --printshellcmds


# Getting read length
FQ="input/utricle_w8/fastq/FTL_855_S5_L003_R1_001.fastq.gz"
zcat $FQ | awk 'NR==2 {print length($0); exit}'


for fq in input/utricle_w8/fastq/*_R1_*.fastq.gz; do
    echo -n "$fq: "
    zcat "$fq" | awk 'NR==2 {print length($0); exit}'
done


### Run Phenotyping on real data
# Utricle_w22
snakemake --configfile config.yml --config dataset=utricle_w22 --cores 24 -n
# Submit as a job
sbatch run_phenotyping_utricle_w22.slurm
1941978
1966090

# Utricle_w8
snakemake --configfile config.yml --config dataset=utricle_w8 --cores 24 -n
# Submit as a job
sbatch run_phenotyping_utricle_w8.slurm
1946813

# Saccule_w8
snakemake --configfile config.yml --config dataset=saccule_w8 --cores 24 -n
# Submit as a job
sbatch run_phenotyping_saccule_w8.slurm
1960731

# Saccule_w22
snakemake --configfile config.yml --config dataset=saccule_w22 --cores 8 -n
# Submit as a job
sbatch run_phenotyping_saccule_w22.slurm
1965876




###----------- Run pheast
cd ../pheast
mamba env create -n pantry_pheast --file environment.yml
mamba activate pantry_pheast

wget -O test_input_pheast.tar.gz "https://www.dropbox.com/scl/fi/5k2ybqyrtmjll2ig6nnxv/test_input_pheast.tar.gz?rlkey=yzybnhb60b3yxoebrhsoaoumz&dl=0"
tar -xzvf test_input_pheast.tar.gz
rm test_input_pheast.tar.gz
snakemake -j 24 -n --printshellcmds

r-stringi=1.7.6

# Run Pheast on real data
### Run pheast
cd ../pheast
mamba activate pantry_pheast

### Making cross-modality file for each tissue_age
bash scripts/combine_modalities_updated_060425.sh --path_to_output output/utricle_w8 alt_polyA alt_TSS expression isoforms splicing stability
bash scripts/combine_modalities_updated_060425.sh --path_to_output output/utricle_w22 alt_polyA alt_TSS expression isoforms splicing stability
bash scripts/combine_modalities_updated_060425.sh --path_to_output output/saccule_w8 alt_polyA alt_TSS expression isoforms splicing stability
bash scripts/combine_modalities_updated_060425.sh --path_to_output output/saccule_w22 alt_polyA alt_TSS expression isoforms splicing stability

### Copy phenotypes from phenotype dir to pheast (with checking, sorting and filtering, if needed)
# saccule w8
awk '{ print $2 "\t" $2 }' input/saccule_w8/saccule_w8.fam > input/saccule_w8/saccule_w8.sample_map.txt
awk '{ print $2 }' input/saccule_w8/saccule_w8.fam > input/saccule_w8/saccule_w8.individuals.txt

./scripts/prepare_phenotypes.py \
    --indir ../phenotyping/output/saccule_w8 \
    --outdir input/saccule_w8/phenotypes/ \
    --map input/saccule_w8/saccule_w8.sample_map.txt \
    --individuals input/saccule_w8/saccule_w8.individuals.txt
cp ../phenotyping/input/saccule_w8/saccule_w8.samples.txt input/saccule_w8/


# saccule w22
awk '{ print $2 "\t" $2 }' input/saccule_w22/saccule_w22.fam > input/saccule_w22/saccule_w22.sample_map.txt
awk '{ print $2 }' input/saccule_w22/saccule_w22.fam > input/saccule_w22/saccule_w22.individuals.txt

./scripts/prepare_phenotypes.py \
    --indir ../phenotyping/output/saccule_w22 \
    --outdir input/saccule_w22/phenotypes/ \
    --map input/saccule_w22/saccule_w22.sample_map.txt \
    --individuals input/saccule_w22/saccule_w22.individuals.txt
cp ../phenotyping/input/saccule_w22/saccule_w22.samples.txt input/saccule_w22/


# utricle w8
awk '{ print $2 "\t" $2 }' input/utricle_w8/utricle_w8.fam > input/utricle_w8/utricle_w8.sample_map.txt
awk '{ print $2 }' input/utricle_w8/utricle_w8.fam > input/utricle_w8/utricle_w8.individuals.txt

./scripts/prepare_phenotypes.py \
    --indir ../phenotyping/output/utricle_w8 \
    --outdir input/utricle_w8/phenotypes/ \
    --map input/utricle_w8/utricle_w8.sample_map.txt \
    --individuals input/utricle_w8/utricle_w8.individuals.txt
cp ../phenotyping/input/utricle_w8/utricle_w8.samples.txt input/utricle_w8/


# utricle w22
awk '{ print $2 "\t" $2 }' input/utricle_w22.fam > input/utricle_w22.sample_map.txt
awk '{ print $2 }' input/utricle_w22.fam > input/utricle_w22.individuals.txt

./scripts/prepare_phenotypes.py \
    --indir ../phenotyping/output/utricle_w22 \
    --outdir input/utricle_w22/phenotypes/ \
    --map input/utricle_w22/utricle_w22.sample_map.txt \
    --individuals input/utricle_w22/utricle_w22.individuals.txt
cp ../phenotyping/input/utricle_w22/utricle_w8.samples.txt input/utricle_w22/


### GPU setup
srun --account=gpu-core-sponsored --partition=gpu-core-sponsored --nodes=1 --ntasks=1 --cpus-per-task=1 --gres=gpu:1 --mem-per-cpu=64G --time=0-00:10:00 --pty /bin/bash
mamba activate pantry_pheast # Just to use apptainer

# Activating docker container
CONTAINER="/path/to/posit-base-20250410-ppm-testing.sif"
apptainer shell --nv --bind /path/to/data $CONTAINER

# Activating conda env in the docker container
MAMBA="$(conda info --base)/etc/profile.d/mamba.sh"
source $MAMBA
mamba activate pantry_pheast

# Check - in container
whoami
nvidia-smi
nvcc --version
rpm -qa | grep cuda

# Check - in python
python
import torch
print(torch.cuda.is_available()) # Should print true
print(torch.cuda.get_device_name(0))



### Submit slurm job using gpu
# First dry run
cd pheast/
snakemake --configfile config.yml --config dataset=utricle_w22 --jobs 16 --printshellcmds -n
snakemake --configfile config.yml --config dataset=utricle_w8 --jobs 16 --printshellcmds -n 


### Main run
# utricle_w8
sbatch --export=ALL,DATASET=utricle_w8 --job-name=pheast_utricle_w8 run_pheast_gpu_wrapper.slurm
2058813
# On cpus
sbatch --export=ALL,DATASET=utricle_w8 --job-name=pheast_utricle_w8 run_pheast_cpu.slurm
2058867


# utricle_w22
sbatch --export=ALL,DATASET=utricle_w22 --job-name=pheast_utricle_w22 run_pheast_gpu_wrapper.slurm
2058817
# On cpus
sbatch --export=ALL,DATASET=utricle_w22 --job-name=pheast_utricle_w22 run_pheast_cpu.slurm
2058819


# saccule_w8
sbatch --export=ALL,DATASET=saccule_w8 --job-name=pheast_saccule_w8 run_pheast_gpu_wrapper.slurm
2058141
# Run on cpu
sbatch --export=ALL,DATASET=saccule_w8 --job-name=pheast_saccule_w8 run_pheast_cpu.slurm
2058183

# saccule_w22
sbatch --export=ALL,DATASET=saccule_w22 --job-name=pheast_saccule_w22 run_pheast_gpu_wrapper.slurm
2057616
# Run on cpu
sbatch --export=ALL,DATASET=saccule_w22 --job-name=pheast_saccule_w22 run_pheast_cpu.slurm
2058152


### Run TWAS Association Test
cd pheast/twas_assoc/
# need to install  bioconda::r-coloc, but it downgrads lots of packages
snakemake -s run_twas_with_tissue.smk --jobs 16 --printshellcmds -n
# Submit as a job
sbatch ./run_twas_assoc.slurm
2060067




###---------------- Dl mouse RNA edit sites ref file
wget -P Pantry_RNA_editing/phenotyping/input/ref/ http://srv00.recas.ba.infn.it/webshare/ATLAS/download/TABLE1_mm10_v3.txt.gz
Rscript prep.R

cd Pantry_RNA_editing/phenotyping/
mamba activate pantry
mamba create -n apptainer_env -c conda-forge apptainer

### utricle_w22
# Dry-run
snakemake --configfile config.yml --config dataset=utricle_w22 --cores 16 -n
# Submit the job - on cpu
sbatch --export=ALL,DATASET=utricle_w22 --job-name=pheno_utricle_w22 run_phenotyping.slurm
2060266

### utricle_w8
# Dry-run
snakemake --configfile config.yml --config dataset=utricle_w8 --cores 16 -n
# Submit the job - on gpu
sbatch --export=ALL,DATASET=utricle_w8 --job-name=pheno_utricle_w8 run_phenotyping.slurm
2060314

### saccule_w8
# Dry-run
snakemake --configfile config.yml --config dataset=saccule_w8 --cores 16 -n
# Submit the job - on gpu
sbatch --export=ALL,DATASET=saccule_w8 --job-name=pheno_saccule_w8 run_phenotyping.slurm
2060315

### saccule_w22
# Dry-run
snakemake --configfile config.yml --config dataset=saccule_w22 --cores 16 -n
# Submit the job - on gpu
sbatch --export=ALL,DATASET=saccule_w22 --job-name=pheno_saccule_w22 run_phenotyping.slurm
2060316



###-------------------------- Pheast
mamba activate pantry_pheast

## utricle w22
# Making cross-modality
cd phenotyping
bash scripts/combine_modalities_updated_060425.sh --path_to_output output/utricle_w22 alt_polyA alt_TSS expression RNA_editing isoforms splicing stability

# Copy phenotypes to pheast input/
cd pheast
cp ../../pheast/input/utricle_w22/utricle_w22.* input/utricle_w22/

./scripts/prepare_phenotypes.py \
    --indir ../phenotyping/output/utricle_w22 \
    --outdir input/utricle_w22/phenotypes/ \
    --map input/utricle_w22/utricle_w22.sample_map.txt \
    --individuals input/utricle_w22/utricle_w22.individuals.txt

# Dry run
snakemake --configfile config.yml --config dataset=utricle_w22 --cores 16 --printshellcmds -n
# Submit on gpu
sbatch --export=ALL,DATASET=utricle_w22 --job-name=pheast_utricle_w22 run_pheast_gpu_wrapper.slurm
2361072

# Run on cpu to complete the trans qtl mapping
sbatch --export=ALL,DATASET=utricle_w22 --job-name=pheast_utricle_w22 run_pheast_cpu.slurm




## utricle w8
# Making cross-modality
cd phenotyping
bash scripts/combine_modalities_updated_060425.sh --path_to_output output/utricle_w8 alt_polyA alt_TSS expression RNA_editing isoforms splicing stability

# Copy phenotypes to pheast input/
cd pheast
cp ../../pheast/input/utricle_w8/utricle_w8.* input/utricle_w8/

./scripts/prepare_phenotypes.py \
    --indir ../phenotyping/output/utricle_w8 \
    --outdir input/utricle_w8/phenotypes/ \
    --map input/utricle_w8/utricle_w8.sample_map.txt \
    --individuals input/utricle_w8/utricle_w8.individuals.txt

# Dry run
snakemake --configfile config.yml --config dataset=utricle_w8 --cores 16 --printshellcmds -n
# Submit on gpu
sbatch --export=ALL,DATASET=utricle_w8 --job-name=pheast_utricle_w8 run_pheast_gpu_wrapper.slurm
2361110

# Run on cpu to complete the trans qtl mapping
sbatch --export=ALL,DATASET=utricle_w8 --job-name=pheast_utricle_w8 run_pheast_cpu.slurm



## saccule_w22
# Making cross-modality
cd phenotyping
bash scripts/combine_modalities_updated_060425.sh --path_to_output output/saccule_w22 alt_polyA alt_TSS expression RNA_editing isoforms splicing stability

# Copy phenotypes to pheast input/
cd pheast
cp ../../pheast/input/saccule_w22/saccule_w22.* input/saccule_w22/
./scripts/prepare_phenotypes.py \
    --indir ../phenotyping/output/saccule_w22 \
    --outdir input/saccule_w22/phenotypes/ \
    --map input/saccule_w22/saccule_w22.sample_map.txt \
    --individuals input/saccule_w22/saccule_w22.individuals.txt

# Dry run
snakemake --configfile config.yml --config dataset=saccule_w22 --cores 16 --printshellcmds -n
# Submit on gpu
sbatch --export=ALL,DATASET=saccule_w22 --job-name=pheast_saccule_w22 run_pheast_gpu_wrapper.slurm
2361121

# Run on cpu to complete the trans qtl mapping
sbatch --export=ALL,DATASET=saccule_w22 --job-name=pheast_saccule_w22 run_pheast_cpu.slurm



## saccule_w8
# Making cross-modality
cd phenotyping
bash scripts/combine_modalities_updated_060425.sh --path_to_output output/saccule_w8 alt_polyA alt_TSS expression RNA_editing isoforms splicing stability

# Copy phenotypes to pheast input/
cd pheast
cp ../../pheast/input/saccule_w8/saccule_w8.* input/saccule_w8/
./scripts/prepare_phenotypes.py \
    --indir ../phenotyping/output/saccule_w8 \
    --outdir input/saccule_w8/phenotypes/ \
    --map input/saccule_w8/saccule_w8.sample_map.txt \
    --individuals input/saccule_w8/saccule_w8.individuals.txt

# Dry run
snakemake --configfile config.yml --config dataset=saccule_w8 --cores 16 --printshellcmds -n
# Submit on gpu
sbatch --export=ALL,DATASET=saccule_w8 --job-name=pheast_saccule_w8 run_pheast_gpu_wrapper.slurm
2361131

# Run on cpu to complete the trans qtl mapping
sbatch --export=ALL,DATASET=saccule_w8 --job-name=pheast_saccule_w8 run_pheast_cpu.slurm



### TWAS Association
cd pheast/twas_assoc/
# need to install  bioconda::r-coloc, but it downgrads lots of packages
snakemake -s run_twas_with_tissue.smk --jobs 16 --printshellcmds -n
# Submit as a job
sbatch ./run_twas_assoc.slurm
2378171

PATH_TO_FILE="/path/to/project/analyses/01_rna_genotype_matching/figures/02.3_expr_pc1_vs_pc2_DEseq2_normalized_scatter_mislabeled_tagged.jpeg"
scp user@login_node:$PATH_TO_FILE .
