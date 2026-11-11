#!/bin/bash

srun -p <cpu_partition> -N 1 -n 1 -c 4 --mem=16g --time=8:00:00 -A <account> --pty bash -i
mamba activate Env1

# Making dirs
mkdir -p results/genotype data/bam data/rnaseq intermediate config env logs results rules scripts src 
mkdir -p results/genotype/utricle_w8 results/genotype/utricle_w22 results/genotype/saccule_w8 results/genotype/saccule_w22
mkdir -p intermediate/count_rsem intermediate/count_matrix

REF="../../resources/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa"

# Copy GWAS results
cp -f ../../data/gwas_results/*.gwas results/gwas_results/

# geno
plink2  --tfile results/genotype/mouse_geno_all \
        --mouse \
        --make-pgen \
        --not-chr Y \
        --sort-vars \
        --nonfounders \
        --out results/genotype/geno_all

plink2  --pfile results/genotype/geno_all \
        --fa $REF \
        --ref-from-fa force \
        --mouse \
        --make-pgen \
        --not-chr Y \
        --sort-vars \
        --nonfounders \
        --out results/genotype/geno_all_refalt_ok


# Making ld ref file
for chr in {1..20}; do
plink2  --pfile results/genotype/geno_all_refalt_ok \
        --make-bed \
        --nonfounders \
        --chr $chr \
        --out results/genotype/LDREF/mouse_mm10_${chr}
done

# Copy to Pantry dir
cp -r results/genotype/LDREF  ../../analyses/03_pantry/pheast/input/

### Convert to sumstat format
Rscript scripts/clean_gwas_to_make_sumstat.R

# Getting LDSC tool and creating the env
git clone https://github.com/bulik/ldsc.git
cd ldsc
mamba env create --file environment.yml
source activate ldsc
mamba activate ldsc

./ldsc.py -h
./munge_sumstats.py -h

cd .. # (analyses/02_gwas_sumstats)
for file in results/gwas_results/*.gwas.cln; do
        base_name=$(basename "$file" | sed -E 's/(\.[0-9]{8})?\.gwas\.cln$//')
        src/ldsc/munge_sumstats.py --sumstats ${file} \
                                --out results/gwas_results/sumstats_format/${base_name} \
                                --N 100
done

for file in results/gwas_results/sumstats_format/*.gz; do
        gunzip $file
done

# Copy to Pantry dir
cp results/gwas_results/sumstats_format/*.sumstats ../../analyses/03_pantry/pheast/input/gwas_sumstat/

