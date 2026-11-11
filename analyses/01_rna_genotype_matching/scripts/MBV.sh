#!/bin/bash
### Running MBV to ensure good match between the sequence and genotype data

# Path to QTLtools
QTLtools="QTLtools"

# Preparing genotype file, first, sorting and basic filtering
REF="resources/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa"
plink2  --tfile results/genotype/mouse_geno_all \
        --mouse \
        --fa $REF \
        --maf 0.05 \
        --geno 0.05 \
        --nonfounders \
        --not-chr X Y \
        --chr 1-19 \
        --snps-only 'just-acgt' \
        --sort-vars \
        --make-pgen \
        --out results/genotype/mouse_geno_all_cln_sorted_maf0.05_geno0.05

# Second, correcting REF and ALT alleles
plink2  --pfile results/genotype/mouse_geno_all_cln_sorted_maf0.05_geno0.05 \
        --fa $REF \
        --maf 0.05 \
        --ref-from-fa force \
        --recode vcf id-paste=iid \
        --out results/genotype/mouse_geno_all_cln_sorted_maf0.05_geno0.05_RefAltOk

# --export vcf bgz id-paste=iid \
bgzip results/genotype/mouse_geno_all_cln_sorted_maf0.05_geno0.05_RefAltOk.vcf
tabix -p vcf results/genotype/mouse_geno_all_cln_sorted_maf0.05_geno0.05_RefAltOk.vcf.gz

# We have an non-standard line in heder (##chrSet=<ID=1,autosomePairCt=19,X,Y>), removing it,
zcat $VCF | grep -v "##chrSet=" | bgzip > results/genotype/mouse_geno_all_cln_sorted_maf0.05_geno0.05_RefAltOk_fixed.vcf.gz
tabix -p vcf results/genotype/mouse_geno_all_cln_sorted_maf0.05_geno0.05_RefAltOk_fixed.vcf.gz

### Run the MBV
BAM="intermediate/bams/FTL_602_Aligned.sortedByCoord.out.bam"
VCF="results/genotype/mouse_geno_all_cln_sorted_maf0.05_geno0.05_RefAltOk.vcf.gz"

# Index bam
samtools index $BAM

# Run the mbv
$QTLtools       mbv \
                --bam $BAM \
                --vcf $VCF \
                --filter-mapping-quality 150 \
                --filter-base-quality 20 \
                --filter-minimal-coverage 10 \
                --out results/mbv_result/Mm_FTL_602_BXD45_RwwJ.bamstat.txt

### FTL_872 (BXD77 in given metadata but matched by BXD75 in our ASE approach analysis)
BAM_MM_872="intermediate/bams/FTL_872_Aligned.sortedByCoord.out.bam"
VCF_MM="results/genotype/mouse_geno_all_cln_sorted_maf0.05_geno0.05_RefAltOk.vcf.gz"

# Index bam
samtools index $BAM_MM_872

# Run the mbv
$QTLtools   mbv \
                --bam $BAM_MM_872 \
                --vcf $VCF_MM \
                --filter-mapping-quality 150 \
                --filter-base-quality 20 \
                --filter-minimal-coverage 10 \
                --out results/mbv_result/Mm_FTL_872_BXD77_RwwJ.bamstat.txt



###----------- Test on human data
BAM_HS="/path/to/cybertron_data/project01_vestibular/vestibular_project_data/human/RNAseq/ProcessedData/BamFiles_RGok/FTL_817_Aligned.sortedByCoord.out.bam"
VCF_HS="/path/to/cybertron_data/project01_vestibular/human/vestibular_project_human/mixup_analysis/geno_phased_imputed_clean_for_ASE_noDup.vcf.gz"

# Index bam
samtools index $BAM_HS

# Run the mbv
$QTLtools   mbv \
                --bam $BAM_HS \
                --vcf $VCF_HS \
                --filter-mapping-quality 150 \
                --filter-base-quality 20 \
                --filter-minimal-coverage 10 \
                --out results/mbv_result/Hs_FTL_817_FTL_32287621.bamstat.txt

### FTL_844 (FTL_23884570 matched by FTL_31819005 in our ASE approach analysis)
BAM_HS_844="/path/to/cybertron_data/project01_vestibular/vestibular_project_data/human/RNAseq/ProcessedData/BamFiles_RGok/FTL_844_Aligned.sortedByCoord.out.bam"
samtools index $BAM_HS_844 # Index bam
$QTLtools   mbv \
                --bam $BAM_HS_844 \
                --vcf $VCF_HS \
                --filter-mapping-quality 150 \
                --filter-base-quality 20 \
                --filter-minimal-coverage 10 \
                --out results/mbv_result/Hs_FTL_844_FTL_23884570.bamstat.txt



### FTL_834 (FTL_17465303 matched by FTL_10380376 in our ASE approach analysis)
BAM_HS_834="/path/to/cybertron_data/project01_vestibular/vestibular_project_data/human/RNAseq/ProcessedData/BamFiles_RGok/FTL_834_Aligned.sortedByCoord.out.bam"
samtools index $BAM_HS_834 # Index bam
$QTLtools   mbv \
                --bam $BAM_HS_834 \
                --vcf $VCF_HS \
                --filter-mapping-quality 150 \
                --filter-base-quality 20 \
                --filter-minimal-coverage 10 \
                --out results/mbv_result/Hs_FTL_834_FTL_17465303.bamstat.txt


### FTL_551 (FTL_31907773 matched by FTL_31854925 in our ASE approach analysis)
BAM_HS_551="/path/to/cybertron_data/project01_vestibular/vestibular_project_data/human/RNAseq/ProcessedData/BamFiles_RGok/FTL_551_Aligned.sortedByCoord.out.bam"
samtools index $BAM_HS_551 # Index bam
$QTLtools    mbv \
                --bam $BAM_HS_551 \
                --vcf $VCF_HS \
                --filter-mapping-quality 150 \
                --filter-base-quality 20 \
                --filter-minimal-coverage 10 \
                --out results/mbv_result/Hs_FTL_551_FTL_31907773.bamstat.txt


PATH_TO_DEST="/path/to/project/analyses/01_rna_genotype_matching/results/mbv_result"
rsync -avz user@login_node:${PATH_TO_DEST}/* .




ls -1 intermediate/bams/FTL_*_Aligned.sortedByCoord.out.bam | grep -v -E "FTL_752|FTL_775" > src/bam_aligned_to_genome.list
sbatch scripts/05_job3_mbv.sh
629879
