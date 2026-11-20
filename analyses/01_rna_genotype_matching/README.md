# 01 Sample matching

Checks that every bulk RNA-seq library comes from the BXD strain it is labelled with, then writes the per-group
inputs of Pantry (03). Two independent approaches agree on the matches: QTLtools MBV and GATK ASEReadCounter at
exonic SNPs (share of RNA alleles that match each strain's genotype).

```
analyses/01_rna_genotype_matching/
  Snakefile   metadata -> STAR -> RSEM -> expression PCA ; genotypes -> MBV / ASE -> similarity -> mixup summary
              -> genotype subsets -> Pantry inputs
  scripts/    01_ ... 11.4_ in run order (R scripts set the working directory themselves); MBV.sh (genotype prep and
              first MBV test), MBV.R (per-library MBV plots), get_genotypes_from_tped.R, qc_rna_to_geno_similarity.py,
              assemble_expression_modified.py, manplot.R; run_script_notes.sh = the interactive notes of the 2024-25 run
  src/        library sheets (mouse_samples_info.xlsx, 4 batches), library lists, lookup tables, exon regions
  results/    genotype/ (BXD genotypes, MBV / ASE VCFs, per-group subsets), count_rsem/, mbv_result/,
              ase_counter_result/ (rna_to_geno_summary.tsv, mixup_table_forcheck.csv), pantry_inputs/ (FASTQ maps,
              sample lists, batches, PLINK genotypes per group = what 03 reads)
  figures/    genotype and expression PCA, MBV plots
  intermediate/  PCA, count tables, covariates; bams/ (not kept: rerunning re-aligns)
```

Libraries kept for Pantry: RNA-to-genotype similarity >= 0.9, one library per strain per group
(scripts 11.1-11.4). The BAMs are temporary: they are removed once RSEM, MBV and ASE are done.

Run: `snakemake -s analyses/01_rna_genotype_matching/Snakefile --profile profiles/slurm` (from the project root).
