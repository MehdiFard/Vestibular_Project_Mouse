# 10 GWAS of the young-mice vestibular traits (FaST-LMM)

In-house GWAS of four young-mice traits (VsEP threshold, P2 latency, P1 amplitude, beam test; 85-86 inbred strains:
66 BXD + 19 other inbred strains incl. the BXD parents; 444-544 mice) with FaST-LMM: linear mixed model,
leave-one-chromosome-out kinship, without (`noPCs`) and with (`withPCs`) the first `N_PCS` genotype PCs as covariates.
It is the 2023 pipeline (old `VestProj/GWAS`, run by hand; kept in `BACKUP/kept_from_old_folders/`) as a Snakefile.
The GWAS that 02 and 04 use are the ones received from UCLA (`../../data/gwas_results`; the four traits are copied to
`input/gwas_ucla_results/`); `compare_with_received` checks ours against them, and nothing downstream reads these results.

```
analyses/10_gwas/
  Snakefile   phenotypes -> genotypes -> genotype_qc -> pca -> fastlmm (noPCs, withPCs) -> manhattan_qq;
              compare_with_received
  scripts/    01_make_phenotypes.R     workbook sheets 1-4 -> intermediate/phenotype/<trait>.phen (IID = <strain>_<n>)
              02_make_genotypes.R      each mouse gets its strain's genotype -> intermediate/genotype/<trait>.{tped,tfam}
              03_pca.R                 strain-level PCA, % variance, covariates, scree + PC1/PC2 plots
              04_fastlmm.py            FaST-LMM LOCO, one trait, one model (env gwas_fastlmm)
              05_manhattan_qq.R        Manhattan + QQ (genomic inflation in the title)
              06_compare_with_received.R  same chr:bp SNPs as the received GWAS: correlation of -log10 P, top hits
  input/      raw inputs: links to the workbook and the genotypes (../../data), gwas_ucla_results/ (copy); see input/README.md
  intermediate/  phenotype/ (<trait>.phen, .pheno.txt), genotype/ (per-trait tped/tfam, QC'd bed sets),
              genotype_pca/ (one mouse per strain, pruned SNPs, .raw for the PCA)
  results/    genotype_pca/ (<trait>.eigenvec.tsv, .variance.tsv, .covar.pcs), noPCs/<trait>.gwas.loco.tsv,
              withPCs/<trait>.gwas.loco.tsv, compare_with_received.tsv
  logs/       one log per rule and trait (logs/<rule>/)
  figures/    pca/<trait>_{scree,PC1_PC2}.jpeg, noPCs/ and withPCs/<trait>.{manhattan,QQ}.jpeg,
              compare_with_received.jpeg
```

- Genotype QC: SNPs only (ACGT), call rate >= 0.9, MAF >= 0.05, REF/ALT from GRCm38, no Y; chr X coded 20 for
  FaST-LMM. Gives exactly the SNP sets of the received GWAS.
- PCA: one mouse per strain (mice of a strain have the same genotype), SNPs LD-pruned among the strains (window 50,
  step 5, r2 0.5), standardised; `prcomp` in R, so % variance is of the total. Each mouse gets its strain's scores.
- Changes from the 2023 run: the PCA was on mice (strains weighted by their number of mice, 1-22) with PC1-2 and %
  variance of the 10 computed PCs only; files were paired with traits by sorted order; plink2 now gets `--mouse`.

Run: `snakemake -s analyses/10_gwas/Snakefile --profile profiles/slurm` (from the project root).
