# 05 Saccule vs utricle: PCA and differential expression

Expression comparison of the two mouse vestibular otolith organs (saccule, utricle) at 8 and 22 weeks, built on the
kallisto quantifications of the Pantry phenotyping run (`../03_pantry`).

Run (from the project root; the DE steps go to Slurm):

```bash
snakemake -s analyses/05_dge/Snakefile --profile profiles/slurm
```

## Inputs

| | |
|---|---|
| Libraries | the Pantry QTL set: one library per strain per group (`saccule_w8`, `saccule_w22`, `utricle_w8`, `utricle_w22`) |
| Counts | kallisto `est_counts` from `../03_pantry/phenotyping/intermediate/<group>/expression/`, summed transcript -> gene with **tximport** |
| Annotation | `Mus_musculus.GRCm38.102.chr.gtf` (Ensembl 102) |
| Organ / age | the Pantry group of each library |
| Batch | `../03_pantry/pheast/input/<group>/<group>.batches.tsv` (sequencing plates 1-4) |
| Strain | RNA-to-genotype concordance from `../01_rna_genotype_matching` |
| Animal, sex | submission workbook (`data/fastq_bulk/bulk_rnaseq_metadata_mouse_b1_to_b4.xlsx`); sex also inferred from expression |

Pantry's own `expression.bed.gz` is quantile-normalised TPM and not a valid DESeq2 input, which is why counts are
re-imported from the kallisto runs.

## Design

Saccule and utricle come from the **same animal**. The animal is taken from the submission workbook (the two organs of
an animal have consecutive library IDs, the same strain and age); where the workbook has no such pair, strain + age +
sex is used. Within any one organ/age group each strain appears once.

| Key | Contrast | Model |
|---|---|---|
| A | Utricle vs saccule at w8 | `~ strain + organ` (paired) |
| B | Utricle vs saccule at w22 | `~ strain + organ` (paired) |
| C | w22 vs w8 in saccule | `~ batch + age` |
| D | w22 vs w8 in utricle | `~ batch + age` |
| E | Organ x age interaction | `~ batch + organ + age + organ:age` |
| S | Unpaired organ contrasts | `~ batch + organ` (sensitivity check) |

Positive `log2FoldChange` = higher in **utricle** (A, B), higher at **w22** (C, D), or the organ gap **widening** with
age (E). Significance: `padj < 0.05` **and** `|log2FC| > 1`. Fold changes are shrunken with `ashr`; the unshrunken
value is kept as `log2FC_raw`.

### Gene filtering

Both filters are set in `scripts/00_config.R` and re-applied by each contrast to its own subset:

1. **Biotype:** `KEEP_BIOTYPES <- c("protein_coding", "lincRNA")`, applied before testing so other biotypes never enter
   the FDR correction. In Ensembl 102 `lincRNA` is the intergenic class only; add `antisense`, `processed_transcript`,
   `bidirectional_promoter_lncRNA` and `sense_intronic` for lncRNA in the broad sense.
2. **Expression:** >= 10 counts in at least as many libraries as the smallest group.

### Sex

Inferred from expression: **Xist** (female) against the Y genes **Ddx3y / Uty / Eif2s3y / Kdm5d** (male), in
`scripts/01_build_inputs.R`; `scripts/13_infer_sex.R` checks the calls against the recorded sex and draws the figure.

## Checks and sensitivity analyses

| Script | Question |
|---|---|
| `06_batch_correction_check.R` | why batch is in the design rather than corrected out first |
| `07_biotype_separation_check.R` | what the biotype filter costs |
| `08_limma_crosscheck.R` | limma-voom rerun of every contrast |
| `09_vst_and_ntop_check.R` | what VST does; PCA top-N sensitivity |
| `10_pairing_and_paired_age.R` | what pairing buys; age contrast within strain (`~ strain + batch + age`) |
| `11_model_diagram.R` | what each model formula compares |
| `14_sex_effect.R` | sex effect and a female-only age contrast |
| `15_paired_vs_unpaired_check.R` | paired vs unpaired organ contrasts |
| `16_interaction_illustration.R` | what the organ x age interaction tests |
| `17_validate_against_backup_sheet.py` | every annotation checked against the complete submission workbook |

`rebuild_291/` repeats the analysis as a sensitivity check, adding the libraries that passed the genotype match but
were left out by the one-library-per-strain-per-group rule (quantified with the same kallisto index and parameters).

## Layout

```
analyses/05_dge/
├── Snakefile
├── src/                               sample metadata (strain, sex, plate)
├── rebuild_291/                       sensitivity re-analysis (own scripts/)
├── scripts/
│   ├── 00_config.R                    paths, thresholds, colours
│   ├── 01_build_inputs.R              sample table, annotation, tximport matrix, sex calls
│   ├── 02_pca.R                       PCA, scree, covariate R², batch and sex checks
│   ├── 03_differential_expression.R   the DESeq2 analyses
│   ├── 04_figures.R                   volcano, concordance, top-gene figures
│   ├── 06-17_*                        checks (table above); 12_build_master_metadata.R = master metadata file
│   ├── run_de.slurm, run_pca.slurm, run_paired_check.slurm   Slurm jobs (set --partition / --account)
│   └── run_all.sh                     the original launcher
├── data/                              count matrix, sample table, cached objects
├── results/tables/                    full and significant-only results (tsv)
├── results/figures/                   300 dpi PNGs
└── logs/
```

## Environment

`Env1` (`envs/Env1.yml`) has DESeq2, tximport and the Python packages. DESeq2 fits one small GLM per gene, so the
BLAS is pinned to a single thread and parallelism is over genes with `BiocParallel`. Run the DE step on a compute
node (`run_de.slurm`).
