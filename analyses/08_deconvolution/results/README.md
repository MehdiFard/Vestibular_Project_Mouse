# Cell-type deconvolution of bulk vestibular RNA-seq

Goal: estimate the cell-type composition of every bulk RNA-seq sample for cell-type **interaction eQTL** mapping
(GTEx / Kim-Hellmuth 2020 style; `../../09_interaction_qtl_mapping`). This folder covers the deconvolution step.

## Design
- **Pooling**: ages pooled within each organ, organs kept separate
  - saccule = saccule_w8 + saccule_w22
  - utricle = utricle_w8 + utricle_w22
- **Cell types (6, abundant)**: Hair_cell (I/II collapsed), Supporting_cell,
  Transitional_epith, Glia_Schwann, Fibrocyte_mesench, Roof_dark_cell
- **Reference**: cleaned snRNA atlas (`all_integrated_clean.rds`), 4 donors
  (FTL_1146/1147/1148/1149), downsampled to <=2000 cells/type
- **Methods**: Bisque, MuSiC (cell-type fractions) and ssGSEA marker
  enrichment (GTEx-style). Genes restricted to protein-coding, minus mt-/ribosomal.

## Scripts (`../scripts/`)
| step | script | env | what |
|---|---|---|---|
| 01 | `01_build_bulk_tpm.py`     | Env1       | pool kallisto gene TPM per organ, ENSMUSG→symbol |
| 02 | `02_build_sc_reference.R`  | seurat_env | counts + metadata + markers from the snRNA atlas |
| 03 | `03_run_deconv.R`          | deconv_env | Bisque + MuSiC + ssGSEA per organ |
| 03b| `03b_pseudobulk_check.R`   | deconv_env | validate mechanics on pseudobulk (known truth) |
| 04 | `04_validate_compare.R`    | deconv_env | comparison + validation plots, correlation table |

Run all: `snakemake -s analyses/08_deconvolution/Snakefile --profile profiles/slurm` (from the project root)
(`deconv_env` = isolated conda env: BisqueRNA[github cozygene/bisque], MuSiC[github
xuranw/MuSiC], GSVA, Biobase, SingleCellExperiment.)

## Outputs
- `inputs/`  bulk_tpm_<tissue>.csv, sc_counts.rds, sc_meta.rds, sc_markers.*, sc_composition.csv
- `tables/`  fractions_bisque_<tissue>.csv, fractions_music_<tissue>.csv,
  enrichment_ssgsea_<tissue>.csv, method_correlation_summary.csv, pseudobulk_validation.csv
- `plots/`   00_pseudobulk_validation, 01_mean_composition_by_method, 02_bisque_vs_music_scatter,
  03_bisque_vs_ssgsea_scatter, 04_fraction_boxplots, 05_hair_vs_supporting, 06_ssgsea_enrichment

## Validation
`03b_pseudobulk_check.R` deconvolves pseudobulk built from the reference cells (known composition) to check the
setup; `04_validate_compare.R` compares the three methods on the real bulk. The ssGSEA enrichment scores are rank
based within each sample, so they are less sensitive than the fractions to the difference between the
single-nucleus reference and the whole-tissue bulk RNA; they are the covariate used downstream.

## Used by the interaction eQTL mapping (09)
- Use `enrichment_ssgsea_<tissue>.csv` as the interaction covariate, **inverse-normal
  transformed** across samples (per GTEx).
- Apply the MAF≥0.05-in-both-halves-of-the-enrichment-distribution safeguard.
- Correct with eigenMT (effective tests per gene) + Benjamini-Hochberg across genes; FDR<0.05.
