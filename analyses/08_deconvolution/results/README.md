# Cell-type deconvolution of bulk vestibular RNA-seq

Goal: estimate the cell-type composition of every bulk RNA-seq sample so that we
can later run cell-type **interaction eQTL** mapping (GTEx / Kim-Hellmuth 2020 style).
This folder covers the **deconvolution step only**; QTL mapping is planned next.

## Design
- **Pooling**: ages pooled within each organ, organs kept separate
  - saccule = saccule_w8 + saccule_w22 (128 bulk samples)
  - utricle = utricle_w8 + utricle_w22 (132 bulk samples)
- **Cell types (6, abundant)**: Hair_cell (I/II collapsed), Supporting_cell,
  Transitional_epith, Glia_Schwann, Fibrocyte_mesench, Roof_dark_cell
- **Reference**: cleaned snRNA atlas (`all_integrated_clean.rds`), 4 donors
  (FTL_1146/1147/1148/1149), downsampled to <=2000 cells/type
- **Methods**: Bisque (planned primary), MuSiC (sensitivity), ssGSEA marker
  enrichment (GTEx-style). Genes restricted to protein-coding, minus mt-/ribosomal.

## Scripts (scripts/deconv_scripts/, run from project root)
| step | script | env | what |
|---|---|---|---|
| 01 | `01_build_bulk_tpm.py`     | Env1       | pool kallisto gene TPM per organ, ENSMUSG→symbol |
| 02 | `02_build_sc_reference.R`  | seurat_env | counts + metadata + markers from the snRNA atlas |
| 03 | `03_run_deconv.R`          | deconv_env | Bisque + MuSiC + ssGSEA per organ |
| 03b| `03b_pseudobulk_check.R`   | deconv_env | validate mechanics on pseudobulk (known truth) |
| 04 | `04_validate_compare.R`    | deconv_env | comparison + validation plots, correlation table |

Run all: `sbatch scripts/deconv_scripts/run_deconv.sbatch`
(`deconv_env` = isolated conda env: BisqueRNA[github cozygene/bisque], MuSiC[github
xuranw/MuSiC], GSVA, Biobase, SingleCellExperiment.)

## Outputs
- `inputs/`  bulk_tpm_<tissue>.csv, sc_counts.rds, sc_meta.rds, sc_markers.*, sc_composition.csv
- `tables/`  fractions_bisque_<tissue>.csv, fractions_music_<tissue>.csv,
  enrichment_ssgsea_<tissue>.csv, method_correlation_summary.csv, pseudobulk_validation.csv
- `plots/`   00_pseudobulk_validation, 01_mean_composition_by_method, 02_bisque_vs_music_scatter,
  03_bisque_vs_ssgsea_scatter, 04_fraction_boxplots, 05_hair_vs_supporting, 06_ssgsea_enrichment

## KEY FINDINGS (read before using the numbers)
1. **Mechanics validated on pseudobulk** (`00_pseudobulk_validation.png`): deconvolving
   pseudobulk built from the single cells themselves (composition known) gives
   **MuSiC r = 0.95** and **Bisque r = 0.45** vs truth. So the setup is correct, and
   MuSiC — not Bisque — is the stronger method in-domain (opposite to the a-priori plan).
2. **Absolute fractions are NOT reliable on the real bulk.** On real bulk, MuSiC
   collapses (~0.85–0.93 Supporting_cell, little variation) and Bisque returns near-
   uniform proportions; the two barely agree (median Spearman ≈ 0.18). Marker-restricting
   the genes makes MuSiC worse (~1.0 Supporting), not better. Cause: the reference is
   **single-NUCLEUS** RNA while the bulk is **whole-cell / total** RNA (kallisto TPM) —
   a strong domain shift — compounded by an **incomplete 6-type reference** for a
   whole-organ dissection (surrounding bone/blood/nerve are unmodeled, so NNLS dumps
   residual signal into one type). Do **not** interpret the Bisque/MuSiC fractions as
   real composition.
3. **Use the ssGSEA marker-enrichment scores** (`06_ssgsea_enrichment.png`) as the
   per-sample cell-type covariate. Enrichment is rank/within-sample based, so it is
   robust to the nuclei-vs-bulk shift (this is exactly why GTEx used xCell enrichment,
   not fractions). The scores show real cross-sample variation (e.g. Hair_cell 0.01–0.57)
   and are the consensus signal both fraction methods partially track (Bisque–ssGSEA 0.49,
   MuSiC–ssGSEA 0.41). Best-behaved types: Hair_cell, Transitional_epith, Glia_Schwann,
   Roof_dark_cell.

## For the next step (interaction eQTL)
- Use `enrichment_ssgsea_<tissue>.csv` as the interaction covariate, **inverse-normal
  transformed** across samples (per GTEx).
- Apply the MAF≥0.05-in-both-halves-of-the-enrichment-distribution safeguard.
- Correct with eigenMT (effective tests per gene) + Benjamini-Hochberg across genes; FDR<0.05.
