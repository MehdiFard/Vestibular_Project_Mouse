# 03 Pantry

Multimodal RNA phenotyping and xQTL / heritability / TWAS-weight mapping with Pantry (Munro et al. 2024) for the four
groups `saccule_w8`, `saccule_w22`, `utricle_w8`, `utricle_w22`. Pantry's own pipelines run inside `phenotyping/` and
`pheast/`; `Snakefile` here runs them in order. TWAS association with the GWAS traits is `../04_twas`.
Pantry's own README: `README_Pantry.md`.

```
analyses/03_pantry/
  Snakefile        phenotyping -> cross-modality -> pheast inputs (+ editing regression) -> pheast (GPU), per group
  phenotyping/     7 modalities (alt_polyA, alt_TSS, expression, isoforms, RNA_editing, splicing, stability);
                   config.yml (datasets = groups), input/<group>/ (FASTQ map + samples from 01; FASTQ links),
                   input/ref/ (GTF, genome -> ../../../resources/genome; mouse editing sites mm10)
  pheast/          cis-QTL (tensorQTL), heritability, TWAS weights (FUSION); input/<group>/ (genotypes, batches from
                   01; phenotypes/ from phenotyping/output), input/{gwas_sumstat,LDREF} from 02
  develop_notes.sh interactive notes of the runs (environments, commands per group)
  prep.R           (Pantry upstream)
```

This is the final run (folder `Pantry_RNA_editing` before October 2026). The earlier run without RNA editing (six
modalities) is superseded and was moved out of the project (`vestibulart_project_mouse_superseded/Pantry_plain_run`,
listed for deletion).

## RNA editing and TWAS weights (fixed October 2026, as in the human project)

1. **0/0 fix** (`phenotyping/scripts/assemble_bed.py`): samples below the coverage cutoff (`0/0` in the edit matrix,
   coverage < 10; about 16% of the values) were scored (0+0.5)/(0+0.5) = 100% edited; they are now missing, so the
   >40%-missing filter and the mean imputation apply.
2. **Host-gene expression regressed out** (`pheast/scripts/regress_rna_editing_host_expression.py`, Li et al. 2022):
   each site's normalised editing level is regressed on log2(TPM+1) of its host gene(s); residuals are rank-INT
   normalised. Cross-modality is rebuilt with the regressed editing (`pheast/scripts/build_cross_modality.py`).
   The unregressed (0/0-fixed) inputs stay in `pheast/input/<group>/phenotypes_unregressed/`.
3. **FUSION sign fix** (`pheast/scripts/fusion_twas/FUSION.compute_weights.R`, `bed <- 2 - bed`): snpStats gives A2
   allele counts, FUSION expects A1 dosages. The TWAS weights of **all** modalities were recomputed with it.

Rerun (October 2026): RNA_editing and cross_modality fully (QTL, heritability, TWAS weights); the other six
modalities' TWAS weights only, so their QTL and heritability results are unchanged. The previous weights
(`pheast/intermediate/<group>/twas_pre_fix_2026-10/`, 226 GB) were deleted on 2026-10-29; key tables of the previous
results are in `BACKUP/earlier_versions/pre_fix_2026-10/`.

pheast runs in the Apptainer image of config `pheast_container` when set (as in the July 2025 run).

Run: `snakemake -s analyses/03_pantry/Snakefile --profile profiles/slurm` (from the project root).
