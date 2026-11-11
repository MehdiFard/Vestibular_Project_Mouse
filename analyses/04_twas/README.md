# 04 TWAS

FUSION TWAS association + coloc of the Pantry pheast weights (`../03_pantry/pheast/intermediate/<group>/twas`) with
the six mouse vestibular GWAS traits (`../03_pantry/pheast/input/gwas_sumstat`, from 02), per group, modality, trait
and chromosome (1-20 as in the LD reference).

```
analyses/04_twas/
  Snakefile                 runs run_twas_with_tissue.smk inside this folder
  run_twas_with_tissue.smk  assoc_test_chr -> assemble_assoc_tests -> assemble_hits_all_traits_modalities
  run_twas_assoc.slurm      the original launcher
  <group>/intermediate/<modality>/<trait>/   one FUSION table per chromosome (correct, per group)
  <group>/output/<modality>/                 <modality>.<trait>.tsv and twas_hits.tsv (P < 5e-8 / 7)
```

Was `pheast/twas_assoc/` inside the Pantry run (the human project keeps TWAS as its own analysis, too).

## Fixed October 2026

- **Groups merged in the tables:** `assemble_assoc_tests` and `assemble_hits_all_traits_modalities` expanded their
  inputs over all four groups (`dataset=datasets`), so every group's `<modality>.<trait>.tsv` held the tests of all
  groups and the four `twas_hits.tsv` were identical (the per-chromosome files were right). Both rules now read their
  own group only (`{{dataset}}`); the hits filter splits on tabs.
- **TWAS.Z sign:** the weights are recomputed with the FUSION sign fix of `../03_pantry` (all modalities), so the
  enet / top1 TWAS.Z have the right direction; no correction afterwards is needed.

Hits after the October 2026 rerun (P < 5e-8 / 7): saccule w8 38, saccule w22 21, utricle w8 36, utricle w22 31
(126 associations, 103 gene-trait pairs, 91 genes). Each `<group>/output/<modality>/twas_hits.tsv` holds the hits of all
seven modalities of that group (the same table in every modality folder).

Counts quoted from the old tables (e.g. the 2026 funder report's 293 associations) were re-derived from the fixed run;
the old hit tables are in `BACKUP/earlier_versions/pre_fix_2026-10/analyses/04_twas/`.

Run: `snakemake -s analyses/04_twas/Snakefile --profile profiles/slurm` (from the project root).
