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

## Notes

- `assemble_assoc_tests` and `assemble_hits_all_traits_modalities` read their own group only (`{{dataset}}`), so each
  group's tables hold that group's tests; the hits filter splits on tabs.
- The weights come from `../03_pantry` with the FUSION sign fix (A1 dosages), so the enet / top1 TWAS.Z have the
  right direction without any correction afterwards.
- Each `<group>/output/<modality>/twas_hits.tsv` holds the hits of all seven modalities of that group (the same table
  in every modality folder).

Run: `snakemake -s analyses/04_twas/Snakefile --profile profiles/slurm` (from the project root).
