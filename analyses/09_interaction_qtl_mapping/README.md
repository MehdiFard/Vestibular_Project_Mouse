# 09 Cell-type interaction eQTL

tensorQTL cis-interaction mapping of bulk expression with the ssGSEA cell-type enrichment (`../08_deconvolution`)
as the interaction term: per group with Pantry's covariates (`ieqtl/`) or genotype PCs only (`ieqtl_core/`), and per
organ with both ages pooled (`ieqtl_pooled/<organ>_{core,full}`); compared with Pantry's marginal eQTLs and across
approaches (`mapping_report/`).

```
analyses/09_interaction_qtl_mapping/
  Snakefile   ieqtl_group / ieqtl_group_core (GPU) -> compare_marginal ; pool_genotypes -> pool_expression ->
              ieqtl_pooled (GPU) -> compare_pooled ; master_compare ; reports
  scripts/    05_interaction_eqtl.py, 06_ieqtl_vs_eqtl_compare.py, 07a_pool_genotypes.sh, 07b_pool_expression_covar.py,
              08_interaction_eqtl_pooled.py, 09_pooled_compare.py, 10_master_compare.py, build_*.py,
              run_ieqtl_*_gpu.sbatch (the original launchers; submit from this folder)
  results/    ieqtl/, ieqtl_core/, ieqtl_pooled/, pooled_inputs/, mapping_report/ (each with its deck / document)
```

The first ieQTL run (`ieqtl_core_v1`, `ieqtl_pooled_v1`) is in `BACKUP/earlier_versions/09_interaction_qtl_mapping/`.

Rerun October 2026 from `../03_pantry` and the rerun `../08_deconvolution` (the earlier results came from the
superseded plain Pantry run). The per-group "full" model uses Pantry's covariates, which in the final run are 3 genotype
PCs + 8 expression PCs (the plain run's also had batch); the pooled models add age and batch
(`pheast/input/<group>/<group>.batches.tsv`). ieGenes (master_ieGene_counts.csv), before -> after: per-group full
42 -> 56, per-group core 13 -> 16, pooled core 79 -> 76, pooled full 120 -> 123.

Run: `snakemake -s analyses/09_interaction_qtl_mapping/Snakefile --profile profiles/slurm` (from the project root; GPU partition in
config.yaml `slurm.gpu_partition`).
