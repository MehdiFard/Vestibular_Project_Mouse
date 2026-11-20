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
              08_interaction_eqtl_pooled.py, 09_pooled_compare.py, 10_master_compare.py,
              run_ieqtl_*_gpu.sbatch (the original launchers; submit from this folder)
  results/    ieqtl/, ieqtl_core/, ieqtl_pooled/, pooled_inputs/, mapping_report/
```

Covariates: the per-group "full" model uses Pantry's covariates (genotype PCs + expression PCs), the "core" model
genotype PCs only; the pooled models add age and batch (`../03_pantry/pheast/input/<group>/<group>.batches.tsv`).

Run: `snakemake -s analyses/09_interaction_qtl_mapping/Snakefile --profile profiles/slurm` (from the project root; GPU partition in
config.yaml `slurm.gpu_partition`).
