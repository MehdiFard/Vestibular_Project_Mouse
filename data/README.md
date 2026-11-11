# data/ (not in git)

Raw inputs: links to the canonical copies in `projects_data/vestibular_project_data/mouse/` (config.yaml `data:`).

```
data/
  fastq_bulk/          -> rnaseq_fastq/bulk: batch1-4 paired-end FASTQs (<library>_S<n>_L<lane>_R{1,2}_001.fastq.gz)
                          and bulk_rnaseq_metadata_mouse_b1_to_b4.xlsx (submission workbook: library, strain, organ,
                          age, sex, animal)
  fastq_scrna/         -> rnaseq_fastq/single_cell: batch1 snRNA-seq FASTQs + single_cell_rnaseq_metadata_mouse_b1.xlsx
  genotype/            -> genotype: BXD-panel genotypes (all_strains_mm10.tped, mouse_genotype_mm10.{tped,tfam})
  gwas_results/        -> gwas_results: mouse vestibular GWAS results as received (6 traits, .gwas + plots)
  phenotype/           -> phenotype: young-mice VsEP + beam test workbook, one row per mouse
                          (Young_mice_VsEP_Beam_test_010923_strain_names_fixed_20230110_CP.xlsx; 10_gwas)
```
