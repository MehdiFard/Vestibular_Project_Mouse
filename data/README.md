# data/ (not in git)

Raw inputs of every analysis (config.yaml `data:`). Analyses link to these files; nothing here is modified.

```
data/
  fastq_bulk/          batch1-4 paired-end FASTQs (<library>_S<n>_L<lane>_R{1,2}_001.fastq.gz) and
                       bulk_rnaseq_metadata_mouse_b1_to_b4.xlsx (submission workbook: library, strain, organ, age, sex,
                       animal)
    run_info/          sequencing-run files: vendor md5sum.txt, demultiplexing Stats/ and Reports/ (batch3, batch4),
                       undetermined reads (batch4), the original sample sheet (mouse_samples_metadata.xlsx)
  fastq_scrna/         batch1 snRNA-seq FASTQs + single_cell_rnaseq_metadata_mouse_b1.xlsx
    run_info/batch1/   vendor md5sum.txt, demultiplexing Stats/ and Reports/
  genotype/            BXD-panel genotypes (all_strains_mm10.tped, mouse_genotype_mm10.{tped,tfam})
  gwas_results/        mouse vestibular GWAS results as received (6 traits, .gwas + plots)
  phenotype/           young-mice VsEP + beam test workbook, one row per mouse
                       (Young_mice_VsEP_Beam_test_010923_strain_names_fixed_20230110_CP.xlsx; 10_gwas)
```

Some FASTQs were renamed from the vendor names (barcode names -> FTL_1017-1021, FTL_1046 S52); md5 values in
`run_info/*/md5sum.txt` are those of the vendor files.
