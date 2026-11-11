# 08 Cell-type deconvolution

Bulk saccule / utricle expression deconvolved with the clean snRNA-seq reference (`../07_single_cell_rna`): Bisque and MuSiC
(fractions) and ssGSEA marker enrichment, validated on pseudobulk made from the reference. The ssGSEA enrichment is
the interaction covariate of `../09_interaction_qtl_mapping`.

```
analyses/08_deconvolution/
  Snakefile   bulk_tpm -> sc_reference -> deconvolution (3 methods) -> pseudobulk_check -> compare_methods -> report
  scripts/    01_build_bulk_tpm.py, 02_build_sc_reference.R (seurat_env), 03_run_deconv.R, 03b_pseudobulk_check.R,
              04_validate_compare.R (deconv_env), make_intro_assets.py, build_deconv_deck.py / _doc.py,
              run_deconv.sbatch (the original launcher; submit from this folder)
  results/    inputs/ (bulk TPM, sc counts / meta / markers), tables/ (fractions, enrichment, validation), plots/,
              intro_assets/, README.md (description of the outputs)
  report/     deconvolution_report.pptx, Deconvolution_Report.docx
```

Rerun October 2026 from the expression of `../03_pantry` (the earlier results came from the superseded plain Pantry
run; same samples and genes, r = 0.9999). The cell-type colours come from `../07_single_cell_rna/scripts/celltype_palette.R`.

Run: `snakemake -s analyses/08_deconvolution/Snakefile --profile profiles/slurm` (from the project root).
