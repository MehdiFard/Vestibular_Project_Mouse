# 07 Single-nucleus RNA-seq

Four 10x snRNA-seq libraries of saccule and utricle (FTL_1146-FTL_1149), GRCm39 / Cell Ranger 9.0.1 with introns.

```
analyses/07_single_cell_rna/
  Snakefile   FastQC + Cell Ranger, then annotation -> integration -> per organ -> organ DE -> clean re-run -> reports
  scripts/    annot2 (FTL_1147 marker annotation) -> annot_hc_subcluster (hair-cell types) -> annot3 (label transfer,
              Burns 2015) -> annot4 (all samples integrated) ; annot5 (per organ) -> annot6 (cell-level organ DE) ;
              annot7/8 (pseudobulk organ DE, edgeR) ; clean re-run: annot9 (SoupX + doublets) -> annot10 (integrate)
              -> annot11 (compare with the original) -> annot12/13 (clean pseudobulk DE) ; recolor_*.R (figures),
              vestibular_markers.R, celltype_palette.R, make_metadata.R ; run_*.sbatch (the original launchers) ;
              run_script_notes.sh (interactive notes)
  src/        Cell Ranger, refdata-gex-GRCm39-2024-A, sample sheet, rna_metadata.tsv, annot_ref/ (marker references,
              gene biotypes)
  data/       links to the raw FASTQs (../../data/fastq_scrna)
  intermediate/  cell_ranger_output/, fastqc/
  results/    cell_annotation/, per_tissue/, tissue_DE/, pseudobulk_DE/, cell_annotation_clean/ (clean objects,
              compositions, recovery table, clean pseudobulk DE), report_assets/
```

The clean integration (`results/cell_annotation_clean/all_integrated_clean.rds`) is the reference of
`../08_deconvolution`.

Run: `snakemake -s analyses/07_single_cell_rna/Snakefile --profile profiles/slurm` (from the project root).
