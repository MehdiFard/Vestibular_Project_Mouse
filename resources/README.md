# resources/ (not in git)

References and tools shared by several analyses (config.yaml `refs:` and `tools:`); analysis-only references are in
`analyses/*/src/`.

```
resources/
  genome/       Mus_musculus.GRCm38.dna_sm.primary_assembly.fa (+ .fai, .dict), Mus_musculus.GRCm38.102.chr.gtf
                (Ensembl 102), Mus_musculus.GRCm38.102.chr.collapsed.genes.gtf (GTEx-style gene collapse; 01);
                used by 01, 02, 05 and by Pantry, which links to them from analyses/03_pantry/phenotyping/input/ref/
  index_star/   STAR 2.7.11b index of the genome (01)
  index_rsem/   RSEM index rsem_indx_m38 (01)
  tools/        STAR 2.7.11b and QTLtools binaries (01)
```

The snRNA-seq reference (GRCm39, Cell Ranger refdata-gex-GRCm39-2024-A) and Cell Ranger 9.0.1 are in
`analyses/07_single_cell_rna/src/`.
