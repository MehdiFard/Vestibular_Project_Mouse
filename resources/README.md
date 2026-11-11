# resources/ (not in git)

References shared by several analyses (config.yaml `refs:`); analysis-only references are in `analyses/*/src/`.

```
resources/
  genome/   Mus_musculus.GRCm38.dna_sm.primary_assembly.fa (+ .fai), Mus_musculus.GRCm38.102.chr.gtf (Ensembl 102);
            used by 01, 02, 05 and by Pantry, which links to them from analyses/03_pantry/phenotyping/input/ref/
```

The STAR and RSEM indexes and the collapsed GTF of 01 stay outside the project (config_lab.yaml `refs:`).
The snRNA-seq reference (GRCm39, Cell Ranger refdata-gex-GRCm39-2024-A) is in `analyses/07_single_cell_rna/src/`.
