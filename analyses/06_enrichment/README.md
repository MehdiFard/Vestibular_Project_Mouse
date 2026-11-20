# 06 Pathway analysis

Two pathway analyses of the same five differential-expression contrasts (A–E) from `../05_dge`.
Run: `snakemake -s analyses/06_enrichment/Snakefile --profile profiles/slurm` (from the project root).

| | What it does |
|---|---|
| [**`gsea/`**](gsea) | GSEA — ranks all ~17,000 genes by Wald statistic, asks whether a set clusters at one end |
| [**`over_rep/`**](over_rep) | ORA — takes only the DE genes, hypergeometric test |

**Start with `gsea/`**, the primary analysis; `over_rep/` is there for comparison.

Both use the same mouse-native MSigDB 2026.1.Mm gene sets (Hallmark, GO:BP, Reactome),
keyed by Ensembl ID, with the same FDR < 0.05 rule — so their results are directly
comparable. `over_rep/scripts/02_figures.R` reads the GSEA output to build that comparison,
so **run `gsea/` first**.

## Shared infrastructure

```
analyses/06_enrichment/
├── Snakefile
├── gsea/        GSEA: scripts, tables, figures, README
├── over_rep/    ORA:  scripts, tables, figures, README
└── logs/        SLURM logs and the package-install log
```

The analysis runs in the **`enrich_r`** conda environment (R 4.5.3, `org.Mm.eg.db` 3.22.0, `msigdbr` 26.1.1).

`msigdbr` downloads its database on first use and caches it under `~/.cache/R/msigdbr`;
the scripts set `TMPDIR` inside the project so that cache write succeeds — the default
`/tmp` is on a different device and the rename fails.
