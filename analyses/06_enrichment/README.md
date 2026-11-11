# Pathway analysis

Two pathway analyses of the same five differential-expression contrasts (A–E) from `../05_dge`.
Run: `snakemake -s analyses/06_enrichment/Snakefile --profile profiles/slurm` (from the project root).
Before October 2026 this folder was `Pantry/sac_vs_utr/pathway_analysis/`.

| | What it does | Report |
|---|---|---|
| [**`gsea/`**](gsea) | GSEA — ranks all ~17,000 genes by Wald statistic, asks whether a set clusters at one end | 24 slides |
| [**`over_rep/`**](over_rep) | ORA — takes only the DE genes, hypergeometric test | 18 slides |

**Start with `gsea/`.** It is the better tool for this dataset; `over_rep/` is there for
comparison and documents why.

Both use the same mouse-native MSigDB 2026.1.Mm gene sets (Hallmark, GO:BP, Reactome),
keyed by Ensembl ID, with the same FDR < 0.05 rule — so their results are directly
comparable. `over_rep/scripts/02_figures.R` reads the GSEA output to build that comparison,
so **run `gsea/` first**.

## Headline biology

1. At 8 weeks the organs differ **metabolically** — utricle higher in oxidative
   phosphorylation and mTORC1; saccule higher in epithelial-mesenchymal transition.
2. By 22 weeks the difference is **inflammatory** — TNFA/NFkB, complement, interferon
   gamma all higher in saccule at w22 but not at w8.
3. Ageing is inflammatory in both organs, **far more so in saccule** (IFN-γ NES +3.12 vs +2.42).
4. The utricle ages differently in kind — it **loses** oxidative phosphorylation (NES −3.28).
5. That asymmetry is what the organ × age interaction detects, and why 36 of its 37
   significant Hallmark sets are negative.

## Shared infrastructure

```
analyses/06_enrichment/
├── Snakefile
├── gsea/        GSEA: scripts, tables, figures, report, README
├── over_rep/    ORA:  scripts, tables, figures, report, README
└── logs/        SLURM logs and the package-install log
```

Mouse annotation (`org.Mm.eg.db` 3.22.0, `msigdbr` 26.1.1) is installed in the **`enrich_r`** conda environment
(October 2026; until then it sat in a project-local `rlib/`, same versions).

Two interpreters: the analysis runs in **`enrich_r`** (R 4.5.3), the decks build in
**`Env1`** (the only env with `python-pptx`).

`msigdbr` downloads its database on first use and caches it under `~/.cache/R/msigdbr`;
the scripts set `TMPDIR` inside the project so that cache write succeeds — the default
`/tmp` is on a different device and the rename fails.
