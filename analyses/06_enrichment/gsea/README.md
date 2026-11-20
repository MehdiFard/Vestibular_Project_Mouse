# Pathway analysis — GSEA across the five saccule/utricle contrasts

Gene-set enrichment for every contrast in `../` (A–E), run with **fgsea** on the
DESeq2 Wald statistics.

```
gsea/
├── scripts/
│   ├── 01_run_gsea.R          fgsea, 5 contrasts x 3 collections
│   ├── 02_figures.R           all 28 figures
│   └── run_gsea.slurm         step 1 as a SLURM job
├── tables/
│   ├── gsea_summary.tsv              counts per contrast x collection
│   ├── gsea_all_contrasts.tsv.gz     every set, with leading-edge genes
│   ├── contrast_nes_correlation.tsv  contrast similarity matrix
│   ├── gsea_objects.rds              combined results + gene sets
│   └── hallmark/ go_bp/ reactome/    one table per contrast
├── figures/
│   ├── summary/     3 cross-contrast figures
│   ├── hallmark/    NES bars, dot plots, enrichment curves per contrast
│   └── go_bp/ reactome/  NES bars per contrast
```

## Method

| Choice | What was done | Why |
|---|---|---|
| Ranking | DESeq2 **Wald statistic**, one per gene | carries direction *and* precision; ranking on raw log2FC puts noisy low-count genes at the extremes |
| Pre-filtering | **none** — the whole ranked list | GSEA exists to find sets of genes moving together by a little; thresholding first discards exactly that |
| Gene sets | mouse-native **MSigDB 2026.1.Mm**, keyed by Ensembl ID | no human→mouse ortholog mapping, no ID conversion |
| Collections | Hallmark (50), GO:BP (7,781), Reactome (1,333) | Hallmark for the headline, GO/Reactome for detail |
| Set size | 15–500 genes | standard; avoids tiny and uselessly broad sets |
| Engine | `fgsea` multilevel, `eps = 0`, seed 42 | exact small p-values |
| FDR | Benjamini–Hochberg within each contrast × collection | |

Positive NES always means the set sits at the **top** of that contrast's ranking.
What the top means differs by contrast, so every figure spells it out.

## Caveats

- **GO:BP is redundant.** Many of its sets share most of their genes, so significant GO:BP sets are not independent
  findings; Hallmark is the headline collection.
- **Contrast E is fitted without strain**, unlike A and B, so its estimates are noisier than the organ contrasts.
- **C and D compare mostly different animals across ages**, with batch partly confounded with age.
- **GSEA tests placement in a ranking**, not pathway activity or causation.
- **Enrichment in bulk tissue can reflect cell composition** as much as regulation within a cell type.

## Environment

The GSEA runs in the **`enrich_r`** conda env (R 4.5.3, fgsea, clusterProfiler, org.Mm.eg.db 3.22.0,
msigdbr 26.1.1).

```bash
ENRICH="$(conda info --base)/envs/enrich_r"

sbatch scripts/run_gsea.slurm          # or: $ENRICH/bin/Rscript scripts/01_run_gsea.R
$ENRICH/bin/Rscript scripts/02_figures.R
```

`msigdbr` downloads the gene-set database on first use and caches it under
`~/.cache/R/msigdbr`; the scripts set `TMPDIR` inside the project so that cache write
succeeds (the default `/tmp` is on a different device and the rename fails).
