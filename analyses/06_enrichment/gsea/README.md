# Pathway analysis — GSEA across the five saccule/utricle contrasts

Gene-set enrichment for every contrast in `../` (A–E), run with **fgsea** on the
DESeq2 Wald statistics.

```
pathway_analysis/
├── scripts/
│   ├── 01_run_gsea.R          fgsea, 5 contrasts x 3 collections
│   ├── 02_figures.R           all 28 figures
│   ├── 03_build_report.py     the deck
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
├── report/pathway_analysis_report.pptx    24 slides
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

## Results

| Contrast | Hallmark | GO:BP | Reactome | Strongest Hallmark set |
|---|---|---|---|---|
| A — Utricle vs Saccule, w8 | 10 | 127 | 28 | Epithelial-mesenchymal transition (−2.23) |
| B — Utricle vs Saccule, w22 | 20 | 605 | 98 | Epithelial-mesenchymal transition (−2.53) |
| C — w22 vs w8, Saccule | 21 | 660 | 86 | Interferon gamma response (**+3.12**) |
| D — w22 vs w8, Utricle | 20 | 364 | 183 | Oxidative phosphorylation (**−3.28**) |
| E — organ × age interaction | 37 | 1,203 | 343 | Oxidative phosphorylation (−3.07) |

Counts are gene sets at FDR < 0.05.

### The biology, in five lines

1. **At 8 weeks the organs differ metabolically.** Utricle is higher in oxidative
   phosphorylation, mTORC1 signalling and cholesterol homeostasis; saccule is higher in
   epithelial-mesenchymal transition and apical junction.
2. **By 22 weeks the difference has become inflammatory.** TNFA/NFkB, complement,
   interferon gamma and allograft rejection are all significantly higher in saccule at
   w22 but *not* at w8.
3. **Ageing is inflammatory in both organs, far more so in saccule.** Interferon gamma is
   the strongest result in the analysis (NES +3.12, FDR 6.5e-30) against +2.42 in utricle.
4. **The utricle ages differently in kind** — it loses oxidative phosphorylation
   (NES −3.28) and MYC targets, which the saccule does not.
5. **That asymmetry is what the interaction detects**, and why 36 of its 37 significant
   Hallmark sets are negative: saccule gains inflammation faster *and* utricle loses
   metabolism faster, and both push the organ gap the same way.

Epithelial-mesenchymal transition is the only Hallmark set significant in all five contrasts.

### Contrast similarity (GO:BP NES, Spearman)

Three clean blocks: the organ contrasts agree with each other (A–B **0.73**), the ageing
contrasts agree with each other (C–D **0.71**), and organ and age are independent
(0.04 to −0.16). E anti-correlates with C (**−0.57**) because the saccule's ageing
response is what drives the interaction.

## Caveats

- **GO:BP is redundant.** Many of its 7,781 sets share most of their genes, so 1,203
  significant sets is not 1,203 independent findings. **Quote Hallmark.**
- **Contrast E is fitted without strain**, unlike A and B, so its estimates are noisier
  than the organ contrasts and its gene-level list is the least reproducible result in
  the project (37% agreement with limma).
- **C and D compare mostly different animals across ages**, with batch partly confounded
  with age. The inflammatory signal survives that, but effect sizes are not precise.
- **GSEA tests placement in a ranking**, not pathway activity or causation.
- **Immune/interferon enrichment in bulk tissue usually reflects cell composition**
  (infiltrating immune cells) as much as regulation within a cell type. The mast-cell
  proteases at the top of the interaction DEG list point the same way.

## Environment

The GSEA runs in the **`enrich_r`** conda env (R 4.5.3, fgsea, clusterProfiler, org.Mm.eg.db 3.22.0,
msigdbr 26.1.1).

The deck is built with **`Env1`** (the only env with `python-pptx`). Two interpreters:

```bash
ENRICH="$(conda info --base)/envs/enrich_r"
ENV1="$(conda info --base)/envs/Env1"

sbatch scripts/run_gsea.slurm          # or: $ENRICH/bin/Rscript scripts/01_run_gsea.R
$ENRICH/bin/Rscript scripts/02_figures.R
$ENV1/bin/python  scripts/03_build_report.py
```

`msigdbr` downloads the gene-set database on first use and caches it under
`~/.cache/R/msigdbr`; the scripts set `TMPDIR` inside the project so that cache write
succeeds (the default `/tmp` is on a different device and the rename fails).
