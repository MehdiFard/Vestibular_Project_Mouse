# Over-representation analysis (ORA)

The classical pathway test, run on the DE gene lists — the counterpart to [`../gsea`](../gsea).
Included for comparison; **GSEA is the better tool for this dataset** (see below).

## The difference in one line

| | Input | Question |
|---|---|---|
| **ORA** (here) | only the DE genes (FDR < 0.05 **and** \|log2FC\| > 1) | is this set over-represented in my list? (hypergeometric) |
| **GSEA** (`../gsea`) | all ~17,000 genes, ranked by Wald statistic | do this set's genes cluster at one end of the ranking? |

## Two choices that matter

1. **Universe** = the genes actually **tested in that contrast**, not the genome and not all
   of MSigDB. A bigger universe than was tested inflates every p-value. This is the most
   common way ORA is got wrong.
2. **Direction** — up and down genes are tested separately as well as together. A set that
   is half up and half down looks enriched in the combined list while being incoherent.

Everything else is identical to the GSEA run (same MSigDB 2026.1.Mm sets, 15–500 genes,
FDR < 0.05), so the two are directly comparable.

## What ORA had to work with

| Contrast | Universe | DE genes | Up | Down |
|---|---|---|---|---|
| A — Utricle vs Saccule, w8 | 17,066 | 756 | 263 | 493 |
| B — Utricle vs Saccule, w22 | 17,143 | 860 | 261 | 599 |
| C — w22 vs w8, Saccule | 17,366 | 168 | 100 | 68 |
| D — w22 vs w8, Utricle | 17,287 | 146 | 94 | 52 |
| E — organ × age interaction | 18,272 | 122 | 41 | 81 |

## Head to head (significant sets, FDR < 0.05)

| Contrast | Hallmark GSEA / ORA | GO:BP GSEA / ORA | Reactome GSEA / ORA |
|---|---|---|---|
| A | 10 / 6 | 127 / **510** | 28 / 26 |
| B | 20 / 6 | 605 / 539 | 98 / 32 |
| C | 21 / 3 | 660 / 80 | 86 / 5 |
| D | 20 / **0** | 364 / 12 | 183 / 1 |
| E | 37 / 4 | 1,203 / 80 | 343 / 12 |

GSEA wins in Hallmark and Reactome in all five. **GO:BP is the exception** — ORA returns
510 terms for contrast A against GSEA's 127, because a 756-gene list hits many small,
heavily overlapping GO terms. That is more *terms*, not more biology.

Per Hallmark set across all contrasts: **14** called by both, **94** by GSEA only,
**5** by ORA only. ORA is close to a subset of GSEA here.

## The worked example — why the gap exists

`figures/summary/03_worked_example_oxphos.png`. Oxidative phosphorylation in contrast D is
the strongest single result in the whole GSEA run (NES −3.28, FDR 2e-37). Of its **194**
genes tested:

- **132 pass FDR < 0.05**
- **88% move in the same direction** (down), median log2FC **−0.13**
- **not one passes |log2FC| > 1**

So none of them ever enter the DE list, and ORA cannot see the set at all. A coherent,
highly significant, biologically real shift is invisible to a fold-change cut.

## Which to use

- **GSEA** when you want to know what is happening biologically — uses all the data,
  detects coordinated small changes, no threshold dependence.
- **ORA** when you already have a gene list you care about and want to characterise it —
  a cluster, a module, a marker panel, a list from a paper.

For this dataset GSEA is clearly better: the DE lists are small because |log2FC| > 1 is
strict, and the real biology lives in coordinated sub-two-fold changes. ORA does reproduce
the strongest gene-level signals (EMT in the organ contrasts, interferon/inflammatory
response in saccule ageing) — where effects are large, both methods agree.

## Running it

```bash
ENRICH="$(conda info --base)/envs/enrich_r"
ENV1="$(conda info --base)/envs/Env1"

$ENRICH/bin/Rscript scripts/01_run_ora.R      # ~3 min, login node is fine
$ENRICH/bin/Rscript scripts/02_figures.R      # also reads ../gsea for the comparison
$ENV1/bin/python  scripts/03_build_report.py  # 18 slides
```

`02_figures.R` depends on `../gsea/tables/gsea_all_contrasts.tsv.gz`, so run the GSEA first.
