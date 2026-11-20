# Over-representation analysis (ORA)

The classical pathway test, run on the DE gene lists — the counterpart to [`../gsea`](../gsea).
Included for comparison with the GSEA (the primary analysis).

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

## Which to use

- **GSEA** when you want to know what is happening biologically — uses all the data,
  detects coordinated small changes, no threshold dependence.
- **ORA** when you already have a gene list you care about and want to characterise it —
  a cluster, a module, a marker panel, a list from a paper.

ORA only sees genes that pass both the FDR and the |log2FC| > 1 cut, so a set whose genes all shift consistently but
by less than two-fold is invisible to it; GSEA, which uses the whole ranking, detects such coordinated changes.
`scripts/02_figures.R` draws the comparison of the two methods.

## Running it

```bash
ENRICH="$(conda info --base)/envs/enrich_r"

$ENRICH/bin/Rscript scripts/01_run_ora.R      # a few minutes
$ENRICH/bin/Rscript scripts/02_figures.R      # also reads ../gsea for the comparison
```

`02_figures.R` depends on `../gsea/tables/gsea_all_contrasts.tsv.gz`, so run the GSEA first.
