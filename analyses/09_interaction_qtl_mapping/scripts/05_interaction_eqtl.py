#!/usr/bin/env python3
"""
05_interaction_eqtl.py  --  cell-type interaction eQTL mapping (GTEx / Kim-Hellmuth style).

For a given Pantry group (e.g. saccule_w8) and each cell type, tests the
genotype x cell-type-enrichment interaction on gene expression with tensorQTL:

    expression ~ genotype + covariates + enrichment + genotype:enrichment

The interaction covariate is the ssGSEA marker-enrichment score (deconvolution
step), inverse-normal transformed across the group's samples. tensorQTL applies
the GTEx safeguard maf_threshold_interaction=0.05 (MAF >=5% in BOTH halves of the
interaction distribution) and eigenMT + Benjamini-Hochberg (run_eigenmt=True),
so <prefix>.cis_qtl_top_assoc.txt.gz gives one interaction test per gene with
eigenMT p (pval_emt) and BH-adjusted p (pval_adj_bh) for ieGene discovery.

Run on a GPU node (see run_ieqtl_gpu.sbatch). Env: Env1 (tensorqtl 1.0.10, torch+cuda).

  python3 05_interaction_eqtl.py --group saccule_w8 --celltypes all
"""
import argparse, os
import numpy as np, pandas as pd, torch
import tensorqtl
from tensorqtl import genotypeio, cis
from scipy.stats import rankdata, norm

PANTRY = "../../analyses/03_pantry"
ENRICH = "../08_deconvolution/results/tables"          # enrichment_ssgsea_<tissue>.csv
OUTROOT = "results/ieqtl"
CELLTYPES = ["Hair_cell","Supporting_cell","Transitional_epith",
             "Glia_Schwann","Fibrocyte_mesench","Roof_dark_cell"]

def inverse_normal(x):
    """rank-based inverse-normal transform of a 1-D vector"""
    return norm.ppf((rankdata(x) - 0.5) / len(x))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", required=True, help="e.g. saccule_w8")
    ap.add_argument("--celltypes", default="all", help="'all' or comma-separated")
    ap.add_argument("--covar", default="full", choices=["full", "core"],
                    help="full = Pantry covariates (geno PCs + RNA PCs); core = geno PCs only")
    ap.add_argument("--outroot", default=OUTROOT)
    ap.add_argument("--window", type=int, default=1_000_000)
    args = ap.parse_args()
    tissue = args.group.split("_")[0]
    cts = CELLTYPES if args.celltypes == "all" else args.celltypes.split(",")
    out_dir = os.path.join(args.outroot, args.group)
    os.makedirs(out_dir, exist_ok=True)

    print("* GPU:", torch.cuda.is_available(),
          torch.cuda.get_device_name(0) if torch.cuda.is_available() else "(CPU!)")

    # ---- load shared inputs once -------------------------------------------
    expr_bed = f"{PANTRY}/phenotyping/output/{args.group}/expression.bed.gz"
    covar_f  = f"{PANTRY}/pheast/intermediate/{args.group}/covar/expression.covar.tsv"
    geno_pfx = f"{PANTRY}/pheast/input/{args.group}/{args.group}"

    pheno, pheno_pos = tensorqtl.read_phenotype_bed(expr_bed)
    covar = pd.read_csv(covar_f, sep="\t", index_col=0).T              # samples x covars
    if args.covar == "core":
        # keep the technical/biological covariates (genotype PCs + sequencing batch)
        # and drop only the RNA/expression PCs -- so core vs full isolates the RNA PCs.
        # (no age term: each per-group analysis is a single age)
        covar = covar[[c for c in covar.columns
                       if c.startswith("geno_PC") or c == "batch"]]
    pr = genotypeio.PlinkReader(geno_pfx)
    geno = pr.load_genotypes()                                        # variants x samples
    variant_df = pr.bim.set_index("snp")[["chrom", "pos"]]

    # drop phenotypes on chroms absent from genotypes (matches Pantry)
    chroms = set(variant_df["chrom"].unique())
    keep = pheno_pos["chr"].isin(chroms)
    pheno, pheno_pos = pheno.loc[keep], pheno_pos.loc[keep]

    enr = pd.read_csv(os.path.join(ENRICH, f"enrichment_ssgsea_{tissue}.csv"), index_col=0)

    # common samples, consistent order across every input
    samples = [s for s in pheno.columns if s in covar.index and s in geno.columns and s in enr.index]
    print(f"* {args.group}: {len(samples)} samples, {pheno.shape[0]} genes, {geno.shape[0]} variants")
    pheno = pheno[samples]; covar = covar.loc[samples]; geno = geno[samples]

    # ---- one interaction scan per cell type --------------------------------
    for ct in cts:
        if ct not in enr.columns:
            print(f"  [skip] {ct}: no enrichment column"); continue
        inter = pd.DataFrame({ct: inverse_normal(enr.loc[samples, ct].values)}, index=samples)
        print(f"  [run] {ct}: interaction eQTL")
        cis.map_nominal(geno, variant_df, pheno, pheno_pos, prefix=ct,
                        covariates_df=covar, interaction_df=inter,
                        maf_threshold_interaction=0.05, window=args.window,
                        run_eigenmt=True, output_dir=out_dir,
                        write_top=True, write_stats=False, verbose=False)
    print("DONE", args.group)

if __name__ == "__main__":
    main()
