#!/usr/bin/env python3
"""
08_interaction_eqtl_pooled.py -- interaction eQTL on POOLED (w8+w22) data.

Same model/safeguards as 05, but on the pooled inputs built by 07a/07b, and with
a selectable covariate set for the RNA-PC sensitivity analysis:
  --covar core : geno PC1-3 + age            (no RNA PCs)
  --covar full : core + RNA(expr) PC1-8

  python3 08_interaction_eqtl_pooled.py --tissue saccule --covar core
Run on a GPU node (Env1: tensorqtl + torch/cuda).
"""
import argparse, os
import numpy as np, pandas as pd, torch
import tensorqtl
from tensorqtl import genotypeio, cis
from scipy.stats import rankdata, norm

POOL="results/pooled_inputs"
ENR="../08_deconvolution/results/tables"
OUTROOT="results/ieqtl_pooled"
CELLTYPES=["Hair_cell","Supporting_cell","Transitional_epith",
           "Glia_Schwann","Fibrocyte_mesench","Roof_dark_cell"]

def inv_norm(x): return norm.ppf((rankdata(x)-0.5)/len(x))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--tissue", required=True, choices=["saccule","utricle"])
    ap.add_argument("--covar", required=True, choices=["core","full"])
    ap.add_argument("--celltypes", default="all")
    ap.add_argument("--window", type=int, default=1_000_000)
    a=ap.parse_args()
    cts=CELLTYPES if a.celltypes=="all" else a.celltypes.split(",")
    out_dir=os.path.join(OUTROOT, f"{a.tissue}_{a.covar}"); os.makedirs(out_dir, exist_ok=True)
    print("* GPU:", torch.cuda.is_available(),
          torch.cuda.get_device_name(0) if torch.cuda.is_available() else "(CPU!)")

    pheno,pheno_pos=tensorqtl.read_phenotype_bed(f"{POOL}/{a.tissue}_expression.bed.gz")
    covar=pd.read_csv(f"{POOL}/{a.tissue}_covar_{a.covar}.tsv", sep="\t", index_col=0).T
    pr=genotypeio.PlinkReader(f"{POOL}/{a.tissue}_pooled")
    geno=pr.load_genotypes(); variant_df=pr.bim.set_index("snp")[["chrom","pos"]]
    chroms=set(variant_df["chrom"].unique())
    keep=pheno_pos["chr"].isin(chroms); pheno,pheno_pos=pheno.loc[keep],pheno_pos.loc[keep]
    enr=pd.read_csv(f"{ENR}/enrichment_ssgsea_{a.tissue}.csv", index_col=0)

    samples=[s for s in pheno.columns if s in covar.index and s in geno.columns and s in enr.index]
    print(f"* {a.tissue}/{a.covar}: {len(samples)} samples, {pheno.shape[0]} genes, "
          f"{geno.shape[0]} variants, {covar.shape[1]} covariates")
    pheno=pheno[samples]; covar=covar.loc[samples]; geno=geno[samples]

    for ct in cts:
        if ct not in enr.columns: continue
        inter=pd.DataFrame({ct: inv_norm(enr.loc[samples,ct].values)}, index=samples)
        print(f"  [run] {ct}")
        cis.map_nominal(geno, variant_df, pheno, pheno_pos, prefix=ct,
                        covariates_df=covar, interaction_df=inter,
                        maf_threshold_interaction=0.05, window=a.window,
                        run_eigenmt=True, output_dir=out_dir,
                        write_top=True, write_stats=False, verbose=False)
    print("DONE", a.tissue, a.covar)

if __name__=="__main__":
    main()
