#!/usr/bin/env python3
"""
07b_pool_expression_covar.py -- pooled (w8+w22) expression BED + covariate sets.

Expression: concatenate the two ages' unnormalized gene TPM (same genes), then
normalize with Pantry's OWN normalize_phenotypes.py (quantile-normalize ->
inverse-normal, dropping >50%-zero / no-variation genes) so the processing is
identical to the per-group pipeline. bgzip + tabix the result.

Covariates (two sets, for the sensitivity analysis):
  core = geno PC1-3 + age(w8=0/w22=1)                       [no RNA PCs]
  full = core + expression(RNA) PC1-8
Also writes the exprPC x enrichment correlation, to flag composition-correlated PCs.

Run from project root (Env1 python). Genotype PCA (07a) must already be done.
"""
import os, subprocess, sys
import numpy as np, pandas as pd

ENV1=os.path.dirname(sys.executable)   # bin/ of the env running this script (Env1)
PANTRY="../../analyses/03_pantry"
NORM=f"{PANTRY}/phenotyping/scripts/normalize_phenotypes.py"
OUT="results/pooled_inputs"
ENR="../08_deconvolution/results/tables"

def main():
    for tis in ["saccule","utricle"]:
        print(f"=== {tis} ===")
        w8 = pd.read_csv(f"{PANTRY}/phenotyping/output/{tis}_w8/unnorm/expression.bed", sep="\t",
                         dtype={"#chr":str}, low_memory=False)
        w22= pd.read_csv(f"{PANTRY}/phenotyping/output/{tis}_w22/unnorm/expression.bed", sep="\t",
                         dtype={"#chr":str}, low_memory=False)
        meta=list(w8.columns[:4])                       # #chr start end phenotype_id
        s8=list(w8.columns[4:]); s22=list(w22.columns[4:])
        pooled=w8.merge(w22, on=meta)                   # same genes -> inner join on positions+id
        pooled=pooled.sort_values(["#chr","start"])
        pooled.to_csv(f"{OUT}/{tis}_unnorm.bed", sep="\t", index=False)
        samples=s8+s22
        with open(f"{OUT}/{tis}.samples.txt","w") as f: f.write("\n".join(samples)+"\n")
        print(f"  genes={len(pooled)}  samples={len(samples)} (w8={len(s8)}, w22={len(s22)})")

        # Pantry normalization (identical to per-group)
        subprocess.run([f"{ENV1}/python", NORM, "-i", f"{OUT}/{tis}_unnorm.bed",
                        "-s", f"{OUT}/{tis}.samples.txt", "-o", f"{OUT}/{tis}_expression.bed"], check=True)
        subprocess.run([f"{ENV1}/bgzip","-f",f"{OUT}/{tis}_expression.bed"], check=True)
        subprocess.run([f"{ENV1}/tabix","-p","bed",f"{OUT}/{tis}_expression.bed.gz"], check=True)

        # expression (RNA) PCs from the normalized matrix
        expr=pd.read_csv(f"{OUT}/{tis}_expression.bed.gz", sep="\t")
        X=expr[samples].T.values.astype(float)          # samples x genes
        Xc=X-X.mean(0)
        U,S,_=np.linalg.svd(Xc, full_matrices=False)
        pcs=U[:,:10]*S[:10]
        exprPC=pd.DataFrame(pcs, index=samples, columns=[f"exprPC{i+1}" for i in range(10)])

        # genotype PCs (plink2 eigenvec: #FID IID PC1..PC5)
        ev=pd.read_csv(f"{OUT}/{tis}_pca.eigenvec", sep="\t").rename(columns={"IID":"sample"})
        ev=ev.set_index("sample")[[f"PC{i}" for i in range(1,6)]]
        ev.columns=[f"geno_PC{i}" for i in range(1,6)]

        age=pd.Series([0]*len(s8)+[1]*len(s22), index=samples, name="age")

        # sequencing batch (plate 1-4) from Pantry's batch files (pheast/input/<group>/<group>.batches.tsv; the
        # final Pantry run's covariate files have genotype + phenotype PCs only, no batch row). Codes are the same
        # across tissues/ages, so they concatenate directly.
        bat={}
        for grp in (f"{tis}_w8", f"{tis}_w22"):
            b=pd.read_csv(f"{PANTRY}/pheast/input/{grp}/{grp}.batches.tsv", sep="\t", header=None, index_col=0)[1]
            bat.update(b.to_dict())
        batch=pd.Series([bat[s] for s in samples], index=samples, name="batch")

        # core = technical/biological covariates (NO RNA PCs); full = core + RNA PC1-8.
        # Per-group Pantry recipe is 3 geno PCs + 8 expression PCs; age and batch are added because ages are pooled.
        core=pd.concat([ev.loc[samples,["geno_PC1","geno_PC2","geno_PC3"]], age, batch], axis=1)
        full=pd.concat([core, exprPC.loc[samples,[f"exprPC{i}" for i in range(1,9)]]], axis=1)
        for name,cov in [("core",core),("full",full)]:
            t=cov.T.copy(); t.index.name="ID"
            t.to_csv(f"{OUT}/{tis}_covar_{name}.tsv", sep="\t")
        print(f"  covar core={core.shape[1]} full={full.shape[1]} covariates")

        # RNA PC vs cell-type enrichment correlation (flag composition-correlated PCs)
        enr=pd.read_csv(f"{ENR}/enrichment_ssgsea_{tis}.csv", index_col=0).loc[samples]
        cors=pd.DataFrame({ct:[np.corrcoef(exprPC[f"exprPC{i+1}"],enr[ct])[0,1] for i in range(8)]
                           for ct in enr.columns}, index=[f"exprPC{i+1}" for i in range(8)])
        cors.to_csv(f"{OUT}/{tis}_exprPC_vs_enrichment_corr.csv")
        top=cors.abs().max(axis=1)
        print(f"  |corr| exprPC~enrichment (max over cell types): "
              + ", ".join(f"{k}={v:.2f}" for k,v in top.items()))
    print("DONE")

if __name__=="__main__":
    main()
