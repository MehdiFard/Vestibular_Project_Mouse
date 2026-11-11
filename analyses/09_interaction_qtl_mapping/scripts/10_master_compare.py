#!/usr/bin/env python3
"""
10_master_compare.py -- unify ALL interaction-eQTL approaches + standard eQTL.

Approaches:
  per-group full  : results/ieqtl/<group>/            (Pantry covars incl RNA PCs)
  per-group core  : results/ieqtl_core/<group>/       (geno PCs only)
  pooled core     : results/ieqtl_pooled/<tis>_core/  (geno PCs + age)
  pooled full     : results/ieqtl_pooled/<tis>_full/  (+ RNA PC1-8)
  standard eQTL   : Pantry expression.cis_qtl (avg eGenes/group)

Outputs (results/mapping_report/): master comparison figures + tables.
Run: Env1 python.
"""
import os, glob
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

D="results"; OUT=os.path.join(D,"mapping_report"); os.makedirs(OUT,exist_ok=True)
PANTRY="../../analyses/03_pantry"
CTS=["Hair_cell","Supporting_cell","Transitional_epith","Glia_Schwann","Fibrocyte_mesench","Roof_dark_cell"]

def sig(path):
    d=pd.read_csv(path,sep="\t"); return int((d.pval_adj_bh<0.05).sum())

def counts(pattern):
    """sum ieGenes per cell type over all files matching pattern -> dict ct->n, plus total"""
    c={ct:0 for ct in CTS}
    for f in glob.glob(pattern):
        ct=os.path.basename(f).split(".")[0]
        if ct in c: c[ct]+=sig(f)
    return c

app={
 "per-group full":  counts(f"{D}/ieqtl/*/*.cis_qtl_top_assoc.txt.gz"),
 "per-group core":  counts(f"{D}/ieqtl_core/*/*.cis_qtl_top_assoc.txt.gz"),
 "pooled core":     counts(f"{D}/ieqtl_pooled/*_core/*.cis_qtl_top_assoc.txt.gz"),
 "pooled full":     counts(f"{D}/ieqtl_pooled/*_full/*.cis_qtl_top_assoc.txt.gz"),
}
tot={k:sum(v.values()) for k,v in app.items()}
egene_avg=float(np.mean([(pd.read_csv(f"{PANTRY}/pheast/output/{g}/qtl/expression.cis_qtl.txt.gz",sep="\t").qval<0.05).sum()
                         for g in ["saccule_w8","saccule_w22","utricle_w8","utricle_w22"]]))

# master table
tab=pd.DataFrame(app).T[CTS]; tab["TOTAL"]=tab.sum(axis=1)
tab.to_csv(f"{OUT}/master_ieGene_counts.csv")
print(tab.to_string()); print("avg eGenes/group:",round(egene_avg))

# M1: total ieGenes by approach (+ eQTL reference)
fig,ax=plt.subplots(figsize=(8,4.6))
labels=[k.replace(" ","\n") for k in app]+["standard\neQTL (avg)"]; vals=[tot[k] for k in app]+[egene_avg]
cols=["#8497B0","#5B7DB1","#2E86AB","#E4572E","#5B8C5A"]
ax.bar(labels,vals,color=cols)
for i,v in enumerate(vals): ax.text(i,v,f"{v:.0f}",ha="center",va="bottom",fontsize=10)
ax.set_yscale("log"); ax.set_ylabel("ieGene hits / eGenes (log)")
ax.set_title("Interaction-eQTL discovery across all approaches")
fig.tight_layout(); fig.savefig(f"{OUT}/M1_approach_totals.png",dpi=200); plt.close()

# M2: ieGenes by cell type across the four ieQTL approaches
fig,ax=plt.subplots(figsize=(11,4.8)); x=np.arange(len(CTS)); w=0.2
for i,(k,c) in enumerate(zip(app.keys(),["#8497B0","#5B7DB1","#2E86AB","#E4572E"])):
    vals=[app[k][ct] for ct in CTS]
    ax.bar(x+(i-1.5)*w,vals,w,label=k,color=c)
    for j,v in enumerate(vals):
        if v: ax.text(x[j]+(i-1.5)*w,v,v,ha="center",va="bottom",fontsize=7)
ax.set_xticks(x); ax.set_xticklabels(CTS,rotation=25,ha="right")
ax.set_ylabel("ieGenes (summed over groups/tissues)")
ax.set_title("ieGenes by cell type across approaches")
ax.legend(fontsize=9,ncol=4); fig.tight_layout(); fig.savefig(f"{OUT}/M2_by_celltype.png",dpi=200); plt.close()

print("saved master figures + table ->",OUT)
