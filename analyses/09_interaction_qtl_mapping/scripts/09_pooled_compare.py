#!/usr/bin/env python3
"""
09_pooled_compare.py -- pooled ieQTL: covariate sensitivity + cross-approach comparison.

Produces:
  S1  pooled core vs full ieGene counts per cell type (RNA-PC sensitivity)
  S2  |corr| RNA-PC x cell-type enrichment heatmap (which PCs capture composition)
  S3  approach comparison: per-group ieGenes vs pooled(core) vs pooled(full) vs eGenes
  robust_ieGenes.csv  : ieGenes significant in BOTH covariate models (confident set)
Run: Env1 python.
"""
import os, glob
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

POOL="results/ieqtl_pooled"; INP="results/pooled_inputs"
PG="results/ieqtl"; PANTRY="../../analyses/03_pantry"
PLT=os.path.join(POOL,"plots"); os.makedirs(PLT,exist_ok=True)
CTS=["Hair_cell","Supporting_cell","Transitional_epith","Glia_Schwann","Fibrocyte_mesench","Roof_dark_cell"]
TIS=["saccule","utricle"]

def sig_genes(path):
    d=pd.read_csv(path,sep="\t"); return set(d.loc[d.pval_adj_bh<0.05,"phenotype_id"])

# collect pooled sig sets per tissue/covar/celltype
pooled={}   # (tis,cov,ct)->set
for f in glob.glob(f"{POOL}/*/*.cis_qtl_top_assoc.txt.gz"):
    cfg=f.split("/")[-2]; ct=os.path.basename(f).split(".")[0]; tis,cov=cfg.rsplit("_",1)
    pooled[(tis,cov,ct)]=sig_genes(f)
pd.DataFrame([dict(tissue=t,covar=c,cell_type=ct,ieGenes=len(v)) for (t,c,ct),v in sorted(pooled.items())]) \
  .to_csv(f"{POOL}/pooled_ieGene_summary.csv",index=False)

# ---- S1: core vs full per cell type ----------------------------------------
fig,axes=plt.subplots(1,2,figsize=(12,4.6),sharey=True)
for ax,tis in zip(axes,TIS):
    core=[len(pooled.get((tis,"core",ct),set())) for ct in CTS]
    full=[len(pooled.get((tis,"full",ct),set())) for ct in CTS]
    x=np.arange(len(CTS)); w=0.38
    ax.bar(x-w/2,core,w,label="core (no RNA PCs)",color="#2E86AB")
    ax.bar(x+w/2,full,w,label="full (+RNA PC1-8)",color="#E4572E")
    for i,v in enumerate(core): ax.text(i-w/2,v,v,ha="center",va="bottom",fontsize=8)
    for i,v in enumerate(full): ax.text(i+w/2,v,v,ha="center",va="bottom",fontsize=8)
    ax.set_title(tis); ax.set_xticks(x); ax.set_xticklabels(CTS,rotation=40,ha="right",fontsize=8)
axes[0].set_ylabel("ieGenes (FDR<0.05)"); axes[0].legend(fontsize=9)
fig.suptitle("Pooled ieGenes: covariate sensitivity (core vs full)")
fig.tight_layout(); fig.savefig(f"{PLT}/S1_core_vs_full.png",dpi=200); plt.close()

# ---- S2: RNA-PC x enrichment correlation heatmap ---------------------------
fig,axes=plt.subplots(1,2,figsize=(12,4.2))
for ax,tis in zip(axes,TIS):
    c=pd.read_csv(f"{INP}/{tis}_exprPC_vs_enrichment_corr.csv",index_col=0)[CTS]
    im=ax.imshow(c.values.astype(float),cmap="RdBu_r",vmin=-0.8,vmax=0.8,aspect="auto")
    ax.set_xticks(range(len(CTS))); ax.set_xticklabels(CTS,rotation=40,ha="right",fontsize=8)
    ax.set_yticks(range(len(c.index))); ax.set_yticklabels(c.index,fontsize=8)
    for i in range(c.shape[0]):
        for j in range(c.shape[1]):
            ax.text(j,i,f"{c.values[i,j]:.2f}",ha="center",va="center",fontsize=7)
    ax.set_title(f"{tis}: corr(RNA PC, enrichment)")
fig.colorbar(im,ax=axes,shrink=0.7,label="Pearson r")
fig.savefig(f"{PLT}/S2_rnaPC_enrichment_heatmap.png",dpi=200,bbox_inches="tight"); plt.close()

# ---- S3: approach comparison (total ieGenes) -------------------------------
# per-group union of ieGenes
pg=set()
for f in glob.glob(f"{PG}/*/*.cis_qtl_top_assoc.txt.gz"): pg|= {(f.split('/')[-2],g) for g in sig_genes(f)}
pg_total=sum(1 for f in glob.glob(f"{PG}/*/*.cis_qtl_top_assoc.txt.gz") for _ in sig_genes(f))
core_total=sum(len(v) for k,v in pooled.items() if k[1]=="core")
full_total=sum(len(v) for k,v in pooled.items() if k[1]=="full")
egene_avg=np.mean([ (pd.read_csv(f"{PANTRY}/pheast/output/{g}/qtl/expression.cis_qtl.txt.gz",sep="\t").qval<0.05).sum()
                    for g in ["saccule_w8","saccule_w22","utricle_w8","utricle_w22"]])
fig,ax=plt.subplots(figsize=(7.5,4.5))
labels=["per-group\nieQTL","pooled core\nieQTL","pooled full\nieQTL","standard eQTL\n(avg/group)"]
vals=[pg_total,core_total,full_total,egene_avg]
cols=["#8497B0","#2E86AB","#E4572E","#5B8C5A"]
ax.bar(labels,vals,color=cols)
for i,v in enumerate(vals): ax.text(i,v,f"{v:.0f}",ha="center",va="bottom",fontsize=10)
ax.set_yscale("log"); ax.set_ylabel("ieGene hits / eGenes (log)")
ax.set_title("Interaction eQTL discovery across approaches")
fig.tight_layout(); fig.savefig(f"{PLT}/S3_approach_comparison.png",dpi=200); plt.close()

# ---- robust ieGenes (sig in BOTH core and full) ----------------------------
rows=[]
for tis in TIS:
    for ct in CTS:
        both=pooled.get((tis,"core",ct),set()) & pooled.get((tis,"full",ct),set())
        for g in both: rows.append(dict(tissue=tis,cell_type=ct,phenotype_id=g))
rob=pd.DataFrame(rows); rob.to_csv(f"{POOL}/robust_ieGenes.csv",index=False)

print(f"per-group ieGene hits: {pg_total}")
print(f"pooled core: {core_total} | pooled full: {full_total} | robust(core∩full): {len(rob)}")
print(f"avg eGenes/group: {egene_avg:.0f}")
print("saved plots + robust_ieGenes.csv ->",POOL)
