#!/usr/bin/env python3
"""
06_ieqtl_vs_eqtl_compare.py  --  compare cell-type interaction eQTLs (this work)
with the standard cis-eQTLs mapped earlier in Pantry.

Standard eQTL (Pantry, expression.cis_qtl) finds marginal effects (thousands of
eGenes); the interaction test isolates the cell-type-specific subset. We show:
  A. eGene vs ieGene counts per group (log scale)
  B. ieGene counts heatmap (group x cell type)
  C. how many ieGenes are also standard eGenes vs interaction-only
  D. for ieGenes, marginal genotype p (pval_g) vs interaction p (pval_gi)

Outputs: results/ieqtl/plots/*.png + comparison table.
Run: Env1 python (pandas, matplotlib).
"""
import os, glob
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

PANTRY="../../analyses/03_pantry"
IEQTL="results/ieqtl"
PLT=os.path.join(IEQTL,"plots"); os.makedirs(PLT,exist_ok=True)
GROUPS=["saccule_w8","saccule_w22","utricle_w8","utricle_w22"]
CTS=["Hair_cell","Supporting_cell","Transitional_epith","Glia_Schwann","Fibrocyte_mesench","Roof_dark_cell"]

# ---- standard eGenes per group (qval<0.05) ---------------------------------
egenes={}
for g in GROUPS:
    d=pd.read_csv(f"{PANTRY}/pheast/output/{g}/qtl/expression.cis_qtl.txt.gz",sep="\t")
    egenes[g]=set(d.loc[d["qval"]<0.05,"phenotype_id"])

# ---- ieGenes per group/cell type -------------------------------------------
ie_rows=[]           # sig ieGene rows (with pval_g, pval_gi)
count=pd.DataFrame(0,index=CTS,columns=GROUPS)
ie_by_group={g:set() for g in GROUPS}
for g in GROUPS:
    for ct in CTS:
        f=f"{IEQTL}/{g}/{ct}.cis_qtl_top_assoc.txt.gz"
        if not os.path.exists(f): continue
        d=pd.read_csv(f,sep="\t")
        s=d[d["pval_adj_bh"]<0.05]
        count.loc[ct,g]=len(s)
        ie_by_group[g]|=set(s["phenotype_id"])
        for _,r in s.iterrows():
            ie_rows.append(dict(group=g,cell_type=ct,phenotype_id=r["phenotype_id"],variant_id=r["variant_id"],
                                pval_g=r["pval_g"],pval_gi=r["pval_gi"],pval_adj_bh=r["pval_adj_bh"],
                                is_eGene=r["phenotype_id"] in egenes[g]))
ie=pd.DataFrame(ie_rows)
ie.drop(columns=["variant_id","pval_adj_bh"]).to_csv(f"{IEQTL}/ieqtl_vs_eqtl_table.csv",index=False)
# full list of significant ieGenes with gene symbols (named in the reports)
sym=pd.read_csv(f"{PANTRY}/../05_dge/data/gene_annotation.tsv",sep="\t").set_index("gene_id")["gene_name"]
ie.assign(symbol=ie["phenotype_id"].map(sym))[["group","cell_type","phenotype_id","variant_id","pval_gi","pval_adj_bh","symbol"]] \
  .to_csv(f"{IEQTL}/ieGenes_significant.csv",index=False)

# ---- A. counts eGenes vs ieGenes per group ---------------------------------
eg=[len(egenes[g]) for g in GROUPS]; ig=[len(ie_by_group[g]) for g in GROUPS]
x=np.arange(len(GROUPS)); w=0.38
fig,ax=plt.subplots(figsize=(8,4.5))
ax.bar(x-w/2,eg,w,label="standard eGenes (cis-eQTL)",color="#2E86AB")
ax.bar(x+w/2,[max(v,0.3) for v in ig],w,label="ieGenes (interaction, any cell type)",color="#E4572E")
ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels(GROUPS,rotation=15)
for i,v in enumerate(ig): ax.text(i+w/2,max(v,0.3),str(v),ha="center",va="bottom",fontsize=9)
for i,v in enumerate(eg): ax.text(i-w/2,v,str(v),ha="center",va="bottom",fontsize=9)
ax.set_ylabel("number of genes (log scale)")
ax.set_title("Standard cis-eQTL vs cell-type interaction eQTL (per group)")
ax.legend(fontsize=9); fig.tight_layout(); fig.savefig(f"{PLT}/A_egene_vs_iegene_counts.png",dpi=200); plt.close()

# ---- B. ieGene heatmap group x cell type -----------------------------------
fig,ax=plt.subplots(figsize=(7,4.5))
im=ax.imshow(count.values,cmap="OrRd",aspect="auto")
ax.set_xticks(range(len(GROUPS))); ax.set_xticklabels(GROUPS,rotation=20)
ax.set_yticks(range(len(CTS))); ax.set_yticklabels(CTS)
for i in range(len(CTS)):
    for j in range(len(GROUPS)):
        ax.text(j,i,int(count.values[i,j]),ha="center",va="center",
                color="black" if count.values[i,j]<count.values.max()*0.6 else "white",fontsize=10)
ax.set_title("ieGenes (FDR<0.05) by group and cell type")
fig.colorbar(im,ax=ax,shrink=0.8,label="ieGenes"); fig.tight_layout()
fig.savefig(f"{PLT}/B_iegene_heatmap.png",dpi=200); plt.close()

# ---- C. ieGenes: also a standard eGene, or interaction-only? ----------------
if len(ie):
    shared=int(ie["is_eGene"].sum()); only=int((~ie["is_eGene"]).sum())
    fig,ax=plt.subplots(figsize=(5.5,4.5))
    ax.bar(["also a standard\neGene","interaction-only\n(not an eGene)"],[shared,only],
           color=["#2E86AB","#E4572E"])
    for i,v in enumerate([shared,only]): ax.text(i,v,str(v),ha="center",va="bottom",fontsize=11)
    ax.set_ylabel("ieGene hits"); ax.set_title("Are the ieGenes captured by standard eQTL?")
    fig.tight_layout(); fig.savefig(f"{PLT}/C_iegene_overlap.png",dpi=200); plt.close()

    # ---- D. marginal vs interaction p ------------------------------------
    fig,ax=plt.subplots(figsize=(5.8,5))
    col=ie["is_eGene"].map({True:"#2E86AB",False:"#E4572E"})
    ax.scatter(-np.log10(ie["pval_g"]),-np.log10(ie["pval_gi"]),c=col,s=45,edgecolor="k",lw=0.3)
    ax.axhline(-np.log10(0.05),ls=":",c="grey")
    ax.set_xlabel("-log10 marginal genotype p (pval_g)")
    ax.set_ylabel("-log10 interaction p (pval_gi)")
    ax.set_title("ieGenes: interaction vs marginal signal")
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([0],[0],marker='o',ls='',mfc="#2E86AB",mec='k',label='also eGene'),
                       Line2D([0],[0],marker='o',ls='',mfc="#E4572E",mec='k',label='interaction-only')],fontsize=9)
    fig.tight_layout(); fig.savefig(f"{PLT}/D_marginal_vs_interaction.png",dpi=200); plt.close()

    print(f"ieGenes total hits: {len(ie)} | also standard eGene: {shared} | interaction-only: {only}")
print("eGene counts:",{g:len(egenes[g]) for g in GROUPS})
print("ieGene (union/group):",{g:len(ie_by_group[g]) for g in GROUPS})
print("saved plots ->",PLT)
