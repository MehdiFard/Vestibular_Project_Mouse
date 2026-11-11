#!/usr/bin/env python
"""Illustration + ANNOTATED formula images for the deconvolution intro slides.
Every formula term is labelled, and 'signature' / 'fraction' are defined visually."""
import os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch, Rectangle

OUT="results/intro_assets"; os.makedirs(OUT,exist_ok=True)
INK="#1a1a2e"; ACC="#2e86ab"; GREY="#555b6e"; RED="#E4572E"; GRN="#1b7a3d"; PUR="#6a4c93"
PAL=["#4363D8","#469990","#BFEF45","#3CB44B","#E6194B","#000075","#F58231"]
plt.rcParams.update({"mathtext.fontset":"cm","font.family":"DejaVu Sans"})

def callout(ax, term_x, term_y, label, color, lx, ly):
    ax.annotate(label, xy=(term_x, term_y), xytext=(lx, ly), ha="center", va="top",
                fontsize=10, color=color,
                arrowprops=dict(arrowstyle="-", color=color, lw=1.1),
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=color, lw=1.1))

# ================= intro 1: mixture model, annotated =================
fig,ax=plt.subplots(figsize=(11,5.6)); ax.axis("off"); ax.set_xlim(0,11); ax.set_ylim(0,10)
np.random.seed(3)
cx,cy=1.9,8.2
for i in range(16):
    a=np.random.rand()*2*np.pi; r=np.random.rand()*1.3
    ax.add_patch(Circle((cx+r*np.cos(a), cy+r*np.sin(a)),0.24,color=PAL[i%len(PAL)],ec="white",lw=0.7))
ax.text(cx,cy+1.7,"cell types in the tissue",ha="center",fontsize=11,color=INK,weight="bold")
ax.add_patch(FancyArrowPatch((3.5,8.2),(4.6,8.2),arrowstyle="-|>",mutation_scale=18,color=GREY,lw=2))
props=[0.30,0.24,0.16,0.10,0.09,0.06,0.05]; bot=6.7; x0,w=5.0,1.1
for p,c in zip(props,PAL): ax.add_patch(Rectangle((x0,bot),w,p*3.0,color=c,ec="white",lw=0.6)); bot+=p*3.0
ax.text(x0+w/2,6.35,"one bulk sample",ha="center",fontsize=11,color=INK,weight="bold")
ax.text(8.4,8.2,"a bulk sample = the sum of\nits cell types' expression,\nweighted by how many of\neach are present",
        ha="center",va="center",fontsize=10.5,color=GREY)

# equation (coloured terms) at y=4.5
yq=4.6
ax.text(2.7,yq,"=",fontsize=24,ha="center",color=INK); ax.text(7.05,yq,"+",fontsize=22,ha="center",color=INK)
ax.text(1.55,yq,r"$B_{g,s}$",ha="center",va="center",fontsize=26,color=ACC)
ax.text(3.5,yq,r"$\sum_{c}$",ha="center",va="center",fontsize=26,color=GREY)
ax.text(4.55,yq,r"$w_{c,s}$",ha="center",va="center",fontsize=26,color=GRN)
ax.text(5.95,yq,r"$S_{g,c}$",ha="center",va="center",fontsize=26,color=PUR)
ax.text(7.75,yq,r"$\varepsilon_{g,s}$",ha="center",va="center",fontsize=24,color="#8a8a8a")
callout(ax,1.55,yq-0.55,"bulk expression\n(gene g, sample s)\nwhat we measure",ACC,1.5,3.4)
callout(ax,4.55,yq-0.55,"cell-type FRACTION\n(cell type c in sample s)\nwhat we estimate",GRN,4.35,3.4)
callout(ax,5.95,yq-0.55,"cell-type SIGNATURE\n(gene g in cell type c)\nfrom single-cell",PUR,7.2,3.4)
callout(ax,7.75,yq-0.5,"noise","#8a8a8a",9.4,3.6)
ax.text(3.4,0.35,r"sum over cell types $c$",ha="center",fontsize=9.5,color=GREY)
fig.tight_layout(); fig.savefig(f"{OUT}/intro1_mixture.png",dpi=200,bbox_inches="tight"); plt.close()

# ================= intro 2: signatures & fractions =================
fig=plt.figure(figsize=(11,5.6)); ax=fig.add_axes([0,0,1,1]); ax.axis("off"); ax.set_xlim(0,12); ax.set_ylim(0,10)
def mat(x,y,w,h,label,dim,color):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.02,rounding_size=0.08",fc=color,ec=INK,lw=1.3,alpha=0.85))
    ax.text(x+w/2,y+h/2,label,ha="center",va="center",fontsize=20,color="white",weight="bold")
    ax.text(x+w/2,y-0.45,dim,ha="center",fontsize=9.5,color=GREY)
ax.text(6.0,9.4,"Deconvolution splits bulk into SIGNATURES x FRACTIONS",ha="center",fontsize=13,color=INK,weight="bold")
mat(0.5,4.6,2.3,3.4,"B",r"genes $\times$ samples",ACC)
ax.text(3.15,6.3,r"$\approx$",fontsize=24,color=INK,ha="center")
mat(3.6,4.6,1.9,3.4,"S",r"genes $\times$ cell types",PUR)
ax.text(5.85,6.3,r"$\times$",fontsize=22,color=INK,ha="center")
mat(6.3,5.9,4.9,2.1,"W",r"cell types $\times$ samples",GRN)
# annotate S = signatures, W = fractions
ax.annotate("cell-type SIGNATURES\neach column = one cell type's\ntypical expression profile",
            xy=(4.55,4.6), xytext=(4.55,3.1), ha="center", va="top", fontsize=10, color=PUR,
            arrowprops=dict(arrowstyle="-",color=PUR,lw=1.1),
            bbox=dict(boxstyle="round,pad=0.3",fc="white",ec=PUR,lw=1.1))
ax.annotate("cell-type FRACTIONS\neach column = one sample's\nproportions (sum to 1)",
            xy=(8.7,5.9), xytext=(7.9,3.9), ha="center", va="top", fontsize=10, color=GRN,
            arrowprops=dict(arrowstyle="-",color=GRN,lw=1.1),
            bbox=dict(boxstyle="round,pad=0.3",fc="white",ec=GRN,lw=1.1))
# objective (NNLS) at bottom-centre, clear of both callouts
ax.text(6.0,0.75,r"$\widehat{W}=\mathrm{argmin}_{\,W\geq 0}\;\| B - S\,W \|^{2}$",
        ha="center",va="center",fontsize=15,color=INK)
# mini insets: a signature (gene profile) and a fraction (stacked bar)
sig=fig.add_axes([0.055,0.06,0.15,0.19]); sig.bar(range(8),[3,1,5,2,4,1,2,3],color=PUR); sig.set_title("one signature",fontsize=8,color=PUR)
sig.set_xticks([]); sig.set_yticks([]); sig.set_xlabel("genes",fontsize=7,color=GREY)
frac=fig.add_axes([0.905,0.07,0.045,0.20]); b=0
for p,c in zip([0.34,0.26,0.18,0.12,0.10],PAL): frac.bar(0,p,bottom=b,color=c,width=0.8); b+=p
frac.set_title("one\nsample",fontsize=8,color=GRN); frac.set_xticks([]); frac.set_yticks([]); frac.set_xlim(-0.6,0.6)
fig.savefig(f"{OUT}/intro2_deconv.png",dpi=200,bbox_inches="tight"); plt.close()

# ================= intro 3: interaction eQTL, annotated =================
fig,ax=plt.subplots(figsize=(9,5.2))
G=np.array([0,1,2])
ax.plot(G,2.0+0.25*G,"-o",color=ACC,lw=2.5,ms=8,label="low cell-type abundance")
ax.plot(G,2.0+1.35*G,"-o",color=RED,lw=2.5,ms=8,label="high cell-type abundance")
ax.set_xticks([0,1,2]); ax.set_xticklabels(["AA","Aa","aa"])
ax.set_xlabel("genotype",fontsize=12); ax.set_ylabel("expression",fontsize=12)
ax.set_title("The genotype effect depends on cell-type abundance = a cell-type interaction eQTL",
             fontsize=12,color=INK,weight="bold")
ax.legend(fontsize=10,loc="upper left"); ax.grid(alpha=0.25)
# annotated formula below
ax.text(0.5,-0.24,
        r"$y = \beta_0 + \beta_g\,G + \beta_w\,w_c + \boldsymbol{\beta_{gw}\,(G\times w_c)} + C\gamma$",
        transform=ax.transAxes, ha="center", fontsize=16, color=INK)
for tx,txt,col in [(0.30,"genotype\neffect",ACC),(0.44,"abundance\neffect",GRN),
                   (0.62,"INTERACTION =\ncell-type-specific effect",RED),(0.80,"other\ncovariates",GREY)]:
    ax.text(tx,-0.46,txt,transform=ax.transAxes,ha="center",va="top",fontsize=8.5,color=col)
ax.text(0.5,-0.60,r"$y$ = expression,  $G$ = genotype dosage,  $w_c$ = cell-type fraction/enrichment",
        transform=ax.transAxes,ha="center",fontsize=9.5,color=GREY)
fig.subplots_adjust(bottom=0.40)
fig.savefig(f"{OUT}/intro3_interaction.png",dpi=200,bbox_inches="tight"); plt.close()

print("saved annotated intro assets ->",OUT)
