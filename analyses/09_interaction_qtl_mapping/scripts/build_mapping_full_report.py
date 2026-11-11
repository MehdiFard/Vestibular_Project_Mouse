#!/usr/bin/env python
"""FULL report of the mapping step: deck (pptx) + detailed Word (docx). Env1 python.
Data-driven: reads master_ieGene_counts.csv + robust_ieGenes.csv for the numbers."""
import os
import pandas as pd
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from docx import Document
from docx.shared import Pt as DPt, RGBColor as DRGB, Inches as DIn
from docx.enum.text import WD_ALIGN_PARAGRAPH as AL
from docx.enum.table import WD_TABLE_ALIGNMENT

D="../../analyses/09_interaction_qtl_mapping/results"
IEP=os.path.join(D,"ieqtl","plots"); SP=os.path.join(D,"ieqtl_pooled","plots")
MR=os.path.join(D,"mapping_report"); os.makedirs(MR,exist_ok=True)
INK=RGBColor(0x1a,0x1a,0x2e); ACC=RGBColor(0x2e,0x86,0xab); GREY=RGBColor(0x55,0x5b,0x6e)
RED=RGBColor(0xb0,0x2a,0x2a); GRN=RGBColor(0x1b,0x7a,0x3d)

tab=pd.read_csv(os.path.join(MR,"master_ieGene_counts.csv"),index_col=0)
def tot(a): return int(tab.loc[a,"TOTAL"])
pg_full,pg_core=tot("per-group full"),tot("per-group core")
pl_core,pl_full=tot("pooled core"),tot("pooled full")
robust=pd.read_csv(os.path.join(D,"ieqtl_pooled","robust_ieGenes.csv"))
nrob=len(robust)

# ============ DECK ============
prs=Presentation(); prs.slide_width=Inches(13.333); prs.slide_height=Inches(7.5)
BLANK=prs.slide_layouts[6]; SW=prs.slide_width
def slide(): return prs.slides.add_slide(BLANK)
def box(s,l,t,w,h):
    tb=s.shapes.add_textbox(Inches(l),Inches(t),Inches(w),Inches(h)); tb.text_frame.word_wrap=True; return tb.text_frame
def setp(p,t,sz=18,b=False,c=INK,it=False,after=6):
    p.text=t
    for r in p.runs: r.font.size=Pt(sz); r.font.bold=b; r.font.color.rgb=c; r.font.italic=it
    p.space_after=Pt(after); return p
def title(s,t,sub=None):
    tf=box(s,0.6,0.35,12.1,1.1); setp(tf.paragraphs[0],t,29,True,INK)
    bar=s.shapes.add_shape(1,Inches(0.62),Inches(1.12),Inches(2.4),Pt(3)); bar.fill.solid(); bar.fill.fore_color.rgb=ACC; bar.line.fill.background()
    if sub: setp(tf.add_paragraph(),sub,14,False,GREY,it=True)
def bull(s,items,l=0.7,t=1.5,w=6.0,sz=15):
    tf=box(s,l,t,w,5.6)
    for i,(txt,lvl,col,bd) in enumerate(items):
        p=tf.paragraphs[0] if i==0 else tf.add_paragraph()
        setp(p,("• "+txt) if lvl==0 else ("   – "+txt),sz-(1 if lvl else 0),bd,col,after=6)
def pic(s,path,l,t,w=None,h=None):
    if os.path.exists(path): s.shapes.add_picture(path,Inches(l),Inches(t),Inches(w) if w else None,Inches(h) if h else None)
def foot(s,n): setp(box(s,11.6,7.05,1.6,0.4).paragraphs[0],str(n),10,False,GREY,after=0)

# 1 title
s=slide()
r=s.shapes.add_shape(1,0,0,SW,Inches(2.8)); r.fill.solid(); r.fill.fore_color.rgb=INK; r.line.fill.background()
setp(box(s,0.8,0.85,11.7,1.5).paragraphs[0],"Cell-type interaction-eQTL mapping",34,True,RGBColor(255,255,255))
setp(box(s,0.8,2.0,11.7,0.7).paragraphs[0],"Full report — vestibular maculae (mouse), GTEx / Kim-Hellmuth approach on GPU",16,False,RGBColor(0xC8,0xE6,0xF0),it=True)
tf=box(s,0.8,3.25,11.7,3.3)
for i,(t,c) in enumerate([
 ("Interaction: expression ~ genotype + covariates + cell-type score + genotype×score.",INK),
 ("Cell-type score = ssGSEA marker enrichment (INT); tensorQTL on A100; MAF≥0.05 both halves; eigenMT+BH.",INK),
 (f"Four analyses: per-group (full/core) and pooled (core/full).  {nrob} robust ieGenes.",GREY)]):
    setp((tf.paragraphs[0] if i==0 else tf.add_paragraph()),"• "+t,16,False,c,after=10)

# 2 method / design matrix
s=slide(); title(s,"Design","Same model, four covariate/sample designs")
bull(s,[
 ("Inputs per group: Pantry genotypes, quantile-normalized expression, covariates.",0,INK,False),
 ("Interaction term: ssGSEA cell-type enrichment, inverse-normal transformed.",0,INK,False),
 ("Safeguards: MAF ≥ 0.05 in both halves of the enrichment; eigenMT + BH; FDR<0.05.",0,INK,False),
 ("Analyses compared:",0,ACC,True),
 ("per-group full — single age, Pantry covariates (incl. 8 RNA PCs),",1,INK,False),
 ("per-group core — single age, geno PCs + batch (no RNA PCs),",1,INK,False),
 ("pooled core — ages pooled (n≈130), geno PCs + age + batch,",1,INK,False),
 ("pooled full — pooled + RNA PC1-8 (same recipe as per-group full).",1,INK,False),
 ("6 cell types each; ssGSEA is the robust cell-type signal (fractions were unreliable).",0,GREY,False),
],w=12.3)

# 3 standard eQTL vs interaction
s=slide(); title(s,"Standard eQTL vs. interaction eQTL","Different questions; interaction isolates a cell-type-specific subset")
pic(s,os.path.join(IEP,"A_egene_vs_iegene_counts.png"),6.4,1.6,w=6.7)
bull(s,[
 ("Standard cis-eQTL: thousands of eGenes/group (marginal, cell-type-averaged).",0,INK,False),
 ("Interaction eQTL: a focused set whose effect scales with cell-type abundance.",0,INK,False),
 ("37/42 per-group ieGenes are also eGenes; 5 are interaction-only (novel).",0,GRN,True),
 ("Interaction answers: in which cell type does the effect act?",0,ACC,True),
],l=0.55,w=5.7,t=1.8,sz=14)

# 4 per-group results
s=slide(); title(s,"Per-group results","Hair & supporting cells dominate; discovery tracks n")
pic(s,os.path.join(IEP,"B_iegene_heatmap.png"),6.6,1.6,w=6.5)
bull(s,[
 ("Full covariate model, per age-group: 42 ieGenes.",0,INK,False),
 ("Hair_cell strongest in BOTH organs (saccule_w8=10, utricle_w8=8).",0,GRN,True),
 ("Week-22 groups → ~0 (fewer samples, less power).",0,RED,True),
 ("Interaction-only vs marginal signal shown in overlap panels (C, D).",0,GREY,False),
],l=0.55,w=5.8,t=1.8,sz=14)

# 5 discovery across approaches
s=slide(); title(s,"Discovery across approaches","Pooling and covariate choice both matter")
pic(s,os.path.join(MR,"M1_approach_totals.png"),6.3,1.6,w=6.8)
bull(s,[
 (f"per-group full: {pg_full}    per-group core: {pg_core}",0,INK,False),
 (f"pooled core: {pl_core}    pooled full: {pl_full}",0,GRN,True),
 ("Pooling ages (n≈130) roughly doubles discovery vs per-group.",0,INK,False),
 ("Standard eQTL (~6,600/group) shown for scale — a different test.",0,GREY,False),
],l=0.55,w=5.6,t=1.8,sz=14)

# 6 by cell type across approaches
s=slide(); title(s,"By cell type, across approaches","Hair & supporting cells are consistently the strongest")
pic(s,os.path.join(MR,"M2_by_celltype.png"),1.1,1.7,w=11.1)
bull(s,[
 ("Hair and supporting cells carry signal in every design → the robust, best-powered cell types.",0,GRN,True),
 ("Rare types appear only sporadically and are covariate-sensitive.",0,INK,False),
],l=0.7,w=12.0,t=6.2,sz=13)

# 7 covariate sensitivity
s=slide(); title(s,"Covariate sensitivity (RNA PCs ↔ composition)","Over-adjustment check")
pic(s,os.path.join(SP,"S1_core_vs_full.png"),0.5,1.6,w=7.9)
pic(s,os.path.join(SP,"S2_rnaPC_enrichment_heatmap.png"),8.3,1.7,w=4.8)
bull(s,[
 ("RNA PCs do capture composition (saccule PC1↔Supporting r=-0.75; utricle PC2↔Hair +0.66).",0,INK,False),
 ("Net: full>core (denoising gain); hair & supporting GAIN → robust, not under-adjustment.",0,GRN,True),
 ("Low-count types still shift (saccule Transitional 3→1, utricle Glia 1→0) → covariate-sensitive, lower confidence.",0,RED,True),
],l=0.6,w=12.2,t=5.9,sz=12)

# 8 robust
s=slide(); title(s,f"Robust ieGenes ({nrob})","Significant under BOTH pooled covariate models")
bull(s,[
 ("Confident set (robust_ieGenes.csv), dominated by hair & supporting cells.",0,GREY,False),
 ("Hair (saccule): Cadps, Acpp, Gm15816 … ;  (utricle): Kcnf1, Gm15816.",0,INK,False),
 ("Supporting (saccule): Cadps, Filip1, Golm1 … ;  Roof/dark (utricle): Slc6a4, Add2.",0,INK,False),
 ("Cadps recurs in hair + supporting — Ca2+-dependent secretion, plausible in sensory epithelium.",0,INK,False),
 ("Caveat: pooled ≠ per-group (age covariate + re-normalization) → limited exact-gene overlap;",0,RED,False),
 ("treat cell-type assignments, not individual variants, as the durable result.",1,RED,False),
],w=12.4)

# 9 summary
s=slide(); title(s,"Summary & next steps")
bull(s,[
 ("A working GPU interaction-eQTL pipeline recovers cell-type-specific eQTLs from bulk vestibular RNA-seq.",0,INK,False),
 ("Hair and supporting cells give robust signal across every design; rarer types are underpowered.",0,INK,False),
 ("Covariate choice matters: report the robust (covariate-insensitive) set; flag covariate-sensitive hits.",0,INK,False),
 ("Pooling ages boosts power but changes the model — complementary to per-group.",0,INK,False),
 ("NEXT: colocalize robust ieGenes with VsEP / beam-walk GWAS to assign balance loci to cell types.",0,ACC,True),
],w=12.4)

for i,s in enumerate(prs.slides,1): foot(s,i)
DECK=os.path.join(MR,"Mapping_Full_Report.pptx"); prs.save(DECK)

# ============ WORD ============
doc=Document(); st=doc.styles["Normal"]; st.font.name="Calibri"; st.font.size=DPt(11)
st.paragraph_format.space_after=DPt(6); st.paragraph_format.line_spacing=1.08
def para(t="",size=11,bold=False,italic=False,color=None,align=None,after=6):
    p=doc.add_paragraph(); p.paragraph_format.space_after=DPt(after)
    if align: p.alignment=align
    r=p.add_run(t); r.bold=bold; r.italic=italic; r.font.size=DPt(size)
    if color is not None: r.font.color.rgb=color
    return p
def h1(t):
    p=doc.add_heading(t,1)
    for r in p.runs: r.font.color.rgb=DRGB(0x1a,0x1a,0x2e)
def h2(t):
    p=doc.add_heading(t,2)
    for r in p.runs: r.font.color.rgb=DRGB(0x2e,0x86,0xab)
def bl(t,lead=None):
    p=doc.add_paragraph(style="List Bullet"); p.paragraph_format.space_after=DPt(3)
    if lead:
        r=p.add_run(lead); r.bold=True; p.add_run(t)
    else: p.add_run(t)
def fig(path,cap,w=6.0):
    if os.path.exists(path):
        doc.add_picture(path,width=DIn(w)); doc.paragraphs[-1].alignment=AL.CENTER
    para(cap,size=9.5,italic=True,color=DRGB(0x55,0x5b,0x6e),align=AL.CENTER,after=12)

t=doc.add_heading("Cell-type interaction-eQTL mapping of the vestibular maculae",0)
for r in t.runs: r.font.color.rgb=DRGB(0x1a,0x1a,0x2e)
para("Full report of the mapping step",size=13,italic=True,color=DRGB(0x55,0x5b,0x6e),align=AL.CENTER)
para("Mouse vestibular project · 2026",size=11,color=DRGB(0x55,0x5b,0x6e),align=AL.CENTER,after=12)

h1("1. Overview")
para("Standard cis-eQTL mapping reports the average genetic effect on expression across all cell types in a "
     "tissue. This analysis instead asks in which cell type a regulatory effect acts, by testing whether a "
     "variant's effect scales with a cell type's abundance — a cell-type interaction eQTL (ieQTL), following "
     "Kim-Hellmuth et al. (GTEx, 2020). Cell-type abundance is the ssGSEA marker-enrichment score from the "
     "deconvolution step (absolute fractions from Bisque/MuSiC were shown to be unreliable for this "
     "nucleus-reference/whole-cell-bulk setting). We ran four analyses and compare them with each other and "
     "with the standard cis-eQTLs mapped in Pantry.")

h1("2. Methods")
para("For each analysis and cell type we fitted, per gene–cis-variant pair within 1 Mb:")
para("expression ~ genotype + covariates + enrichment + genotype × enrichment,",italic=True,align=AL.CENTER)
para("and tested the genotype × enrichment term with tensorQTL 1.0.10 on an NVIDIA A100 GPU.")
bl("MAF ≥ 0.05 in both the upper and lower halves of the enrichment distribution "
   "(so an interaction cannot be driven by a genotype seen only at one extreme of abundance); "
   "eigenMT for the effective number of cis-variants per gene, then Benjamini-Hochberg across genes; "
   "ieGenes called at FDR < 0.05.","Safeguards: ")
bl("per-group full (single age, Pantry covariates incl. 8 RNA PCs); per-group core (genotype PCs + batch); "
   "pooled core (ages pooled within organ, n≈130; genotype PCs + age + batch); pooled full (+ RNA PC1–8, the same recipe as per-group full). "
   "Pooled genotypes were merged on shared SNPs and pooled expression re-normalized with Pantry's own "
   "script.","Four analyses: ")

h1("3. Standard eQTL vs. interaction eQTL")
para("The standard cis-eQTL scan identifies thousands of eGenes per group (~6,600 on average), capturing every "
     "marginal genetic effect on expression. The interaction test is a harder, different question and returns a "
     "focused set. Of the 42 per-group ieGenes, 37 are also standard eGenes (now localized to a cell type) and "
     "5 are interaction-only — regulation invisible to a bulk eQTL scan.")
fig(os.path.join(IEP,"A_egene_vs_iegene_counts.png"),
    "Figure 1. Standard eGenes vs. interaction ieGenes per group (log scale).",5.6)

h1("4. Results by analysis")
h2("4.1 Per-group")
para(f"Under the full covariate model, the four age-groups yielded {pg_full} ieGenes; under the core model "
     f"(genotype PCs + batch) they yielded {pg_core}. In both, hair and supporting cells in the week-8 groups "
     "carry almost all of the signal, and the week-22 groups (smaller n) yield essentially none.")
fig(os.path.join(IEP,"B_iegene_heatmap.png"),
    "Figure 2. Per-group ieGenes (FDR<0.05) by group and cell type (full model).",5.4)
h2("4.2 Pooled")
para(f"Pooling the two ages within each organ (saccule n=128, utricle n=132) increased discovery to {pl_core} "
     f"ieGenes under the core model and {pl_full} under the full model, roughly doubling the per-group yield "
     "as expected from the larger sample size.")
fig(os.path.join(MR,"M1_approach_totals.png"),
    "Figure 3. ieGene discovery across all four analyses and standard eQTL (log scale).",5.4)
para("Across every analysis, hair and supporting cells are the consistently strongest cell types; rarer types "
     "appear only sporadically (Figure 4).")
fig(os.path.join(MR,"M2_by_celltype.png"),
    "Figure 4. ieGenes by cell type across the four analyses.",6.2)

h1("5. Covariate sensitivity (RNA PCs and cell composition)")
para("Expression (RNA) PCs used as covariates can themselves capture cell composition, so including them risks "
     "regressing out the interaction being tested. We confirmed this — for example saccule expression-PC1 "
     "correlates r = −0.75 with the supporting-cell score and utricle expression-PC2 r = +0.66 with the "
     "hair-cell score. The effect on discovery is not uniform. Overall the full model finds more ieGenes than "
     "core (pooled 120 vs 79): at n≈130 the denoising gain outweighs the composition absorbed, and crucially the "
     "abundant hair and supporting cells GAIN ieGenes under the full model — so their signals are strengthened, "
     "not weakened, and are not artifacts of under-adjustment. Conversely, several low-count signals vanish once "
     "RNA PCs are added (transitional-epithelial in saccule 3→1; glia 1→0 in utricle); these "
     "are covariate-sensitive and lower-confidence, the over-adjustment pattern to watch for.")
fig(os.path.join(SP,"S1_core_vs_full.png"),
    "Figure 5. Pooled ieGenes per cell type: core (no RNA PCs) vs full (+RNA PC1–8).",6.2)
fig(os.path.join(SP,"S2_rnaPC_enrichment_heatmap.png"),
    "Figure 6. Correlation of each RNA PC with each cell type's enrichment score.",6.2)

h1(f"6. Robust ieGenes ({nrob})")
para(f"We report the {nrob} ieGenes significant under BOTH pooled covariate models as the robust, "
     "covariate-insensitive set (robust_ieGenes.csv), dominated by hair and supporting cells. Examples include "
     "Cadps (a Ca²⁺-dependent secretion activator, recurring in both hair and supporting cells), Acpp in hair cells and Ralyl in fibrocytes "
     "in hair cells, Kcnf1 (a potassium channel) in utricle hair cells, and Slc6a4 in utricle roof/dark cells.")

h1("7. Caveats and next steps")
bl("robust signal is concentrated in hair and supporting cells; rarer cell types remain underpowered.")
bl("pooled and per-group analyses share few exact genes because pooling changes the model (age covariate, "
   "re-normalization); treat cell-type assignments, not individual variants, as the durable result.")
bl("report the covariate-insensitive robust set and flag covariate-sensitive hits explicitly.")
bl("next: colocalize robust ieGenes with the vestibular-function GWAS (VsEP, beam-walk) to assign balance loci "
   "to the cell type in which they act — the payoff of the deconvolution → interaction pipeline.")
para("Scripts: scripts/05,06,07,08,09,10 (+ sbatch runners). "
     "Outputs: results/ieqtl/, ieqtl_core/, ieqtl_pooled/, mapping_report/.",
     size=10, italic=True, color=DRGB(0x55,0x5b,0x6e))
DOC=os.path.join(MR,"Mapping_Full_Report.docx"); doc.save(DOC)
print("SAVED",DECK,"|",len(prs.slides._sldIdLst),"slides")
print("SAVED",DOC)
