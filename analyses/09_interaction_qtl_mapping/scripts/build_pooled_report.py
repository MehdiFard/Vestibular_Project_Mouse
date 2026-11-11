#!/usr/bin/env python
"""Pooled-ages ieQTL report: deck (pptx) + detailed Word (docx). Env1 python."""
import os
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from docx import Document
from docx.shared import Pt as DPt, RGBColor as DRGB, Inches as DIn
from docx.enum.text import WD_ALIGN_PARAGRAPH as AL

RES="../../analyses/09_interaction_qtl_mapping/results"
POOL=os.path.join(RES,"ieqtl_pooled"); PLT=os.path.join(POOL,"plots")
INK=RGBColor(0x1a,0x1a,0x2e); ACC=RGBColor(0x2e,0x86,0xab); GREY=RGBColor(0x55,0x5b,0x6e)
RED=RGBColor(0xb0,0x2a,0x2a); GRN=RGBColor(0x1b,0x7a,0x3d)

# ===================== DECK =====================
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
    tf=box(s,0.6,0.35,12.1,1.1); setp(tf.paragraphs[0],t,30,True,INK)
    bar=s.shapes.add_shape(1,Inches(0.62),Inches(1.15),Inches(2.4),Pt(3)); bar.fill.solid(); bar.fill.fore_color.rgb=ACC; bar.line.fill.background()
    if sub: setp(tf.add_paragraph(),sub,15,False,GREY,it=True)
def bull(s,items,l=0.7,t=1.5,w=6.0,sz=16):
    tf=box(s,l,t,w,5.5)
    for i,(txt,lvl,col,bd) in enumerate(items):
        p=tf.paragraphs[0] if i==0 else tf.add_paragraph()
        setp(p,("• "+txt) if lvl==0 else ("   – "+txt),sz-(1 if lvl else 0),bd,col,after=7)
def pic(s,path,l,t,w=None,h=None):
    if os.path.exists(path): s.shapes.add_picture(path,Inches(l),Inches(t),Inches(w) if w else None,Inches(h) if h else None)
def foot(s,n): setp(box(s,11.6,7.05,1.6,0.4).paragraphs[0],str(n),10,False,GREY,after=0)

s=slide()
r=s.shapes.add_shape(1,0,0,SW,Inches(2.7)); r.fill.solid(); r.fill.fore_color.rgb=INK; r.line.fill.background()
setp(box(s,0.8,0.9,11.7,1.4).paragraphs[0],"Pooled-ages interaction-eQTL mapping",34,True,RGBColor(255,255,255))
setp(box(s,0.8,2.0,11.7,0.7).paragraphs[0],"Higher power + covariate (RNA-PC) sensitivity analysis",17,False,RGBColor(0xC8,0xE6,0xF0),it=True)
tf=box(s,0.8,3.2,11.7,3.2)
for i,(t,c) in enumerate([
 ("Ages pooled within organ: saccule n=128, utricle n=132 (vs 55–74 per age-group).",INK),
 ("Merged genotypes (shared SNPs) + Pantry-identical re-normalized expression + recomputed covariates.",INK),
 ("Two covariate models: core (geno PCs + age + batch) vs full (+ RNA PC1-8).",GREY)]):
    setp((tf.paragraphs[0] if i==0 else tf.add_paragraph()),"• "+t,17,False,c,after=10)

s=slide(); title(s,"Pooling boosts discovery","More samples → more ieGenes")
pic(s,os.path.join(PLT,"S3_approach_comparison.png"),6.5,1.7,w=6.4)
bull(s,[
 ("per-group ieGenes (4 age-groups): 42",0,INK,False),
 ("pooled, core covariates: 79",0,GRN,True),
 ("pooled, full covariates: 120",0,GRN,True),
 ("robust set (significant in BOTH models): 49",0,ACC,True),
 ("Standard eQTL still finds ~6,600 eGenes/group (different, marginal test).",0,GREY,False),
],l=0.55,w=5.7,t=1.8,sz=15)

s=slide(); title(s,"Covariate sensitivity (your caution)","RNA PCs can capture cell composition")
pic(s,os.path.join(PLT,"S1_core_vs_full.png"),0.5,1.55,w=8.2)
bull(s,[
 ("Overall full>core: at n≈130 the RNA-PC denoising gain outweighs the composition they absorb.",0,INK,False),
 ("Hair & supporting GAIN strongly under full (sac hair 23→40, supp 33→48) → robust.",0,GRN,True),
 ("Saccule Transitional 3→1, utricle Glia 1→0 → covariate-sensitive, low confidence.",0,RED,True),
 ("→ low-count, covariate-sensitive = lower confidence.",0,RED,False),
],l=8.6,w=4.4,t=1.6,sz=13)

s=slide(); title(s,"Which RNA PCs capture composition?","Correlation of RNA PCs with cell-type enrichment")
pic(s,os.path.join(PLT,"S2_rnaPC_enrichment_heatmap.png"),1.2,1.7,w=10.9)
bull(s,[
 ("saccule exprPC1 ↔ Supporting r=-0.75; exprPC2 ↔ Hair +0.45.  utricle exprPC2 ↔ Hair +0.66; exprPC4 ↔ Roof +0.57.",0,INK,False),
 ("A signal that vanishes only when a composition-correlated PC is added = over-adjustment, not a false interaction.",0,ACC,True),
],l=0.7,w=12.0,t=6.1,sz=13)

s=slide(); title(s,"Robust ieGenes (49) & takeaways")
bull(s,[
 ("Confident set = significant under BOTH covariate models (result/.../ieqtl_pooled/robust_ieGenes.csv).",0,GREY,False),
 ("Hair cell (saccule): Cadps, Acpp, Gm15816 …  (utricle): Kcnf1, Gm15816.",0,INK,False),
 ("Supporting cell (saccule): Cadps, Filip1, Golm1 … ;  Roof/dark (utricle): Slc6a4, Add2.",0,INK,False),
 ("Cadps recurs in hair + supporting — Ca2+-dependent secretion, plausible in sensory epithelium.",0,INK,False),
 ("Comparison caveat:",0,ACC,True),
 ("pooled ≠ per-group analysis (age covariate + re-normalization) → limited exact-gene overlap;",1,RED,False),
 ("mixing ages can itself induce/erase interactions, so pooled and per-group are complementary.",1,RED,False),
 ("Next: colocalize robust ieGenes with VsEP / beam-walk GWAS to assign balance loci to cell types.",0,ACC,True),
],w=12.4)

for i,s in enumerate(prs.slides,1): foot(s,i)
DECK=os.path.join(POOL,"pooled_ieqtl_report.pptx"); prs.save(DECK)

# ===================== WORD =====================
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

t=doc.add_heading("Pooled-ages cell-type interaction-eQTL mapping",0)
for r in t.runs: r.font.color.rgb=DRGB(0x1a,0x1a,0x2e)
para("Higher-powered mapping and covariate-sensitivity analysis",size=13,italic=True,color=DRGB(0x55,0x5b,0x6e),align=AL.CENTER)
para("Mouse vestibular project · 2026",size=11,color=DRGB(0x55,0x5b,0x6e),align=AL.CENTER,after=12)

h1("1. Motivation")
para("The per-age-group interaction scans (saccule/utricle at weeks 8 and 22) were underpowered: the week-22 "
     "groups (n≈55–58) returned essentially no ieGenes. Pooling the two ages within each organ roughly doubles "
     "the sample size (saccule n=128, utricle n=132) and should recover more cell-type-specific signal. Because "
     "pooling mixes ages, we add an age covariate and treat the pooled and per-group analyses as complementary "
     "rather than identical.")

h1("2. Building the pooled inputs")
bl("the two ages' genotypes were merged on their shared SNPs (saccule 103,345; utricle 103,372) and "
   "genotype PCs recomputed on the pooled sample.", "Genotypes: ")
bl("the two ages' gene TPM were concatenated and re-normalized with Pantry's own normalize_phenotypes.py "
   "(quantile-normalize → inverse-normal, dropping >50%-zero / invariant genes), so processing matches the "
   "original pipeline exactly.", "Expression: ")
bl("two sets were built — core (genotype PC1–3 + age + batch) and full (core + expression/RNA PC1–8) — to test "
   "covariate sensitivity.", "Covariates: ")
para("Mapping used the same GPU tensorQTL interaction model and safeguards as before (MAF≥0.05 in both halves "
     "of the enrichment distribution; eigenMT + Benjamini–Hochberg; FDR<0.05), with the ssGSEA enrichment score "
     "(inverse-normal transformed) as the interaction term.")

h1("3. Pooling increases discovery")
para("Pooling raised the ieGene count from 42 (summed across the four age-groups) to 79 under the core model "
     "and 120 under the full model. Standard cis-eQTL still identifies far more eGenes (~6,600 per group) because "
     "it is a different, marginal test; the interaction scan isolates the cell-type-specific subset.")
fig(os.path.join(PLT,"S3_approach_comparison.png"),
    "Figure 1. ieGene discovery across approaches (log scale): per-group vs pooled (core, full) vs standard eGenes.",5.6)

h1("4. Covariate sensitivity (RNA PCs and cell composition)")
para("A concern with this design is that expression (RNA) PCs used as covariates can themselves capture cell "
     "composition, so including them may regress out the very interaction being tested. We confirmed the RNA PCs "
     "do track composition — for example saccule expression-PC1 correlates r = −0.75 with the supporting-cell "
     "enrichment score, and utricle expression-PC2 correlates r = +0.66 with the hair-cell score (Figure 3).")
fig(os.path.join(PLT,"S1_core_vs_full.png"),
    "Figure 2. Pooled ieGenes per cell type under core (no RNA PCs) vs full (+RNA PC1–8).",6.4)
para("The effect is not uniform. Overall the full model finds more ieGenes than core (120 vs 79): at n≈130 the "
     "denoising gain from the RNA PCs outweighs the composition they absorb. Importantly, the two abundant "
     "cell types — hair and supporting cells — GAIN ieGenes under the full model (hair 23→40, supporting 33→48 "
     "in saccule), so their signals are strengthened, not weakened, by adding RNA PCs; they are not artifacts of "
     "under-adjustment. Conversely, several low-count signals disappear once RNA PCs are added — transitional-"
     "epithelial cells in saccule (3→1) and glia cells in utricle (1→0). "
     "These are covariate-sensitive and lower-confidence, exactly the over-adjustment pattern to be wary of, since "
     "the RNA PCs partly encode those cell types' abundance.")
fig(os.path.join(PLT,"S2_rnaPC_enrichment_heatmap.png"),
    "Figure 3. Correlation of each RNA PC with each cell type's enrichment score.",6.4)
para("We therefore report the 49 ieGenes that are significant under BOTH covariate models as the robust, "
     "covariate-insensitive set (robust_ieGenes.csv). It is dominated by hair and supporting cells and includes, "
     "for example, Cadps (a Ca²⁺-dependent secretion activator, recurring in both hair and supporting cells), "
     "Acpp and Ralyl in hair cells, and Slc6a4 in utricle roof/dark cells.")

h1("5. Comparison with the previous approaches")
bl("thousands of eGenes per group (marginal effects); the interaction scans annotate a small subset with a "
   "cell type. Complementary, not competing.", "vs standard cis-eQTL (Pantry): ")
bl("pooling recovers more ieGenes (66–99 vs 42) but the specific genes differ substantially, because pooling "
   "changes the model (age covariate, re-normalization across 128–132 samples). Limited exact-gene overlap is "
   "expected; mixing ages can itself create or mask interactions. The two analyses are complementary.",
   "vs per-age-group ieQTL: ")

h1("6. Caveats and next steps")
bl("robust ieGenes are concentrated in hair and supporting cells; rarer cell types remain underpowered.")
bl("pooled and per-group hits overlap little at the individual-gene level; treat cell-type assignments, not "
   "specific variants, as the durable result for now.")
bl("next: colocalize the robust ieGenes with the vestibular-function GWAS (VsEP, beam-walk) to assign balance "
   "loci to the cell type in which they act — the payoff of the whole deconvolution→interaction pipeline.")
para("Scripts: scripts/07a,07b,08,09. Outputs: results/ieqtl_pooled/ "
     "(pooled_ieGene_summary.csv, robust_ieGenes.csv, plots/).", size=10, italic=True, color=DRGB(0x55,0x5b,0x6e))
DOC=os.path.join(POOL,"Pooled_ieQTL_Report.docx"); doc.save(DOC)
print("SAVED", DECK, "|", len(prs.slides._sldIdLst),"slides")
print("SAVED", DOC)
