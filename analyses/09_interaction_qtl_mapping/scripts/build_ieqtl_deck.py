#!/usr/bin/env python
"""Interaction-eQTL mapping report deck (PowerPoint). Env1 python."""
import os
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
RES="../../analyses/09_interaction_qtl_mapping/results"
IE=os.path.join(RES,"ieqtl"); PLT=os.path.join(IE,"plots")
OUT=os.path.join(IE,"ieqtl_mapping_report.pptx")
INK=RGBColor(0x1a,0x1a,0x2e); ACC=RGBColor(0x2e,0x86,0xab); GREY=RGBColor(0x55,0x5b,0x6e)
RED=RGBColor(0xb0,0x2a,0x2a); GRN=RGBColor(0x1b,0x7a,0x3d)
prs=Presentation(); prs.slide_width=Inches(13.333); prs.slide_height=Inches(7.5)
BLANK=prs.slide_layouts[6]; SW,SH=prs.slide_width,prs.slide_height
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
def bullets(s,items,l=0.7,t=1.5,w=6.0,sz=16):
    tf=box(s,l,t,w,5.5)
    for i,(txt,lvl,col,bd) in enumerate(items):
        p=tf.paragraphs[0] if i==0 else tf.add_paragraph()
        setp(p,("• "+txt) if lvl==0 else ("   – "+txt),sz-(1 if lvl else 0),bd,col,after=7)
def pic(s,path,l,t,w=None,h=None):
    if os.path.exists(path): s.shapes.add_picture(path,Inches(l),Inches(t),Inches(w) if w else None,Inches(h) if h else None)
def foot(s,n): setp(box(s,11.6,7.05,1.6,0.4).paragraphs[0],str(n),10,False,GREY,after=0)

# 1 title
s=slide()
r=s.shapes.add_shape(1,0,0,SW,Inches(2.7)); r.fill.solid(); r.fill.fore_color.rgb=INK; r.line.fill.background()
setp(box(s,0.8,0.9,11.7,1.4).paragraphs[0],"Cell-type interaction-eQTL mapping",34,True,RGBColor(255,255,255))
setp(box(s,0.8,2.0,11.7,0.7).paragraphs[0],"Vestibular maculae (mouse) — GTEx / Kim-Hellmuth 2020 approach, tensorQTL on GPU",17,False,RGBColor(0xC8,0xE6,0xF0),it=True)
tf=box(s,0.8,3.2,11.7,3.2)
for i,(t,c) in enumerate([
 ("Interaction test: expression ~ genotype + covariates + cell-type score + genotype×score.",INK),
 ("Cell-type score = ssGSEA marker enrichment (inverse-normal transformed); 6 cell types.",INK),
 ("Per Pantry age-group; MAF≥0.05 both halves; eigenMT + BH, FDR<0.05. 42 ieGenes.",GREY)]):
    setp((tf.paragraphs[0] if i==0 else tf.add_paragraph()),"• "+t,17,False,c,after=10)

# 2 method
s=slide(); title(s,"Method","One interaction scan per group × cell type")
bullets(s,[
 ("Genotypes, quantile-normalized expression, and covariates: the exact Pantry per-group inputs.",0,INK,False),
 ("Interaction covariate: ssGSEA cell-type enrichment, inverse-normal transformed across samples.",0,INK,False),
 ("tensorQTL cis.map_nominal (v1.0.10) on an A100 GPU, 1 Mb cis-window.",0,INK,False),
 ("Safeguards (as in the paper):",0,ACC,True),
 ("MAF ≥ 0.05 in BOTH halves of the enrichment distribution (maf_threshold_interaction),",1,INK,False),
 ("eigenMT (effective tests/gene) + Benjamini-Hochberg across genes, FDR < 0.05.",1,INK,False),
 ("4 groups × 6 cell types = 24 scans; ~2.5 min each on GPU.",0,GREY,False),
],w=12.2)

# 3 results heatmap
s=slide(); title(s,"Results — 42 ieGenes (FDR<0.05)","Hair & supporting cells dominate; discovery tracks sample size")
pic(s,os.path.join(PLT,"B_iegene_heatmap.png"),6.5,1.5,w=6.6)
bullets(s,[
 ("Hair_cell is the strongest type in BOTH organs (saccule_w8=10, utricle_w8=8).",0,GRN,True),
 ("Supporting_cell next (9 / 3).",0,INK,False),
 ("Both week-22 groups → ~0 ieGenes.",0,RED,True),
 ("→ power tracks n (w8 n≈73–74 vs w22 n≈55–58), as in the bulk xQTLs.",0,INK,False),
],l=0.55,w=5.8,t=1.7,sz=15)

# 4 vs standard eQTL
s=slide(); title(s,"vs. standard cis-eQTL","Interaction isolates a focused, cell-type-specific subset")
pic(s,os.path.join(PLT,"A_egene_vs_iegene_counts.png"),6.3,1.6,w=6.8)
bullets(s,[
 ("Standard eQTL: thousands of eGenes/group (marginal effects, cell-type-averaged).",0,INK,False),
 ("Interaction eQTL: a focused set whose effect depends on cell-type abundance.",0,INK,False),
 ("Different questions — counts are not directly comparable.",0,GREY,False),
 ("The interaction test answers 'in which cell type does the effect act?'",0,ACC,True),
],l=0.55,w=5.6,t=1.7,sz=15)

# 5 overlap
s=slide(); title(s,"Are the ieGenes new?","37/42 also standard eGenes; 5 interaction-only")
pic(s,os.path.join(PLT,"C_iegene_overlap.png"),0.7,1.7,w=5.2)
pic(s,os.path.join(PLT,"D_marginal_vs_interaction.png"),6.4,1.6,w=5.0)
bullets(s,[
 ("Most ieGenes also have a marginal eQTL — but now attributed to a cell type.",0,INK,False),
 ("5 are interaction-only (no standard eGene): the genuinely novel cell-type-specific signals.",0,GRN,True),
],l=0.6,w=12.0,t=6.2,sz=14)

# 6 examples
s=slide(); title(s,"Example ieGenes")
bullets(s,[
 ("Hair_cell (saccule w8): Slc7a11, Fam124a, Nedd4, Tlr3, Fstl4.",0,INK,False),
 ("Hair_cell (utricle w8): Slk, Suclg2, Jade2, Zfp658.",0,INK,False),
 ("Supporting_cell (saccule w8): Znrf3 (strongest, BH=3e-4), Col20a1, Mrap, Limk2.",0,INK,False),
 ("Roof/dark (utricle w8): Pla2g7, Gdpd5.  Transitional (utricle w8): Lrrc6.",0,INK,False),
 ("Caveats:",0,ACC,True),
 ("no exact-gene replication across age/organ (low power, different animals),",1,RED,False),
 ("within-organ recurrence (Znrf3, Slc7a11) partly reflects correlated enrichment scores,",1,RED,False),
 ("a few technical/low-confidence hits (Rps26, Olfr1193, Gm*/Rik).",1,RED,False),
],w=12.2)

# 7 next
s=slide(); title(s,"Caveats & next steps")
bullets(s,[
 ("Full list: results/ieqtl/ieGenes_significant.csv (+ summary, comparison table).",0,GREY,False),
 ("Underpowered at per-age-group n — w22 yields ~0.",0,INK,False),
 ("NEXT — pool ages within organ (n≈130) to roughly double power:",0,ACC,True),
 ("merge genotypes, re-normalize expression, recompute covariates, re-map.",1,INK,False),
 ("Covariate sensitivity: Pantry covariates include 8 RNA PCs that may capture cell composition;",0,ACC,True),
 ("re-run core model (no RNA PCs) vs +RNA PC1-8 — a hit that vanishes only with a PC correlated to the",1,INK,False),
 ("tested cell abundance may be overadjustment, not a false interaction.",1,INK,False),
],w=12.4)

for i,s in enumerate(prs.slides,1): foot(s,i)
prs.save(OUT); print("SAVED",OUT,"slides",len(prs.slides._sldIdLst))
