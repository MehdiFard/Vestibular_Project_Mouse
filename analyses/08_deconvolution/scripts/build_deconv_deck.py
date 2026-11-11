#!/usr/bin/env python
"""Deconvolution report deck (PowerPoint). Run with Env1 python (has python-pptx)."""
import os
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN

RES = "../../analyses/08_deconvolution/results"
PLT = os.path.join(RES, "plots")
INTRO = os.path.join(RES, "intro_assets")
OUT = os.path.join(os.path.dirname(RES), "report", "deconvolution_report.pptx")

INK=RGBColor(0x1a,0x1a,0x2e); ACC=RGBColor(0x2e,0x86,0xab); GREY=RGBColor(0x55,0x5b,0x6e)
RED=RGBColor(0xb0,0x2a,0x2a); GRN=RGBColor(0x1b,0x7a,0x3d)
prs=Presentation(); prs.slide_width=Inches(13.333); prs.slide_height=Inches(7.5)
BLANK=prs.slide_layouts[6]
SW, SH = prs.slide_width, prs.slide_height

def slide():
    return prs.slides.add_slide(BLANK)
def box(s, l, t, w, h):
    tb=s.shapes.add_textbox(Inches(l),Inches(t),Inches(w),Inches(h)); tb.text_frame.word_wrap=True
    return tb.text_frame
def setpar(p, text, size=18, bold=False, color=INK, italic=False, bullet=False, lvl=0, after=6):
    p.text=text; p.level=lvl
    for r in p.runs: r.font.size=Pt(size); r.font.bold=bold; r.font.color.rgb=color; r.font.italic=italic
    p.space_after=Pt(after)
    return p
def title(s, text, sub=None):
    tf=box(s,0.6,0.35,12.1,1.1); setpar(tf.paragraphs[0],text,30,True,INK)
    bar=s.shapes.add_shape(1, Inches(0.62), Inches(1.15), Inches(2.4), Pt(3))
    bar.fill.solid(); bar.fill.fore_color.rgb=ACC; bar.line.fill.background()
    if sub:
        p=tf.add_paragraph(); setpar(p,sub,15,False,GREY,italic=True)
def bullets(s, items, l=0.7, t=1.5, w=6.0, h=5.5, size=16):
    tf=box(s,l,t,w,h)
    for i,(txt,lvl,col,bold) in enumerate(items):
        p=tf.paragraphs[0] if i==0 else tf.add_paragraph()
        setpar(p, ("• "+txt) if lvl==0 else ("   – "+txt), size-(1 if lvl else 0), bold, col, after=7); p.level=0
def picture(s, path, l, t, w=None, h=None):
    if os.path.exists(path):
        return s.shapes.add_picture(path, Inches(l), Inches(t),
                                    Inches(w) if w else None, Inches(h) if h else None)
def footer(s, n):
    tf=box(s,11.4,7.05,1.7,0.4); setpar(tf.paragraphs[0],f"{n}",10,False,GREY,after=0)

# 1 — title
s=slide()
r=s.shapes.add_shape(1, 0, 0, SW, Inches(2.7)); r.fill.solid(); r.fill.fore_color.rgb=INK; r.line.fill.background()
tf=box(s,0.8,0.9,11.7,1.6); setpar(tf.paragraphs[0],"Cell-type deconvolution of bulk vestibular RNA-seq",34,True,RGBColor(255,255,255))
p=tf.add_paragraph(); setpar(p,"Toward cell-type interaction-eQTL mapping (GTEx / Kim-Hellmuth, Science 2020)",18,False,RGBColor(0xC8,0xE6,0xF0),italic=True)
tf=box(s,0.8,3.2,11.7,3.2)
for i,(t,c) in enumerate([
    ("Bulk saccule (n=128) & utricle (n=132), ages pooled within organ; snRNA reference (4 donors, 6 cell types).",INK),
    ("Three methods compared: Bisque, MuSiC, and ssGSEA marker enrichment.",INK),
    ("Validated on pseudobulk; characterised on real bulk. Deconvolution step only — QTL mapping is next.",GREY)]):
    p=tf.paragraphs[0] if i==0 else tf.add_paragraph(); setpar(p,"• "+t,17,False,c,after=10)

# --- intro slides (concept, formulas, illustrations) ---
s=slide(); title(s,"Bulk RNA-seq is a mixture of cell types",
                 "Standard eQTLs average genetic effects over every cell type")
picture(s, os.path.join(INTRO,"intro1_mixture.png"), 2.42, 1.45, w=8.5)
bullets(s,[("A variant that acts in only one cell type is diluted by all the others — its effect is hidden in bulk.",0,INK,False)],
        l=1.0, t=6.55, w=11.3, size=15)

s=slide(); title(s,"Deconvolution: estimating cell-type composition",
                 "Factorize bulk into cell-type signatures × per-sample fractions")
picture(s, os.path.join(INTRO,"intro2_deconv.png"), 2.57, 1.45, w=8.2)
bullets(s,[("Reference signatures S come from single-cell data; three methods are compared here (Bisque, MuSiC, ssGSEA).",0,INK,False)],
        l=1.0, t=6.55, w=11.3, size=15)

s=slide(); title(s,"Why it matters: cell-type interaction eQTL",
                 "Recover cell-type-specific genetic regulation from bulk (GTEx / Kim-Hellmuth 2020)")
picture(s, os.path.join(INTRO,"intro3_interaction.png"), 2.77, 1.45, w=7.8)
bullets(s,[("The genotype × abundance interaction term is the cell-type-specific effect — the goal this deconvolution feeds.",0,INK,False)],
        l=1.0, t=6.55, w=11.3, size=15)

# 2 — background
s=slide(); title(s,"Why deconvolution?","Bulk eQTLs average over cell types and hide cell-type-specific regulation")
bullets(s,[
 ("A bulk sample is a mixture: its expression is a weighted average across all cell types present.",0,INK,False),
 ("A variant that regulates a gene only in hair cells is diluted by every other cell type → missed or averaged-out in a standard eQTL scan.",0,INK,False),
 ("GTEx/Kim-Hellmuth solution — cell-type interaction eQTLs (ieQTLs):",0,ACC,True),
 ("estimate each cell type's abundance per bulk sample (deconvolution),",1,INK,False),
 ("test genotype × cell-type-abundance interaction on expression,",1,INK,False),
 ("the interaction term = a cell-type-specific regulatory effect, recovered from bulk.",1,INK,False),
 ("This report delivers the deconvolution (the per-sample cell-type estimates); mapping follows.",0,GREY,False),
], w=12.0)

# 3 — design
s=slide(); title(s,"Design & data","Ages pooled within organ; organs separate")
bullets(s,[
 ("Bulk (kallisto TPM): saccule = w8+w22 (128), utricle = w8+w22 (132).",0,INK,False),
 ("Reference: cleaned snRNA atlas, 4 donors (FTL_1146–1149), ≤2000 cells/type → 11,094 cells.",0,INK,False),
 ("6 abundant cell types: Hair_cell (I/II collapsed), Supporting_cell, Transitional_epith, Glia_Schwann, Fibrocyte_mesench, Roof_dark_cell.",0,INK,False),
 ("Genes: protein-coding only, mitochondrial & ribosomal removed (16,840 common genes).",0,INK,False),
 ("Methods:",0,ACC,True),
 ("Bisque — learns a single-cell→bulk transform (planned primary),",1,INK,False),
 ("MuSiC — weights genes by cross-donor consistency (sensitivity),",1,INK,False),
 ("ssGSEA — GTEx-style marker-enrichment score.",1,INK,False),
], w=12.2)

# 4 — pseudobulk validation
s=slide(); title(s,"Validation: pseudobulk self-recovery","Deconvolve pseudobulk built from the single cells (composition known)")
picture(s, os.path.join(PLT,"00_pseudobulk_validation.png"), 6.7, 1.5, w=6.2)
bullets(s,[
 ("Mechanics are correct — both methods run and recover known composition.",0,INK,False),
 ("MuSiC r = 0.95 (points hug the diagonal).",0,GRN,True),
 ("Bisque r = 0.45 (compressed toward uniform).",0,RED,True),
 ("MuSiC out-performs Bisque in-domain — the opposite of the a-priori plan.",0,INK,False),
], l=0.6, w=6.0, t=1.7, size=16)

# 5 — real bulk fractions fail
s=slide(); title(s,"On real bulk, absolute fractions are unreliable","Not a bug — a nuclei-vs-whole-cell domain shift")
picture(s, os.path.join(PLT,"01_mean_composition_by_method.png"), 6.6, 1.6, w=6.4)
bullets(s,[
 ("MuSiC collapses (~0.85–0.93 Supporting_cell, little variation).",0,RED,False),
 ("Bisque returns near-uniform proportions (~0.18 each).",0,RED,False),
 ("Marker-restricting genes made MuSiC worse (~1.0 Supporting).",0,INK,False),
 ("Cause:",0,ACC,True),
 ("single-NUCLEUS reference vs whole-cell / total-RNA bulk,",1,INK,False),
 ("incomplete 6-type reference for a whole-organ dissection.",1,INK,False),
 ("→ Do not read Bisque/MuSiC fractions as real composition.",0,RED,True),
], l=0.6, w=5.9, t=1.6, size=15)

# 6 — method disagreement
s=slide(); title(s,"Methods disagree on real bulk","Per-sample estimates barely correlate")
picture(s, os.path.join(PLT,"02_bisque_vs_music_scatter.png"), 5.9, 1.5, w=7.0)
bullets(s,[
 ("Median Spearman across cell types:",0,ACC,True),
 ("Bisque vs MuSiC: 0.18",1,RED,False),
 ("Bisque vs ssGSEA: 0.49",1,INK,False),
 ("MuSiC vs ssGSEA: 0.41",1,INK,False),
 ("ssGSEA is the consensus signal both fraction methods partially track.",0,INK,False),
], l=0.55, w=5.2, t=1.6, size=15)

# 7 — ssGSEA usable
s=slide(); title(s,"The usable signal: ssGSEA marker enrichment","Rank-based → robust to the nuclei-vs-bulk shift (what GTEx used)")
picture(s, os.path.join(PLT,"06_ssgsea_enrichment.png"), 6.5, 1.6, w=6.5)
bullets(s,[
 ("Per-sample cell-type abundance proxy; not a fraction, but varies with true abundance.",0,INK,False),
 ("Real cross-sample variation (Hair_cell 0.01–0.57).",0,GRN,True),
 ("Best-behaved types: Hair_cell, Transitional_epith, Glia_Schwann, Roof_dark_cell.",0,INK,False),
 ("This is the per-sample covariate that feeds interaction-eQTL mapping.",0,ACC,True),
], l=0.6, w=5.8, t=1.7, size=15)

# 8 — summary / next
s=slide(); title(s,"Summary & next steps")
bullets(s,[
 ("Deconvolution pipeline built (scripts/), 3 methods, validated on pseudobulk.",0,INK,False),
 ("Absolute-fraction deconvolution is not reliable here (nuclei-vs-bulk domain shift + incomplete reference).",0,INK,False),
 ("MuSiC > Bisque in-domain (r 0.95 vs 0.45); neither trustworthy on real bulk.",0,INK,False),
 ("USE the ssGSEA enrichment scores as the cell-type covariate for mapping.",0,GRN,True),
 ("Next — cell-type interaction-eQTL mapping (tensorQTL):",0,ACC,True),
 ("ssGSEA score inverse-normal transformed as the interaction term,",1,INK,False),
 ("MAF ≥ 0.05 in both halves of the enrichment distribution,",1,INK,False),
 ("eigenMT + Benjamini-Hochberg, FDR < 0.05.",1,INK,False),
], w=12.2)

for i,s in enumerate(prs.slides,1): footer(s,i)
prs.save(OUT); print("SAVED", OUT, "| slides:", len(prs.slides._sldIdLst))
