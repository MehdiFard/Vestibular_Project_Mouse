#!/usr/bin/env python
"""Lab-meeting deck for the saccule-vs-utricle per-tissue integration + DE analysis."""
import os, sys
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE
from PIL import Image

FIG_PT = "../../analyses/07_single_cell_rna/results/per_tissue"
FIG_DE = "../../analyses/07_single_cell_rna/results/tissue_DE"
OUT = sys.argv[1] if len(sys.argv) > 1 else \
    "../../analyses/07_single_cell_rna/report/saccule_vs_utricle_labmeeting.pptx"

INK  = RGBColor(0x1a, 0x1a, 0x2e); ACC = RGBColor(0x2e, 0x86, 0xab)
ACC2 = RGBColor(0xc0, 0x39, 0x2b); GREY = RGBColor(0x55, 0x5b, 0x6e)
LIGHT = RGBColor(0xf4, 0xf6, 0xf8); WHITE = RGBColor(0xff, 0xff, 0xff)
GREEN = RGBColor(0x2e, 0x7d, 0x32)

prs = Presentation(); prs.slide_width = Inches(13.333); prs.slide_height = Inches(7.5)
SW, SH = prs.slide_width, prs.slide_height
BLANK = prs.slide_layouts[6]

def slide(): return prs.slides.add_slide(BLANK)
def rect(s,x,y,w,h,c):
    sh=s.shapes.add_shape(MSO_SHAPE.RECTANGLE,x,y,w,h); sh.fill.solid()
    sh.fill.fore_color.rgb=c; sh.line.fill.background(); sh.shadow.inherit=False; return sh
def tbox(s,x,y,w,h,anchor=MSO_ANCHOR.TOP):
    tb=s.shapes.add_textbox(x,y,w,h); tf=tb.text_frame; tf.word_wrap=True
    tf.vertical_anchor=anchor; tf.margin_left=tf.margin_right=Inches(0.05)
    tf.margin_top=tf.margin_bottom=Inches(0.02); return tb,tf
def run(p,t,sz,c,bold=False,italic=False):
    r=p.add_run(); r.text=t; f=r.font; f.size=Pt(sz); f.color.rgb=c
    f.bold=bold; f.italic=italic; f.name="Calibri"; return r
def header(s,kick,title):
    rect(s,0,0,SW,Inches(1.12),INK); rect(s,0,Inches(1.12),SW,Inches(0.06),ACC)
    _,tf=tbox(s,Inches(0.5),Inches(0.10),Inches(12.3),Inches(0.98),MSO_ANCHOR.MIDDLE)
    run(tf.paragraphs[0],kick.upper(),12,ACC,bold=True)
    run(tf.add_paragraph(),title,25,WHITE,bold=True)
def bullets(s,x,y,w,h,items):
    _,tf=tbox(s,x,y,w,h)
    for i,it in enumerate(items):
        text,lvl,c,bold=(it if isinstance(it,tuple) else (it,0,INK,False))
        p=tf.paragraphs[0] if i==0 else tf.add_paragraph()
        p.space_after=Pt(6)
        run(p,("   – " if lvl else "• ") if lvl==0 or lvl==1 else "",13.5 if lvl==0 else 12,
            ACC if lvl==0 else GREY, bold=True)
        run(p,text,13.5 if lvl==0 else 12,c,bold=bold)
    return tf
def pic(s,path,x,y,w,h,cap=None):
    if not os.path.exists(path):
        b=rect(s,x,y,w,h,LIGHT); b.line.color.rgb=ACC2
        _,tf=tbox(s,x,y,w,h,MSO_ANCHOR.MIDDLE); pp=tf.paragraphs[0]; pp.alignment=PP_ALIGN.CENTER
        run(pp,"missing: "+os.path.basename(path),11,ACC2); return
    iw,ih=Image.open(path).size; br=w/h; ir=iw/ih
    dw,dh=(w,int(w/ir)) if ir>br else (int(h*ir),h)
    s.shapes.add_picture(path,x+(w-dw)//2,y+(h-dh)//2,dw,dh)
    if cap:
        _,tf=tbox(s,x,y+h+Inches(0.02),w,Inches(0.3)); pp=tf.paragraphs[0]
        pp.alignment=PP_ALIGN.CENTER; run(pp,cap,10,GREY,italic=True)
def footer(s,n):
    _,tf=tbox(s,Inches(0.4),SH-Inches(0.4),Inches(11),Inches(0.3))
    run(tf.paragraphs[0],"Saccule vs utricle — per-tissue integration & differential expression",9,GREY)
    _,tf2=tbox(s,SW-Inches(1.1),SH-Inches(0.4),Inches(0.7),Inches(0.3))
    pp=tf2.paragraphs[0]; pp.alignment=PP_ALIGN.RIGHT; run(pp,str(n),10,GREY)
def pt(n): return os.path.join(FIG_PT,n)
def de(n): return os.path.join(FIG_DE,n)

# ---- 1 title
s=slide(); rect(s,0,0,SW,SH,INK); rect(s,0,Inches(4.2),SW,Inches(0.08),ACC)
_,tf=tbox(s,Inches(0.9),Inches(1.6),Inches(11.5),Inches(2.4))
run(tf.paragraphs[0],"Saccule vs utricle: are there molecular differences?",36,WHITE,bold=True)
run(tf.add_paragraph(),"Per-tissue integration, marker annotation, and differential expression",19,ACC)
_,tf=tbox(s,Inches(0.9),Inches(4.5),Inches(11.5),Inches(1.8))
for t,sz in [("Saccule (FTL_1146 + 1148, 25,172 cells)  ·  Utricle (FTL_1147 + 1149, 24,037 cells)",15),
             ("10x 3′ · mouse · integrated per tissue, then compared cell-type-by-cell-type",13),
             ("Run on the cluster (SLURM batch job)",12)]:
    run(tf.add_paragraph(),t,sz,WHITE if sz>13 else GREY)
run(tbox(s,Inches(0.9),Inches(6.6),Inches(11),Inches(0.4))[1].paragraphs[0],
    "Lab meeting · Mehdi Esmaeili Fard",12,GREY)

# ---- 2 workflow
s=slide(); header(s,"Workflow","Four steps, run as a cluster job")
bullets(s,Inches(0.55),Inches(1.5),Inches(12.2),Inches(5.4),[
    ("1 · Integrate the two SACCULE replicates (CCA batch correction)",0,INK,True),
    ("2 · Integrate the two UTRICLE replicates (separately)",0,INK,True),
    ("3 · Marker-based annotation of each tissue (14 canonical panels → cell types), plus per-cell-type marker genes",0,INK,True),
    ("4 · Differential expression utricle vs saccule WITHIN each shared cell type",0,INK,True),
    ("volcano per cell type · targeted regional genes (Otx1/Otx2/Emx2/Tbx1…) · pseudobulk replicate check",1,GREY,False),
    ("Why per-tissue first? clean within-tissue maps; the cross-tissue DE then uses the combined object so cell types are matched",0,INK,False),
    ("Compute: SLURM batch job (16 cores / 128 G); scripts annot5 (per tissue) + annot6 (DE)",0,GREY,False),
])
footer(s,2)

# ---- 3 per-tissue annotation
s=slide(); header(s,"Steps 1–3 · annotation","Both organs recover the same cell types")
pic(s,pt("saccule_UMAP_celltype.png"),Inches(0.3),Inches(1.3),Inches(6.4),Inches(4.6),
    cap="Saccule — 25,172 cells (2 reps)")
pic(s,pt("utricle_UMAP_celltype.png"),Inches(6.7),Inches(1.3),Inches(6.4),Inches(4.6),
    cap="Utricle — 24,037 cells (2 reps)")
bullets(s,Inches(0.5),Inches(6.05),Inches(12.4),Inches(1.2),[
    ("Same cell types in both: hair cells, supporting, transitional epithelial, dark/roof, glia, fibrocyte, immune, endothelial, melanocyte",0,INK,False),
    ("Hair cells labelled 'HC_type_I' in saccule vs 'Hair_cell' in utricle is an argmax quirk (I/II not resolvable at 3′) — treat as one hair-cell compartment",0,ACC2,False),
])
footer(s,3)

# ---- 4 DE method + summary
s=slide(); header(s,"Step 4 · differential expression","utricle vs saccule, per cell type")
bullets(s,Inches(0.5),Inches(1.4),Inches(6.1),Inches(5.3),[
    ("For each shared cell type (≥40 cells/tissue): FindMarkers utricle vs saccule on RNA (log-normalized) data",0,INK,False),
    ("Significant = |log2FC| > 0.5 and adj.p < 0.05",1,GREY,False),
    ("Consistent asymmetry: many MORE genes 'up in saccule' than 'up in utricle' in the epithelial types",0,INK,True),
    ("Caveat — this uses Wilcoxon with CELLS as units, which massively inflates significance (24k cells → tiny p-values). Read the counts as relative, not absolute.",0,ACC2,False),
    ("(DESeq2/edgeR not installed → no proper pseudobulk stats; see slide 7 for the replicate-level check)",1,GREY,False),
])
rows=[("cell type","↑saccule","↑utricle"),("HC_type_I","161","63"),("Supporting_cell","164","22"),
      ("Hair_cell","78","33"),("Fibrocyte_mesench","49","42"),("Glia_Schwann","43","22"),
      ("Neuron","45","1"),("Roof_dark_cell","4","14"),("Pericyte_SMC","15","1"),
      ("Transitional_epith","0","5")]
t=s.shapes.add_table(len(rows),3,Inches(7.0),Inches(1.4),Inches(5.9),Inches(5.1)).table
for c,w in zip(range(3),[Inches(3.3),Inches(1.3),Inches(1.3)]): t.columns[c].width=w
for r in range(len(rows)):
    for c in range(3):
        cell=t.cell(r,c); cell.vertical_anchor=MSO_ANCHOR.MIDDLE
        cell.margin_top=cell.margin_bottom=Inches(0.02)
        pp=cell.text_frame.paragraphs[0]
        if c: pp.alignment=PP_ALIGN.CENTER
        hd=(r==0); run(pp,rows[r][c],11 if hd else 10.5,WHITE if hd else INK,bold=(hd or c==0))
        cell.fill.solid(); cell.fill.fore_color.rgb=INK if hd else (LIGHT if r%2 else WHITE)
footer(s,4)

# ---- 5 regional markers negative result
s=slide(); header(s,"Key result 1","The textbook regional genes do NOT separate the organs")
pic(s,de("regional_markers_DotPlot_by_tissue.png"),Inches(0.3),Inches(1.35),Inches(7.4),Inches(5.2))
bullets(s,Inches(7.9),Inches(1.5),Inches(5.1),Inches(5.0),[
    ("Otx1, Otx2, Emx2, Tbx1, Lfng, Bmp4 are the classic saccule/utricle patterning genes",0,INK,False),
    ("Here they are expressed in BOTH organs at similar levels — none is tissue-specific",0,ACC2,True),
    ("Otx2/Emx2/Tbx1 mark supporting cells in utricle AND saccule alike",1,GREY,False),
    ("None reaches significance in the DE (no gold points on the volcanoes)",0,INK,False),
    ("Interpretation: these are EMBRYONIC patterning genes; by this postnatal stage they no longer distinguish the maculae",0,INK,False),
])
footer(s,5)

# ---- 6 volcano technical signature
s=slide(); header(s,"Key result 2","DE hits are dominated by a technical signature")
pic(s,de("volcano_HC_type_I.png"),Inches(0.3),Inches(1.35),Inches(7.2),Inches(5.2),
    cap="HC_type_I (hair cells): utricle vs saccule")
bullets(s,Inches(7.7),Inches(1.5),Inches(5.3),Inches(5.0),[
    ("The 'up in saccule' side is overwhelmingly mitochondrial (mt-Co1, mt-Nd1, mt-Cytb…) and ribosomal (Rplp1, Rps21, Rps24) genes",0,ACC2,True),
    ("Plus ambient/low-specificity genes (Slco1c1, Apod)",1,GREY,False),
    ("This is a classic QUALITY / dissociation / batch signature between the two tissue preps — not cell-intrinsic biology",0,INK,False),
    ("Same pattern in Supporting_cell and Hair_cell volcanoes",0,GREY,False),
    ("→ the large 'saccule-up' counts are inflated by these technical genes",0,INK,True),
])
footer(s,6)

# ---- 7 pseudobulk near-identity (headline)
s=slide(); header(s,"Key result 3 (headline)","At the replicate level, the two organs are nearly identical")
pic(s,de("pseudobulk_scatter_HC_type_I.png"),Inches(0.3),Inches(1.35),Inches(5.9),Inches(5.2),
    cap="HC_type_I pseudobulk (mean of replicates)")
pic(s,de("pseudobulk_scatter_Supporting_cell.png"),Inches(6.3),Inches(1.35),Inches(5.9),Inches(5.2),
    cap="Supporting cells pseudobulk")
bullets(s,Inches(0.5),Inches(6.55),Inches(12.4),Inches(0.8),[
    ("Averaging expression across the 2 replicates per tissue, utricle vs saccule falls almost perfectly on the diagonal (y = x) → no robust transcriptional difference; the per-cell 'significance' was driven by cell number + technical variance",0,GREEN,True),
])
footer(s,7)

# ---- 8 candidates + caveats
s=slide(); header(s,"Candidates & caveats","What (if anything) might be real")
bullets(s,Inches(0.5),Inches(1.4),Inches(6.2),Inches(5.3),[
    ("A few non-technical candidates survive (follow-up, not conclusions):",0,INK,True),
    ("↑ saccule: Otos (otospiralin), Coch (cochlin), Corin, Slc4a11, Ecrg4, Acan",1,INK,False),
    ("↑ utricle: Ranbp3l, Rspo3, Kcnq1, Enpep, Dct, Stac",1,INK,False),
    ("Otos and Coch (secreted matrix / otoconial) are plausibly real regional differences",1,GREY,False),
])
bullets(s,Inches(6.9),Inches(1.4),Inches(6.1),Inches(5.3),[
    ("Why to stay cautious:",0,ACC2,True),
    ("mito/ribosomal/ambient genes contaminate the saccule-up side (quality/batch)",1,INK,False),
    ("only 2 replicates per tissue; Wilcoxon on cells inflates p-values",1,INK,False),
    ("pseudobulk shows near-identity (slide 7)",1,INK,False),
    ("To settle it: pseudobulk DESeq2 (2v2), regress out mito/ribo, validate Otos/Coch by ISH/qPCR",0,INK,True),
])
footer(s,8)

# ---- 9 summary
s=slide(); header(s,"Summary","Saccule vs utricle — a shared molecular program")
bullets(s,Inches(0.5),Inches(1.5),Inches(12.3),Inches(5.3),[
    ("Both organs were integrated and annotated independently and recover the SAME cell types in similar proportions",0,INK,False),
    ("The classic regional patterning genes (Otx1/Otx2/Emx2/Tbx1) do not distinguish the organs postnatally",0,INK,False),
    ("Per-cell DE returns many 'hits', but they are dominated by mitochondrial/ribosomal technical genes, and pseudobulk (replicate-level) expression is nearly identical",0,INK,False),
    ("→ No robust, cell-intrinsic saccule-vs-utricle transcriptional signature at this depth",0,GREEN,True),
    ("A short list of candidates (Otos, Coch, Corin, Rspo3, Kcnq1) is worth a targeted, properly-powered follow-up",0,INK,False),
    ("",0,INK,False),
    ("Next steps: pseudobulk DESeq2 with mito/ribo regressed · ISH/qPCR validation of top candidates · consider that the maculae may simply be transcriptionally very similar",0,ACC2,True),
])
footer(s,9)

prs.save(OUT); print("SAVED:",OUT,"| slides:",len(prs.slides._sldIdLst))
