#!/usr/bin/env python
"""Deck: saccule-vs-utricle pseudobulk PCA + edgeR DE (protein-coding + lncRNA)."""
import os, sys
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE
from PIL import Image

FIG = "../../analyses/07_single_cell_rna/results/pseudobulk_DE"
OUT = sys.argv[1] if len(sys.argv) > 1 else \
    "../../analyses/07_single_cell_rna/report/saccule_vs_utricle_pseudobulk_DE.pptx"

INK=RGBColor(0x1a,0x1a,0x2e); ACC=RGBColor(0x2e,0x86,0xab); ACC2=RGBColor(0xc0,0x39,0x2b)
GREY=RGBColor(0x55,0x5b,0x6e); LIGHT=RGBColor(0xf4,0xf6,0xf8); WHITE=RGBColor(0xff,0xff,0xff)
GREEN=RGBColor(0x2e,0x7d,0x32)

prs=Presentation(); prs.slide_width=Inches(13.333); prs.slide_height=Inches(7.5)
SW,SH=prs.slide_width,prs.slide_height; BLANK=prs.slide_layouts[6]
def slide(): return prs.slides.add_slide(BLANK)
def rect(s,x,y,w,h,c):
    sh=s.shapes.add_shape(MSO_SHAPE.RECTANGLE,x,y,w,h); sh.fill.solid()
    sh.fill.fore_color.rgb=c; sh.line.fill.background(); sh.shadow.inherit=False; return sh
def tbox(s,x,y,w,h,a=MSO_ANCHOR.TOP):
    tb=s.shapes.add_textbox(x,y,w,h); tf=tb.text_frame; tf.word_wrap=True; tf.vertical_anchor=a
    tf.margin_left=tf.margin_right=Inches(0.05); tf.margin_top=tf.margin_bottom=Inches(0.02); return tb,tf
def run(p,t,sz,c,bold=False,italic=False):
    r=p.add_run(); r.text=t; f=r.font; f.size=Pt(sz); f.color.rgb=c; f.bold=bold; f.italic=italic; f.name="Calibri"; return r
def header(s,k,t):
    rect(s,0,0,SW,Inches(1.12),INK); rect(s,0,Inches(1.12),SW,Inches(0.06),ACC)
    _,tf=tbox(s,Inches(0.5),Inches(0.10),Inches(12.3),Inches(0.98),MSO_ANCHOR.MIDDLE)
    run(tf.paragraphs[0],k.upper(),12,ACC,bold=True); run(tf.add_paragraph(),t,25,WHITE,bold=True)
def bullets(s,x,y,w,h,items):
    _,tf=tbox(s,x,y,w,h)
    for i,it in enumerate(items):
        text,lvl,c,bold=(it if isinstance(it,tuple) else (it,0,INK,False))
        p=tf.paragraphs[0] if i==0 else tf.add_paragraph(); p.space_after=Pt(6)
        run(p,("   – " if lvl else "• "),13.5 if lvl==0 else 12,ACC if lvl==0 else GREY,bold=True)
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
    run(tf.paragraphs[0],"Saccule vs utricle — pseudobulk PCA & edgeR DE (protein-coding + lncRNA)",9,GREY)
    _,tf2=tbox(s,SW-Inches(1.1),SH-Inches(0.4),Inches(0.7),Inches(0.3)); pp=tf2.paragraphs[0]
    pp.alignment=PP_ALIGN.RIGHT; run(pp,str(n),10,GREY)
def F(n): return os.path.join(FIG,n)

# 1 title
s=slide(); rect(s,0,0,SW,SH,INK); rect(s,0,Inches(4.2),SW,Inches(0.08),ACC)
_,tf=tbox(s,Inches(0.9),Inches(1.5),Inches(11.5),Inches(2.6))
run(tf.paragraphs[0],"Saccule vs utricle: PCA & proper (edgeR) DE",35,WHITE,bold=True)
run(tf.add_paragraph(),"Pseudobulk analysis, protein-coding + lncRNA genes only",19,ACC)
_,tf=tbox(s,Inches(0.9),Inches(4.5),Inches(11.5),Inches(1.9))
for t,sz in [("Reconciling the single-cell data with the bulk PCA separation",15),
             ("4 samples · 2 organs × 2 replicates · all female · edgeR pseudobulk (Env1)",13),
             ("Run on the cluster (SLURM); scripts annot7 (prep) + annot8 (DE)",12)]:
    run(tf.add_paragraph(),t,sz,WHITE if sz>13 else GREY)
run(tbox(s,Inches(0.9),Inches(6.6),Inches(11),Inches(0.4))[1].paragraphs[0],"Lab meeting · Mehdi Esmaeili Fard",12,GREY)

# 2 what/why
s=slide(); header(s,"Motivation & approach","Bulk shows tissue separation — does single-cell?")
bullets(s,Inches(0.55),Inches(1.5),Inches(12.2),Inches(5.4),[
    ("Observation: the BULK RNA-seq PCA cleanly separates saccule from utricle — we expected the same here",0,INK,True),
    ("Earlier single-cell DE (Wilcoxon on cells) was unreliable — cells are not replicates, and hits were mito/ribosomal",0,ACC2,False),
    ("So we redid it properly:",0,INK,True),
    ("Keep only protein-coding + lncRNA genes (remove pseudogenes / miRNA / IG-TR)",1,GREY,False),
    ("Build PSEUDOBULK (sum counts per sample, and per cell-type × sample) — replicate-aware",1,GREY,False),
    ("PCA of the samples; edgeR DE utricle vs saccule at organ level AND per cell type",1,GREY,False),
    ("edgeR/DESeq2 available in Env1; pseudobulk built from the combined 11-cell-type object (matched labels)",0,GREY,False),
])
footer(s,2)

# 3 gene filter
s=slide(); header(s,"Gene filtering","protein-coding + lncRNA")
bullets(s,Inches(0.55),Inches(1.5),Inches(12.2),Inches(5.4),[
    ("Biotypes from the 10x GRCm39 GTF: kept gene_type ∈ {protein_coding, lncRNA}",0,INK,False),
    ("25,203 of 25,234 detected genes retained (removed: pseudogenes, IG/TR genes, miRNA, rRNA, snoRNA, TEC)",1,GREY,False),
    ("Note: GENCODE calls it 'lncRNA' (subsumes the old 'lincRNA')",1,GREY,False),
    ("Caveat — mitochondrial and ribosomal-protein genes are protein_coding, so they SURVIVE this filter",0,ACC2,True),
    ("i.e. the filter cleans pseudogenes/predicted genes but not the mito/ribo technical signal",1,GREY,False),
    ("(An optional mito/ribosomal blocklist can be added if wanted)",1,GREY,False),
])
footer(s,3)

# 4 PCA raw — batch dominates
s=slide(); header(s,"Result 1 — PCA (raw)","The dominant axis is BATCH, not tissue")
pic(s,F("PCA_organ_pseudobulk.png"),Inches(0.3),Inches(1.35),Inches(7.4),Inches(5.3))
bullets(s,Inches(7.9),Inches(1.6),Inches(5.1),Inches(5.0),[
    ("PC1 (54%): the two REPLICATES separate — rep1 (1146+1147) left, rep2 (1148+1149) right",0,ACC2,True),
    ("Same-batch saccule + utricle cluster together; tissue does NOT separate on PC1 or PC2",0,INK,False),
    ("→ a strong batch effect between the two 10x runs dominates the variation",0,INK,True),
    ("This is why a naive analysis shows no tissue signal — and a clue to the bulk vs single-cell difference",0,GREY,False),
])
footer(s,4)

# 5 PCA batch-corrected — tissue separates
s=slide(); header(s,"Result 2 — PCA (batch-corrected)","Once batch is removed, the organs DO separate")
pic(s,F("PCA_organ_batch_corrected.png"),Inches(0.3),Inches(1.35),Inches(7.4),Inches(5.3))
bullets(s,Inches(7.9),Inches(1.6),Inches(5.1),Inches(5.0),[
    ("After limma removeBatchEffect(replicate):",0,INK,False),
    ("PC2 separates saccule (both positive) from utricle (both negative)",0,GREEN,True),
    ("So a consistent tissue signal IS present — matching your bulk PCA",0,INK,False),
    ("but it is a MODEST secondary axis (18.8%), previously masked by the batch effect",0,INK,False),
    ("Reconciliation: the tissue difference is real but small; the batch was hiding it",0,ACC,True),
])
footer(s,5)

# 6 organ DE
s=slide(); header(s,"Result 3 — organ-level DE (edgeR, batch-blocked)","No gene reaches significance")
pic(s,F("volcano_organ.png"),Inches(0.3),Inches(1.35),Inches(7.2),Inches(5.3))
bullets(s,Inches(7.7),Inches(1.6),Inches(5.3),Inches(5.0),[
    ("Design ~ replicate + tissue (paired); 0 genes at FDR<0.05, |log2FC|>1",0,ACC2,True),
    ("With only 2 replicates per organ, power is very low — a real but small effect won't reach significance",0,INK,False),
    ("Top nominal candidates (not significant): Acan, Bglap (matrix / otoconial), Vcam1, Hes1 (Notch)",0,INK,False),
    ("Y-chromosome genes appeared in the noise — but all 4 samples are FEMALE (Xist+, Y≈0), so no sex confound",0,GREY,False),
])
footer(s,6)

# 7 cell-level DE
s=slide(); header(s,"Result 4 — cell-level DE (edgeR, per cell type)","Zero significant genes in every cell type")
pic(s,F("cellDE_summary_bar.png"),Inches(0.3),Inches(1.4),Inches(7.0),Inches(5.2),
    cap="Significant genes per cell type (all zero)")
bullets(s,Inches(7.5),Inches(1.6),Inches(5.4),Inches(5.0),[
    ("Pseudobulk 2v2 within each of 9 testable cell types, ~replicate + tissue",0,INK,False),
    ("Every cell type: 0 significant genes (FDR<0.05, |log2FC|>1)",0,ACC2,True),
    ("Consistent with the pseudobulk near-identity seen earlier — within a matched cell type, the two organs are transcriptionally almost the same",0,INK,False),
    ("Again a power issue at n=2, but the effect sizes are genuinely small",0,GREY,False),
])
footer(s,7)

# 8 reconciliation
s=slide(); header(s,"Why bulk separates but single-cell 'doesn't'","Reconciling the two results")
bullets(s,Inches(0.55),Inches(1.5),Inches(12.2),Inches(5.4),[
    ("The tissue difference IS there — it shows up on PC2 after batch correction",0,GREEN,True),
    ("But three things weaken it in this single-cell dataset:",0,INK,True),
    ("Batch effect: a large rep1-vs-rep2 technical difference dominates PC1 (54%) and inflates dispersion",1,INK,False),
    ("Power: only 2 biological replicates per organ → edgeR can't reach FDR significance for a small effect",1,INK,False),
    ("Composition-controlled: per-cell-type pseudobulk removes the cell-mix differences that bulk 'sees'",1,INK,False),
    ("Your bulk PCA likely had more replicates and/or less batch structure, so the same modest signal separated cleanly",0,INK,False),
    ("Net: not a contradiction — same direction, but underpowered and batch-confounded here",0,ACC,True),
])
footer(s,8)

# 9 summary
s=slide(); header(s,"Summary & recommendations","")
bullets(s,Inches(0.55),Inches(1.5),Inches(12.2),Inches(3.5),[
    ("Saccule and utricle share cell types and are transcriptionally very similar per cell type",0,INK,False),
    ("A consistent tissue signal exists (PC2 after batch correction) but is modest",0,INK,False),
    ("No gene passes FDR at organ or cell level with 2 replicates + a dominant batch effect",0,INK,False),
    ("All samples female → no sex confound; nominal candidates point to matrix/otoconial genes (Acan, Bglap)",0,GREY,False),
])
rect(s,Inches(0.55),Inches(5.05),Inches(12.2),Inches(0.045),ACC2)
bullets(s,Inches(0.55),Inches(5.2),Inches(12.2),Inches(1.9),[
    ("To detect the tissue difference cleanly",0,ACC2,True),
    ("Add biological replicates (≥3–4 per organ) — the single biggest fix for power",1,INK,False),
    ("Control the batch (balanced design + ~batch in the model, as done here)",1,INK,False),
    ("Validate the top candidates (Acan/Bglap/otoconial genes) directly by ISH/qPCR",1,INK,False),
])
footer(s,9)

prs.save(OUT); print("SAVED:",OUT,"| slides:",len(prs.slides._sldIdLst))
