#!/usr/bin/env python
"""Report deck: per-organ (saccule / utricle) integration + marker annotation."""
import os, sys
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE
from PIL import Image

FIG = "../../analyses/07_single_cell_rna/results/per_tissue"
OUT = sys.argv[1] if len(sys.argv) > 1 else \
    "../../analyses/07_single_cell_rna/report/per_organ_annotation_report.pptx"

INK=RGBColor(0x1a,0x1a,0x2e); ACC=RGBColor(0x2e,0x86,0xab); ACC2=RGBColor(0xc0,0x39,0x2b)
GREY=RGBColor(0x55,0x5b,0x6e); LIGHT=RGBColor(0xf4,0xf6,0xf8); WHITE=RGBColor(0xff,0xff,0xff)
SAC=RGBColor(0xc0,0x39,0x2b); UTR=RGBColor(0x2e,0x86,0xab)

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
def header(s,k,t,kc=ACC):
    rect(s,0,0,SW,Inches(1.12),INK); rect(s,0,Inches(1.12),SW,Inches(0.06),kc)
    _,tf=tbox(s,Inches(0.5),Inches(0.10),Inches(12.3),Inches(0.98),MSO_ANCHOR.MIDDLE)
    run(tf.paragraphs[0],k.upper(),12,kc,bold=True); run(tf.add_paragraph(),t,25,WHITE,bold=True)
def bullets(s,x,y,w,h,items):
    _,tf=tbox(s,x,y,w,h)
    for i,it in enumerate(items):
        text,lvl,c,bold=(it if isinstance(it,tuple) else (it,0,INK,False))
        p=tf.paragraphs[0] if i==0 else tf.add_paragraph(); p.space_after=Pt(6)
        run(p,("   – " if lvl else "• "),13 if lvl==0 else 11.5,ACC if lvl==0 else GREY,bold=True)
        run(p,text,13 if lvl==0 else 11.5,c,bold=bold)
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
def table(s,rows,x,y,w,h,hdr=INK):
    t=s.shapes.add_table(len(rows),len(rows[0]),x,y,w,h).table
    for r in range(len(rows)):
        for c in range(len(rows[0])):
            cell=t.cell(r,c); cell.vertical_anchor=MSO_ANCHOR.MIDDLE
            cell.margin_top=cell.margin_bottom=Inches(0.015); cell.margin_left=Inches(0.06)
            pp=cell.text_frame.paragraphs[0]
            if c: pp.alignment=PP_ALIGN.CENTER
            h0=(r==0); run(pp,rows[r][c],10.5 if h0 else 10,WHITE if h0 else INK,bold=(h0 or c==0))
            cell.fill.solid(); cell.fill.fore_color.rgb=hdr if h0 else (LIGHT if r%2 else WHITE)
    return t
def footer(s,n):
    _,tf=tbox(s,Inches(0.4),SH-Inches(0.4),Inches(11),Inches(0.3))
    run(tf.paragraphs[0],"Per-organ integration & marker annotation — saccule vs utricle",9,GREY)
    _,tf2=tbox(s,SW-Inches(1.1),SH-Inches(0.4),Inches(0.7),Inches(0.3)); pp=tf2.paragraphs[0]
    pp.alignment=PP_ALIGN.RIGHT; run(pp,str(n),10,GREY)
def F(n): return os.path.join(FIG,n)

# 1 title
s=slide(); rect(s,0,0,SW,SH,INK); rect(s,0,Inches(4.2),SW,Inches(0.08),ACC)
_,tf=tbox(s,Inches(0.9),Inches(1.6),Inches(11.5),Inches(2.4))
run(tf.paragraphs[0],"Per-organ integration & annotation",37,WHITE,bold=True)
run(tf.add_paragraph(),"Saccule and utricle integrated & annotated separately",19,ACC)
_,tf=tbox(s,Inches(0.9),Inches(4.5),Inches(11.5),Inches(1.8))
for t,sz in [("Each organ = its own 2 replicates, integrated on its own (not pooled with the other organ)",15),
             ("Saccule: FTL_1146 + 1148 (25,172 cells) · Utricle: FTL_1147 + 1149 (24,037 cells)",13),
             ("CCA integration → clustering → 14-panel marker annotation",12)]:
    run(tf.add_paragraph(),t,sz,WHITE if sz>13 else GREY)
run(tbox(s,Inches(0.9),Inches(6.6),Inches(11),Inches(0.4))[1].paragraphs[0],"Lab meeting · Mehdi Esmaeili Fard",12,GREY)

# 2 approach
s=slide(); header(s,"Approach","Integrate each organ on its own, then annotate")
bullets(s,Inches(0.55),Inches(1.5),Inches(12.2),Inches(5.4),[
    ("For EACH organ separately (this replaces the all-4-sample integration):",0,INK,True),
    ("Read the 2 replicates → QC (mt- < 5%, nCount > 800, nFeature > 500) → normalize",1,GREY,False),
    ("CCA integration of the 2 replicates (batch-correct within the organ)",1,GREY,False),
    ("Cluster (res 0.3) → UMAP → annotate each cluster by 14 canonical marker panels (module-score argmax)",1,GREY,False),
    ("FindAllMarkers per cell type (marker evidence)",1,GREY,False),
    ("Why per-organ? a clean, organ-specific map and composition, without the other organ influencing clustering",0,INK,False),
    ("Same marker panels as before (hair / supporting / transitional / dark / glia / neuron / immune / endothelial / melanocyte / fibrocyte / pericyte)",0,GREY,False),
])
footer(s,2)

# 3 saccule UMAP + counts
s=slide(); header(s,"Saccule — annotation","25,172 cells · 9 cell types",kc=SAC)
pic(s,F("saccule_UMAP_celltype.png"),Inches(0.3),Inches(1.35),Inches(7.3),Inches(5.3))
table(s,[["cell type","cells"],["HC_type_I *","13,287"],["Supporting_cell","5,211"],
         ["Fibrocyte_mesench","1,926"],["Transitional_epith","1,423"],["Hair_cell *","1,302"],
         ["Glia_Schwann","1,226"],["Roof_dark_cell","570"],["Immune_macrophage","167"],
         ["Melanocyte","60"]],
      Inches(7.8),Inches(1.4),Inches(5.1),Inches(4.8),hdr=SAC)
run(tbox(s,Inches(7.8),Inches(6.35),Inches(5.1),Inches(0.7))[1].paragraphs[0],
    "* HC_type_I + Hair_cell are both the hair-cell compartment (~14,600) — the I/II label is not reliable (see slide 7)",10,ACC2,italic=True)
footer(s,3)

# 4 saccule quality + markers
s=slide(); header(s,"Saccule — integration & markers","Replicates mix well; canonical markers are clean",kc=SAC)
pic(s,F("saccule_UMAP_replicate.png"),Inches(0.3),Inches(1.35),Inches(6.3),Inches(5.3),
    cap="rep1/rep2 overlap within clusters → integration OK")
pic(s,F("saccule_DotPlot_markers.png"),Inches(6.7),Inches(1.35),Inches(6.4),Inches(5.3),
    cap="One representative marker per cell type")
footer(s,4)

# 5 utricle UMAP + counts
s=slide(); header(s,"Utricle — annotation","24,037 cells · 8 cell types",kc=UTR)
pic(s,F("utricle_UMAP_celltype.png"),Inches(0.3),Inches(1.35),Inches(7.3),Inches(5.3))
table(s,[["cell type","cells"],["Hair_cell *","8,062"],["Supporting_cell","6,703"],
         ["Transitional_epith","3,601"],["Fibrocyte_mesench","2,595"],["Glia_Schwann","1,537"],
         ["Roof_dark_cell","1,157"],["HC_type_II *","272"],["Immune_macrophage","110"]],
      Inches(7.8),Inches(1.4),Inches(5.1),Inches(4.4),hdr=UTR)
run(tbox(s,Inches(7.8),Inches(6.0),Inches(5.1),Inches(1.0))[1].paragraphs[0],
    "* Hair_cell + HC_type_II are the hair-cell compartment (~8,300). No Melanocyte cluster resolved here (present in saccule) — likely a resolution/rarity effect, not a real absence.",10,ACC2,italic=True)
footer(s,5)

# 6 utricle quality + markers
s=slide(); header(s,"Utricle — integration & markers","Replicates mix well; markers clean",kc=UTR)
pic(s,F("utricle_UMAP_replicate.png"),Inches(0.3),Inches(1.35),Inches(6.3),Inches(5.3),
    cap="rep1/rep2 overlap within clusters")
pic(s,F("utricle_DotPlot_markers.png"),Inches(6.7),Inches(1.35),Inches(6.4),Inches(5.3),
    cap="Representative markers by cell type")
footer(s,6)

# 7 comparison + caveats
s=slide(); header(s,"Saccule vs utricle","Same cell types; differences are labelling/batch, not biology")
pic(s,F("saccule_composition_bar.png"),Inches(0.3),Inches(1.4),Inches(4.3),Inches(4.2),
    cap="Saccule composition")
pic(s,F("utricle_composition_bar.png"),Inches(4.7),Inches(1.4),Inches(4.3),Inches(4.2),
    cap="Utricle composition")
bullets(s,Inches(9.2),Inches(1.45),Inches(3.9),Inches(5.4),[
    ("Both organs recover the same core types: hair, supporting, transitional, dark/roof, glia, fibrocyte, immune",0,INK,False),
    ("Hair-cell label differs (saccule→HC_type_I, utricle→HC_type_II) — an ARGMAX artifact; I/II is not resolvable on 3′ data. Both are just hair cells.",0,ACC2,True),
    ("rep2 (FTL_1148/1149) is enriched in supporting / dark / transitional cells in BOTH organs → a batch/dissection effect",0,INK,False),
    ("Per-cell-type, the two organs are transcriptionally near-identical (edgeR: 0 DE genes)",0,GREY,False),
])
footer(s,7)

# 8 summary
s=slide(); header(s,"Summary","")
bullets(s,Inches(0.55),Inches(1.6),Inches(12.2),Inches(5.0),[
    ("Each organ was integrated from its own 2 replicates and annotated independently",0,INK,False),
    ("Saccule: 9 types / 25,172 cells; Utricle: 8 types / 24,037 cells — replicates integrate well, markers are clean",0,INK,False),
    ("Both organs share the same cell-type inventory; apparent differences are:",0,INK,True),
    ("hair-cell subtype labels (HC_type_I vs HC_type_II) — an argmax artifact, not biology",1,ACC2,False),
    ("a rare Melanocyte cluster resolved only in saccule (resolution/rarity)",1,GREY,False),
    ("rep2 enrichment of non-sensory cells — a batch/dissection effect (rep1 vs rep2)",1,GREY,False),
    ("Consistent with the DE analysis: within matched cell types the organs are transcriptionally very similar",0,INK,False),
    ("Saved objects: results/per_tissue/{saccule,utricle}_integrated_annotated.rds",0,GREY,False),
])
footer(s,8)

prs.save(OUT); print("SAVED:",OUT,"| slides:",len(prs.slides._sldIdLst))
