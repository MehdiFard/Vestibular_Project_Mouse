#!/usr/bin/env python
"""Consolidated full report: annotation (1 sample -> per-organ -> all 4) + PCA/DE."""
import os, sys
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE
from PIL import Image

RES = "../../analyses/07_single_cell_rna/results"
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(RES), "report", "vestibular_FULL_report.pptx")

INK=RGBColor(0x1a,0x1a,0x2e); ACC=RGBColor(0x2e,0x86,0xab); ACC2=RGBColor(0xc0,0x39,0x2b)
GREY=RGBColor(0x55,0x5b,0x6e); LIGHT=RGBColor(0xf4,0xf6,0xf8); WHITE=RGBColor(0xff,0xff,0xff)
GREEN=RGBColor(0x2e,0x7d,0x32); SAC=ACC2; UTR=ACC

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
    run(tf.paragraphs[0],k.upper(),12,kc,bold=True); run(tf.add_paragraph(),t,24,WHITE,bold=True)
def bullets(s,x,y,w,h,items):
    _,tf=tbox(s,x,y,w,h)
    for i,it in enumerate(items):
        text,lvl,c,bold=(it if isinstance(it,tuple) else (it,0,INK,False))
        p=tf.paragraphs[0] if i==0 else tf.add_paragraph(); p.space_after=Pt(5)
        run(p,("   – " if lvl else "• "),12.5 if lvl==0 else 11,ACC if lvl==0 else GREY,bold=True)
        run(p,text,12.5 if lvl==0 else 11,c,bold=bold)
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
        pp.alignment=PP_ALIGN.CENTER; run(pp,cap,9.5,GREY,italic=True)
def footer(s,n):
    _,tf=tbox(s,Inches(0.4),SH-Inches(0.38),Inches(11),Inches(0.3))
    run(tf.paragraphs[0],"Mouse vestibular maculae — full scRNA-seq report",9,GREY)
    _,tf2=tbox(s,SW-Inches(1.0),SH-Inches(0.38),Inches(0.7),Inches(0.3)); pp=tf2.paragraphs[0]
    pp.alignment=PP_ALIGN.RIGHT; run(pp,str(n),10,GREY)
def CA(n): return os.path.join(RES,"cell_annotation",n)
def PT(n): return os.path.join(RES,"per_tissue",n)
def PB(n): return os.path.join(RES,"pseudobulk_DE",n)
N=[0]
def num(): N[0]+=1; return N[0]

# 1 title
s=slide(); rect(s,0,0,SW,SH,INK); rect(s,0,Inches(4.3),SW,Inches(0.08),ACC)
_,tf=tbox(s,Inches(0.9),Inches(1.5),Inches(11.6),Inches(2.6))
run(tf.paragraphs[0],"Cell-type annotation & organ comparison",34,WHITE,bold=True)
run(tf.add_paragraph(),"of the mouse vestibular maculae (utricle & saccule)",26,WHITE,bold=True)
run(tf.add_paragraph(),"Single-cell RNA-seq — a full report",18,ACC)
_,tf=tbox(s,Inches(0.9),Inches(4.7),Inches(11.6),Inches(1.7))
for t,sz in [("From one sample → per-organ → all samples → PCA & differential expression",15),
             ("4 samples · 2 organs × 2 replicates · 10x 3′ v4 · Mus musculus (GRCm39) · ~49,000 cells",13)]:
    run(tf.add_paragraph(),t,sz,WHITE if sz>13 else GREY)
run(tbox(s,Inches(0.9),Inches(6.6),Inches(11),Inches(0.4))[1].paragraphs[0],"Lab meeting · Mehdi Esmaeili Fard",12,GREY)

# 2 dataset
s=slide(); header(s,"Dataset & question","Utricle vs saccule, two replicates each")
bullets(s,Inches(0.5),Inches(1.4),Inches(6.4),Inches(5.2),[
    ("10x 3′ v4 single-cell RNA-seq of the mouse vestibular maculae",0,INK,False),
    ("QC: nCount>800, nFeature>500, mt- <5% (mouse mito prefix)",1,GREY,False),
    ("~49,000 cells after QC; all samples female",1,GREY,False),
    ("Two goals:",0,INK,True),
    ("(1) a trustworthy cell-type map of the sensory maculae",1,INK,False),
    ("(2) compare molecular profile of utricle vs saccule",1,INK,False),
    ("Challenge: specialized tissue — generic references lack hair/supporting cells",0,ACC2,False),
])
gx,gy,cw,ch=Inches(7.15),Inches(1.7),Inches(1.87),Inches(1.25)
grid=[["","Rep 1","Rep 2"],["Utricle","FTL_1147","FTL_1149"],["Saccule","FTL_1146","FTL_1148"]]
for r in range(3):
    for c in range(3):
        hd=(r==0 or c==0)
        b=rect(s,gx+c*cw,gy+r*ch,cw-Inches(0.05),ch-Inches(0.05),ACC if hd else LIGHT)
        tf=b.text_frame; tf.vertical_anchor=MSO_ANCHOR.MIDDLE; pp=tf.paragraphs[0]; pp.alignment=PP_ALIGN.CENTER
        run(pp,grid[r][c],13,WHITE if hd else INK,bold=hd)
footer(s,num())

# 3 roadmap
s=slide(); header(s,"Roadmap","How the analysis built up")
steps=[("1","One utricle sample","marker-based annotation of FTL_1147 (proof of concept)"),
       ("2","Validation","independent label transfer (Burns 2015) + automatic bulk-SingleR"),
       ("3","Per-organ integration","saccule and utricle each integrated from its 2 replicates"),
       ("4","All-samples integration","all 4 co-embedded → one 11-type map + composition"),
       ("5","Organ comparison","pseudobulk PCA + edgeR DE: does utricle differ from saccule?")]
y=Inches(1.5)
for i,(nn,t,d) in enumerate(steps):
    yy=y+i*Inches(1.02)
    c=rect(s,Inches(0.6),yy,Inches(0.7),Inches(0.8),INK)
    tf=c.text_frame; tf.vertical_anchor=MSO_ANCHOR.MIDDLE; pp=tf.paragraphs[0]; pp.alignment=PP_ALIGN.CENTER
    run(pp,nn,22,WHITE,bold=True)
    _,tf=tbox(s,Inches(1.5),yy,Inches(11.2),Inches(0.8),MSO_ANCHOR.MIDDLE)
    run(tf.paragraphs[0],t+"  —  ",15,ACC,bold=True); run(tf.paragraphs[0],d,13.5,INK)
footer(s,num())

# 4 marker method (Part 1)
s=slide(); header(s,"Part 1 · annotation method","Marker-based cell typing")
bullets(s,Inches(0.5),Inches(1.4),Inches(6.4),Inches(5.3),[
    ("Standard Seurat pipeline: QC → normalize → HVG → PCA → cluster → UMAP",0,INK,False),
    ("Score 14 canonical marker panels (AddModuleScore); assign each cluster by argmax; verify with DE genes",0,INK,False),
    ("Panels from the inner-ear literature + standard lineages:",0,INK,True),
    ("Sensory: hair cells (Myo7a, Pou4f3), Type I (Spp1, Ocm), Type II (Anxa4, Mapt), supporting (Sox2, Gjb2), transitional epith. (Krt8), dark/roof (Kcnq1)",1,GREY,False),
    ("Non-sensory: glia (Plp1), neuron (Snap25), immune (Ptprc), endothelial (Pecam1), melanocyte (Dct/Mlana), fibrocyte (Col1a1), pericyte (Acta2)",1,GREY,False),
    ("Sources: Burns et al. 2015; McInturff et al. 2018; CellMarker / PanglaoDB",0,GREY,False),
])
pic(s,CA("DotPlot_markers_FTL_1147.png"),Inches(7.05),Inches(1.4),Inches(6.0),Inches(5.2),
    cap="Marker DotPlot validates the panels (FTL_1147)")
footer(s,num())

# 5 single utricle sample
s=slide(); header(s,"Part 1 · step 1","One utricle sample (FTL_1147) — a clean first map")
pic(s,CA("UMAP_markers_refined_FTL_1147.png"),Inches(0.4),Inches(1.35),Inches(7.4),Inches(5.4))
bullets(s,Inches(7.9),Inches(1.6),Inches(5.1),Inches(5.0),[
    ("~8,650 cells; six well-separated types",0,INK,True),
    ("Hair cells ~3,700 · Supporting ~3,440",1,INK,False),
    ("Glia 663 · Endothelial 662 · Pericyte 184",1,INK,False),
    ("Sensory epithelium dominates — expected for a macula dissection",0,INK,False),
    ("Every call backed by DE genes + module-score heatmap",0,GREY,False),
])
footer(s,num())

# 6 HC I/II
s=slide(); header(s,"Part 1 · caveat","Hair-cell Type I vs II — a limit of 3′ data")
pic(s,CA("UMAP_HC_subtype_FTL_1147.png"),Inches(0.4),Inches(1.35),Inches(6.3),Inches(5.4))
bullets(s,Inches(6.9),Inches(1.6),Inches(6.1),Inches(5.0),[
    ("Sub-clustered hair cells; scored Type I (Spp1/Ocm) vs Type II (Anxa4/Mapt)",0,INK,False),
    ("Confident-only: ~1,110 Type II; the rest 'unresolved'; NO confident Type I",0,ACC2,True),
    ("Why: 10x 3′ captures Ocm/Spp1 too sparsely to separate the subtypes",0,INK,False),
    ("Consequence: treat hair cells as ONE compartment; the I/II labels elsewhere are argmax artifacts",0,INK,True),
])
footer(s,num())

# 7 validation
s=slide(); header(s,"Part 1 · step 2","Validation — independent reference agrees")
pic(s,CA("UMAP_marker_vs_burns_FTL_1147.png"),Inches(0.4),Inches(1.4),Inches(8.7),Inches(4.4),
    cap="Left: our marker labels.  Right: Burns 2015 labels via SingleR")
bullets(s,Inches(9.3),Inches(1.5),Inches(3.7),Inches(5.2),[
    ("Label transfer from Burns 2015 (mouse utricle) via SingleR + Seurat",0,INK,False),
    ("86.9% agreement on epithelial cells (HC→HC 0.86, SC→SC 0.88)",0,GREEN,True),
    ("Independent confirmation of the marker calls",0,INK,False),
    ("(An earlier automatic bulk-SingleR pass named only broad lineages — no hair/supporting cells)",0,GREY,False),
])
footer(s,num())

# 8 per-organ
s=slide(); header(s,"Part 1 · step 3","Per-organ integration — each organ on its own")
pic(s,PT("saccule_UMAP_celltype.png"),Inches(0.3),Inches(1.4),Inches(6.4),Inches(4.7),
    cap="Saccule — 25,172 cells, 9 types")
pic(s,PT("utricle_UMAP_celltype.png"),Inches(6.7),Inches(1.4),Inches(6.4),Inches(4.7),
    cap="Utricle — 24,037 cells, 8 types")
bullets(s,Inches(0.5),Inches(6.15),Inches(12.4),Inches(1.2),[
    ("Each organ integrated from its 2 replicates (CCA); replicates mix well; same core cell types recovered",0,INK,False),
    ("Hair-cell label differs (saccule→HC_type_I, utricle→HC_type_II) = argmax artifact, not biology — both are hair cells",0,ACC2,False),
])
footer(s,num())

# 9 all-4 integration
s=slide(); header(s,"Part 1 · step 4","All samples integrated — one map, reproducible")
pic(s,CA("UMAP_integrated_celltype_FTL_all.png"),Inches(0.35),Inches(1.4),Inches(7.0),Inches(5.4))
pic(s,CA("UMAP_integrated_tissue_sample_FTL_all.png"),Inches(7.45),Inches(1.5),Inches(5.5),Inches(2.7),
    cap="Cells co-embed by biology, not sample → integrated")
bullets(s,Inches(7.55),Inches(4.35),Inches(5.4),Inches(2.4),[
    ("CCA integration of all 4 samples → 11 cell types",0,INK,False),
    ("All 11 types present in BOTH replicates of BOTH organs → reproducible, not batch artifacts",0,GREEN,True),
])
footer(s,num())

# 10 composition
s=slide(); header(s,"Part 1 · composition","Utricle vs saccule — shared inventory")
pic(s,CA("barplot_celltype_composition_FTL_all.png"),Inches(0.4),Inches(1.4),Inches(7.4),Inches(5.4))
bullets(s,Inches(7.9),Inches(1.6),Inches(5.1),Inches(5.0),[
    ("Both organs dominated by sensory epithelium (hair + supporting)",0,INK,False),
    ("All 11 cell types recovered in both replicates of both tissues",0,INK,False),
    ("Rep-2 samples carry more transitional-epithelial & dark/roof cells (both organs)",0,ACC2,True),
    ("→ a batch/dissection difference (rep1 vs rep2), not utricle-vs-saccule",1,GREY,False),
    ("'HC_type_I' = hair-cell compartment (I/II not reliable)",0,GREY,False),
])
footer(s,num())

# 11 Part 2 PCA
s=slide(); header(s,"Part 2 · organ comparison","Bulk PCA separates the organs — does single-cell?")
pic(s,PB("PCA_organ_pseudobulk.png"),Inches(0.3),Inches(1.4),Inches(6.4),Inches(4.5),
    cap="Raw pseudobulk PCA")
pic(s,PB("PCA_organ_batch_corrected.png"),Inches(6.7),Inches(1.4),Inches(6.4),Inches(4.5),
    cap="After removing the batch effect")
bullets(s,Inches(0.5),Inches(6.0),Inches(12.4),Inches(1.3),[
    ("Raw PCA: PC1 (54%) separates by REPLICATE/BATCH (rep1 vs rep2), not tissue — a strong technical effect dominates",0,ACC2,True),
    ("After removing batch: saccule vs utricle DO separate (PC2) — a real but modest tissue signal, matching your bulk PCA",0,GREEN,True),
])
footer(s,num())

# 12 DE
s=slide(); header(s,"Part 2 · differential expression","Proper edgeR pseudobulk: no gene is significant")
pic(s,PB("volcano_organ.png"),Inches(0.3),Inches(1.4),Inches(6.3),Inches(5.3),
    cap="Organ-level edgeR (batch-blocked)")
pic(s,PB("cellDE_summary_bar.png"),Inches(6.7),Inches(1.5),Inches(6.3),Inches(3.4),
    cap="Cell-level DE per type — all zero")
bullets(s,Inches(6.7),Inches(5.05),Inches(6.3),Inches(2.0),[
    ("edgeR (~replicate + tissue), protein-coding + lncRNA only",0,INK,False),
    ("0 significant genes at organ level AND every cell type (FDR<0.05, |log2FC|>1)",0,ACC2,True),
    ("Top nominal candidates: Acan, Bglap (matrix/otoconial) — need validation",0,GREY,False),
])
footer(s,num())

# 13 reconciliation
s=slide(); header(s,"Part 2 · reconciliation","Why bulk separates but single-cell 'doesn't'")
bullets(s,Inches(0.55),Inches(1.5),Inches(12.2),Inches(5.4),[
    ("It's not a contradiction — the tissue signal IS there (PC2 after batch correction), same direction as your bulk",0,GREEN,True),
    ("Three reasons it's weak in this single-cell dataset:",0,INK,True),
    ("Batch effect — a large rep1-vs-rep2 technical difference dominates the variation (PC1, 54%)",1,INK,False),
    ("Power — only 2 biological replicates per organ → edgeR can't reach FDR for a small effect",1,INK,False),
    ("Composition-controlled — per-cell-type pseudobulk removes the cell-mix differences that bulk 'sees'",1,INK,False),
    ("Your bulk likely had more replicates / less batch structure, so the same modest signal separated cleanly",0,INK,False),
])
footer(s,num())

# 14 summary
s=slide(); header(s,"Summary & next steps","")
bullets(s,Inches(0.55),Inches(1.45),Inches(12.2),Inches(3.6),[
    ("Built a robust cell-type map: one sample → per-organ → all 4, cross-validated by an independent reference",0,INK,False),
    ("Both organs share the same cell types in similar proportions; hair-cell Type I/II is not resolvable at 3′ depth",0,INK,False),
    ("A consistent utricle-vs-saccule signal exists (PC2 after batch correction) but is modest and not significant per-gene",0,INK,False),
    ("Dominant technical factor: a batch effect between the two 10x runs (rep1 vs rep2)",0,INK,False),
])
rect(s,Inches(0.55),Inches(5.15),Inches(12.2),Inches(0.045),ACC2)
bullets(s,Inches(0.55),Inches(5.3),Inches(12.2),Inches(1.8),[
    ("Next steps",0,ACC2,True),
    ("Add biological replicates (≥3–4 per organ) — the single biggest fix for power",1,INK,False),
    ("Validate candidates (Acan/Bglap/otoconial) and hair-cell subtypes by ISH / qPCR",1,INK,False),
    ("Optional: mito/ribosomal blocklist; higher clustering resolution to resolve rare types (neurons, endothelial)",1,INK,False),
])
footer(s,num())

prs.save(OUT); print("SAVED:",OUT,"| slides:",len(prs.slides._sldIdLst))
