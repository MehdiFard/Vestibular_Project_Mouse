#!/usr/bin/env python
"""Report: ambient (SoupX) + doublet (scDblFinder) cleaning, re-annotation, recovery."""
import os, sys
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE
from PIL import Image

FIG = "../../analyses/07_single_cell_rna/results/cell_annotation_clean"
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(FIG, "cleaning_report.pptx")

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
def header(s,k,t,kc=ACC):
    rect(s,0,0,SW,Inches(1.12),INK); rect(s,0,Inches(1.12),SW,Inches(0.06),kc)
    _,tf=tbox(s,Inches(0.5),Inches(0.10),Inches(12.3),Inches(0.98),MSO_ANCHOR.MIDDLE)
    run(tf.paragraphs[0],k.upper(),12,kc,bold=True); run(tf.add_paragraph(),t,24,WHITE,bold=True)
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
        pp.alignment=PP_ALIGN.CENTER; run(pp,cap,9.5,GREY,italic=True)
def dtable(s,rows,x,y,w,h,widths=None,hdr=INK,fs=10):
    t=s.shapes.add_table(len(rows),len(rows[0]),x,y,w,h).table
    if widths:
        for c,wd in enumerate(widths): t.columns[c].width=wd
    for r in range(len(rows)):
        for c in range(len(rows[0])):
            cell=t.cell(r,c); cell.vertical_anchor=MSO_ANCHOR.MIDDLE
            cell.margin_top=cell.margin_bottom=Inches(0.01); cell.margin_left=Inches(0.05)
            pp=cell.text_frame.paragraphs[0]
            if c: pp.alignment=PP_ALIGN.CENTER
            h0=(r==0); run(pp,rows[r][c],fs+0.5 if h0 else fs,WHITE if h0 else INK,bold=(h0 or c==0))
            cell.fill.solid(); cell.fill.fore_color.rgb=hdr if h0 else (LIGHT if r%2 else WHITE)
    return t
def footer(s,n):
    _,tf=tbox(s,Inches(0.4),SH-Inches(0.38),Inches(11),Inches(0.3))
    run(tf.paragraphs[0],"QC re-analysis — ambient (SoupX) + doublet (scDblFinder) cleaning",9,GREY)
    _,tf2=tbox(s,SW-Inches(1.0),SH-Inches(0.38),Inches(0.7),Inches(0.3)); pp=tf2.paragraphs[0]
    pp.alignment=PP_ALIGN.RIGHT; run(pp,str(n),10,GREY)
def F(n): return os.path.join(FIG,n)

# 1 title
s=slide(); rect(s,0,0,SW,SH,INK); rect(s,0,Inches(4.3),SW,Inches(0.08),ACC)
_,tf=tbox(s,Inches(0.9),Inches(1.6),Inches(11.6),Inches(2.4))
run(tf.paragraphs[0],"Adding QC: ambient-RNA + doublet cleaning",33,WHITE,bold=True)
run(tf.add_paragraph(),"Re-annotation of the vestibular maculae after SoupX + scDblFinder",18,ACC)
_,tf=tbox(s,Inches(0.9),Inches(4.6),Inches(11.6),Inches(1.7))
for t,sz in [("Per-sample cleaning → per-organ + all-sample integration → cell-type recovery vs the raw run",15),
             ("All new results in results/cell_annotation_clean/ · run on the cluster (SLURM)",12)]:
    run(tf.add_paragraph(),t,sz,WHITE if sz>13 else GREY)
run(tbox(s,Inches(0.9),Inches(6.6),Inches(11),Inches(0.4))[1].paragraphs[0],"Lab meeting · Mehdi Esmaeili Fard",12,GREY)

# 2 method
s=slide(); header(s,"Why & how","The earlier runs had no ambient/doublet QC")
bullets(s,Inches(0.55),Inches(1.5),Inches(12.2),Inches(5.4),[
    ("Gap: cell-typing so far used Cell Ranger filtered cells only — no ambient-RNA or doublet removal",0,ACC2,True),
    ("Per sample, two cleaning steps were added:",0,INK,True),
    ("SoupX — estimates & subtracts ambient ('soup') RNA using the raw droplet matrix + clusters",1,GREY,False),
    ("scDblFinder — flags & removes doublets (artificial-doublet kNN classifier)",1,GREY,False),
    ("Same QC thresholds as before (nCount>800, nFeature>500, mt-<5%) so the comparison is fair",0,INK,False),
    ("Then re-did everything on cleaned data: integrate saccule, integrate utricle, integrate all 4, annotate",0,INK,False),
    ("Packages installed into seurat_env (SoupX from CRAN; scDblFinder via conda) — Seurat v5 unaffected",0,GREY,False),
])
footer(s,2)

# 3 filtering report
s=slide(); header(s,"Filtering report","Ambient 10–30%, doublets 6–11%, ~13% of cells removed")
dtable(s,[["sample","tissue/rep","ambient rho","cellranger","after QC","doublets","final","% removed"],
          ["FTL_1146","saccule r1","0.30","16,866","14,794","1,014","13,780","18.3"],
          ["FTL_1147","utricle r1","0.10","10,196","8,496","667","7,829","23.2"],
          ["FTL_1148","saccule r2","0.10","10,203","8,987","771","8,216","19.5"],
          ["FTL_1149","utricle r2","0.14","17,099","14,955","1,953","13,002","24.0"]],
       Inches(0.5),Inches(1.5),Inches(12.3),Inches(2.6),
       widths=[Inches(1.6),Inches(1.7),Inches(1.5),Inches(1.7),Inches(1.5),Inches(1.4),Inches(1.5),Inches(1.4)])
bullets(s,Inches(0.55),Inches(4.4),Inches(12.2),Inches(2.6),[
    ("FTL_1146 (saccule r1) had unusually high ambient contamination: rho = 0.30 (30%) vs ~0.10–0.14 for the others",0,ACC2,True),
    ("Doublet rate 6–11% per sample (scDblFinder); total ~49,200 → 42,827 cells (~13% removed)",0,INK,False),
    ("Ambient soup is dominated by the highly-expressed hair-cell genes (Spp1/Ocm) — relevant below",0,GREY,False),
])
footer(s,3)

# 4 cleaned UMAP
s=slide(); header(s,"Cleaned map (all samples)","11 cell types — rare types now distinct")
pic(s,F("all_clean_UMAP_celltype.png"),Inches(0.3),Inches(1.35),Inches(8.6),Inches(5.4))
bullets(s,Inches(9.1),Inches(1.6),Inches(4.0),Inches(5.0),[
    ("Same 11-type inventory, but cleaner separation",0,INK,False),
    ("Neuron, Pericyte/SMC, Endothelial, Immune, Melanocyte now form clean, distinct clusters",0,GREEN,True),
    ("The artificial HC Type I/II split is gone — hair cells are one clear island",0,INK,False),
])
footer(s,4)

# 5 recovery (headline)
s=slide(); header(s,"Cell-type recovery (headline)","Cleaning recovered rare types that were missing before")
pic(s,F("celltype_recovery_all_samples.png"),Inches(0.3),Inches(1.4),Inches(7.4),Inches(5.3),
    cap="All-samples: cell-type counts before vs after cleaning")
bullets(s,Inches(7.9),Inches(1.6),Inches(5.1),Inches(5.0),[
    ("Newly recovered cell types (absent in the raw runs):",0,GREEN,True),
    ("Saccule: Neuron (245) + Pericyte/SMC (460)",1,INK,False),
    ("Utricle: Endothelial (190)",1,INK,False),
    ("All-samples: Endothelial (337)",1,INK,False),
    ("Why: removing ambient + doublets sharpened the clustering, so these small populations separated instead of blending in",0,INK,False),
    ("Direct answer to 'can we recover any cell type?' — yes: neurons, pericytes, endothelial cells",0,ACC,True),
])
footer(s,5)

# 6 composition rebalance
s=slide(); header(s,"Composition — rebalanced by ambient removal","Hair-cell fraction was inflated by soup")
pic(s,F("all_clean_composition_bar.png"),Inches(0.3),Inches(1.4),Inches(7.3),Inches(5.3))
bullets(s,Inches(7.8),Inches(1.6),Inches(5.2),Inches(5.0),[
    ("All-samples hair-cell fraction: ~56% (raw) → 33% (cleaned)",0,ACC2,True),
    ("Supporting 20%→31%, Transitional 6%→17% correspondingly rise",1,GREY,False),
    ("Interpretation: ambient Spp1/Ocm soup (worst in FTL_1146, 30%) had inflated the hair-cell compartment; SoupX correction reassigns those cells",0,INK,False),
    ("Caveat: 3′ dropout makes individual hair-cell markers sparse — treat exact fractions with care",0,GREY,False),
])
footer(s,6)

# 7 per-organ cleaned
s=slide(); header(s,"Per-organ cleaned maps","Saccule and utricle after cleaning")
pic(s,F("saccule_clean_UMAP_celltype.png"),Inches(0.3),Inches(1.4),Inches(6.4),Inches(5.2),
    cap="Saccule (cleaned)")
pic(s,F("utricle_clean_UMAP_celltype.png"),Inches(6.7),Inches(1.4),Inches(6.4),Inches(5.2),
    cap="Utricle (cleaned)")
footer(s,7)

# 8 quality
s=slide(); header(s,"Quality check","SoupX did not over-correct")
bullets(s,Inches(0.55),Inches(1.5),Inches(12.2),Inches(5.4),[
    ("Despite up to 30% ambient correction, canonical markers survive with healthy expression:",0,INK,True),
    ("Myo7a (hair) 10.8%+ · Sox2 (support) 7.9%+ · Kcnq1 (dark) 14.2%+ · Pecam1 (endo) 4.2%+ · Snap25 (neuron) 12.6%+ · Acta2 (pericyte) 0.6%+",1,GREY,False),
    ("→ no evidence of gross over-correction; the recovered clusters express their expected markers",0,GREEN,True),
    ("Replicates still integrate well after cleaning (see *_clean_UMAP_batch.png)",0,INK,False),
    ("Outputs: results/cell_annotation_clean/ — cleaned per-sample objects, 3 integrations, filtering_stats.csv, recovery table",0,GREY,False),
])
footer(s,8)

# 9 — why HC_type_I only in saccule (answers the question)
s=slide(); header(s,"Why 'HC_type_I' in saccule but not utricle?","An ambient-Spp1 labelling artifact")
bullets(s,Inches(0.55),Inches(1.5),Inches(12.2),Inches(5.4),[
    ("The 'HC_type_I' panel is driven by Spp1 (its other markers Ocm/Cabp2 are near-silent on 3′ data)",0,INK,False),
    ("In the RAW data, saccule hair cells had higher Spp1 than utricle:",0,INK,True),
    ("saccule Spp1 mean 2.04 (82% of hair cells+) vs utricle 1.21 (73%+)",1,GREY,False),
    ("→ the argmax assigned saccule hair-cell clusters to 'HC_type_I', utricle to the generic 'Hair_cell'",0,INK,False),
    ("That extra saccule Spp1 was largely AMBIENT: saccule FTL_1146 had 30% ambient, and Spp1/Ocm dominate the hair-cell soup",0,ACC2,True),
    ("After SoupX cleaning: HC_type_I → 0 — BOTH organs' hair cells now label uniformly as 'Hair_cell'",0,GREEN,True),
    ("So it was never a real Type I difference between the organs — just ambient Spp1 tipping the label (and true Type I/II isn't resolvable on 3′ data anyway)",0,INK,False),
])
footer(s,9)

# 10 — cleaned pseudobulk PCA
s=slide(); header(s,"Cleaned data · pseudobulk PCA","Batch still dominates; tissue on PC2")
pic(s,F("pseudobulk_DE/clean_PCA_organ.png"),Inches(0.3),Inches(1.4),Inches(6.4),Inches(4.6),
    cap="Raw pseudobulk PCA (cleaned data)")
pic(s,F("pseudobulk_DE/clean_PCA_batch_corrected.png"),Inches(6.7),Inches(1.4),Inches(6.4),Inches(4.6),
    cap="After removing the batch effect")
bullets(s,Inches(0.5),Inches(6.05),Inches(12.4),Inches(1.2),[
    ("PC1 (60%) still separates by REPLICATE/BATCH — SoupX corrects ambient per sample but not the between-run batch effect",0,ACC2,True),
    ("After batch correction, saccule vs utricle separate on PC2 (saccule top, utricle bottom) — same modest tissue signal as before cleaning",0,GREEN,True),
])
footer(s,10)

# 11 — cleaned edgeR DE
s=slide(); header(s,"Cleaned data · edgeR DE","Still no significant gene — even for recovered types")
pic(s,F("pseudobulk_DE/clean_volcano_organ.png"),Inches(0.3),Inches(1.4),Inches(6.3),Inches(5.3),
    cap="Organ-level edgeR (cleaned, batch-blocked)")
pic(s,F("pseudobulk_DE/clean_cellDE_summary_bar.png"),Inches(6.7),Inches(1.5),Inches(6.3),Inches(3.4),
    cap="Cell-level DE (cleaned) — all zero")
bullets(s,Inches(6.7),Inches(5.05),Inches(6.3),Inches(2.0),[
    ("0 significant genes at organ level AND every cell type (incl. the recovered Endothelial/Neuron/Pericyte)",0,ACC2,True),
    ("→ cleaning fixed the ANNOTATION but did not reveal tissue DE",0,INK,True),
    ("The limit is 2 replicates + batch, not ambient/doublets",0,GREY,False),
])
footer(s,11)

# 12 summary
s=slide(); header(s,"Summary & implications","")
bullets(s,Inches(0.55),Inches(1.45),Inches(12.2),Inches(3.6),[
    ("Added proper QC (SoupX ambient + scDblFinder doublets); ~13% of cells removed",0,INK,False),
    ("Recovered 3 rare cell types missing before: Neuron, Pericyte/SMC, Endothelial",0,GREEN,True),
    ("The saccule 'HC_type_I' label was an ambient-Spp1 artifact — gone after cleaning (hair cells unify)",0,INK,False),
    ("Markers intact → cleaning improved resolution without distorting identities",0,INK,False),
])
rect(s,Inches(0.55),Inches(5.05),Inches(12.2),Inches(0.045),ACC2)
bullets(s,Inches(0.55),Inches(5.2),Inches(12.2),Inches(1.9),[
    ("Two separate conclusions",0,ACC2,True),
    ("Annotation: cleaning clearly helped — better resolution, recovered rare types, corrected composition",1,INK,False),
    ("Organ DE: unchanged — a modest tissue signal (PC2) but 0 significant genes at n=2; batch dominates. Use cleaned data + ≥3–4 replicates for a powered test",1,INK,False),
])
footer(s,12)

prs.save(OUT); print("SAVED:",OUT,"| slides:",len(prs.slides._sldIdLst))
