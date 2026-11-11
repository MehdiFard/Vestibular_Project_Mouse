#!/usr/bin/env python
"""One-slide deck: cell types recovered after QC cleaning (SoupX + scDblFinder)."""
import os
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN

OUT="../../analyses/07_single_cell_rna/results/cell_annotation_clean/recovered_celltypes.pptx"
INK=RGBColor(0x1a,0x1a,0x2e); ACC=RGBColor(0x2e,0x86,0xab); GREY=RGBColor(0x55,0x5b,0x6e)
GRN=RGBColor(0x1b,0x7a,0x3d); WHITE=RGBColor(0xff,0xff,0xff); LT=RGBColor(0xED,0xF3,0xF6)

prs=Presentation(); prs.slide_width=Inches(13.333); prs.slide_height=Inches(7.5)
s=prs.slides.add_slide(prs.slide_layouts[6])

def txt(l,t,w,h,text,size,bold=False,color=INK,align=PP_ALIGN.LEFT,italic=False):
    tf=s.shapes.add_textbox(Inches(l),Inches(t),Inches(w),Inches(h)).text_frame; tf.word_wrap=True
    p=tf.paragraphs[0]; p.alignment=align; r=p.add_run(); r.text=text
    r.font.size=Pt(size); r.font.bold=bold; r.font.color.rgb=color; r.font.italic=italic
    return tf

# title + accent bar
txt(0.6,0.4,12.1,1.0,"Cell types recovered after QC cleaning",28,True,INK)
bar=s.shapes.add_shape(1,Inches(0.62),Inches(1.18),Inches(2.6),Pt(3)); bar.fill.solid(); bar.fill.fore_color.rgb=ACC; bar.line.fill.background()
txt(0.62,1.28,12.0,0.5,"SoupX (ambient-RNA removal) + scDblFinder (doublet removal) uncovered a rare vascular/neural compartment",14,False,GREY,italic=True)

# table
rows=[["Cell type","Analysis","Cells before","Cells after","Status"],
      ["Neuron","saccule","0","245","recovered (new)"],
      ["Pericyte_SMC  (pericyte / smooth-muscle)","saccule","0","460","recovered (new)"],
      ["Endothelial","utricle","0","190","recovered (new)"],
      ["Endothelial","all samples","0","337","recovered (new)"]]
nr,nc=len(rows),len(rows[0])
gt=s.shapes.add_table(nr,nc,Inches(0.8),Inches(2.3),Inches(11.7),Inches(2.9)).table
gt.columns[0].width=Inches(4.6); gt.columns[1].width=Inches(2.0)
gt.columns[2].width=Inches(1.7); gt.columns[3].width=Inches(1.5); gt.columns[4].width=Inches(1.9)
for i,row in enumerate(rows):
    gt.rows[i].height=Inches(0.55)
    for j,val in enumerate(row):
        c=gt.cell(i,j); c.text=""
        para=c.text_frame.paragraphs[0]; para.alignment=PP_ALIGN.CENTER if j>0 else PP_ALIGN.LEFT
        r=para.add_run(); r.text=val; r.font.size=Pt(13)
        if i==0:
            c.fill.solid(); c.fill.fore_color.rgb=INK; r.font.bold=True; r.font.color.rgb=WHITE
        else:
            c.fill.solid(); c.fill.fore_color.rgb=LT if i%2 else WHITE
            r.font.color.rgb=GRN if j==4 else INK
            if j==4: r.font.bold=True

# note
txt(0.8,5.5,11.7,1.4,
    "All four had 0 cells before cleaning — they were masked by ambient RNA and doublets. "
    "Cleaning also corrected a hair-cell artifact: the ambient-Spp1–inflated “HC_type_I” cluster "
    "collapsed (24,262 → 0 in the combined run) while true Hair_cell counts rose (3,286 → 14,096).",
    13,False,GREY)
txt(0.8,6.9,11.7,0.4,"Source: results/cell_annotation_clean/celltype_recovery_comparison.csv",10,False,GREY,italic=True)

prs.save(OUT); print("SAVED",OUT)
