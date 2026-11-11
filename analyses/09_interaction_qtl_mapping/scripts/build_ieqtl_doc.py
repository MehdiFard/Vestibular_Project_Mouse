#!/usr/bin/env python
"""Detailed interaction-eQTL mapping write-up (Word). Env1 python."""
import os
from docx import Document
from docx.shared import Pt, RGBColor, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH as AL
from docx.enum.table import WD_TABLE_ALIGNMENT
RES="../../analyses/09_interaction_qtl_mapping/results"
IE=os.path.join(RES,"ieqtl"); PLT=os.path.join(IE,"plots")
OUT=os.path.join(IE,"ieQTL_Mapping_Report.docx")
INK=RGBColor(0x1a,0x1a,0x2e); ACC=RGBColor(0x2e,0x86,0xab); GREY=RGBColor(0x55,0x5b,0x6e)
doc=Document(); st=doc.styles["Normal"]; st.font.name="Calibri"; st.font.size=Pt(11)
st.paragraph_format.space_after=Pt(6); st.paragraph_format.line_spacing=1.08
def para(t="",size=11,bold=False,italic=False,color=None,align=None,after=6):
    p=doc.add_paragraph(); p.paragraph_format.space_after=Pt(after)
    if align: p.alignment=align
    r=p.add_run(t); r.bold=bold; r.italic=italic; r.font.size=Pt(size)
    if color is not None: r.font.color.rgb=color
    return p
def h1(t):
    p=doc.add_heading(t,1)
    for r in p.runs: r.font.color.rgb=INK
def h2(t):
    p=doc.add_heading(t,2)
    for r in p.runs: r.font.color.rgb=ACC
def bullet(t,lead=None):
    p=doc.add_paragraph(style="List Bullet"); p.paragraph_format.space_after=Pt(3)
    if lead:
        r=p.add_run(lead); r.bold=True; p.add_run(t)
    else: p.add_run(t)
def figure(path,cap,w=6.0):
    if os.path.exists(path):
        doc.add_picture(path,width=Inches(w)); doc.paragraphs[-1].alignment=AL.CENTER
    para(cap,size=9.5,italic=True,color=GREY,align=AL.CENTER,after=12)

t=doc.add_heading("Cell-type interaction-eQTL mapping of the vestibular maculae",0)
for r in t.runs: r.font.color.rgb=INK
para("Recovering cell-type-specific genetic regulation from bulk RNA-seq (GTEx / Kim-Hellmuth 2020)",
     size=13,italic=True,color=GREY,align=AL.CENTER)
para("Mouse vestibular project · 2026",size=11,color=GREY,align=AL.CENTER,after=12)

h1("1. Overview")
para("Standard cis-eQTL mapping treats each bulk sample as a single unit and reports the average genetic "
     "effect on expression across all the cell types it contains. This analysis instead asks in which cell "
     "type a regulatory effect acts, by testing whether a variant's effect on expression scales with the "
     "abundance of a given cell type — a cell-type interaction eQTL (ieQTL). The cell-type abundance proxy is "
     "the ssGSEA marker-enrichment score produced in the deconvolution step (fractions from Bisque/MuSiC were "
     "shown to be unreliable here; the enrichment score is the robust signal). This document reports the "
     "mapping and compares it with the standard cis-eQTLs mapped earlier in Pantry.")

h1("2. Methods")
para("For each Pantry group (saccule and utricle, at weeks 8 and 22) and each of six abundant cell types "
     "(hair, supporting, transitional-epithelial, glia/Schwann, fibrocyte, roof/dark), we fitted:")
para("expression ~ genotype + covariates + enrichment + genotype × enrichment,", italic=True, align=AL.CENTER)
para("testing the genotype × enrichment interaction term for every gene–cis-variant pair within a 1 Mb window.")
bullet("genotypes, quantile-normalized expression, and covariates were the exact per-group inputs already "
       "prepared in Pantry (genotype PCs + expression PCs + sex).", "Inputs: ")
bullet("the ssGSEA cell-type enrichment score, inverse-normal transformed across the group's samples, added "
       "as a single interaction term per scan.", "Interaction covariate: ")
bullet("tensorQTL 1.0.10 (cis.map_nominal) on an NVIDIA A100 GPU; ~2.5 minutes per scan.", "Engine: ")
bullet("MAF ≥ 0.05 was required in both the upper and lower halves of the enrichment distribution "
       "(maf_threshold_interaction=0.05), so an interaction cannot be driven by a genotype seen only at one "
       "extreme of cell abundance; significance used eigenMT (effective number of cis-variants per gene) with "
       "Benjamini-Hochberg across genes, calling ieGenes at FDR < 0.05.", "Safeguards: ")

h1("3. Results")
h2("3.1 ieGene discovery")
para("Across the four groups and six cell types we identified 42 significant ieGenes (FDR < 0.05). Discovery "
     "was concentrated in the two abundant sensory cell types — hair cells (10 in saccule week 8, 8 in utricle "
     "week 8) and supporting cells (9 and 3) — the cell types whose enrichment scores are best estimated. Both "
     "week-22 groups yielded essentially no ieGenes, mirroring the power difference seen in the bulk xQTLs: the "
     "week-8 groups have more samples (n≈73–74 vs 55–58) and therefore more power for the harder interaction "
     "test.")
figure(os.path.join(PLT,"B_iegene_heatmap.png"),
       "Figure 1. Number of ieGenes (FDR<0.05) by group and cell type. Hair and supporting cells in the "
       "week-8 groups carry almost all of the signal.", 5.8)

h2("3.2 Comparison with standard cis-eQTL")
para("The standard cis-eQTL scan (Pantry) identifies thousands of eGenes per group (7,559 in saccule week 8; "
     "8,073 in utricle week 8), because it captures every marginal genetic effect on expression regardless of "
     "which cell type produces it. The interaction test is a different, harder question and returns a much "
     "smaller, focused set (Figure 2). The two are complementary rather than competing: the eQTL scan "
     "catalogues regulatory variants, and the interaction scan annotates a subset of them with the cell type "
     "in which they act.")
figure(os.path.join(PLT,"A_egene_vs_iegene_counts.png"),
       "Figure 2. Standard eGenes vs. interaction ieGenes per group (log scale).", 5.8)
para("Of the 42 ieGene hits, 37 are also standard eGenes — i.e. they have a detectable marginal effect that "
     "the interaction now localizes to a cell type — while 5 are interaction-only, with no significant marginal "
     "eQTL. These interaction-only genes are the clearest illustration of the method's value: regulation that "
     "is invisible to a bulk eQTL scan but emerges once cell-type abundance is modelled.")
figure(os.path.join(PLT,"C_iegene_overlap.png"),
       "Figure 3. ieGenes that are also standard eGenes vs. interaction-only.", 4.3)

h2("3.3 Example ieGenes")
para("Representative hits include, in hair cells, Slc7a11, Fam124a, Nedd4, Tlr3 (saccule week 8) and Slk, "
     "Suclg2, Jade2 (utricle week 8); and in supporting cells, Znrf3 (the strongest single hit, BH p = 3×10⁻⁴), "
     "Col20a1 and Mrap. Roof/dark cells contributed Pla2g7 and Gdpd5, and transitional-epithelial cells Lrrc6. "
     "The full annotated list is in ieGenes_significant.csv.")

h1("4. Caveats")
bullet("no ieGene reached significance in more than one group (age or organ); with modest per-group sample "
       "sizes and different animals across groups, individual hits should be treated as candidates.",
       "No cross-group replication: ")
bullet("some genes appear for several cell types within one group (e.g. Znrf3, Slc7a11). Because the abundant "
       "cell types' enrichment scores are positively correlated, a single interaction can surface under more "
       "than one cell type and does not by itself prove independent cell-type specificity.",
       "Correlated enrichment: ")
bullet("a few hits are ribosomal or predicted genes (Rps26, several Gm*/Rik) or biologically unexpected "
       "(Olfr1193) and are lower confidence.", "Technical hits: ")
bullet("the per-group covariates include eight expression (RNA) PCs, which can themselves capture cell "
       "composition; an interaction that survives only without a composition-correlated PC (or vanishes only "
       "when one is added) should be checked for over-adjustment. A covariate-sensitivity analysis is planned.",
       "Covariate over-adjustment: ")

h1("5. Next steps")
bullet("pool the two ages within each organ (n≈128–132) to roughly double the sample size and power for the "
       "interaction test; this requires merged genotypes, re-normalized expression and recomputed covariates.")
bullet("run the covariate-sensitivity models (core covariates ± RNA PCs 1–6) and report the correlation "
       "between each RNA PC and the tested cell-type enrichment.")
bullet("for robust ieGenes, colocalize with the vestibular-function GWAS (VsEP, beam-walk) to assign balance "
       "loci to specific cell types.")
para("Scripts: scripts/05_interaction_eqtl.py, 06_ieqtl_vs_eqtl_compare.py, run_ieqtl_gpu.sbatch. "
     "Outputs: results/ieqtl/.", size=10, italic=True, color=GREY)
doc.save(OUT); print("SAVED",OUT)
