#!/usr/bin/env python
"""Detailed deconvolution write-up (Word). Run with Env1 python (has python-docx)."""
import os
from docx import Document
from docx.shared import Pt, RGBColor, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH

RES = "../../analyses/08_deconvolution/results"
PLT = os.path.join(RES, "plots")
OUT = os.path.join(os.path.dirname(RES), "report", "Deconvolution_Report.docx")
INK=RGBColor(0x1a,0x1a,0x2e); ACC=RGBColor(0x2e,0x86,0xab); GREY=RGBColor(0x55,0x5b,0x6e)

doc=Document()
st=doc.styles["Normal"]; st.font.name="Calibri"; st.font.size=Pt(11)
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
    return p
def h2(t):
    p=doc.add_heading(t,2)
    for r in p.runs: r.font.color.rgb=ACC
    return p
def bullet(t,lead=None):
    p=doc.add_paragraph(style="List Bullet"); p.paragraph_format.space_after=Pt(3)
    if lead:
        r=p.add_run(lead); r.bold=True; p.add_run(t)
    else: p.add_run(t)
def figure(path,cap,w=6.0):
    if os.path.exists(path):
        doc.add_picture(path,width=Inches(w)); doc.paragraphs[-1].alignment=WD_ALIGN_PARAGRAPH.CENTER
    para(cap,size=9.5,italic=True,color=GREY,align=WD_ALIGN_PARAGRAPH.CENTER,after=12)

# cover
t=doc.add_heading("Cell-type deconvolution of bulk vestibular RNA-seq",0)
for r in t.runs: r.font.color.rgb=INK
para("Detailed methods and results — toward cell-type interaction-eQTL mapping",
     size=13,italic=True,color=GREY,align=WD_ALIGN_PARAGRAPH.CENTER)
para("Mouse vestibular project · 2026",size=11,color=GREY,align=WD_ALIGN_PARAGRAPH.CENTER,after=12)

h1("1. Aim and rationale")
para("A bulk RNA-seq sample is a mixture of many cell types, and its measured expression is the "
     "weighted average of the expression of those cell types. Standard eQTL mapping therefore averages "
     "genetic-regulatory effects over the whole tissue: a variant that regulates a gene only in one cell "
     "type is diluted by every other cell type and is easily missed or reported as a washed-out effect. "
     "The vestibular maculae are a clear example — hair cells and supporting cells dominate the tissue, so "
     "a regulatory effect specific to a rare cell type is essentially invisible in bulk.")
para("The strategy of Kim-Hellmuth et al. (GTEx, Science 2020) recovers this signal without needing "
     "single-cell data from every genotyped individual. First, computational deconvolution estimates how "
     "much of each cell type is present in each bulk sample. Then a genotype × cell-type-abundance "
     "interaction term is added to the eQTL model; where a variant's effect scales with a cell type's "
     "abundance, that is a cell-type interaction eQTL (ieQTL) — a cell-type-specific effect recovered from "
     "bulk. This document covers the first half only: producing a reliable per-sample cell-type signal. "
     "The QTL mapping itself is the next step.")

h1("2. Data and design")
bullet("saccule = saccule_w8 + saccule_w22 (128 samples); utricle = utricle_w8 + utricle_w22 (132 samples). "
       "Bulk expression is kallisto gene-level TPM (linear scale), gene IDs mapped ENSMUSG → symbol.",
       "Bulk (ages pooled within organ, organs separate): ")
bullet("the cleaned, jointly integrated single-nucleus atlas (all_integrated_clean.rds), 4 donors "
       "(FTL_1146–1149), downsampled to ≤2000 cells per type (11,094 cells).", "Single-cell reference: ")
bullet("Hair_cell (HC_type_I/II collapsed — that split is an annotation-threshold artefact), Supporting_cell, "
       "Transitional_epith, Glia_Schwann, Fibrocyte_mesench, Roof_dark_cell.", "Six abundant cell types: ")
bullet("restricted to protein-coding genes with mitochondrial and ribosomal genes removed (16,840 genes "
       "common to bulk and reference). This is essential: the bulk TPM is otherwise dominated by technical "
       "genes (e.g. the miRNA Mir6236 at ~59,000 TPM, several Gm* small RNAs, ~10% mitochondrial) that are "
       "absent or divergent in nuclei and distort the decomposition.", "Gene filtering: ")

h1("3. Methods")
h2("3.1 Deconvolution methods")
para("Three complementary methods were run for each organ:")
bullet("learns a gene-specific transformation between the single-cell reference and bulk expression, then "
       "decomposes. Because the bulk and single-cell come from different animals, it was run with "
       "use.overlap = FALSE (no shared samples). Installed from source (cozygene/bisque; the CRAN package "
       "is archived).", "Bisque (planned primary): ")
bullet("weights genes by their consistency across reference donors, up-weighting stable cell-type markers. "
       "Run via music_prop with the four donors as subjects.", "MuSiC (sensitivity): ")
bullet("a marker-enrichment score: for each bulk sample it scores how enriched each cell type's marker set "
       "is (rank-based, within-sample), using the top 50 markers per type from the reference. This is the "
       "GTEx-style approach (GTEx used xCell, an enrichment method, rather than absolute fractions).",
       "ssGSEA (enrichment): ")
para("All analyses were run on the compute cluster (Slurm) in an isolated conda environment "
     "(deconv_env), never on the login node.")

h2("3.2 Validation strategy")
para("Absolute fractions cannot be checked against a ground truth in the real bulk (the true composition of "
     "each HMDP animal is unknown). We therefore validated the method mechanics on pseudobulk: summing the "
     "single-cell counts within each donor yields four pseudobulk samples whose true composition is known "
     "exactly. Deconvolving these isolates method mechanics from the biological domain shift between nuclei "
     "and whole cells: if a method cannot even recover the composition of pseudobulk made from the same "
     "cells, the setup is at fault; if it can, then failure on real bulk points to the domain shift.")

h1("4. Results")
h2("4.1 Method mechanics are correct (pseudobulk)")
para("On pseudobulk with known composition, MuSiC recovered the true proportions almost perfectly "
     "(Pearson r = 0.95; e.g. donor FTL_1146 hair cells estimated 0.31 vs. true 0.31, supporting cells "
     "0.42 vs. 0.43). Bisque was markedly weaker (r = 0.45), compressing all types toward a uniform ~0.15–0.24. "
     "Two conclusions follow: the pipeline is implemented correctly, and — contrary to the a-priori plan — "
     "MuSiC is the stronger method in-domain, not Bisque.")
figure(os.path.join(PLT,"00_pseudobulk_validation.png"),
       "Figure 1. Estimated vs. true proportion on pseudobulk. MuSiC (r = 0.95) tracks the diagonal; "
       "Bisque (r = 0.45) is compressed toward uniform.", 6.3)

h2("4.2 Absolute fractions are not reliable on the real bulk")
para("On the real bulk, both fraction methods degrade badly. MuSiC collapses, assigning ~0.85–0.93 of every "
     "sample to supporting cells with almost no cross-sample variation; Bisque returns near-uniform "
     "proportions (~0.18 for most types). The two methods barely agree across samples (median Spearman ≈ 0.18). "
     "Restricting the genes to the cell-type markers made MuSiC worse (~1.0 supporting), not better.")
para("The cause is a domain shift: the reference is single-NUCLEUS RNA whereas the bulk is whole-cell / "
     "total RNA (kallisto TPM). Nuclear and whole-cell transcriptomes differ systematically (nuclear-retained "
     "transcripts, intronic content, gene-length and cytoplasmic-abundance biases), which breaks the linear-"
     "mixing assumption both methods rely on. This is compounded by an incomplete reference: a whole-organ "
     "dissection includes compartments not represented by the six sensory/stromal types (bone, blood, nerve), "
     "so the non-negative least-squares fit dumps unexplained signal into whichever type best absorbs it "
     "(here, supporting cells). The Bisque and MuSiC fractions should therefore not be interpreted as real "
     "cell-type composition.")
figure(os.path.join(PLT,"01_mean_composition_by_method.png"),
       "Figure 2. Mean composition by method vs. observed single-cell composition. On real bulk MuSiC "
       "collapses to supporting cells and Bisque is near-uniform, unlike the observed scRNA composition.", 6.0)

h2("4.3 The usable signal: ssGSEA marker enrichment")
para("The ssGSEA enrichment score is robust to the domain shift because it is rank-based and computed within "
     "each sample: it asks only whether a cell type's marker genes are relatively high in a given sample, not "
     "how the sample decomposes as a linear mixture. This is precisely why GTEx used an enrichment method "
     "rather than absolute fractions. The enrichment scores show real cross-sample variation (e.g. hair-cell "
     "score ranges 0.01–0.57) and biologically sensible structure (supporting and hair cells score highest, "
     "rare types lowest). They are also the consensus signal that both fraction methods partially track "
     "(Bisque vs ssGSEA Spearman 0.49; MuSiC vs ssGSEA 0.41). The best-behaved types — with the clearest "
     "cross-sample spread — are Hair_cell, Transitional_epith, Glia_Schwann and Roof_dark_cell.")
figure(os.path.join(PLT,"06_ssgsea_enrichment.png"),
       "Figure 3. ssGSEA marker-enrichment scores across bulk samples, per organ. These per-sample scores are "
       "the cell-type covariate used for interaction-eQTL mapping.", 6.0)

h1("5. Conclusion and next steps")
para("Deconvolution was implemented and validated. The key result is a methodological one: absolute-fraction "
     "deconvolution (Bisque, MuSiC) is not trustworthy for this dataset because the single-nucleus reference "
     "and the whole-cell bulk occupy different measurement domains, and the six-type reference does not span "
     "the whole-organ bulk. The ssGSEA marker-enrichment score is the reliable per-sample cell-type signal and "
     "is what we will carry forward.")
para("For the interaction-eQTL mapping (next step):")
bullet("use the ssGSEA enrichment score as the interaction covariate, inverse-normal transformed across the "
       "mapping samples;")
bullet("require MAF ≥ 0.05 separately in both the upper and lower halves of the enrichment distribution, so an "
       "interaction cannot be driven by a genotype seen only at one extreme of cell abundance;")
bullet("correct with eigenMT (effective number of cis-variants per gene) and Benjamini-Hochberg across genes, "
       "calling cell-type ieGenes at FDR < 0.05;")
bullet("focus on the best-powered cell types (hair, supporting, transitional-epithelial, glia, roof/dark).")

para("All scripts are in scripts/ and all outputs (tables and plots) in results/; "
     "the full pipeline reruns with scripts/run_deconv.sbatch.", size=10, italic=True, color=GREY)

doc.save(OUT); print("SAVED", OUT)
