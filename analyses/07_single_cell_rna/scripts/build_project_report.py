#!/usr/bin/env python
"""Funder project report (Word): bulk + Pantry xQTL/TWAS + single-nucleus RNA-seq.
Mouse arm only. Official/progress-report tone."""
import os, sys
from docx import Document
from docx.shared import Pt, RGBColor, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT

ROOT = "../../analyses/07_single_cell_rna"
RES  = os.path.join(ROOT, "results")
CLN  = os.path.join(RES, "cell_annotation_clean")
ASS  = os.path.join(RES, "report_assets")
OUT  = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "report", "Vestibular_Project_Report_2026.docx")
INK=RGBColor(0x1a,0x1a,0x2e); ACC=RGBColor(0x2e,0x86,0xab); GREY=RGBColor(0x55,0x5b,0x6e)

# ---------------- Pantry numbers, read from the final run (03 pheast output, 04 TWAS hits) ----------------
import gzip, csv
ANA = os.path.dirname(ROOT)
GROUPS = ["saccule_w8", "saccule_w22", "utricle_w8", "utricle_w22"]
MODS = [("expression", "Gene expression"), ("isoforms", "Isoform ratio"), ("splicing", "Splicing (intron-excision)"),
        ("stability", "RNA stability"), ("alt_TSS", "Alternative TSS"), ("alt_polyA", "Alternative polyA"),
        ("RNA_editing", "RNA editing")]
def n_rows(path):
    with gzip.open(path, "rt") as f: return sum(1 for _ in f) - 1
QTL = {(m, g): n_rows(f"{ANA}/03_pantry/pheast/output/{g}/qtl/{m}.cis_independent_qtl.txt.gz") for m, _ in MODS for g in GROUPS}
QTL_GROUP = {g: sum(QTL[(m, g)] for m, _ in MODS) for g in GROUPS}
QTL_TOTAL = sum(QTL_GROUP.values())
# TWAS hits per group (P < 5e-8 / 7 modalities; 04_twas writes the same all-modality table into every modality folder)
TWAS = {}
for g in GROUPS:
    with open(f"{ANA}/04_twas/{g}/output/expression/twas_hits.tsv") as f:
        TWAS[g] = [r for r in csv.DictReader(f, delimiter="\t")]
TWAS_N = {g: len(TWAS[g]) for g in GROUPS}
TWAS_TOTAL = sum(TWAS_N.values())
gene_of = lambda r: r["ID"].split("__")[0]
TWAS_GENES = len({gene_of(r) for g in GROUPS for r in TWAS[g]})
TWAS_PAIRS = len({(gene_of(r), r["TRAIT"]) for g in GROUPS for r in TWAS[g]})
fmt = lambda n: f"{n:,}"

doc = Document()
st = doc.styles["Normal"]; st.font.name="Calibri"; st.font.size=Pt(11)
st.paragraph_format.space_after=Pt(6); st.paragraph_format.line_spacing=1.08

def para(text="", size=11, bold=False, italic=False, color=None, align=None, after=6):
    p=doc.add_paragraph(); p.paragraph_format.space_after=Pt(after)
    if align: p.alignment=align
    r=p.add_run(text); r.bold=bold; r.italic=italic; r.font.size=Pt(size)
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
def bullet(text, bold_lead=None):
    p=doc.add_paragraph(style="List Bullet"); p.paragraph_format.space_after=Pt(3)
    if bold_lead:
        r=p.add_run(bold_lead); r.bold=True; p.add_run(text)
    else: p.add_run(text)
    return p
def figure(path, caption, width=6.0):
    if os.path.exists(path):
        doc.add_picture(path, width=Inches(width))
        doc.paragraphs[-1].alignment=WD_ALIGN_PARAGRAPH.CENTER
    else:
        para("[missing figure: %s]" % os.path.basename(path), italic=True, color=GREY,
             align=WD_ALIGN_PARAGRAPH.CENTER)
    return para(caption, size=9.5, italic=True, color=GREY,
                align=WD_ALIGN_PARAGRAPH.CENTER, after=12)

# ---------------- cover ----------------
t=doc.add_heading("Vestibular Project — Mouse Data Report",0)
for r in t.runs: r.font.color.rgb=INK
para("Regulatory genetics and single-cell characterization of the inner-ear vestibular maculae",
     size=13, italic=True, color=GREY, align=WD_ALIGN_PARAGRAPH.CENTER)
para("Progress report to funders · 2026", size=11, color=GREY,
     align=WD_ALIGN_PARAGRAPH.CENTER, after=14)

# ---------------- executive summary ----------------
h1("Executive summary")
para("The mouse arm of the vestibular project is now data-complete and has moved from data generation into "
     "integrative analysis. Bulk RNA-seq of the saccule and utricle across two ages in the Hybrid Mouse "
     "Diversity Panel (HMDP) has been processed end-to-end, and we have extended it well beyond a conventional "
     "expression-QTL scan: using the laboratory's Pantry framework we mapped quantitative-trait loci for seven "
     "distinct molecular phenotypes and, by integrating these with existing vestibular-function GWAS, performed "
     "a transcriptome-wide association study (TWAS) that nominates candidate genes for balance phenotypes. In "
     "parallel, we generated and analyzed single-nucleus RNA-seq (snRNA-seq) of both organs, producing a "
     "quality-controlled, cell-type-resolved atlas of the vestibular maculae. Together these two layers — a "
     "population-scale map of regulatory variation and a cell-type-resolved map of the tissue — position the "
     "project for cell-type-aware regulatory mapping, the central goal of the next phase.")

# ---------------- overview ----------------
h1("Background and objectives")
para("The vestibular maculae — the utricle and the saccule — are the inner-ear sensory organs that detect "
     "linear acceleration and head position and are central to balance. The goal of this project is to "
     "understand how genetic variation shapes gene regulation in these organs and how that regulation relates "
     "to vestibular function and its decline with age. Building on the bulk RNA-seq and QTL mapping summarized "
     "in the previous report (March 2025), the current period delivered two advances: (i) completion of "
     "multi-phenotype QTL mapping and a first TWAS across seven molecular phenotypes, and (ii) generation and "
     "analysis of single-nucleus RNA-seq that resolves the cell types of both organs. This report covers the "
     "mouse arm only; the parallel human study has been completed and submitted for publication, and will be "
     "reported separately.")

# ---------------- data ----------------
h1("Data")
h2("Bulk RNA-seq (Hybrid Mouse Diversity Panel)")
para("Saccule and utricle were profiled from the HMDP at two ages — week 8 (young) and week 22 (aged). In "
     "total 321 libraries were sequenced across four batches (saccule: 160; utricle: 161), split evenly "
     "between the two ages. After sequence quality control and matching to genotypes, the usable sample sizes "
     "per group were: saccule week 8 = 73, saccule week 22 = 55, utricle week 8 = 74, utricle week 22 = 58 "
     "(~30,700 genes quantified per group). Genotypes were available for the full panel, enabling cis-QTL "
     "mapping in every tissue-by-age group.")
h2("Single-nucleus RNA-seq")
para("Four snRNA-seq libraries were generated — two biological replicates each of utricle and saccule — on "
     "the 10x Genomics Chromium 3′ v4 platform (mouse, GRCm39). Approximately 49,000 nuclei passed initial "
     "quality control and entered the analysis described below. With these datasets in hand, the mouse "
     "resource is complete.")

# ---------------- bulk + pantry ----------------
h1("Bulk RNA-seq and multi-phenotype QTL mapping (Pantry)")
h2("Expression profiling and differential expression")
para("Reads were quality-controlled and trimmed (FastQC, Trimmomatic), aligned to the mouse reference genome "
     "(GRCm38) with STAR, and gene expression was quantified with RSEM. Principal-component analysis of the "
     "expression matrix cleanly separated the two organs along the first principal component (~25% of variance), "
     "confirming that the saccule and utricle carry distinct bulk transcriptomes (Fig. 1). Differential-"
     "expression analysis contrasting young and aged animals identified 23 up-regulated and 40 down-regulated "
     "genes (|log2 fold-change| > 1, FDR < 5%), highlighting age-associated changes in both organs (Fig. 2).")
figure(os.path.join(ASS,"bulk_expression_PCA.jpeg"),
       "Fig. 1. Bulk expression PCA. Samples separate by organ on PC1 (saccule vs. utricle); "
       "shape denotes age (week 8 / week 22).", 5.2)
figure(os.path.join(ASS,"bulk_volcano.png"),
       "Fig. 2. Differential expression between young and aged animals: 23 up- (red) and 40 down-regulated "
       "(blue) genes at |log2FC| > 1 and FDR < 5%.", 5.2)

h2("Regulatory mapping across seven molecular phenotypes")
para("To move beyond a standard expression-QTL scan, we applied Pantry (PAN-TRanscriptomic phenotYping; Munro "
     "et al., Nature Communications 2024) to extract seven molecular phenotypes from the same RNA-seq data — gene "
     "expression, isoform ratio, intron-excision (splicing) ratio, RNA stability, alternative transcription-"
     "start-site (TSS) usage, alternative polyadenylation (polyA) usage, and A-to-I RNA editing — and mapped "
     "cis-QTLs for each, in every tissue-by-age group, with tensorQTL. RNA editing levels were measured at known "
     "mouse editing sites; sites covered by too few reads in a sample were treated as missing, and the expression "
     "of each site's host gene was regressed out so that editing QTLs are not simply expression QTLs.")
para(f"Across the four groups we identified {fmt(QTL_TOTAL)} significant independent cis-QTLs spanning the "
     f"seven phenotypes (Table 1). Gene expression alone yielded {fmt(min(QTL[('expression', g)] for g in GROUPS))}–"
     f"{fmt(max(QTL[('expression', g)] for g in GROUPS))} independent QTLs per group, and the "
     "additional phenotypes contributed thousands of further regulatory signals — substantially "
     "expanding the regulatory landscape relative to a conventional eQTL analysis and capturing variation in "
     "splicing, RNA stability, and transcript-boundary usage that expression alone would miss. The young "
     "(week 8) groups yielded more QTLs than the aged (week 22) groups; this partly reflects their larger "
     "sample sizes, and hence greater statistical power, in addition to any biological age effect. Cis-"
     "heritability was estimated for all phenotypes to support the association analyses below.")

# Table 1
para("Table 1. Significant independent cis-QTLs per molecular phenotype and group.",
     size=10, bold=True, after=4)
rows=[["Molecular phenotype","Saccule w8","Saccule w22","Utricle w8","Utricle w22"]] + \
     [[lab] + [fmt(QTL[(m, g)]) for g in GROUPS] for m, lab in MODS] + \
     [["Total"] + [fmt(QTL_GROUP[g]) for g in GROUPS]]
tb=doc.add_table(rows=len(rows), cols=5); tb.alignment=WD_TABLE_ALIGNMENT.CENTER
try: tb.style="Light Grid Accent 1"
except Exception: pass
for i,row in enumerate(rows):
    for j,val in enumerate(row):
        c=tb.cell(i,j); c.text=""
        rr=c.paragraphs[0].add_run(val); rr.font.size=Pt(9.5)
        rr.bold = (i==0 or j==0 or i==len(rows)-1)
        if j>0: c.paragraphs[0].alignment=WD_ALIGN_PARAGRAPH.CENTER
para("", after=8)

h2("Transcriptome-wide association with vestibular function (TWAS)")
para("Finally, we integrated the molecular-QTL models with existing HMDP vestibular-function GWAS — vestibular "
     "sensory-evoked potential (VsEP) amplitude, latency and threshold, and a beam-walking balance assay in "
     "young and aged animals — in a transcriptome-wide association study (FUSION framework, with statistical "
     f"colocalization), one molecular phenotype at a time. At P < 5×10⁻⁸/7 this analysis nominated {fmt(TWAS_TOTAL)} "
     f"significant associations across the seven molecular phenotypes and four tissue-by-age groups "
     f"(saccule week 8: {TWAS_N['saccule_w8']}, week 22: {TWAS_N['saccule_w22']}; utricle week 8: "
     f"{TWAS_N['utricle_w8']}, week 22: {TWAS_N['utricle_w22']}), involving {fmt(TWAS_PAIRS)} distinct gene–trait "
     f"pairs and {fmt(TWAS_GENES)} genes. TWAS prioritizes genes whose genetically regulated expression, splicing, or RNA processing is "
     "statistically associated with vestibular function, and provides a direct, mechanism-anchored shortlist of "
     "candidate genes for follow-up.")

# ---------------- snRNA ----------------
h1("Single-nucleus RNA-seq analysis")

h2("Read mapping and count generation")
para("Sequencing reads for the four libraries were processed with the 10x Genomics Cell Ranger pipeline "
     "(v9.0.1). Reads were demultiplexed, barcode- and UMI-corrected, and aligned to the mouse reference "
     "genome (GRCm39) with the Cell Ranger STAR-based aligner; confidently mapped, exonic reads were then "
     "collapsed by UMI to produce a nucleus-by-gene count matrix per library. These filtered count matrices "
     "were the input to the quality control, integration, and annotation steps described below.")

h2("Quality control (ambient RNA and doublets)")
para("Each library was cleaned before analysis. Ambient (“soup”) RNA — cell-free transcripts co-captured with "
     "nuclei — was estimated and subtracted with SoupX, and doublets (droplets containing more than one "
     "nucleus) were removed with scDblFinder. Ambient contamination was substantial — up to ~30% in one "
     "saccule library and ~10–14% in the others — and doublets accounted for 6–11% of nuclei; overall ~13% of "
     "nuclei were removed. Critically, this cleaning recovered three cell populations that had been masked in "
     "the raw data — neurons, pericytes/smooth-muscle, and endothelial cells — which now form distinct "
     "clusters. Canonical markers (e.g. Myo7a for hair cells, Sox2 for supporting cells, Pecam1 for "
     "endothelium) were preserved, confirming that the correction sharpened rather than distorted cell "
     "identities.")

h2("Cell-type annotation")
para("Nuclei were integrated across replicates (per organ and jointly) and annotated with curated panels of "
     "canonical inner-ear marker genes. We resolved 11 cell types spanning the sensory epithelium (hair cells, "
     "supporting cells, transitional epithelial cells, and roof/dark cells) and the surrounding non-sensory "
     "compartment (glia/Schwann cells, neurons, immune/macrophage, endothelial cells, melanocytes, fibrocytes, "
     "and pericytes/smooth-muscle). An independent published mouse-utricle reference (Burns et al., 2015), "
     "transferred with SingleR, agreed with our marker-based labels for 86.9% of epithelial cells, providing "
     "external validation. Both organs recovered essentially the same cell types (Fig. 3). Within the sensory "
     "epithelium, per-organ integration resolved a hair-cell subpopulation that is labeled HC_type_I in the "
     "saccule and HC_type_II in the utricle; when all four samples are integrated jointly these merge into a "
     "single hair-cell cluster, indicating that the split reflects an annotation-threshold effect on a "
     "continuous population rather than a strictly organ-specific hair-cell type.")
figure(os.path.join(CLN,"sac_utr_UMAP_sidebyside.png"),
       "Fig. 3. Cell-type annotation of the saccule (left) and utricle (right), shown side by side "
       "(cleaned data, shared colour scheme). Both organs recover the same major cell types.", 6.4)

h2("Cell-type composition of utricle vs. saccule")
para("The two organs showed broadly similar cell-type compositions, dominated by the sensory epithelium (hair "
     "cells and supporting cells) with smaller non-sensory populations (Fig. 4). Replicate-2 libraries carried "
     "more transitional-epithelial and roof/dark cells in both organs, consistent with a dissection/batch "
     "effect rather than a true organ difference.")
figure(os.path.join(CLN,"all_clean_composition_bar.png"),
       "Fig. 4. Cell-type composition per sample, faceted by organ (jointly integrated, cleaned data).", 5.8)

h2("Molecular comparison of the two organs at single-cell depth")
para("We next tested whether the utricle and saccule differ transcriptionally at single-cell depth, using "
     "pseudobulk profiles and edgeR (restricted to protein-coding and lncRNA genes, blocking for batch). The "
     "dominant axis of variation was a technical batch effect between the two sequencing runs; once this was "
     "removed, the organs separated on a secondary axis (Fig. 5), indicating a real but modest tissue signal. "
     "No individual gene reached significance at the organ level or within any cell type, which — with only two "
     "replicates per organ — reflects limited statistical power rather than transcriptional identity between "
     "the organs. This is consistent with the clear tissue separation seen in the much larger bulk RNA-seq "
     "cohort (Fig. 1) and motivates the cell-type-resolved regulatory work outlined below.")
figure(os.path.join(CLN,"pseudobulk_DE","clean_PCA_batch_corrected.png"),
       "Fig. 5. Pseudobulk PCA (cleaned) after removing the batch effect: saccule and utricle separate "
       "on the second principal component.", 4.8)

# ---------------- summary ----------------
h1("Summary and next steps")
bullet(f"mapped {fmt(QTL_TOTAL)} significant independent cis-QTLs across seven molecular phenotypes in the four "
       f"tissue-by-age groups, and a first TWAS nominated {fmt(TWAS_PAIRS)} gene–trait pairs for vestibular function — greatly extending the "
       "regulatory map beyond standard eQTLs.", "Bulk / Pantry: ")
bullet("built a quality-controlled, externally validated 11-cell-type atlas of the vestibular maculae; QC "
       "(SoupX + scDblFinder) recovered neurons, pericytes and endothelial cells and corrected an ambient-RNA "
       "inflation of hair-cell counts.", "Single-nucleus: ")
bullet("the two organs share cell types and are transcriptionally similar at single-cell depth, while the "
       "larger bulk cohort resolves a clear organ difference; the single-cell contrast is modest and "
       "underpowered at two replicates per organ.", "Organ comparison: ")
para("Planned next steps: finalize the TWAS and heritability summaries across all seven phenotypes and "
     "consolidate the candidate-gene shortlist; and use the single-cell atlas to perform cell-type-aware "
     "(deconvolution / interaction) QTL mapping, so that regulatory variants can be assigned to the specific "
     "cell types in which they act. This will link the population-scale regulatory map to the cellular "
     "architecture of the organ and sharpen the biological interpretation of the vestibular-function "
     "associations reported here.", after=4)

doc.save(OUT)
print("SAVED:", OUT)
print("paragraphs:", len(doc.paragraphs), " tables:", len(doc.tables))
