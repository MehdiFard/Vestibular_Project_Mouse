#!/usr/bin/env python3
"""05 - Assemble the results into a PowerPoint deck.

Reads the figures and tables written by scripts 02-04 and produces
report/saccule_vs_utricle_report.pptx.
"""

from pathlib import Path

import pandas as pd
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Emu, Inches, Pt

PROJ = Path("../..")
BASE = PROJ / "analyses" / "05_dge"
FIGS = BASE / "results" / "figures"
TABLES = BASE / "results" / "tables"
DATA = BASE / "data"
OUT = BASE / "report" / "saccule_vs_utricle_report.pptx"

# 16:9
SLIDE_W, SLIDE_H = Inches(13.333), Inches(7.5)

DARK = RGBColor(0x1F, 0x2A, 0x37)
MUTED = RGBColor(0x5B, 0x66, 0x73)
ACCENT = RGBColor(0x76, 0x2A, 0x83)
GREEN = RGBColor(0x1B, 0x78, 0x37)

prs = Presentation()
prs.slide_width, prs.slide_height = SLIDE_W, SLIDE_H
BLANK = prs.slide_layouts[6]


# ----------------------------------------------------------------- primitives
def add_slide():
    return prs.slides.add_slide(BLANK)


def add_text(slide, text, left, top, width, height, size=18, bold=False,
             color=DARK, align=PP_ALIGN.LEFT, italic=False, spacing=1.0):
    box = slide.shapes.add_textbox(left, top, width, height)
    tf = box.text_frame
    tf.word_wrap = True
    lines = text.split("\n") if isinstance(text, str) else text
    for i, line in enumerate(lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.text = line
        p.alignment = align
        p.line_spacing = spacing
        for run in p.runs:
            run.font.size = Pt(size)
            run.font.bold = bold
            run.font.italic = italic
            run.font.color.rgb = color
            run.font.name = "Calibri"
    return box


def add_title(slide, title, subtitle=None):
    add_text(slide, title, Inches(0.55), Inches(0.28), Inches(12.2), Inches(0.7),
             size=28, bold=True)
    if subtitle:
        add_text(slide, subtitle, Inches(0.55), Inches(0.95), Inches(12.2), Inches(0.5),
                 size=14, color=MUTED)


def add_picture_fit(slide, path, left, top, max_w, max_h):
    """Insert a picture scaled to fit inside a box, keeping its aspect ratio."""
    from PIL import Image
    with Image.open(path) as im:
        iw, ih = im.size
    scale = min(max_w / iw, max_h / ih)
    w, h = int(iw * scale), int(ih * scale)
    return slide.shapes.add_picture(str(path), left + int((max_w - w) / 2),
                                    top + int((max_h - h) / 2), width=w, height=h)


def add_table(slide, df, left, top, width, height, font_size=11,
              header_fill=RGBColor(0x2C, 0x3E, 0x50), col_widths=None):
    rows, cols = df.shape[0] + 1, df.shape[1]
    shape = slide.shapes.add_table(rows, cols, left, top, width, height)
    table = shape.table

    if col_widths:
        total = sum(col_widths)
        for j, cw in enumerate(col_widths):
            table.columns[j].width = Emu(int(width * cw / total))

    for j, name in enumerate(df.columns):
        cell = table.cell(0, j)
        cell.text = str(name)
        cell.fill.solid()
        cell.fill.fore_color.rgb = header_fill
        for p in cell.text_frame.paragraphs:
            p.alignment = PP_ALIGN.CENTER
            for r in p.runs:
                r.font.size = Pt(font_size)
                r.font.bold = True
                r.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
                r.font.name = "Calibri"

    for i in range(df.shape[0]):
        for j in range(cols):
            cell = table.cell(i + 1, j)
            cell.text = str(df.iat[i, j])
            cell.fill.solid()
            cell.fill.fore_color.rgb = (RGBColor(0xFF, 0xFF, 0xFF) if i % 2 == 0
                                        else RGBColor(0xF2, 0xF4, 0xF7))
            for p in cell.text_frame.paragraphs:
                p.alignment = PP_ALIGN.CENTER if j > 0 else PP_ALIGN.LEFT
                for r in p.runs:
                    r.font.size = Pt(font_size)
                    r.font.color.rgb = DARK
                    r.font.name = "Calibri"
    return table


def bullets(slide, items, left, top, width, height, size=15, spacing=1.35):
    box = slide.shapes.add_textbox(left, top, width, height)
    tf = box.text_frame
    tf.word_wrap = True
    for i, item in enumerate(items):
        lvl = 0
        if isinstance(item, tuple):
            item, lvl = item
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.text = ("• " if lvl == 0 else "– ") + item
        p.level = lvl
        p.line_spacing = spacing
        p.space_after = Pt(6)
        for r in p.runs:
            r.font.size = Pt(size if lvl == 0 else size - 2)
            r.font.color.rgb = DARK if lvl == 0 else MUTED
            r.font.name = "Calibri"
    return box


def fig_slide(title, image, subtitle=None, notes=None):
    s = add_slide()
    add_title(s, title, subtitle)
    top = Inches(1.5)
    max_h = Inches(5.25) if notes else Inches(5.6)
    add_picture_fit(s, FIGS / image, Inches(0.7), top, Inches(11.9), max_h)
    if notes:
        add_text(s, notes, Inches(0.7), Inches(6.85), Inches(11.9), Inches(0.5),
                 size=12, color=MUTED, italic=True)
    return s


def two_fig_slide(title, img_left, img_right, subtitle=None, notes=None):
    s = add_slide()
    add_title(s, title, subtitle)
    add_picture_fit(s, FIGS / img_left, Inches(0.45), Inches(1.5), Inches(6.1), Inches(5.1))
    add_picture_fit(s, FIGS / img_right, Inches(6.75), Inches(1.5), Inches(6.1), Inches(5.1))
    if notes:
        add_text(s, notes, Inches(0.55), Inches(6.75), Inches(12.2), Inches(0.5),
                 size=12, color=MUTED, italic=True)
    return s


# ----------------------------------------------------------------- load results
samples = pd.read_csv(DATA / "sample_table.tsv", sep="\t")
# master_metadata carries the same 260 libraries plus sequencing QC.
meta = pd.read_csv(TABLES / "master_metadata.tsv", sep="\t")
summ = pd.read_csv(TABLES / "DE_summary.tsv", sep="\t")
overlap = pd.read_csv(TABLES / "contrast_overlap_summary.tsv", sep="\t")
pca_r2 = pd.read_csv(TABLES / "pca_covariate_r2.tsv", sep="\t")
pca_var = pd.read_csv(TABLES / "pca_variance_explained.tsv", sep="\t")
limma_cmp = pd.read_csv(TABLES / "limma_vs_deseq2_summary.tsv", sep="\t")
biotype_sep = pd.read_csv(TABLES / "biotype_separation_comparison.tsv", sep="\t")
best_pc = pd.read_csv(TABLES / "biotype_best_pc_per_factor.tsv", sep="\t")
paired_age = pd.read_csv(TABLES / "paired_age_vs_batch_adjusted.tsv", sep="\t")
sex_eff = pd.read_csv(TABLES / "sex_effect_summary.tsv", sep="\t")
int_ex = pd.read_csv(TABLES / "interaction_examples.tsv", sep="\t")
int_top9 = pd.read_csv(TABLES / "interaction_top9.tsv", sep="\t")
pvu = pd.read_csv(TABLES / "paired_vs_unpaired_summary.tsv", sep="\t")

NICE = {
    "A_utricle_vs_saccule_w8": "Utricle vs Saccule, w8",
    "B_utricle_vs_saccule_w22": "Utricle vs Saccule, w22",
    "C_w22_vs_w8_saccule": "w22 vs w8, Saccule",
    "D_w22_vs_w8_utricle": "w22 vs w8, Utricle",
    "E_organ_by_age_interaction": "Organ x age interaction",
    "S_unpaired_utricle_vs_saccule_w8": "Utricle vs Saccule w8 (unpaired)",
    "S_unpaired_utricle_vs_saccule_w22": "Utricle vs Saccule w22 (unpaired)",
}


def top_table(key, n=12):
    df = pd.read_csv(TABLES / f"DE_{key}.tsv.gz", sep="\t")
    df = df[df["padj"].notna() & (df["padj"] < 0.05) & (df["log2FoldChange"].abs() > 1)]
    df = df.sort_values("padj").head(n)
    out = pd.DataFrame({
        "Gene": df["gene_name"].fillna(df["gene_id"]),
        "Chr": df["chr"].astype(str),
        "log2FC": df["log2FoldChange"].round(2),
        "baseMean": df["baseMean"].round(0).astype(int),
        "padj": df["padj"].apply(lambda v: f"{v:.2e}"),
        "Biotype": df["biotype"].fillna("").str.replace("_", " "),
    })
    return out.reset_index(drop=True)


# ================================================================= 1 title
s = add_slide()
bar = s.shapes.add_textbox(Inches(0), Inches(2.35), SLIDE_W, Inches(0.06))
bar.fill.solid()
add_text(s, "Saccule vs Utricle", Inches(0.9), Inches(2.0), Inches(11.5), Inches(1.0),
         size=44, bold=True, color=ACCENT)
add_text(s, "Expression PCA and differential expression across organ and age",
         Inches(0.9), Inches(3.05), Inches(11.5), Inches(0.6), size=22, color=DARK)
add_text(s,
         "Mouse vestibular otolith organs  |  260 RNA-seq libraries  |  "
         "8 and 22 weeks  |  BXD-style strain panel",
         Inches(0.9), Inches(3.75), Inches(11.5), Inches(0.6), size=15, color=MUTED)
add_text(s, "Counts: kallisto (Pantry phenotyping run) summarised to genes with tximport\n"
            "Analysis: DESeq2  |  Generated from Pantry/sac_vs_utr",
         Inches(0.9), Inches(5.6), Inches(11.5), Inches(0.9), size=12, color=MUTED)

# ================================================================= 2 design
s = add_slide()
add_title(s, "Study design", "What went into the analysis, and how the factors are arranged")

# --- libraries, pairing and sex, per organ/age group -------------------------
# "In a pair" = the animal this library came from contributed both organs, so
# the library can be used in the paired organ contrast.
def group_row(label_organ, label_age, d):
    paired = int(d["has_organ_pair"].sum())
    return [label_organ, label_age, len(d), paired, len(d) - paired,
            int((d["sex_final"] == "F").sum()), int((d["sex_final"] == "M").sum())]

rows = [group_row(o, a, meta[(meta["organ"] == o) & (meta["age"] == a)])
        for o in ["Saccule", "Utricle"] for a in ["w8", "w22"]]
rows.append(group_row("All", "", meta))
counts_df = pd.DataFrame(rows, columns=["Organ", "Age", "Libraries", "In a pair",
                                        "Unpaired", "Female", "Male"])

add_text(s, "Libraries, pairing and sex per group", Inches(0.55), Inches(1.38),
         Inches(7.2), Inches(0.35), size=14, bold=True)
add_table(s, counts_df, Inches(0.55), Inches(1.78), Inches(7.05), Inches(1.85),
          font_size=11, col_widths=[15, 9, 15, 14, 13, 12, 10])
add_text(s,
         "A pair is one animal (submission workbook), both organs present. "
         "70 animals pair completely at w8 and 51 at w22,\n"
         "which is the 242 libraries the organ contrasts use. Within any one group "
         "each strain appears exactly once.",
         Inches(0.55), Inches(3.72), Inches(7.2), Inches(0.5), size=10.5, color=MUTED,
         italic=True, spacing=1.15)

# --- sequencing batch, with sex and the QC that separates batch 4 ------------
brows = []
for b in sorted(meta["batch"].unique()):
    d = meta[meta["batch"] == b]
    brows.append([int(b), int((d["age"] == "w8").sum()), int((d["age"] == "w22").sum()),
                  int((d["sex_final"] == "F").sum()), int((d["sex_final"] == "M").sum()),
                  round(d["n_reads_processed"].median() / 1e6, 1),
                  round(d["pct_pseudoaligned"].median(), 1)])
batch_tab = pd.DataFrame(brows, columns=["Batch", "w8", "w22", "Female", "Male",
                                         "Median reads (M)", "Mapped %"])

add_text(s, "Sequencing batch: age, sex and quality", Inches(0.55), Inches(4.35),
         Inches(7.2), Inches(0.35), size=14, bold=True)
add_table(s, batch_tab, Inches(0.55), Inches(4.75), Inches(7.05), Inches(1.55),
          font_size=11, col_widths=[11, 9, 9, 13, 11, 22, 15])
add_text(s,
         "Organs are balanced inside every batch (44/44, 38/38, 35/37, 11/13), so no "
         "batch favours one organ.\n"
         "Batch 4 is the later 'oldmouse' run: it adds 10 w22 pairs but is the "
         "shallowest and worst-mapping data in the set.",
         Inches(0.55), Inches(6.38), Inches(7.2), Inches(0.5), size=10.5, color=MUTED,
         italic=True, spacing=1.15)

bullets(s, [
    "Both organs came from the same animal, so strain + age + sex identifies a "
    "saccule/utricle pair. That is why strain goes into the organ models.",
    "37 strains were sampled at both ages (32 saccule, 36 utricle), so the age "
    "contrast is mostly between animals rather than within. A within-strain "
    "version is shown later as a check.",
    "Batch and age overlap: batch 1 is w8-only, batch 4 is w22-only. The age "
    "effect is identified from batches 2 and 3, which span both.",
    "Every w8 library is female, so sex cannot be separated from age across all "
    "libraries. Sex is known for all 260 - recovered from Xist and the Y genes "
    "where the sheet had none - so a female-only age contrast is possible.",
    "Strain comes from RNA-to-genotype concordance, not the metadata sheet - the "
    "same source Pantry used to attach genotypes. It covers all 260 libraries and "
    "is 99.1% self-consistent across organ pairs.",
], Inches(7.95), Inches(1.38), Inches(4.95), Inches(5.6), size=13, spacing=1.25)

# ================================================================= 3 methods
s = add_slide()
add_title(s, "How the numbers were produced", "Each step and why it was chosen")
bullets(s, [
    "Counts: kallisto est_counts from the existing Pantry phenotyping run, "
    "summarised transcript to gene with tximport (average transcript length "
    "carried as a DESeq2 offset).",
    ("Uses the same quantification as everything else in the Pantry project, so "
     "the DE results and the QTL results describe the same numbers.", 1),
    "Biotype filter: only protein-coding and lincRNA genes are analysed. "
    "55,401 annotated genes drop to 27,488.",
    ("Pseudogenes, TEC, antisense, miRNA and snoRNA are removed before testing, so "
     "they never enter the multiple-testing burden.", 1),
    "Expression filter: a gene is kept if it reaches 10 counts in at least as many "
    "libraries as the smallest group. 18,273 genes pass across all 260; each "
    "contrast re-applies both filters on its own subset and tests ~17,000.",
    ("Stricter than a non-zero count filter, which lets through genes that are "
     "near zero in every sample and only inflates the multiple-testing burden.", 1),
    "Organ contrasts use ~ strain + organ, fitted separately at each age.",
    ("Pairing removes between-animal and genetic-background variation, which is "
     "the largest nuisance term in a strain panel.", 1),
    "Age contrasts use ~ batch + age, fitted separately in each organ.",
    "Interaction model ~ batch + organ + age + organ:age asks directly whether "
    "the organ difference itself changes with age.",
    "Fold changes are shrunken with ashr; significance is FDR < 0.05 together "
    "with |log2FC| > 1.",
], Inches(0.7), Inches(1.5), Inches(11.9), Inches(5.5), size=14.5)

# ================================================================= 4-8 PCA
# Pull the numbers quoted in the notes straight from the tables, so they cannot
# drift away from the figures they describe.
def var_of(pc):
    return float(pca_var.loc[pca_var["PC"] == pc, "percent_variance"].iloc[0])


def r2_of(pc, cov):
    m = (pca_r2["PC"] == pc) & (pca_r2["covariate"] == cov)
    return float(pca_r2.loc[m, "r2"].iloc[0])


fig_slide("Expression PCA: organ separates on PC1",
          "01_pca_pc1_pc2_organ_age.png",
          "Colour = organ, shape = age, as requested. Variance-stabilised expression, top 2000 most variable genes.",
          f"PC1 ({var_of('PC1')}%) splits saccule from utricle - organ alone explains "
          f"R2 = {r2_of('PC1', 'organ'):.2f} of it. PC2 ({var_of('PC2')}%) is mostly strain "
          f"(R2 = {r2_of('PC2', 'strain'):.2f}); age sits on PC3, not PC2. A handful of points "
          "fall inside the opposite organ's cluster - see the sample-check slide.")

two_fig_slide("The same PCA by group, and by sex",
              "02_pca_pc1_pc2_ellipses.png", "08_pca_by_sex.png",
              "Left: colour = organ, shape = age, with 95% ellipses per organ-age group. "
              "Right: the same coordinates coloured by sex, shape = organ.",
              f"Sex is known for all 260 libraries, recovered from Xist and the Y genes. PC1 is organ "
              f"and nothing else - sex explains R2 = {r2_of('PC1', 'sex'):.3f} of it - so the organ "
              "contrasts cannot be confounded by sex. The right panel also shows the design problem: "
              "every 8-week library is female, so across the whole set 'sex' and 'age' colour the same points.")

fig_slide("Where sex does show up: PC3", "09_pca_pc1_pc3_by_sex.png",
          "Same libraries, plotted on PC1 against PC3. Colour = sex, shape = organ.",
          f"PC3 carries age (R2 = {r2_of('PC3', 'age'):.2f}), batch ({r2_of('PC3', 'batch'):.2f}) and "
          f"sex ({r2_of('PC3', 'sex'):.2f}) together. Restricting to 22 weeks removes age but not batch: "
          "there, batch 4 is entirely female and batch 3 is 48 male to 13 female. This is why sex is "
          "handled by restriction - a female-only age contrast - rather than as a covariate.")

two_fig_slide("How much variance, and what drives it",
              "03_pca_scree.png", "06_pca_covariate_r2.png",
              "Left: variance per PC. Right: R-squared of each PC against each known covariate.",
              f"Organ explains R2 = {r2_of('PC1', 'organ'):.2f} of PC1 and almost nothing else. Age, batch "
              f"and sex all load on PC3 ({r2_of('PC3', 'age'):.2f}, {r2_of('PC3', 'batch'):.2f}, "
              f"{r2_of('PC3', 'sex'):.2f}) and cannot be fully separated in this design, though age is the "
              "strongest of the three. Strain drives PC2 and PC4-PC6, as expected for a strain panel.")

two_fig_slide("Batch effect on the same coordinates",
              "04_pca_by_batch.png", "05_pca_batch_corrected.png",
              "Left: coloured by sequencing batch. Right: after limma::removeBatchEffect, organ and age protected.",
              "The organ split on PC1 is unaffected by batch correction, so it is a real biological signal. "
              "The age separation is partly batch.")

# ================================================================= 8 DE summary
s = add_slide()
add_title(s, "Differential expression: headline numbers",
          "FDR < 0.05 and |log2FC| > 1")
main = summ[summ["analysis"].isin(NICE)].copy()
main["Contrast"] = main["analysis"].map(NICE)
tab = main[["Contrast", "model", "n_libraries", "genes_tested",
            "sig_fdr", "sig_fdr_lfc", "up", "down"]].copy()
tab.columns = ["Contrast", "Model", "Libraries", "Genes tested",
               "FDR<0.05", "FDR + |LFC|>1", "Positive", "Negative"]
tab["Model"] = tab["Model"].str.replace("~ batch + organ + age + organ:age",
                                        "~ batch + organ*age", regex=False)
add_table(s, tab, Inches(0.5), Inches(1.55), Inches(12.3), Inches(2.4), font_size=11,
          col_widths=[3.0, 2.3, 1.0, 1.2, 1.1, 1.3, 1.0, 1.0])
add_text(s,
         "Positive log2FC means higher in Utricle (organ contrasts) or higher at w22 (age contrasts). "
         "For the interaction, positive means the organ gap widens with age.",
         Inches(0.5), Inches(4.25), Inches(12.3), Inches(0.5), size=12, color=MUTED, italic=True)
add_picture_fit(s, FIGS / "15_significant_gene_counts.png",
                Inches(1.9), Inches(4.8), Inches(9.5), Inches(2.5))

# ================================================================= 9-10 organ
fig_slide("Utricle vs Saccule at 8 weeks", "10_volcano_utricle_vs_saccule_w8.png",
          "Paired within strain (~ strain + organ), 70 strain pairs, 140 libraries.")
n_w22 = int(summ.loc[summ["analysis"] == "B_utricle_vs_saccule_w22", "n_libraries"].iloc[0])
fig_slide("Utricle vs Saccule at 22 weeks", "11_volcano_utricle_vs_saccule_w22.png",
          f"Paired within strain (~ strain + organ), {n_w22 // 2} animal pairs, {n_w22} libraries.")

s = add_slide()
add_title(s, "Top organ-discriminating genes",
          "Ranked by adjusted p in the w8 contrast")
add_table(s, top_table("A_utricle_vs_saccule_w8", 11), Inches(0.5), Inches(1.5),
          Inches(6.1), Inches(4.5), font_size=10.5,
          col_widths=[1.7, 0.6, 0.9, 1.0, 1.1, 1.5])
add_text(s, "w8", Inches(0.5), Inches(1.15), Inches(3.0), Inches(0.3), size=14,
         bold=True, color=GREEN)
add_table(s, top_table("B_utricle_vs_saccule_w22", 11), Inches(6.8), Inches(1.5),
          Inches(6.1), Inches(4.5), font_size=10.5,
          col_widths=[1.7, 0.6, 0.9, 1.0, 1.1, 1.5])
add_text(s, "w22", Inches(6.8), Inches(1.15), Inches(3.0), Inches(0.3), size=14,
         bold=True, color=ACCENT)
add_text(s, "Positive log2FC = higher in Utricle. Full tables in results/tables/.",
         Inches(0.5), Inches(6.3), Inches(12.3), Inches(0.4), size=12,
         color=MUTED, italic=True)

fig_slide("Is the organ difference the same at both ages?",
          "16_concordance_organ_w8_vs_w22.png",
          "log2FC at w8 against log2FC at w22, over the 16,834 genes tested in both.",
          "Answer: yes. r = 0.92 with a slope of 0.98, 593 genes significant in both and 100% of them "
          "in the same direction - not one sign flip. Every mouse is either w8 or w22, so these are two "
          "disjoint sets of animals in different batches: this is a genuine replication of the organ "
          "signature, not the same data twice. The off-diagonal genes are mostly threshold effects, and "
          "contrast E puts a number on the real exceptions: 122 of 18,273.")

# ================================================================= 11-13 age
fig_slide("22 weeks vs 8 weeks in Saccule", "12_volcano_w22_vs_w8_saccule.png",
          "~ batch + age over all 128 saccule libraries.")
fig_slide("22 weeks vs 8 weeks in Utricle", "13_volcano_w22_vs_w8_utricle.png",
          "~ batch + age over all 132 utricle libraries.")
fig_slide("Is the ageing response the same in both organs?",
          "17_concordance_age_saccule_vs_utricle.png",
          "log2FC in saccule against log2FC in utricle, over genes tested in both.")

s = add_slide()
add_title(s, "Top age-associated genes", "Ranked by adjusted p within each organ")
add_table(s, top_table("C_w22_vs_w8_saccule", 11), Inches(0.5), Inches(1.5),
          Inches(6.1), Inches(4.5), font_size=10.5,
          col_widths=[1.7, 0.6, 0.9, 1.0, 1.1, 1.5])
add_text(s, "Saccule", Inches(0.5), Inches(1.15), Inches(3.0), Inches(0.3), size=14,
         bold=True, color=GREEN)
add_table(s, top_table("D_w22_vs_w8_utricle", 11), Inches(6.8), Inches(1.5),
          Inches(6.1), Inches(4.5), font_size=10.5,
          col_widths=[1.7, 0.6, 0.9, 1.0, 1.1, 1.5])
add_text(s, "Utricle", Inches(6.8), Inches(1.15), Inches(3.0), Inches(0.3), size=14,
         bold=True, color=ACCENT)
add_text(s, "Positive log2FC = higher at 22 weeks. Read these alongside the batch caveat.",
         Inches(0.5), Inches(6.3), Inches(12.3), Inches(0.4), size=12,
         color=MUTED, italic=True)

# ================================================================= 14 overlap
s = add_slide()
add_title(s, "Shared and condition-specific signal",
          "How much the two organ contrasts, and the two age contrasts, agree")
add_picture_fit(s, FIGS / "18_shared_vs_specific.png",
                Inches(0.6), Inches(1.45), Inches(12.1), Inches(3.4))
ov = overlap.copy()
ov.columns = ["Comparison", "r (all genes)", "r (shared sig.)",
              "Sig. in both", "Sig. only in first", "Sig. only in second"]
for c in ["r (all genes)", "r (shared sig.)"]:
    ov[c] = ov[c].round(3)
add_table(s, ov, Inches(0.9), Inches(5.15), Inches(11.5), Inches(1.2), font_size=11.5,
          col_widths=[2.6, 1.3, 1.3, 1.2, 1.4, 1.4])

# ================================================================= 15 interaction
fig_slide("What the organ x age interaction asks", "31_interaction_illustration.png",
          "Six genes, one for each shape the interaction can take. Lines are the model's fitted means.",
          "A, B, C and D ask 'is this gene different?'. The interaction asks whether the ORGAN DIFFERENCE "
          "itself changes with age. Gpnmb is the control - a large organ difference that does not move "
          "with age, so it scores zero. The other five are significant.")

s = add_slide()
add_title(s, "The same six genes as numbers",
          "interaction = gap(w22) - gap(w8), where gap = log2(Utricle / Saccule)")
add_table(s, int_ex, Inches(0.45), Inches(1.5), Inches(12.45), Inches(2.6),
          font_size=11, col_widths=[1.1, 1.5, 1.5, 1.7, 1.0, 4.4])
bullets(s, [
    "Kcne1 worked through: at w8 the utricle carries 2^3.99 = 16x more Kcne1 than the saccule. "
    "By w22 it carries 2^7.18 = 145x more. The gap grew by 3.19 log2 units, which is the interaction.",
    ("Both organs lose Kcne1 with age - the saccule just loses far more of it. The interaction says "
     "nothing about the direction of the age change, only about the gap between the organs.", 1),
    "Gpnmb is the control: a 9-fold organ difference (gap +3.23) that is the same at both ages "
    "(+3.09), so the interaction is -0.14 and nowhere near significant.",
    "Only 122 of 18,273 genes have a significant interaction, against 756 and 863 for the organ "
    "contrasts themselves. The organ identity of these two tissues is mostly age-stable.",
], Inches(0.5), Inches(4.35), Inches(12.3), Inches(2.6), size=13, spacing=1.2)

s = add_slide()
add_title(s, "The nine most significant interaction genes",
          "Ranked by adjusted p among the 122 genes passing FDR < 0.05 and |log2FC| > 1")
add_picture_fit(s, FIGS / "32_interaction_top9.png",
                Inches(0.4), Inches(1.45), Inches(7.5), Inches(5.3))
add_table(s, int_top9, Inches(8.15), Inches(1.6), Inches(4.75), Inches(3.1),
          font_size=10, col_widths=[1.3, 1.1, 0.9, 1.0, 1.1, 1.1])
bullets(s, [
    "Cma1, Mcpt4 and Tpsb2 are mast-cell proteases and move together: utricle-enriched at "
    "w8, saccule-enriched at w22.",
    "Read those three with care - their mean counts are 9-11, so they top the ranking on a "
    "large fold change in very few reads.",
    "Trpm1, Gchfr, Fam163b and Jchain are the best-supported entries here.",
], Inches(8.15), Inches(4.95), Inches(4.8), Inches(2.2), size=11.5, spacing=1.15)

fig_slide("Organ x age interaction",
          "14_volcano_interaction.png",
          "Genes where the Utricle-Saccule difference itself changes between w8 and w22.")
fig_slide("Genes with the strongest organ x age interaction",
          "20_top_interaction_genes.png",
          "Normalised expression by organ and age.")

s = add_slide()
add_title(s, "Top interaction genes",
          "The organ difference in these genes is not constant across age")
add_table(s, top_table("E_organ_by_age_interaction", 12), Inches(3.2), Inches(1.55),
          Inches(6.9), Inches(4.8), font_size=11,
          col_widths=[1.7, 0.6, 0.9, 1.0, 1.1, 1.5])
add_text(s, "log2FC here is the change in the Utricle-minus-Saccule difference "
            "going from w8 to w22.",
         Inches(0.6), Inches(6.5), Inches(12.1), Inches(0.4), size=12,
         color=MUTED, italic=True)

# ================================================================= 16 top genes
fig_slide("Expression of the top organ-marker genes",
          "19_top_organ_genes.png",
          "Normalised counts split by organ and age.")

# ================================================================= 17 paired
fig_slide("What pairing by strain bought", "21_paired_vs_unpaired.png",
          "Same w8 organ contrast fitted with and without strain in the model.",
          "Pairing adds about 1,400 genes at FDR < 0.05 (10,923 vs 9,537). It does NOT add genes past "
          "|log2FC| > 1 - the unpaired run reports slightly more of those. The next slide explains why "
          "that is a symptom of noise, not of extra signal.")

s = add_slide()
add_title(s, "Why the unpaired run reports MORE genes",
          "The headline count is an effect-size count, not a power count")
pvu_tab = pvu.copy()
pvu_tab["model"] = pvu_tab["model"].str.replace("~ ", "", regex=False)
pvu_tab = pvu_tab[["comparison", "model", "n_libraries", "residual_df",
                   "sig_fdr", "lfc_over_1", "both", "median_se"]]
pvu_tab.columns = ["Age", "Model", "Libraries", "Resid. df", "FDR < 0.05",
                   "|log2FC| > 1", "Both", "Median SE"]
add_table(s, pvu_tab, Inches(0.45), Inches(1.45), Inches(12.45), Inches(2.3),
          font_size=10.5, col_widths=[0.8, 3.4, 1.2, 1.1, 1.3, 1.4, 1.0, 1.2])
bullets(s, [
    "Read the last two columns first: 'Both' equals '|log2FC| > 1' on every row. Every gene that "
    "clears two-fold also clears FDR, in both models - so the headline number is counting effect "
    "sizes and the FDR test is doing no filtering at all here.",
    "Pairing's real benefit is in the FDR column. On the SAME libraries, strain adds 1,475 genes at "
    "w8 (10,923 vs 9,448) and 1,503 at w22 (10,969 vs 9,466), despite spending 70 degrees of freedom "
    "on strain coefficients.",
    "Dropping strain makes each fold change about 20% less precise (median SE ratio 1.20 at w8, "
    "1.17 at w22). A noisier estimate overshoots and undershoots equally often, but far more genes "
    "sit just below |log2FC| = 1 than just above, so noise promotes more than it demotes: at w8, "
    "76 genes cross the line and only 36 fall back; at w22, 95 in and 44 out.",
    "So the extra unpaired 'hits' are the threshold catching noise. They are not extra discoveries.",
], Inches(0.5), Inches(4.0), Inches(12.3), Inches(3.0), size=13, spacing=1.2)

# ================================================================= biotype filter
s = add_slide()
add_title(s, "Why only protein-coding and lincRNA genes",
          "Keeping all biotypes looks like it separates age better - but what separates is batch")
add_picture_fit(s, FIGS / "22_biotype_pca_comparison.png",
                Inches(0.4), Inches(1.35), Inches(12.5), Inches(3.3))

bt = biotype_sep.copy()
bt_tab = pd.DataFrame({
    "Gene set": bt["gene_set"],
    "Genes": bt["n_genes"],
    "PC1 % var": bt["PC1_pct"].round(1),
    "Organ separation (PC1)": bt["organ_sep_PC1"].round(2),
    "PC1 R² organ": bt["PC1_r2_organ"].round(3),
    "PC2 R² age": bt["PC2_r2_age"].round(3),
    "PC2 R² batch": bt["PC2_r2_batch"].round(3),
})
add_table(s, bt_tab, Inches(0.6), Inches(4.7), Inches(12.1), Inches(1.0), font_size=11,
          col_widths=[2.6, 0.9, 1.1, 1.7, 1.2, 1.1, 1.2])
add_text(s,
         "Organ separation is unchanged (standardised distance 4.50 → 4.71). With all biotypes PC2 tracks "
         "BATCH (R²=0.78) more than age (R²=0.65), and batch overlaps age here, so the batch axis reads as age "
         "separation. After the filter age has not gone - it moves to PC3, where it is both stronger (R²=0.70) "
         "and now beats batch (0.55) instead of losing to it. The non-coding biotypes track batch (mean R²=0.19) "
         "far more than organ (0.035): they are low-count, multi-mapping features sensitive to library prep.",
         Inches(0.6), Inches(5.95), Inches(12.1), Inches(1.2), size=12, color=MUTED, italic=True)

# ================================================================= limma cross-check
s = add_slide()
add_title(s, "Cross-check: DESeq2 vs limma-voom",
          "Same libraries, same filters, same designs - a different statistical engine")
add_picture_fit(s, FIGS / "23_limma_vs_deseq2_lfc.png",
                Inches(0.55), Inches(1.4), Inches(12.2), Inches(3.9))

lc = limma_cmp.copy()
lc_tab = pd.DataFrame({
    "Contrast": lc["contrast"],
    "DESeq2": lc["deseq2_sig"],
    "limma": lc["limma_sig"],
    "Shared": lc["shared"],
    "% of DESeq2 recovered": lc["pct_of_deseq2_recovered"],
    "Jaccard": lc["jaccard"],
    "r (log2FC)": lc["r_lfc"],
})
add_table(s, lc_tab, Inches(1.4), Inches(5.5), Inches(10.5), Inches(1.6), font_size=10.5,
          col_widths=[2.6, 0.9, 0.9, 0.9, 1.7, 1.0, 1.1])

fig_slide("Which genes each engine finds", "24_limma_vs_deseq2_overlap.png",
          "Significant at FDR < 0.05 and |log2FC| > 1 in DESeq2 only, limma only, or both.",
          "The organ contrasts agree almost completely (96-97% of DESeq2 hits recovered, r = 0.98). The age "
          "contrasts agree on 90-93%, with limma calling roughly twice as many overall. The interaction is the "
          "one place they diverge (37% shared) - interaction terms carry the least information, so the two "
          "shrinkage schemes pull apart. Treat the interaction list as the least settled result here.")

# ================================================================= what pairing means
s = add_slide()
add_title(s, "What \u201cstrain identifies a pair\u201d actually means",
          "Each strain contributed one saccule and one utricle, from the same animal, at the same age")
add_picture_fit(s, FIGS / "27_pairing_illustration.png",
                Inches(0.5), Inches(1.4), Inches(12.3), Inches(4.0))
add_text(s,
         "Every grey line is one strain. The vertical spread of the green points is genetic background: strains "
         "differ from each other by several-fold. The line slopes are the organ effect. Averaged per-gene, strain "
         "explains 62% of the variation among the w8 libraries while organ explains only 10% \u2013 so the thing we "
         "are measuring is smaller than the thing we are blocking out (adjusted R2 0.10 against 0.24; the raw "
         "0.10-against-0.62 comparison is inflated by the animal factor's 69 parameters). An unpaired model has to treat "
         "that vertical spread as noise; putting strain in the model cancels it, leaving only the slopes.",
         Inches(0.5), Inches(5.6), Inches(12.3), Inches(1.4), size=12.5, color=MUTED, italic=True)

# ================================================================= paired age
s = add_slide()
add_title(s, "Putting the both-age strains to work",
          "32-36 strains were sampled at 8 and 22 weeks, so their age effect can be measured within strain")
add_picture_fit(s, FIGS / "28_paired_age_check.png",
                Inches(0.55), Inches(1.35), Inches(7.4), Inches(4.6))

pa = paired_age.copy()
pa_tab = pd.DataFrame({
    "Organ": pa["organ"],
    "Strains": pa["n_strains"],
    "Main analysis": pa["batch_adjusted_sig"],
    "Within-strain": pa["paired_sig"],
    "Shared": pa["shared"],
    "Same direction": pa["same_direction_pct"].astype(str) + "%",
})
add_table(s, pa_tab, Inches(8.15), Inches(1.9), Inches(4.7), Inches(1.1), font_size=10.5,
          col_widths=[1.2, 0.9, 1.3, 1.3, 0.9, 1.2])
bullets(s, [
    "Model is ~ strain + batch + age, and batch IS estimable here because "
    "strains differ in which batches their two libraries fell into.",
    "64-72 libraries beat 128-132: more hits from half the data.",
    "99-100% of the main analysis hits move the same direction here.",
], Inches(8.15), Inches(3.3), Inches(4.8), Inches(2.6), size=12)
add_text(s,
         "This is the strongest evidence the age results are real: within the same genetic line, adjusted for "
         "batch, the same genes move the same way.",
         Inches(0.55), Inches(6.2), Inches(12.3), Inches(0.8), size=12, color=MUTED, italic=True)

# ================================================================= sex
s = add_slide()
add_title(s, "Sex: recovered, and what it cost",
          "52 libraries had no sex recorded; the expression data settles it directly")
add_picture_fit(s, FIGS / "30_sex_inference.png",
                Inches(0.4), Inches(1.3), Inches(6.4), Inches(4.3))
add_picture_fit(s, FIGS / "31_age_effect_females_only.png",
                Inches(6.9), Inches(1.3), Inches(6.1), Inches(4.3))

se = sex_eff.copy()
se_tab = pd.DataFrame({
    "Organ": se["organ"],
    "Age hits, all libs": se["age_hits_all_samples"],
    "Age hits, females only": se["age_hits_females_only"],
    "Shared": se["shared"],
    "% recovered": se["pct_recovered"],
    "Sex hits in w22": se["sex_hits_w22"],
    "Age hits that are also sex hits": se["age_hits_that_are_also_sex_hits"]})
add_table(s, se_tab, Inches(0.8), Inches(5.75), Inches(11.8), Inches(1.0), font_size=10.5,
          col_widths=[1.1, 1.5, 1.7, 0.9, 1.2, 1.4, 2.2])
add_text(s,
         "Xist against the Y genes separates the sexes completely - 100% concordant with all 208 libraries whose "
         "sex WAS recorded. Only 31 and 8 genes respond to sex at all, and just 11 and 3 of the age hits are among "
         "them. The female-only age contrast recovers 90-95% of the original hits and finds more besides, because "
         "dropping the males removes heterogeneity from the w22 group.",
         Inches(0.5), Inches(6.85), Inches(12.4), Inches(0.6), size=11, color=MUTED, italic=True)

# ================================================================= PCA method notes
s = add_slide()
add_title(s, "Two notes on how the PCA was built",
          "What the transform does, and why the gene count for the PCA does not matter")
add_picture_fit(s, FIGS / "25_vst_variance_stabilisation.png",
                Inches(0.4), Inches(1.4), Inches(6.2), Inches(4.5))
add_picture_fit(s, FIGS / "26_pca_ntop_sensitivity.png",
                Inches(6.75), Inches(1.4), Inches(6.2), Inches(4.5))
add_text(s,
         "Left - VST (variance-stabilising transform) puts every gene on a log2-like scale where spread no longer "
         "depends on expression level. Raw counts let the few highest-expressed genes dominate; log2(count+1) "
         "overcorrects and makes the lowest-expressed genes the noisiest. VST is the flat one.\n"
         "Right - the PCA uses the 2000 most variable genes, the usual convention. PC1 tracks organ just as well "
         "at 250 genes (R²=0.82) or all 18,273 (R²=0.81) as at 2000 (R²=0.85), so nothing here hinges on that choice. "
         "This applies to the PCA only: every differential expression test uses all ~17,000 filtered genes.",
         Inches(0.5), Inches(6.0), Inches(12.4), Inches(1.3), size=11.5, color=MUTED, italic=True)

# ================================================================= 18 caveats
s = add_slide()
add_title(s, "Caveats and what to do next", "Read the results with these in mind")
bullets(s, [
    "Age and batch are partly confounded. Batch 1 is w8-only and batch 4 is "
    "w22-only, so every age result carries some batch signal. The organ results "
    "do not have this problem - organ is balanced inside every batch.",
    "Sex is now known for all 260 libraries. All w8 libraries are female, so the "
    "all-library age contrast still cannot separate age from sex - but the "
    "female-only version can, and it reproduces 90-95% of those hits.",
    "A few libraries sit inside the opposite organ's PC1 cluster. These are worth "
    "checking against genotype before the results are finalised.",
    "The organ and age results reproduce under limma-voom; the organ x age "
    "interaction list is the least reproducible and should be treated as provisional.",
    "Worth adding next:",
    ("Pathway / GO enrichment on each significant list - the old project already "
     "has a GSEA script that can be pointed at these tables.", 1),
    ("The same four contrasts on the other Pantry modalities (splicing, RNA "
     "editing, alt TSS/polyA, stability) - differential splicing between the two "
     "organs is a natural companion to this.", 1),
    ("Overlap between the organ-specific genes here and the cis-QTL genes from "
     "Pheast, to see whether organ-biased genes are also more often regulated.", 1),
    ("Recovering strain IDs for the 24 later-run w22 libraries - they are "
     "currently dropped from the paired w22 organ contrast, which loses 27 of "
     "113 libraries.", 1),
], Inches(0.7), Inches(1.45), Inches(11.9), Inches(5.7), size=14)

# ================================================================= 19 files
s = add_slide()
add_title(s, "Where everything lives", "Pantry/sac_vs_utr/")
files_df = pd.DataFrame({
    "Path": ["scripts/00_config.R", "scripts/01_build_inputs.R", "scripts/02_pca.R",
             "scripts/03_differential_expression.R", "scripts/04_figures.R",
             "scripts/05_build_report.py", "scripts/06_batch_correction_check.R",
             "scripts/07_biotype_separation_check.R", "scripts/08_limma_crosscheck.R",
             "scripts/run_all.sh",
             "data/", "results/tables/", "results/figures/", "report/"],
    "Contents": ["Shared paths, thresholds and colours",
                 "Sample table, gene annotation, tximport count matrix",
                 "PCA, scree, covariate R-squared, batch check",
                 "The five DESeq2 analyses plus sensitivity checks",
                 "Volcano, concordance and top-gene figures",
                 "Builds this deck",
                 "Why batch stays in the design instead of being corrected out",
                 "Whether dropping the other biotypes costs any real separation",
                 "limma-voom rerun of every contrast, and the agreement tables",
                 "Runs the whole thing end to end",
                 "Count matrix, sample table, cached DESeq2 objects",
                 "Full and significant-only result tables (tsv)",
                 "All figures as 300 dpi PNG",
                 "This PowerPoint"],
})
add_table(s, files_df, Inches(1.1), Inches(1.6), Inches(11.1), Inches(5.0),
          font_size=12, col_widths=[3.0, 6.0])

# ----------------------------------------------------------------- save
OUT.parent.mkdir(parents=True, exist_ok=True)
prs.save(OUT)
print(f"Wrote {OUT}")
print(f"{len(prs.slides.__iter__.__self__._sldIdLst)} slides")
