#!/usr/bin/env python3
"""06 - Deck for the 291-library rebuild, against the 260-library original.
"""

from pathlib import Path

import pandas as pd
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Emu, Inches, Pt

PA = Path("../../analyses/05_dge/rebuild_291")
FIGS = PA / "results" / "figures"
TABLES = PA / "results" / "tables"
OUT = PA / "report" / "rebuild_291_report.pptx"

SLIDE_W, SLIDE_H = Inches(13.333), Inches(7.5)
DARK = RGBColor(0x1F, 0x2A, 0x37)
MUTED = RGBColor(0x5B, 0x66, 0x73)
ACCENT = RGBColor(0x76, 0x2A, 0x83)

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
    for i, line in enumerate(text.split("\n") if isinstance(text, str) else text):
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
             size=27, bold=True)
    if subtitle:
        add_text(slide, subtitle, Inches(0.55), Inches(0.95), Inches(12.2), Inches(0.5),
                 size=14, color=MUTED)


def add_picture_fit(slide, path, left, top, max_w, max_h):
    from PIL import Image
    with Image.open(path) as im:
        iw, ih = im.size
    scale = min(max_w / iw, max_h / ih)
    w, h = int(iw * scale), int(ih * scale)
    return slide.shapes.add_picture(str(path), left + int((max_w - w) / 2),
                                    top + int((max_h - h) / 2), width=w, height=h)


def add_table(slide, df, left, top, width, height, font_size=11, col_widths=None):
    rows, cols = df.shape[0] + 1, df.shape[1]
    table = slide.shapes.add_table(rows, cols, left, top, width, height).table
    if col_widths:
        total = sum(col_widths)
        for j, cw in enumerate(col_widths):
            table.columns[j].width = Emu(int(width * cw / total))
    for j, name in enumerate(df.columns):
        cell = table.cell(0, j)
        cell.text = str(name)
        cell.fill.solid()
        cell.fill.fore_color.rgb = RGBColor(0x2C, 0x3E, 0x50)
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


def bullets(slide, items, left, top, width, height, size=14, spacing=1.25):
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
    add_picture_fit(s, FIGS / image, Inches(0.7), Inches(1.5), Inches(11.9),
                    Inches(5.25) if notes else Inches(5.6))
    if notes:
        add_text(s, notes, Inches(0.7), Inches(6.85), Inches(11.9), Inches(0.5),
                 size=12, color=MUTED, italic=True)
    return s


def two_fig_slide(title, left_img, right_img, subtitle=None, notes=None):
    s = add_slide()
    add_title(s, title, subtitle)
    add_picture_fit(s, FIGS / left_img, Inches(0.45), Inches(1.5), Inches(6.1), Inches(5.1))
    add_picture_fit(s, FIGS / right_img, Inches(6.75), Inches(1.5), Inches(6.1), Inches(5.1))
    if notes:
        add_text(s, notes, Inches(0.55), Inches(6.75), Inches(12.2), Inches(0.5),
                 size=12, color=MUTED, italic=True)
    return s



# ----------------------------------------------------------------- data
V1 = Path("../../analyses/05_dge")
T1 = V1 / "results" / "tables"
samples = pd.read_csv(PA / "data" / "sample_table.tsv", sep="\t")
s2 = pd.read_csv(TABLES / "DE_summary.tsv", sep="\t")
s1 = pd.read_csv(T1 / "DE_summary.tsv", sep="\t")
cmp = pd.read_csv(TABLES / "v1_vs_v2_summary.tsv", sep="\t")
conc = pd.read_csv(TABLES / "v1_vs_v2_concordance.tsv", sep="\t")
g2 = pd.read_csv(TABLES / "gsea_summary.tsv", sep="\t")
hm = pd.read_csv(TABLES / "gsea_v1_vs_v2_hallmark.tsv", sep="\t")
r2 = pd.read_csv(TABLES / "pca_covariate_r2.tsv", sep="\t")
var = pd.read_csv(TABLES / "pca_variance_explained.tsv", sep="\t")

NICE = {"A_utricle_vs_saccule_w8": "Utricle vs Saccule, w8",
        "B_utricle_vs_saccule_w22": "Utricle vs Saccule, w22",
        "C_w22_vs_w8_saccule": "w22 vs w8, Saccule",
        "D_w22_vs_w8_utricle": "w22 vs w8, Utricle",
        "E_organ_by_age_interaction": "Organ x age interaction"}


def r2_of(pc, cov):
    m = (r2["PC"] == pc) & (r2["covariate"] == cov)
    return float(r2.loc[m, "r2"].iloc[0]) if m.any() else float("nan")


def var_of(pc):
    return float(var.loc[var["PC"] == pc, "percent_variance"].iloc[0])


# ================================================================= 1 title
s = add_slide()
add_text(s, "Rebuild at 291 libraries", Inches(0.9), Inches(2.0), Inches(11.5), Inches(1.0),
         size=42, bold=True, color=ACCENT)
add_text(s, "The same five contrasts, 31 more libraries, animal identity from the submission record",
         Inches(0.9), Inches(3.05), Inches(11.5), Inches(0.8), size=20)
add_text(s, "Side by side with the original 260-library analysis  |  nothing in ../ was modified",
         Inches(0.9), Inches(3.95), Inches(11.5), Inches(0.6), size=15, color=MUTED)
add_text(s, "Generated from Pantry/sac_vs_utr/rebuild_291",
         Inches(0.9), Inches(5.8), Inches(11.5), Inches(0.5), size=12, color=MUTED)

# ================================================================= 2 what changed
s = add_slide()
add_title(s, "What changed, and why", "Two independent changes, both enabled by the complete submission workbook")
bullets(s, [
    "31 libraries added. Everything with RNA-to-genotype concordance >= 0.90 is now included. "
    "The original run additionally kept only ONE library per strain per organ-age group.",
    ("That rule exists because QTL mapping needs one genotype per sample. Differential expression "
     "has no such requirement - a second animal of the same strain is a valid biological replicate, "
     "and strain stays in the model as a blocking factor either way.", 1),
    ("The 30 libraries below the 0.90 concordance threshold remain excluded; their identity is "
     "not established.", 1),
    "Animal identity now comes from the submission record. Consecutive library IDs with the same "
    "recorded strain and the two organs identify one animal directly.",
    ("The original run inferred it from strain + age + sex, which cannot separate two mice of the "
     "same inbred strain. That reproduced 115 of 116 pairs - and the one disagreement was a genuine "
     "error, merging the saccule of one BXD128 animal with the utricle of another.", 1),
    "The 31 new libraries were quantified here with the same kallisto 0.48.0 index and the same "
    "parameters as the original Pantry run, so the count matrices are directly comparable.",
    "Everything else is unchanged: same gene filters, same models, same thresholds, same gene sets.",
], Inches(0.7), Inches(1.5), Inches(12.0), Inches(5.6), size=13.5)

# ================================================================= 3 design
s = add_slide()
add_title(s, "Sample counts", "291 libraries against the original 260")
g = samples.groupby(["organ", "age"]).size().unstack(fill_value=0)
rows = []
for organ in ["Saccule", "Utricle"]:
    for age in ["w8", "w22"]:
        rows.append([organ, age, int(g.loc[organ, age])])
t = pd.DataFrame(rows, columns=["Organ", "Age", "Libraries (291 run)"])
old = {"Saccule_w8": 73, "Saccule_w22": 55, "Utricle_w8": 74, "Utricle_w22": 58}
t["Libraries (260 run)"] = [old[f"{r.Organ}_{r.Age}"] for r in t.itertuples()]
t["Change"] = t["Libraries (291 run)"] - t["Libraries (260 run)"]
t = t[["Organ", "Age", "Libraries (260 run)", "Libraries (291 run)", "Change"]]
add_table(s, t, Inches(0.7), Inches(1.6), Inches(6.0), Inches(2.0), font_size=12)

paired = samples[samples["has_organ_pair"]]
pr = paired.groupby("age")["pair_id"].nunique()
pt = pd.DataFrame({"Age": ["w8", "w22", "Total"],
                   "Animal pairs (260 run)": [70, 52, 122],
                   "Animal pairs (291 run)": [int(pr.get("w8", 0)), int(pr.get("w22", 0)),
                                              int(paired["pair_id"].nunique())]})
pt["Change"] = pt["Animal pairs (291 run)"] - pt["Animal pairs (260 run)"]
add_table(s, pt, Inches(7.1), Inches(1.6), Inches(5.6), Inches(1.7), font_size=12)
src = samples["pair_source"].value_counts()
bullets(s, [
    f"Animal identity comes from the submission record for {int(src.get('submission record', 0))} "
    f"libraries and falls back to strain + age + sex for {int(src.get('strain + age + sex', 0))} "
    "(six w22 pairs in batch 3 whose two organs were not submitted adjacently).",
    "Libraries with no partner organ are still dropped from the organ contrasts - they carry no "
    "information once animal identity is in the model.",
], Inches(0.7), Inches(4.2), Inches(12.0), Inches(2.2), size=13)

# ================================================================= 4 PCA
fig_slide("Global structure is unchanged", "01_pca_pc1_pc2.png",
          f"PCA of all 291 libraries, top 2000 variable genes.",
          f"PC1 ({var_of('PC1')}%) still separates the organs, and organ still explains "
          f"R2 = {r2_of('PC1','organ'):.2f} of it. Adding 31 libraries did not disturb the structure.")
fig_slide("Covariate structure is unchanged", "02_pca_covariate_r2.png",
          "R-squared of each component against each covariate.",
          f"Organ on PC1 ({r2_of('PC1','organ'):.2f}), strain on PC2 ({r2_of('PC2','strain'):.2f}), "
          f"age on PC3 ({r2_of('PC3','age'):.2f}) with batch ({r2_of('PC3','batch'):.2f}) and sex "
          f"({r2_of('PC3','sex'):.2f}) alongside it - the same picture as the 260-library run.")

# ================================================================= 5 DE headline
s = add_slide()
add_title(s, "Differential expression: 260 vs 291", "FDR < 0.05, and FDR < 0.05 with |log2FC| > 1")
tab = cmp.copy()
tab = tab[["label", "n1", "n2", "fdr1", "fdr2", "sig1", "sig2", "d_sig"]]
tab.columns = ["Contrast", "Libs 260", "Libs 291", "FDR only (260)", "FDR only (291)",
               "DEGs (260)", "DEGs (291)", "Change"]
add_table(s, tab, Inches(0.45), Inches(1.5), Inches(12.45), Inches(2.3), font_size=10.5,
          col_widths=[3.0, 1.0, 1.0, 1.5, 1.5, 1.3, 1.3, 1.0])
bullets(s, [
    "The FDR-only columns are where extra libraries buy power - that number is a sensitivity "
    "measure. The DEG columns are gated by the two-fold cut and move much less.",
    "Fold changes themselves are essentially unchanged between runs (next slide), so this is added "
    "sensitivity, not a different answer.",
], Inches(0.5), Inches(4.1), Inches(12.3), Inches(2.4), size=13)

two_fig_slide("Where the extra libraries land",
              "21_v1_vs_v2_fdr_counts.png", "20_v1_vs_v2_counts.png",
              "Left: genes at FDR < 0.05. Right: genes also passing |log2FC| > 1.",
              "The age contrasts (C and D) gain the most, which is what was expected - they are "
              "library-limited, not pair-limited.")

# ================================================================= 6 agreement
s = add_slide()
add_title(s, "Do the two runs agree?", "Gene-level concordance between the 260- and 291-library results")
ct = conc[["label", "genes", "r", "sig_v1", "sig_v2", "shared", "pct_v1_kept"]].copy()
ct["r"] = ct["r"].round(3)
ct.columns = ["Contrast", "Genes compared", "r (log2FC)", "DEGs 260", "DEGs 291",
              "Shared", "% of 260 kept"]
add_table(s, ct, Inches(0.6), Inches(1.5), Inches(12.1), Inches(2.3), font_size=10.5,
          col_widths=[3.0, 1.5, 1.3, 1.2, 1.2, 1.1, 1.5])
add_picture_fit(s, FIGS / "22_v1_vs_v2_lfc.png", Inches(1.2), Inches(4.0), Inches(11.0), Inches(3.1))

# ================================================================= 7 GSEA
s = add_slide()
add_title(s, "Pathway results", "GSEA on the 291-library rankings")
gt = g2.copy()
CN = {"hallmark": "Hallmark", "go_bp": "GO:BP", "reactome": "Reactome"}
gt["collection"] = gt["collection"].map(CN)
gt = gt[["contrast_label", "collection", "sig", "top", "top_NES"]]
gt.columns = ["Contrast", "Collection", "FDR < 0.05", "Strongest set", "NES"]
gt["Strongest set"] = gt["Strongest set"].str.slice(0, 38)
add_table(s, gt, Inches(0.5), Inches(1.45), Inches(12.3), Inches(5.2), font_size=9.5,
          col_widths=[2.6, 1.6, 1.2, 3.6, 0.9])

two_fig_slide("Pathways: 260 vs 291",
              "30_gsea_v1_vs_v2.png", "31_gsea_nes_agreement.png",
              "Left: significant sets per contrast. Right: Hallmark NES, run against run.",
              f"Hallmark enrichment scores correlate at r = "
              f"{hm['nes1'].corr(hm['nes2']):.3f}. "
              f"{int((hm['call']=='both').sum())} sets are called by both runs, "
              f"{int((hm['call']=='291 only').sum())} only by the 291 run, "
              f"{int((hm['call']=='260 only').sum())} only by the 260 run.")

# ================================================================= 8 verdict
s = add_slide()
add_title(s, "What this means for the paper", "Read the numbers, then decide")
bullets(s, [
    "The biology does not change. Fold changes correlate near-perfectly between runs, the PCA "
    "structure is identical, and the pathway story is the same. Nothing published from the "
    "260-library run would need retracting.",
    "What the rebuild buys is sensitivity, concentrated in the age contrasts, plus one corrected "
    "animal pair.",
    "Reasons to adopt it: more libraries, animal identity taken from the record rather than "
    "inferred, and a pairing error fixed. Reasons to keep the 260: it is already written up, "
    "and it is the set that matches the QTL mapping exactly.",
    "A middle option: keep the 260 as the primary analysis and cite the 291 rebuild as a "
    "sensitivity analysis showing the results are stable to sample selection.",
    "Whichever is chosen, the BXD128 w22 pairing error should be corrected - it is a correctness "
    "issue, not a power one.",
], Inches(0.7), Inches(1.5), Inches(12.0), Inches(5.6), size=13.5)

# ================================================================= 9 files
s = add_slide()
add_title(s, "Where everything lives", "rebuild_291/  (the 260-library analysis in ../ is untouched)")
add_text(s,
         "scripts/\n"
         "  00_config.R                  settings for this rebuild only\n"
         "  run_kallisto.slurm           quantifies the 31 libraries Pantry never ran\n"
         "  01_build_inputs.R            tximport over all 291, sample table, animal identity\n"
         "  02_pca.R                     PCA and covariate R-squared\n"
         "  03_differential_expression.R the same five contrasts\n"
         "  04_figures_and_compare.R     volcanoes plus every 260-vs-291 comparison\n"
         "  05_gsea.R                    GSEA and its own 260-vs-291 comparison\n"
         "  06_build_report.py           this deck\n\n"
         "data/        sample table, count matrix, cached DESeq2 objects\n"
         "quant/       kallisto output for the 31 new libraries\n"
         "results/tables/   DE_summary.tsv, v1_vs_v2_summary.tsv, v1_vs_v2_concordance.tsv,\n"
         "                  gsea_summary.tsv, gsea_v1_vs_v2_hallmark.tsv, per-contrast results\n"
         "results/figures/  PCA, volcanoes, and the comparison figures\n",
         Inches(0.7), Inches(1.5), Inches(12.0), Inches(5.6), size=12.5, spacing=1.15)

OUT.parent.mkdir(parents=True, exist_ok=True)
prs.save(OUT)
print("Wrote", OUT)
print(len(prs.slides._sldIdLst), "slides")
