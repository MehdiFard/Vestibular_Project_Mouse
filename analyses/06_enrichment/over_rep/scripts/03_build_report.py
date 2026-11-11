#!/usr/bin/env python3
"""03 - Assemble the over-representation results into a PowerPoint deck.

Reads the tables and figures written by scripts 01-02 and produces
report/over_representation_report.pptx.
"""

from pathlib import Path

import pandas as pd
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Emu, Inches, Pt

PA = Path("../../../analyses/06_enrichment/over_rep")
FIGS = PA / "figures"
TABLES = PA / "tables"
OUT = PA / "report" / "over_representation_report.pptx"

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
GSEA = PA.parent / "gsea"
summ = pd.read_csv(PA / "tables" / "ora_summary.tsv", sep="\t")
counts = pd.read_csv(PA / "tables" / "ora_input_gene_counts.tsv", sep="\t")
comb = pd.read_csv(PA / "tables" / "ora_all_contrasts.tsv.gz", sep="\t")
vs = pd.read_csv(PA / "tables" / "ora_vs_gsea_counts.tsv", sep="\t")
agree = pd.read_csv(PA / "tables" / "hallmark_ora_vs_gsea.tsv", sep="\t")

CONTRASTS = [
    ("A_organ_w8", "A - Utricle vs Saccule at 8 weeks"),
    ("B_organ_w22", "B - Utricle vs Saccule at 22 weeks"),
    ("C_age_sac", "C - 22 vs 8 weeks in Saccule"),
    ("D_age_utr", "D - 22 vs 8 weeks in Utricle"),
    ("E_interaction", "E - Organ x age interaction"),
]
COLL = {"hallmark": "Hallmark", "go_bp": "GO Biological Process", "reactome": "Reactome"}


def top_terms(cid, cn, n=10):
    d = comb[(comb["contrast"] == cid) & (comb["collection"] == cn)
             & (comb["gene_list"] != "all") & (comb["p.adjust"] < 0.05)]
    d = d.reindex(d["p.adjust"].sort_values().index).head(n)
    return pd.DataFrame({
        "Pathway": d["pathway_label"].str.slice(0, 40),
        "List": d["gene_list"],
        "Genes": d["k"].astype(int).astype(str) + "/" + d["n"].astype(int).astype(str),
        "Fold": d["fold_enrichment"].round(1),
        "FDR": d["p.adjust"].apply(lambda v: f"{v:.1e}"),
    }).reset_index(drop=True)


# ================================================================= 1 title
s = add_slide()
add_text(s, "Over-representation analysis", Inches(0.9), Inches(2.0), Inches(11.5),
         Inches(1.0), size=42, bold=True, color=ACCENT)
add_text(s, "The classical pathway test, run on the DE gene lists",
         Inches(0.9), Inches(3.05), Inches(11.5), Inches(0.6), size=22)
add_text(s, "Hypergeometric test  |  same gene sets and same contrasts as ../gsea  |  "
            "included for comparison",
         Inches(0.9), Inches(3.75), Inches(11.5), Inches(0.6), size=15, color=MUTED)
add_text(s, "Generated from Pantry/sac_vs_utr/pathway_analysis/over_rep",
         Inches(0.9), Inches(5.7), Inches(11.5), Inches(0.5), size=12, color=MUTED)

# ================================================================= 2 method
s = add_slide()
add_title(s, "ORA and GSEA are different questions",
          "Both are valid; they just use different inputs")
bullets(s, [
    "ORA takes only the DE genes - FDR < 0.05 and |log2FC| > 1 - and asks, by a "
    "hypergeometric test, whether a gene set is over-represented in that list. A gene is "
    "in or out; nothing below the cut counts.",
    "GSEA keeps all ~17,000 genes, ranks them by the Wald statistic, and asks whether a "
    "set clusters at one end. Nothing is discarded.",
    "Universe: the background here is the genes TESTED in that contrast, not the whole "
    "genome and not all of MSigDB. Using a larger universe than was actually tested makes "
    "every p-value look better than it is - the most common way ORA is got wrong.",
    "Direction: up and down genes are tested separately as well as together. A set that is "
    "half up and half down can look enriched in the combined list while being biologically "
    "incoherent.",
    "Everything else is held identical to the GSEA run - same MSigDB 2026.1.Mm sets, same "
    "15-500 size limits, same FDR < 0.05 - so the two are directly comparable.",
], Inches(0.7), Inches(1.5), Inches(12.0), Inches(5.6), size=14)

# ================================================================= 3 inputs
s = add_slide()
add_title(s, "What ORA had to work with",
          "The DE lists are all that ORA sees, out of ~17,000-18,000 genes tested")
add_picture_fit(s, FIGS / "summary/01_de_gene_counts.png",
                Inches(0.4), Inches(1.4), Inches(6.5), Inches(5.6))
ct = counts.copy()
ct = ct[["contrast_label", "universe", "n_all", "n_up", "n_down"]]
ct.columns = ["Contrast", "Universe", "DE genes", "Up", "Down"]
add_table(s, ct, Inches(7.2), Inches(2.0), Inches(5.7), Inches(2.2), font_size=11,
          col_widths=[2.6, 1.1, 1.1, 0.7, 0.8])
add_text(s, "Contrasts C, D and E give ORA only 122-168 genes to work with. That is the "
            "whole input. GSEA uses every tested gene in all five.",
         Inches(7.2), Inches(4.5), Inches(5.7), Inches(1.4), size=12.5, color=MUTED,
         italic=True, spacing=1.2)

# ================================================================= 4-6 comparison
fig_slide("GSEA finds more, except where the DE list is large",
          "summary/02_ora_vs_gsea_counts.png",
          "Significant gene sets at FDR < 0.05. ORA run on the combined DE list.",
          "GSEA wins in Hallmark and Reactome in all five contrasts. GO:BP is the exception: ORA "
          "returns 510 terms for contrast A against GSEA's 128, because a 756-gene list hits many "
          "small overlapping GO terms - more terms, not more biology. Where the DE list is small the "
          "gap goes the other way and is total: D gives ORA 146 genes and no Hallmark set at all.")

fig_slide("Why ORA misses what GSEA finds",
          "summary/03_worked_example_oxphos.png",
          "Oxidative phosphorylation in contrast D - the strongest single result in the GSEA run.",
          "132 of its 194 genes are significantly changed, 88% in the same direction, but not one "
          "reaches two-fold. The set is coherently shifted and entirely invisible to a fold-change cut.")

n_both = int((agree["call"] == "both").sum())
n_g = int((agree["call"] == "GSEA only").sum())
n_o = int((agree["call"] == "ORA only").sum())
fig_slide("Where the two methods agree",
          "summary/04_hallmark_agreement.png",
          "Every Hallmark set in every contrast: GSEA NES against ORA significance.",
          f"{n_both} sets are called by both, {n_g} by GSEA only, {n_o} by ORA only. ORA is close to a "
          "subset of GSEA here - it rarely finds anything GSEA missed, and misses a great deal that "
          "GSEA found.")

# ================================================================= per contrast
for cid, label in CONTRASTS:
    left = FIGS / f"hallmark/{cid}_ora_bar.png"
    right = FIGS / f"go_bp/{cid}_ora_bar.png"
    if not left.exists() and not right.exists():
        continue
    s = add_slide()
    add_title(s, f"{label}: ORA", "Left: Hallmark. Right: GO Biological Process. "
                                  "Bars point right for the up list, left for the down list.")
    if left.exists():
        add_picture_fit(s, left, Inches(0.45), Inches(1.5), Inches(6.1), Inches(5.1))
    else:
        add_text(s, "No Hallmark set reached FDR < 0.05\nin either direction.",
                 Inches(0.45), Inches(3.4), Inches(6.1), Inches(1.0), size=16,
                 color=MUTED, italic=True, align=PP_ALIGN.CENTER)
    if right.exists():
        add_picture_fit(s, right, Inches(6.75), Inches(1.5), Inches(6.1), Inches(5.1))
    tt = top_terms(cid, "hallmark", 8)
    if len(tt):
        s2 = add_slide()
        add_title(s2, f"{label}: top ORA terms", "Hallmark, from the direction-split lists")
        add_table(s2, tt, Inches(1.4), Inches(1.6), Inches(10.5), Inches(3.2),
                  font_size=11, col_widths=[4.0, 0.9, 1.1, 0.8, 1.2])
        add_text(s2, "Genes = how many of the DE list fall in that set, out of the DE genes that "
                     "are in any set. Fold = how much more often than expected by chance.",
                 Inches(1.4), Inches(5.3), Inches(10.5), Inches(0.8), size=12,
                 color=MUTED, italic=True, spacing=1.2)

# ================================================================= verdict
s = add_slide()
add_title(s, "Which to use", "They answer different questions, so the answer depends on the question")
bullets(s, [
    "Use GSEA when you want to know what is happening biologically. It uses all the data, "
    "detects coordinated small changes, and does not depend on where you put a threshold.",
    "Use ORA when you already have a gene list you care about and want to characterise it - "
    "a cluster, a module, a validated marker panel, a list from someone else's paper.",
    "In this dataset GSEA is clearly the better tool. The DE lists are small (122-863 genes) "
    "because |log2FC| > 1 is strict, and the real biology - oxidative phosphorylation, "
    "interferon response - lives in coordinated sub-two-fold changes that never enter the list.",
    "ORA does reproduce the strongest signals: epithelial-mesenchymal transition in the organ "
    "contrasts, interferon and inflammatory response in the saccule ageing contrast. Where the "
    "effect is large and gene-level, both methods agree.",
    "ORA is also sensitive to a choice GSEA does not have to make - the universe. Swapping the "
    "tested-gene background for 'all genes in MSigDB' would inflate every result here.",
], Inches(0.7), Inches(1.5), Inches(12.0), Inches(5.6), size=13.5)

# ================================================================= files
s = add_slide()
add_title(s, "Where everything lives", "pathway_analysis/over_rep/")
add_text(s,
         "scripts/\n"
         "  01_run_ora.R            hypergeometric test, 5 contrasts x 3 collections x 3 lists\n"
         "  02_figures.R            all 19 figures, including the GSEA comparison\n"
         "  03_build_report.py      this deck\n\n"
         "tables/\n"
         "  ora_summary.tsv               counts per contrast x collection x list\n"
         "  ora_all_contrasts.tsv.gz      every set tested, with gene counts and fold enrichment\n"
         "  ora_input_gene_counts.tsv     DE list sizes and universes\n"
         "  ora_vs_gsea_counts.tsv        the head-to-head counts\n"
         "  hallmark_ora_vs_gsea.tsv      per-set agreement between the two methods\n"
         "  hallmark/ go_bp/ reactome/    one table per contrast per direction\n\n"
         "figures/\n"
         "  summary/    DE list sizes, ORA vs GSEA counts, the worked example, agreement\n"
         "  hallmark/ go_bp/ reactome/    ORA bars per contrast\n\n"
         "../gsea/      the GSEA run on the same contrasts",
         Inches(0.7), Inches(1.5), Inches(12.0), Inches(5.6), size=12.5, spacing=1.15)

OUT.parent.mkdir(parents=True, exist_ok=True)
prs.save(OUT)
print("Wrote", OUT)
print(len(prs.slides._sldIdLst), "slides")
