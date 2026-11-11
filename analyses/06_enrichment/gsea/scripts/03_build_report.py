#!/usr/bin/env python3
"""03 - Assemble the GSEA results into a PowerPoint deck.

Reads the tables and figures written by scripts 01-02 and produces
report/pathway_analysis_report.pptx.
"""

from pathlib import Path

import pandas as pd
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Emu, Inches, Pt

PA = Path("../../../analyses/06_enrichment/gsea")
FIGS = PA / "figures"
TABLES = PA / "tables"
OUT = PA / "report" / "pathway_analysis_report.pptx"

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
summ = pd.read_csv(TABLES / "gsea_summary.tsv", sep="\t")
comb = pd.read_csv(TABLES / "gsea_all_contrasts.tsv.gz", sep="\t")

CONTRASTS = [
    ("A_organ_w8", "A - Utricle vs Saccule at 8 weeks"),
    ("B_organ_w22", "B - Utricle vs Saccule at 22 weeks"),
    ("C_age_sac", "C - 22 vs 8 weeks in Saccule"),
    ("D_age_utr", "D - 22 vs 8 weeks in Utricle"),
    ("E_interaction", "E - Organ x age interaction"),
]
COLL = {"hallmark": "Hallmark", "go_bp": "GO Biological Process", "reactome": "Reactome"}


def top_paths(cid, cn, n=10):
    d = comb[(comb["contrast"] == cid) & (comb["collection"] == cn) & (comb["padj"] < 0.05)]
    d = d.reindex(d["padj"].sort_values().index).head(n)
    return pd.DataFrame({
        "Pathway": d["pathway_label"].str.slice(0, 44),
        "NES": d["NES"].round(2),
        "FDR": d["padj"].apply(lambda v: f"{v:.1e}"),
        "Size": d["size"],
        "Direction": d["direction"].str.slice(0, 32),
    }).reset_index(drop=True)


# ================================================================= 1 title
s = add_slide()
add_text(s, "Pathway analysis", Inches(0.9), Inches(2.0), Inches(11.5), Inches(1.0),
         size=44, bold=True, color=ACCENT)
add_text(s, "GSEA across all five saccule / utricle contrasts",
         Inches(0.9), Inches(3.05), Inches(11.5), Inches(0.6), size=22)
add_text(s, "fgsea on DESeq2 Wald statistics  |  mouse-native MSigDB 2026.1.Mm  |  "
            "Hallmark, GO Biological Process, Reactome",
         Inches(0.9), Inches(3.75), Inches(11.5), Inches(0.6), size=15, color=MUTED)
add_text(s, "Generated from Pantry/sac_vs_utr/pathway_analysis/gsea",
         Inches(0.9), Inches(5.7), Inches(11.5), Inches(0.5), size=12, color=MUTED)

# ================================================================= 2 method
s = add_slide()
add_title(s, "How the analysis was run", "Choices that matter for interpreting the rest")
bullets(s, [
    "Ranking metric: the DESeq2 Wald statistic, one value per gene per contrast. It "
    "carries both the direction of the effect and its precision.",
    ("Ranking on raw log2FC instead would put noisy low-count genes at the extremes - "
     "the same problem that puts Cma1 (mean count 11) at the top of the interaction DEG list.", 1),
    "No pre-filtering. The whole ranked list goes into GSEA, not just the genes passing "
    "FDR < 0.05 and |log2FC| > 1.",
    ("This is the point of GSEA: it detects a set of genes moving together by a small "
     "amount, which a per-gene threshold discards. Thresholding first would waste the method.", 1),
    "Gene sets: mouse-native MSigDB 2026.1.Mm, keyed by Ensembl gene ID, so no human-to-mouse "
    "ortholog mapping and no ID conversion anywhere in the pipeline.",
    "Sets restricted to 15-500 genes. fgsea multilevel, eps = 0, seed 42. FDR is "
    "Benjamini-Hochberg within each contrast x collection.",
    "Positive NES always means the set sits at the top of that contrast's ranking - the "
    "direction is spelled out on every figure, since it differs between contrasts.",
], Inches(0.7), Inches(1.5), Inches(12.0), Inches(5.6), size=14)

# ================================================================= 3 summary table
s = add_slide()
add_title(s, "What each contrast turned up",
          "Gene sets significant at FDR < 0.05, of those tested")
tab = summ.copy()
tab["collection"] = tab["collection"].map(COLL)
tab = tab[["contrast_label", "collection", "pathways_tested", "sig", "sig_up", "sig_down", "top"]]
tab.columns = ["Contrast", "Collection", "Tested", "Significant", "NES +", "NES -", "Strongest set"]
tab["Strongest set"] = tab["Strongest set"].str.slice(0, 34)
add_table(s, tab, Inches(0.4), Inches(1.45), Inches(12.5), Inches(5.0), font_size=10,
          col_widths=[2.6, 2.2, 1.0, 1.2, 0.8, 0.8, 3.0])

# ================================================================= 4-6 summary figs
fig_slide("Significantly enriched sets per contrast",
          "summary/01_significant_per_contrast.png",
          "FDR < 0.05. The y scale differs by collection.",
          "The interaction (E) returns the most sets in every collection, and almost all of them "
          "negative - the two organs diverge with age far more than either differs at a fixed age.")

fig_slide("Hallmark enrichment across all five contrasts",
          "summary/02_hallmark_nes_heatmap.png",
          "NES per set per contrast; stars mark FDR. 42 of 50 Hallmark sets are significant somewhere.",
          "Read it in columns: A and B are the organ axis, C and D the ageing axis, E the difference "
          "between them. Epithelial-mesenchymal transition is the one set significant in all five.")

fig_slide("Do the contrasts agree at the pathway level?",
          "summary/03_contrast_similarity.png",
          "Spearman correlation of GO:BP normalised enrichment scores between contrasts.",
          "Three blocks: the two organ contrasts agree (0.73), the two ageing contrasts agree (0.71), "
          "and organ vs age are independent (0.04 to -0.16). E anti-correlates with C (-0.57) because "
          "the saccule is what drives the interaction.")

# ================================================================= per contrast
for cid, label in CONTRASTS:
    d = comb[(comb["contrast"] == cid) & (comb["padj"] < 0.05)]
    up_dir = d[d["NES"] > 0]["direction"].iloc[0] if (d["NES"] > 0).any() else "-"
    dn_dir = d[d["NES"] < 0]["direction"].iloc[0] if (d["NES"] < 0).any() else "-"
    dirline = f"Positive NES = {up_dir}   |   negative NES = {dn_dir}"

    two_fig_slide(f"{label}: Hallmark",
                  f"hallmark/{cid}_nes_bar.png", f"hallmark/{cid}_dotplot.png",
                  dirline,
                  "Left: NES for the strongest sets. Right: the same sets with gene-set size and "
                  "significance.")

    s = add_slide()
    add_title(s, f"{label}: strongest Hallmark sets", dirline)
    add_picture_fit(s, FIGS / f"hallmark/{cid}_enrichment_curves.png",
                    Inches(0.4), Inches(1.45), Inches(6.4), Inches(5.5))
    add_table(s, top_paths(cid, "hallmark", 10), Inches(7.0), Inches(1.6),
              Inches(5.9), Inches(3.6), font_size=9.5,
              col_widths=[3.0, 0.8, 1.0, 0.7, 2.0])
    add_text(s, "The running-score curves show where each set's genes fall in the ranked list. "
                "A peak near the left edge means the set is concentrated at the positive end.",
             Inches(7.0), Inches(5.45), Inches(5.9), Inches(1.2), size=11.5, color=MUTED,
             italic=True, spacing=1.2)

    two_fig_slide(f"{label}: GO Biological Process and Reactome",
                  f"go_bp/{cid}_nes_bar.png", f"reactome/{cid}_nes_bar.png",
                  dirline,
                  "Left: GO Biological Process. Right: Reactome. Both show the top 10 each way by "
                  "adjusted p.")

# ================================================================= biology
s = add_slide()
add_title(s, "What the pathways say", "Reading the five contrasts together")
bullets(s, [
    "The two organs differ metabolically at 8 weeks. Utricle runs higher oxidative "
    "phosphorylation, mTORC1 signalling and cholesterol homeostasis; saccule runs higher "
    "epithelial-mesenchymal transition and apical junction.",
    "By 22 weeks that difference has become inflammatory. TNFA/NFkB, complement, interferon "
    "gamma and allograft rejection are all significantly higher in saccule at w22 but not at w8.",
    "Ageing is inflammatory in both organs, but much more so in saccule. Interferon gamma is "
    "the single strongest result in the whole analysis (NES +3.12, FDR 6.5e-30 in saccule) "
    "against NES +2.42 in utricle.",
    "The utricle's ageing signature is different in kind: it loses oxidative phosphorylation "
    "(NES -3.28) and MYC targets, which the saccule does not.",
    "That asymmetry is exactly what the interaction picks up, and why 36 of its 37 significant "
    "Hallmark sets are negative: the saccule gains inflammation faster, and the utricle loses "
    "metabolism faster, and both push the organ gap the same way.",
    "Epithelial-mesenchymal transition is higher in saccule at both ages and is the only "
    "Hallmark set significant in all five contrasts.",
], Inches(0.7), Inches(1.5), Inches(12.0), Inches(5.6), size=13.5)

# ================================================================= caveats
s = add_slide()
add_title(s, "Caveats", "What not to over-read")
bullets(s, [
    "GO Biological Process is highly redundant - many of its 7,781 sets share most of their "
    "genes, so 1,203 significant sets for the interaction is not 1,203 independent findings. "
    "Hallmark is the non-redundant view and the one to quote.",
    "Contrast E is fitted without strain, unlike A and B. The interaction estimates are "
    "therefore noisier than the organ contrasts, and its gene-level list is the least "
    "reproducible result in the project (37% agreement with limma).",
    "C and D compare mostly different animals across ages, and batch is partly confounded "
    "with age. The inflammatory ageing signal is strong enough to survive that, but the "
    "effect sizes should not be taken as precise.",
    "GSEA tests whether a set is non-randomly placed in the ranking. It does not establish "
    "that a pathway is active, nor direction of causation.",
    "Enrichment of immune and interferon sets in bulk tissue usually reflects a change in "
    "cell composition (infiltrating immune cells) as much as regulation within a cell type. "
    "The mast-cell protease genes in the interaction DEG list point the same way.",
], Inches(0.7), Inches(1.5), Inches(12.0), Inches(5.6), size=13.5)

# ================================================================= files
s = add_slide()
add_title(s, "Where everything lives", "pathway_analysis/")
add_text(s,
         "scripts/\n"
         "  01_run_gsea.R           fgsea over 5 contrasts x 3 collections\n"
         "  02_figures.R            all 28 figures\n"
         "  03_build_report.py      this deck\n"
         "  run_gsea.slurm          step 1 as a SLURM job\n\n"
         "tables/\n"
         "  gsea_summary.tsv               counts per contrast x collection\n"
         "  gsea_all_contrasts.tsv.gz      every set, every contrast, with leading-edge genes\n"
         "  contrast_nes_correlation.tsv   the similarity matrix\n"
         "  hallmark/ go_bp/ reactome/     one table per contrast\n\n"
         "figures/\n"
         "  summary/    the three cross-contrast figures\n"
         "  hallmark/   NES bars, dot plots, enrichment curves per contrast\n"
         "  go_bp/ reactome/   NES bars per contrast\n\n"
         "../over_rep/  the same contrasts run as over-representation analysis",
         Inches(0.7), Inches(1.5), Inches(12.0), Inches(5.6), size=12.5, spacing=1.15)

OUT.parent.mkdir(parents=True, exist_ok=True)
prs.save(OUT)
print("Wrote", OUT)
print(len(prs.slides.__iter__.__self__._sldIdLst), "slides")
