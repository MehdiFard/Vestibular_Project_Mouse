#!/usr/bin/env python3
"""17 - Validate the data-derived sample annotation against the COMPLETE
submission spreadsheet found in the backup archive.

The sheet used when this project was built
(metadata_mouse_all_plates_with_matched_geno_ids.csv) covers FTL_588-FTL_945
only, so every batch-4 library was missing from it. The backup archive holds
the full four-batch workbook (FTL_588-FTL_1057), which lets the strain, organ,
age and sex assignments used here be checked against an independent record.

Writes results/tables/backup_sheet_validation.tsv (per library) and
backup_sheet_validation_summary.tsv.
"""

import re
from pathlib import Path

import pandas as pd

PROJ = Path("../..")
BASE = PROJ / "analyses" / "05_dge"
# submission workbook (byte-identical to the /data/rss backup copy, which compute nodes cannot read)
XLSX = PROJ / "data" / "fastq_bulk" / "bulk_rnaseq_metadata_mouse_b1_to_b4.xlsx"
T = BASE / "results" / "tables"


def norm(s):
    """Strain names to a comparable form: drop the sub-strain suffix after '/',
    upper-case, strip punctuation. BXD125/RwwJ -> BXD125."""
    if pd.isna(s):
        return None
    return re.sub(r"[^A-Z0-9]", "", str(s).split("/")[0].upper().strip())


# ------------------------------------------------------------------ the sheet
frames = []
xl = pd.ExcelFile(XLSX)
for s in xl.sheet_names:
    d = xl.parse(s).rename(columns=lambda c: str(c).strip())
    d["batch_sheet"] = s
    frames.append(d[["Library ID", "mouse strain", "tissue", "Sex", "Age", "batch_sheet"]])
sh = pd.concat(frames, ignore_index=True)

# batch-4 IDs carry a suffix: FTL_1022_RNAseq_oldmouse
sh["library_id"] = sh["Library ID"].astype(str).str.extract(r"(FTL_\d+)", expand=False)
sh = sh.dropna(subset=["library_id"])
sh["sheet_strain"] = sh["mouse strain"]
sh["sheet_strain_norm"] = sh["mouse strain"].map(norm)
sh["sheet_organ"] = sh["tissue"].astype(str).str.strip().str.capitalize()
sh["sheet_sex"] = sh["Sex"].astype(str).str.strip().str[0].str.upper()
sh["sheet_age"] = sh["Age"].astype(str).str.extract(r"(\d+)")[0].map(
    lambda v: f"w{v}" if pd.notna(v) else None)

print(f"backup workbook: {len(sh)} libraries, "
      f"FTL_{sh.library_id.str.replace('FTL_','').astype(int).min()}"
      f"-FTL_{sh.library_id.str.replace('FTL_','').astype(int).max()}")

# ------------------------------------------------------------------ compare
# master_metadata already carries PC1-PC5, so no second join is needed
md = pd.read_csv(T / "master_metadata.tsv", sep="\t")
m = md.merge(sh.drop(columns=["Library ID", "mouse strain", "tissue", "Sex", "Age"]),
             on="library_id", how="left")
m["geno_strain_norm"] = m["strain"].map(norm)
m["strain_agrees"] = m.geno_strain_norm == m.sheet_strain_norm
m["organ_agrees"] = m.organ == m.sheet_organ
m["age_agrees"] = m.age == m.sheet_age
m["sex_agrees"] = m.sex_final == m.sheet_sex


def classify(r):
    """Why does the genotype call differ from the sheet?"""
    if r.strain_agrees or pd.isna(r.sheet_strain_norm):
        return ""
    g, s = str(r.geno_strain_norm), str(r.sheet_strain_norm)
    # BXD131a / BXD147a / BXD183a: sub-line letter suffix on the sheet
    if s.rstrip("ABC") == g:
        return "spelling: sub-line suffix"
    # BXF213 for BXD213, C57BI6J for C57BL6J: transcription typos
    if s.replace("BXF", "BXD") == g or s.replace("I", "L") == g:
        return "spelling: typo"
    if r.strain_gap < 0.01:
        return "near-tie, ambiguous"
    return "genuine relabel"


m["disagreement_type"] = m.apply(classify, axis=1)

cols = ["library_id", "organ", "age", "batch", "strain", "sheet_strain",
        "strain_agrees", "disagreement_type", "strain_simil", "strain_runnerup",
        "strain_runnerup_simil", "strain_gap", "sheet_organ", "organ_agrees",
        "PC1", "sheet_age", "age_agrees", "sex_final", "sex_source", "sheet_sex",
        "sex_agrees", "batch_sheet"]
m[cols].to_csv(T / "backup_sheet_validation.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ summary
in_sheet = m.sheet_strain_norm.notna()
sexed = m.sheet_sex.isin(["F", "M"])
inferred = m.sex_source.str.contains("absent", na=False)
rows = [
    ("Libraries analysed", len(m)),
    ("Found in the complete backup sheet", int(in_sheet.sum())),
    ("Strain agrees", int(m.strain_agrees.sum())),
    ("Strain disagrees", int((~m.strain_agrees & in_sheet).sum())),
    ("  of which spelling only", int((m.disagreement_type.str.startswith("spelling")).sum())),
    ("  of which near-tie (gap < 0.01)", int((m.disagreement_type == "near-tie, ambiguous").sum())),
    ("  of which genuine relabel", int((m.disagreement_type == "genuine relabel").sum())),
    ("Organ agrees", int(m.organ_agrees.sum())),
    ("Age agrees", int(m.age_agrees.sum())),
    ("Sex recorded in the backup sheet", int(sexed.sum())),
    ("Sex agrees where recorded", int((m.sex_agrees & sexed).sum())),
    ("Sex inferred from expression", int(inferred.sum())),
    ("  of those, now checkable against the sheet", int((inferred & sexed).sum())),
    ("  of those, inference confirmed", int((inferred & sexed & m.sex_agrees).sum())),
    ("Animal pairs", int(m[m.has_organ_pair].pair_id.nunique())),
    ("  with a different strain call for the two organs",
     int((m[m.has_organ_pair].groupby("pair_id")["strain"].nunique() > 1).sum())),
]
summ = pd.DataFrame(rows, columns=["check", "n"])
summ.to_csv(T / "backup_sheet_validation_summary.tsv", sep="\t", index=False)
print()
print(summ.to_string(index=False))

print("\n=== strain disagreements ===")
d = m[~m.strain_agrees & in_sheet]
print(d[["library_id", "batch", "organ", "strain", "sheet_strain", "strain_simil",
         "strain_runnerup", "strain_gap", "disagreement_type"]]
      .sort_values(["disagreement_type", "library_id"]).to_string(index=False))

print("\n=== organ disagreements (PC1 adjudicates: saccule positive, utricle negative) ===")
print(m[~m.organ_agrees][["library_id", "strain", "organ", "sheet_organ", "PC1"]]
      .to_string(index=False))
