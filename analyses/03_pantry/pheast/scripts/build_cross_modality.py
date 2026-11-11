#!/usr/bin/env python3
# Cross-modality phenotypes with the expression-regressed RNA editing rows.
# Copied from the human project (vestibulart_project_human_pantry, 05_pantry), with the dataset (group) as argument.
# cross_modality is a plain concatenation of the 7 modality tables with "modality:" prefixes
# (phenotyping/scripts/combine_modalities_updated_060425.sh, no re-normalization), so swapping the RNA editing
# rows of the fixed cross_modality input for the regressed values is the same as rebuilding it.
# Check first: the fixed cross_modality RNA editing rows equal the fixed RNA_editing input.
# Run from analyses/03_pantry/pheast (rule pheast_inputs in ../Snakefile):  python3 scripts/build_cross_modality.py <group>
import os, shutil, subprocess, sys
import numpy as np
import pandas as pd

DS = sys.argv[1]
FIX = f"input/{DS}/phenotypes_unregressed"
REG = f"input/{DS}/phenotypes"
read = lambda p: pd.read_csv(p, sep="\t", dtype={"#chr": str}, low_memory=False)

cross = read(f"{FIX}/cross_modality.bed.gz")
ed_fix = read(f"{FIX}/RNA_editing.bed.gz").set_index("phenotype_id")
ed_reg = read(f"{REG}/RNA_editing.bed.gz").set_index("phenotype_id")
subj = list(cross.columns[4:])
is_ed = cross.phenotype_id.str.startswith("RNA_editing:")
ids = cross.loc[is_ed, "phenotype_id"].str.replace("RNA_editing:", "", n=1)

same = np.allclose(cross.loc[is_ed, subj].to_numpy(float), ed_fix.loc[ids, subj].to_numpy(float), atol=1e-5)
print(f"cross_modality rows: {len(cross):,} | RNA editing rows: {is_ed.sum():,} | identical to fixed RNA_editing input: {same}")
assert same and is_ed.sum() == len(ed_reg)

cross.loc[is_ed, subj] = ed_reg.loc[ids, subj].to_numpy()
out = f"{REG}/cross_modality.bed"
cross.to_csv(out, sep="\t", index=False, float_format="%g")
subprocess.run(["bgzip", "-f", out], check=True)
subprocess.run(["tabix", "-f", "-p", "bed", out + ".gz"], check=True)
shutil.copy2(f"{FIX}/cross_modality.phenotype_groups.txt", f"{REG}/cross_modality.phenotype_groups.txt")
print(f"wrote {out}.gz (+ .tbi, phenotype_groups)")
