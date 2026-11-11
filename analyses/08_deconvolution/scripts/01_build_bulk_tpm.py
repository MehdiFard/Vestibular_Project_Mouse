#!/usr/bin/env python3
"""
01_build_bulk_tpm.py  --  Build linear-scale bulk TPM matrices for deconvolution.

Deconvolution (Bisque / MuSiC / ssGSEA) needs LINEAR expression (TPM), not the
inverse-normal-transformed values used for QTL mapping. Pantry's gene-level TPM
lives in <group>/unnorm/expression.bed (genes = ENSMUSG, GRCm38).

We pool the two ages within each tissue (ages kept together, tissues separate):
    saccule = saccule_w8 + saccule_w22
    utricle = utricle_w8 + utricle_w22
and relabel genes ENSMUSG -> symbol (to match the single-cell reference, which
uses gene symbols). Output: one gene(symbol) x sample TPM CSV per tissue.

Run: python3 scripts/01_build_bulk_tpm.py
(Light I/O; safe to run without a compute allocation.)
"""
import os
import pandas as pd

PANTRY = "../../analyses/03_pantry"
GTF    = os.path.join(PANTRY, "phenotyping/input/ref/Mus_musculus.GRCm38.102.chr.gtf")
OUTDIR = "results/inputs"

# tissue -> list of Pantry groups to pool (ages together)
TISSUES = {
    "saccule": ["saccule_w8", "saccule_w22"],
    "utricle": ["utricle_w8", "utricle_w22"],
}

def unnorm_bed(group):
    return os.path.join(PANTRY, "phenotyping/output", group, "unnorm", "expression.bed")

def build_ensembl_to_symbol(gtf_path):
    """Stream the GTF 'gene' lines -> {ENSMUSG: symbol}."""
    ens2sym = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attr = f[8]
            gid = sym = None
            for field in attr.strip().split(";"):
                field = field.strip()
                if field.startswith("gene_id "):
                    gid = field.split('"')[1]
                elif field.startswith("gene_name "):
                    sym = field.split('"')[1]
            if gid and sym:
                ens2sym[gid] = sym
    return ens2sym

def load_group_tpm(group):
    """Read one group's unnorm gene TPM bed -> DataFrame (index=ENSMUSG, cols=samples)."""
    df = pd.read_csv(unnorm_bed(group), sep="\t")
    df = df.set_index("phenotype_id").drop(columns=["#chr", "start", "end"])
    return df

def main():
    os.makedirs(OUTDIR, exist_ok=True)
    print("[1] Parsing GTF for ENSMUSG -> symbol ...")
    ens2sym = build_ensembl_to_symbol(GTF)
    print(f"    mapped {len(ens2sym):,} genes")

    for tissue, groups in TISSUES.items():
        print(f"[2] {tissue}: pooling {groups}")
        mats = [load_group_tpm(g) for g in groups]
        # same gene set across ages; align on gene index, concatenate sample columns
        tpm = pd.concat(mats, axis=1, join="inner")
        print(f"    genes={tpm.shape[0]:,}  samples={tpm.shape[1]}")

        # ENSMUSG -> symbol; drop unmapped; collapse duplicate symbols by summing TPM
        tpm.index = tpm.index.map(lambda g: ens2sym.get(g))
        tpm = tpm[tpm.index.notna()]
        n_before = tpm.shape[0]
        tpm = tpm.groupby(level=0).sum()
        print(f"    symbol-mapped genes={tpm.shape[0]:,} (collapsed from {n_before:,})")

        out = os.path.join(OUTDIR, f"bulk_tpm_{tissue}.csv")
        tpm.to_csv(out)
        print(f"    saved {out}")

    print("DONE")

if __name__ == "__main__":
    main()
