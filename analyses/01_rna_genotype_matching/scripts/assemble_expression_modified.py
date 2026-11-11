import argparse
from gtfparse import read_gtf
import numpy as np
import pandas as pd
from pathlib import Path
import rnaseqnorm

parser = argparse.ArgumentParser(description="Combine RSEM output files into BED files with log2 counts, TPM, and IQN units")
parser.add_argument("rsem_dir", type=Path, help="Directory with RSEM output files to be combined. Each sample is named from its file prefix.")
parser.add_argument("samples", type=Path, help="File containing sample IDs, one per line, for which to load RSEM outputs")
parser.add_argument("annotations", help="GTF annotation file.")
parser.add_argument("out_prefix", help="Prefix (including path) of output files, which will be appended with '.log2.bed', '.tpm.bed', and '.iqn.bed'")
args = parser.parse_args()

with args.samples.open() as f:
    ids = f.read().splitlines()

log2_list = []  # List to store log2 DataFrames
tpm_list = []  # List to store tpm DataFrames
sample_names = []  # List to store sample names

for i, sample in enumerate(ids):
    fname = args.rsem_dir / f"{sample}.genes.results.gz"
    d = pd.read_csv(fname, sep="\t", index_col="gene_id")
    if i == 0:
        common_index = d.index
    else:
        assert d.index.equals(common_index)
    # Store sample names
    sample_names.append(sample)
    log2_list.append(np.log2(d["expected_count"] + 1))
    tpm_list.append(d["TPM"])

log2 = pd.concat(log2_list, axis=1)
tpm = pd.concat(tpm_list, axis=1)

# Assign correct sample names to columns
log2.columns = sample_names
tpm.columns = sample_names

iqn = rnaseqnorm.normalize_quantiles(tpm)
iqn = rnaseqnorm.inverse_normal_transform(iqn)

anno = read_gtf(args.annotations, result_type='pandas')
anno = anno.loc[anno["feature"] == "gene", :]
anno["tss"] = np.where(anno["strand"] == "+", anno["start"], anno["end"])
anno["start"] = anno["tss"] - 1  # BED coordinates are 0-based
anno["end"] = anno["tss"]
anno["#chr"] = anno["seqname"]
anno = anno.sort_values(["#chr", "start"])
anno = anno[["#chr", "start", "end", "gene_id"]]

# Save filtered IQN file to avoid tensorQTL error on phenotypes with 1 nonzero value:
iqnfilt = iqn.copy()
iqnfilt = iqnfilt.loc[np.sum(iqnfilt > 0, axis=1) > 1, :]

log2 = anno.merge(log2.reset_index(), on="gene_id", how="inner")
tpm = anno.merge(tpm.reset_index(), on="gene_id", how="inner")
iqn = anno.merge(iqn.reset_index(), on="gene_id", how="inner")
iqnfilt = anno.merge(iqnfilt.reset_index(), on="gene_id", how="inner")

# Filtering and handling chromosome names
iqnfilt = iqnfilt.loc[iqnfilt["#chr"].isin([str(i) for i in range(1, 20)]), :]

def chr_to_int(chr):
    try:
        return int(chr)
    except ValueError:
        return chr

iqnfilt["#chr"] = iqnfilt["#chr"].apply(chr_to_int)  # Handle non-numeric values
iqnfilt = iqnfilt.sort_values(["#chr", "start"])

log2.to_csv(f"{args.out_prefix}.log2.bed", sep="\t", index=False, float_format="%g")
tpm.to_csv(f"{args.out_prefix}.tpm.bed", sep="\t", index=False, float_format="%g")
iqn.to_csv(f"{args.out_prefix}.iqn.bed", sep="\t", index=False, float_format="%g")
iqnfilt.to_csv(f"{args.out_prefix}.iqn.filtered.bed", sep="\t", index=False, float_format="%g")

