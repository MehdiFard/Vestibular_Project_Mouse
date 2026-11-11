#!/usr/bin/env python

### 04 GWAS of one trait with FaST-LMM (https://fastlmm.github.io/FaST-LMM/), leave-one-chromosome-out kinship
# from the tested SNPs; noPCs = no covariates, withPCs = the strain genotype PCs in results/genotype_pca/<trait>.covar.pcs
# Usage (from analyses/10_gwas): python scripts/04_fastlmm.py <trait> <noPCs|withPCs>

import os
import sys
from fastlmm.association import single_snp

trait, model = sys.argv[1], sys.argv[2]
assert model in ("noPCs", "withPCs")
genofile = f"intermediate/genotype/{trait}_clean_RefAltOk"
phenofile = f"intermediate/phenotype/{trait}.pheno.txt"
covarfile = f"results/genotype_pca/{trait}.covar.pcs" if model == "withPCs" else None
os.makedirs(f"results/{model}", exist_ok=True)

print("##### Performing GWAS on:", trait, "data -", model, "#####")
single_snp(test_snps=genofile, pheno=phenofile, covar=covarfile, count_A1=False, leave_out_one_chrom=True,
           output_file_name=f"results/{model}/{trait}.gwas.loco.tsv")

print("End of Analyses")

### END
