#!/bin/bash
# 07a_pool_genotypes.sh -- merge w8+w22 genotypes per tissue (shared SNPs) + PCA.
# Run in analyses/09_interaction_qtl_mapping (Snakefile; plink and plink2 on PATH). Outputs to results/pooled_inputs/.
set -euo pipefail
PLINK=plink
PLINK2=plink2
PANTRY=../../analyses/03_pantry
OUT=results/pooled_inputs
mkdir -p "$OUT"

for tis in saccule utricle; do
    echo "=== $tis ==="
    w8=$PANTRY/pheast/input/${tis}_w8/${tis}_w8
    w22=$PANTRY/pheast/input/${tis}_w22/${tis}_w22
    comm -12 <(cut -f2 ${w8}.bim | sort) <(cut -f2 ${w22}.bim | sort) > $OUT/${tis}.shared_snps.txt
    echo "shared SNPs: $(wc -l < $OUT/${tis}.shared_snps.txt)"
    $PLINK --bfile $w8  --extract $OUT/${tis}.shared_snps.txt --make-bed --out $OUT/${tis}_w8s  >/dev/null 2>&1
    $PLINK --bfile $w22 --extract $OUT/${tis}.shared_snps.txt --make-bed --out $OUT/${tis}_w22s >/dev/null 2>&1
    # merge; on allele-mismatch, drop the offending SNPs and retry
    if ! $PLINK --bfile $OUT/${tis}_w8s --bmerge $OUT/${tis}_w22s --make-bed --out $OUT/${tis}_pooled >/dev/null 2>&1; then
        echo "  merge conflict -> excluding $(wc -l < $OUT/${tis}_pooled-merge.missnp) SNPs and retrying"
        $PLINK --bfile $OUT/${tis}_w8s  --exclude $OUT/${tis}_pooled-merge.missnp --make-bed --out $OUT/${tis}_w8s2  >/dev/null 2>&1
        $PLINK --bfile $OUT/${tis}_w22s --exclude $OUT/${tis}_pooled-merge.missnp --make-bed --out $OUT/${tis}_w22s2 >/dev/null 2>&1
        $PLINK --bfile $OUT/${tis}_w8s2 --bmerge $OUT/${tis}_w22s2 --make-bed --out $OUT/${tis}_pooled >/dev/null 2>&1
    fi
    $PLINK2 --bfile $OUT/${tis}_pooled --pca 5 --out $OUT/${tis}_pca >/dev/null 2>&1
    echo "  pooled: $(wc -l < $OUT/${tis}_pooled.fam) samples, $(wc -l < $OUT/${tis}_pooled.bim) SNPs, PCA done"
    rm -f $OUT/${tis}_w8s* $OUT/${tis}_w22s*   # tidy intermediates
done
echo "DONE"
