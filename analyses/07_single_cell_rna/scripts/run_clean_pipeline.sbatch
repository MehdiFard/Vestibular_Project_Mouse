#!/bin/bash
#SBATCH -p <cpu_partition>
#SBATCH -A <account>
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=128g
#SBATCH --time=8:00:00
#SBATCH -J vestib_clean
#SBATCH -o logs/clean_pipeline_%j.log
#SBATCH -e logs/clean_pipeline_%j.log
# Ambient (SoupX) + doublet (scDblFinder) cleaning -> per-organ + all-sample
# integration/annotation -> comparison with previous runs.
# Submit with:  sbatch scripts/run_clean_pipeline.sbatch

set -o pipefail
cd analyses/07_single_cell_rna || exit 1
R="$(conda info --base)/envs/seurat_env/bin/Rscript"

echo "=== START $(date) on $(hostname) ==="
echo ">>> [1/3] per-sample cleaning (SoupX + scDblFinder)"
$R scripts/annot9_clean_qc.R       || { echo "CLEAN FAILED"; exit 1; }
echo ">>> [2/3] integration + annotation (cleaned)"
$R scripts/annot10_integrate_clean.R || { echo "INTEGRATE FAILED"; exit 1; }
echo ">>> [3/3] compare with previous runs"
$R scripts/annot11_compare_clean.R || { echo "COMPARE FAILED"; exit 1; }
echo "=== ALL DONE $(date) ==="
