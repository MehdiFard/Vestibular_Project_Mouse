#!/bin/bash
#SBATCH -p <cpu_partition>
#SBATCH -A <account>
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=128g
#SBATCH --time=8:00:00
#SBATCH -J vestib_tissue
#SBATCH -o logs/tissue_analysis_%j.log
#SBATCH -e logs/tissue_analysis_%j.log
# Per-tissue integration + annotation (saccule, utricle) then utricle-vs-saccule DE.
# Submit with:  sbatch scripts/run_tissue_analysis.sbatch

set -o pipefail
cd analyses/07_single_cell_rna || exit 1
R="$(conda info --base)/envs/seurat_env/bin/Rscript"

echo "=== START $(date) on $(hostname) ==="

echo ">>> [1/3] integrate + annotate SACCULE"
$R scripts/annot5_integrate_tissue.R saccule || { echo "SACCULE FAILED"; exit 1; }

echo ">>> [2/3] integrate + annotate UTRICLE"
$R scripts/annot5_integrate_tissue.R utricle || { echo "UTRICLE FAILED"; exit 1; }

echo ">>> [3/3] utricle-vs-saccule differential expression"
$R scripts/annot6_tissue_DE.R || { echo "DE FAILED"; exit 1; }

echo "=== ALL DONE $(date) ==="
