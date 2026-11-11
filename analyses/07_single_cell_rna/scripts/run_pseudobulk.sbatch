#!/bin/bash
#SBATCH -p <cpu_partition>
#SBATCH -A <account>
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=64g
#SBATCH --time=2:00:00
#SBATCH -J vestib_pb
#SBATCH -o logs/pseudobulk_%j.log
#SBATCH -e logs/pseudobulk_%j.log
# Pseudobulk PCA + edgeR DE (organ + cell level), protein_coding + lncRNA only.
# annot7 needs Seurat (seurat_env); annot8 needs edgeR/DESeq2 (Env1).
# Submit with:  sbatch scripts/run_pseudobulk.sbatch

set -o pipefail
cd analyses/07_single_cell_rna || exit 1
RSEURAT="$(conda info --base)/envs/seurat_env/bin/Rscript"
RENV1="$(conda info --base)/envs/Env1/bin/Rscript"

echo "=== START $(date) on $(hostname) ==="
echo ">>> [1/2] pseudobulk prep (seurat_env)"
$RSEURAT scripts/annot7_pseudobulk_prep.R || { echo "PREP FAILED"; exit 1; }
echo ">>> [2/2] PCA + edgeR DE (Env1)"
$RENV1 scripts/annot8_pseudobulk_DE.R || { echo "DE FAILED"; exit 1; }
echo "=== ALL DONE $(date) ==="
