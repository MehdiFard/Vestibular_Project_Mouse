#!/bin/bash
#SBATCH --job-name=sac_utr_gsea
#SBATCH --output=../logs/%x-%j.log
#SBATCH --error=../logs/%x-%j.err
#SBATCH --partition=<cpu_partition>
#SBATCH --account=<account>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
set -euo pipefail
ROOT="$PWD/analyses/06_enrichment"
E="$(conda info --base)/envs/enrich_r"
cd "$ROOT/gsea"
export TMPDIR="$ROOT/logs"
echo "GSEA start $(date)"
"$E/bin/Rscript" scripts/01_run_gsea.R
echo "GSEA end $(date)"
