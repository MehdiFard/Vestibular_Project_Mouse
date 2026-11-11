#!/bin/bash
### Quantify the 31 libraries that Pantry never ran, using the same kallisto
### index and the same parameters as the original phenotyping run.
#SBATCH --job-name=v2_kallisto
#SBATCH --output=logs/%x-%A_%a.log
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --partition=<cpu_partition>
#SBATCH --account=<account>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=4:00:00
#SBATCH --array=1-31%8

set -euo pipefail
RB="$PWD/analyses/05_dge/rebuild_291"
PH="$PWD/analyses/03_pantry/phenotyping"
KAL="$(conda info --base)/envs/pantry/bin/kallisto"
cd "$RB"

read -r LIB DATASET R1 R2 < <(sed -n "${SLURM_ARRAY_TASK_ID}p" data/quant_jobs.tsv)
IDX="$PH/intermediate/${DATASET}/reference/kallisto.idx"
OUT="$RB/quant/${LIB}"

echo "[$SLURM_ARRAY_TASK_ID] $LIB  dataset=$DATASET  $(date)"
mkdir -p "$OUT"
"$KAL" quant --index "$IDX" --output-dir "$OUT" --threads "$SLURM_CPUS_PER_TASK" "$R1" "$R2"
echo "done $LIB $(date)"
