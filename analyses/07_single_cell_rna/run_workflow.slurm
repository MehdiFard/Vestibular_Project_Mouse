#!/bin/bash
# Run workflow on SLURM

#SBATCH --job-name=scRNA_workflow_cellranger
#SBATCH --output=logs/%x-%j.log
#SBATCH --error=logs/%x-%j.err
#SBATCH --partition=cpu-core-sponsored
#SBATCH --account=cpu-core-sponsored
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=34
#SBATCH --mem=256G
#SBATCH --time=36:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<email>


## Go to WD
cd analyses/07_single_cell_rna

## Activate environment
MAMBA="$(conda info --base)/etc/profile.d/mamba.sh"
source $MAMBA
mamba activate Env1


# Run Snakemake
snakemake   --cores 34 \
            --default-resources mem_mb=256000 cpus=34 \
            --printshellcmds \
            --latency-wait 60 \
            --keep-going \
            --rerun-incomplete \
            --scheduler greedy

### END
