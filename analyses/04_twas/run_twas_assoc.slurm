#!/bin/bash
### Running TWAS association test

#SBATCH --job-name=twas_assoc
#SBATCH --output=logs/%x-%j.log
#SBATCH --error=logs/%x-%j.err
#SBATCH --partition=cpu-core-sponsored
#SBATCH --account=cpu-core-sponsored
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<email>


## Go to WD
# cd analyses/04_twas


## Activate environment
MAMBA="$(conda info --base)/etc/profile.d/mamba.sh"
source $MAMBA
mamba activate pantry_pheast


# Run Snakemake
snakemake -s run_twas_with_tissue.smk \
            --default-resources mem_mb=64000 cpus=24 \
            --jobs 24 \
            --printshellcmds \
            --latency-wait 60

### END