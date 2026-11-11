#!/bin/bash
### Running Pheast by pantry

#SBATCH --output=logs/%x-%j.log
#SBATCH --error=logs/%x-%j.err
#SBATCH --partition=cpu-core-sponsored
#SBATCH --account=cpu-core-sponsored
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=34
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<email>


## Go to WD
# cd analyses/03_pantry/pheast/


## Activate environment
MAMBA="$(conda info --base)/etc/profile.d/mamba.sh"
source $MAMBA
mamba activate pantry_pheast

# ========== Check DATASET ==========
: "${DATASET:?Environment variable DATASET not set. Please submit with --export=ALL,DATASET=your_dataset_name}"

snakemake   --configfile config.yml \
            --config dataset=$DATASET \
            --jobs 34 \
            --default-resources mem_mb=256000 cpus=34 \
            --printshellcmds \
            --latency-wait 60


# -R tensorqtl_trans
# --use-conda \
# --rerun-incomplete \
# --keep-going \

### END