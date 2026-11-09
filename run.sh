#!/bin/bash
#SBATCH --job-name=med_surv
#SBATCH --time=1:00:00
#SBATCH --account=XXXX
#SBATCH --mail-type=end,fail
#SBATCH --mem=1g
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1000

module load R
R CMD BATCH --no-save --no-restore simu0.R