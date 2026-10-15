#!/bin/bash
#SBATCH --job-name=brms_sor
#SBATCH --account=laurenpetrullo
#SBATCH --mail-type=ALL
#SBATCH --mail-user=schiro@arizona.edu
#SBATCH --partition=standard  
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120gb
#SBATCH --time=32:00:00
#SBATCH --output=sor_%A_%a.out
#SBATCH --error=ser_%A_%a.err

# Go to working directory
cd /xdisk/laurenpetrullo/schiro/sky_island/brms_models
echo "directory set"


# Load R module
module load R
echo "R module loaded at $(date)"
which Rscript
Rscript --version



# Run the R script
Rscript code_HPC.R 