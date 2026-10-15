#!/bin/bash
#SBATCH --job-name=brms_sor_an
#SBATCH --account=rc_petrullo_pi
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=200g
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gabriesc@unc.edu
#SBATCH --output=gtdbtk_classify_%j.out
#SBATCH --error=gtdbtk_classify_%j.err


# Go to working directory
cd /proj/petrullolab/schiro/brms_sky_islands
echo "directory set"


# Load R module
module load anaconda
conda activate /proj/petrullolab/software/conda_envs/rbayes

echo "R module loaded at $(date)"
which Rscript
Rscript --version


# Run the R script
Rscript code_HPC_an_diet.R 