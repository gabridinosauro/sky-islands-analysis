#!/usr/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=dada2_16s
#SBATCH --output=outputr%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --account=laurenpetrullo
#SBATCH --mail-type=ALL
#SBATCH --mail-user=schiro@arizona.edu
#SBATCH --partition=standard  
#SBATCH --time=100:00:00
#SBATCH --ntasks=94
#SBATCH --nodes=1
#SBATCH --mem=470gb

# --------------------------------------------------
# Load modules here
# --------------------------------------------------
module load contrib
module load micromamba
eval "$(micromamba shell hook --shell bash)"  # Initialize micromamba for bash
micromamba activate /contrib/laurenpetrullo/dada2

# --------------------------------------------------
# Execute commands here
# --------------------------------------------------
cd /xdisk/laurenpetrullo/schiro/vsearch_GTDB
Rscript /xdisk/laurenpetrullo/schiro/vsearch_GTDB/tax_ass.R 