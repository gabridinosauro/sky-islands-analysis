#!/bin/bash -l

# --------------------------------------------------------------
### PART 1: Requests resources to run your job.
# --------------------------------------------------------------
### Optional. Set the job name
#SBATCH --job-name=sina
### Optional. Set the output filename.
### SLURM reads %x as the job name and %j as the job ID
#SBATCH --output=%x-%j.out
### REQUIRED. Specify the PI group for this job
#SBATCH --account=laubitz
### Optional. Request email when job begins and ends
### SBATCH --mail-type=ALL
### Optional. Specify email address to use for notification
### SBATCH --mail-user=schiro@arizona.edu
### REQUIRED. Set the partition for your job.
#SBATCH --partition=standard
### REQUIRED. Set the number of cores that will be used for this job. 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
### REQUIRED. Set the number of nodes
#SBATCH --nodes=1
### REQUIRED. Set the memory required for this job.
#SBATCH --mem=240gb
### REQUIRED. Specify the time required for this job, hhh:mm:ss
#SBATCH --time=30:00:00


# --------------------------------------------------------------
### PART 2: Executes bash commands to run your job
# --------------------------------------------------------------
module load micromamba
module load fasttree
cd /xdisk/laurenpetrullo/schiro/SILVA

# First run allignment
micromamba activate /contrib/laurenpetrullo/sina
sina -i filtered_ASVs_and_arch.fasta  -o allignment.fasta --db SILVA_138.2_SSURef_NR99_03_07_24_opt.arb -p 48


# Now  trim gaps
 micromamba activate /contrib/laurenpetrullo/trimal
 trimal -in allignment.fasta -out allignment_trimmed.fasta  -gt 0.2


# make a tree
module load fasttree
export OMP_NUM_THREADS=48

FastTreeMP -nt -gtr -log fasttree_no_fastest.log < allignment_trimmed.fasta > asv_tree_no_fastest.nwk