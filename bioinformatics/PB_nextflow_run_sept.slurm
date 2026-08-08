#!/bin/bash
#SBATCH --job-name=hifi16s
#SBATCH --time=2-00:00:00        # adjust as needed (1 day here)
#SBATCH --cpus-per-task=32        # only for the Nextflow head job
#SBATCH --mem=128G
#SBATCH -p general               # your partition
#SBATCH -q public                # your QOS
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user="%u@asu.edu"
#SBATCH --export=NONE            # start with a clean environment

# --- software env ---
module load nextflow-25.04.3-gcc-12.1.0
module load mamba
module load openjdk-17.0.3_7-gcc-12.1.0            # or anaconda/mamba on your system

mkdir -p /scratch/gschiro/tmp
export TMPDIR=/scratch/gschiro/tmp
export R_MAX_VSIZE='128Gb'


cp -r /home/gschiro/PB_16s_workflow/HiFi-16S-workflow /scratch/gschiro/HiFi-16S-workflow_4_new
cd /scratch/gschiro/HiFi-16S-workflow_4_new


# --- launch pipeline ---
# Assumes you already made samples.tsv with absolute paths to your 5 FASTQs
nextflow run main.nf  \
  --input /data/tsuzuki7/04_wild_rodent/data_2025_pb_log_run/file_list_2025.tsv \
  --metadata /data/tsuzuki7/04_wild_rodent/data_2025_pb_log_run/condition_2025.tsv \
  --outdir /data/tsuzuki7/04_wild_rodent/data_2025_pb_log_run/results_final_4_err \
  --process.queue general \
  --publish_dir_mode copy \
  --dada2_cpu 32  \
  --vsearch_cpu 32 \
  --max_ee 4
