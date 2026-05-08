#!/bin/bash
#SBATCH --job-name=hdlss_sim
#SBATCH --output=logs/hdlss_%A_%a.out
#SBATCH --error=logs/hdlss_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0-12:00:00
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --array=1-192
#SBATCH --mail-user=rag334@pitt.edu
#SBATCH --mail-type=END,FAIL

# Load R modules
module purge
module load gcc r

# Run the R script — it reads SLURM_ARRAY_TASK_ID automatically
Rscript run_config.R
