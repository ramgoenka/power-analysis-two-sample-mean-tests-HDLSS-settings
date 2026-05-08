#!/bin/bash
#SBATCH --job-name=hdlss_hotelling
#SBATCH --output=logs_hotelling/hotelling_%A_%a.out
#SBATCH --error=logs_hotelling/hotelling_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --array=1-144
#SBATCH --mail-user=rag334@pitt.edu
#SBATCH --mail-type=END,FAIL

# Load R modules
module purge
module load gcc r

# Create output directories if needed
mkdir -p results_hotelling
mkdir -p logs_hotelling

# Run the R script
Rscript run_hotelling.R
