#!/bin/bash
#SBATCH --job-name=hdlss_l2_b1
#SBATCH --output=logs/hdlss_%A_%a.out
#SBATCH --error=logs/hdlss_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0-02:00:00
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --array=1-900
#SBATCH --mail-user=rag334@pitt.edu
#SBATCH --mail-type=END,FAIL

# ---------------------------------------------------------------------------
# Batch 1: p in {200, 500, 1000}
# Grid size: 900 configs (60 identity + 420 cs + 420 ar1)
#
# For batch 2 (p = 2000):
#   - Edit run_config.R: set p_grid <- c(2000)
#   - Edit this file:    --array=1-300, --job-name=hdlss_l2_b2
#                        consider --time=0-04:00:00 and --mem=16G
#
# For batch 3 (p = 10000):
#   - Edit run_config.R: set p_grid <- c(10000)
#   - Edit this file:    --array=1-300, --job-name=hdlss_l2_b3
#                        --time=1-00:00:00 and --mem=32G recommended
# ---------------------------------------------------------------------------

# Load R module
module purge
module load gcc r

# Prevent HDNRA's OpenMP from oversubscribing CPUs
# (each SLURM task only gets 1 CPU, so threading must match)
export OMP_NUM_THREADS=1

# Run the R script -- it reads SLURM_ARRAY_TASK_ID automatically
Rscript run_config.R
