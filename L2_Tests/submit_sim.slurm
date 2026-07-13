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

# Load R module
module purge
module load gcc r
export OMP_NUM_THREADS=1
Rscript run_config.R
