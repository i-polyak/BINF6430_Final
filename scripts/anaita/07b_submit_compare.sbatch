#!/bin/bash
#SBATCH --job-name=compare_supplement
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=2
#SBATCH --partition=courses

cd $HOME/binf6430_repro_env

apptainer exec --pwd $HOME/binf6430_repro_env $HOME/bioc315.sif \
  Rscript 07b_compare_supplement.R

echo "Job finished with exit code $?"
