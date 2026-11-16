#!/bin/bash
#SBATCH --job-name=deseq2
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=00:45:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=short

# BINF 6430 reproducibility project
# Runs 05_deseq2.R inside the Bioconductor 3.15 / renv environment

cd $HOME/binf6430_repro_env

apptainer exec --pwd $HOME/binf6430_repro_env $HOME/bioc315.sif \
  Rscript 05_deseq2.R

echo "Job finished with exit code $?"
