#!/bin/bash
#SBATCH --job-name=tximport
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=short

# BINF 6430 reproducibility project
# Runs 03_import_and_coldata.R (tx2gene + tximport) inside the
# Bioconductor 3.15 / renv environment

cd $HOME/binf6430_repro_env

apptainer exec --pwd $HOME/binf6430_repro_env $HOME/bioc315.sif \
  Rscript 03_import_and_coldata.R

echo "Job finished with exit code $?"
