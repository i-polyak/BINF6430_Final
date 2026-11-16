#!/bin/bash
#SBATCH --job-name=heatmap
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=00:15:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --partition=courses

cd $HOME/binf6430_repro_env

apptainer exec --pwd $HOME/binf6430_repro_env $HOME/bioc315.sif \
  Rscript 08_heatmap.R

echo "Job finished with exit code $?"
