#!/bin/bash
#SBATCH --job-name=tximport
#SBATCH --partition=courses
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=tximport_%j.out

conda deactivate 2>/dev/null || true

source /scratch/caldwell.ja/final_project/scripts/run_R421.sh

mkdir -p /scratch/caldwell.ja/final_project/results

rscript_r421 /scratch/caldwell.ja/final_project/scripts/tximport_mouse.R

