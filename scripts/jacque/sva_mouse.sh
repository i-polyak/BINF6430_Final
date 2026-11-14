#!/bin/bash
#SBATCH --job-name=sva_mouse
#SBATCH --partition=courses
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=sva_mouse_%j.out

conda deactivate 2>/dev/null || true

source /scratch/caldwell.ja/final_project/scripts/run_R421.sh

rscript_r421 /scratch/caldwell.ja/final_project/scripts/sva_mouse.R matched aswritten
rscript_r421 /scratch/caldwell.ja/final_project/scripts/sva_mouse.R matched standard
rscript_r421 /scratch/caldwell.ja/final_project/scripts/sva_mouse.R allGEO aswritten
rscript_r421 /scratch/caldwell.ja/final_project/scripts/sva_mouse.R allGEO standard
