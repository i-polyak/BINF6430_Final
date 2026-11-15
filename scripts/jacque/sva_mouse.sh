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

SCRIPT_PATH=/scratch/caldwell.ja/final_project/src/BINF6430_Final/scripts/jacque

source $SCRIPT_PATH/run_R421.sh

rscript_r421 $SCRIPT_PATH/sva_mouse.R matched
rscript_r421 $SCRIPT_PATH/sva_mouse.R allGEO
