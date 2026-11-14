#!/bin/bash
#SBATCH --job-name=fgsea_mouse
#SBATCH --partition=courses
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/caldwell.ja/final_project/results/logs/fgsea_mouse_%j.out

conda deactivate 2>/dev/null || true

SCRIPT_PATH=/scratch/caldwell.ja/final_project/scripts

source $SCRIPT_PATH/run_R421.sh

rscript_r421 $SCRIPT_PATH/fgsea_mouse.R paperS2

for set in matched allGEO; do
  for design in sv1 none batch; do
    rscript_r421 /full/path/to/fgsea_mouse.R $set $design
  done
done
