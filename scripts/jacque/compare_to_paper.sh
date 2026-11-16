#!/bin/bash
#SBATCH --job-name=compare_paper
#SBATCH --partition=courses
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/caldwell.ja/final_project/results/logs/compare_paper_%j.out

conda deactivate 2>/dev/null || true

SCRIPT_PATH=/scratch/caldwell.ja/final_project/src/BINF6430_Final/scripts/jacque
source $SCRIPT_PATH/run_R421.sh

for set in matched allGEO; do
  for design in sv1 none batch; do
    rscript_r421 $SCRIPT_PATH/compare_to_paper.R $set $design
  done
done

rscript_r421 $SCRIPT_PATH/compare_to_paper.R authorsSV1 sv1

