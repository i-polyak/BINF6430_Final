#!/bin/bash
#SBATCH --job-name=deseq2_mouse
#SBATCH --partition=courses
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=deseq2_mouse_%j.out

conda deactivate 2>/dev/null || true

SCRIPT_PATH=/scratch/caldwell.ja/final_project/src/BINF6430_Final/scripts/jacque

source $SCRIPT_PATH/run_R421.sh

rscript_r421 $SCRIPT_PATH/deseq2_mouse.R matched sv1
rscript_r421 $SCRIPT_PATH/deseq2_mouse.R matched none
rscript_r421 $SCRIPT_PATH/deseq2_mouse.R matched batch
rscript_r421 $SCRIPT_PATH/deseq2_mouse.R allGEO sv1
rscript_r421 $SCRIPT_PATH/deseq2_mouse.R allGEO none
rscript_r421 $SCRIPT_PATH/deseq2_mouse.R allGEO batch
rscript_r421 $SCRIPT_PATH/deseq2_mouse.R authorsSV1 sv1
