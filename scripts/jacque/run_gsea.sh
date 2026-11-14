#!/bin/bash
SCRIPT_PATH=/scratch/caldwell.ja/final_project/scripts
source $SCRIPT_PATH/run_R421.sh

for set in matched allGEO; do
  for design in sv1 none batch; do
    for r in stat lfc; do
      rscript_r421 $SCRIPT_PATH/fgsea_mouse.R $set $design $r
    done
  done
done
