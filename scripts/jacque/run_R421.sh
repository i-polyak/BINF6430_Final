#!/bin/bash

export FINAL_PROJECT=/scratch/caldwell.ja/final_project
export CONTAINERS=$FINAL_PROJECT/sif
export R_LIBS_USER=/scratch/caldwell.ja/final_project/Rlib_bioc315
export R_LIBS=$R_LIBS_USER:/usr/local/lib/R/site-library/
export R_421=$CONTAINERS/bioc315.sif

export R_LIBS=/scratch/caldwell.ja/final_project/Rlib_bioc315:/usr/local/lib/R/site-library

mkdir -p $R_LIBS_USER

run_r421(){
	apptainer exec --bind /scratch/caldwell.ja $R_421 R "$@"
}

rscript_r421(){
	apptainer exec --bind /scratch/caldwell.ja $R_421 Rscript "$@"
}

