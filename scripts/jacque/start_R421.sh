#!/bin/bash

export R_LIBS_USER=/scratch/caldwell.ja/final_project/Rlib_4.2.1
export R_421=/scratch/caldwell.ja/final_project/sif/r421.sif

mkdir -p $R_LIBS_USER

run_r421(){

    apptainer exec --bind /scratch/caldwell.ja $R_421 R "$@"

}
