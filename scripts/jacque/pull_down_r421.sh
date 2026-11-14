#!/bin/bash

export APPTAINERS_CACHEDIR=/scratch/caldwell.ja/apptainer_cache
mkdir -p /scratch/caldwell.ja/final_project/sif
cd /scratch/caldwell.ja/final_project/sif

apptainer pull r421.sif docker://rocker/r-ver:4.2.1
apptainer exec r421.sif R --version

apptainer pull bioc315.sif docker:://bioconductor/bioconductor_docker:RELEASE_3_15
apptainer exec bioc315.sif R --version

