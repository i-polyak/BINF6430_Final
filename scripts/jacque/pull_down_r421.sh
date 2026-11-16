#!/bin/bash

export APPTAINERS_CACHEDIR=/scratch/caldwell.ja/apptainer_cache
mkdir -p /scratch/caldwell.ja/final_project/sif
cd /scratch/caldwell.ja/final_project/sif

# Needed R version 4.2.1; but the r421.sif file didn't work for this pipeline,
# replaced by bioc315.sif

#  apptainer pull r421.sif docker://rocker/r-ver:4.2.1
#  apptainer exec r421.sif R --version

apptainer pull bioc315.sif docker:://bioconductor/bioconductor_docker:RELEASE_3_15
apptainer exec bioc315.sif R --version

