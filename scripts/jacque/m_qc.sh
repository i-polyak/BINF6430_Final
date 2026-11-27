#!/bin/bash
# Use `source <scriptname>` to run this script so you can use the MULTIQC alias!
#    ENVIRONMENT

# Variable to the path of image/container
MULTIQC_PATH=/shared/container_repository/explorer/multiqc/1.26/multiqc_latest.sif

# temp alias for the program using signularity
alias MULTIQC="singularity exec -B /scratch:/scratch,/courses:/courses,/home:/home $MULTIQC_PATH multiqc"
