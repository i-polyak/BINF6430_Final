#!/bin/bash

# script file to move the files off scratch so we save the important ones!
#
# FASTQ files can be downloaded from SRA so not needed
# FASTQC files can be calculated from FASTQ (can be easily regenerated)
#    multiqc*.html summarizes so keep this instead of the fastqc files.
# ENA/Gencode files -- can be downloaded any time.
# will need to keep the abdunance.tsv files from kallisto
#
# will need the refs/gencode.vM30.tx2gene.tsv
# will need kallisto {index.batch,quant.sbatch,samples.txt,qc_mapping.tsv,meta_b
loom.tsv}
# will also need the multiqc file which is in kallisto/multiqc_report/multiqc_re
port.html

# Directories with data
SCRATCH=/scratch/caldwell.ja/final_project
KALLISTO=${SCRATCH}/kallisto
MULTIQC=${KALLISTO}/multiqc_report

# back up to: (note ~/courses is a link to /courses)
# back up to: (note ~/courses is a link to /courses)
#
BACKUP=/home/caldwell.ja/courses/6430/final_proj
SAVE_KALLISTO=${BACKUP}/kallisto

mkdir -p ${SAVE_KALLISTO}
cd ${SAVE_KALLISTO}
cp ${SCRATCH}/refs/gencode.vM30.tx2gene.tsv .
cp ${KALLISTO}/multiqc_report/multiqc_report.html kallisto_multiqc_report.html 
cp ${KALLISTO}/index.sbatch .
cp ${KALLISTO}/quant.sbatch .
cp ${KALLISTO}/samples.txt .
cp ${KALLISTO}/qc_mapping.tsv .
cp ${KALLISTO}/meta_bloom.tsv .
cp ${KALLISTO}/*.yml .
cp ${KALLISTO}/after_kallisto_save.sh .
