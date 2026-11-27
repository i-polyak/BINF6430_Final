#!/bin/bash
#SBATCH --job-name=fqdump
#SBATCH --partition=short
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=fqdump_%j.log

eval "$(conda shell.bash hook)"
conda activate rnaseq_download

cd /scratch/caldwell.ja/final_project/data/GSE288245
mkdir -p fastq tmp

for srr in $(cat SRR_Acc_List.txt); do
  if [[ -s fastq/${srr}_1.fastq.gz && -s fastq/${srr}_2.fastq.gz && ! -e fastq/${srr}_1.fastq && ! -e fastq/${srr}_2.fastq ]]; then
    echo "skip $srr (done)"; continue
  fi
  echo "processing $srr  $(date)"
  rm -f fastq/${srr}_*.fastq fastq/${srr}_*.fastq.gz
  fasterq-dump --split-files --threads 8 sra/"$srr" -O fastq/ --temp tmp/
  pigz -p 8 fastq/${srr}_*.fastq
done
echo "all done $(date)"
