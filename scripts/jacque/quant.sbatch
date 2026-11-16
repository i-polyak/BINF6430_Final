#!/bin/bash
#SBATCH --job-name=kal_quant
#SBATCH --partition=courses
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=1-33%8
#SBATCH --output=/scratch/%u/final_project/kallisto/logs/%x_%A_%a.out

set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate kallisto_48

BASE=/scratch/$USER/final_project
cd $BASE/kallisto
IDX=$BASE/refs/gencode.vM30.k31.idx
S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

kallisto version
kallisto quant -i $IDX -o quant/$S -b 100 -t 4 \
  fastq/${S}_1.fastq.gz fastq/${S}_2.fastq.gz
