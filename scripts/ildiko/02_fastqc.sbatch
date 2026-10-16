#!/bin/bash
#SBATCH --job-name=fp_fastqc
#SBATCH --partition=courses
#SBATCH --array=1-33%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/courses/BINF6430.202710/students/polyak.i/FinalProject/logs/ildiko/fastqc_%A_%a.out

# Run FastQC on both reads of one sample per array task.
set -eo pipefail

# ---- Paths ----
PROJ=/courses/BINF6430.202710/students/polyak.i/FinalProject
RAW=/scratch/polyak.i/binf6430/FinalProject/raw
OUT=/scratch/polyak.i/binf6430/FinalProject/fastqc
ACC_LIST=$PROJ/metadata/SRR_Acc_List.txt

# ---- Pick this task's SRR ID ----
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ACC_LIST" | tr -d '\r')
echo "Task ${SLURM_ARRAY_TASK_ID}: ${SRR} on $(hostname) at $(date)"

# ---- Skip if both reports already exist ----
if [[ -f "$OUT/${SRR}_1_fastqc.html" && -f "$OUT/${SRR}_2_fastqc.html" ]]; then
    echo "${SRR} already has FastQC reports, skipping."
    exit 0
fi

# ---- Activate environment ----
eval "$(/home/polyak.i/.local/bin/micromamba shell hook --shell bash)"
micromamba activate fp-qc

# ---- Run FastQC on both reads ----
fastqc -t "$SLURM_CPUS_PER_TASK" -o "$OUT" \
    "$RAW/${SRR}_1.fastq.gz" "$RAW/${SRR}_2.fastq.gz"

echo "${SRR} FastQC finished at $(date)"
