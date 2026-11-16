#!/bin/bash
#SBATCH --job-name=fp_download
#SBATCH --partition=courses
#SBATCH --array=1-33%6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/courses/BINF6430.202710/students/polyak.i/FinalProject/logs/ildiko/download_%A_%a.out

# Download one SRA run per array task, convert to FASTQ, compress, verify.
set -eo pipefail

# ---- Paths ----
PROJ=/courses/BINF6430.202710/students/polyak.i/FinalProject
RAW=/scratch/polyak.i/binf6430/FinalProject/raw
ACC_LIST=$PROJ/metadata/SRR_Acc_List.txt

# ---- Pick this task's SRR ID (line N of the list; strip any Windows line endings) ----
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ACC_LIST" | tr -d '\r')
TMP=$RAW/tmp/$SRR

echo "Task ${SLURM_ARRAY_TASK_ID}: ${SRR} on $(hostname) at $(date)"

# ---- Skip if this sample already finished ----
if [[ -f "$RAW/${SRR}.done" ]]; then
    echo "${SRR} already done, skipping."
    exit 0
fi

# ---- Activate environment (aliases and shell functions don't exist in batch jobs) ----
eval "$(/home/polyak.i/.local/bin/micromamba shell hook --shell bash)"
micromamba activate fp-download

# ---- 1. Download .sra ----
prefetch "$SRR" -O "$RAW" --max-size 50G

# ---- 2. Convert to paired FASTQ (spot/read counts print to this log) ----
mkdir -p "$TMP"
fasterq-dump "$RAW/$SRR" -O "$RAW" -t "$TMP" -e "$SLURM_CPUS_PER_TASK"

# ---- 3. Compress (replaces .fastq with .fastq.gz) ----
pigz -p "$SLURM_CPUS_PER_TASK" "$RAW/${SRR}_1.fastq" "$RAW/${SRR}_2.fastq"

# ---- 4. Verify compressed files are intact ----
pigz -t "$RAW/${SRR}_1.fastq.gz" "$RAW/${SRR}_2.fastq.gz"

# ---- 5. Clean up only after everything above succeeded ----
rm -rf "$RAW/$SRR" "$TMP"
touch "$RAW/${SRR}.done"

echo "${SRR} finished at $(date)"
EOF
