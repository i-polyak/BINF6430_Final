#!/bin/bash
#SBATCH --job-name=fp_kal_quant
#SBATCH --partition=courses
#SBATCH --array=1-33%6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/courses/BINF6430.202710/students/polyak.i/FinalProject/logs/ildiko/quant_%A_%a.out

# kallisto 0.48.0 quant, unstranded baseline, one sample per array task.
set -eo pipefail

# ---- Paths ----
PROJ=/courses/BINF6430.202710/students/polyak.i/FinalProject
RAW=/scratch/polyak.i/binf6430/FinalProject/raw
IDX=/scratch/polyak.i/binf6430/FinalProject/reference/gencode.vM30.kallisto_0.48.0.idx
OUTDIR=/scratch/polyak.i/binf6430/FinalProject/kallisto/baseline_unstranded
ACC_LIST=$PROJ/metadata/SRR_Acc_List.txt

# ---- Pick this task's SRR ID ----
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ACC_LIST" | tr -d '\r')
echo "Task ${SLURM_ARRAY_TASK_ID}: ${SRR} on $(hostname) at $(date)"

# ---- Skip if already done ----
if [[ -f "$OUTDIR/${SRR}.done" ]]; then
    echo "${SRR} already quantified, skipping."
    exit 0
fi

# ---- Activate environment ----
eval "$(/home/polyak.i/.local/bin/micromamba shell hook --shell bash)"
micromamba activate fp-quant
kallisto version

# ---- kallisto only creates the last folder level, so make the parent first ----
mkdir -p "$OUTDIR"

# ---- Quantify (no strand flag = unstranded; no -b = no bootstraps) ----
kallisto quant -i "$IDX" -o "$OUTDIR/$SRR" -t "$SLURM_CPUS_PER_TASK" \
    "$RAW/${SRR}_1.fastq.gz" "$RAW/${SRR}_2.fastq.gz"

touch "$OUTDIR/${SRR}.done"
echo "${SRR} finished at $(date)"
