#!/bin/bash
#SBATCH --job-name=fp_kal_index
#SBATCH --partition=courses
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/courses/BINF6430.202710/students/polyak.i/FinalProject/logs/ildiko/kallisto_index_%j.out

# Build kallisto 0.48.0 index from unmodified GENCODE vM30 transcripts (default k=31).
set -eo pipefail

REF=/scratch/polyak.i/binf6430/FinalProject/reference
FASTA=$REF/gencode.vM30.transcripts.fa.gz
IDX=$REF/gencode.vM30.kallisto_0.48.0.idx

eval "$(/home/polyak.i/.local/bin/micromamba shell hook --shell bash)"
micromamba activate fp-quant

echo "Started on $(hostname) at $(date)"
kallisto version
md5sum "$FASTA"

kallisto index -i "$IDX" "$FASTA"

ls -lh "$IDX"
echo "Finished at $(date)"
