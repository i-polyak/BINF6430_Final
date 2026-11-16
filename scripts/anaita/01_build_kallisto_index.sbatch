#!/bin/bash
# =============================================================================
# build_kallisto_index.sbatch
# Build the kallisto index for GENCODE vM30 (mm39) used by kallisto_quant.sbatch.
# Submit from the LOGIN node:  sbatch scripts/anaita/build_kallisto_index.sbatch
#
# =============================================================================

#SBATCH --job-name=kallisto_index
#SBATCH --partition=courses
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G                 # index build is memory-heavier than quant per-sample
#SBATCH --output=/scratch/lashley.an/final_project/logs/index_%j.out
#SBATCH --error=/scratch/lashley.an/final_project/logs/index_%j.err

set -euo pipefail

eval "$(conda shell.bash hook)"
conda activate BINF6430-ovx-quant        # kallisto 0.48.0 — same env as quant step

REF_DIR=/scratch/lashley.an/final_project/reference
FASTA="$REF_DIR/gencode.vM30.transcripts.fa.gz"
IDX="$REF_DIR/gencode.vM30.kallisto.idx"

mkdir -p "$REF_DIR"

# ASSUMPTION (paper does not state kallisto version or index parameters):
# default k-mer length (k=31), full transcriptome fasta, no masking.
echo "kallisto: $(kallisto version)"
echo "Building index from: $FASTA"

kallisto index -i "$IDX" "$FASTA"

echo "Index built at: $IDX"
echo "DONE"
