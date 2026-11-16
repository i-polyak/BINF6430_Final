#!/bin/bash
# =============================================================================
# kallisto_quant.sbatch
# Quantify all GSE288245 samples (Bloom et al. 2026 reproduction) with kallisto.
# Submit from the LOGIN node:  sbatch scripts/anaita/kallisto_quant.sbatch
# =============================================================================

# --- SLURM resource request (read by sbatch; %j expands to the job ID) --------
#SBATCH --job-name=kallisto_quant
#SBATCH --partition=courses
#SBATCH --time=06:00:00          # wall-time ceiling; generous headroom for 33 samples
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16       # cores available to kallisto's threads
#SBATCH --mem=24G                # index (~2.6 GB) + overhead; memory is not the bottleneck
#SBATCH --output=/scratch/lashley.an/final_project/logs/quant_%j.out
#SBATCH --error=/scratch/lashley.an/final_project/logs/quant_%j.err

# Fail fast: exit on first error (-e), error on unset variables (-u), and make a
# pipeline fail if any stage fails (-o pipefail). A broken sample stops the job
# instead of silently continuing.
set -euo pipefail

# A batch job runs in a non-interactive shell where `conda activate` is not set
# up by default; this initializes conda so activate works inside sbatch.
eval "$(conda shell.bash hook)"
conda activate BINF6430-ovx-quant        # environment holding kallisto 0.48.0

# --- Paths --------------------------------------------------------------------
IDX=/scratch/lashley.an/final_project/reference/gencode.vM30.kallisto.idx   # index (GENCODE vM30)
READS=/scratch/caldwell.ja/final_project/data/GSE288245/fastq              # paired FASTQs (teammate's scratch, read-only)
OUT=/scratch/lashley.an/final_project/quant/GSE288245                      # my outputs, one folder per sample
mkdir -p "$OUT"

# Record what actually ran, into the job log.
echo "kallisto: $(kallisto version)"
echo "samples: $(ls "$READS"/*_1.fastq.gz | wc -l)"

# --- Quantify every sample ----------------------------------------------------
# Loop over each read-1 file; for each one, derive the sample ID and hand
# kallisto the matching read-1 / read-2 pair.
for r1 in "$READS"/*_1.fastq.gz; do
  # basename strips the folder and the _1.fastq.gz suffix, leaving the SRR
  # accession, e.g. /.../SRR32153993_1.fastq.gz -> SRR32153993
  sample=$(basename "$r1" _1.fastq.gz)
  echo "=== $sample ==="

  # kallisto quant arguments:
  #   -i  transcriptome index to pseudoalign against
  #   -o  per-sample output folder (writes abundance.tsv, abundance.h5, run_info.json)
  #   -t  threads = the cores SLURM allocated to this job
  #   then the two paired-end read files, R1 then R2
  # No strand flag  -> matches the paper's unstated kallisto default (faithful reproduction)
  # No -b bootstraps -> DESeq2 uses raw estimated counts, so bootstraps aren't needed
  kallisto quant -i "$IDX" -o "$OUT/$sample" -t "$SLURM_CPUS_PER_TASK" \
    "$r1" "$READS/${sample}_2.fastq.gz"
done

echo "DONE"   # printed only if every sample finished (because of set -e)
