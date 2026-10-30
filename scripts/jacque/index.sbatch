#!/bin/bash
#SBATCH --job-name=kal_index
#SBATCH --partition=courses
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/%u/final_project/kallisto/logs/%x_%j.out

set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate kallisto_48

cd /scratch/$USER/final_project/refs
kallisto version

# tx -> gene table (tx_id, gene_id, gene_name) for tximport afterwards....
zcat gencode.vM30.transcripts.fa.gz | grep '^>' | sed 's/^>//' \
	| awk -F'|' 'BEGIN{OFS="\t"}{print $1,$2,$6}' > gencode.vM30.tx2gene.tsv

#strip headers to transcript ID only
zcat gencode.vM30.transcripts.fa.gz | sed '/^>/s/|.*//' | gzip > gencode.vM30.transcripts.clean.fa.gz

kallisto index -i gencode.vM30.k31.idx gencode.vM30.transcripts.clean.fa.gz
