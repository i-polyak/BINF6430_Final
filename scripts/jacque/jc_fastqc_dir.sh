#!/bin/bash
#SBATCH --job-name=JC_First_fastqc                          # Name of the job
#SBATCH -N 1                                                # number of nodes
#SBATCH -c 4                                # number of cpus-per-task (threads)
#SBATCH --mem 32G                           # memory pool for all cores
#SBATCH -t 8:00:00                                          # time (HH:MM:SS)
#SBATCH --mail-type=END,FAIL                # Get an email when the program completes or fails
#SBATCH --mail-user=caldwell.ja@northeastern.edu      # where to send the email
#SBATCH --out=/scratch/caldwell.ja/final_project/logs/%x_%j.log   # captured stdout
#SBATCH --error=/scratch/caldwell.ja/final_project/logs/logs/%x_%j.err # captured stdin
#
#
### Environment ###
module load fastqc         # Now we have access to fastqc


### Variables ##
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

# Full path to the FASTQ file directory. You have to also replace <interesting cell type> with a cell type from the data
FASTQ_DIR=/scratch/caldwell.ja/final_project/data/GSE288245/fastq

# Full path to the FastQC output directory.
OUTPUT_DIR=/scratch/caldwell.ja/final_project/data/GSE288245/qc/fastqc
mkdir -p ${OUTPUT_DIR}

# Create variable for the FastQC flags
FASTQC_FLAGS="--extract --delete -t 4"   # Could add threads to this!

#find subdirectories in FASTQ_DIR and store in local_dirs[]
mapfile -d '' local_dirs < <(find -L ${FASTQ_DIR} -mindepth 1 -type d -printf "%f\0"0)

if [ "${#local_dirs[@]}" -eq 0 ]; then # no subdirectories found
	echo "${FASTQ_DIR} found, running fastqc"
	fastqc ${FASTQC_FLAGS} ${FASTQ_DIR}/*fastq.gz -o ${OUTPUT_DIR}
	echo "all files in ${FASTQ_DIR} processed."
	exit 0
fi

#Loop only runs if local_dirs[] has elements
# For each sub-directory ; 1. create the output subdirectory 
#                          2. find all the files in the subdirectory -type f *.fastqgz
#			   3. Run fastqc saving files in output subdirectory

for l_dir in "${local_dirs[@]}"; do
	echo "Creating output directory $l_dir"
	mkdir -p ${OUTPUT_DIR}/${l_dir}
	# run the fastqc program
	fastqc ${FASTQ_FLAGS} ${FASTQ_DIR}/${l_dir}/*fastq.gz -o ${OUTPUT_DIR}/${l_dir}

done
