Readme file for BINF6430_Final/scripts/jacque

# contains:

## after_kallisto_backup.sh -- script for moving a bunch of files that would take more time to calculate if lost.  Saved to /courses partition

## index.sbatch -- SLURM script for creating kallisto index file on the M30 mouse genome.

## quant.sbatch -- SLURM script for creating kallisto files showing the quantification of expressed genes when given a reference genome

## mapping_rate.sh -- takes the output from kallisto (after quant.sbatch is run) and processes it into a file that is filtered and sorted. File that is created is called qc_mapping.tsv.

