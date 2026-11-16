# Environment (Jacque)

## Upstream (QC and quantification): conda environments
- `fastqc_env.yml`: FastQC
- `kallisto_48.yml`: kallisto 0.48.0
- `multiqc.yml`: MultiQC

Created by exporting each environment, for example:

    conda activate multiqc
    conda export -n multiqc --file multiqc.yml --format yaml

Recreate with `conda env create -f <file>.yml`.

## Downstream (tximport, SVA, DESeq2, fgsea): container
- R 4.2.1, Bioconductor 3.15, matching the versions reported in the paper:
  DESeq2 1.36.0, tximport 1.24.0, fgsea 1.22.0
- Container image: `final_project/sif/bioc315.sif`, pulled with `pull_down_r421.sh` (scripts/jacque)
- R library: `final_project/Rlib_bioc315/`
- `run_R421.sh` (scripts/jacque) sets both and defines `rscript_r421`.
  Run `conda deactivate` before sourcing it.
