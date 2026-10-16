# RNA-seq Reproduction — Environment, QC, Index, and Quantification

**Author:** Anaita Lashley
**Project:** BINF 6430 Final Project (Group 1) — reproducibility study
**Paper:** Bloom et al. (2026), *JCI Insight* 11(17):e207270 — progressive hypothalamic neuroinflammation in ovariectomized mice
**Dataset:** GEO GSE288245 — mouse bulk RNA-seq, posterior hypothalamus (PH) + preoptic area (POA), 33 paired-end samples from 17 mice
**Scope reproduced:** mouse bulk RNA-seq pipeline (Fig. 1D, Fig. 2)
**HPC:** Northeastern Explorer
**Last updated:** 2026-10-15

---

## Overview

This log records my portion of the pipeline so far: HPC setup, conda environments,
quality-control review, building the kallisto index, and quantifying all 33 samples. Raw
reads are read from a teammate's scratch directory (not re-downloaded); all of my outputs
are written to my own scratch. `tximport`, DESeq2, and GSEA follow in later phases.

A recurring operational rule on Explorer: the **login node** has internet access but a
tight per-user memory cap, while **compute nodes** have memory and CPU but no outbound
internet. So anything that downloads runs on the login node, and anything memory-heavy
(conda solves, the kallisto index build) runs on a compute node. Memory-heavy jobs run
on the login node are killed (exit code 137). A related gotcha: an interactive node
request (`srun`) cannot be nested inside another, and a non-interactive batch job needs
conda initialized explicitly before `conda activate` works (see §6).

---

## Open questions and assumptions

These are the decisions made in the absence of complete information from the paper, and
the questions still open. They are collected here because they are the points most likely
to affect whether our results match the authors' — and the first thing a reproducibility
review should see. Each is revisited in the relevant section below.

- **kallisto version — assumed 0.48.0.** The paper does not report it. The downstream
  package versions it does give (`tximport` 1.24.0, DESeq2 1.36.0, `fgsea` 1.22.0) belong
  to the Bioconductor 3.15 cycle (spring 2022), and kallisto 0.50.0 (mid-2023) changed the
  index format and dropped bulk RNA-seq features — so a pre-0.50 version is appropriate and
  0.48.0 is the best-supported choice. Documented as an assumption.
- **kallisto index settings — default k-mer (31) assumed.** The paper states no index
  parameters; defaults are used. Build fingerprint for cross-checking: 897,905 contigs /
  123,081,114 k-mers (2.6 GB index) from `gencode.vM30.transcripts.fa.gz`.
- **Strandedness — decided for the primary run; alternative kept as a check.** Not stated
  in the paper. The primary quantification uses **no strand flag**, matching the paper's
  unstated kallisto default (the faithful reproduction). The library is reverse-stranded
  (KAPA mRNA HyperPrep Stranded kit), so a `--rf-stranded` run is retained as a documented
  sensitivity check rather than used as the primary setting — adding a flag the authors
  never specified would itself be a deviation.
- **Sample-to-condition mapping — unresolved.** The GEO labels (4 PH sham / 6 PH 14-day)
  conflict with the analysis code (5 / 5), and GEO uses "sham" where the paper uses
  "intact." Each SRR run must be mapped to its condition via the authors' sample key, and
  any unreported sample exclusions checked, before differential expression. (Quantification
  does not need this — it runs per sample — so quant proceeds in parallel with resolving it.)
- **Trimming — assumed none.** The paper describes no trimming step, and the deposited
  reads look untrimmed (uniform 50 bp length, low-complexity runs retained). We reproduce
  without trimming; a trim-vs-no-trim comparison is a candidate sensitivity check.
- **Downstream modeling — to confirm in the DE phase.** We intend to use the published
  design (`~ SV1 + condition`) with the authors' surrogate-variable values; the paper's
  limited batch modeling — and the observation that SV1 largely separates the two brain
  regions — is a concern to revisit when reproducing Fig. 1D / Fig. 2.

---

## 1. HPC setup — interactive compute-node aliases

Two aliases were added to `~/.bashrc` to request an interactive compute node, so work
does not run on the login node:

```bash
# light work (FastQC, index build)
alias getnode='srun --time=08:00:00 --ntasks=1 --nodes=1 --mem=16G --cpus-per-task=16 --partition=courses --pty bash'

# heavy work (quantifying many samples)
alias getbignode='srun --time=16:00:00 --ntasks=1 --nodes=1 --mem=128G --cpus-per-task=46 --partition=courses --pty bash'
```

After editing `~/.bashrc`, load them into the current shell:

```bash
source ~/.bashrc
```

Note: these request *interactive* nodes and cannot be nested (running `getbignode` from
inside a `getnode` session fails). Long unattended work is submitted with `sbatch` instead
(see §6).

---

## 2. Conda environments

Pinning exact tool versions is the point of this study: the paper does not provide an
environment file, which is why the kallisto version had to be reconstructed (see below).
Environment files are exported and committed to the repo so each group member rebuilds an
identical environment.

**Channel configuration** (per Bioconda's documented setup): channels ordered
`conda-forge` (highest priority), then `bioconda`, with strict channel priority. The
`defaults` channel is not used — Bioconda removed it from the recommended set (Aug 2024),
and avoiding it also sidesteps Anaconda's licensing terms for institutional use.

### QC environment (reused)

An existing general-purpose QC environment was reused: **FastQC 0.12.1**, **MultiQC 1.19**
(bioconda). On Explorer, FastQC is also available via `module load fastqc`; both routes
are recorded so any difference is traceable. QC-tool versions affect the report's
presentation, not the counts that flow downstream.

### Quantification environment (built)

**kallisto 0.48.0.**

*Version rationale (assumption):* the paper does not state the kallisto version. The
downstream package versions it does report — `tximport` 1.24.0, DESeq2 1.36.0,
`fgsea` 1.22.0 — belong to the Bioconductor 3.15 cycle (spring 2022), and kallisto 0.48.0
was the release current at that time. kallisto 0.50.0 (mid-2023) introduced a new index
format and dropped bulk RNA-seq features, so a pre-0.50 version is appropriate for
reproducing a bulk RNA-seq paper. We therefore assume **0.48.0** and record it as a
stated assumption.

Build (run on a compute node; strict channel priority set first):

```bash
conda config --set channel_priority strict
getnode
conda create -n BINF6430-ovx-quant -c conda-forge -c bioconda kallisto=0.48.0
conda activate BINF6430-ovx-quant
kallisto version          # confirm 0.48.0
```

Export for the repo:

```bash
conda env export --no-builds > environment/BINF6430-ovx-quant.yml
```

*(A static kallisto 0.48.0 Linux binary from the pachterlab v0.48.0 release is an
equivalent fallback if conda cannot solve the environment on a given machine.)*

---

## 3. Data

Raw reads: GSE288245, 33 paired-end FASTQ files (50 bp, NovaSeq), obtained from SRA.
The download was performed by a teammate (Jacqueline); I read the reads directly from her
scratch rather than re-downloading, to avoid duplicating a large (~117 GB) transfer. The
FASTQs sit in a `fastq/` subfolder:

```
/scratch/caldwell.ja/final_project/data/GSE288245/fastq/
```

Files are named `SRR<run>_1.fastq.gz` / `SRR<run>_2.fastq.gz` (fasterq-dump split pairs),
spanning SRR32153993–SRR32154025. Confirmed the full set before quantifying:

```bash
ls /scratch/caldwell.ja/final_project/data/GSE288245/fastq/*_1.fastq.gz | wc -l   # -> 33
```

Note: downloading from SRA to Explorer is done one accession at a time because of the
cluster's egress restrictions.

---

## 4. Quality control

FastQC + MultiQC were run by a teammate (Jacqueline); I reviewed the MultiQC report.
The report is copied into the repo at `results/qc/multiqc_report.html` and linked from
`results/README.md` via raw.githack, because GitHub does not render HTML inline.

Observations:

- **GC content:** 7 samples warned. The warnings fall on the **R1 (read-1)** files, which
  is consistent with read-1 composition bias (e.g. random-hexamer priming) rather than a
  problem with specific libraries.
- **Overrepresented sequences:** 8 samples failed. Most top sequences are low-complexity
  homopolymer runs (poly-A, poly-C/T) — expected artifacts of mRNA selection, and the
  reason such sequences return "No Hit" on BLAST.
- The highest-abundance *complex* sequence
  (`GGATTTAGCTCAGTGGTAGAGCGCTTGCCTAGCAAGCGCAAGGCCCTGGG`) carries a conserved tRNA sequence
  motif and appears across all samples — most likely a tRNA fragment (to be confirmed
  against GtRNAdb/Rfam rather than general NCBI BLAST, which is unreliable for short,
  structured sequences).
- Warnings/failures in these two modules are expected for RNA-seq, whose libraries are
  non-uniform by design; they do not indicate failed sequencing runs.

**Trimming decision:** the paper describes no trimming step before kallisto, and the
deposited reads show the hallmarks of untrimmed data (uniform 50 bp length; low-complexity
runs retained). To reproduce faithfully, no trimming is performed. A trim-vs-no-trim
sensitivity check is a candidate Phase 3 experiment.

---

## 5. Reference and kallisto index

**Reference:** GENCODE **vM30** mouse transcripts, full transcript set (all biotypes), on
**GRCm39/mm39** — matching the paper's annotation and build. The full transcripts file is
used (not the protein-coding-only file) because the paper pseudoaligned on "the entire
transcriptome annotation" and only restricted to protein-coding genes downstream.

Download (login node — needs internet):

```bash
mkdir -p /scratch/lashley.an/final_project/reference
cd /scratch/lashley.an/final_project/reference
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M30/gencode.vM30.transcripts.fa.gz
```

Build the index (compute node — memory):

```bash
getnode
conda activate BINF6430-ovx-quant
cd /scratch/lashley.an/final_project/reference
kallisto index -i gencode.vM30.kallisto.idx gencode.vM30.transcripts.fa.gz
ls -lh gencode.vM30.kallisto.idx
```

k-mer length: default (**31**), appropriate for 50 bp reads. The same kallisto (0.48.0) is
used for the index and the quantification so the index format matches.

Build result (2026-10-15): 897,905 contigs / 123,081,114 k-mers; index file 2.6 GB. These
figures double as a cross-check — a teammate building from the same FASTA with kallisto
0.48.0 should reproduce them.

---

## 6. Quantification

All 33 samples were quantified with `kallisto quant` (paired-end), reading the index from
my scratch and the FASTQs from the teammate's `fastq/` folder, and writing one output
folder per sample to my own scratch:

```
/scratch/lashley.an/final_project/quant/GSE288245/
```

**Decisions:**

- **No strand flag** for the primary run — matches the paper's unstated kallisto default
  (see Open questions). A `--rf-stranded` run is kept as a separate sensitivity check.
- **No bootstraps** (`-b` omitted) — DESeq2 uses raw estimated counts, which do not need
  bootstrap replicates (those matter for transcript-level uncertainty methods such as
  sleuth).
- Threads tied to the job's allocated cores.

Each sample folder contains `abundance.tsv`, `abundance.h5`, and `run_info.json`.

The run is submitted as a **batch job** rather than run interactively, so it survives a
dropped session and avoids the compute-node nesting issue. The script is
`scripts/anaita/kallisto_quant.sbatch`. One batch-specific detail: a non-interactive SLURM
shell needs conda initialized explicitly (`eval "$(conda shell.bash hook)"`) before
`conda activate` will work.

```bash
#!/bin/bash
# --- SLURM resource request (read by sbatch; %j is the job ID) ---------------
#SBATCH --job-name=kallisto_quant
#SBATCH --partition=courses
#SBATCH --time=06:00:00          # wall-time ceiling; generous headroom for 33 samples
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16       # cores for kallisto's threads
#SBATCH --mem=24G                # index (~2.6 GB) + overhead; memory is not the constraint
#SBATCH --output=/scratch/lashley.an/final_project/logs/quant_%j.out
#SBATCH --error=/scratch/lashley.an/final_project/logs/quant_%j.err

# Fail fast: exit on first error (-e), error on unset variables (-u), and make
# a pipeline fail if any stage fails (-o pipefail) — so a broken sample stops
# the job instead of silently continuing.
set -euo pipefail

# A batch job runs in a non-interactive shell where `conda activate` is not set
# up by default; this initializes conda so activate works inside sbatch.
eval "$(conda shell.bash hook)"
conda activate BINF6430-ovx-quant        # environment holding kallisto 0.48.0

# --- Paths -------------------------------------------------------------------
IDX=/scratch/lashley.an/final_project/reference/gencode.vM30.kallisto.idx   # index (GENCODE vM30)
READS=/scratch/caldwell.ja/final_project/data/GSE288245/fastq              # paired FASTQs (teammate's scratch, read-only)
OUT=/scratch/lashley.an/final_project/quant/GSE288245                      # my outputs, one folder per sample
mkdir -p "$OUT"

# Record what actually ran, into the job log.
echo "kallisto: $(kallisto version)"
echo "samples: $(ls "$READS"/*_1.fastq.gz | wc -l)"

# --- Quantify every sample ---------------------------------------------------
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
```

Submitted from the login node (no node needed to queue a job):

```bash
mkdir -p /scratch/lashley.an/final_project/logs
sbatch scripts/anaita/kallisto_quant.sbatch
squeue -u lashley.an                                          # watch the queue
tail -f /scratch/lashley.an/final_project/logs/quant_*.out    # watch progress
```

**Verification (on completion):** check the pseudoalignment rate for every sample —

```bash
grep p_pseudoaligned /scratch/lashley.an/final_project/quant/GSE288245/*/run_info.json
```

Consistent, reasonable rates with no outliers indicate a clean run.

**Status (2026-10-15): complete.** Job 11115023 finished `COMPLETED` (exit `0:0`) in
1h30m, and all 33 samples produced output. Pseudoalignment rates ranged **85.2%–92.1%**
(most ~88–92%), tightly clustered with no outliers — the two lowest (SRR32153997 at 85.2%,
SRR32154004 at 86.5%) are well within normal spread, not anomalies. The uniformly high
rates also corroborate that the reference (GENCODE vM30 / mm39) and kallisto 0.48.0 are
well matched to the data; a wrong build or annotation would show up as depressed, scattered
rates. The abundance outputs are the artifact to preserve (copied out of scratch before the
purge window).

---

## 7. File locations / repo policy

- **Repo (GitHub):** scripts, environment `.yml` files, small results (MultiQC report, this
  log).
- **Scratch (not in repo — large and purgeable):** raw reads (read from teammate's
  scratch), reference FASTA, kallisto index, quantification outputs.
- **Pattern:** read inputs from the teammate's scratch; write all outputs to my own scratch.

---

## 8. Recorded versions and parameters

| Component | Value |
|---|---|
| Genome build | GRCm39 / mm39 |
| Annotation | GENCODE vM30 (full transcripts, all biotypes) |
| kallisto | 0.48.0 (assumed — see §2) |
| kallisto index k-mer | 31 (default) |
| kallisto quant strand | none / unstranded (primary); `--rf-stranded` kept as sensitivity check |
| kallisto quant bootstraps | none (`-b` omitted) |
| FastQC | 0.12.1 (conda); module version recorded separately |
| MultiQC | 1.19 |
| Reads | 50 bp paired-end (NovaSeq), GSE288245 |

Downstream versions reported by the paper (for the quant/DE phases): `tximport` 1.24.0,
DESeq2 1.36.0, `smartSVA` 0.1.3, `fgsea` 1.22.0, R 4.2.1.

---

## 9. Next steps

- Check pseudoalignment rates across all 33 samples once the batch job finishes (quant QC).
- Resolve the GEO sample-label vs analysis-code discrepancy (GEO lists 4 PH sham / 6 PH
  14-day; the code uses 5 / 5), mapping each SRR run to its condition via the authors'
  sample key.
- Build a `tx2gene` map from GENCODE vM30 and import counts with `tximport`.
- Run DESeq2 (`~ SV1 + condition`) using the authors' surrogate-variable values; compare
  protein-coding DEG counts to Fig. 2.
- Run `fgsea` on the Hallmark gene sets; compare to Fig. 1D.
- Optional sensitivity checks: `--rf-stranded` quantification; trim-vs-no-trim.
