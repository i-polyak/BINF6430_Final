Decision log — BINF6430 Final (Bloom et al. 2026 reanalysis)

Shared by Ildiko, Anaita, and Jacque. One entry per decision, newest at the bottom.

Record the decision and its evidence, not who caught what.
Don't edit old entries; if a decision changes, add a new entry and fill in "Supersedes".

Personal notebooks: docs/<name>/notebook.md

## Entry template

### <YYYY-MM-DD>  <short title>
- **Status:**       proposed | agreed (date, where)
- **Decision:**
- **Evidence:**     (file, line, log, or table that supports it)
- **Who:**
- **Affects:**      (scripts, results, tables)
- **Supersedes:**   (date + title of earlier entry, or "none")

---

## Entries

### 2026-11-27 Shared documentation locations:
- **Status:**       proposed
- **Decision:**     docs/decision_log.md for group decisions
	      docs/<name>/notebook.md for personal ones
- **Evidence:**     Marcus's 11/16 review, item 10 (dated decision log)
- **Who:**	      Jacque
- **Affects:**      Team
- **Supersedes:**   none

### 2026-11-27 All runs submitted from repo clone
- **Status:** decided
- **Decision:** Every job is submitted from the repo clone and reads scripts only from it.
- **Evidence:** Marcus's item 2
- **Who:** Jacque
- **Affects:** Jacque 
#### R scripts: authors_sv1.R, check_basemean.R, collect_gsea.R, compare_to_paper.R, deseq2_mouse.R, fgsea_mouse.R, nes_heatmap.R, sva_mouse.R, tximport_mouse.R
#### Wrappers that send Slurm logs or output to final_project/results: compare_to_paper.sh, fgsea_mouse.sh, tximport_mouse.sh, plus run_gsea.sh, which points at the old folder
#### README.md
- **Supersedes:** None. Retires scripts/jacque/copy_to_repo.sh

### 2026-11-27 Paths from PROJECT_ROOT and DATA_ROOT; starting with as_ran records 
- **Status:** decided
- **Decision:** No hard coded personal paths in scripts.  PROJECT_ROOT is the clone (submit from repo root); DATA_ROOT holds large files outside git (FASTQs; refs, kallisto output) and each person overrides it with their own path
- **Evidence:** Marcus's item 2
- **Who:** Jacque
- **Affects:** All wrapper scripts,R scripts that build paths.
	Modify:  in scripts/jacque: all files EXCEPT fdump.sbatch.
- **Supersedes:** none
- **Note:** Download of data files will not be run, but the rest of the pipeline will be run end to end (inlcuding the kallisto runs).

### 2026-11-27 Switch from the R 4.2.1 container to Explorer's standard R
- **Status:** proposed
- **Decision:** Use Explorer's standard R module (version to be filled in after
	module avail R) instead of the R 4.2.1 container with a /scratch 
	directory ; R version will be 4.4.1
- **Evidence:** Marcus' recommendation of standard installs over individually
	built stacks; his review notes that the /scratch library is purged
	storage and the container was recorded by tag, not digest.
- **Who:** Jacque
- **Affects:**  All R wrappers and scripts in Jacque's directory; every
	downstream DESeq2/fgsea table; possibly compare with Ildiko's (R 4.2.3) and Anaita's (R 4.2.1) results
- **Supersedes:** original plan to match the paper's R 4.2.1.  Retires scripts in directory: scripts/jacque/{pull_down_r421.sh, run_R421.sh, start_R421.sh}

- **NOTE:** No longer a version-matched reproduction. The existing 4.2.1 seed-6430 results (tag in github pre-marcus-review) stay as the baseline for a continuous comparison (Wald stat, NES) after rerun.
