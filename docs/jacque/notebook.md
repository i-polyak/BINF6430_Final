# Jacque -- lab notebook, BINF6430 Final Project
# Bloom et al. 2026 reanalysis
#
#

## 2026-11-27

Context: working through comments from Marcus' 11/16 review ("What I'd do before the final report).

- synced fork with upstream (i-polyak/BINF6430_Final); fast-foraward to 379e549

- Disabled push to upstraem (git remote set-url --push upstream no_phs).

- Deleted untracked scripts/jacque/Rstudio_version/ (Never committed).

- Tagged 78bc39, not 379e549: Marcus' review quotes reports/09_findings_summary.md which 379e549 (11/16) removed.  Confirmed the file is present in the tag.

- History:  my scripts were built and run outside the repo, then copied in later; not all were rerun from the clone.  Originals are scattered:
	/scratch/caldwell.ja/final_project/my_scripts/
	/scratch/caldwell.ja/final_project/kallisto/       (quant.sbatch)
	/scratch/caldwell.ja/final_project/samplekey.txt   (mine: SRR accessions)
  Wrappers tximport_mouse.sh and run_gsea.sh source
  	/scratch/caldwell.ja/final_project/scripts/, which no longer exists.


- Author's repo used for comparison:
	https://github.com/jcbloom/hypo_menopausNavarro.git
	commit commit c220d3566603e50cc9fb91304766115874bd16e7 (2026-06-22)
    local clone: /scratch/caldwell.ja/final_project/src/hypo_menopauseNavarro/
  Their Mouse/samplekey.txt (library names) is not the same file as mine.

- Found in my_scripts/ but not in rep: qdump.sbatch, jc_fastqc_dir.sh,
  m_qc.sh, session_info.txt, R_pkg_versions, install_in_R, start_R421.sh,
  tximport_*.out, earlier sva/deseq2/compare logs (11868xxx).

- Still to locate: samples.txt, index.sbatch, mapping_rate.sh originals.

- diff of kallisto/ originals vs repo:
  index.sbatch, quant.sbatch, fastqc_env.yml, multiqc.yml, kallisto_48.yml identical.
  mapping_rate.sh: repo copy has two lines broken by wrapping during copy-in
  (comment wrapped onto ")l"; jq "@tsv" split). Original parses; it is what
  produced qc_mapping.tsv. Repo copy to be replaced with original.
  => QC version mismatch (yml vs MultiQC report) predates the copy-in.

- found some old notes from building kallisto the first time in the kallisto
  directory, and copied them to docs/jacque/early_notes {log.txt, learning.txt}

- Chose at the time of copying scripts into the repo, not to copy some files 
  from my_Scripts; reasons not recorded

- Renamed original scripts/ folder (now my_scripts/) deliverately, so that
  any wrapper still reference the old location with fail rather than silently
  use an out-of-repo copy.  Stages rerun after the rename (sva, deseq2, compare,
  fgsea) therefore ran from within the repo.  kallisto, tximport and run_gsea
  were not rerun after the rename and still reference the old path.

- File mtimes vs sacct submit times: sva_mouse.R/.sh, authors_sv1.R, fgsea_mouse.sh,
  compare_to_paper.sh, deseq2_mouse.sh all modified before the job that used them
  => committed versions (097c102) are what ran. deseq2_mouse.R mtime is
  2026-11-15 12:25 (after the run) but content unchanged since 097c102 (git status clean).

- 097c102 change to deseq2_mouse.R: input file renamed from
  sampleTable_SVs_<set>_standard.txt to sampleTable_SVs_<set>.txt, following
  removal of the as-written SVA variant (sva_mouse.R, same commit). sva ran
  19:34 and wrote the new names before deseq2 ran at 19:51.

- results/sampleTable_SVs_{matched,allGEO,authorsSV1}.txt written 2026-11-14
  19:35–19:43; no _standard files remain => deseq2 (19:51) read these.
  Provenance of seed-6430 results confirmed: repo scripts at 097c102 are what ran
  for sva, deseq2, compare, fgsea; tximport ran from old folder with scripts
  identical to the repo copies.

- Committed as-ran-recods (16d509f, 89 files) and pushed; opened PR to i-polyak/BINF6430_Final.  mapping_rate.sh mode now 755, and broken copy corrected, matching the orginal run. 

- Merged PR; and notified team that it is there.

### Summary: If you don't document something, it didn't happen. So I have copied all of the log files into the github directory. Addtionally - fixed a bug in a copied script that had an error in it: mapping_rate.sh ; Confirmed provenance of seed-6430 R stages (SVA, DESeq2, comparsion, fgsea), runs with scripts (scripts creation dates were all before runs were done.  for kallisto runs, quant and index files were present, and identical to what is in repo copies, all logs are present and do not report an error.  Tximport; ran from older folder (pre copy into the repo) Confirmed that final run (all items after tximport) were done after the other.   Files involved pushed to github; PR merged; team notified.

## 2026-11-27 (cont.)  Phase B: decision log

- PR #5 (as-ran-records) merged as 482cf4f; started branch decision-log.

- Session timed out while writing decision_log.md; no swap file recovered; rewrote it inside tmux.

- Wrote 4 decision-log entries: docs locations; runs from clone (retire
  copy_to_repo.sh); PROJECT_ROOT/DATA_ROOT; switch to Explorer's standard R.

- Plan: rerun everything from QC onward from the clone; download not rerun.

- tximport 11840148 produced txi.rds (2026-11-12 16:50). Wrapper has a typo'd
  path (final_proect) yet ran; how is unresolved and will be superseded by the rerun.

- fgsea_mouse.sh has the same typo in its --output line; fix in Phase C.

- QC version mismatch (Marcus): multiqc.yml says MultiQC 1.35, fastqc_env.yml says
  FastQC 0.13.0; MultiQC report says 1.26 / 0.12.1. yml files identical to the
  kallisto/ originals, so the mismatch predates the copy-in. m_qc.sh references a
  container — likely QC did not run from these conda envs (not yet confirmed).
  QC rerun from the clone will record versions in the job log.

 
