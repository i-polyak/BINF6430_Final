[View the MultiQC report](https://raw.githack.com/i-polyak/BINF6430_Final/main/results/qc/multiqc_report.html)


# Results overview (Jacque)

Reproduction of the Bloom et al. (JCI Insight) mouse OVX hypothalamus analysis:
kallisto → tximport → SmartSVA → DESeq2 → Hallmark GSEA, compared with the paper's
Table S2 (DE results) and Table S1 (GSEA results). All runs use seed 6430.
Scripts are in `scripts/jacque/`; environment in `environment/jacque/`; logs in `logs/jacque/`.

## Sample sets
- **matched (30):** the authors' samples that could be matched to a GEO SRR.
  Three of the authors' intact controls (L16_5443_S1, L16_5452_S10, L534_3_S3) have no SRR.
- **allGEO (33):** every sample in GSE288245, including three GEO-only samples
  (L534_1, L534_5, L534_10) not found in the authors' sample table.
- **authorsSV1 (30):** the matched set, using the authors' published SV1 in place of ours.

## Runs
| Run | Samples | Design | SV1 |
|---|---|---|---|
| matched_sv1 / _none / _batch | 30 | ~ SV1 + condition / ~ condition / ~ batch + condition | ours |
| allGEO_sv1 / _none / _batch | 33 | same | ours |
| authorsSV1_sv1 | 30 | ~ SV1 + condition | authors' |
| paperS2 | — | authors' DE table (Table S2), GSEA step only | — |

## Where things are
| What | Location |
|---|---|
| SV tables | `deseq2/jacque/sva/` |
| Full DE tables (not in repo, ~230 MB) | Explorer: `/scratch/caldwell.ja/final_project/results/DE_<set>_<design>/`; regenerate with `deseq2_mouse.sh` |
| Gene-level comparison to Table S2 | `comparisons/jacque/<set>_<design>/comparison_to_paper.txt` |
| Hallmark GSEA vs Table S1 | `fgsea/jacque/<set>_<design>/GSEA_stat/`, `GSEA_lfc/` |
| GSEA on the authors' Table S2 | `fgsea/jacque/paperS2/<rank>/` |
| All GSEA runs in one table | `fgsea/jacque/GSEA_all_runs_stat.txt` |
| NES heatmaps (paper vs ours) | `fgsea/jacque/figures/` |
| SV1 diagnostics | `logs/jacque/sva_mouse_11875470.out` |

Each folder has its own README with column definitions and settings.

## Start here
- `fgsea/jacque/GSEA_all_runs_stat.txt`: pathway-level agreement with the paper, every run
- `fgsea/jacque/figures/NES_heatmap_matched_sv1_stat.png`: main reproduction vs the paper
- `comparisons/jacque/matched_sv1/comparison_to_paper.txt`: gene-level agreement, main run

