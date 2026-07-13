# Python Scripts Visualizer for Host-Parasite Metabolomics Data of Rafflesia

Every PNG in this package was verified byte-identical (MD5) against the image actually embedded
in `Rafflesia_negLCMS_REVISED_trackchanges.docx`. 

## Updates from the developmental version

| Folder | Manuscript figure | Script | What changed vs. original |
|---|---|---|---|
| `figure1/` | Figure 1 (all-sample PCA) | `1_pca_3d_sample_groups_fixed.py` | Dropped 22 exact-duplicate rows from raw data; added `random_state=0` for exact reproducibility |
| `figure2/` | Figure 2 (host/non-host PCA) | `2_FINAL_pooled_68CI_fig1orient.py` | Same dedup fix; ellipsoid radius fixed to true 68% coverage (was ~20%); ellipsoids pooled across localities (was one per locality); non-host ellipsoid removed; equal-aspect-ratio rendering fix; camera angle matched to Figure 1 |
| `figure3/` | Figure 3 (supergroup PCA) | `3_pca_3d_supergroups_fixed.py` | Same dedup fix |
| `figure4/` | Figure 4 (infection-enriched scatter) | `4_infection_enriched_metabolite_scatter_FDR.py` | Added Benjamini–Hochberg FDR correction (previously absent). 0/6 metabolites survive q<0.05 |
| `figure5/` | Figure 5 (non-host vs. infected bar chart) | `5_nonhost_vs_infected_log2fc_barchart_FDR.py` | Rebuilt: recomputed Welch's t-tests across the full feature family per locality from deduplicated raw data, applied FDR correction, added hatching for non-significant bars, p-values shown as decimals (not scientific notation) |
| `figure6/` | Figure 6 (top-20 stacked bar) | `6_stacked_barplot_top20_FIXED.py` | Removed an unexplained exclusion of 12-Hydroxyoctadecanoic acid, which is the #1 most abundant seed compound and is named in the Discussion |
| `figure7/` | Figure 7 (metabolite heatmap) | `7_metabolite_heatmap_FINAL.py` | Rebuilt 33/34 rows from deduplicated data; removed "Naringenin" (row could not be reconstructed/verified from any raw data under any tested matching rule — unrelated to the duplicate-row bug) |
| `raw_data/` | — | — | The 5 source Excel files (`rafflesia_dataset1.xlsx`–`5.xlsx`) needed to run any script above from scratch |

## Regenerating a figure
Each script expects its corresponding `rafflesia_datasetN.xlsx` in the same working directory.
Copy the relevant file from `raw_data/` alongside the script and run with `python3`.

## CSV/data outputs included alongside each figure
Figure 1 → rafflesia_dataset1.xlsx
Figure 2 → rafflesia_dataset1.xlsx
Figure 3 → rafflesia_dataset1.xlsx
Figure 4 → rafflesia_dataset2.xlsx
Figure 5 → rafflesia_dataset1.xlsx + rafflesia_dataset3.xlsx
Figure 6 → rafflesia_dataset4.xlsx
Figure 7 → metabolite_heatmap_matrix_FINAL.csv, plus possibly rafflesia_dataset5.xlsx as source context

## Project
These Python scripts are generated for Rafflesia host-parasite metabolomics project lead by Dr. Jeanmaire Molina.
Manuscript to cite (please cite once it's published):

Molina J, Azbalimov R, Yin P, Wicaksono A, Bernier F, Hill J, Burger M, Wen J, Pell S. 2026. Untargeted metabolomics reveals host responses and metabolites linked to host compatibility in Rafflesiaceae parasitism

## Declaration of AI usage
The Python scripts were generated and checked by assistance of ChatGPT 5.2 and Claude Sonnet 5. All authors have checked, tested, and validated the scripts. 
