# Rafflesiaceae–Tetrastigma LC–MS Figure Reproducibility Package

This repository contains the Python scripts, input datasets, and figure outputs used to regenerate the LC–MS metabolomics figures for the manuscript:

**Untargeted metabolomics reveals host responses and metabolites linked to host compatibility in Rafflesiaceae parasitism**

Every PNG figure in this package was checked against the corresponding figure embedded in the manuscript.

## Purpose of this package

These scripts are provided to support figure reproducibility for a specific LC–MS metabolomics study of the *Rafflesia–Tetrastigma* host–parasite system. They are not intended as a general-purpose metabolomics analysis pipeline or standalone software tool.

The package allows readers to inspect the input datasets, rerun the figure-generation scripts, and reproduce the plotted outputs used in the manuscript.

## Package structure

```text
.
├── raw_data/            rafflesia_dataset1.xlsx (main feature matrix = Supplementary Dataset S1; source of all figures)
│                        Supplementary_Dataset_S2_sample_metadata.xlsx (sample metadata = Supplementary Dataset S2)
│                        heatmap_metabolites.csv (Fig. 7 metabolite list, selected by pathway relevance)
├── supplementary_datasets/  Supplementary Dataset S1 (feature matrix = raw_data/rafflesia_dataset1.xlsx),
│                        Dataset S2 (sample metadata for every run), Table S1 (copy of figure5 output)
├── figure1/ … figure7/  script, figure output, and CSV outputs for each manuscript figure
├── supplementary/       Supplementary Figs. S1–S4 (2D PCA projections, scree plots, PCA without Ampelopsis/seeds)
├── pubstyle.py          shared print style and group labels
├── common.py            shared raw-data loading, name cleaning, group means
├── verify_figures.py    independent check of figure CSVs against raw data
├── statistics/          PERMANOVA on biological-sample means (results in permanova_results.txt)
├── run_scripts.sh
├── LICENSE
└── README.md
```

## Figure-to-script map (revision 2)

| Folder | Manuscript item | Script | Main input |
|---|---|---|---|
| `figure1/` | Fig. 1, all-sample 3D PCA | `1_pca_3d_sample_groups_v3.py` | dataset1 |
| `figure2/` | Fig. 2, host/non-host 3D PCA with pooled 68% ellipsoids | `2_PCA_pooled_v4.py` | dataset1 |
| `figure3/` | Fig. 3, supergroup 3D PCA | `3_pca_3d_supergroups_v3.py` | dataset1 |
| `figure4/` | Fig. 4, infected vs apparently uninfected (biological samples) | `4_infection_candidates_scatter_bio.py` (+ `bio_stats.py`) | dataset1 |
| `figure5/` | Fig. 5 and Table S1, non-host vs infected (biological samples) | `5_nonhost_vs_infected_dotplots_bio.py` (+ `bio_stats.py`) | dataset1 |
| `figure6/` | Fig. 6, top-20 compound stacked bars | `6_stacked_barplot_top20_v5.py` | dataset1 |
| `figure7/` | Fig. 7, metabolite heatmap | `7_metabolite_heatmap_v4.py` | dataset1 (+ `heatmap_metabolites.csv`) |
| `supplementary/` | Figs. S1–S4 | `pca_2d_supplementary.py` | dataset1 |
| `statistics/` | PERMANOVA (Results 3.2–3.3) | `permanova.py` (+ `bio_stats.py`) | dataset1 |

## Changes in revision 2

- **Biological replication (Figs. 4–5, PERMANOVA):** the two technical LC-MS/MS injections of each sample (adjacent columns, e.g. `infectedILO` / `infectedILO.1`) are averaged per biological sample before testing. Earlier versions tested technical runs as independent observations. Welch's t-tests; log2FC = log2((mean A + 1)/(mean B + 1)); negative intensities set to 0; Benjamini–Hochberg FDR across all testable metabolites/features per contrast (606 name-consolidated metabolites, pooled Fig. 4; 1,369 features CAM and 1,476 features ILO, Fig. 5). No metabolite or feature has q < 0.05.
- **Fig. 4** is computed directly from dataset1 (the earlier version read precomputed run-level statistics from a separate summary table). All metabolites meeting p < 0.05 and log2FC > 1 at the biological-sample level are labeled.
- **Fig. 5** is selected by a stated rule: log2FC > 1 (non-host/infected) at both localities and p < 0.05 in at least one, natural-product annotations only. All qualifying features are in `figure5/fig5_rule_based_candidates.csv` (Supplementary Table S1); the nine with the lowest p are plotted as biological-sample dot plots, labeled by retention time.
- **PERMANOVA:** Euclidean distances on autoscaled log10(x + 1) values, 9,999 permutations, restricted within locality for pooled contrasts; Bray–Curtis as sensitivity check. Non-host vs infected p = 0.002 (R² = 0.11); infected vs apparently uninfected p = 0.18.
- **Fig. 2:** no ellipsoid is drawn for the two-sample aerial stem/leaf group (it was not visible in the figure).
- **Figs. 6–7 data correction:** group means are now computed directly from the individual runs in dataset1. The precomputed `aveuninfecraffspec-stemleaf` column (used via dataset4/dataset5 in earlier versions) did not equal the mean of its runs; this affected the stem/leaf bar in Fig. 6 and all Z-scores in Fig. 7. The heatmap metabolite list is in `raw_data/heatmap_metabolites.csv`; the earlier intermediate tables (datasets 2–5) are no longer used and have been removed.
- **Name consolidation:** annotations differing only in capitalization or a ".mol" suffix (e.g., "EPICATECHIN"/"Epicatechin", "CITRIC ACID"/"Citric acid") are merged in all name-consolidated analyses (Figs. 4, 6, 7); the pooled Fig. 4 family is 606 metabolites.
- **Verification:** `run_scripts.sh` ends by running `verify_figures.py`, which recomputes Figs. 4–7 values independently from dataset1 and checks that every metabolite meeting the Fig. 4 thresholds is labelled.
- **Fig. 6** ranks the 20 most abundant compounds across all annotated features in dataset1 (summed group means), excluding non-natural or unresolved annotations (`common.is_natural`); earlier versions ranked only a hand-curated subset (former dataset4) that omitted abundant natural products such as isovitexin and L-malic acid.
- **Fig. 7:** exported at 600 dpi.
- **Supplementary Figs. S1–S4:** PC1–PC3 and PC2–PC3 projections and scree plots for Figs. 1–3, and Fig. 1 recomputed without Ampelopsis and R. speciosa seeds.
- All figures are drawn at final print width (174 mm) with 6.5–8 pt lettering and plain-language group labels, using the shared `pubstyle.py`; PNGs are exported at 600 dpi.
- All scripts read data from `../raw_data/` using platform-independent paths.

## Regenerating figures

Run each script from inside its own folder, e.g.

```bash
cd figure4
python3 4_infection_candidates_scatter_bio.py
```

or run every script with `run_scripts.sh` from the repository root (Linux: `./run_scripts.sh`; Mac: `sh run_scripts.sh`; Windows: `bash run_scripts.sh` in Git Bash).

### Automating script execution

Additionally, you can run every script at once by running the Bash script `run_scripts.bh`.
To do so:

1. Run your terminal or console in the parent folder (folder where `run_scripts.bh` is located).
2. (Linux) Run the command `./run_scripts.bh`.
3. (Mac) Run the command `sh run_scripts.bh`.
4. (Windows) Install an application to run Unix commands (such as Git Bash) and run the command `bash run_scripts.bh`.

## Software requirements

The scripts were written for Python 3 and require the following packages:

```text
pandas
numpy
matplotlib
scikit-learn
openpyxl
```

and `scipy`. A minimal installation command is:

```bash
pip install pandas numpy matplotlib scikit-learn openpyxl scipy
```

## Supplementary datasets

- **Dataset S1** (`raw_data/rafflesia_dataset1.xlsx`): the full feature matrix (1,540 rows; 114 sample runs and 4 procedural blanks). Technical duplicate injections are adjacent columns (e.g., `infectedILO`, `infectedILO.1`).
- **Dataset S2** (`raw_data/Supplementary_Dataset_S2_sample_metadata.xlsx`): one row per run in Dataset S1, with biological sample ID, technical replicate, group, locality, tissue, species, infection status, and collection date. *Tetrastigma* species are difficult to identify without leaves or reproductive structures; where several species are listed for a group, the species of individual samples could not be confirmed.
- **Table S1** (`figure5/Table_S1_nonhost_vs_infected.xlsx`): generated by the Fig. 5 script.

## Supplementary datasets

- **Dataset S1** (`supplementary_datasets/Supplementary_Dataset_S1_feature_matrix.xlsx`) is identical to `raw_data/rafflesia_dataset1.xlsx`.
- **Dataset S2** (`supplementary_datasets/Supplementary_Dataset_S2_sample_metadata.xlsx`) describes every run column in Dataset S1: biological sample ID, technical replicate, group, locality, tissue, species, infection status, and collection date. Technical replicates are adjacent columns in Dataset S1. *Tetrastigma* species are difficult to identify without leaves or reproductive structures; where several species are listed for a group, the exact species of individual samples could not be confirmed.
- **Table S1** is generated by `figure5/5_nonhost_vs_infected_dotplots_bio.py`; the copy in `supplementary_datasets/` is for convenience.

## Note on rafflesia_dataset1.xlsx

The raw matrix includes precomputed group-average columns (`ave…`). These are ignored by every script; all means are computed from the individual runs. The `aveuninfecraffspec-stemleaf` column does not equal the mean of its four runs and should not be used.

## Statistical notes

PCA (Figs. 1–3, S1–S4) uses autoscaled feature matrices and plots individual technical runs for analytical reproducibility; points are not independent biological replicates. All hypothesis tests (Figs. 4–5, PERMANOVA) use biological-sample means. Nominal p-values and FDR-adjusted q-values should be interpreted separately; highlighted metabolites are exploratory candidates, not confirmed markers. CAM comparisons rest on two biological samples per group.

Statistical tests are run on all annotated features. Annotations judged unlikely to be natural products were excluded only from figure labels and interpretation (screening with ChatGPT followed by manual author review).

## Figure 7 heatmap note

The Figure 7 heatmap is computed from `rafflesia_dataset1.xlsx` for the 33 metabolites listed in `raw_data/heatmap_metabolites.csv` (Naringenin, which could not be verified in the raw matrix, is excluded). The script exports the final 33-row heatmap matrix as:

```text
metabolite_heatmap_matrix.csv
```

## Recommended citation

Please cite the associated manuscript once published:

Molina, J., Abzalimov, R., Yin, P., Wicaksono, A., Bürger, M., Hill, J., Bernier, F., Wen, J., & Pell, S. 2026. *Untargeted metabolomics reveals host responses and metabolites linked to host compatibility in Rafflesiaceae parasitism.*

## AI assistance disclosure

Generative AI tools, including ChatGPT and Claude, were used during different stages of this work to assist with Python code drafting, code checking, natural-product screening of annotations for display, figure review, and manuscript consistency review. Early script development used ChatGPT GPT-5.2, while later code review, figure checking, and manuscript consistency checks used ChatGPT GPT-5.5 and Claude Sonnet 5.

All scripts, statistical outputs, figure files, biological interpretations, and manuscript conclusions were reviewed and approved by the authors. The authors retain responsibility for the analyses, interpretations, and final reported conclusions.
