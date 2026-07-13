# Rafflesiaceae–Tetrastigma LC–MS Figure Reproducibility Package

This archive contains the Python scripts, input datasets, and figure outputs used to regenerate the LC–MS metabolomics figures for the manuscript:

**Untargeted metabolomics reveals host responses and metabolites linked to host compatibility in Rafflesiaceae parasitism**

Every PNG figure in this package was checked against the corresponding figure embedded in the manuscript.

## Purpose of this package

These scripts are provided to support figure reproducibility for a specific LC–MS metabolomics study of the *Rafflesia–Tetrastigma* host–parasite system. They are not intended as a general-purpose metabolomics analysis pipeline or standalone software tool.

The package allows readers to inspect the input datasets, rerun the figure-generation scripts, and reproduce the plotted outputs used in the manuscript.

## Package structure

```text
.
├── raw_data/
│   ├── rafflesia_dataset1.xlsx
│   ├── rafflesia_dataset2.xlsx
│   ├── rafflesia_dataset3.xlsx
│   ├── rafflesia_dataset4.xlsx
│   └── rafflesia_dataset5.xlsx
├── figure1/
├── figure2/
├── figure3/
├── figure4/
├── figure5/
├── figure6/
├── figure7/
├── LICENSE
└── README.md
```

Each `figureN/` folder contains the script, figure output, and associated CSV or intermediate output files for the corresponding manuscript figure.

## Figure-to-dataset map

| Folder | Manuscript figure | Main input data | Notes |
|---|---|---|---|
| `figure1/` | Figure 1, all-sample PCA | `rafflesia_dataset1.xlsx` | Uses deduplicated LC–MS feature table |
| `figure2/` | Figure 2, selected host/non-host PCA | `rafflesia_dataset1.xlsx` | Uses selected sample groups and recalculated PCA |
| `figure3/` | Figure 3, supergroup PCA | `rafflesia_dataset1.xlsx` | Uses pooled biological supergroups |
| `figure4/` | Figure 4, infection-enrichment significance scatter | `rafflesia_dataset2.xlsx` | Includes nominal p-values and Benjamini–Hochberg FDR values |
| `figure5/` | Figure 5, non-host versus infected comparison | `rafflesia_dataset1.xlsx`; `rafflesia_dataset3.xlsx` | Uses deduplicated raw data and locality-specific comparisons |
| `figure6/` | Figure 6, top-compound stacked bar plot | `rafflesia_dataset4.xlsx` | Uses top-abundance compound summaries |
| `figure7/` | Figure 7, metabolite heatmap | `rafflesia_dataset5.xlsx` | Heatmap generated directly from dataset 5 after removal of `Naringenin` |

## Updates from the developmental version

| Folder | Manuscript figure | Main changes from developmental version |
|---|---|---|
| `figure1/` | Figure 1, all-sample PCA | Removed exact duplicate rows from the feature table and fixed reproducibility settings |
| `figure2/` | Figure 2, selected host/non-host PCA | Recalculated PCA after selected-group filtering; corrected ellipsoid calculation and rendering; removed unsupported non-host ellipsoid |
| `figure3/` | Figure 3, supergroup PCA | Removed exact duplicate rows and corrected group assignment logic |
| `figure4/` | Figure 4, infection-enrichment significance scatter | Added Benjamini–Hochberg FDR values and retained nominal and corrected statistics in the output |
| `figure5/` | Figure 5, non-host versus infected comparison | Recomputed locality-specific Welch tests from deduplicated data, applied FDR correction, and marked non-FDR-supported bars visually |
| `figure6/` | Figure 6, top-compound stacked bar plot | Restored 12-hydroxyoctadecanoic acid, a highly abundant seed-associated compound discussed in the manuscript |
| `figure7/` | Figure 7, metabolite heatmap | Rebuilt the heatmap from `rafflesia_dataset5.xlsx` and removed `Naringenin` because the row could not be verified for the final heatmap |

## Regenerating figures

Each figure folder should be run independently.

To regenerate a figure:

1. Copy the required input dataset from `raw_data/` into the corresponding `figureN/` folder, unless it is already present.
2. Open a terminal in the `figureN/` folder.
3. Run the relevant Python script.

Example:

```bash
cd figure1
python3 figure1_script_name.py
```

The script will regenerate the PNG figure and any associated CSV or summary output files.

## Software requirements

The scripts were written for Python 3 and require the following packages:

```text
pandas
numpy
matplotlib
scikit-learn
openpyxl
```

Some scripts may also require:

```text
scipy
statsmodels
```

A minimal installation command is:

```bash
pip install pandas numpy matplotlib scikit-learn openpyxl scipy statsmodels
```

## Statistical notes

The PCA figures use autoscaled metabolite-feature matrices. In Figures 1–3, duplicate feature rows were removed before PCA. Because PCA variance percentages depend on the exact input matrix and sample set, the percentages reported in each figure correspond to the specific dataset and filtering rules used for that figure.

Figures 4 and 5 include Benjamini–Hochberg FDR correction to account for multiple testing. Nominal p-values and FDR-adjusted values should be interpreted separately. Metabolites or features that do not survive FDR correction should be treated as exploratory candidates rather than confirmed statistically significant markers.

For Figure 5, locality-specific comparisons are shown for selected metabolites. Some comparisons rely on limited biological replication and should therefore be interpreted as exploratory feature-level patterns requiring further validation.

## Figure 7 heatmap note

The Figure 7 heatmap is generated directly from `rafflesia_dataset5.xlsx`. The compound `Naringenin` is removed within the script because it could not be verified for the final heatmap. The script exports the final 33-row heatmap matrix as:

```text
metabolite_heatmap_matrix.csv
```

## Recommended citation

Please cite the associated manuscript once published:

Molina, J., Abzalimov, R., Yin, P., Wicaksono, A., Bürger, M., Hill, J., Bernier, F., Wen, J., & Pell, S. 2026. *Untargeted metabolomics reveals host responses and metabolites linked to host compatibility in Rafflesiaceae parasitism.*

## AI assistance disclosure

Generative AI tools, including ChatGPT and Claude, were used during different stages of this work to assist with Python code drafting, code checking, figure review, and manuscript consistency review. Early script development used ChatGPT GPT-5.2, while later code review, figure checking, and manuscript consistency checks used ChatGPT GPT-5.5 and Claude Sonnet 5.

All scripts, statistical outputs, figure files, biological interpretations, and manuscript conclusions were reviewed and approved by the authors. The authors retain responsibility for the analyses, interpretations, and final reported conclusions.
