from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# Upload file
BASE_DIR = Path(__file__).resolve().parent
RAW_FILE = BASE_DIR / "rafflesia_dataset5.xlsx"
SHEET_NAME = "Sheet1"

PNG_OUT = BASE_DIR / "metabolite_heatmap.png"
SVG_OUT = BASE_DIR / "metabolite_heatmap.svg"
CSV_OUT = BASE_DIR / "metabolite_heatmap_matrix.csv"

full = pd.read_excel(RAW_FILE, sheet_name=SHEET_NAME)

# In dataset5, metabolite names are stored in the last unnamed column
if "Unnamed: 4" in full.columns:
    full = full.rename(columns={"Unnamed: 4": "Metabolite"})

value_columns = ["Infected ILO", "Uninfected ILO", "Aerial stem/leaf", "Nonhost ILO"]

required_columns = value_columns + ["Metabolite"]
missing = [c for c in required_columns if c not in full.columns]
if missing:
    raise ValueError(f"Missing required columns: {missing}")

full = full.set_index("Metabolite")

# Naringenin could not be verified from the raw feature matrix
full = full.drop(index="Naringenin", errors="ignore")

matrix = full[value_columns].copy()
matrix = matrix.apply(pd.to_numeric, errors="coerce")

if matrix.isna().any().any():
    bad_rows = matrix[matrix.isna().any(axis=1)]
    raise ValueError(f"Non-numeric or missing values found:\n{bad_rows}")

matrix.to_csv(CSV_OUT)

# ---------- Plot ----------
fig_height = max(8, 0.30 * len(matrix) + 1.6)
fig, ax = plt.subplots(figsize=(10.5, fig_height))
norm = TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2)
image = ax.imshow(matrix.values, aspect="auto", cmap="RdBu_r", norm=norm)
ax.set_xticks(np.arange(len(matrix.columns)))
ax.set_xticklabels(matrix.columns, rotation=45, ha="right", rotation_mode="anchor", fontsize=10)
ax.set_yticks(np.arange(len(matrix.index)))
ax.set_yticklabels(matrix.index, fontsize=8)
ax.set_xticks(np.arange(0.5, len(matrix.columns), 1), minor=True)
ax.set_yticks(np.arange(0.5, len(matrix.index), 1), minor=True)
ax.grid(which="minor", color="white", linestyle="-", linewidth=0.75)
ax.tick_params(which="minor", bottom=False, left=False)
ax.tick_params(axis="both", length=0)

for spine in ax.spines.values():
    spine.set_visible(False)

colorbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.05)
colorbar.set_label("Row-wise Z-score", fontsize=10)
colorbar.ax.tick_params(labelsize=9)
colorbar.set_ticks([-2, -1, 0, 1, 2])

ax.set_title("Metabolite Heatmap", fontsize=14, pad=12)
fig.tight_layout()
fig.savefig(PNG_OUT, dpi=300, bbox_inches="tight", pad_inches=0.2)
fig.savefig(SVG_OUT, bbox_inches="tight", pad_inches=0.2)
print("Saved:", PNG_OUT.name, SVG_OUT.name)
print(f"Rows: {len(matrix)}")
plt.close(fig)
