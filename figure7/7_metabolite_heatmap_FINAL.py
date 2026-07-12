from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

BASE_DIR = Path(__file__).resolve().parent
MATRIX_FILE = BASE_DIR / "corrected_heatmap_matrix_with_status.csv"

PNG_OUT = BASE_DIR / "metabolite_heatmap_FINAL.png"
SVG_OUT = BASE_DIR / "metabolite_heatmap_FINAL.svg"
CSV_OUT = BASE_DIR / "metabolite_heatmap_matrix_FINAL.csv"

full = pd.read_csv(MATRIX_FILE, index_col=0)

# Naringenin could not be reconstructed/verified from the raw feature matrix
# under any tested compound-matching rule; removed rather than left unverified.
full = full.drop(index="Naringenin")

value_columns = ["Infected ILO", "Uninfected ILO", "Aerial stem/leaf", "Nonhost ILO"]
matrix = full[value_columns].copy()
matrix.to_csv(CSV_OUT)

# ---------- Plot (same logic as the original script) ----------
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

dominant_group = matrix.idxmax(axis=1)
for row in range(1, len(matrix)):
    if dominant_group.iloc[row] != dominant_group.iloc[row - 1]:
        ax.axhline(row - 0.5, color="black", linewidth=0.8)

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

print("Saved:", PNG_OUT.name, SVG_OUT.name, CSV_OUT.name)
print(f"Rows: {len(matrix)} (Naringenin removed)")
plt.close(fig)
