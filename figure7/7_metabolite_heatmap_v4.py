"""Fig. 7: row-wise Z-score heatmap of selected metabolites in the R. speciosa system (ILO).
Group means are computed directly from dataset1 (raw_data/heatmap_metabolites.csv gives the list and order of
metabolites, selected by pathway relevance). Each metabolite's mean is taken over all
features carrying that annotation (case-insensitive; '.mol' suffix removed)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from pubstyle import apply_style, WIDTH
from common import load_raw, group_mean
apply_style()
BASE = Path(__file__).resolve().parent
GROUPS = {"infectedILO": "Infected host\ntissue (n = 9)", "UNinfectedILO": "Uninfected host\ntissue (n = 10)",
          "uninfecraffspec-stemleaf": "Uninfected\nstem/leaf (n = 2)", "nonhostILO": "Non-host\nroot (n = 4)"}

raw = load_raw()
d5 = pd.read_csv(BASE.parent / "raw_data" / "heatmap_metabolites.csv")
names = d5["Metabolite"].dropna().astype(str).str.strip().str.replace(r"\.mol$", "", regex=True)
names = [n for n in names if n != "Naringenin"]          # could not be verified in the raw matrix
means = pd.DataFrame({g: group_mean(raw, g) for g in GROUPS}); means["key"] = raw["name_key"]
means = means.groupby("key").mean()
missing = [n for n in names if n.lower() not in means.index]
assert not missing, f"metabolites not found in dataset1: {missing}"
M = means.loc[[n.lower() for n in names]]; M.index = names
Z = M.sub(M.mean(axis=1), axis=0).div(M.std(axis=1, ddof=0), axis=0)
Z.to_csv(BASE / "metabolite_heatmap_matrix.csv"); M.to_csv(BASE / "metabolite_heatmap_group_means.csv")

fig, ax = plt.subplots(figsize=(WIDTH, 0.16 * len(Z) + 1.4))
im = ax.imshow(Z.values, aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2))
ax.set_xticks(range(len(GROUPS))); ax.set_xticklabels(GROUPS.values(), fontsize=7)
ax.set_yticks(range(len(Z))); ax.set_yticklabels([n.replace("?-D", "β-D") for n in Z.index], fontsize=7)
ax.set_xticks(np.arange(0.5, len(GROUPS), 1), minor=True); ax.set_yticks(np.arange(0.5, len(Z), 1), minor=True)
ax.grid(which="minor", color="white", lw=0.75); ax.tick_params(which="both", length=0)
for s in ax.spines.values(): s.set_visible(False)
cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.05); cb.set_label("Row-wise Z-score", fontsize=7.5); cb.ax.tick_params(labelsize=7); cb.set_ticks([-2, -1, 0, 1, 2])
fig.tight_layout(); fig.savefig(BASE / "metabolite_heatmap.png", dpi=600, bbox_inches="tight", pad_inches=0.1)
print(f"Rows: {len(Z)}"); print(Z.round(2).to_string())
