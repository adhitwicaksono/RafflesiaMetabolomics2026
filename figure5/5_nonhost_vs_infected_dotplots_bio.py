"""Fig. 5: non-host vs infected Tetrastigma per locality, biological-sample level.
Selection rule (same nominal criteria as Fig. 4, applied per locality): log2FC > 1 (non-host/infected)
at BOTH localities and p < 0.05 in at least one, natural-product annotations only. All qualifying
features are written to fig5_rule_based_candidates.csv (Supplementary Table S1); the nine with the
lowest p-value (minimum across localities) are plotted as dot plots of biological-sample intensities."""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy import stats
from bio_stats import *
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from pubstyle import apply_style, WIDTH
apply_style()
import textwrap

raw = load()
LOC = {"CAM": ("NonhostCAM", "infectedCAM"), "ILO": ("nonhostILO", "infectedILO")}
out = pd.DataFrame({"Feature_Label": raw.Feature_Label, "Name": raw.Name, "RT": raw.RT})
mats = {}
for loc, (n, i) in LOC.items():
    N, I = bio_matrix(raw, n), bio_matrix(raw, i); mats[loc] = (N, I)
    p, fc, q = contrast(N, I)
    out[f"{loc}_log2FC"], out[f"{loc}_p"], out[f"{loc}_q"] = fc, p, q
    print(f"{loc}: nonhost n={N.shape[1]}, infected n={I.shape[1]}, family={np.isfinite(p).sum()}, "
          f"p<0.05={np.sum(p<0.05)}, q<0.05={np.sum(q<0.05)}, min q={np.nanmin(q):.3f}")
out.to_csv(BASE / "fig5_bio_all_features.csv", index=False)

from common import is_natural
out["min_p"] = out[["CAM_p", "ILO_p"]].min(axis=1)
rule = (out.CAM_log2FC > 1) & (out.ILO_log2FC > 1) & ((out.CAM_p < 0.05) | (out.ILO_p < 0.05))
cand = out[rule].copy(); cand["natural_product"] = cand.Name.map(is_natural)
cand.sort_values("min_p").to_csv(BASE / "fig5_rule_based_candidates.csv", index=False)
cand = cand[cand.natural_product].sort_values("min_p")
print(f"rule-based features: {rule.sum()} ({len(cand)} natural-product annotations); plotting top 9")
tab = cand.head(9).reset_index(drop=True)
idx = [int(out.index[out.Feature_Label == fl][0]) for fl in tab.Feature_Label]
tab.to_csv(BASE / "fig5_bio_displayed.csv", index=False)
print(tab[["Name","RT","CAM_log2FC","CAM_p","CAM_q","ILO_log2FC","ILO_p","ILO_q"]].round(4).to_string(index=False))

# dot plots (drawn at print width)
fig, axes = plt.subplots(3, 3, figsize=(WIDTH, 7.6)); axes = axes.ravel()
grp = [("CAM","infected"),("CAM","non-host"),("ILO","infected"),("ILO","non-host")]
col = {"infected":"#1f77b4","non-host":"#ff7f0e"}
rng = np.random.default_rng(0)
for k, (ax, j) in enumerate(zip(axes, idx)):
    r = tab.iloc[k]
    for x, (loc, g) in enumerate(grp):
        N, I = mats[loc]; v = (I if g=="infected" else N)[j]
        y = np.log10(v + 1)
        ax.scatter(x + rng.uniform(-0.08, 0.08, len(y)), y, s=12, c=col[g], edgecolors="k", lw=0.3, zorder=3)
        ax.hlines(y.mean(), x-0.25, x+0.25, color="k", lw=1)
    ax.set_xticks(range(4)); ax.set_xticklabels(["Inf.", "Non-\nhost", "Inf.", "Non-\nhost"], fontsize=7)
    ax.set_title(textwrap.fill(r.Name[0].upper() + r.Name[1:], 30) + f"\n[RT {r.RT:.2f} min]", fontsize=7)
    ymax = ax.get_ylim()[1]; ax.set_ylim(top=ymax * 1.38 + 0.5)
    for xc, loc in [(0.5, "CAM"), (2.5, "ILO")]:
        ax.text(xc, 0.97, f"{loc}\np = {r[f'{loc}_p']:.3f}\nq = {r[f'{loc}_q']:.2f}", transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=6.5)
    ax.axvline(1.5, c="#bbbbbb", lw=0.6)
    if k % 3 == 0: ax.set_ylabel("log$_{10}$(intensity + 1)", fontsize=7)
    ax.tick_params(labelsize=7)
fig.tight_layout(h_pad=1.0, w_pad=0.6)
fig.savefig(BASE / "nonhost_vs_infected_dotplots_bio.png", dpi=600); plt.close(fig)

# ---- Supplementary Table S1 (all natural-product features meeting the Fig. 5 rule)
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment
def clean_name(n):
    n = str(n).split(";")[0].strip()          # drop spectral-library metadata (e.g. "; LC-tDDA; CE30")
    return n.capitalize() if n.isupper() else n
info = raw[["Feature_Label", "m/z", "Ions"]].drop_duplicates("Feature_Label")
t1 = cand.merge(info, on="Feature_Label", how="left")
shown = set(tab.Feature_Label)
wb = Workbook(); ws = wb.active; ws.title = "Table S1"
ws["A1"] = ("Table S1. Features with higher mean abundance in non-host than infected Tetrastigma roots at both localities "
            "(log2FC > 1 in CAM and ILO) and p < 0.05 in at least one locality (Welch's t-test on biological-sample means; "
            "CAM n = 2 vs 2, ILO n = 4 vs 9). q, Benjamini-Hochberg FDR across all features tested per locality "
            f"({int(out.CAM_p.notna().sum()):,} CAM; {int(out.ILO_p.notna().sum()):,} ILO). Annotations are putative (MSI level 2-3). "
            "Rows are ordered by the lower of the two p-values; features shown in Fig. 5 are marked.")
ws["A1"].alignment = Alignment(wrap_text=True, vertical="top"); ws.merge_cells("A1:K1"); ws.row_dimensions[1].height = 75
ws.append(["Annotation", "RT (min)", "m/z", "Ion", "CAM log2FC", "CAM p", "CAM q", "ILO log2FC", "ILO p", "ILO q", "In Fig. 5"])
for _, r in t1.iterrows():
    ws.append([clean_name(r.Name), round(r.RT, 2), round(r["m/z"], 4), r.Ions, round(r.CAM_log2FC, 2), round(r.CAM_p, 4),
               round(r.CAM_q, 2), round(r.ILO_log2FC, 2), round(r.ILO_p, 4), round(r.ILO_q, 2), "yes" if r.Feature_Label in shown else ""])
excluded = sorted(out[rule & ~out.Name.map(is_natural)].Name)
ws.append([]); ws.append(["Features meeting the rule but excluded as not plausible plant natural products: " + "; ".join(excluded) + "."])
for row in ws.iter_rows():
    for cell in row: cell.font = Font(name="Arial", size=10, bold=(cell.row == 2))
ws.column_dimensions["A"].width = 55
for col in "BCDEFGHIJK": ws.column_dimensions[col].width = 11
wb.save(BASE / "Table_S1_nonhost_vs_infected.xlsx")
print("Table S1 rows:", len(t1))
