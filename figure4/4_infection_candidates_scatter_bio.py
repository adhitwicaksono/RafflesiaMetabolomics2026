"""Fig. 4 (revised): pooled infected vs apparently uninfected, biological-sample level.
Features sharing an annotation are consolidated by name (case-insensitive; '.mol' suffix removed).
Labels: every metabolite meeting p < 0.05 and log2FC > 1 at the biological-sample level."""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from bio_stats import *
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from pubstyle import apply_style, label as glabel, WIDTH
apply_style()
from common import load_raw
raw = load_raw()
num = [c for c in raw.columns if c not in ("RT","m/z","M meas","Ions","MS/MS","Name","Molecular_Formula","Feature_Label","name_key")]
d = raw.copy(); d[num] = d[num].apply(pd.to_numeric, errors="coerce").clip(lower=0)
d = d.groupby("Name", as_index=False)[num].mean()
IG = ["infectedTHAI","infectedCAM","infectedILO"]; UG = ["UNinfectedTHAI","UNinfectedCAM","UNinfectedILO"]
def bm(g):
    X = d[run_cols(d, g)].values; return (X[:,0::2] + X[:,1::2]) / 2
A = np.hstack([bm(g) for g in IG]); B = np.hstack([bm(g) for g in UG])
p, fc, q = contrast(A, B)
res = pd.DataFrame({"Name": d.Name, "log2FC": fc, "p": p, "q_BH": q}).dropna(subset=["p"])
res.sort_values("p").to_csv(BASE/"fig4_bio_consolidated.csv", index=False)
# Labelled metabolites: all metabolites meeting the nominal thresholds (p < 0.05, log2FC > 1)
# at the biological-sample level, each assigned a chemical class.
CLS = {"Glycoside":"#1f77b4", "Phytohormone/signaling":"#2ca02c", "Flavonoid":"#9467bd"}
CLASS_OF = {"Nb-trans-Feruloylserotonin glucoside":"Glycoside",
            "Phloroacetophenone 6'-[xylosyl-(1->6)-glucoside]":"Glycoside",
            "Gibberellin A39":"Phytohormone/signaling", "3,5-Digalloylepicatechin":"Flavonoid"}
hits = res[(res.p < 0.05) & (res.log2FC > 1)]
missing = set(hits.Name) - set(CLASS_OF)
assert not missing, f"unclassified nominal hits: {missing}"
c = hits.set_index("Name")
print("family", len(res), "min q %.2f" % res.q_BH.min(), "p<0.05:", (res.p<0.05).sum(), "nominal hits:", len(hits)); print(c.round(4))
plot = res[res.log2FC >= -0.5]
fig, ax = plt.subplots(figsize=(WIDTH, 4.6))
ax.scatter(plot.log2FC, -np.log10(plot.p), s=6, c="#d0d0d0", edgecolors="none", alpha=0.7, zorder=1)
off = {"Nb-trans-Feruloylserotonin glucoside":(-128,-26), "Phloroacetophenone 6'-[xylosyl-(1->6)-glucoside]":(-40,8),
       "Gibberellin A39":(-82,4), "3,5-Digalloylepicatechin":(8,2)}
names = {"Phloroacetophenone 6'-[xylosyl-(1->6)-glucoside]":"Phloroacetophenone 6'-[xylosyl-(1→6)-glucoside]"}
for n, r in c.iterrows():
    y = -np.log10(r.p); col = CLS[CLASS_OF[n]]
    ax.scatter(r.log2FC, y, s=40, c=col, edgecolors="#404040", lw=0.6, zorder=3)
    ax.annotate(f"{names.get(n, n)}\n(p={r.p:.3f}, q={r.q_BH:.2f})", (r.log2FC, y), xytext=off.get(n, (6, 4)), textcoords="offset points", fontsize=7, color=col)
ax.axhline(-np.log10(0.05), ls="--", c="#303030", lw=1); ax.axvline(1, ls="--", c="#303030", lw=1)
ax.set_xlabel("log₂FC (infected / uninfected)"); ax.set_ylabel("−log₁₀(p-value)")
ax.set_ylim(top=2.05)
from matplotlib.lines import Line2D
h = [Line2D([0],[0],marker="o",color="none",markerfacecolor="#d0d0d0",markersize=5,label="Other metabolites")]
for cls,col in CLS.items(): h.append(Line2D([0],[0],marker="o",color="none",markerfacecolor=col,markeredgecolor="#404040",markersize=6,label=cls))
ax.legend(handles=h, frameon=False, loc="upper right", fontsize=7, title="Chemical class", title_fontsize=7.5)
fig.text(0.01, 0.005, f"Welch's t-test on biological-sample means (infected n = 13, uninfected n = 14);\nBH FDR across {len(res)} metabolites: none has q < 0.05. log₂FC < −0.5 not shown.",
        fontsize=7, style="italic", va="bottom")
fig.tight_layout(rect=(0,0.06,1,1)); fig.savefig(BASE/"infection_candidates_scatter_bio.png", dpi=600); plt.close(fig)
