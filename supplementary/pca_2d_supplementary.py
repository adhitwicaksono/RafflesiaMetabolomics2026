"""Supplementary PC1–PC3 / PC2–PC3 panels with scree
plots (Figs. S1–S3), and Fig. 1 recomputed without Ampelopsis and R. speciosa seeds (Fig. S4).
Same preprocessing as the submitted scripts: exact duplicate rows dropped, missing -> 0,
autoscaling, PCA (random_state=0). Points are technical runs (two per biological sample)."""
from pathlib import Path
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.stats import chi2
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from pubstyle import apply_style, label as glabel, WIDTH
apply_style()
BASE = Path(__file__).resolve().parent
raw = pd.read_excel(BASE.parent / "raw_data" / "rafflesia_dataset1.xlsx").dropna(how="all").drop_duplicates().set_index("Feature_Label")
raw = raw.loc[:, ~raw.columns.str.contains("ave", case=False)]
GROUPS = ["sapbud","infectedTHAI","UNinfectedTHAI","Raffbudlag","infectedCAM","UNinfectedCAM","NonhostCAM","raffseed",
          "Raffbudspec","infectedILO","UNinfectedILO","uninfecraffspec-stemleaf","nonhostILO","Ampelopsis"]
def grp(c):
    n = c.lower()
    for g in ["UNinfectedTHAI","UNinfectedCAM","UNinfectedILO"] + [g for g in GROUPS if not g.startswith("UNinfected")]:
        if g.lower() in n: return g
    return None
col2g = {c: grp(c) for c in raw.columns if grp(c)}
STYLE = {"sapbud":("#6baed6","o"),"Raffbudspec":("#4292c6","o"),"Raffbudlag":("#2171b5","o"),
  "infectedTHAI":("#31a354","s"),"infectedCAM":("#74c476","s"),"infectedILO":("#006d2c","s"),
  "UNinfectedTHAI":("#fd8d3c","D"),"UNinfectedCAM":("#fdae6b","D"),"UNinfectedILO":("#fdd0a2","D"),
  "uninfecraffspec-stemleaf":("#9e9ac8","1"),"nonhostILO":("#a05a2c","X"),"NonhostCAM":("#8c6d31","X"),
  "raffseed":("#252525","^"),"Ampelopsis":("#f768a1","v"),
  "BUD":("#3182bd","o"),"INFECTED":("#31a354","s"),"UNINFECTED":("#fd8d3c","D"),"nonhostTET":("#8c6d31","X"),
  "UNINFRAFFSPEC":("#9e9ac8","1"),"RAFFSEED":("#252525","^")}
SUPER = {"BUD":["sapbud","Raffbudlag","Raffbudspec"],"INFECTED":["infectedTHAI","infectedCAM","infectedILO"],
  "UNINFECTED":["UNinfectedTHAI","UNinfectedCAM","UNinfectedILO"],"nonhostTET":["NonhostCAM","nonhostILO"],
  "UNINFRAFFSPEC":["uninfecraffspec-stemleaf"],"RAFFSEED":["raffseed"]}

def run_pca(exclude=(), keep=None, n=10):
    cols = [c for c, g in col2g.items() if g not in exclude and (keep is None or g in keep)]
    X = raw[cols].T.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    p = PCA(n_components=n, random_state=0); S = p.fit_transform(StandardScaler().fit_transform(X.to_numpy()))
    return S, np.array([col2g[c] for c in cols]), p.explained_variance_ratio_ * 100

def ellipse(ax, P, color, cov=0.68):
    if len(P) < 3: return
    C = np.cov(P, rowvar=False); w, v = np.linalg.eigh(C); r = np.sqrt(chi2.ppf(cov, 2) * w)
    ang = np.degrees(np.arctan2(v[1, 1], v[0, 1]))
    ax.add_patch(Ellipse(P.mean(0), 2*r[1], 2*r[0], angle=ang, color=color, alpha=0.18, lw=0, zorder=0))

def scatter(ax, S, labels, order, a, b, ve, ells=None, legend=True):
    for g in order:
        m = labels == g
        if not m.any(): continue
        c, mk = STYLE[g]
        ax.scatter(S[m, a], S[m, b], s=12 if mk != "1" else 22, c=c, marker=mk, edgecolors="none" if mk != "1" else None, lw=1.0 if mk == "1" else None, label=glabel(g), zorder=2)
    for members, col in (ells or []):
        ellipse(ax, S[np.isin(labels, members)][:, [a, b]], col)
    ax.set_xlabel(f"PC{a+1} ({ve[a]:.1f}%)"); ax.set_ylabel(f"PC{b+1} ({ve[b]:.1f}%)")
    ax.axhline(0, c="#999", lw=0.5, zorder=0); ax.axvline(0, c="#999", lw=0.5, zorder=0)
    if legend: ax.legend(title="Group", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)

def scree(ax, ve):
    k = np.arange(1, len(ve)+1)
    ax.bar(k, ve, color="#9ecae1", edgecolor="#3182bd"); ax.plot(k, np.cumsum(ve), "o-", c="#08519c", ms=3)
    ax.set_xticks(k); ax.set_xlabel("Principal component"); ax.set_ylabel("Variance explained (%)")
    ax.legend(["Cumulative", "Individual"], frameon=False)

ORDER1 = ["sapbud","Raffbudspec","Raffbudlag","infectedTHAI","infectedCAM","infectedILO","UNinfectedTHAI","UNinfectedCAM",
          "UNinfectedILO","uninfecraffspec-stemleaf","nonhostILO","NonhostCAM","raffseed","Ampelopsis"]
ELL2 = [(["infectedTHAI","infectedCAM","infectedILO"], "#006d2c"), (["UNinfectedTHAI","UNinfectedCAM","UNinfectedILO"], "#fd8d3c")]
sets = {
  "1": dict(kw={}, order=ORDER1, ells=None),
  "2": dict(kw={"exclude": {"sapbud","Raffbudspec","Raffbudlag","raffseed"}}, order=[g for g in ORDER1 if g not in {"sapbud","Raffbudspec","Raffbudlag","raffseed"}], ells=ELL2),
  "S4": dict(kw={"exclude": {"Ampelopsis","raffseed"}}, order=[g for g in ORDER1 if g not in {"Ampelopsis","raffseed"}], ells=None),
}
summary = {}
for key, d in sets.items():
    S, lab, ve = run_pca(**d["kw"])
    if key == "2": S[:, 2] *= -1  # same PC3 orientation as submitted Fig. 2
    summary[key] = ve[:3]
    if key != "S4":
        fig, axs = plt.subplots(2, 2, figsize=(WIDTH, 5.8)); axs = axs.ravel()
        scatter(axs[0], S, lab, d["order"], 0, 2, ve, d["ells"], legend=False); scatter(axs[1], S, lab, d["order"], 1, 2, ve, d["ells"], legend=False)
        scree(axs[2], ve); h, l = axs[0].get_legend_handles_labels()
        axs[3].axis("off"); axs[3].legend(h, l, title="Group", loc="center", frameon=False); fig.tight_layout()
        fig.savefig(BASE / f"FigS{key}_PC13_PC23_scree.png", dpi=600); plt.close(fig)
    else:
        fig, axs = plt.subplots(2, 2, figsize=(WIDTH, 5.8)); axs = axs.ravel()
        scatter(axs[0], S, lab, d["order"], 0, 1, ve, None, legend=False); scatter(axs[1], S, lab, d["order"], 0, 2, ve, None, legend=False); scree(axs[2], ve)
        h, l = axs[0].get_legend_handles_labels(); axs[3].axis("off"); axs[3].legend(h, l, title="Group", loc="center", frameon=False)
        fig.tight_layout(); fig.savefig(BASE / "FigS4_no_Ampelopsis_seeds.png", dpi=600); plt.close(fig)

# Fig. 3 supergroups (Ampelopsis not included, as in submitted Fig. 3)
cols = [c for c, g in col2g.items() if any(g in m for m in SUPER.values())]
lab = np.array([next(k for k, m in SUPER.items() if col2g[c] in m) for c in cols])
X = raw[cols].T.apply(pd.to_numeric, errors="coerce").fillna(0.0)
p = PCA(n_components=10, random_state=0); S = p.fit_transform(StandardScaler().fit_transform(X.to_numpy())); ve = p.explained_variance_ratio_ * 100
summary["3"] = ve[:3]; o3 = ["BUD","INFECTED","UNINFECTED","nonhostTET","UNINFRAFFSPEC","RAFFSEED"]
fig, axs = plt.subplots(2, 2, figsize=(WIDTH, 5.8)); axs = axs.ravel(); scatter(axs[0], S, lab, o3, 0, 2, ve, legend=False); scatter(axs[1], S, lab, o3, 1, 2, ve, legend=False); scree(axs[2], ve)
h, l = axs[0].get_legend_handles_labels(); axs[3].axis("off"); axs[3].legend(h, l, title="Group", loc="center", frameon=False); fig.tight_layout()
fig.savefig(BASE / "FigS3_PC13_PC23_scree.png", dpi=600); plt.close(fig)
for k, v in summary.items(): print(k, np.round(v, 1), "cum3 %.1f" % v.sum(), "PC1+2 %.1f" % v[:2].sum())
