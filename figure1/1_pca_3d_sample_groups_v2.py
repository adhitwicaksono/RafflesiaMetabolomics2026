import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa

# --- 1) Load sheet ---
file_path = "rafflesia_dataset1.xlsx"
df = pd.read_excel(file_path, sheet_name="Sheet1")

# --- 2) Keep only sample-intensity columns ---
# Use Feature_Label as feature ID when available
df = df.dropna(how="all")

# Drop exact duplicate rows (data-entry copy/paste duplicates in source spreadsheet)
n_before = len(df)
df = df.drop_duplicates()
n_after = len(df)
print(f"Dropped {n_before - n_after} exact-duplicate rows ({n_before} -> {n_after})")

if "Feature_Label" in df.columns:
    df = df.set_index("Feature_Label")
else:
    df = df.set_index(df.columns[0])

# drop columns containing 'ave'
df = df.loc[:, ~df.columns.str.contains("ave", case=False)]

# --- 3) Supergroup assignment ---
SUPERGROUPS = [
    "sapbud", "infectedTHAI", "UNinfectedTHAI", "Raffbudlag",
    "infectedCAM", "UNinfectedCAM", "NonhostCAM", "raffseed",
    "Raffbudspec", "infectedILO", "UNinfectedILO",
    "uninfecraffspec-stemleaf", "nonhostILO", "Ampelopsis"
]

def to_supergroup(col_name: str) -> str:
    name = col_name.lower()
    # ensure “UNinfected” doesn’t fall into “infected”
    ordered = ["UNinfectedTHAI","UNinfectedCAM","UNinfectedILO"] + \
              [g for g in SUPERGROUPS if not g.startswith("UNinfected")]
    for g in ordered:
        if g.lower() in name:
            return g
    return "Unknown"

col2grp = {c: to_supergroup(c) for c in df.columns}
valid_cols = [c for c,g in col2grp.items() if g != "Unknown"]
df = df[valid_cols]                      # keep only mapped samples
groups = [col2grp[c] for c in valid_cols]

# --- 4) Matrix orientation: rows = samples, cols = features ---
X = df.T   # shape: (n_samples, n_features)

# Ensure all intensity values are numeric
X = X.apply(pd.to_numeric, errors="coerce").fillna(0.0)

# --- 5) Autoscaling (mean-centered, unit variance) ---
X_scaled = StandardScaler().fit_transform(X)

# --- 6) PCA ---
pca = PCA(n_components=3, random_state=0)
scores = pca.fit_transform(X_scaled)
ve = pca.explained_variance_ratio_ * 100

# --- 7) Plot (publication-clean) ---
# marker/color style similar to your reference
style = {
    "sapbud": ("#6baed6","o"),
    "Raffbudspec": ("#4292c6","o"),
    "Raffbudlag": ("#2171b5","o"),
    "infectedTHAI": ("#31a354","s"),
    "infectedCAM": ("#74c476","s"),
    "infectedILO": ("#006d2c","s"),
    "UNinfectedTHAI": ("#fd8d3c","D"),
    "UNinfectedCAM": ("#fdae6b","D"),
    "UNinfectedILO": ("#fdd0a2","D"),
    "uninfecraffspec-stemleaf": ("#9e9ac8","1"),
    "nonhostILO": ("#a05a2c","X"),
    "NonhostCAM": ("#8c6d31","X"),
    "raffseed": ("#252525","^"),
    "Ampelopsis": ("#f768a1","v"),
}

legend_order = [
    "sapbud","Raffbudspec","Raffbudlag",
    "infectedTHAI","infectedCAM","infectedILO",
    "UNinfectedTHAI","UNinfectedCAM","UNinfectedILO",
    "uninfecraffspec-stemleaf","nonhostILO","NonhostCAM",
    "raffseed","Ampelopsis"
]

fig = plt.figure(figsize=(11, 8))
ax = fig.add_subplot(111, projection="3d")

for g in legend_order:
    idx = [i for i, gg in enumerate(groups) if gg == g]
    if not idx: 
        continue
    col, mk = style[g]
    ax.scatter(scores[idx,0], scores[idx,1], scores[idx,2],
               s=70, c=col, marker=mk, edgecolor="none", alpha=1.0, depthshade=False, label=g)

ax.set_xlabel(f"PC1 ({ve[0]:.1f}%)")
ax.set_ylabel(f"PC2 ({ve[1]:.1f}%)")
ax.set_zlabel(f"PC3 ({ve[2]:.1f}%)")
# rotate to a reference-like view (tweak as you wish)
ax.view_init(elev=22, azim=-35)
ax.grid(True, linestyle="--", alpha=0.4)

# legend in the same order as legend_order
handles, labels = ax.get_legend_handles_labels()
ordered = [(h,l) for g in legend_order for h,l in zip(handles, labels) if l == g]
if ordered:
    h_ord, l_ord = zip(*ordered)
    ax.legend(h_ord, l_ord, title="Group", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
else:
    ax.legend(title="Group", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)

plt.tight_layout()

# --- 8) Export PCA scores & loadings ---
scores_out = "PCA_scores_autoscaled_main-data.csv"
loadings_out = "PCA_loadings_autoscaled_main-data.csv"

scores_df = pd.DataFrame(scores, index=X.index, columns=["PC1", "PC2", "PC3"])
scores_df["Group"] = groups
loadings_df = pd.DataFrame(pca.components_.T, index=df.index, columns=["PC1", "PC2", "PC3"])

scores_df.to_csv(scores_out)
loadings_df.to_csv(loadings_out)

# --- 9) Save figure files ---
figure_png = "PCA_3D_sample_groups.png"
figure_svg = "PCA_3D_sample_groups.svg"

fig.savefig(figure_png, dpi=300, bbox_inches="tight")
fig.savefig(figure_svg, bbox_inches="tight")

print("Saved:", scores_out, loadings_out, figure_png, figure_svg)

plt.close(fig)
