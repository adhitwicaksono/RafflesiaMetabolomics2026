import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa

# --- 1) Load data ---
file_path = "rafflesia_dataset1.xlsx"
df = pd.read_excel(file_path, sheet_name="Sheet1")

# --- 2) Feature IDs and sample columns ---
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

# --- 3) Remove average columns ---
df = df.loc[:, ~df.columns.str.contains("ave", case=False)]

# --- 4) Supergroup mapping ---
supergroups = {
    "BUD": ["sapbud", "Raffbudlag", "Raffbudspec"],
    "UNINFECTED": ["UNinfectedTHAI", "UNinfectedCAM", "UNinfectedILO"],
    "INFECTED": ["infectedTHAI", "infectedCAM", "infectedILO"],
    "RAFFSEED": ["raffseed"],
    "UNINFRAFFSPEC": ["uninfecraffspec-stemleaf"],
    "nonhostTET": ["NonhostCAM", "nonhostILO"],
}

def assign_supergroup(col_name: str) -> str:
    low = col_name.lower()
    for sg, members in supergroups.items():
        for m in members:
            if m.lower() in low:
                return sg
    return "Unknown"

col2grp = {c: assign_supergroup(c) for c in df.columns}
valid_cols = [c for c, g in col2grp.items() if g != "Unknown"]
df = df[valid_cols]
sample_groups = pd.Series([col2grp[c] for c in valid_cols], index=valid_cols, name="Group")

# --- 5) Transpose: rows = samples ---
X = df.T  # samples x features

# --- 6) Fill missing with 0 ---
X = X.apply(pd.to_numeric, errors="coerce").fillna(0.0)

# --- 7) Autoscaling (z-score) ---
X_scaled = StandardScaler().fit_transform(X)

# --- 8) PCA (3 components) ---
pca = PCA(n_components=3, random_state=0)
scores = pca.fit_transform(X_scaled)
ve = pca.explained_variance_ratio_ * 100

# --- 9) Plot (3D) ---
style = {
    "BUD": ("#3182bd", "o"),
    "INFECTED": ("#31a354", "s"),
    "RAFFSEED": ("#252525", "^"),
    "UNINFECTED": ("#fd8d3c", "D"),
    "UNINFRAFFSPEC": ("#9e9ac8", "P"),
    "nonhostTET": ("#8c6d31", "X"),
}
order = ["BUD", "INFECTED", "UNINFECTED", "nonhostTET", "UNINFRAFFSPEC", "RAFFSEED"]

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection="3d")

scores_df = pd.DataFrame(scores, index=X.index, columns=["PC1", "PC2", "PC3"])
scores_df["Group"] = sample_groups.values

for g in order:
    sub = scores_df[scores_df["Group"] == g]
    if sub.empty:
        continue
    color, marker = style[g]
    ax.scatter(
        sub["PC1"], sub["PC2"], sub["PC3"],
        c=color, marker=marker, s=70,
        edgecolor="none", alpha=1.0, depthshade=False, label=g
    )

ax.set_xlabel(f"PC1 ({ve[0]:.1f}%)")
ax.set_ylabel(f"PC2 ({ve[1]:.1f}%)")
ax.set_zlabel(f"PC3 ({ve[2]:.1f}%)")
ax.view_init(elev=22, azim=-35)
ax.legend(title="Supergroup", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
plt.tight_layout()

# --- Export scores and loadings before opening the interactive plot window ---
scores_out = "PCA_scores_supergroups.csv"
loadings_out = "PCA_loadings_supergroups.csv"
scores_df.to_csv(scores_out)
pd.DataFrame(pca.components_.T, index=X.columns, columns=["PC1", "PC2", "PC3"]).to_csv(loadings_out)

# --- Save figure files ---
fig_out_png = "PCA_3D_supergroups.png"
fig_out_svg = "PCA_3D_supergroups.svg"
fig.savefig(fig_out_png, dpi=300, bbox_inches="tight")
fig.savefig(fig_out_svg, bbox_inches="tight")

print("Saved:", scores_out, loadings_out, fig_out_png, fig_out_svg)

plt.close(fig)
