import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa
from scipy.stats import chi2

# For a 3D multivariate-normal covariance ellipsoid, the fraction of points
# enclosed at radius n_std (in SD units) follows a chi-squared distribution
# with 3 degrees of freedom: coverage = chi2.cdf(n_std**2, df=3).
# n_std=1.0 only encloses ~19.9% of points in 3D (vs. 68% in 1D).
# Solve for the n_std that gives true 95% coverage:
TARGET_COVERAGE = 0.68
N_STD = chi2.ppf(TARGET_COVERAGE, df=3) ** 0.5  # ~1.87

# --- 1) Load sheet ---
file_path = r'..\raw_data\rafflesia_dataset1.xlsx'
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
    ordered = ["UNinfectedTHAI", "UNinfectedCAM", "UNinfectedILO"] + \
              [g for g in SUPERGROUPS if not g.startswith("UNinfected")]
    for g in ordered:
        if g.lower() in name:
            return g
    return "Unknown"

col2grp = {c: to_supergroup(c) for c in df.columns}

# Groups excluded from this PCA
EXCLUDED_GROUPS = {
    "sapbud",
    "Raffbudspec",
    "Raffbudlag",
    "raffseed"
}

valid_cols = [
    c for c, g in col2grp.items()
    if g != "Unknown" and g not in EXCLUDED_GROUPS
]

df = df[valid_cols]
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
    "sapbud": ("#6baed6", "o"),
    "Raffbudspec": ("#4292c6", "o"),
    "Raffbudlag": ("#2171b5", "o"),
    "infectedTHAI": ("#31a354", "s"),
    "infectedCAM": ("#74c476", "s"),
    "infectedILO": ("#006d2c", "s"),
    "UNinfectedTHAI": ("#fd8d3c", "D"),
    "UNinfectedCAM": ("#fdae6b", "D"),
    "UNinfectedILO": ("#fdd0a2", "D"),
    "uninfecraffspec-stemleaf": ("#9e9ac8", "1"),
    "nonhostILO": ("#a05a2c", "X"),
    "NonhostCAM": ("#8c6d31", "X"),
    "raffseed": ("#252525", "^"),
    "Ampelopsis": ("#f768a1", "v"),
}

legend_order = [
    "infectedTHAI",
    "infectedCAM",
    "infectedILO",
    "UNinfectedTHAI",
    "UNinfectedCAM",
    "UNinfectedILO",
    "uninfecraffspec-stemleaf",
    "nonhostILO",
    "NonhostCAM",
    "Ampelopsis"
]

fig = plt.figure(figsize=(11, 8))
ax = fig.add_subplot(111, projection="3d")

def plot_group_ellipsoid(
    ax,
    points,
    color,
    n_std=N_STD,
    alpha=0.15,
    resolution=40
):
    """
    Draw a transparent 3D covariance ellipsoid around one PCA group.

    Parameters
    ----------
    ax : matplotlib 3D axis
        Axis on which the ellipsoid is drawn.
    points : array-like, shape (n_samples, 3)
        PC1, PC2, and PC3 coordinates for one group.
    color : str
        Matplotlib-compatible color.
    n_std : float
        Ellipsoid radius in standard-deviation units.
    alpha : float
        Transparency of the ellipsoid.
    resolution : int
        Surface smoothness.
    """

    points = np.asarray(points, dtype=float)

    # Fewer than three observations cannot define a useful 3D covariance shape
    if points.shape[0] < 3:
        return

    center = points.mean(axis=0)
    covariance = np.cov(points, rowvar=False)

    # Ellipsoid directions and axis lengths
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)

    # Prevent tiny negative eigenvalues caused by numerical rounding
    eigenvalues = np.clip(eigenvalues, 0, None)
    radii = n_std * np.sqrt(eigenvalues)

    # Parametric unit sphere
    u = np.linspace(0, 2 * np.pi, resolution)
    v = np.linspace(0, np.pi, resolution)

    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones_like(u), np.cos(v))

    sphere = np.stack((x, y, z), axis=-1)

    # Scale, rotate, and move the sphere to the group center
    ellipsoid = sphere @ np.diag(radii) @ eigenvectors.T
    ellipsoid += center

    ax.plot_surface(
        ellipsoid[:, :, 0],
        ellipsoid[:, :, 1],
        ellipsoid[:, :, 2],
        color=color,
        alpha=alpha,
        edgecolor="none",
        linewidth=0,
        shade=False
    )

for g in legend_order:
    idx = [i for i, gg in enumerate(groups) if gg == g]

    if not idx:
        continue

    col, mk = style[g]

    ax.scatter(
        scores[idx, 0],
        scores[idx, 1],
        scores[idx, 2],
        s=70,
        c=col,
        marker=mk,
        edgecolors="none",
        alpha=1.0,
        depthshade=False,
        label=g
    )

# --- Pooled supergroup ellipsoids ---
# Draw ONE ellipsoid per biological supergroup (infected / uninfected / non-host),
# pooling across localities (THAI, CAM, ILO), rather than one ellipsoid per
# locality subgroup. Locality-level markers/colors above are unchanged.
ELLIPSOID_SUPERGROUPS = {
    "infected":   (["infectedTHAI", "infectedCAM", "infectedILO"],   "#006d2c"),  # infectedILO dark green
    "uninfected": (["UNinfectedTHAI", "UNinfectedCAM", "UNinfectedILO"], "#fd8d3c"),  # mid orange
}

for label, (members, ell_color) in ELLIPSOID_SUPERGROUPS.items():
    idx = [i for i, gg in enumerate(groups) if gg in members]
    if not idx:
        continue
    plot_group_ellipsoid(
        ax,
        scores[idx, :3],
        color=ell_color,
        n_std=N_STD,
        alpha=0.15
    )

# uninfecraffspec-stemleaf has no locality split, so it can still get its own
# ellipsoid if it has enough points; not part of ELLIPSOID_SUPERGROUPS above.
idx_stem = [i for i, gg in enumerate(groups) if gg == "uninfecraffspec-stemleaf"]
if idx_stem:
    col, _ = style["uninfecraffspec-stemleaf"]
    plot_group_ellipsoid(
        ax,
        scores[idx_stem, :3],
        color=col,
        n_std=N_STD,
        alpha=0.15
        )

ax.set_xlabel(f"PC1 ({ve[0]:.1f}%)")
ax.set_ylabel(f"PC2 ({ve[1]:.1f}%)")
ax.set_zlabel(f"PC3 ({ve[2]:.1f}%)")
# rotate to a reference-like view (tweak as you wish)
ax.view_init(elev=22, azim=-35)
ax.grid(True, linestyle="--", alpha=0.4)

# Equal aspect ratio: scale each axis box by its true data range so that
# distances (and the covariance ellipsoids) are rendered geometrically
# correctly rather than independently stretched to fill the box.
ranges = np.ptp(scores[:, :3], axis=0)  # PC1, PC2, PC3 point-to-point range
ax.set_box_aspect(tuple(ranges))

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
