import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import textwrap

# ---------- 1) Load data ----------
file_path = "rafflesia_dataset4.xlsx"
df = pd.ExcelFile(file_path).parse(0)   # first sheet

# ---------- 2) Keep group-average columns + Feature_Label ----------
df.columns = df.columns.str.strip()
ave_cols = [c for c in df.columns if c.lower().startswith("ave")]
df = df[["Feature_Label"] + ave_cols].copy()

# ---------- 3) Replace negatives with 0 ----------
df[ave_cols] = df[ave_cols].clip(lower=0)

# ---------- 4) Extract compound names (after the '|') ----------
df["Compound"] = df["Feature_Label"].astype(str).str.split("|").str[-1].str.strip()

# Standardize selected compound names
COMPOUND_RENAMES = {
    "CITRIC ACID": "Citric acid",
    "Palmitic Acid": "Palmitic acid",
}

df["Compound"] = df["Compound"].replace(COMPOUND_RENAMES)

# Exclude manually flagged compounds before ranking the top 20
EXCLUDED_COMPOUNDS = {
    # (previously excluded "12-Hydroxyoctadecanoic acid" here, but that
    # compound is the #1 most abundant feature in raffseed and is cited by
    # name in the manuscript Discussion as evidence of seed fatty-acid/
    # oxylipin enrichment -- excluding it directly contradicted the text.)
}

df = df[~df["Compound"].isin(EXCLUDED_COMPOUNDS)].copy()

# ---------- 5) Collapse duplicates by averaging intensities ----------
# (mean across rows for each compound name)
compound_avg = df.groupby("Compound", as_index=True)[ave_cols].mean()

# ---------- 6) Rank compounds and keep top 20 ----------
# Default: rank by overall abundance across all groups
overall = compound_avg.sum(axis=1)
top20_names = overall.nlargest(20).index

# If you instead want to rank by RAFFSEED only, uncomment:
# raff_col = [c for c in ave_cols if "raffseed" in c.lower()][0]
# top20_names = compound_avg[raff_col].nlargest(20).index

top20 = compound_avg.loc[top20_names]

# ---------- 7) Total Sum Scaling (normalize each group column to sum=1) ----------
tss = top20.div(top20.sum(axis=0), axis=1)

# (Optional) choose a specific order for the x-axis:
# order = ["averaffseed","aveRaffbudspec","aveinfectedILO","aveUNinfectedILO",
#          "aveuninfecraffspec-stemleaf","avenonhostILO"]
# cols_to_plot = [c for c in order if c in tss.columns]
# Otherwise, keep whatever columns exist in the file:
cols_to_plot = list(tss.columns)

# ---------- 8) Plot stacked bar ----------
fig, ax = plt.subplots(figsize=(15, 7))
tss[cols_to_plot].T.plot(kind="bar", stacked=True, ax=ax, width=0.85, colormap="tab20")

ax.set_ylabel("Normalized Abundance (Total Sum Scaling)")
ax.set_xlabel("")
handles, labels = ax.get_legend_handles_labels()
wrapped_labels = [textwrap.fill(label, width=42) for label in labels]

ax.legend(
    handles,
    wrapped_labels,
    title="Compound Name",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    fontsize=8,
    frameon=False
)

fig.subplots_adjust(
    left=0.09,
    right=0.66,
    bottom=0.30,
    top=0.92
)

plt.setp(
    ax.get_xticklabels(),
    rotation=45,
    ha="right",
    rotation_mode="anchor"
)

# Save figure (PNG/SVG)
fig.savefig(
    "top20_compounds_stacked_TSS_FIXED.png",
    dpi=300,
    bbox_inches="tight",
    pad_inches=0.2
)

fig.savefig(
    "top20_compounds_stacked_TSS_FIXED.svg",
    bbox_inches="tight",
    pad_inches=0.2
)

# ---------- 9) Export the numbers used in the figure ----------
# Raw means (averaged duplicates) for the top 20:
top20.to_csv("top20_compounds_raw_means_FIXED.csv")
# TSS-normalized values used for plotting:
tss.to_csv("top20_compounds_TSS_FIXED.csv")

print("Saved: top20_compounds_stacked_TSS_FIXED.png/.svg, top20_compounds_raw_means_FIXED.csv, top20_compounds_TSS_FIXED.csv")
