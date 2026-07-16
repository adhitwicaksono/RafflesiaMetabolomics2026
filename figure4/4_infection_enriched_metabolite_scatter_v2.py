#!/usr/bin/env python3
"""
Generate a one-sided significance scatterplot for metabolites enriched in
infected versus uninfected tissues.

This is intentionally not labelled as a conventional volcano plot because
only the infection-enriched side is emphasized.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# ---------- Configuration ----------
BASE_DIR = Path(__file__).resolve().parent
INPUT_FILE = BASE_DIR.parent / "raw_data" / "rafflesia_dataset2.xlsx"
SHEET_NAME = "Sheet1"

MIN_DISPLAY_LOG2FC = -0.5
SIGNIFICANT_LOG2FC = 1.0
P_VALUE_THRESHOLD = 0.05
FDR_THRESHOLD = 0.10  # relaxed FDR cutoff, shown for transparency only

OUTPUT_PNG = BASE_DIR / "infection_enriched_metabolite_scatter.png"
OUTPUT_SVG = BASE_DIR / "infection_enriched_metabolite_scatter.svg"
OUTPUT_CSV = BASE_DIR / "infection_enriched_significant_metabolites.csv"

CLASS_COLORS = {
    "Glycoside": "#1f77b4",
    "Phytohormone/signaling": "#2ca02c",
    "Phenolic acid/phenylpropanoid": "#d62728",
    "Flavonoid": "#9467bd",
}
OTHER_COLOR = "#d9d9d9"

# Manual label offsets, in display points, to reproduce the intended layout.
LABEL_OFFSETS = {
    "Nb-trans-Feruloylserotonin glucoside": (8, 8),
    "Gibberellin A93": (-60, 18),
    "ellagic acid": (-68, 10),
    "4-Methoxycinnamic acid": (-105, 18),
    "quercetin-3-O-deoxyhexosyl(1-2)pentoside": (-150, 18),
}


def main() -> None:
    if not INPUT_FILE.exists():
        raise FileNotFoundError(
            f"Input file not found: {INPUT_FILE.name}\n"
            f"Place the Excel file in the directory {BASE_DIR.parent / 'raw_data'}."
        )

    df = pd.read_excel(INPUT_FILE, sheet_name=SHEET_NAME)

    required_columns = {
        "Metabolite",
        "log2FC_pooled",
        "p_value_Welch",
        "FDR_BH",
        "Class",
    }
    missing = sorted(required_columns.difference(df.columns))
    if missing:
        raise ValueError(
            "Missing required Excel column(s): " + ", ".join(missing)
        )

    # Ensure numeric values and remove rows that cannot be plotted.
    df["log2FC_pooled"] = pd.to_numeric(df["log2FC_pooled"], errors="coerce")
    df["p_value_Welch"] = pd.to_numeric(df["p_value_Welch"], errors="coerce")
    df["FDR_BH"] = pd.to_numeric(df["FDR_BH"], errors="coerce")
    df = df.dropna(subset=["log2FC_pooled", "p_value_Welch"]).copy()

    # Prevent log10(0) while preserving extremely small p-values.
    smallest_positive = np.nextafter(0, 1)
    df["neglog10_p"] = -np.log10(df["p_value_Welch"].clip(lower=smallest_positive))

    # Keep the same displayed x-range logic as the reference figure.
    plotted = df[df["log2FC_pooled"] >= MIN_DISPLAY_LOG2FC].copy()

    significant = plotted[
        (plotted["log2FC_pooled"] >= SIGNIFICANT_LOG2FC)
        & (plotted["p_value_Welch"] < P_VALUE_THRESHOLD)
    ].copy()

    # Only these four classes are highlighted in the reference figure.
    highlighted = significant[significant["Class"].isin(CLASS_COLORS)].copy()

    fig, ax = plt.subplots(figsize=(11, 8))

    # Background points: all displayed metabolites, including significant
    # metabolites whose classes are not among the four highlighted categories.
    ax.scatter(
        plotted["log2FC_pooled"],
        plotted["neglog10_p"],
        s=20,
        c=OTHER_COLOR,
        alpha=0.28,
        edgecolors="none",
        zorder=1,
    )

    # Highlight selected chemical classes.
    for chemical_class, color in CLASS_COLORS.items():
        subset = highlighted[highlighted["Class"] == chemical_class]
        if subset.empty:
            continue

        fdr_pass = subset["FDR_BH"] < FDR_THRESHOLD
        ax.scatter(
            subset.loc[~fdr_pass, "log2FC_pooled"],
            subset.loc[~fdr_pass, "neglog10_p"],
            s=85,
            c=color,
            edgecolors="#404040",
            linewidths=0.7,
            alpha=0.95,
            zorder=3,
        )
        # FDR-surviving points get a heavier black outline to visually
        # distinguish them from raw-p-only "significant" points.
        ax.scatter(
            subset.loc[fdr_pass, "log2FC_pooled"],
            subset.loc[fdr_pass, "neglog10_p"],
            s=85,
            c=color,
            edgecolors="black",
            linewidths=2.0,
            alpha=0.95,
            zorder=3,
        )

        for _, row in subset.iterrows():
            metabolite = str(row["Metabolite"])
            offset = LABEL_OFFSETS.get(metabolite, (8, 8))
            horizontal_alignment = "left" if offset[0] >= 0 else "left"
            label_text = f"{metabolite}\n(p={row['p_value_Welch']:.3f}, FDR={row['FDR_BH']:.2f})"

            ax.annotate(
                label_text,
                xy=(row["log2FC_pooled"], row["neglog10_p"]),
                xytext=offset,
                textcoords="offset points",
                fontsize=8,
                color=color,
                ha=horizontal_alignment,
                va="bottom",
                arrowprops={
                    "arrowstyle": "-",
                    "color": color,
                    "linewidth": 0.8,
                    "alpha": 0.55,
                },
                zorder=4,
            )

    # Threshold lines.
    y_threshold = -np.log10(P_VALUE_THRESHOLD)
    ax.axvline(
        SIGNIFICANT_LOG2FC,
        color="#303030",
        linestyle="--",
        linewidth=1.2,
        zorder=2,
    )
    ax.axhline(
        y_threshold,
        color="#303030",
        linestyle="--",
        linewidth=1.2,
        zorder=2,
    )

    ax.text(
        SIGNIFICANT_LOG2FC,
        ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 1,
        f"log\u2082FC = {SIGNIFICANT_LOG2FC:g}",
        ha="center",
        va="top",
        fontsize=9,
    )
    ax.text(
        MIN_DISPLAY_LOG2FC + 0.03,
        y_threshold + 0.03,
        f"p = {P_VALUE_THRESHOLD:g}",
        ha="left",
        va="bottom",
        fontsize=9,
    )

    ax.set_xlim(
        MIN_DISPLAY_LOG2FC,
        max(3.25, float(plotted["log2FC_pooled"].max()) + 0.15),
    )
    ax.set_ylim(
        -0.05,
        max(2.62, float(plotted["neglog10_p"].max()) + 0.12),
    )

    ax.set_xlabel("log\u2082FC (infected / uninfected)")
    ax.set_ylabel("\u2212log\u2081\u2080(p-value)")

    # Deliberately no title: this is a one-sided significance scatterplot,
    # not a conventional two-sided volcano plot.

    n_fdr_pass = int((highlighted["FDR_BH"] < FDR_THRESHOLD).sum())
    n_highlighted = len(highlighted)
    fdr_note = (
        f"{n_fdr_pass}/{n_highlighted} highlighted metabolites remain below "
        f"FDR (Benjamini\u2013Hochberg) = {FDR_THRESHOLD:g};\n"
        f"none remain below FDR = 0.05. Raw p-value threshold shown is uncorrected."
    )
    ax.text(
        0.02, 0.02, fdr_note,
        transform=ax.transAxes,
        fontsize=8,
        color="#303030",
        ha="left",
        va="bottom",
        style="italic",
        bbox={"facecolor": "white", "edgecolor": "#cccccc", "boxstyle": "round,pad=0.4", "alpha": 0.9},
    )

    legend_handles = [
        Line2D(
            [0], [0],
            marker="o",
            color="none",
            markerfacecolor=OTHER_COLOR,
            markeredgecolor="none",
            markersize=8,
            label="Other metabolites",
            alpha=0.55,
        )
    ]
    for chemical_class, color in CLASS_COLORS.items():
        legend_handles.append(
            Line2D(
                [0], [0],
                marker="o",
                color="none",
                markerfacecolor=color,
                markeredgecolor="#404040",
                markersize=8,
                label=chemical_class,
            )
        )
    legend_handles.append(
        Line2D(
            [0], [0],
            marker="o",
            color="none",
            markerfacecolor="white",
            markeredgecolor="black",
            markeredgewidth=2.0,
            markersize=8,
            label=f"FDR < {FDR_THRESHOLD:g} (thick outline)",
        )
    )

    ax.legend(
        handles=legend_handles,
        title="Chemical class",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
    )


    fig.subplots_adjust(left=0.10, right=0.70, bottom=0.10, top=0.96)

    significant[
        ["Metabolite", "log2FC_pooled", "p_value_Welch", "neglog10_p", "FDR_BH", "Class"]
    ].sort_values(
        ["p_value_Welch", "log2FC_pooled"],
        ascending=[True, False],
    ).to_csv(OUTPUT_CSV, index=False)

    fig.savefig(OUTPUT_PNG, dpi=600, bbox_inches="tight")
    fig.savefig(OUTPUT_SVG, bbox_inches="tight")

    print(
        "Saved:",
        OUTPUT_PNG.name,
        OUTPUT_SVG.name,
        OUTPUT_CSV.name,
    )
    print(
        f"Threshold-significant metabolites: {len(significant)}; "
        f"highlighted metabolites: {len(highlighted)}"
    )


    plt.close(fig)


if __name__ == "__main__":
    main()
