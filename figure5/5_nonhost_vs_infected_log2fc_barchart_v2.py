from pathlib import Path

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

# ---------- Configuration ----------
BASE_DIR = Path(__file__).resolve().parent
RAW_FILE = BASE_DIR / "rafflesia_dataset1.xlsx"     # full untargeted feature matrix
TOP_FILE = BASE_DIR / "rafflesia_dataset3.xlsx"     # pre-selected top nonhost-vs-infected compounds
SHEET_NAME = "Sheet1"

FDR_ALPHA = 0.05
MIN_N_PER_GROUP = 2  # CAM groups only have n=2 biological replicates (n=4 technical runs)


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction (no external stats dependency)."""
    p = np.asarray(pvals, dtype=float)
    n = np.sum(~np.isnan(p))
    q = np.full_like(p, np.nan)
    valid_idx = np.where(~np.isnan(p))[0]
    order = valid_idx[np.argsort(p[valid_idx])]
    ranked_p = p[order]
    ranks = np.arange(1, n + 1)
    raw_q = ranked_p * n / ranks
    q_sorted = np.minimum.accumulate(raw_q[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0, 1)
    q[order] = q_sorted
    return q


def cols_for(prefix: str, df: pd.DataFrame) -> list:
    """Exact technical-replicate columns for a group, e.g. 'infectedCAM', 'infectedCAM.1', ..."""
    return [c for c in df.columns if c == prefix or c.startswith(prefix + ".")]


def main() -> None:
    # ---------- 1) Recompute the FULL family of tests from raw data ----------
    # A proper FDR correction requires the p-value distribution across every
    # metabolite actually tested for this contrast, not just the pre-selected
    # "top hits." We therefore recompute Welch's t-tests for every feature in
    # the deduplicated raw matrix, per locality, then correct across that
    # full family before looking up q-values for the highlighted compounds.
    raw = pd.read_excel(RAW_FILE, sheet_name=SHEET_NAME)
    n_before = len(raw)
    raw = raw.drop_duplicates().reset_index(drop=True)
    print(f"Dropped {n_before - len(raw)} exact-duplicate rows from raw matrix "
          f"({n_before} -> {len(raw)})")

    group_cols = {
        "CAM": {"infected": cols_for("infectedCAM", raw), "nonhost": cols_for("NonhostCAM", raw)},
        "ILO": {"infected": cols_for("infectedILO", raw), "nonhost": cols_for("nonhostILO", raw)},
    }

    family = pd.DataFrame({"Feature_Label": raw["Feature_Label"]})
    for loc, g in group_cols.items():
        n_inf, n_non = len(g["infected"]), len(g["nonhost"])
        print(f"{loc}: infected n={n_inf} technical runs, nonhost n={n_non} technical runs")
        inf = raw[g["infected"]].apply(pd.to_numeric, errors="coerce")
        non = raw[g["nonhost"]].apply(pd.to_numeric, errors="coerce")
        log2fc = np.log2((non.mean(axis=1) + 1) / (inf.mean(axis=1) + 1))
        pvals = np.full(len(raw), np.nan)
        for i in range(len(raw)):
            a = non.iloc[i].dropna().values
            b = inf.iloc[i].dropna().values
            if len(a) < MIN_N_PER_GROUP or len(b) < MIN_N_PER_GROUP:
                continue
            _, p = stats.ttest_ind(a, b, equal_var=False)
            pvals[i] = p
        family[f"{loc}_log2FC"] = log2fc.values
        family[f"{loc}_p"] = pvals
        family[f"{loc}_FDR_BH"] = bh_fdr(pvals)

    n_tested_cam = family["CAM_p"].notna().sum()
    n_tested_ilo = family["ILO_p"].notna().sum()
    print(f"Full test family sizes: CAM n={n_tested_cam}, ILO n={n_tested_ilo} features")

    # ---------- 2) Load the pre-selected top compounds and attach corrected stats ----------
    df = pd.read_excel(TOP_FILE, sheet_name=SHEET_NAME)
    df.columns = df.columns.str.strip()
    df = df[
        ~df["Compound"].astype(str).str.contains("NCGC00380170", case=False, na=False)
    ].copy()

    # Match each displayed compound to its exact source row via its original
    # (log2FC, p) fingerprint, since compound names are not unique keys
    # (multiple ion features can share a compound name).
    matched_rows = []
    for _, row in df.iterrows():
        dist = (
            (family["CAM_log2FC"] - row["CAM_log2FC"]).abs()
            + (family["ILO_log2FC"] - row["ILO_log2FC"]).abs()
            + (family["CAM_p"] - row["CAM_p"]).abs() * 5
            + (family["ILO_p"] - row["ILO_p"]).abs() * 5
        )
        matched_rows.append(family.loc[dist.idxmin()])
    matched = pd.DataFrame(matched_rows).reset_index(drop=True)

    df = df.reset_index(drop=True)
    df["CAM_FDR"] = matched["CAM_FDR_BH"]
    df["ILO_FDR"] = matched["ILO_FDR_BH"]

    # ---------- 3) Plot ----------
    fig, ax = plt.subplots(figsize=(14, 7.5))

    y = np.arange(len(df))
    h = 0.35

    def bar_style(fdr_series):
        # Solid, fully opaque bars for FDR < alpha; hatched, faded bars otherwise.
        alphas = np.where(fdr_series < FDR_ALPHA, 1.0, 0.35)
        hatches = np.where(fdr_series < FDR_ALPHA, "", "//")
        return alphas, hatches

    cam_alpha, cam_hatch = bar_style(df["CAM_FDR"])
    ilo_alpha, ilo_hatch = bar_style(df["ILO_FDR"])

    for i in range(len(df)):
        ax.barh(
            y[i] - h / 2, df["CAM_log2FC"].iloc[i], height=h,
            color="tab:blue", alpha=cam_alpha[i], hatch=cam_hatch[i],
            edgecolor="black", linewidth=0.4,
            label="NonhostCAM/infectedCAM" if i == 0 else None,
        )
        ax.barh(
            y[i] + h / 2, df["ILO_log2FC"].iloc[i], height=h,
            color="tab:orange", alpha=ilo_alpha[i], hatch=ilo_hatch[i],
            edgecolor="black", linewidth=0.4,
            label="NonhostILO/infectedILO" if i == 0 else None,
        )

    ax.set_yticks(y)
    ax.set_yticklabels(df["Compound"], fontsize=8)
    ax.set_xlabel("log2FC")

    def format_p(p: float) -> str:
        # Decimal (not scientific) notation, with enough precision that very
        # small p-values don't collapse to "0.0000".
        if p < 0.0001:
            return f"{p:.6f}"
        return f"{p:.4f}"

    for i, r in enumerate(df.itertuples()):
        ax.text(
            r.CAM_log2FC + 0.15, i - h / 2,
            f"p={format_p(r.CAM_p)}, q={r.CAM_FDR:.2f}",
            va="center", fontsize=6.5,
        )
        ax.text(
            r.ILO_log2FC + 0.15, i + h / 2,
            f"p={format_p(r.ILO_p)}, q={r.ILO_FDR:.2f}",
            va="center", fontsize=6.5,
        )

    ax.set_xlim(0, ax.get_xlim()[1] + 2.0)

    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor="tab:blue", edgecolor="black", label="NonhostCAM/infectedCAM"),
        Patch(facecolor="tab:orange", edgecolor="black", label="NonhostILO/infectedILO"),
        Patch(facecolor="white", edgecolor="black", hatch="//", label=f"FDR \u2265 {FDR_ALPHA:g} (not significant)"),
    ]
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)

    note = (
        f"q-values: Benjamini\u2013Hochberg FDR across the full tested feature family "
        f"(CAM n={n_tested_cam}, ILO n={n_tested_ilo} deduplicated metabolite features),\n"
        f"not just the {len(df)} compounds shown. CAM comparisons use only "
        f"n=2 biological replicates per group (4 technical runs) \u2014 treat as exploratory."
    )
    ax.text(
        0.01, -0.13, note,
        transform=ax.transAxes,
        fontsize=8, style="italic", color="#303030", ha="left", va="top",
    )

    fig.tight_layout()

    # ---------- 4) Save outputs ----------
    png_out = BASE_DIR / "nonhost_vs_infected_log2fc_barchart_FDR_decimal.png"
    svg_out = BASE_DIR / "nonhost_vs_infected_log2fc_barchart_FDR_decimal.svg"
    csv_out = BASE_DIR / "nonhost_vs_infected_metabolites_FDR_decimal.csv"

    fig.savefig(png_out, dpi=300, bbox_inches="tight")
    fig.savefig(svg_out, bbox_inches="tight")
    df.to_csv(csv_out, index=False)

    print("Saved:", png_out.name, svg_out.name, csv_out.name)
    print(
        f"Bars surviving FDR<{FDR_ALPHA:g}: "
        f"{int((df['CAM_FDR'] < FDR_ALPHA).sum() + (df['ILO_FDR'] < FDR_ALPHA).sum())} / {2 * len(df)}"
    )

    plt.close(fig)


if __name__ == "__main__":
    main()
