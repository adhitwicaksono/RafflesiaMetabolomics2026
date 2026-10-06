"""Fig. 6: top-20 compounds (total-sum scaled) across R. speciosa life stages and ILO host tissues.
Selection rule: group means computed from the individual runs in dataset1 for ALL annotated features;
features sharing an annotation (case-insensitive) averaged; negative means set to zero; compounds
ranked by summed group means after excluding non-natural or unresolved annotations (common.is_natural);
the top 20 are plotted after total-sum scaling within each group."""
import sys, textwrap
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import pandas as pd, matplotlib.pyplot as plt
from pubstyle import apply_style, label as glabel, WIDTH
from common import load_raw, group_mean, is_natural
apply_style()
BASE = Path(__file__).resolve().parent
GROUPS = ["raffseed", "Raffbudspec", "infectedILO", "UNinfectedILO", "uninfecraffspec-stemleaf", "nonhostILO"]

raw = load_raw()
means = pd.DataFrame({"ave" + g: group_mean(raw, g) for g in GROUPS})
means["Compound"] = raw["Name"].values
means = means.groupby("Compound").mean().clip(lower=0)
ranked = means.sum(axis=1).sort_values(ascending=False)
excluded_above = [n for n in ranked.index[:30] if not is_natural(n)]
ranked = ranked[[is_natural(n) for n in ranked.index]]
top20 = means.loc[ranked.index[:20]]
tss = top20.div(top20.sum(axis=0), axis=1)
top20.to_csv(BASE / "top20_compounds_raw_means.csv"); tss.to_csv(BASE / "top20_compounds_TSS.csv")
ranked.head(30).to_csv(BASE / "top30_ranking_after_exclusion.csv", header=["summed_group_means"])
print("excluded from the top ranks:", excluded_above)

fig, ax = plt.subplots(figsize=(WIDTH, 4.4))
tss.T.plot(kind="bar", stacked=True, ax=ax, width=0.85, colormap="tab20")
ax.set_ylabel("Relative abundance (total-sum scaled)"); ax.set_xlabel("")
ax.set_xticklabels([glabel(c) for c in tss.columns], rotation=35, ha="right", rotation_mode="anchor")
h, l = ax.get_legend_handles_labels()
ax.legend(h, [textwrap.fill(x, 34) for x in l], title="Compound", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=6.5, frameon=False)
fig.subplots_adjust(left=0.09, right=0.56, bottom=0.22, top=0.97)
fig.savefig(BASE / "top20_compounds_stacked_TSS.png", dpi=600, bbox_inches="tight", pad_inches=0.05)
print(tss.round(3).to_string())
