"""Independent check that figure CSVs match values recomputed from raw_data/rafflesia_dataset1.xlsx.
Run after the figure scripts. Prints PASS/FAIL per check."""
import re, numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
B = Path(__file__).resolve().parent
raw = pd.read_excel(B / "raw_data/rafflesia_dataset1.xlsx").drop_duplicates().reset_index(drop=True)
raw["key"] = raw.Name.astype(str).str.strip().str.replace(r"\.mol$", "", regex=True).str.lower()
cols = lambda g: [c for c in raw.columns if c == g or re.fullmatch(re.escape(g) + r"\.\d+", str(c))]
ok = lambda c, msg: print(("PASS " if c else "FAIL ") + msg)
num = lambda g: raw[cols(g)].apply(pd.to_numeric, errors="coerce")

# Fig. 4: recompute Welch p on biological means for the labelled compounds and all nominal hits
f4 = pd.read_csv(B / "figure4/fig4_bio_consolidated.csv")
d = pd.concat([num(g).clip(lower=0) for g in ["infectedTHAI","infectedCAM","infectedILO","UNinfectedTHAI","UNinfectedCAM","UNinfectedILO"]], axis=1)
d["key"] = raw.key; d = d.groupby("key").mean()
bio = lambda gs: np.hstack([(lambda X: (X[:, 0::2] + X[:, 1::2]) / 2)(d[[c for c in d.columns if re.fullmatch(re.escape(g) + r"(\.\d+)?", c)]].values) for g in gs])
A, Bm = bio(["infectedTHAI","infectedCAM","infectedILO"]), bio(["UNinfectedTHAI","UNinfectedCAM","UNinfectedILO"])
p = stats.ttest_ind(A, Bm, axis=1, equal_var=False).pvalue
chk = pd.Series(p, index=d.index).dropna()
f4["key"] = f4.Name.str.lower()
m = f4.set_index("key").join(chk.rename("p_check"), how="inner")
ok(len(f4) == len(chk), f"Fig 4 tested-metabolite count {len(f4)} (independent {len(chk)})")
ok(np.allclose(m.p, m.p_check, rtol=1e-6), "Fig 4 p-values match independent recomputation")
q = stats.false_discovery_control(chk.values, method="bh")
ok(np.isclose(f4.q_BH.min(), q.min()), f"Fig 4 minimum q {f4.q_BH.min():.3f}")
hits = sorted(f4[(f4.p < 0.05) & (f4.log2FC > 1)].Name)
print("     Fig 4 nominal hits (all must be labelled in the figure):", hits)

# Fig. 5: recompute locality tests; check the selection rule and that the plotted nine are the top nine
f5 = pd.read_csv(B / "figure5/fig5_bio_displayed.csv"); allf = pd.read_csv(B / "figure5/fig5_bio_all_features.csv")
cand = pd.read_csv(B / "figure5/fig5_rule_based_candidates.csv")
ok(len(allf) == len(raw), f"Fig 5 feature table covers all {len(raw)} deduplicated features")
P, FC = {}, {}
for loc, (n, i) in {"CAM": ("NonhostCAM", "infectedCAM"), "ILO": ("nonhostILO", "infectedILO")}.items():
    N = num(n).clip(lower=0).values; I = num(i).clip(lower=0).values
    N = (N[:, 0::2] + N[:, 1::2]) / 2; I = (I[:, 0::2] + I[:, 1::2]) / 2
    P[loc] = stats.ttest_ind(N, I, axis=1, equal_var=False).pvalue
    FC[loc] = np.log2((N.mean(1) + 1) / (I.mean(1) + 1))
rule = (FC["CAM"] > 1) & (FC["ILO"] > 1) & ((P["CAM"] < 0.05) | (P["ILO"] < 0.05))
ok(set(raw.Feature_Label[rule]) == set(cand.Feature_Label), f"Fig 5 rule-based list complete ({rule.sum()} features)")
nat = cand[cand.natural_product].sort_values("min_p")
ok(list(nat.Feature_Label[:9]) == list(f5.Feature_Label), "Fig 5 plots the nine natural-product features with lowest p")
idx = [raw.index[raw.Feature_Label == fl][0] for fl in f5.Feature_Label]
ok(all(np.allclose(f5[f"{l}_p"], P[l][idx], rtol=1e-6) for l in P), "Fig 5 p-values match independent recomputation")
print("     excluded as non-natural:", sorted(cand[~cand.natural_product].Name))

# Fig. 6: top 20 over ALL annotated features (excluding non-natural/unresolved), means from raw runs
import sys; sys.path.insert(0, str(B)); from common import is_natural
t = pd.read_csv(B / "figure6/top20_compounds_raw_means.csv", index_col=0); tss = pd.read_csv(B / "figure6/top20_compounds_TSS.csv", index_col=0)
G6 = ["raffseed","Raffbudspec","infectedILO","UNinfectedILO","uninfecraffspec-stemleaf","nonhostILO"]
mm = pd.DataFrame({"ave" + g: num(g).mean(axis=1) for g in G6}); mm["key"] = raw.key
mm = mm.groupby("key").mean().clip(lower=0)
ok(np.allclose(mm.loc[t.index.str.lower()].values, t.values, rtol=1e-6), "Fig 6 group means match raw data")
rank = mm.sum(axis=1).sort_values(ascending=False)
natural = [k for k in rank.index if is_natural(raw.Name[raw.key == k].iloc[0])]
ok(set(natural[:20]) == set(t.index.str.lower()), "Fig 6 shows the 20 most abundant natural-product annotations")
ok(np.allclose(tss.sum(axis=0), 1), "Fig 6 TSS columns sum to 1")
print("     Fig 6 ranks 1-25 (all features; review for non-natural names):", [raw.Name[raw.key == k].iloc[0][:40] for k in rank.index[:25]])

# Fig. 7: Z-scores from raw group means
z = pd.read_csv(B / "figure7/metabolite_heatmap_matrix.csv", index_col=0)
G7 = ["infectedILO","UNinfectedILO","uninfecraffspec-stemleaf","nonhostILO"]
mg = pd.DataFrame({g: num(g).mean(axis=1) for g in G7}); mg["key"] = raw.key; mg = mg.groupby("key").mean().loc[z.index.str.lower()]
zz = mg.sub(mg.mean(axis=1), axis=0).div(mg.std(axis=1, ddof=0), axis=0)
ok(np.allclose(zz.values, z.values, atol=1e-6), "Fig 7 Z-scores match raw group means")
stored = raw[[c for c in raw.columns if str(c).startswith("ave")]]
bad = [c for c in stored.columns if cols(c[3:]) and not np.allclose(stored[c], num(c[3:]).mean(axis=1), atol=1)]
print("     note: precomputed 'ave' columns in dataset1 that do NOT equal the mean of their runs:", bad, "(not used by any script)")
