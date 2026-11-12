"""Blank-filter sensitivity analysis (revision 2; Methods 2.3, Results 3.1).

Dataset S1 intensities are blank-subtracted: the mean of the four procedural blank runs
was subtracted from every run, so the blank columns themselves are centred on zero.
For each feature the blank mean is recovered as -min(all runs), which is exact whenever
at least one run had zero raw signal (true for 1,475/1,518 features; the rest give a
lower bound). Which features had blank signal does not depend on this reconstruction:
a feature had signal in at least one blank if and only if its four blank-subtracted blank
values differ (1,012 features; the only exception would be identical non-zero values in all
four blanks). The strict filter (ii) therefore uses the blank columns directly, whereas the
3x ratio in filter (i) uses the reconstructed blank mean. Analyses are then repeated after excluding
  (i)  features whose highest group mean (raw) is < 3x the blank mean, and
  (ii) all features with any blank signal.
Also reports duplicate-injection agreement in total feature intensity, and Ampelopsis
distances to host vs non-host Tetrastigma (Discussion 4.1).
Run from this folder:  python3 blank_sensitivity.py
"""
import re, numpy as np, pandas as pd
from scipy.spatial import procrustes
from scipy.spatial.distance import cdist
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from bio_stats import load, run_cols, bio_matrix
import types
from pathlib import Path

raw = load()
# reuse pseudo-F / permutation functions from permanova.py without running its tests
_pm = types.SimpleNamespace()
_src = open(Path(__file__).with_name('permanova.py')).read().split('tests=')[0].replace('raw = load()', '')
_ns = {}
exec(_src, _ns)
class _PM:
    def __setattr__(self, k, v):
        _ns[k] = v
    def __getattr__(self, k):
        return _ns[k]
_pm = _PM()
META = ['RT', 'm/z', 'M meas', 'Ions', 'MS/MS', 'Name', 'Molecular_Formula', 'Feature_Label']
BL = ['Blank', 'Blank.1', 'Blank.2', 'Blank.3']
runs = [c for c in raw.columns if c not in META and not c.startswith('ave') and c not in BL]
groups = sorted({re.sub(r'\.\d+$', '', c) for c in runs})
V = raw[runs + BL].apply(pd.to_numeric, errors='coerce')
bmean = -V.min(axis=1).clip(upper=0)
R = V.add(bmean, axis=0)
gmax = pd.DataFrame({g: R[run_cols(raw, g)].mean(axis=1) for g in groups}).max(axis=1)
ratio = (gmax / bmean.replace(0, np.nan)).fillna(np.inf)
keep3 = (ratio >= 3).values
strict = (V[BL].max(axis=1) == V[BL].min(axis=1)).values   # no blank signal: four blank values identical (read directly from blank columns)
print(f'features {len(raw)}; with blank signal {(bmean > 0).sum()}; ratio<3 {(~keep3).sum()}; no blank signal {strict.sum()}')

def pca3(mask):
    X = raw.loc[mask, runs].apply(pd.to_numeric, errors='coerce').fillna(0).T.values
    X = X[:, X.std(0) > 0]
    p = PCA(3, random_state=0); s = p.fit_transform(StandardScaler().fit_transform(X))
    return s, p.explained_variance_ratio_ * 100
s0, v0 = pca3(np.ones(len(raw), bool))
for lab, m in (('ratio>=3', keep3), ('no blank signal', strict)):
    s, v = pca3(m)
    print(f'PCA all samples [{lab}]: PC1-3 = {v.round(1)} (sum {v.sum():.1f}); Procrustes m2 vs full = {procrustes(s0, s)[2]:.3f}')

TESTS = {'infected vs uninfected (pooled, strata)': (['infectedTHAI', 'infectedCAM', 'infectedILO', 'UNinfectedTHAI', 'UNinfectedCAM', 'UNinfectedILO'], True),
         'non-host vs infected (CAM+ILO, strata)': (['NonhostCAM', 'infectedCAM', 'nonhostILO', 'infectedILO'], True),
         'non-host vs infected ILO': (['nonhostILO', 'infectedILO'], False)}
for lab, m in (('all', np.ones(len(raw), bool)), ('ratio>=3', keep3), ('no blank signal', strict)):
    _pm.raw = raw.loc[m].reset_index(drop=True)
    for k, (gs, st) in TESTS.items():
        X, labs, loc = _pm.mat(gs); strata = loc if st else np.zeros(len(labs))
        for met in ('euc', 'bc'):
            F, R2, p = _pm.perm(X, labs, strata, met, nperm=9999)
            print(f'PERMANOVA [{lab}] {k} {met}: F={F:.2f} R2={R2:.3f} p={p:.4f}')

# candidates
key = raw.Name.astype(str).str.strip().str.replace(r'\.mol$', '', regex=True).str.lower()
for n in ["Nb-trans-Feruloylserotonin glucoside", "Phloroacetophenone 6'-[xylosyl-(1->6)-glucoside]", "Gibberellin A39", "3,5-Digalloylepicatechin"]:
    print(f'Fig. 4 candidate {n}: min sample/blank ratio {ratio[key == n.lower()].min():.0f}')
t1 = pd.read_csv(Path('../figure5/fig5_rule_based_candidates.csv'))
r1 = [ratio[(raw.RT.round(2) == round(r.RT, 2)) & (key == str(r.Name).lower().strip())].min() for _, r in t1.iterrows()]
print(f'Table S1: {len(t1)} features; min sample/blank ratio {np.nanmin(r1):.1f}')

# duplicate-injection agreement in total intensity
rd = []
for g in groups:
    s = raw[run_cols(raw, g)].apply(pd.to_numeric, errors='coerce').clip(lower=0).sum().values
    a, b = s[0::2], s[1::2]; rd += list(np.abs(a - b) / ((a + b) / 2))
rd = np.array(rd) * 100
print(f'duplicate injections, total intensity difference: median {np.median(rd):.1f}%, 90th pct {np.percentile(rd, 90):.1f}%')

# Ampelopsis vs host / non-host distances (biological samples, autoscaled log10 as PERMANOVA)
_pm.raw = raw
G = {'amp': ['Ampelopsis'], 'host': TESTS['infected vs uninfected (pooled, strata)'][0], 'non': ['NonhostCAM', 'nonhostILO']}
M = {k: np.hstack([bio_matrix(raw, g) for g in v]).T for k, v in G.items()}
X = np.vstack(list(M.values())); keep = np.log10(X + 1).std(0) > 0
Xl = np.log10(X + 1)[:, keep]; Z = (Xl - Xl.mean(0)) / Xl.std(0)
i = np.cumsum([0] + [len(M[k]) for k in G])
A, H, N = (Z[i[j]:i[j + 1]] for j in range(3)); Ar, Hr, Nr = (X[:, keep][i[j]:i[j + 1]] for j in range(3))
print(f'Ampelopsis mean distance to hosts vs non-hosts: Euclidean {cdist(A, H).mean():.1f} vs {cdist(A, N).mean():.1f}; '
      f'Bray-Curtis {cdist(Ar, Hr, "braycurtis").mean():.2f} vs {cdist(Ar, Nr, "braycurtis").mean():.2f}')
