"""Shared: biological-replicate-level statistics for Figs. 4-5.
Technical duplicate injections (adjacent columns, e.g. infectedILO / infectedILO.1)
are averaged per biological sample BEFORE testing; Welch's t-test; BH FDR across
all testable features in each contrast."""
from pathlib import Path
import re, numpy as np, pandas as pd
from scipy import stats

BASE = Path(__file__).resolve().parent
RAW = BASE.parent / "raw_data" / "rafflesia_dataset1.xlsx"

def load():
    raw = pd.read_excel(RAW).drop_duplicates().reset_index(drop=True)
    return raw

def run_cols(raw, g):
    return [c for c in raw.columns if c == g or re.fullmatch(re.escape(g) + r"\.\d+", c)]

def bio_matrix(raw, g, clip=True):
    """n_features x n_biological; adjacent technical runs averaged."""
    X = raw[run_cols(raw, g)].apply(pd.to_numeric, errors="coerce")
    if clip:
        X = X.clip(lower=0)   # negative values = baseline-subtraction artifacts
    X = X.values
    assert X.shape[1] % 2 == 0, g
    return (X[:, 0::2] + X[:, 1::2]) / 2.0

def bh(p):
    p = np.asarray(p, float); q = np.full_like(p, np.nan); ok = ~np.isnan(p)
    if ok.any():
        q[ok] = stats.false_discovery_control(p[ok], method="bh")
    return q

def contrast(A, B):
    """A vs B (rows=features). Welch p, log2((meanA+1)/(meanB+1))."""
    with np.errstate(all="ignore"):
        p = stats.ttest_ind(A, B, axis=1, equal_var=False).pvalue
    fc = np.log2((A.mean(1) + 1) / (B.mean(1) + 1))
    return p, fc, bh(p)
