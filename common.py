"""Shared data helpers: group means computed directly from raw_data/rafflesia_dataset1.xlsx
(never from precomputed 'ave' columns), and annotation-name cleaning so that names differing
only in letter case or a '.mol' suffix are consolidated together."""
from pathlib import Path
import re
import numpy as np
import pandas as pd

RAW = Path(__file__).resolve().parent / "raw_data" / "rafflesia_dataset1.xlsx"

def load_raw():
    raw = pd.read_excel(RAW).drop_duplicates().reset_index(drop=True)
    raw = raw.loc[:, [c for c in raw.columns if not str(c).lower().startswith("ave")]]
    raw["Name"] = raw["Name"].astype(str).str.strip().str.replace(r"\.mol$", "", regex=True)
    raw["name_key"] = raw["Name"].str.lower()
    # display name per key: prefer a form that is not ALL CAPS
    disp = {}
    for k, grp in raw.groupby("name_key")["Name"]:
        forms = sorted(set(grp), key=lambda s: (s.isupper(), sum(ch.isupper() for ch in s), s))
        disp[k] = forms[0]
    raw["Name"] = raw["name_key"].map(disp)
    return raw

def run_cols(df, g):
    return [c for c in df.columns if c == g or re.fullmatch(re.escape(g) + r"\.\d+", str(c))]

def group_mean(df, g, clip=False):
    X = df[run_cols(df, g)].apply(pd.to_numeric, errors="coerce")
    if clip:
        X = X.clip(lower=0)
    return X.mean(axis=1)


# Annotations judged not plausible plant natural products or unresolved library identifiers
# (flagged during screening, confirmed on manual review). Used by Figs. 5 and 6.
NOT_NATURAL = {
    "1-(4-(5-maleimidopentyl)aminobenzyl)ethylenediaminetetraacetic acid",  # synthetic chelator
    "c.i. acid green 3",                                                   # synthetic dye
}

def is_unresolved(name):
    """Library identifiers rather than compound names (e.g. 'NCGC00380149-01_C22H26O11_...')."""
    return bool(re.match(r"^NCGC\d", str(name))) or bool(re.search(r"_C\d+H\d+", str(name)))

def is_natural(name):
    return str(name).strip().lower() not in NOT_NATURAL and not is_unresolved(name)
