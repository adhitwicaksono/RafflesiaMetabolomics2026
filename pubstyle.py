"""Shared print style: figures drawn at final print width (174 mm, double column),
lettering >= 7 pt, readable group labels (species names italic)."""
import matplotlib as mpl

WIDTH = 6.85  # inches (174 mm)

def apply_style():
    mpl.rcParams.update({
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7,
        "legend.fontsize": 7, "legend.title_fontsize": 7.5,
        "savefig.dpi": 600, "mathtext.default": "regular",
    })

def it(s):
    return r"$\it{" + s.replace(" ", r"\ ") + "}$"

LABELS = {
    "sapbud": it("Sapria") + " buds (THAI)",
    "Raffbudspec": it("R. speciosa") + " buds (ILO)",
    "Raffbudlag": it("R. lagascae") + " buds (CAM)",
    "infectedTHAI": "Infected host (THAI)", "infectedCAM": "Infected host (CAM)", "infectedILO": "Infected host (ILO)",
    "UNinfectedTHAI": "Uninfected host (THAI)", "UNinfectedCAM": "Uninfected host (CAM)", "UNinfectedILO": "Uninfected host (ILO)",
    "uninfecraffspec-stemleaf": "Uninfected host stem/leaf (ILO)",
    "nonhostILO": "Non-host (ILO)", "NonhostCAM": "Non-host (CAM)",
    "raffseed": it("R. speciosa") + " seeds",
    "Ampelopsis": it("Ampelopsis"),
    "BUD": "Buds", "INFECTED": "Infected hosts", "UNINFECTED": "Uninfected hosts",
    "nonhostTET": "Non-hosts", "UNINFRAFFSPEC": "Uninfected host stem/leaf", "RAFFSEED": "Seeds",
    "averaffseed": "Seeds", "aveRaffbudspec": "Buds", "aveinfectedILO": "Infected host tissue",
    "aveUNinfectedILO": "Uninfected host tissue", "aveuninfecraffspec-stemleaf": "Uninfected stem/leaf",
    "avenonhostILO": "Non-host root",
}

def label(g):
    return LABELS.get(g, g)
