"""
plot_covid_sensitivity.py — Three-period bar chart:
  Pre-ShotSpotter (2009-2016) / Post-SS excl COVID (2017-19, 22-23) / COVID (2020-21)
Shows annualized homicide rates to make periods of different lengths comparable.
"""
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_style import apply_academic_style, COLORS, add_source_note, SOURCE_DEFAULT

P = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CLEAN = os.path.join(P, "data", "processed")
PLOTS = os.path.join(P, "plots")

apply_academic_style()

C_BLUE   = COLORS["control"]    # pre-ShotSpotter
C_ORANGE = COLORS["orange"]     # post-SS excl COVID
C_RED    = COLORS["treated"]    # COVID period

SELECTED = ["Austin", "Humboldt Park", "North Lawndale",
            "Englewood", "West Englewood", "South Shore"]

# Load data
hom = pd.read_csv(os.path.join(CLEAN, "Homicide_six_hoods.csv"))
hom["Date"] = pd.to_datetime(hom["Date"])
hom["Year"] = hom["Date"].dt.year

# Compute annualized rates for each period
PRE_YEARS = list(range(2009, 2017))       # 8 years
POST_YEARS = [2017, 2018, 2019, 2022, 2023]  # 5 years
COVID_YEARS = [2020, 2021]                # 2 years

results = {}
gmax = 0
for hood in SELECTED:
    h = hom[hom["Community"] == hood]
    yearly = h.groupby("Year").size().reindex(range(2009, 2024), fill_value=0)

    pre   = yearly[PRE_YEARS].sum() / len(PRE_YEARS)
    post  = yearly[POST_YEARS].sum() / len(POST_YEARS)
    covid = yearly[COVID_YEARS].sum() / len(COVID_YEARS)
    results[hood] = (pre, post, covid)
    gmax = max(gmax, pre, post, covid)

# Plot
fig, axes = plt.subplots(2, 3, figsize=(14, 10))
axes = axes.flatten()

periods = ["Pre-SS\n(2009-16)", "Post-SS\nexcl. COVID\n(2017-19, 22-23)", "COVID\n(2020-21)"]
colors = [C_BLUE, C_ORANGE, C_RED]

for idx, hood in enumerate(SELECTED):
    ax = axes[idx]
    pre, post, covid = results[hood]
    vals = [pre, post, covid]

    bars = ax.bar(periods, vals, color=colors, width=0.55,
                  edgecolor="black", linewidth=0.5)
    ax.set_ylim(0, gmax * 1.18)
    ax.set_ylabel("Gun Homicides / Year", fontsize=9)
    ax.yaxis.grid(True, alpha=0.3, linewidth=0.5)
    ax.set_axisbelow(True)
    ax.set_title(hood, fontsize=12, fontweight="bold", pad=8)

    # Value labels on bars
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + gmax * 0.015,
                f"{val:.1f}", ha="center", va="bottom",
                fontsize=9, fontweight="bold")

    # Percent change annotations
    if pre > 0:
        pct_post = (post - pre) / pre * 100
        pct_covid = (covid - pre) / pre * 100
        ax.text(0.5, -0.13,
                f"vs Pre-SS:  Post {pct_post:+.0f}%  |  COVID {pct_covid:+.0f}%",
                fontsize=7, color="#555555", ha="center",
                transform=ax.transAxes)

# Legend
patches = [
    mpatches.Patch(color=C_BLUE, label="Pre-ShotSpotter (2009-2016)"),
    mpatches.Patch(color=C_ORANGE, label="Post-SS excl. COVID (2017-19, 22-23)"),
    mpatches.Patch(color=C_RED, label="COVID Period (2020-2021)"),
]
fig.legend(handles=patches, loc="lower center", ncol=3,
           frameon=False, fontsize=10, bbox_to_anchor=(0.5, -0.02))

fig.suptitle("Gun Homicide Rates: Pre-ShotSpotter vs. Post-SS vs. COVID\n"
             "(Annualized counts per period)",
             fontsize=15, fontweight="bold", y=0.99)
add_source_note(fig, SOURCE_DEFAULT,
                sample_n="Annualized: total / years in period",
                y=0.005, x=0.5)
for txt in fig.texts:
    if txt.get_text().startswith("Source:"):
        txt.set_ha("center")

fig.tight_layout(rect=[0, 0.05, 1, 0.95])

out = os.path.join(PLOTS, "covid_sensitivity.png")
fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
plt.close(fig)
sz = os.path.getsize(out) / 1024
print(f"OK: covid_sensitivity.png  ({sz:.0f} KB)")
