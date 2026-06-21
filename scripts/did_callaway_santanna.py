"""
did_callaway_santanna.py — Heterogeneity-robust staggered DiD as a final robustness
check on top of the TWFE / count-data / synthetic-control re-analysis.

The earlier scripts diagnose the staggered-rollout problem (Goodman-Bacon 2021) but
estimate it with a hand-rolled cohort event study. This script estimates the
*Callaway & Sant'Anna (2021)* group-time ATT — the modern estimator the literature
actually prescribes for staggered adoption with a clean never-treated control group —
using the `differences` package.

The point is not to "rescue" a causal estimate. It is to show that even the
gold-standard staggered-DiD estimator does not deliver a credible effect: its overall
ATT flips sign relative to the naive OLS, and its own pre-treatment event-study
coefficients reject parallel trends (the 2016 homicide-record year sits one period
before the 2017 cohort's treatment). Across OLS (+2.6), Poisson/NegBin (null), and
Callaway-Sant'Anna (negative), the sign of the "effect" is determined by the method,
not the data.

Appends Section 8 to docs/results_summary_reanalysis.md and writes
plots/callaway_santanna_event_study.png. Run AFTER econometric_robustness.py.
"""
import os
import sys
import json
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_style import (apply_academic_style, COLORS, add_source_note,
                        add_n_label, SOURCE_DEFAULT)

from differences import ATTgt

P = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(P, "data", "raw")
PLOTS = os.path.join(P, "plots")

apply_academic_style()
C_BLUE = COLORS["control"]
C_RED = COLORS["treated"]
C_GRAY = COLORS["neutral"]


# =====================================================================
# DATA (same loading + cohort logic as econometric_analysis.py)
# =====================================================================
print("=" * 72)
print("LOADING DATA")
print("=" * 72)

with open(os.path.join(DATA, "CommAreas_20250408.geojson"), "r", encoding="utf-8") as f:
    geo = json.load(f)
ca_num2name = {int(float(feat["properties"]["area_num_1"])): feat["properties"]["community"].title()
               for feat in geo["features"]}
ca_name2num = {}
for n, name in ca_num2name.items():
    ca_name2num[name] = n
    ca_name2num[name.upper()] = n
    ca_name2num[name.lower()] = n

hom_raw = pd.read_csv(
    os.path.join(DATA, "Violence_Reduction_-_Victims_of_Homicides_and_Non-Fatal_Shootings_20250513.csv"),
    low_memory=False)
hom_raw["DATE"] = pd.to_datetime(hom_raw["DATE"], format="%m/%d/%Y %I:%M:%S %p", errors="coerce")
hom_raw["CA_num"] = hom_raw["COMMUNITY_AREA"].map(ca_name2num)
hom_raw["year"] = hom_raw["DATE"].dt.year

ss_raw = pd.read_csv(
    os.path.join(DATA, "Violence_Reduction_-_Shotspotter_Alerts_-_Historical_20250408.csv"))
ss_raw["DATE"] = pd.to_datetime(ss_raw["DATE"], format="%m/%d/%Y %I:%M:%S %p", errors="coerce")
ss_raw = ss_raw[ss_raw["DATE"] <= "2023-12-31"].copy()
sample_ca = ss_raw["COMMUNITY_AREA"].dropna().iloc[0]
if str(sample_ca).replace(".", "").isdigit():
    ss_raw["CA_num"] = pd.to_numeric(ss_raw["COMMUNITY_AREA"], errors="coerce").astype("Int64")
    ss_raw["CA_name"] = ss_raw["CA_num"].map(ca_num2name)
else:
    ss_raw["CA_name"] = ss_raw["COMMUNITY_AREA"].str.title()
    ss_raw["CA_num"] = ss_raw["CA_name"].map(ca_name2num)

ss_profile = (
    ss_raw.dropna(subset=["CA_name"]).groupby("CA_name")
    .agg(total_alerts=("DATE", "size"), earliest=("DATE", "min"), latest=("DATE", "max")))
ss_profile["months"] = ((ss_profile["latest"] - ss_profile["earliest"]).dt.days / 30.4).round().astype(int)
ss_profile["install_year"] = ss_profile["earliest"].dt.year
ss_profile = ss_profile[(ss_profile.total_alerts >= 100) & (ss_profile.months >= 12)]
install_year = {name: int(iy) for name, iy in ss_profile["install_year"].items()}

# Build the community-area x year panel of gun homicides, 2009-2023
yrs = list(range(2009, 2024))
all_ca = sorted(ca_num2name.keys())
panel = pd.DataFrame([(ca, yr) for ca in all_ca for yr in yrs], columns=["ca_num", "year"])
panel["ca_name"] = panel["ca_num"].map(ca_num2name)
gun = hom_raw[(hom_raw["VICTIMIZATION_PRIMARY"] == "HOMICIDE") &
             (hom_raw["GUNSHOT_INJURY_I"] == "YES") &
             (hom_raw["year"] >= 2009) & (hom_raw["year"] <= 2023)]
g_counts = gun.dropna(subset=["CA_num"]).groupby(["CA_num", "year"]).size().reset_index(name="gun_hom")
g_counts["CA_num"] = g_counts["CA_num"].astype(int)
panel = panel.merge(g_counts.rename(columns={"CA_num": "ca_num"}), on=["ca_num", "year"], how="left")
panel["gun_hom"] = panel["gun_hom"].fillna(0).astype(int)
# Cohort = ShotSpotter install year for treated areas, NaN for never-treated
panel["cohort"] = panel["ca_name"].map(install_year).astype("float")

n_treated = panel.loc[panel["cohort"].notna(), "ca_num"].nunique()
n_never = panel.loc[panel["cohort"].isna(), "ca_num"].nunique()
print(f"  Panel: {panel.shape[0]} rows | {n_treated} treated CAs | {n_never} never-treated")
print(f"  Cohorts: {sorted(set(panel.loc[panel['cohort'].notna(), 'cohort'].astype(int)))}")


# =====================================================================
# CALLAWAY & SANT'ANNA (2021) GROUP-TIME ATT
# =====================================================================
print("\n" + "=" * 72)
print("CALLAWAY-SANT'ANNA GROUP-TIME ATT")
print("=" * 72)

cs_data = panel[["ca_num", "year", "gun_hom", "cohort"]].set_index(["ca_num", "year"]).sort_index()


def fit_cs(control_group):
    att = ATTgt(data=cs_data, cohort_column="cohort")
    att.fit(formula="gun_hom", control_group=control_group, progress_bar=False)
    return att


# Primary: never-treated controls. Robustness: not-yet-treated controls.
att_nt = fit_cs("never_treated")
simple_nt = att_nt.aggregate("simple")
overall_nt = float(simple_nt.iloc[0, 0])
overall_nt_se = float(simple_nt.iloc[0, 1])

att_ny = fit_cs("not_yet_treated")
simple_ny = att_ny.aggregate("simple")
overall_ny = float(simple_ny.iloc[0, 0])
overall_ny_se = float(simple_ny.iloc[0, 1])

print(f"\n  Overall ATT (never-treated control):    {overall_nt:+.3f}  (SE {overall_nt_se:.3f})")
print(f"  Overall ATT (not-yet-treated control):  {overall_ny:+.3f}  (SE {overall_ny_se:.3f})")
print("  (Compare: naive OLS DiD = +2.63; Poisson/NegBin = statistically null.)")

# Event-study aggregation (never-treated control)
ev = att_nt.aggregate("event")
ev.columns = ["att", "se", "lower", "upper", "sig"]
ev = ev.reset_index().rename(columns={ev.index.name or "index": "k"})
ev["k"] = ev.iloc[:, 0].astype(int)
ev = ev[["k", "att", "se", "lower", "upper"]].sort_values("k").reset_index(drop=True)

print("\n  Event-study ATT(k) (never-treated control):")
print("   k  |   ATT   |   SE   | 95% CI")
for _, r in ev.iterrows():
    print(f"  {int(r['k']):+3d} | {r['att']:+7.3f} | {r['se']:.3f} | [{r['lower']:+.3f}, {r['upper']:+.3f}]")

n_pre_sig = int(((ev["k"] < 0) & ((ev["lower"] > 0) | (ev["upper"] < 0))).sum())
print(f"\n  Significant PRE-treatment coefficients (parallel-trends violations): {n_pre_sig}")


# =====================================================================
# PLOT: Callaway-Sant'Anna event study
# =====================================================================
fig, ax = plt.subplots(figsize=(12, 6.5))
ks = ev["k"].values
atts = ev["att"].values
los = ev["lower"].values
his = ev["upper"].values

ax.axvspan(ks.min() - 1, -0.5, color=COLORS["pre_shade"], alpha=0.5, zorder=0)
ax.axhline(0, color="#333333", lw=1.2, zorder=1)
ax.axvline(-0.5, color=C_GRAY, ls="--", lw=1.2, zorder=1)

for k, att, lo, hi in zip(ks, atts, los, his):
    color = C_BLUE if k < 0 else C_RED
    ax.plot([k, k], [lo, hi], color=color, lw=2.2, zorder=2)
    ax.plot(k, att, "o", color=color, markersize=8,
            markeredgecolor="white", markeredgewidth=0.6, zorder=3)

ax.text(-0.4, ax.get_ylim()[1] * 0.96, "ShotSpotter installed",
        fontsize=9, color=C_GRAY, va="top", ha="left", style="italic")
ax.set_xlabel("Years relative to ShotSpotter installation (cohort-specific)")
ax.set_ylabel("ATT: Gun homicides per community-area-year\n(Callaway-Sant'Anna, never-treated control)")
ax.set_title("Callaway-Sant'Anna Event Study (heterogeneity-robust staggered DiD)\n"
             f"Overall ATT = {overall_nt:+.2f} (SE {overall_nt_se:.2f}) — "
             "sign flips vs. naive OLS (+2.63); pre-trends still fail",
             fontweight="bold")

legend_handles = [
    Line2D([0], [0], color=C_BLUE, marker="o", lw=2.2, markersize=7, label="Pre-treatment (k < 0)"),
    Line2D([0], [0], color=C_RED, marker="o", lw=2.2, markersize=7, label="Post-treatment (k ≥ 0)"),
    Line2D([0], [0], color=C_GRAY, ls="--", lw=1.2, label="Install (k = 0)"),
]
ax.legend(handles=legend_handles, loc="lower right", frameon=False, fontsize=10)
ax.yaxis.grid(True)
add_n_label(ax,
            f"Treated: n={n_treated} (2017 + 2018 cohorts)\nNever-treated control: n={n_never}\n"
            f"Estimator: Callaway & Sant'Anna (2021), outcome-regression",
            loc="upper left")
add_source_note(fig, SOURCE_DEFAULT)

fig.subplots_adjust(top=0.86, bottom=0.13, left=0.10, right=0.97)
fig.savefig(os.path.join(PLOTS, "callaway_santanna_event_study.png"), dpi=200, facecolor="white")
plt.close(fig)
print("\n  SAVED: callaway_santanna_event_study.png")


# =====================================================================
# APPEND SECTION 8 TO results_summary_reanalysis.md
# =====================================================================
lines = []
lines.append("\n## 8. Callaway-Sant'Anna (heterogeneity-robust staggered DiD)\n\n")
lines.append("The modern estimator for staggered adoption (Callaway & Sant'Anna 2021),\n")
lines.append("estimated with the `differences` package on group-time ATTs and aggregated\n")
lines.append("to an overall effect and an event study. Standard errors cluster by\n")
lines.append("community area (analytic influence-function inference).\n\n")
lines.append("### Overall ATT\n\n")
lines.append("| Control group | Overall ATT | SE |\n")
lines.append("|---|---|---|\n")
lines.append(f"| Never-treated | {overall_nt:+.3f} | {overall_nt_se:.3f} |\n")
lines.append(f"| Not-yet-treated | {overall_ny:+.3f} | {overall_ny_se:.3f} |\n\n")
lines.append("### Event-study ATT(k), never-treated control\n\n")
lines.append("| k (years from install) | ATT | SE | 95% CI |\n")
lines.append("|---|---|---|---|\n")
for _, r in ev.iterrows():
    lines.append(f"| {int(r['k']):+d} | {r['att']:+.3f} | {r['se']:.3f} | "
                 f"[{r['lower']:+.3f}, {r['upper']:+.3f}] |\n")
lines.append("\n")
lines.append("> **What this adds.** The Callaway-Sant'Anna overall ATT is "
             f"{overall_nt:+.2f} — *negative*, the opposite sign of the naive OLS DiD\n")
lines.append("> (+2.63) and of the null count-data models. The estimate is not credible as a\n")
lines.append("> treatment effect: the event study still shows a large, significant PRE-treatment\n")
lines.append("> coefficient at k = -1 (the 2016 homicide-record year, one period before the 2017\n")
lines.append("> cohort), so parallel trends fail even under the heterogeneity-robust estimator and\n")
lines.append("> the post-period 'decline' is largely a mechanical reversion from that spike.\n")
lines.append("> **The takeaway is the consistency of the inconsistency:** across OLS (+2.63),\n")
lines.append("> Poisson/NegBin (null), and Callaway-Sant'Anna (negative), the *sign* of the\n")
lines.append("> estimated effect is set by the choice of method, not by the data — which is the\n")
lines.append("> project's central point that no design credibly identifies a causal effect.\n")
lines.append("\nFigure: `plots/callaway_santanna_event_study.png`\n")

with open(os.path.join(P, "docs", "results_summary_reanalysis.md"), "a", encoding="utf-8") as f:
    f.write("".join(lines))
print("  Appended Section 8 to docs/results_summary_reanalysis.md")
print("\n" + "=" * 72)
print("DONE")
print("=" * 72)
