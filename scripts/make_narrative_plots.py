"""
make_narrative_plots.py — Two figures that visualize the central findings
of the econometric re-analysis.

  1. plots/specification_comparison.png — translates the OLS DiD into the
     same multiplicative (IRR) scale as the Poisson and Negative Binomial
     estimates and plots all three with confidence intervals. The OLS
     "headline" appears as significant; the proper count-data models do
     not. This is the single-figure visualization of the central finding.

  2. plots/falsification_plot.png — shows the actual post-2017 DiD coefficient
     alongside the four pre-period placebo coefficients (fake treatment
     dates in 1999, 2003, 2007, 2011). When the same design produces
     similar-magnitude effects in placebo windows, the design fails its
     falsification test.
"""
import os
import sys
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import statsmodels.api as sm
import statsmodels.formula.api as smf

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_style import (apply_academic_style, COLORS, format_pvalue,
                        sig_stars, add_source_note, SOURCE_DEFAULT)

P = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(P, "data", "raw")
PLOTS = os.path.join(P, "plots")

apply_academic_style()

C_BLUE = COLORS["control"]
C_RED = COLORS["treated"]
C_GRAY = COLORS["neutral"]
C_LIGHTGRAY = COLORS["grid"]
C_GREEN = COLORS["green"]


# =====================================================================
# DATA LOADING (re-using logic from prior scripts)
# =====================================================================
with open(os.path.join(DATA, "CommAreas_20250408.geojson"), "r", encoding="utf-8") as f:
    geo = json.load(f)
ca_num2name = {int(float(feat["properties"]["area_num_1"])): feat["properties"]["community"].title()
               for feat in geo["features"]}
ca_name2num = {}
for n, name in ca_num2name.items():
    ca_name2num[name] = n
    ca_name2num[name.upper()] = n
    ca_name2num[name.lower()] = n

hom_raw = pd.read_csv(os.path.join(DATA,
    "Violence_Reduction_-_Victims_of_Homicides_and_Non-Fatal_Shootings_20250513.csv"),
    low_memory=False)
hom_raw["DATE"] = pd.to_datetime(hom_raw["DATE"], format="%m/%d/%Y %I:%M:%S %p", errors="coerce")
hom_raw["CA_num"] = hom_raw["COMMUNITY_AREA"].map(ca_name2num)
hom_raw["year"] = hom_raw["DATE"].dt.year

ss_raw = pd.read_csv(os.path.join(DATA,
    "Violence_Reduction_-_Shotspotter_Alerts_-_Historical_20250408.csv"))
ss_raw["DATE"] = pd.to_datetime(ss_raw["DATE"], format="%m/%d/%Y %I:%M:%S %p", errors="coerce")
ss_raw = ss_raw[ss_raw["DATE"] <= "2023-12-31"].copy()
sample = ss_raw["COMMUNITY_AREA"].dropna().iloc[0]
if str(sample).replace(".", "").isdigit():
    ss_raw["CA_name"] = pd.to_numeric(ss_raw["COMMUNITY_AREA"], errors="coerce").astype("Int64").map(ca_num2name)
else:
    ss_raw["CA_name"] = ss_raw["COMMUNITY_AREA"].str.title()

ss_profile = (ss_raw.dropna(subset=["CA_name"]).groupby("CA_name")
              .agg(total=("DATE", "size"), earliest=("DATE", "min"), latest=("DATE", "max")))
ss_profile["months"] = ((ss_profile["latest"] - ss_profile["earliest"]).dt.days / 30.4).round().astype(int)
ALL_SS = sorted(ss_profile[(ss_profile.total >= 100) & (ss_profile.months >= 12)].index.tolist())


def build_panel(year_start, year_end, treat_set, post_year=None):
    gun = hom_raw[
        (hom_raw["VICTIMIZATION_PRIMARY"] == "HOMICIDE") &
        (hom_raw["GUNSHOT_INJURY_I"] == "YES") &
        (hom_raw["DATE"] >= f"{year_start}-01-01") &
        (hom_raw["DATE"] <= f"{year_end}-12-31")
    ]
    panel = pd.DataFrame([(ca, yr) for ca in sorted(ca_num2name.keys())
                          for yr in range(year_start, year_end + 1)],
                         columns=["ca_num", "year"])
    panel["ca_name"] = panel["ca_num"].map(ca_num2name)
    g = gun.dropna(subset=["CA_num"]).groupby(["CA_num", "year"]).size().reset_index(name="gun_hom")
    g["CA_num"] = g["CA_num"].astype(int)
    panel = panel.merge(g.rename(columns={"CA_num": "ca_num"}), on=["ca_num", "year"], how="left")
    panel["gun_hom"] = panel["gun_hom"].fillna(0).astype(int)
    panel["treat"] = panel["ca_name"].isin(treat_set).astype(int)
    panel["post"] = (panel["year"] >= (post_year or 2017)).astype(int)
    panel["did"] = panel["treat"] * panel["post"]
    return panel


# =====================================================================
# FIGURE 1: SPECIFICATION COMPARISON (OLS vs Poisson vs NegBin)
# =====================================================================
print("Building specification comparison plot...")

panel_main = build_panel(2009, 2023, ALL_SS)

# Estimate the three models
ols = smf.ols("gun_hom ~ did + C(ca_num) + C(year)", data=panel_main).fit(
    cov_type="cluster", cov_kwds={"groups": panel_main["ca_num"]})
ols_b = ols.params["did"]
ols_se = ols.bse["did"]
ols_ci = ols.conf_int().loc["did"].values
ols_p = ols.pvalues["did"]

# Translate OLS effect to IRR-equivalent: divide by treated-area pre-treatment mean
treated_pre_mean = panel_main[(panel_main["treat"] == 1) & (panel_main["post"] == 0)]["gun_hom"].mean()
ols_irr = (treated_pre_mean + ols_b) / treated_pre_mean
ols_irr_lo = (treated_pre_mean + ols_ci[0]) / treated_pre_mean
ols_irr_hi = (treated_pre_mean + ols_ci[1]) / treated_pre_mean

# Poisson
pois = smf.glm("gun_hom ~ did + C(ca_num) + C(year)", data=panel_main,
               family=sm.families.Poisson()).fit(
    cov_type="cluster", cov_kwds={"groups": panel_main["ca_num"]})
pois_b = pois.params["did"]
pois_se = pois.bse["did"]
pois_ci = pois.conf_int().loc["did"].values
pois_p = pois.pvalues["did"]
pois_irr = float(np.exp(pois_b))
pois_irr_lo = float(np.exp(pois_ci[0]))
pois_irr_hi = float(np.exp(pois_ci[1]))

# Negative Binomial — dispersion estimated (Cameron-Trivedi NB2 auxiliary regression),
# not fixed at 1, to match scripts/econometric_robustness.py.
_mu = pois.fittedvalues
_z = ((panel_main["gun_hom"].astype(float) - _mu) ** 2 - panel_main["gun_hom"]) / _mu
nb_alpha = max(float(sm.OLS(_z, _mu).fit().params.iloc[0]), 1e-6)
nb = smf.glm("gun_hom ~ did + C(ca_num) + C(year)", data=panel_main,
             family=sm.families.NegativeBinomial(alpha=nb_alpha)).fit(
    cov_type="cluster", cov_kwds={"groups": panel_main["ca_num"]})
nb_b = nb.params["did"]
nb_se = nb.bse["did"]
nb_ci = nb.conf_int().loc["did"].values
nb_p = nb.pvalues["did"]
nb_irr = float(np.exp(nb_b))
nb_irr_lo = float(np.exp(nb_ci[0]))
nb_irr_hi = float(np.exp(nb_ci[1]))


fig, ax = plt.subplots(figsize=(11, 7))

# Three rows
specs = [
    ("OLS DiD\n(implied multiplicative effect)", ols_irr, ols_irr_lo, ols_irr_hi, ols_p, C_RED, "headline as published"),
    ("Poisson DiD with two-way FE", pois_irr, pois_irr_lo, pois_irr_hi, pois_p, C_BLUE, "proper count-data model"),
    ("Negative Binomial DiD with two-way FE", nb_irr, nb_irr_lo, nb_irr_hi, nb_p, C_BLUE, "robust to overdispersion"),
]
y_pos = list(range(len(specs)))[::-1]

ax.axvline(1.0, color="#333333", lw=1, zorder=0)
ax.axvspan(0.95, 1.05, color=COLORS["grid"], alpha=0.6, zorder=0)

for y, (label, irr, lo, hi, p, color, note) in zip(y_pos, specs):
    ax.plot([lo, hi], [y, y], color=color, lw=3, zorder=2)
    ax.plot(irr, y, "o", color=color, markersize=12, zorder=3,
            markeredgecolor="white", markeredgewidth=1.2)
    ax.text(hi + 0.02, y + 0.18, f"IRR = {irr:.3f}  ({format_pvalue(p)}{sig_stars(p)})",
            va="center", fontsize=10, color=color, fontweight="bold")
    ax.text(hi + 0.02, y - 0.18, note, va="center", fontsize=9,
            color=COLORS["text_muted"], style="italic")

ax.set_yticks(y_pos)
ax.set_yticklabels([s[0] for s in specs], fontsize=11)
ax.set_xlabel("Incidence-rate ratio (IRR) on gun homicides\n"
              "Values > 1 mean ShotSpotter areas had relatively more homicides; values < 1 mean fewer",
              fontsize=10)
ax.set_title("ShotSpotter DiD: OLS Headline vs. Count-Data Specifications\n"
             "Same panel (51 treated areas, 26 controls, 2009–2023, two-way FE, clustered SEs)",
             fontweight="bold", fontsize=13, pad=14)
ax.set_xlim(0.6, 1.55)
ax.set_ylim(-0.5, len(specs) - 0.4)
ax.grid(True, axis="x", alpha=0.6, lw=0.5)
ax.set_axisbelow(True)
ax.spines["left"].set_visible(False)
ax.tick_params(axis="y", length=0)

# Takeaway moved to figure-level text below the axes (no longer competes for axis space)
_red = round((1 - pois_irr_lo) * 100)
_inc = round((pois_irr_hi - 1) * 100)
fig.text(0.5, 0.07,
         f"Takeaway: the OLS headline implies a large positive effect (IRR = {ols_irr:.2f}, p < 0.001),\n"
         f"but both count-data models are null (IRR ≈ {pois_irr:.2f}) — "
         f"a {_red}% reduction to a {_inc}% increase.",
         fontsize=10, color="#333333", style="italic", ha="center", va="bottom")

add_source_note(fig, SOURCE_DEFAULT,
                sample_n="N = 1,155 community-area-years (77 areas × 15 years)",
                y=0.005, x=0.5)
for txt in fig.texts:
    if txt.get_text().startswith("Source:"):
        txt.set_ha("center")

fig.subplots_adjust(top=0.88, bottom=0.30, left=0.30, right=0.97)
fig.savefig(os.path.join(PLOTS, "specification_comparison.png"), dpi=200,
            facecolor="white")
plt.close(fig)
print("  SAVED: specification_comparison.png")


# =====================================================================
# FIGURE 2: FALSIFICATION PLOT
# =====================================================================
print("\nBuilding falsification plot...")

# Run the same DiD with placebo treatment dates plus the actual one
spec_results = []
for fake_year in [1999, 2003, 2007, 2011]:
    p_panel = build_panel(fake_year - 8, fake_year + 4, ALL_SS, post_year=fake_year)
    m = smf.ols("gun_hom ~ did + C(ca_num) + C(year)", data=p_panel).fit(
        cov_type="cluster", cov_kwds={"groups": p_panel["ca_num"]})
    spec_results.append({"label": f"Placebo {fake_year}",
                         "fake": True,
                         "b": m.params["did"], "se": m.bse["did"], "p": m.pvalues["did"],
                         "ci": m.conf_int().loc["did"].values})

# Actual 2017
spec_results.append({"label": "Actual 2017", "fake": False,
                     "b": ols_b, "se": ols_se, "p": ols_p, "ci": ols_ci})


fig, ax = plt.subplots(figsize=(11, 7.5))
y_pos = list(range(len(spec_results)))[::-1]

ax.axvline(0, color="#333333", lw=1, zorder=0)

for y, r in zip(y_pos, spec_results):
    color = C_GRAY if r["fake"] else C_RED
    ax.plot([r["ci"][0], r["ci"][1]], [y, y], color=color, lw=3, zorder=2)
    ax.plot(r["b"], y, "o", color=color, markersize=11, zorder=3,
            markeredgecolor="white", markeredgewidth=1.2)
    # white backing box so the zero line never strikes through the label text
    ax.text(r["ci"][1] + 0.18, y,
            f"β = {r['b']:+.2f}  ({format_pvalue(r['p'])}{sig_stars(r['p'])})",
            va="center", fontsize=10, color=color, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.22", fc="white", alpha=0.9, ec="none"),
            zorder=4)

ax.set_yticks(y_pos)
ax.set_yticklabels([r["label"] for r in spec_results], fontsize=11)
ax.set_xlabel("OLS DiD coefficient on gun homicides per community-area-year\n"
              "(same specification: gun_hom ~ did + C(ca_num) + C(year), cluster-robust SEs)",
              fontsize=10)
ax.set_title("Falsification Test: Same Design Produces 'Significant' Effects\n"
             "in Pre-ShotSpotter Placebo Windows",
             fontweight="bold", fontsize=13, pad=14)
ax.grid(True, axis="x", alpha=0.6, lw=0.5)
ax.set_axisbelow(True)
ax.spines["left"].set_visible(False)
ax.tick_params(axis="y", length=0)
ax.set_xlim(-4.5, 7)

# Legend (top-left of axes, inside)
legend_elems = [
    mpatches.Patch(color=C_RED, label="Actual treatment (2017)"),
    mpatches.Patch(color=C_GRAY, label="Placebo (fake treatment date)"),
]
ax.legend(handles=legend_elems, loc="upper right", frameon=False, fontsize=10)

# Takeaway as figure-level text BELOW the axes — won't compress the data area
fig.text(0.5, 0.07,
         "Takeaway: if +2.6 reflected a real treatment effect, placebo years should be null. "
         "Two of four placebos are significant at p < 0.01,\n"
         "so the design is capturing pre-existing trend differences — the post-2017 "
         "estimate cannot be cleanly attributed to ShotSpotter.",
         fontsize=10, color="#333333", style="italic", ha="center", va="bottom")

add_source_note(fig, SOURCE_DEFAULT,
                sample_n="Each window: 51 SS areas + 26 controls × 13 years",
                y=0.005, x=0.5)
for txt in fig.texts:
    if txt.get_text().startswith("Source:"):
        txt.set_ha("center")

fig.subplots_adjust(top=0.88, bottom=0.26, left=0.16, right=0.97)
fig.savefig(os.path.join(PLOTS, "falsification_plot.png"), dpi=200,
            facecolor="white")
plt.close(fig)
print("  SAVED: falsification_plot.png")

print("\nDone. Two narrative plots written to plots/.")
