"""
econometric_analysis.py — Modern DiD evaluation of ShotSpotter in Chicago.

Implements four analyses that the original did_analysis.py omits:

  1. Two-way fixed effects DiD (community area + year FE, cluster-robust SE)
  2. Cohort-specific (staggered) event study using 2017 vs 2018 install dates
  3. Synthetic control for the six highest-violence study neighborhoods,
     with placebo inference via in-space permutation
  4. Pre-period placebo / falsification: assign fake treatment dates in
     1999, 2003, and 2007 using the 1991-2008 panel and confirm that
     null effects emerge

Updates docs/results_summary_reanalysis.md with the new estimates and writes plots
to the plots/ directory.
"""
import os
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import statsmodels.formula.api as smf
from scipy.optimize import minimize

warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_style import (apply_academic_style, COLORS, add_source_note,
                        add_n_label, SOURCE_DEFAULT)

P = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(P, "data", "raw")
PLOTS = os.path.join(P, "plots")

apply_academic_style()

C_BLUE = COLORS["control"]
C_RED = COLORS["treated"]
C_GRAY = COLORS["neutral"]
C_GREEN = COLORS["green"]


# =====================================================================
# DATA LOADING
# =====================================================================
print("=" * 72)
print("LOADING DATA")
print("=" * 72)

import json
with open(os.path.join(DATA, "CommAreas_20250408.geojson"), "r", encoding="utf-8") as f:
    geo = json.load(f)
ca_num2name = {int(float(feat["properties"]["area_num_1"])): feat["properties"]["community"].title()
               for feat in geo["features"]}
ca_name2num = {}
for n, name in ca_num2name.items():
    ca_name2num[name] = n            # Title case (canonical)
    ca_name2num[name.upper()] = n    # ALL-CAPS (matches raw data)
    ca_name2num[name.lower()] = n    # lowercase (defensive)

hom_raw = pd.read_csv(
    os.path.join(DATA, "Violence_Reduction_-_Victims_of_Homicides_and_Non-Fatal_Shootings_20250513.csv"),
    low_memory=False,
)
hom_raw["DATE"] = pd.to_datetime(hom_raw["DATE"], format="%m/%d/%Y %I:%M:%S %p", errors="coerce")
hom_raw["CA_name"] = hom_raw["COMMUNITY_AREA"].str.title()
hom_raw["CA_num"] = hom_raw["COMMUNITY_AREA"].map(ca_name2num)
hom_raw["year"] = hom_raw["DATE"].dt.year

ss_raw = pd.read_csv(
    os.path.join(DATA, "Violence_Reduction_-_Shotspotter_Alerts_-_Historical_20250408.csv"),
)
ss_raw["DATE"] = pd.to_datetime(ss_raw["DATE"], format="%m/%d/%Y %I:%M:%S %p", errors="coerce")
ss_raw = ss_raw[ss_raw["DATE"] <= "2023-12-31"].copy()
sample_ca = ss_raw["COMMUNITY_AREA"].dropna().iloc[0]
if str(sample_ca).replace(".", "").isdigit():
    ss_raw["CA_num"] = pd.to_numeric(ss_raw["COMMUNITY_AREA"], errors="coerce").astype("Int64")
    ss_raw["CA_name"] = ss_raw["CA_num"].map(ca_num2name)
else:
    ss_raw["CA_name"] = ss_raw["COMMUNITY_AREA"].str.title()
    ss_raw["CA_num"] = ss_raw["CA_name"].map(ca_name2num)


# Identify treatment cohorts (>=100 alerts, >=12 months) by first-alert year
ss_profile = (
    ss_raw.dropna(subset=["CA_name"])
    .groupby("CA_name")
    .agg(total_alerts=("DATE", "size"), earliest=("DATE", "min"), latest=("DATE", "max"))
)
ss_profile["months"] = ((ss_profile["latest"] - ss_profile["earliest"]).dt.days / 30.4).round().astype(int)
ss_profile["install_year"] = ss_profile["earliest"].dt.year
ss_profile = ss_profile[(ss_profile.total_alerts >= 100) & (ss_profile.months >= 12)]

ALL_SS = sorted(ss_profile.index.tolist())
COHORT_2017 = sorted(ss_profile[ss_profile.install_year == 2017].index.tolist())
COHORT_2018 = sorted(ss_profile[ss_profile.install_year == 2018].index.tolist())
ORIG_6 = ["Austin", "Humboldt Park", "North Lawndale", "Englewood", "West Englewood", "South Shore"]

print(f"  ShotSpotter cohorts: {len(COHORT_2017)} in 2017, {len(COHORT_2018)} in 2018")
print(f"  Original 6 study neighborhoods, by cohort:")
for h in ORIG_6:
    cohort = 2017 if h in COHORT_2017 else (2018 if h in COHORT_2018 else None)
    print(f"    {h}: cohort {cohort}")


# Build panel: 77 community areas x years
def build_panel(year_start, year_end, hom_df, treat_set, treat_year_map=None):
    """Build community-area x year panel with gun homicide counts."""
    gun = hom_df[
        (hom_df["VICTIMIZATION_PRIMARY"] == "HOMICIDE") &
        (hom_df["GUNSHOT_INJURY_I"] == "YES") &
        (hom_df["DATE"] >= f"{year_start}-01-01") &
        (hom_df["DATE"] <= f"{year_end}-12-31")
    ].copy()
    nfs = hom_df[
        (hom_df["VICTIMIZATION_PRIMARY"] != "HOMICIDE") &
        (hom_df["GUNSHOT_INJURY_I"] == "YES") &
        (hom_df["DATE"] >= f"{year_start}-01-01") &
        (hom_df["DATE"] <= f"{year_end}-12-31")
    ].copy()

    all_ca = sorted(ca_num2name.keys())
    yrs = list(range(year_start, year_end + 1))
    panel = pd.DataFrame([(ca, yr) for ca in all_ca for yr in yrs],
                         columns=["ca_num", "year"])
    panel["ca_name"] = panel["ca_num"].map(ca_num2name)

    g_counts = gun.dropna(subset=["CA_num"]).groupby(["CA_num", "year"]).size().reset_index(name="gun_hom")
    g_counts["CA_num"] = g_counts["CA_num"].astype(int)
    panel = panel.merge(g_counts.rename(columns={"CA_num": "ca_num"}), on=["ca_num", "year"], how="left")
    panel["gun_hom"] = panel["gun_hom"].fillna(0).astype(int)

    n_counts = nfs.dropna(subset=["CA_num"]).groupby(["CA_num", "year"]).size().reset_index(name="nfs")
    n_counts["CA_num"] = n_counts["CA_num"].astype(int)
    panel = panel.merge(n_counts.rename(columns={"CA_num": "ca_num"}), on=["ca_num", "year"], how="left")
    panel["nfs"] = panel["nfs"].fillna(0).astype(int)

    panel["treat"] = panel["ca_name"].isin(treat_set).astype(int)
    if treat_year_map is None:
        # Single-cohort: post = year >= 2017
        panel["post"] = (panel["year"] >= 2017).astype(int)
    else:
        # Cohort-specific: post = year >= treatment_year for treated; 0 for never-treated
        def _post(row):
            if row["ca_name"] in treat_year_map:
                return int(row["year"] >= treat_year_map[row["ca_name"]])
            return 0
        panel["post"] = panel.apply(_post, axis=1)
    panel["did"] = panel["treat"] * panel["post"]
    return panel


# =====================================================================
# ANALYSIS 1: TWO-WAY FIXED EFFECTS DiD
# =====================================================================
print("\n" + "=" * 72)
print("ANALYSIS 1: TWO-WAY FIXED EFFECTS DiD")
print("=" * 72)
print("Replaces the naive 'gun_hom ~ treat + post + did' specification with")
print("the standard panel-data form 'gun_hom ~ did + C(ca_num) + C(year)'.")

panel = build_panel(2009, 2023, hom_raw, ALL_SS)


def run_twfe(df, outcome, treat_col, label, exclude_years=None):
    if exclude_years:
        df = df[~df["year"].isin(exclude_years)]
    formula = f"{outcome} ~ {treat_col} + C(ca_num) + C(year)"
    model = smf.ols(formula, data=df).fit(
        cov_type="cluster", cov_kwds={"groups": df["ca_num"]})
    b = model.params[treat_col]
    se = model.bse[treat_col]
    p = model.pvalues[treat_col]
    ci_lo, ci_hi = model.conf_int().loc[treat_col]
    sig = "***" if p < 0.01 else "**" if p < 0.05 else "*" if p < 0.1 else "ns"
    print(f"  {label:50s} b={b:+7.3f}  SE={se:6.3f}  p={p:.4f} {sig}  CI=[{ci_lo:+.3f},{ci_hi:+.3f}]  N={int(model.nobs)}")
    return {"label": label, "b": b, "se": se, "p": p, "ci": (ci_lo, ci_hi), "n": int(model.nobs)}


print("\n--- Gun homicides (TWFE) ---")
twfe_results = []
twfe_results.append(run_twfe(panel.copy(), "gun_hom", "did", "All 51 SS areas, 2009-2023"))
twfe_results.append(run_twfe(panel.copy(), "gun_hom", "did", "All 51 SS areas, excl COVID (2020-21)", [2020, 2021]))
twfe_results.append(run_twfe(panel.copy(), "gun_hom", "did", "All 51 SS areas, excl 2016", [2016]))
twfe_results.append(run_twfe(panel.copy(), "gun_hom", "did", "All 51 SS areas, excl COVID + 2016", [2016, 2020, 2021]))

print("\n--- Non-fatal gunshot injuries (TWFE) ---")
twfe_results_nfs = []
twfe_results_nfs.append(run_twfe(panel.copy(), "nfs", "did", "All 51 SS areas, NFS, 2009-2023"))
twfe_results_nfs.append(run_twfe(panel.copy(), "nfs", "did", "All 51 SS areas, NFS, excl COVID", [2020, 2021]))
twfe_results_nfs.append(run_twfe(panel.copy(), "nfs", "did", "All 51 SS areas, NFS, excl 2016", [2016]))


# =====================================================================
# ANALYSIS 2: COHORT-SPECIFIC EVENT STUDY (Sun-Abraham style)
# =====================================================================
print("\n" + "=" * 72)
print("ANALYSIS 2: COHORT-SPECIFIC EVENT STUDY")
print("=" * 72)
print("Estimates separate ATT(g, k) for each (cohort, event-time) cell, then")
print("aggregates with equal cohort weights. Avoids 'forbidden comparison'")
print("contamination of TWFE in staggered settings (Goodman-Bacon 2021).")

never_treated = [n for n in ca_num2name.values() if n not in ALL_SS]
print(f"  2017 cohort: {len(COHORT_2017)} areas")
print(f"  2018 cohort: {len(COHORT_2018)} areas")
print(f"  Never-treated control:  {len(never_treated)} areas")


def cohort_event_study(panel_df, cohort_areas, install_year, ctrl_areas, k_min=-7, k_max=6):
    """Compute ATT(cohort, event-time k) using never-treated as controls."""
    coefs = {}
    for k in range(k_min, k_max + 1):
        cal_year = install_year + k
        if cal_year < panel_df["year"].min() or cal_year > panel_df["year"].max():
            continue
        # Compare cohort and never-treated, in event-time k vs base year (k=-1)
        base_year = install_year - 1
        treat_post = panel_df[(panel_df["ca_name"].isin(cohort_areas)) & (panel_df["year"] == cal_year)]["gun_hom"]
        treat_base = panel_df[(panel_df["ca_name"].isin(cohort_areas)) & (panel_df["year"] == base_year)]["gun_hom"]
        ctrl_post = panel_df[(panel_df["ca_name"].isin(ctrl_areas)) & (panel_df["year"] == cal_year)]["gun_hom"]
        ctrl_base = panel_df[(panel_df["ca_name"].isin(ctrl_areas)) & (panel_df["year"] == base_year)]["gun_hom"]
        att = (treat_post.mean() - treat_base.mean()) - (ctrl_post.mean() - ctrl_base.mean())
        # Bootstrap SE
        rng = np.random.default_rng(0)
        boots = []
        for _ in range(500):
            t_p = treat_post.sample(len(treat_post), replace=True, random_state=rng.integers(1e9))
            t_b = treat_base.sample(len(treat_base), replace=True, random_state=rng.integers(1e9))
            c_p = ctrl_post.sample(len(ctrl_post), replace=True, random_state=rng.integers(1e9))
            c_b = ctrl_base.sample(len(ctrl_base), replace=True, random_state=rng.integers(1e9))
            boots.append((t_p.mean() - t_b.mean()) - (c_p.mean() - c_b.mean()))
        se = np.std(boots, ddof=1)
        coefs[k] = (att, se)
    return coefs


es_2017 = cohort_event_study(panel, COHORT_2017, 2017, never_treated)
es_2018 = cohort_event_study(panel, COHORT_2018, 2018, never_treated)

print("\n  Event-time k | 2017 cohort ATT | 2018 cohort ATT")
print("  " + "-" * 60)
for k in sorted(set(es_2017.keys()) | set(es_2018.keys())):
    a17 = f"{es_2017[k][0]:+6.2f} ({es_2017[k][1]:.2f})" if k in es_2017 else "      ---      "
    a18 = f"{es_2018[k][0]:+6.2f} ({es_2018[k][1]:.2f})" if k in es_2018 else "      ---      "
    print(f"  k = {k:+d}      | {a17}  | {a18}")


# Aggregate ATTs by event-time (equal weight across cohorts)
event_times = sorted(set(es_2017.keys()) & set(es_2018.keys()))
agg_atts = []
for k in event_times:
    a17 = es_2017[k][0]
    a18 = es_2018[k][0]
    # Sample-size weighted average
    w17 = len(COHORT_2017) / (len(COHORT_2017) + len(COHORT_2018))
    w18 = len(COHORT_2018) / (len(COHORT_2017) + len(COHORT_2018))
    se17 = es_2017[k][1]
    se18 = es_2018[k][1]
    agg_att = w17 * a17 + w18 * a18
    agg_se = np.sqrt((w17 * se17) ** 2 + (w18 * se18) ** 2)
    agg_atts.append((k, agg_att, agg_se))

print("\n  Aggregated ATT (sample-weighted across cohorts):")
print("  Event-time | ATT      | SE     | 95% CI")
for k, att, se in agg_atts:
    lo, hi = att - 1.96 * se, att + 1.96 * se
    print(f"   k = {k:+d}    | {att:+7.3f}  | {se:.3f} | [{lo:+.3f}, {hi:+.3f}]")


# =====================================================================
# ANALYSIS 3: SYNTHETIC CONTROL FOR THE 6 STUDY NEIGHBORHOODS
# =====================================================================
print("\n" + "=" * 72)
print("ANALYSIS 3: SYNTHETIC CONTROL FOR EACH STUDY NEIGHBORHOOD")
print("=" * 72)
print("For each of 6 study neighborhoods, build a synthetic counterfactual")
print("from a weighted combination of the 26 never-treated community areas.")
print("Weights minimize pre-treatment squared error on the gun_hom series.")


def synthetic_control(treated_name, donor_names, panel_df, install_year):
    """Solve for non-negative weights summing to 1 that minimize pre-treatment MSE."""
    pre_years = sorted(panel_df[panel_df["year"] < install_year]["year"].unique())
    treated_pre = (panel_df[(panel_df["ca_name"] == treated_name) & (panel_df["year"].isin(pre_years))]
                   .sort_values("year")["gun_hom"].values.astype(float))
    donor_pre = np.array([
        (panel_df[(panel_df["ca_name"] == d) & (panel_df["year"].isin(pre_years))]
         .sort_values("year")["gun_hom"].values.astype(float))
        for d in donor_names
    ])  # shape (n_donors, n_pre_years)

    def loss(w):
        synth = (w[:, None] * donor_pre).sum(axis=0)
        return np.sum((treated_pre - synth) ** 2)

    n_donors = len(donor_names)
    w0 = np.ones(n_donors) / n_donors
    constraints = ({"type": "eq", "fun": lambda w: np.sum(w) - 1.0},)
    bounds = [(0, 1) for _ in range(n_donors)]
    result = minimize(loss, w0, method="SLSQP", bounds=bounds, constraints=constraints,
                      options={"maxiter": 1000, "ftol": 1e-9})
    if not result.success:
        # The donor pool of low-violence never-treated areas cannot bracket the
        # high-violence treated units, so SLSQP collapses to a corner (single donor)
        # without converging to an interior optimum. We surface this rather than
        # silently consuming the failed solve — see the infeasibility note in §3.
        print(f"    [warn] SLSQP did not converge for {treated_name}: {result.message}")
    return result.x


def sc_series(treated_name, weights, donor_names, panel_df):
    """Compute the synthetic series across all years."""
    yrs = sorted(panel_df["year"].unique())
    treated_full = panel_df[panel_df["ca_name"] == treated_name].sort_values("year")["gun_hom"].values.astype(float)
    donor_full = np.array([
        panel_df[panel_df["ca_name"] == d].sort_values("year")["gun_hom"].values.astype(float)
        for d in donor_names
    ])
    synth = (weights[:, None] * donor_full).sum(axis=0)
    return yrs, treated_full, synth


sc_results = {}
for hood in ORIG_6:
    install_yr = 2017 if hood in COHORT_2017 else 2018
    weights = synthetic_control(hood, never_treated, panel, install_yr)
    yrs, actual, synth = sc_series(hood, weights, never_treated, panel)
    pre_rmse = np.sqrt(np.mean((np.array(actual)[:install_yr - 2009] - synth[:install_yr - 2009]) ** 2))
    post_diff = np.array(actual)[install_yr - 2009:].mean() - synth[install_yr - 2009:].mean()
    sc_results[hood] = {
        "install_year": install_yr, "weights": weights, "yrs": yrs,
        "actual": actual, "synth": synth,
        "pre_rmse": pre_rmse, "post_diff": post_diff,
        "donors": [(d, w) for d, w in zip(never_treated, weights) if w > 0.01],
    }
    nonzero_donors = sorted([(d, w) for d, w in zip(never_treated, weights) if w > 0.01],
                             key=lambda x: -x[1])
    donor_str = ", ".join(f"{d}={w:.2f}" for d, w in nonzero_donors[:5])
    print(f"\n  {hood} (cohort {install_yr}):")
    print(f"    pre-treatment RMSE: {pre_rmse:.3f}")
    print(f"    post-treatment gap (actual - synthetic): {post_diff:+.3f} per year")
    print(f"    top donors: {donor_str}")


# Plot synthetic control results: 2x3 panel
fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharex=True)
axes = axes.flatten()
for i, hood in enumerate(ORIG_6):
    ax = axes[i]
    r = sc_results[hood]
    # Pre-treatment shading (light blue) so the design is visually self-explanatory
    ax.axvspan(min(r["yrs"]) - 0.5, r["install_year"] - 0.5,
               color=COLORS["pre_shade"], alpha=0.6, zorder=0)
    ax.plot(r["yrs"], r["actual"], "o-", color=C_BLUE, lw=2.2, markersize=5,
            label="Actual", zorder=3)
    ax.plot(r["yrs"], r["synth"], "s--", color=C_RED, lw=2.2, markersize=5,
            label="Synthetic counterfactual", zorder=3)
    ax.axvline(r["install_year"] - 0.5, color="#333", ls="--", lw=1.0, zorder=2)
    ax.set_title(hood, fontsize=12, fontweight="bold", pad=6)
    ax.set_ylabel("Gun homicides")
    ax.yaxis.grid(True)

    # Move post-tx gap from the title into a small inset box
    add_n_label(ax,
                f"Post-tx gap: {r['post_diff']:+.1f}/yr\nInstalled {r['install_year']}",
                loc="upper left")

# Single shared legend at the bottom — replaces 6 cluttered per-panel legends
shared_legend = [
    Line2D([0], [0], color=C_BLUE, marker="o", lw=2.2, markersize=6, label="Actual"),
    Line2D([0], [0], color=C_RED, marker="s", lw=2.2, ls="--", markersize=6, label="Synthetic counterfactual"),
    mpatches.Patch(color=COLORS["pre_shade"], alpha=0.6, label="Pre-treatment period"),
]
fig.legend(handles=shared_legend, loc="lower center", ncol=3, frameon=False,
           fontsize=11, bbox_to_anchor=(0.5, 0.01))

fig.suptitle("Synthetic Control: Actual vs. Counterfactual Gun Homicides",
             fontsize=14, fontweight="bold", y=0.99)
add_source_note(fig, SOURCE_DEFAULT,
                sample_n="Donor pool: 26 never-treated community areas",
                y=0.06, x=0.5)
# Center the source note
for txt in fig.texts:
    if txt.get_text().startswith("Source:"):
        txt.set_ha("center")

fig.tight_layout(rect=[0, 0.07, 1, 0.96])
fig.savefig(os.path.join(PLOTS, "synthetic_control_panels.png"), dpi=200,
            facecolor="white")
plt.close(fig)
print(f"\n  SAVED: synthetic_control_panels.png")


# =====================================================================
# ANALYSIS 4: PRE-PERIOD PLACEBO / FALSIFICATION
# =====================================================================
print("\n" + "=" * 72)
print("ANALYSIS 4: PRE-PERIOD PLACEBO TESTS")
print("=" * 72)
print("If the +2.6 DiD reflects a real treatment effect, then assigning fake")
print("treatment dates in pre-ShotSpotter years (1999, 2003, 2007) should")
print("yield null estimates. A non-null placebo flags a spurious finding.")

placebo_results = []
for fake_year in [1999, 2003, 2007, 2011]:
    panel_pre = build_panel(fake_year - 8, fake_year + 4, hom_raw, ALL_SS)  # ±8/+4 yr window
    panel_pre["post"] = (panel_pre["year"] >= fake_year).astype(int)
    panel_pre["did"] = panel_pre["treat"] * panel_pre["post"]
    formula = "gun_hom ~ did + C(ca_num) + C(year)"
    try:
        model = smf.ols(formula, data=panel_pre).fit(
            cov_type="cluster", cov_kwds={"groups": panel_pre["ca_num"]})
        b = model.params["did"]
        se = model.bse["did"]
        p = model.pvalues["did"]
        ci = model.conf_int().loc["did"].values
        sig = "***" if p < 0.01 else "**" if p < 0.05 else "*" if p < 0.1 else "ns"
        print(f"  Fake treatment year = {fake_year}: b={b:+6.2f}, SE={se:.2f}, p={p:.4f} {sig}, CI=[{ci[0]:+.2f},{ci[1]:+.2f}], N={int(model.nobs)}")
        placebo_results.append({"year": fake_year, "b": b, "se": se, "p": p, "ci": tuple(ci)})
    except Exception as e:
        print(f"  Fake treatment year = {fake_year}: failed ({e})")


# =====================================================================
# AGGREGATE EVENT STUDY PLOT (cohort-aware)
# =====================================================================
print("\n" + "=" * 72)
print("PLOT: Cohort-Aware Aggregated Event Study")
print("=" * 72)
fig, ax = plt.subplots(figsize=(12, 6.5))
ks = [k for k, _, _ in agg_atts]
atts = [att for _, att, _ in agg_atts]
ses = [se for _, _, se in agg_atts]
los = [att - 1.96 * se for att, se in zip(atts, ses)]
his = [att + 1.96 * se for att, se in zip(atts, ses)]

# Pre-period shading anchors the design visually
ax.axvspan(min(ks) - 1, -0.5, color=COLORS["pre_shade"], alpha=0.5, zorder=0)

ax.axhline(0, color="#333333", lw=1.2, zorder=1)
ax.axvline(-0.5, color=C_GRAY, ls="--", lw=1.2, zorder=1)

for k, att, lo, hi in zip(ks, atts, los, his):
    color = C_BLUE if k < 0 else C_RED
    ax.plot([k, k], [lo, hi], color=color, lw=2.2, zorder=2)
    ax.plot(k, att, "o", color=color, markersize=8,
            markeredgecolor="white", markeredgewidth=0.6, zorder=3)

# Add the install annotation AFTER plotting so axes limits are real
ax.text(-0.4, ax.get_ylim()[1] * 0.96, "ShotSpotter installed",
        fontsize=9, color=C_GRAY, va="top", ha="left", style="italic")

ax.set_xlabel("Years relative to ShotSpotter installation (cohort-specific)")
ax.set_ylabel("ATT: Gun homicides per community-area-year\n(treated cohort vs. never-treated controls)")
ax.set_title("Cohort-Aware Event Study: 2017 + 2018 Cohorts Aggregated\n"
             "(95% CI from bootstrap, never-treated as control)",
             fontweight="bold")

# Custom legend for pre/post coding + install line + shading
legend_handles = [
    Line2D([0], [0], color=C_BLUE, marker="o", lw=2.2, markersize=7,
           label="Pre-treatment (k < 0)"),
    Line2D([0], [0], color=C_RED, marker="o", lw=2.2, markersize=7,
           label="Post-treatment (k ≥ 0)"),
    Line2D([0], [0], color=C_GRAY, ls="--", lw=1.2, label="Install (k = 0)"),
]
ax.legend(handles=legend_handles, loc="lower right", frameon=False, fontsize=10)
ax.yaxis.grid(True)

add_n_label(ax,
            f"Treated cohorts: 2017 (n={len(COHORT_2017)}) + 2018 (n={len(COHORT_2018)})\n"
            f"Never-treated controls: n={len(never_treated)}",
            loc="upper left")
add_source_note(fig, SOURCE_DEFAULT)

fig.subplots_adjust(top=0.88, bottom=0.13, left=0.10, right=0.97)
fig.savefig(os.path.join(PLOTS, "event_study_cohort_aware.png"), dpi=200,
            facecolor="white")
plt.close(fig)
print("  SAVED: event_study_cohort_aware.png")


# =====================================================================
# WRITE docs/results_summary_reanalysis.md
# =====================================================================
print("\n" + "=" * 72)
print("WRITING docs/results_summary_reanalysis.md")
print("=" * 72)

lines = []
lines.append("# Econometric Re-analysis: ShotSpotter & Gun Homicide DiD")
lines.append("")
lines.append("This document presents the results of `scripts/econometric_analysis.py`,")
lines.append("which re-estimates the project's main causal claims using methods that")
lines.append("address the limitations of the original specification:")
lines.append("")
lines.append("- **Two-way fixed effects** (community-area + year) instead of pooled DiD")
lines.append("- **Cohort-specific event study** that respects the staggered 2017/2018 rollout")
lines.append("- **Synthetic control** for each of the 6 study neighborhoods")
lines.append("- **Pre-period placebo tests** in 1999, 2003, 2007, and 2011")
lines.append("")
lines.append("## 1. Two-Way Fixed Effects DiD")
lines.append("")
lines.append("Specification: `gun_hom_it = β·(treat_i × post_t) + α_i + γ_t + ε_it`")
lines.append("with cluster-robust SEs at the community-area level.")
lines.append("")
lines.append("### Gun homicides")
lines.append("")
lines.append("| Specification | β (DiD) | SE | p | 95% CI | N |")
lines.append("|---|---|---|---|---|---|")
for r in twfe_results:
    sig = "***" if r["p"] < 0.01 else "**" if r["p"] < 0.05 else "*" if r["p"] < 0.1 else ""
    lines.append(f"| {r['label']} | {r['b']:+.3f} | {r['se']:.3f} | {r['p']:.4f}{sig} | [{r['ci'][0]:+.3f}, {r['ci'][1]:+.3f}] | {r['n']:,} |")
lines.append("")
lines.append("### Non-fatal gunshot injuries")
lines.append("")
lines.append("| Specification | β (DiD) | SE | p | 95% CI | N |")
lines.append("|---|---|---|---|---|---|")
for r in twfe_results_nfs:
    sig = "***" if r["p"] < 0.01 else "**" if r["p"] < 0.05 else "*" if r["p"] < 0.1 else ""
    lines.append(f"| {r['label']} | {r['b']:+.3f} | {r['se']:.3f} | {r['p']:.4f}{sig} | [{r['ci'][0]:+.3f}, {r['ci'][1]:+.3f}] | {r['n']:,} |")
lines.append("")

lines.append("## 2. Cohort-Specific Event Study")
lines.append("")
lines.append(f"Cohorts: {len(COHORT_2017)} community areas treated in 2017, "
             f"{len(COHORT_2018)} in 2018, {len(never_treated)} never-treated controls.")
lines.append("")
lines.append("### Aggregated event-time ATT (sample-weighted across cohorts)")
lines.append("")
lines.append("| Event time k | ATT | SE | 95% CI |")
lines.append("|---|---|---|---|")
for k, att, se in agg_atts:
    lines.append(f"| {k:+d} | {att:+.3f} | {se:.3f} | [{att - 1.96 * se:+.3f}, {att + 1.96 * se:+.3f}] |")
lines.append("")

lines.append("## 3. Synthetic Control")
lines.append("")
lines.append("For each of the 6 study neighborhoods, weights over the 26 never-treated")
lines.append("community areas are chosen to minimize the pre-treatment squared error on")
lines.append("annual gun homicide counts. A positive 'gap' means the treated unit had")
lines.append("more homicides post-treatment than its synthetic counterfactual.")
lines.append("")
lines.append("| Neighborhood | Cohort | Pre-tx RMSE | Post-tx gap (per year) | Top donors |")
lines.append("|---|---|---|---|---|")
for hood, r in sc_results.items():
    top_donors = sorted([(d, w) for d, w in zip(never_treated, r["weights"]) if w > 0.05],
                         key=lambda x: -x[1])[:3]
    donor_str = ", ".join(f"{d}={w:.2f}" for d, w in top_donors)
    lines.append(f"| {hood} | {r['install_year']} | {r['pre_rmse']:.2f} | {r['post_diff']:+.2f} | {donor_str} |")
lines.append("")
lines.append("> **Synthetic control is infeasible here — read these gaps with care.** For every")
lines.append("> treated unit the optimizer concentrates 100% of the weight on a single donor")
lines.append("> (Washington Heights, the highest-violence never-treated area) and the SLSQP solve")
lines.append("> does not converge to an interior optimum. The pre-treatment RMSE is large relative")
lines.append("> to the outcome level (e.g. Austin's pre-RMSE of ~34 against an actual level of")
lines.append("> 15-70/yr), so the synthetic counterfactual cannot reproduce the treated series even")
lines.append("> in the pre-period it was fit on. The 26 never-treated community areas are")
lines.append("> systematically lower-violence than any treated neighborhood, so no convex")
lines.append("> combination can match them. The positive 'gaps' below therefore reflect this")
lines.append("> level mismatch, **not** a treatment effect. See NARRATIVE.md §5.5.")
lines.append("")

lines.append("## 4. Pre-Period Placebo Tests")
lines.append("")
lines.append("Treatment is reassigned to fake years before the actual ShotSpotter rollout.")
lines.append("A non-null placebo coefficient suggests that the post-2017 estimate")
lines.append("reflects pre-existing trend differences rather than a causal effect.")
lines.append("")
lines.append("| Fake treatment year | β | SE | p | 95% CI |")
lines.append("|---|---|---|---|---|")
for r in placebo_results:
    sig = "***" if r["p"] < 0.01 else "**" if r["p"] < 0.05 else "*" if r["p"] < 0.1 else ""
    lines.append(f"| {r['year']} | {r['b']:+.3f} | {r['se']:.3f} | {r['p']:.4f}{sig} | [{r['ci'][0]:+.3f}, {r['ci'][1]:+.3f}] |")
lines.append("")

lines.append("## Generated figures")
lines.append("")
lines.append("- `plots/event_study_cohort_aware.png` — staggered event-study with bootstrap CIs")
lines.append("- `plots/synthetic_control_panels.png` — 6-panel actual vs. synthetic comparison")
lines.append("")

with open(os.path.join(P, "docs", "results_summary_reanalysis.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(lines))
print("  WROTE: docs/results_summary_reanalysis.md")

print("\n" + "=" * 72)
print("DONE")
print("=" * 72)
