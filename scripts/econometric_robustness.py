"""
econometric_robustness.py — Final-layer robustness on top of econometric_analysis.py.

Adds three improvements that the prior script left on the table:

  1. Count-data regression (Poisson and Negative Binomial with two-way FE)
     — addresses the OLS misspecification of homicide counts as continuous.
  2. Joint Wald F-test for pre-trends — turns the "pre-period coefficients
     look bad" qualitative reading into a single test statistic.
  3. Permutation-based inference for synthetic control — runs "placebo"
     synthetic controls on each of the 26 never-treated community areas and
     compares the post-treatment gap of each treated unit to the placebo
     distribution.

Appends the new tables to RESULTS_SUMMARY_V2.md.
"""
import os
import sys
import json
import warnings
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_style import (apply_academic_style, COLORS, format_pvalue,
                        add_source_note, SOURCE_DEFAULT)

P = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(P, "data", "raw")
PLOTS = os.path.join(P, "plots")

apply_academic_style()

C_BLUE = COLORS["control"]
C_RED = COLORS["treated"]
C_GRAY = COLORS["neutral"]

# =====================================================================
# DATA (replicates loading logic from econometric_analysis.py)
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
    low_memory=False,
)
hom_raw["DATE"] = pd.to_datetime(hom_raw["DATE"], format="%m/%d/%Y %I:%M:%S %p", errors="coerce")
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

ss_profile = (
    ss_raw.dropna(subset=["CA_name"]).groupby("CA_name")
    .agg(total_alerts=("DATE", "size"), earliest=("DATE", "min"), latest=("DATE", "max"))
)
ss_profile["months"] = ((ss_profile["latest"] - ss_profile["earliest"]).dt.days / 30.4).round().astype(int)
ss_profile["install_year"] = ss_profile["earliest"].dt.year
ss_profile = ss_profile[(ss_profile.total_alerts >= 100) & (ss_profile.months >= 12)]
ALL_SS = sorted(ss_profile.index.tolist())
COHORT_2017 = sorted(ss_profile[ss_profile.install_year == 2017].index.tolist())
COHORT_2018 = sorted(ss_profile[ss_profile.install_year == 2018].index.tolist())
ORIG_6 = ["Austin", "Humboldt Park", "North Lawndale", "Englewood", "West Englewood", "South Shore"]
NEVER_TREATED = [n for n in ca_num2name.values() if n not in ALL_SS]


def build_panel(year_start, year_end, treat_set):
    gun = hom_raw[
        (hom_raw["VICTIMIZATION_PRIMARY"] == "HOMICIDE") &
        (hom_raw["GUNSHOT_INJURY_I"] == "YES") &
        (hom_raw["DATE"] >= f"{year_start}-01-01") &
        (hom_raw["DATE"] <= f"{year_end}-12-31")
    ].copy()
    all_ca = sorted(ca_num2name.keys())
    yrs = list(range(year_start, year_end + 1))
    panel = pd.DataFrame([(ca, yr) for ca in all_ca for yr in yrs],
                         columns=["ca_num", "year"])
    panel["ca_name"] = panel["ca_num"].map(ca_num2name)
    g_counts = gun.dropna(subset=["CA_num"]).groupby(["CA_num", "year"]).size().reset_index(name="gun_hom")
    g_counts["CA_num"] = g_counts["CA_num"].astype(int)
    panel = panel.merge(g_counts.rename(columns={"CA_num": "ca_num"}),
                        on=["ca_num", "year"], how="left")
    panel["gun_hom"] = panel["gun_hom"].fillna(0).astype(int)
    panel["treat"] = panel["ca_name"].isin(treat_set).astype(int)
    panel["post"] = (panel["year"] >= 2017).astype(int)
    panel["did"] = panel["treat"] * panel["post"]
    return panel


panel = build_panel(2009, 2023, ALL_SS)
print(f"  Panel built: {len(panel)} rows, {panel['ca_num'].nunique()} areas, {panel['year'].nunique()} years")


# =====================================================================
# 1. POISSON AND NEGATIVE BINOMIAL REGRESSION WITH TWO-WAY FE
# =====================================================================
print("\n" + "=" * 72)
print("1. COUNT-DATA REGRESSION (POISSON & NEGATIVE BINOMIAL)")
print("=" * 72)
print("Replaces OLS for the count outcome `gun_hom`. Coefficients are on")
print("the log-mean scale; exp(beta) gives the multiplicative effect.")


def run_count_model(family, label, df=panel):
    formula = "gun_hom ~ did + C(ca_num) + C(year)"
    try:
        model = smf.glm(formula, data=df, family=family).fit(
            cov_type="cluster", cov_kwds={"groups": df["ca_num"]})
        b = model.params["did"]
        se = model.bse["did"]
        p = model.pvalues["did"]
        ci = model.conf_int().loc["did"].values
        irr = np.exp(b)
        irr_lo, irr_hi = np.exp(ci)
        sig = "***" if p < 0.01 else "**" if p < 0.05 else "*" if p < 0.1 else ""
        print(f"  {label:35s}  beta={b:+.4f}  SE={se:.4f}  p={p:.4f}{sig}")
        print(f"  {'':35s}  IRR=exp(beta)={irr:.4f}  CI=[{irr_lo:.4f}, {irr_hi:.4f}]")
        return {"label": label, "b": b, "se": se, "p": p, "ci": tuple(ci),
                "irr": irr, "irr_ci": (irr_lo, irr_hi)}
    except Exception as e:
        print(f"  {label}: FAILED ({type(e).__name__}: {e})")
        return None


print("\n--- Poisson ---")
poisson_res = []
poisson_res.append(run_count_model(sm.families.Poisson(), "Full panel 2009-2023"))
poisson_res.append(run_count_model(sm.families.Poisson(), "Excl COVID (2020-21)",
                                    panel[~panel["year"].isin([2020, 2021])]))
poisson_res.append(run_count_model(sm.families.Poisson(), "Excl 2016",
                                    panel[panel["year"] != 2016]))

print("\n--- Negative Binomial (alpha=1) ---")
nb_res = []
nb_res.append(run_count_model(sm.families.NegativeBinomial(alpha=1.0), "Full panel 2009-2023"))
nb_res.append(run_count_model(sm.families.NegativeBinomial(alpha=1.0), "Excl COVID (2020-21)",
                               panel[~panel["year"].isin([2020, 2021])]))
nb_res.append(run_count_model(sm.families.NegativeBinomial(alpha=1.0), "Excl 2016",
                               panel[panel["year"] != 2016]))


# =====================================================================
# 2. JOINT WALD F-TEST FOR PRE-TRENDS
# =====================================================================
print("\n" + "=" * 72)
print("2. JOINT WALD TEST FOR PRE-TRENDS")
print("=" * 72)
print("Estimates the full TWFE event-study and tests H0: all pre-treatment")
print("interaction coefficients = 0 jointly. Failure is informative.")

es_panel = panel.copy()
all_years = sorted(es_panel["year"].unique())
for yr in all_years:
    if yr != 2016:
        es_panel[f"tx{yr}"] = (es_panel["treat"] * (es_panel["year"] == yr)).astype(int)

interact_terms = [f"tx{yr}" for yr in all_years if yr != 2016]
formula_es = "gun_hom ~ C(ca_num) + C(year) + " + " + ".join(interact_terms)
model_es = smf.ols(formula_es, data=es_panel).fit(
    cov_type="cluster", cov_kwds={"groups": es_panel["ca_num"]})

pre_terms = [f"tx{yr}" for yr in all_years if yr < 2016]
post_terms = [f"tx{yr}" for yr in all_years if yr > 2016]

print(f"\n  Pre-treatment coefficients (years 2009-2015, base = 2016):")
for term in pre_terms:
    yr = int(term[2:])
    b = model_es.params[term]
    p = model_es.pvalues[term]
    sig = "***" if p < 0.01 else "**" if p < 0.05 else "*" if p < 0.1 else ""
    print(f"    {yr}: beta = {b:+.3f}  p = {p:.4f}{sig}")

# Wald test: H0 that all pre-treatment coefficients are zero
hypothesis_pre = " = ".join(pre_terms) + " = 0"
wald_pre = model_es.wald_test(hypothesis_pre, scalar=True)
print(f"\n  Joint Wald F-test, H0: all 7 pre-treatment coefficients = 0")
print(f"    F = {wald_pre.statistic:.3f}")
print(f"    p = {wald_pre.pvalue:.6f}")
pre_trend_level = "1%" if wald_pre.pvalue < 0.01 else "5%" if wald_pre.pvalue < 0.05 else None
if pre_trend_level is not None:
    print(f"    -> REJECT parallel trends at the {pre_trend_level} level. Design fails the")
    print(f"       standard pre-trend falsification test.")
else:
    print(f"    -> Fail to reject parallel trends.")

# Also test: are pre-period and post-period coefficients statistically different?
hypothesis_diff = " + ".join(pre_terms) + " = " + " + ".join(post_terms)
wald_diff = model_es.wald_test(hypothesis_diff, scalar=True)
print(f"\n  Wald test, H0: sum of pre-coefs = sum of post-coefs")
print(f"    F = {wald_diff.statistic:.3f}, p = {wald_diff.pvalue:.6f}")


# =====================================================================
# 3. PERMUTATION INFERENCE FOR SYNTHETIC CONTROL
# =====================================================================
print("\n" + "=" * 72)
print("3. PERMUTATION-BASED SYNTHETIC CONTROL INFERENCE")
print("=" * 72)
print("For each treated unit, also runs synthetic control on every never-treated")
print("unit (treating it as if it were 'placebo-treated'). The post-2017 gap of")
print("the actual treated unit is then ranked against the placebo distribution.")
print("Standard SC inference: p-value = fraction of placebos with |gap| >= |treated gap|.")


def synthetic_control(treated_name, donor_names, panel_df, install_year):
    pre_years = sorted(panel_df[panel_df["year"] < install_year]["year"].unique())
    treated_pre = (panel_df[(panel_df["ca_name"] == treated_name) & (panel_df["year"].isin(pre_years))]
                   .sort_values("year")["gun_hom"].values.astype(float))
    donor_pre = np.array([
        (panel_df[(panel_df["ca_name"] == d) & (panel_df["year"].isin(pre_years))]
         .sort_values("year")["gun_hom"].values.astype(float))
        for d in donor_names
    ])
    if len(donor_pre) == 0:
        return None, False
    def loss(w):
        synth = (w[:, None] * donor_pre).sum(axis=0)
        return float(np.sum((treated_pre - synth) ** 2))
    n_donors = len(donor_names)
    w0 = np.ones(n_donors) / n_donors
    constraints = ({"type": "eq", "fun": lambda w: np.sum(w) - 1.0},)
    bounds = [(0.0, 1.0) for _ in range(n_donors)]
    result = minimize(loss, w0, method="SLSQP", bounds=bounds, constraints=constraints,
                      options={"maxiter": 500, "ftol": 1e-9})
    return result.x, bool(result.success)


def gap_metric(unit_name, weights, donor_names, panel_df, install_year):
    """Mean post-treatment difference between unit and synthetic counterfactual."""
    treated_full = panel_df[panel_df["ca_name"] == unit_name].sort_values("year")["gun_hom"].values.astype(float)
    donor_full = np.array([
        panel_df[panel_df["ca_name"] == d].sort_values("year")["gun_hom"].values.astype(float)
        for d in donor_names
    ])
    synth = (weights[:, None] * donor_full).sum(axis=0)
    yrs = sorted(panel_df["year"].unique())
    post_idx = [i for i, y in enumerate(yrs) if y >= install_year]
    pre_idx = [i for i, y in enumerate(yrs) if y < install_year]
    pre_rmse = float(np.sqrt(np.mean((treated_full[pre_idx] - synth[pre_idx]) ** 2)))
    post_gap = float((treated_full[post_idx] - synth[post_idx]).mean())
    return post_gap, pre_rmse


# For each of the 6 treated units, compute the actual gap, then placebo gaps
sc_inference = {}
for hood in ORIG_6:
    install_yr = 2017 if hood in COHORT_2017 else 2018
    # Actual: synthetic control with the 26 never-treated as donors
    w_actual, conv_actual = synthetic_control(hood, NEVER_TREATED, panel, install_yr)
    actual_gap, actual_rmse = gap_metric(hood, w_actual, NEVER_TREATED, panel, install_yr)
    if not conv_actual:
        print(f"    [warn] SLSQP did not converge for {hood} (single-donor corner solution)")

    # Placebo: each never-treated unit, donors = the OTHER 25 never-treated
    placebo_gaps = []
    for placebo in NEVER_TREATED:
        donors = [d for d in NEVER_TREATED if d != placebo]
        try:
            w_p, _ = synthetic_control(placebo, donors, panel, install_yr)
            gap, rmse = gap_metric(placebo, w_p, donors, panel, install_yr)
            placebo_gaps.append((placebo, gap, rmse))
        except Exception:
            continue

    # Restrict placebos to those with reasonable pre-fit (RMSE <= 5x actual)
    valid_placebos = [(p, g, r) for p, g, r in placebo_gaps if r <= 5 * actual_rmse]
    placebo_gap_values = [g for _, g, _ in valid_placebos]
    n_placebo = len(placebo_gap_values)

    # Two-sided p-value: fraction of placebos with |gap| >= |actual|
    if n_placebo > 0:
        p_value = (np.abs(placebo_gap_values) >= np.abs(actual_gap)).sum() / n_placebo
    else:
        p_value = np.nan

    sc_inference[hood] = {
        "install_year": install_yr,
        "actual_gap": actual_gap,
        "actual_rmse": actual_rmse,
        "placebo_gaps": placebo_gap_values,
        "n_placebo": n_placebo,
        "p_value": p_value,
    }
    print(f"\n  {hood} (cohort {install_yr}):")
    print(f"    actual post-tx gap: {actual_gap:+.2f} per year (pre-tx RMSE: {actual_rmse:.2f})")
    if n_placebo > 0:
        pct_lo, pct_hi = np.percentile(placebo_gap_values, [2.5, 97.5])
        print(f"    placebo distribution (n={n_placebo}): "
              f"95% range [{pct_lo:+.2f}, {pct_hi:+.2f}], median {np.median(placebo_gap_values):+.2f}")
        print(f"    permutation p-value: {p_value:.3f}")


# Plot: actual gap vs placebo distribution for each treated unit
fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharex=False, sharey=False)
axes = axes.flatten()
for i, hood in enumerate(ORIG_6):
    ax = axes[i]
    inf = sc_inference[hood]
    if inf["n_placebo"] == 0:
        continue
    ax.hist(inf["placebo_gaps"], bins=15, color=C_GRAY, alpha=0.55,
            edgecolor="white")
    ax.axvline(inf["actual_gap"], color=C_RED, lw=3, zorder=3)
    ax.axvline(0, color="#333333", ls=":", lw=1)
    ax.set_title(f"{hood}\n({format_pvalue(inf['p_value'])} permutation)",
                 fontsize=12, fontweight="bold", pad=6)
    ax.set_xlabel("Post-tx gap (homicides/yr)", fontsize=10)
    ax.set_ylabel("Placebo count" if i % 3 == 0 else "", fontsize=10)
    ax.tick_params(axis="both", labelsize=9)
    ax.yaxis.grid(True)

# Single shared legend at the bottom
shared_legend = [
    mpatches.Patch(color=C_GRAY, alpha=0.55,
                   label=f"Placebo distribution (n=26 never-treated areas)"),
    Line2D([0], [0], color=C_RED, lw=3, label="Actual treated-unit gap"),
    Line2D([0], [0], color="#333", ls=":", lw=1, label="Zero gap"),
]
fig.legend(handles=shared_legend, loc="lower center", ncol=3, frameon=False,
           fontsize=11, bbox_to_anchor=(0.5, 0.01))

fig.suptitle("Synthetic Control: Treated Gap vs. Placebo Distribution",
             fontsize=15, fontweight="bold", y=0.99)
add_source_note(fig, SOURCE_DEFAULT, y=0.06, x=0.5)
for txt in fig.texts:
    if txt.get_text().startswith("Source:"):
        txt.set_ha("center")

fig.tight_layout(rect=[0, 0.07, 1, 0.96])
fig.savefig(os.path.join(PLOTS, "synthetic_control_inference.png"), dpi=200,
            facecolor="white")
plt.close(fig)
print(f"\n  SAVED: synthetic_control_inference.png")


# =====================================================================
# APPEND TO RESULTS_SUMMARY_V2.md
# =====================================================================
print("\n" + "=" * 72)
print("APPENDING TO RESULTS_SUMMARY_V2.md")
print("=" * 72)

lines = []
lines.append("\n---\n")
lines.append("# Robustness Layer (econometric_robustness.py)\n")
lines.append("Three additions on top of the prior re-analysis:\n")
lines.append("- Count-data models (Poisson, Negative Binomial)\n")
lines.append("- Joint Wald test on pre-trends\n")
lines.append("- Permutation inference for synthetic control\n\n")

lines.append("## 5. Count-Data Regression\n")
lines.append("Coefficients are on the log-mean scale; IRR = exp(beta).\n\n")
lines.append("### Poisson with two-way FE\n\n")
lines.append("| Specification | beta | SE | p | IRR | IRR 95% CI |\n")
lines.append("|---|---|---|---|---|---|\n")
for r in poisson_res:
    if r is None:
        continue
    sig = "***" if r["p"] < 0.01 else "**" if r["p"] < 0.05 else "*" if r["p"] < 0.1 else ""
    lines.append(f"| {r['label']} | {r['b']:+.4f} | {r['se']:.4f} | {r['p']:.4f}{sig} | "
                 f"{r['irr']:.4f} | [{r['irr_ci'][0]:.4f}, {r['irr_ci'][1]:.4f}] |\n")
lines.append("\n### Negative Binomial (alpha=1) with two-way FE\n\n")
lines.append("| Specification | beta | SE | p | IRR | IRR 95% CI |\n")
lines.append("|---|---|---|---|---|---|\n")
for r in nb_res:
    if r is None:
        continue
    sig = "***" if r["p"] < 0.01 else "**" if r["p"] < 0.05 else "*" if r["p"] < 0.1 else ""
    lines.append(f"| {r['label']} | {r['b']:+.4f} | {r['se']:.4f} | {r['p']:.4f}{sig} | "
                 f"{r['irr']:.4f} | [{r['irr_ci'][0]:.4f}, {r['irr_ci'][1]:.4f}] |\n")

lines.append("\n## 6. Joint Wald Test for Pre-Trends\n\n")
lines.append("Full TWFE event-study `gun_hom ~ C(ca_num) + C(year) + treat:1[year=k]`\n")
lines.append("for k in 2009..2023 (omitting 2016 as base).\n\n")
lines.append(f"H0: all 7 pre-treatment interaction coefficients (2009-2015) = 0\n\n")
lines.append(f"- F = {wald_pre.statistic:.3f}\n")
lines.append(f"- p = {wald_pre.pvalue:.6f}\n")
if pre_trend_level is not None:
    verdict = f"REJECT parallel trends at the {pre_trend_level} level"
else:
    verdict = "Fail to reject parallel trends"
lines.append(f"- **Verdict:** {verdict}.\n\n")

lines.append("## 7. Permutation Inference for Synthetic Control\n\n")
lines.append("For each treated unit, the post-treatment gap (actual - synthetic) is\n")
lines.append("compared to the distribution of analogous gaps obtained by running\n")
lines.append("synthetic control on each of the 26 never-treated areas. Permutation\n")
lines.append("p-value = fraction of placebos with absolute gap >= absolute treated gap.\n\n")
lines.append("| Neighborhood | Cohort | Actual gap | Pre-tx RMSE | Placebo n | Perm. p |\n")
lines.append("|---|---|---|---|---|---|\n")
for hood, inf in sc_inference.items():
    p_str = f"{inf['p_value']:.3f}" if not np.isnan(inf["p_value"]) else "n/a"
    lines.append(f"| {hood} | {inf['install_year']} | {inf['actual_gap']:+.2f} | "
                 f"{inf['actual_rmse']:.2f} | {inf['n_placebo']} | {p_str} |\n")

lines.append("\n> **These p-values are invalid by construction — do not read them as treatment effects.**\n")
lines.append("> The permutation p of 0.000 for every unit is an artifact of the synthetic-control\n")
lines.append("> infeasibility documented in §3: the treated units' pre-treatment RMSE (10-34) is far\n")
lines.append("> larger than the placebo never-treated units' (which fit themselves well), so the treated\n")
lines.append("> 'gaps' land far outside a tight placebo distribution **because of pre-treatment level\n")
lines.append("> mismatch, not because of any post-treatment effect.** The standard RMSE-ratio guard\n")
lines.append("> cannot rescue this when no donor pool can match the treated units. The honest reading is\n")
lines.append("> that synthetic control is infeasible for these data, not that the effect is significant.\n")
lines.append("> See NARRATIVE.md §5.5.\n")
lines.append("\nFigure: `plots/synthetic_control_inference.png`\n")

with open(os.path.join(P, "docs", "results_summary_reanalysis.md"), "a", encoding="utf-8") as f:
    f.write("".join(lines))

print("  Appended results to docs/results_summary_reanalysis.md")
print("\n" + "=" * 72)
print("DONE")
print("=" * 72)
