"""
placement_model.py — What predicted ShotSpotter placement? (selection / first stage)

The paper's central claim is that ShotSpotter was assigned *because of* pre-existing
violence, not at random, and that this selection is what makes a valid counterfactual
impossible to construct. That claim is asserted from group means. This script estimates
it directly: a cross-sectional model of ShotSpotter placement across all 77 Chicago
community areas on pre-period violence and 2008-2012 Census socioeconomic characteristics.

Three things come out of it:

  1. A balance table (treated vs. never-treated) with standardized differences and the
     univariate discriminating power (AUC) of each characteristic.
  2. A placement model (linear-probability + logit) and a horse-race that asks whether
     anything predicts placement *conditional on* pre-period violence.
  3. A propensity-score overlap plot. Because placement is so well predicted by violence,
     treated and control units barely share common support — the same fact, viewed from
     the selection side, that makes DiD/synthetic-control infeasible in the main paper.

Outputs:
  plots/placement_coefficients.png       — standardized LPM coefficients (forest plot)
  plots/placement_propensity_overlap.png — propensity distributions, treated vs control
  docs/results_summary_placement.md      — balance table + model tables + interpretation
"""
import os
import sys
import json
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import variance_inflation_factor
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_style import (apply_academic_style, COLORS, add_source_note,
                        add_n_label, sig_stars, SOURCE_DEFAULT)

P = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(P, "data", "raw")
PLOTS = os.path.join(P, "plots")
DOCS = os.path.join(P, "docs")

apply_academic_style()
C_BLUE = COLORS["control"]
C_RED = COLORS["treated"]
C_GRAY = COLORS["neutral"]

PRE_START, PRE_END = 2009, 2016   # pre-rollout window used to measure "violence at selection"


# =====================================================================
# DATA LOADING  (mirrors econometric_analysis.py conventions)
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

# --- Homicides (pre-period gun homicide counts) ---
hom_raw = pd.read_csv(
    os.path.join(DATA, "Violence_Reduction_-_Victims_of_Homicides_and_Non-Fatal_Shootings_20250513.csv"),
    low_memory=False,
)
hom_raw["DATE"] = pd.to_datetime(hom_raw["DATE"], format="%m/%d/%Y %I:%M:%S %p", errors="coerce")
hom_raw["CA_num"] = hom_raw["COMMUNITY_AREA"].map(ca_name2num)
hom_raw["year"] = hom_raw["DATE"].dt.year

gun_pre = hom_raw[
    (hom_raw["VICTIMIZATION_PRIMARY"] == "HOMICIDE") &
    (hom_raw["GUNSHOT_INJURY_I"] == "YES") &
    (hom_raw["year"] >= PRE_START) & (hom_raw["year"] <= PRE_END)
].dropna(subset=["CA_num"]).copy()
gun_pre["CA_num"] = gun_pre["CA_num"].astype(int)

pre_counts = gun_pre.groupby("CA_num").size().reindex(sorted(ca_num2name.keys()), fill_value=0)
pre_mean_annual = pre_counts / (PRE_END - PRE_START + 1)

# --- ShotSpotter treatment set (same rule as the main analysis) ---
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
    ss_raw.dropna(subset=["CA_name"])
    .groupby("CA_name")
    .agg(total_alerts=("DATE", "size"), earliest=("DATE", "min"), latest=("DATE", "max"))
)
ss_profile["months"] = ((ss_profile["latest"] - ss_profile["earliest"]).dt.days / 30.4).round().astype(int)
ss_profile = ss_profile[(ss_profile.total_alerts >= 100) & (ss_profile.months >= 12)]
ALL_SS = sorted(ss_profile.index.tolist())   # treated community-area names (title case)

# --- Census socioeconomic indicators (2008-2012) ---
cen = pd.read_csv(os.path.join(
    DATA, "Census_Data_-_Selected_socioeconomic_indicators_in_Chicago__2008___2012_20250408.csv"))
cen.columns = [c.strip() for c in cen.columns]
cen = cen[pd.to_numeric(cen["Community Area Number"], errors="coerce").notna()].copy()  # drop 'Chicago' total
cen["ca_num"] = cen["Community Area Number"].astype(int)

RENAME = {
    "PERCENT OF HOUSING CROWDED":                  "crowded",
    "PERCENT HOUSEHOLDS BELOW POVERTY":            "poverty",
    "PERCENT AGED 16+ UNEMPLOYED":                 "unemp",
    "PERCENT AGED 25+ WITHOUT HIGH SCHOOL DIPLOMA": "no_hs",
    "PERCENT AGED UNDER 18 OR OVER 64":            "dependency",
    "PER CAPITA INCOME":                           "income",
    "HARDSHIP INDEX":                              "hardship",
}
cen = cen.rename(columns=RENAME)

# --- Assemble the cross-section: one row per community area ---
df = pd.DataFrame({"ca_num": sorted(ca_num2name.keys())})
df["ca_name"] = df["ca_num"].map(ca_num2name)
df["treated"] = df["ca_name"].isin(ALL_SS).astype(int)
df["pre_gun_hom"] = df["ca_num"].map(pre_mean_annual.to_dict())
df = df.merge(cen[["ca_num"] + list(RENAME.values())], on="ca_num", how="left")

# income in $000s so its coefficient is on a comparable, readable scale before z-scoring
df["income"] = df["income"] / 1000.0

# --- Racial composition (ACS 2019-2023 5-year, by community area) ---
# Vintage caveat: the only community-area ACS race table the Chicago portal exposes is
# the most-recent 5-year (2019-2023), which post-dates the 2017/2018 rollout. Chicago's
# community-area racial composition is highly persistent decade to decade, so we use it as
# a proxy for the pre-rollout composition and flag the mismatch. Black share uses "Black or
# African American alone" (no Black-not-Hispanic field is published; the Black-Hispanic
# overlap is small in Chicago).
def _norm(s):
    return "".join(ch for ch in str(s).upper() if ch.isalnum())
ca_norm2num = {_norm(name): n for n, name in ca_num2name.items()}
race = pd.read_csv(os.path.join(DATA, "ACS_5yr_race_by_community_area_2023.csv"))
race["ca_num"] = race["community_area"].map(lambda s: ca_norm2num.get(_norm(s)))
_unm = race[race["ca_num"].isna()]["community_area"].tolist()
if _unm:
    print(f"  [warn] race rows unmatched to a community area: {_unm}")
race["pct_black"] = race["black_or_african_american"] / race["total_population"] * 100
race["pct_hispanic"] = race["hispanic_or_latino"] / race["total_population"] * 100
race["pct_white_nh"] = race["white_not_hispanic_or_latino"] / race["total_population"] * 100
race = race.dropna(subset=["ca_num"]).copy()
race["ca_num"] = race["ca_num"].astype(int)
df = df.merge(race[["ca_num", "pct_black", "pct_hispanic", "pct_white_nh"]], on="ca_num", how="left")

# --- Life expectancy (Chicago Public Health Statistics, 2010 by community area) ---
# Treated as a descriptive health/disadvantage summary in the balance table, not a siting
# predictor: the city did not site on life expectancy, but that ShotSpotter areas carry a
# large life-expectancy penalty is a vivid statement of who bears the coverage.
le = pd.read_csv(os.path.join(DATA, "Public_Health_Statistics_Life_Expectancy_By_Community_Area.csv"))
le["ca_num"] = le["community_area"].map(lambda s: ca_norm2num.get(_norm(s)))
le["life_exp"] = pd.to_numeric(le["_2010_life_expectancy"], errors="coerce")
le = le.dropna(subset=["ca_num"]).copy()
le["ca_num"] = le["ca_num"].astype(int)
df = df.merge(le[["ca_num", "life_exp"]], on="ca_num", how="left")

PREDICTORS = ["pre_gun_hom", "poverty", "unemp", "no_hs", "income", "hardship",
              "crowded", "dependency"]
RACE = ["pct_black", "pct_hispanic", "pct_white_nh"]
HEALTH = ["life_exp"]   # descriptive balance only; kept out of the placement regressions
PRETTY = {
    "pre_gun_hom": "Pre-period gun homicides (mean/yr, 2009-16)",
    "poverty":     "Households below poverty (%)",
    "unemp":       "Unemployed, 16+ (%)",
    "no_hs":       "No high-school diploma, 25+ (%)",
    "income":      "Per-capita income ($000s)",
    "hardship":    "Hardship index",
    "crowded":     "Crowded housing (%)",
    "dependency":  "Dependency: under 18 or over 64 (%)",
    "pct_black":    "Black share of population (%)",
    "pct_hispanic": "Hispanic/Latino share (%)",
    "pct_white_nh": "White, non-Hispanic share (%)",
    "life_exp":     "Life expectancy, 2010 (years)",
}

df = df.dropna(subset=PREDICTORS + RACE + HEALTH).reset_index(drop=True)
n_missing = 77 - len(df)
if n_missing:
    print(f"  [note] dropped {n_missing} community area(s) missing Census, ACS, or health indicators")

# z-scored copies so coefficients are comparable "per one SD" effects
for v in PREDICTORS + RACE + HEALTH:
    df[f"z_{v}"] = (df[v] - df[v].mean()) / df[v].std(ddof=0)
ZPRED = [f"z_{v}" for v in PREDICTORS]

n_t = int(df["treated"].sum())
n_c = int((1 - df["treated"]).sum())
print(f"  Community areas analyzed: {len(df)}  (treated={n_t}, never-treated={n_c})")
print(f"  Pre-period gun homicides/yr: treated mean={df.loc[df.treated==1,'pre_gun_hom'].mean():.2f}, "
      f"control mean={df.loc[df.treated==0,'pre_gun_hom'].mean():.2f}")


# =====================================================================
# 1. BALANCE TABLE  (treated vs. never-treated)
# =====================================================================
print("\n" + "=" * 72)
print("1. BALANCE TABLE: what distinguishes treated from never-treated areas")
print("=" * 72)

balance = []
for v in PREDICTORS + HEALTH:
    t = df.loc[df.treated == 1, v]
    c = df.loc[df.treated == 0, v]
    pooled_sd = np.sqrt(((t.var(ddof=1) * (len(t) - 1)) + (c.var(ddof=1) * (len(c) - 1))) /
                        (len(t) + len(c) - 2))
    smd = (t.mean() - c.mean()) / pooled_sd if pooled_sd > 0 else np.nan
    # univariate AUC: how well does this single variable rank placement?
    auc = roc_auc_score(df["treated"], df[v])
    auc = max(auc, 1 - auc)   # direction-free discriminating power
    balance.append({"var": v, "treated_mean": t.mean(), "control_mean": c.mean(),
                    "smd": smd, "auc": auc})
balance = pd.DataFrame(balance).sort_values("auc", ascending=False)

print(f"  {'Characteristic':<44} {'Treated':>9} {'Control':>9} {'Std.diff':>9} {'AUC':>6}")
for _, r in balance.iterrows():
    print(f"  {PRETTY[r['var']]:<44} {r['treated_mean']:>9.2f} {r['control_mean']:>9.2f} "
          f"{r['smd']:>+9.2f} {r['auc']:>6.3f}")

le_t = df.loc[df.treated == 1, "life_exp"].mean()
le_c = df.loc[df.treated == 0, "life_exp"].mean()
print(f"\n  LIFE-EXPECTANCY GAP: ShotSpotter areas {le_t:.1f} yrs vs {le_c:.1f} in never-treated "
      f"areas = {le_t - le_c:+.1f} years (2010).")


# =====================================================================
# 2. PLACEMENT MODEL (LPM + logit) and the violence horse-race
# =====================================================================
print("\n" + "=" * 72)
print("2. PLACEMENT MODEL: treated ~ characteristics (standardized)")
print("=" * 72)

def fit_lpm(cols, label):
    f = "treated ~ " + " + ".join(cols)
    m = smf.ols(f, data=df).fit(cov_type="HC1")
    print(f"\n  [{label}]  LPM, R^2={m.rsquared:.3f}, N={int(m.nobs)}")
    for c in cols:
        print(f"    {c:<16} b={m.params[c]:+.4f}  SE={m.bse[c]:.4f}  p={m.pvalues[c]:.4f} {sig_stars(m.pvalues[c])}")
    return m

def logit_fit(cols):
    """Near-unpenalized L2 logit (sklearn) — robust to quasi-separation. Returns (probs, auc, ps_R2)."""
    X = df[cols].values
    y = df["treated"].values
    clf = LogisticRegression(penalty="l2", C=1e6, solver="lbfgs", max_iter=10000)
    clf.fit(X, y)
    p = clf.predict_proba(X)[:, 1]
    auc = roc_auc_score(y, p)
    # McFadden-style pseudo-R^2 via log-loss against the null (base rate) model
    eps = 1e-12
    ll = np.sum(y * np.log(p + eps) + (1 - y) * np.log(1 - p + eps))
    pbar = y.mean()
    ll0 = np.sum(y * np.log(pbar + eps) + (1 - y) * np.log(1 - pbar + eps))
    ps_r2 = 1 - ll / ll0
    return p, auc, ps_r2, clf

# Horse-race: how much does violence alone explain, and does disadvantage add to it?
models = {
    "Violence only":            ["z_pre_gun_hom"],
    "Disadvantage only (hardship+income)": ["z_hardship", "z_income"],
    "Violence + hardship + income":        ["z_pre_gun_hom", "z_hardship", "z_income"],
    "Full (violence + all SES)":           ZPRED,
}
lpm_results = {}
horse = []
for label, cols in models.items():
    m = fit_lpm(cols, label)
    _, auc, ps_r2, _ = logit_fit(cols)
    lpm_results[label] = m
    horse.append({"model": label, "cols": cols, "r2": m.rsquared, "auc": auc, "ps_r2": ps_r2})
    print(f"    -> logit AUC={auc:.3f}, pseudo-R^2={ps_r2:.3f}")

print("\n  Horse-race summary (does anything beat violence at predicting placement?):")
print(f"  {'Model':<38} {'LPM R^2':>8} {'AUC':>7} {'pseudoR2':>9}")
for h in horse:
    print(f"  {h['model']:<38} {h['r2']:>8.3f} {h['auc']:>7.3f} {h['ps_r2']:>9.3f}")

# VIF diagnostic for the full model (hardship is a composite of several components)
print("\n  Variance inflation factors (full model — collinearity check):")
Xv = df[ZPRED].copy()
Xv.insert(0, "const", 1.0)
for i, c in enumerate(ZPRED):
    vif = variance_inflation_factor(Xv.values, i + 1)
    print(f"    {c:<16} VIF={vif:5.1f}")

# Headline marginal effect from the parsimonious model
head_m = lpm_results["Violence + hardship + income"]
b_v = head_m.params["z_pre_gun_hom"]
p_v = head_m.pvalues["z_pre_gun_hom"]
print(f"\n  HEADLINE: a 1-SD increase in pre-period gun homicides is associated with a "
      f"{b_v*100:+.1f} pp change in P(ShotSpotter), p={p_v:.4f}, "
      f"holding hardship and income fixed.")


# =====================================================================
# 2b. EQUITY: does racial composition predict placement conditional on
#     violence and hardship? (the concentrated-surveillance question)
# =====================================================================
print("\n" + "=" * 72)
print("2b. EQUITY: does race predict placement on top of violence + hardship?")
print("=" * 72)

# Raw racial gaps between treated and never-treated areas
race_balance = []
for v in RACE:
    t = df.loc[df.treated == 1, v]
    c = df.loc[df.treated == 0, v]
    pooled_sd = np.sqrt(((t.var(ddof=1) * (len(t) - 1)) + (c.var(ddof=1) * (len(c) - 1))) /
                        (len(t) + len(c) - 2))
    smd = (t.mean() - c.mean()) / pooled_sd if pooled_sd > 0 else np.nan
    auc = roc_auc_score(df["treated"], df[v]); auc = max(auc, 1 - auc)
    race_balance.append({"var": v, "treated_mean": t.mean(), "control_mean": c.mean(),
                         "smd": smd, "auc": auc})
race_balance = pd.DataFrame(race_balance)
print(f"  {'Characteristic':<32} {'Treated':>9} {'Control':>9} {'Std.diff':>9} {'AUC':>6}")
for _, r in race_balance.iterrows():
    print(f"  {PRETTY[r['var']]:<32} {r['treated_mean']:>9.1f} {r['control_mean']:>9.1f} "
          f"{r['smd']:>+9.2f} {r['auc']:>6.3f}")

# Does race add predictive power on top of violence + hardship + income?
base_cols = ["z_pre_gun_hom", "z_hardship", "z_income"]
race_cols = base_cols + ["z_pct_black", "z_pct_hispanic"]
_, auc_base, psr2_base, _ = logit_fit(base_cols)
_, auc_race, psr2_race, _ = logit_fit(race_cols)
print(f"\n  Predictive gain from adding race:")
print(f"    violence + hardship + income        : logit AUC={auc_base:.3f}, pseudo-R^2={psr2_base:.3f}")
print(f"    + Black share + Hispanic share       : logit AUC={auc_race:.3f}, pseudo-R^2={psr2_race:.3f}")

# Conditional coefficients: race given violence + hardship + income (LPM, HC1)
eq_m = fit_lpm(race_cols, "Placement ~ violence + hardship + income + race")
b_black = eq_m.params["z_pct_black"]; p_black = eq_m.pvalues["z_pct_black"]
b_hisp = eq_m.params["z_pct_hispanic"]; p_hisp = eq_m.pvalues["z_pct_hispanic"]
print(f"\n  CONDITIONAL ON violence + hardship + income:")
print(f"    Black share:    {b_black*100:+.1f} pp per SD, p={p_black:.4f} {sig_stars(p_black)}")
print(f"    Hispanic share: {b_hisp*100:+.1f} pp per SD, p={p_hisp:.4f} {sig_stars(p_hisp)}")
t_black = df.loc[df.treated == 1, "pct_black"].mean()
c_black = df.loc[df.treated == 0, "pct_black"].mean()
t_hisp = df.loc[df.treated == 1, "pct_hispanic"].mean()
c_hisp = df.loc[df.treated == 0, "pct_hispanic"].mean()
print(f"\n  Raw composition: ShotSpotter areas are {t_black:.0f}% Black / {t_hisp:.0f}% Hispanic "
      f"on average, vs {c_black:.0f}% / {c_hisp:.0f}% in never-treated areas.")


# =====================================================================
# 3. PROPENSITY-SCORE OVERLAP  (the selection-side view of 'no counterfactual')
# =====================================================================
print("\n" + "=" * 72)
print("3. PROPENSITY OVERLAP: do treated and control share common support?")
print("=" * 72)

ps, ps_auc, ps_r2, _ = logit_fit(["z_pre_gun_hom", "z_hardship", "z_income"])
df["propensity"] = ps
t_ps = df.loc[df.treated == 1, "propensity"]
c_ps = df.loc[df.treated == 0, "propensity"]
# The six highest-violence study neighborhoods sit at the top of the propensity scale;
# the honest "no counterfactual" statement is that no control lives up there with them.
STUDY6 = ["Austin", "Humboldt Park", "North Lawndale", "Englewood", "West Englewood", "South Shore"]
study6_ps = df[df["ca_name"].isin(STUDY6)]["propensity"]
c_ge_08 = int((c_ps >= 0.8).sum())
c_ge_09 = int((c_ps >= 0.9).sum())
study6_lo = study6_ps.min()
c_above_study6 = int((c_ps >= study6_lo).sum())
print(f"  Propensity model AUC={ps_auc:.3f}")
print(f"  Treated propensity: median {t_ps.median():.2f}  |  Control propensity: median {c_ps.median():.2f}")
print(f"  Never-treated areas with propensity >= 0.8: {c_ge_08} of {n_c};  >= 0.9: {c_ge_09} of {n_c}")
print(f"  Six study neighborhoods' propensity: min={study6_ps.min():.2f}, all >= {study6_lo:.2f}")
print(f"  Never-treated areas as high-propensity as the LEAST-extreme study neighborhood: "
      f"{c_above_study6} of {n_c}")


# =====================================================================
# PLOT A: standardized differences, treated vs. never-treated
# (A multivariate forest plot would be dominated by collinearity — the
#  hardship index has VIF ~55 because it is a composite of the other
#  Census variables. Standardized mean differences show each
#  characteristic's raw gap without that artifact; the conditional
#  "violence adds little once hardship is held fixed" point is made by
#  the horse-race table, not by disentangling collinear coefficients.)
# =====================================================================
border = balance.sort_values("smd", key=lambda s: s.abs(), ascending=True)
labels = [PRETTY[v] for v in border["var"]]
smds = border["smd"].values
aucs = border["auc"].values

fig, ax = plt.subplots(figsize=(10, 6))
ypos = np.arange(len(border))
for y, s, v in zip(ypos, smds, border["var"]):
    color = C_RED if s > 0 else C_BLUE
    ax.plot([0, s], [y, y], color=color, lw=2.4, alpha=0.85, zorder=2)
    ax.plot(s, y, "o", color=color, markersize=10, markeredgecolor="white",
            markeredgewidth=0.8, zorder=3)
for y, s, a in zip(ypos, smds, aucs):
    ax.text(s + (0.06 if s > 0 else -0.06), y, f"AUC {a:.2f}",
            va="center", ha="left" if s > 0 else "right", fontsize=8.5,
            color=COLORS["text_muted"])
ax.axvline(0, color="#333", lw=1.2, zorder=1)
ax.set_yticks(ypos)
ax.set_yticklabels(labels, fontsize=10)
ax.set_xlim(-2.7, 3.1)
ax.set_xlabel("Standardized difference, ShotSpotter minus never-treated (pooled SDs)\n"
              "Positive = higher in ShotSpotter areas.  AUC = single-variable discriminating power.")
ax.set_title("What distinguished ShotSpotter areas from the rest of Chicago",
             fontweight="bold", pad=10)
ax.xaxis.grid(True)
add_n_label(ax, f"N = {len(df)} areas\n{n_t} ShotSpotter / {n_c} never-treated", loc="lower left")
legend = [
    Line2D([0], [0], color=C_RED, marker="o", lw=2.4, markersize=9, label="Higher in ShotSpotter areas"),
    Line2D([0], [0], color=C_BLUE, marker="o", lw=2.4, markersize=9, label="Lower in ShotSpotter areas"),
]
ax.legend(handles=legend, loc="upper left", frameon=False, fontsize=9)
add_source_note(fig, SOURCE_DEFAULT + "  Census 2008-2012 indicators; Public Health life expectancy (2010).")
fig.tight_layout(rect=[0, 0.03, 1, 1])
fig.savefig(os.path.join(PLOTS, "placement_coefficients.png"), dpi=200, facecolor="white")
plt.close(fig)
print("\n  SAVED: placement_coefficients.png")


# =====================================================================
# PLOT B: propensity-score overlap
# =====================================================================
fig, ax = plt.subplots(figsize=(10, 5.5))
bins = np.linspace(0, 1, 21)
ax.hist(c_ps, bins=bins, color=C_BLUE, alpha=0.65, label=f"Never-treated (n={n_c})",
        edgecolor="white", linewidth=0.5)
ax.hist(t_ps, bins=bins, color=C_RED, alpha=0.65, label=f"ShotSpotter (n={n_t})",
        edgecolor="white", linewidth=0.5)
# rug marks
ax.plot(c_ps, np.full_like(c_ps, -0.4), "|", color=C_BLUE, markersize=10, mew=1.5)
ax.plot(t_ps, np.full_like(t_ps, -0.9), "|", color=C_RED, markersize=10, mew=1.5)
ax.axhline(0, color="#999", lw=0.6)
ax.set_xlabel("Estimated propensity to receive ShotSpotter\n"
              "(logit on pre-period violence, hardship, and per-capita income)")
ax.set_ylabel("Number of community areas")
ax.axvspan(study6_lo, 1.02, color=COLORS["accent"], alpha=0.10, zorder=0)
ax.set_title("Placement is highly predictable — and no control matches the study areas\n"
             "The selection-side view of why no valid counterfactual exists",
             fontweight="bold", pad=10)
ax.legend(frameon=False, fontsize=10, loc="upper center")
add_n_label(ax, f"Propensity AUC = {ps_auc:.2f}\n"
                f"Control median {c_ps.median():.2f} vs treated {t_ps.median():.2f}\n"
                f"Six study areas all ≥ {study6_lo:.2f} (shaded);\n"
                f"{c_above_study6} of {n_c} controls reach that band",
            loc="upper right")
add_source_note(fig, SOURCE_DEFAULT + "  Census 2008-2012 indicators; Public Health life expectancy (2010).")
fig.tight_layout(rect=[0, 0.03, 1, 1])
fig.savefig(os.path.join(PLOTS, "placement_propensity_overlap.png"), dpi=200, facecolor="white")
plt.close(fig)
print("  SAVED: placement_propensity_overlap.png")


# =====================================================================
# WRITE docs/results_summary_placement.md
# =====================================================================
print("\n" + "=" * 72)
print("WRITING docs/results_summary_placement.md")
print("=" * 72)

L = []
L.append("# ShotSpotter Placement Model: What Predicted Treatment?")
L.append("")
L.append("Results of `scripts/placement_model.py`. This is the *selection* / first-stage")
L.append("analysis: a cross-section of all 77 Chicago community areas in which ShotSpotter")
L.append("placement is modeled on pre-period (2009-2016) gun-homicide levels and the")
L.append("2008-2012 Census socioeconomic indicators. It converts the paper's asserted")
L.append("selection-on-violence story into a directly estimated fact, and shows the")
L.append("selection-side reason no valid counterfactual exists.")
L.append("")
L.append(f"Sample: {len(df)} community areas ({n_t} ShotSpotter, {n_c} never-treated). "
         "Treatment defined as in the main analysis (≥100 alerts and ≥12 months of coverage).")
L.append("")
L.append("## 1. Balance table (treated vs. never-treated)")
L.append("")
L.append("`Std. diff` is the standardized mean difference (treated − control, in pooled SDs). "
         "`AUC` is the direction-free area under the ROC curve for that single characteristic "
         "as a classifier of placement (0.5 = no discriminating power, 1.0 = perfect).")
L.append("")
L.append("| Characteristic | Treated mean | Control mean | Std. diff | Univariate AUC |")
L.append("|---|---|---|---|---|")
for _, r in balance.iterrows():
    L.append(f"| {PRETTY[r['var']]} | {r['treated_mean']:.2f} | {r['control_mean']:.2f} | "
             f"{r['smd']:+.2f} | {r['auc']:.3f} |")
L.append("")
L.append(f"The health penalty is the most human summary of the gap: ShotSpotter areas averaged "
         f"**{le_t:.1f} years** of life expectancy in 2010 versus **{le_c:.1f}** in never-treated "
         f"areas---a **{le_t - le_c:+.1f}-year** difference. Life expectancy is a descriptive "
         f"disadvantage summary here, not a siting predictor (the city did not site on it); it is "
         f"kept out of the placement regressions below and reported only to characterize who bears "
         f"the coverage. Source: Chicago Public Health Statistics, Life Expectancy by Community "
         f"Area (`qjr3-bm53`).")
L.append("")
L.append("## 2. Placement model and the violence horse-race")
L.append("")
L.append("Linear-probability models (HC1 robust SEs) with standardized predictors, so each")
L.append("coefficient is the change in the probability of placement per one-SD increase.")
L.append("The logit AUC and pseudo-R² come from a near-unpenalized logistic fit on the same")
L.append("predictors. The question: does anything predict placement *on top of* violence?")
L.append("")
L.append("| Model | LPM R² | Logit AUC | Pseudo-R² |")
L.append("|---|---|---|---|")
for h in horse:
    L.append(f"| {h['model']} | {h['r2']:.3f} | {h['auc']:.3f} | {h['ps_r2']:.3f} |")
L.append("")
L.append("Standardized coefficients, parsimonious model (`Violence + hardship + income`):")
L.append("")
L.append("| Predictor | b (Δ prob per SD) | SE | p |")
L.append("|---|---|---|---|")
for c in ["z_pre_gun_hom", "z_hardship", "z_income"]:
    L.append(f"| {c} | {head_m.params[c]:+.4f} | {head_m.bse[c]:.4f} | {head_m.pvalues[c]:.4f} |")
L.append("")
L.append(f"**Headline.** The dominant predictor of placement is structural disadvantage, not "
         f"the gun-homicide count itself. The hardship index alone classifies placement with "
         f"AUC {[h['auc'] for h in horse if h['model']=='Disadvantage only (hardship+income)'][0]:.2f} "
         f"(LPM R² {[h['r2'] for h in horse if h['model']=='Disadvantage only (hardship+income)'][0]:.2f}), "
         f"versus AUC {[h['auc'] for h in horse if h['model']=='Violence only'][0]:.2f} "
         f"(R² {[h['r2'] for h in horse if h['model']=='Violence only'][0]:.2f}) for pre-period violence "
         f"alone. Adding violence on top of hardship barely moves fit "
         f"(R² {[h['r2'] for h in horse if h['model']=='Disadvantage only (hardship+income)'][0]:.2f} "
         f"→ {[h['r2'] for h in horse if h['model']=='Violence + hardship + income'][0]:.2f}), and "
         f"the violence coefficient is not significant once hardship is held fixed "
         f"({b_v*100:+.1f} pp, p = {p_v:.2f}). Read carefully: violence and hardship are entangled "
         f"(the same neighborhoods score high on both), so this does not say violence was irrelevant "
         f"— it says the city's siting tracked the broad map of disinvestment at least as tightly "
         f"as it tracked the specific gun-homicide burden. This refines, rather than overturns, the "
         f"paper's selection-on-violence premise: treatment was assigned on characteristics "
         f"(hardship, violence) that are themselves strong predictors of the outcome, which is what "
         f"makes a clean counterfactual impossible.")
L.append("")
L.append("**Caveat on the violence measure.** Pre-period violence enters as the gun-homicide "
         "*count* (mean/yr, 2009-2016), consistent with the rest of the paper and with the plausible "
         "premise that the city responded to the absolute burden of gunfire. A per-capita *rate* "
         "would be the natural robustness check but requires community-area population, which the "
         "2008-2012 indicators file does not carry. Adding it (one ACS pull) is a listed extension.")
L.append("")
L.append("## 2b. Equity: did placement track racial composition?")
L.append("")
_wnh = race_balance[race_balance["var"] == "pct_white_nh"].iloc[0]
L.append(f"ShotSpotter areas are on average {t_black:.0f}% Black and {t_hisp:.0f}% Hispanic, against "
         f"{c_black:.0f}% and {c_hisp:.0f}% in never-treated areas; the White non-Hispanic share is the "
         f"single sharpest separator in the entire analysis (treated {_wnh['treated_mean']:.0f}% vs "
         f"control {_wnh['control_mean']:.0f}%, AUC {_wnh['auc']:.2f}).")
L.append("")
L.append("| Composition | Treated mean | Control mean | Std. diff | Univariate AUC |")
L.append("|---|---|---|---|---|")
for _, r in race_balance.iterrows():
    L.append(f"| {PRETTY[r['var']]} | {r['treated_mean']:.1f} | {r['control_mean']:.1f} | "
             f"{r['smd']:+.2f} | {r['auc']:.3f} |")
L.append("")
L.append(f"Adding Black and Hispanic share to the placement model lifts logit AUC from "
         f"{auc_base:.3f} to {auc_race:.3f}, and the association survives conditioning: in an LPM on "
         f"pre-period violence, hardship, income, and race, a one-SD increase in the Black share is "
         f"associated with a {b_black*100:+.1f} pp change in placement probability (p = {p_black:.3f}) "
         f"and the Hispanic share {b_hisp*100:+.1f} pp (p = {p_hisp:.3f}), while violence and hardship "
         f"themselves turn insignificant. Read carefully: race, hardship, and violence are deeply "
         f"entangled in a segregated city, and with 77 areas the model cannot cleanly separate them, so "
         f"this is *not* evidence that race drove siting independently of disadvantage. What it shows, "
         f"descriptively, is that ShotSpotter placement tracked Chicago's racial geography at least as "
         f"tightly as it tracked violence or hardship---the concentrated-surveillance pattern advocacy "
         f"groups described.")
L.append("")
L.append("*Vintage caveat:* the community-area race figures are the 2019--2023 ACS 5-year (the only "
         "community-area race table the Chicago open-data portal exposes), which post-dates the "
         "2017/2018 rollout. Chicago's community-area racial composition is highly persistent across "
         "decades, so we use it as a proxy for the pre-rollout composition; a 2008--2012-vintage race "
         "table (via the Census tract API) would remove the mismatch.")
L.append("")
L.append("## 3. Propensity-score overlap")
L.append("")
L.append(f"Propensity model (violence + hardship + income) AUC = {ps_auc:.3f}. Median propensity "
         f"is {t_ps.median():.2f} for ShotSpotter areas versus {c_ps.median():.2f} for never-treated "
         f"areas. Separation is strong but not complete — the two groups mix in the middle of the "
         f"scale, where a handful of lower-hardship areas received ShotSpotter and a handful of "
         f"higher-hardship areas did not.")
L.append("")
L.append(f"The decisive fact is at the top of the scale. The six highest-violence study "
         f"neighborhoods all carry a placement propensity of at least {study6_lo:.2f}, and "
         f"{c_above_study6} of {n_c} never-treated areas reach that band "
         f"({c_ge_09} reach 0.9). These are exactly the units the synthetic control in the main "
         f"paper cannot reproduce. The placement model therefore states the identification problem "
         f"from the selection side: the city installed ShotSpotter in essentially every "
         f"highest-violence, highest-hardship area, leaving no comparably situated untreated unit "
         f"from which to build a counterfactual for them.")
L.append("")
L.append("## Figures")
L.append("")
L.append("- `plots/placement_coefficients.png` — standardized differences (treated vs. never-treated), by characteristic")
L.append("- `plots/placement_propensity_overlap.png` — propensity distributions, treated vs. control")
L.append("")
L.append("## Data note")
L.append("")
L.append("Race data: `data/raw/ACS_5yr_race_by_community_area_2023.csv` (Chicago Data Portal dataset "
         "t68z-cikk, ACS 2019--2023 5-year by community area). All other inputs as in the main pipeline.")

with open(os.path.join(DOCS, "results_summary_placement.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(L))
print("  WROTE: docs/results_summary_placement.md")

print("\n" + "=" * 72)
print("DONE")
print("=" * 72)
