"""
did_analysis.py  --  Difference-in-Differences analysis, statistical testing,
                     and publication-quality plots for ACE 592 Final Project.
"""
import os, sys, warnings
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from scipy import stats
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_style import (apply_academic_style, COLORS, format_pvalue,
                        sig_stars, add_source_note, add_n_label, SOURCE_DEFAULT)

P  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA  = os.path.join(P, "data", "raw")
CLEAN = os.path.join(P, "data", "processed")
PLOTS = os.path.join(P, "plots")

apply_academic_style()

C_BLUE  = COLORS["control"]
C_RED   = COLORS["treated"]
C_GREEN = COLORS["green"]
C_GRAY  = COLORS["neutral"]

ORIG_6 = ["Austin", "Humboldt Park", "North Lawndale",
           "Englewood", "West Englewood", "South Shore"]

def savefig(fig, name):
    p = os.path.join(PLOTS, name)
    fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  SAVED: {name}  ({os.path.getsize(p)/1024:.0f} KB)")

# ==================================================================
# SECTION 0 : DATA LOADING
# ==================================================================
print("=" * 72)
print("LOADING DATA")
print("=" * 72)

# --- GeoJSON  (community area number <-> name mapping) ---
gdf = gpd.read_file(os.path.join(DATA, "CommAreas_20250408.geojson"))
ca_num2name = {}
ca_name2num = {}
for _, r in gdf.iterrows():
    num  = int(float(r["area_num_1"]))
    name = r["community"].title()
    ca_num2name[num]  = name
    ca_name2num[name] = num
    ca_name2num[name.upper()] = num   # for matching ALL-CAPS
print(f"  GeoJSON: {len(gdf)} community areas")

# --- Census ---
census = pd.read_csv(os.path.join(DATA,
    "Census_Data_-_Selected_socioeconomic_indicators_in_Chicago__2008___2012_20250408.csv"))
census["CA_name"] = census["COMMUNITY AREA NAME"].str.title()
census_map = census.set_index("CA_name")
print(f"  Census:  {len(census)} rows")

# --- Raw homicide / shooting data ---
hom_raw = pd.read_csv(os.path.join(DATA,
    "Violence_Reduction_-_Victims_of_Homicides_and_Non-Fatal_Shootings_20250513.csv"),
    low_memory=False)
hom_raw["DATE"] = pd.to_datetime(hom_raw["DATE"],
                                  format="%m/%d/%Y %I:%M:%S %p", errors="coerce")
# COMMUNITY_AREA is ALL-CAPS name in this dataset
hom_raw["CA_name"] = hom_raw["COMMUNITY_AREA"].str.title()
hom_raw["CA_num"]  = hom_raw["COMMUNITY_AREA"].map(ca_name2num)
hom_raw["year"]    = hom_raw["DATE"].dt.year

# Gun homicides 2009-2023
gun = hom_raw[
    (hom_raw["VICTIMIZATION_PRIMARY"] == "HOMICIDE") &
    (hom_raw["GUNSHOT_INJURY_I"] == "YES") &
    (hom_raw["DATE"] >= "2009-01-01") &
    (hom_raw["DATE"] <= "2023-12-31")
].copy()
print(f"  Gun homicides 2009-2023: {len(gun):,}")

# Non-fatal gunshot injuries 2009-2023
nfs = hom_raw[
    (hom_raw["VICTIMIZATION_PRIMARY"] != "HOMICIDE") &
    (hom_raw["GUNSHOT_INJURY_I"] == "YES") &
    (hom_raw["DATE"] >= "2009-01-01") &
    (hom_raw["DATE"] <= "2023-12-31")
].copy()
print(f"  Non-fatal shootings 2009-2023: {len(nfs):,}")

# --- Raw ShotSpotter data ---
ss_raw = pd.read_csv(os.path.join(DATA,
    "Violence_Reduction_-_Shotspotter_Alerts_-_Historical_20250408.csv"))
ss_raw["DATE"] = pd.to_datetime(ss_raw["DATE"],
                                 format="%m/%d/%Y %I:%M:%S %p", errors="coerce")
ss_raw = ss_raw[ss_raw["DATE"] <= "2023-12-31"].copy()
# Determine if COMMUNITY_AREA is name or number
sample_ca = ss_raw["COMMUNITY_AREA"].dropna().iloc[0]
if str(sample_ca).replace(".", "").isdigit():
    # numeric -> map to name
    ss_raw["CA_num"]  = pd.to_numeric(ss_raw["COMMUNITY_AREA"], errors="coerce").astype("Int64")
    ss_raw["CA_name"] = ss_raw["CA_num"].map(ca_num2name)
else:
    ss_raw["CA_name"] = ss_raw["COMMUNITY_AREA"].str.title()
    ss_raw["CA_num"]  = ss_raw["CA_name"].map(ca_name2num)
print(f"  ShotSpotter alerts (through 2023): {len(ss_raw):,}")

# ==================================================================
# SECTION 1 : TASK 1  --  IDENTIFY ALL SHOTSPOTTER NEIGHBORHOODS
# ==================================================================
print("\n" + "=" * 72)
print("TASK 1: IDENTIFY ALL SHOTSPOTTER NEIGHBORHOODS")
print("=" * 72)

ss_profile = (ss_raw.dropna(subset=["CA_name"])
              .groupby("CA_name")
              .agg(total_alerts=("DATE", "size"),
                   earliest=("DATE", "min"),
                   latest=("DATE", "max"))
              .sort_values("total_alerts", ascending=False))
ss_profile["months"] = ((ss_profile["latest"] - ss_profile["earliest"])
                         .dt.days / 30.4).round().astype(int)

print(f"\n{'Community Area':<30s} {'Alerts':>8s} {'Months':>7s}  {'Earliest':>12s}")
print("-" * 70)
for name, row in ss_profile.iterrows():
    print(f"  {name:<28s} {row.total_alerts:>8,d} {row.months:>7d}  "
          f"{row.earliest.strftime('%Y-%m-%d'):>12s}")

# Threshold: >= 100 alerts AND >= 12 months
ALL_SS = sorted(ss_profile[
    (ss_profile.total_alerts >= 100) & (ss_profile.months >= 12)
].index.tolist())
ALL_SS_nums = [ca_name2num.get(n) for n in ALL_SS if ca_name2num.get(n)]

print(f"\nThreshold: >= 100 alerts AND >= 12 months of data")
print(f"ALL ShotSpotter neighborhoods ({len(ALL_SS)}):")
for n in ALL_SS:
    tag = " (original 6)" if n in ORIG_6 else ""
    print(f"  - {n}{tag}")

print(f"\nOriginal 6 study neighborhoods:")
for n in ORIG_6:
    print(f"  - {n}")

# ==================================================================
# SECTION 2 : TASK 2  --  DIFFERENCE-IN-DIFFERENCES
# ==================================================================
print("\n" + "=" * 72)
print("TASK 2: DIFFERENCE-IN-DIFFERENCES ANALYSIS")
print("=" * 72)

# --- Build full panel: 77 community areas x 15 years ---
all_ca_nums = sorted(ca_num2name.keys())
all_years   = list(range(2009, 2024))
panel = pd.DataFrame([(ca, yr) for ca in all_ca_nums for yr in all_years],
                      columns=["ca_num", "year"])
panel["ca_name"] = panel["ca_num"].map(ca_num2name)

# Gun homicide counts
gun_counts = (gun.dropna(subset=["CA_num"])
              .groupby(["CA_num", "year"]).size()
              .reset_index(name="gun_hom"))
gun_counts["CA_num"] = gun_counts["CA_num"].astype(int)
panel = panel.merge(gun_counts.rename(columns={"CA_num": "ca_num"}),
                     on=["ca_num", "year"], how="left")
panel["gun_hom"] = panel["gun_hom"].fillna(0).astype(int)

# Non-fatal shooting counts
nfs_counts = (nfs.dropna(subset=["CA_num"])
              .groupby(["CA_num", "year"]).size()
              .reset_index(name="nfs"))
nfs_counts["CA_num"] = nfs_counts["CA_num"].astype(int)
panel = panel.merge(nfs_counts.rename(columns={"CA_num": "ca_num"}),
                     on=["ca_num", "year"], how="left")
panel["nfs"] = panel["nfs"].fillna(0).astype(int)

# Treatment indicators
panel["treat_6"]   = panel["ca_name"].isin(ORIG_6).astype(int)
panel["treat_all"] = panel["ca_name"].isin(ALL_SS).astype(int)
panel["post"]      = (panel["year"] >= 2017).astype(int)
panel["treat6_post"]  = panel["treat_6"]  * panel["post"]
panel["treatall_post"] = panel["treat_all"] * panel["post"]

print(f"\nPanel: {len(panel)} observations "
      f"({len(all_ca_nums)} areas x {len(all_years)} years)")
print(f"Treatment (6):   {panel['treat_6'].sum()} area-years")
print(f"Treatment (all): {panel['treat_all'].sum()} area-years")

def did_table(treat_col, label, exclude_years=None):
    """Compute simple DiD table."""
    df = panel.copy()
    if exclude_years:
        df = df[~df["year"].isin(exclude_years)]
    pre  = df[df["post"] == 0]
    post = df[df["post"] == 1]
    n_pre_yrs  = pre["year"].nunique()
    n_post_yrs = post["year"].nunique()

    t_pre  = pre[pre[treat_col] == 1]["gun_hom"].sum() / n_pre_yrs
    t_post = post[post[treat_col] == 1]["gun_hom"].sum() / n_post_yrs
    c_pre  = pre[pre[treat_col] == 0]["gun_hom"].sum() / n_pre_yrs
    c_post = post[post[treat_col] == 0]["gun_hom"].sum() / n_post_yrs
    did    = (t_post - t_pre) - (c_post - c_pre)

    print(f"\n--- {label} ---")
    excl = f" (excl {exclude_years})" if exclude_years else ""
    print(f"  Pre years:  {n_pre_yrs},  Post years: {n_post_yrs}{excl}")
    print(f"  {'':20s} {'Pre avg/yr':>12s} {'Post avg/yr':>12s} {'Diff':>10s}")
    print(f"  {'Treatment':<20s} {t_pre:>12.1f} {t_post:>12.1f} {t_post-t_pre:>+10.1f}")
    print(f"  {'Control':<20s} {c_pre:>12.1f} {c_post:>12.1f} {c_post-c_pre:>+10.1f}")
    print(f"  {'DiD estimate':<20s} {'':>12s} {'':>12s} {did:>+10.1f}")
    return {"label": label, "t_pre": t_pre, "t_post": t_post,
            "c_pre": c_pre, "c_post": c_post, "did": did,
            "n_pre": n_pre_yrs, "n_post": n_post_yrs}

# 2A: 6 neighborhoods, full period
r2a = did_table("treat_6", "2A: Original 6 Neighborhoods")
# 2B: All SS neighborhoods, full period
r2b = did_table("treat_all", "2B: All ShotSpotter Neighborhoods")
# 2C: Excluding COVID
r2c_6   = did_table("treat_6",   "2C-6:  6 Neighborhoods, excl COVID",   [2020, 2021])
r2c_all = did_table("treat_all", "2C-all: All SS, excl COVID",           [2020, 2021])
# Extra: Excluding 2016 anomaly
r2_no16_6   = did_table("treat_6",   "Extra: 6 Hoods, excl 2016",        [2016])
r2_no16_all = did_table("treat_all", "Extra: All SS, excl 2016",          [2016])
# Extra: Excluding both COVID and 2016
r2_both_6   = did_table("treat_6",   "Extra: 6 Hoods, excl COVID+2016",  [2016, 2020, 2021])
r2_both_all = did_table("treat_all", "Extra: All SS, excl COVID+2016",    [2016, 2020, 2021])

# --- 2D: DiD Regression ---
print("\n--- 2D: DiD OLS Regression ---")

def run_did_regression(treat_col, interact_col, label, exclude_years=None):
    df = panel.copy()
    if exclude_years:
        df = df[~df["year"].isin(exclude_years)]
    formula = f"gun_hom ~ {treat_col} + post + {interact_col}"
    model = smf.ols(formula, data=df).fit(
        cov_type="cluster", cov_kwds={"groups": df["ca_num"]})
    b3    = model.params[interact_col]
    se3   = model.bse[interact_col]
    t3    = model.tvalues[interact_col]
    p3    = model.pvalues[interact_col]
    ci_lo, ci_hi = model.conf_int().loc[interact_col]
    n     = model.nobs
    r2    = model.rsquared

    print(f"\n  {label}")
    print(f"  Formula: {formula}")
    print(f"  N = {int(n)}, R-sq = {r2:.4f}")
    print(f"  beta3 (DiD) = {b3:+.3f}  (SE={se3:.3f}, t={t3:.3f}, p={p3:.4f})")
    print(f"  95% CI: [{ci_lo:.3f}, {ci_hi:.3f}]")
    sig = "***" if p3 < 0.01 else "**" if p3 < 0.05 else "*" if p3 < 0.1 else "ns"
    print(f"  Significance: {sig}")
    return {"label": label, "b3": b3, "se": se3, "t": t3, "p": p3,
            "ci_lo": ci_lo, "ci_hi": ci_hi, "n": int(n), "r2": r2}

reg_6   = run_did_regression("treat_6",   "treat6_post",
                              "6 Neighborhoods, full period")
reg_all = run_did_regression("treat_all", "treatall_post",
                              "All SS Neighborhoods, full period")
reg_6nc = run_did_regression("treat_6",   "treat6_post",
                              "6 Neighborhoods, excl COVID", [2020, 2021])
reg_allnc = run_did_regression("treat_all", "treatall_post",
                                "All SS Neighborhoods, excl COVID", [2020, 2021])
reg_6_no16 = run_did_regression("treat_6", "treat6_post",
                                 "6 Neighborhoods, excl 2016", [2016])
reg_all_no16 = run_did_regression("treat_all", "treatall_post",
                                   "All SS Neighborhoods, excl 2016", [2016])
reg_6_both = run_did_regression("treat_6", "treat6_post",
                                 "6 Neighborhoods, excl COVID+2016", [2016, 2020, 2021])
reg_all_both = run_did_regression("treat_all", "treatall_post",
                                   "All SS Neighborhoods, excl COVID+2016", [2016, 2020, 2021])

# --- Non-fatal shootings DiD regression ---
print("\n--- Non-Fatal Shootings DiD ---")
panel["treatall_post_nfs"] = panel["treat_all"] * panel["post"]
panel["treat6_post_nfs"]   = panel["treat_6"]   * panel["post"]

def run_nfs_regression(treat_col, interact_col, label, exclude_years=None):
    df = panel.copy()
    if exclude_years:
        df = df[~df["year"].isin(exclude_years)]
    formula = f"nfs ~ {treat_col} + post + {interact_col}"
    model = smf.ols(formula, data=df).fit(
        cov_type="cluster", cov_kwds={"groups": df["ca_num"]})
    b3 = model.params[interact_col]
    se3 = model.bse[interact_col]
    p3 = model.pvalues[interact_col]
    ci_lo, ci_hi = model.conf_int().loc[interact_col]
    print(f"  {label}: b3={b3:+.3f} (SE={se3:.3f}, p={p3:.4f}), "
          f"CI=[{ci_lo:.3f},{ci_hi:.3f}]")
    return {"label": label, "b3": b3, "se": se3, "p": p3,
            "ci_lo": ci_lo, "ci_hi": ci_hi, "n": int(model.nobs)}

nfs_reg_all = run_nfs_regression("treat_all", "treatall_post",
                                  "Non-fatal: All SS, full")
nfs_reg_6   = run_nfs_regression("treat_6", "treat6_post",
                                  "Non-fatal: 6 hoods, full")

# ==================================================================
# SECTION 3 : TASK 3  --  T-TESTS
# ==================================================================
print("\n" + "=" * 72)
print("TASK 3: STATISTICAL TESTING (T-TESTS)")
print("=" * 72)

# Clean data for 6 neighborhoods
hom6 = pd.read_csv(os.path.join(CLEAN, "Homicide_six_hoods.csv"))
hom6["Date"] = pd.to_datetime(hom6["Date"])

# Complete month index
all_months = pd.period_range("2009-01", "2023-12", freq="M")

def ttest_neighborhood(hood, hom_df):
    h = hom_df[hom_df["Community"] == hood].copy()
    h["month"] = h["Date"].dt.to_period("M")
    monthly = h.groupby("month").size()
    monthly = monthly.reindex(all_months, fill_value=0)

    # Jan 2017 (the first ShotSpotter install month) is excluded from both windows
    # as a one-month washout so the "before"/"after" contrast straddles a clean break.
    before = monthly[monthly.index < pd.Period("2017-01", "M")]
    after  = monthly[monthly.index >= pd.Period("2017-02", "M")]

    # Excl COVID
    after_nc = monthly[(monthly.index >= pd.Period("2017-02", "M")) &
                       ~((monthly.index >= pd.Period("2020-01", "M")) &
                         (monthly.index <= pd.Period("2021-12", "M")))]

    t_stat, p_val = stats.ttest_ind(before, after, equal_var=False)
    t_nc, p_nc    = stats.ttest_ind(before, after_nc, equal_var=False)

    se_diff = np.sqrt(before.var()/len(before) + after_nc.var()/len(after_nc))
    diff_nc = after_nc.mean() - before.mean()
    t_crit  = stats.t.ppf(0.975, min(len(before), len(after_nc)) - 1)
    ci_lo_diff = diff_nc - t_crit * se_diff
    ci_hi_diff = diff_nc + t_crit * se_diff

    bm = before.mean()
    pct = (after_nc.mean() - bm) / bm * 100 if bm > 0 else np.nan
    pct_ci_lo = ci_lo_diff / bm * 100 if bm > 0 else np.nan
    pct_ci_hi = ci_hi_diff / bm * 100 if bm > 0 else np.nan

    sig = "***" if p_val < 0.01 else "**" if p_val < 0.05 else "*" if p_val < 0.1 else "ns"

    return {
        "hood": hood,
        "before_mean": bm, "before_std": before.std(),
        "after_mean": after.mean(), "after_std": after.std(),
        "after_nc_mean": after_nc.mean(),
        "pct_full": (after.mean() - bm) / bm * 100 if bm > 0 else np.nan,
        "pct_nc": pct,
        "pct_ci_lo": pct_ci_lo, "pct_ci_hi": pct_ci_hi,
        "t": t_stat, "p": p_val, "sig": sig,
        "t_nc": t_nc, "p_nc": p_nc,
        "n_before": len(before), "n_after": len(after),
    }

# 3A: Per-neighborhood t-tests
print("\n--- 3A: Per-Neighborhood T-Tests ---")
print(f"{'Neighborhood':<20s} {'Before/mo':>10s} {'After/mo':>10s} "
      f"{'Change':>8s} {'t-stat':>8s} {'p-value':>8s} {'Sig':>4s}")
print("-" * 75)

ttest_results = []
for hood in ORIG_6:
    r = ttest_neighborhood(hood, hom6)
    ttest_results.append(r)
    print(f"  {r['hood']:<18s} {r['before_mean']:>10.3f} {r['after_mean']:>10.3f} "
          f"{r['pct_full']:>+7.1f}% {r['t']:>8.3f} {r['p']:>8.4f} {r['sig']:>4s}")

# 3B: Aggregated t-test (all 6 combined)
print("\n--- 3B: All 6 Combined ---")
combined = hom6.copy()
combined["month"] = combined["Date"].dt.to_period("M")
monthly_comb = combined.groupby("month").size().reindex(all_months, fill_value=0)
bef = monthly_comb[monthly_comb.index < pd.Period("2017-01", "M")]
aft = monthly_comb[monthly_comb.index >= pd.Period("2017-02", "M")]
t_c, p_c = stats.ttest_ind(bef, aft, equal_var=False)
pct_c = (aft.mean() - bef.mean()) / bef.mean() * 100
sig_c = "***" if p_c < 0.01 else "**" if p_c < 0.05 else "*" if p_c < 0.1 else "ns"
print(f"  Before: {bef.mean():.3f}/mo | After: {aft.mean():.3f}/mo | "
      f"Change: {pct_c:+.1f}% | t={t_c:.3f}, p={p_c:.4f} {sig_c}")

# ==================================================================
# SECTION 4 : TASK 4  --  PLOTS
# ==================================================================
print("\n" + "=" * 72)
print("TASK 4: PLOTS")
print("=" * 72)

# ----- Aggregate annual data for plotting -----
def annual_group(treat_col):
    """Annual totals for treatment vs. control."""
    t = panel[panel[treat_col] == 1].groupby("year")["gun_hom"].sum()
    c = panel[panel[treat_col] == 0].groupby("year")["gun_hom"].sum()
    return t, c

def annual_group_nfs(treat_col):
    t = panel[panel[treat_col] == 1].groupby("year")["nfs"].sum()
    c = panel[panel[treat_col] == 0].groupby("year")["nfs"].sum()
    return t, c

t_all, c_all = annual_group("treat_all")
t_6, c_6     = annual_group("treat_6")

# ---- PLOT 1 : DiD Parallel Trends (indexed to 100 at 2016) ----
print("\nPlot 1: DiD Parallel Trends")
fig, ax = plt.subplots(figsize=(12, 7))

# Index both series to 100 at 2016
base_t = t_all.loc[2016]
base_c = c_all.loc[2016]
idx_t = t_all / base_t * 100
idx_c = c_all / base_c * 100

ax.axvspan(2020, 2021, alpha=0.5, color=COLORS["covid_shade"],
           label="COVID (2020–21)", zorder=0)
ax.axhline(100, color="#cccccc", ls=":", lw=0.8, zorder=0)

ax.plot(idx_t.index, idx_t.values, "o-", color=C_BLUE, lw=2.2, markersize=6,
        label=f"ShotSpotter neighborhoods (n={len(ALL_SS)})", zorder=3)
ax.plot(idx_c.index, idx_c.values, "s-", color=C_RED, lw=2.2, markersize=6,
        label=f"Rest of Chicago (n={77 - len(ALL_SS)})", zorder=3)

# Single ShotSpotter Installed annotation
ax.axvline(2017, color="#333", ls="--", lw=1.2, zorder=1)
y_top = max(idx_t.max(), idx_c.max())
ax.text(2017.12, y_top * 0.98, "ShotSpotter\ninstalled (2017)",
        fontsize=10, color="#333", va="top", fontweight="bold")

# Annotate the 2016 spike — placed below the curves so it doesn't clash with the legend
ax.annotate("2016 baseline (= 100 by construction):\npre-treatment spike inflates it",
            xy=(2016, 100), xytext=(2012.5, 30),
            fontsize=9.5, color="#333", ha="center",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor=COLORS["grid"], linewidth=0.5),
            arrowprops=dict(arrowstyle="->", color=COLORS["text_muted"],
                            lw=0.8, connectionstyle="arc3,rad=-0.2"))

ax.set_xlabel("Year")
ax.set_ylabel("Gun homicides (indexed: 2016 = 100)")
ax.set_title("Difference-in-Differences: ShotSpotter vs. Rest of Chicago\n"
             "(Gun Homicides, Indexed to 2016 = 100)", fontweight="bold")
ax.set_xticks(range(2009, 2024))
ax.tick_params(axis="x", rotation=45)
ax.legend(loc="upper left", frameon=False, fontsize=10)
ax.yaxis.grid(True)

add_source_note(fig, SOURCE_DEFAULT)
fig.subplots_adjust(top=0.90, bottom=0.16, left=0.08, right=0.97)
savefig(fig, "did_parallel_trends.png")

# ---- PLOT 2 : Event Study ----
print("\nPlot 2: Event Study")

# Build event-study regression: Y_it = a_i + g_t + sum_k b_k*(Treat_i*D_k) + e
es_panel = panel.copy()
# Create interaction dummies (omit 2016 = base year)
for yr in all_years:
    if yr != 2016:
        es_panel[f"tx{yr}"] = (es_panel["treat_all"] * (es_panel["year"] == yr)).astype(int)

interact_terms = [f"tx{yr}" for yr in all_years if yr != 2016]
formula_es = "gun_hom ~ C(ca_num) + C(year) + " + " + ".join(interact_terms)
model_es = smf.ols(formula_es, data=es_panel).fit(
    cov_type="cluster", cov_kwds={"groups": es_panel["ca_num"]})

# Extract event-study coefficients
es_years = [yr for yr in all_years if yr != 2016]
es_coefs = [model_es.params[f"tx{yr}"] for yr in es_years]
es_ci    = model_es.conf_int()
es_lo    = [es_ci.loc[f"tx{yr}", 0] for yr in es_years]
es_hi    = [es_ci.loc[f"tx{yr}", 1] for yr in es_years]
es_rel   = [yr - 2017 for yr in es_years]   # relative to treatment year

# Insert 2016 = 0 (omitted base year)
rel_2016 = 2016 - 2017  # = -1
idx_insert = sum(1 for x in es_rel if x < rel_2016)
es_rel.insert(idx_insert, rel_2016)
es_coefs.insert(idx_insert, 0)
es_lo.insert(idx_insert, 0)
es_hi.insert(idx_insert, 0)

fig, ax = plt.subplots(figsize=(12, 6.5))

# Pre-period shading anchors the design visually
ax.axvspan(min(es_rel) - 1, -0.5, color=COLORS["pre_shade"], alpha=0.5, zorder=0)
ax.axvspan(3, 4, alpha=0.5, color=COLORS["covid_shade"], zorder=0,
           label="COVID (2020–21)")
ax.axhline(0, color="#333333", lw=1.2, zorder=1)
ax.axvline(-0.5, color=C_GRAY, ls="--", lw=1.2, zorder=1)

for i, (yr_r, coef, lo, hi) in enumerate(zip(es_rel, es_coefs, es_lo, es_hi)):
    color = C_BLUE if yr_r < 0 else C_RED
    ax.plot([yr_r, yr_r], [lo, hi], color=color, lw=2.2, zorder=2)
    ax.plot(yr_r, coef, "o", color=color, markersize=8, zorder=3,
            markeredgecolor="white", markeredgewidth=0.6)

ax.set_xlabel("Years relative to ShotSpotter installation  (k = 0 ≡ 2017)")
ax.set_ylabel("Estimated treatment effect\n(Gun homicides per community area)")
ax.set_title("Event Study: Year-by-Year Treatment Effects\n"
             "(Base Year = 2016, Clustered SE at Community Area Level)",
             fontweight="bold")
ax.set_xticks(es_rel)
ax.set_xticklabels([f"{r:+d}" if r != 0 else "0" for r in es_rel])
ax.yaxis.grid(True)

# Custom legend with pre/post coding instead of cramming it elsewhere
legend_handles = [
    Line2D([0], [0], color=C_BLUE, marker="o", lw=2.2, markersize=7,
           label="Pre-treatment (k < 0)"),
    Line2D([0], [0], color=C_RED, marker="o", lw=2.2, markersize=7,
           label="Post-treatment (k ≥ 0)"),
    mpatches.Patch(color=COLORS["covid_shade"], alpha=0.5, label="COVID (2020–21)"),
]
ax.legend(handles=legend_handles, loc="upper left", frameon=False, fontsize=10)

add_source_note(fig, SOURCE_DEFAULT,
                sample_n=f"N = 1,155 community-area-years; treated n={len(ALL_SS)}")

fig.subplots_adjust(top=0.88, bottom=0.13, left=0.10, right=0.97)
savefig(fig, "did_event_study.png")

# ---- PLOT 3 : Annual Trends by Neighborhood ----
print("\nPlot 3: Annual Trends by Neighborhood")
fig, axes = plt.subplots(2, 3, figsize=(15, 10), sharex=True, sharey=True)
axes = axes.flatten()

# Compute global y-max so all panels are visually comparable
yearly_by_hood = {}
for hood in ORIG_6:
    h = hom6[hom6["Community"] == hood]
    yearly_by_hood[hood] = (h.groupby(h["Date"].dt.year).size()
                            .reindex(range(2009, 2024), fill_value=0))
gmax = max(y.max() for y in yearly_by_hood.values())

for idx, hood in enumerate(ORIG_6):
    ax = axes[idx]
    yearly = yearly_by_hood[hood]

    ax.axvspan(2020, 2021, alpha=0.5, color=COLORS["covid_shade"], zorder=0)
    ax.axvline(2017, color="#333", ls="--", lw=1.2, alpha=0.85, zorder=1)

    ax.bar(yearly.index, yearly.values, color=C_BLUE, alpha=0.75, width=0.7, zorder=2)
    ax.plot(yearly.index, yearly.values, "o-", color=C_BLUE, lw=1.6,
            markersize=4, zorder=3)

    # Pre and post means as horizontal reference lines
    pre_mean  = yearly.loc[2009:2016].mean()
    post_mean = yearly.loc[2017:2023].mean()
    ax.axhline(pre_mean, color=C_BLUE, ls=":", lw=1.2, alpha=0.7, xmin=0, xmax=0.53)
    ax.axhline(post_mean, color=C_RED, ls=":", lw=1.2, alpha=0.7, xmin=0.53, xmax=1)

    ax.set_title(f"{hood}    (pre: {pre_mean:.1f}/yr  →  post: {post_mean:.1f}/yr)",
                 fontsize=11, fontweight="bold")
    if idx % 3 == 0:
        ax.set_ylabel("Gun homicides")
    ax.set_xticks(range(2009, 2024, 2))
    ax.tick_params(axis="x", rotation=45, labelsize=9)
    ax.yaxis.grid(True)
    ax.set_ylim(0, gmax * 1.1)

fig.suptitle("Annual Gun Homicide Counts by Neighborhood (2009–2023)",
             fontsize=14, fontweight="bold", y=0.99)

legend_elements = [
    Line2D([0], [0], color="#333", ls="--", lw=1.2, label="ShotSpotter installed (2017)"),
    mpatches.Patch(facecolor=COLORS["covid_shade"], alpha=0.5, label="COVID period (2020–21)"),
    Line2D([0], [0], color=C_BLUE, ls=":", lw=1.2, label="Pre-treatment mean"),
    Line2D([0], [0], color=C_RED, ls=":", lw=1.2, label="Post-treatment mean"),
]
fig.legend(handles=legend_elements, loc="lower center", ncol=4,
           frameon=False, fontsize=10, bbox_to_anchor=(0.5, 0.01))

add_source_note(fig, SOURCE_DEFAULT, y=0.06, x=0.5)
for txt in fig.texts:
    if txt.get_text().startswith("Source:"):
        txt.set_ha("center")

fig.tight_layout(rect=[0, 0.07, 1, 0.96])
savefig(fig, "annual_trends_neighborhoods.png")

# ---- PLOT 4 : Forest Plot ----
print("\nPlot 4: Forest Plot of Treatment Effects")
fig, ax = plt.subplots(figsize=(10, 7))

# Data: pct change excl COVID, with CI
forest_data = []
for r in ttest_results:
    forest_data.append({
        "name": r["hood"], "pct": r["pct_nc"],
        "lo": r["pct_ci_lo"], "hi": r["pct_ci_hi"],
    })

# Add DiD estimate (all SS, excl COVID) — convert to per-area annualized pct
# Use the regression coefficient relative to pre-treatment mean
n_treat = len(ALL_SS)
pre_mean_all = panel[(panel["treat_all"] == 1) & (panel["post"] == 0)]["gun_hom"].mean()
did_pct = reg_allnc["b3"] / pre_mean_all * 100 if pre_mean_all > 0 else 0
did_lo  = reg_allnc["ci_lo"] / pre_mean_all * 100 if pre_mean_all > 0 else 0
did_hi  = reg_allnc["ci_hi"] / pre_mean_all * 100 if pre_mean_all > 0 else 0
forest_data.append({"name": "DiD Estimate\n(All SS, excl COVID)",
                     "pct": did_pct, "lo": did_lo, "hi": did_hi})

n = len(forest_data)
y_pos = list(range(n))

# Compute a sensible right-side x-limit so percentage callouts don't overflow
max_hi = max(d["hi"] for d in forest_data)
xlim_right = max_hi + 18

for i, d in enumerate(forest_data):
    # red = CI excludes zero (significant increase); neutral gray = CI spans zero.
    # Gray rather than green: red/green is the classic colorblind-unsafe pairing.
    color = COLORS["neutral"] if d["lo"] < 0 else C_RED
    if i == n - 1:
        color = "#333333"
        ax.axhline(i, color="#eeeeee", lw=8, zorder=0)

    ax.plot([d["lo"], d["hi"]], [i, i], color=color, lw=2.5, zorder=2)
    ax.plot(d["pct"], i, "o", color=color, markersize=8, zorder=3,
            markeredgecolor="white", markeredgewidth=0.8)
    ax.text(d["hi"] + 2, i, f"{d['pct']:+.1f}%", va="center", fontsize=10,
            fontweight="bold", color=color)

ax.axvline(0, color="#333333", ls="-", lw=1, zorder=1)
ax.set_yticks(y_pos)
ax.set_yticklabels([d["name"] for d in forest_data], fontsize=10)
ax.invert_yaxis()
ax.set_xlim(min(d["lo"] for d in forest_data) - 5, xlim_right)
ax.set_xlabel("Percent change in gun homicide rate\n(Before vs. after ShotSpotter, excl. COVID)")
ax.set_title("Forest Plot: Treatment Effect Estimates by Neighborhood",
             fontweight="bold", pad=12)
ax.xaxis.grid(True)
ax.spines["left"].set_visible(False)
ax.tick_params(axis="y", length=0)
add_source_note(fig, SOURCE_DEFAULT,
                sample_n="t-test on monthly counts; DiD = OLS w/ two-way FE")
fig.subplots_adjust(top=0.92, bottom=0.16, left=0.18, right=0.97)
savefig(fig, "forest_plot_effects.png")

# ---- PLOT 5 : Treatment vs. Control Map ----
print("\nPlot 5: Treatment vs. Control Map")

# Since 51/77 community areas have ShotSpotter, the "control" group is the
# 26 areas WITHOUT ShotSpotter. Find the highest-hardship ones among these.
all_ca_names = [ca_num2name[n] for n in all_ca_nums]
NON_SS = sorted(set(all_ca_names) - set(ALL_SS))
print(f"  Non-ShotSpotter areas ({len(NON_SS)}): {NON_SS}")

# Pick highest-hardship non-SS areas as "highlighted controls"
non_ss_hardship = []
for name in NON_SS:
    crow = census[census["CA_name"].str.lower() == name.lower()]
    hi = crow["HARDSHIP INDEX"].values[0] if len(crow) > 0 else 0
    non_ss_hardship.append((name, hi))
non_ss_hardship.sort(key=lambda x: -x[1])
CONTROLS = [n for n, h in non_ss_hardship[:6]]
print(f"  Highest-hardship non-SS controls: {CONTROLS}")

gdf_plot = gdf.copy()
gdf_plot["ca_name"] = gdf_plot["community"].str.title()

# Highlight: original 6 in dark blue, other SS in light blue, controls in red
def classify(name):
    if name in ORIG_6:
        return "Original 6 (ShotSpotter)"
    elif name in ALL_SS:
        return "Other ShotSpotter"
    elif name in CONTROLS:
        return "Control (highest-hardship non-SS)"
    else:
        return "Other (no ShotSpotter)"

gdf_plot["group"] = gdf_plot["ca_name"].apply(classify)
color_map = {
    "Original 6 (ShotSpotter)": C_BLUE,
    "Other ShotSpotter": "#A6C8E0",
    "Control (highest-hardship non-SS)": C_RED,
    "Other (no ShotSpotter)": "#E8E8E8",
}

fig, ax = plt.subplots(figsize=(9, 11))
# Plot "Other" first so it's the background
for group in ["Other (no ShotSpotter)", "Other ShotSpotter",
              "Control (highest-hardship non-SS)", "Original 6 (ShotSpotter)"]:
    subset = gdf_plot[gdf_plot["group"] == group]
    if len(subset) > 0:
        subset.plot(ax=ax, color=color_map[group],
                    edgecolor="#666666", linewidth=0.4)

# Labels for original 6 and controls only. Adjacent centroids collide for a few
# neighbor pairs (West Englewood/Englewood), so nudge those labels apart.
LABEL_OFFSETS = {
    "West Englewood": (-0.022, -0.006),
    "Englewood": (0.016, 0.006),
    "Irving Park": (0.010, -0.004),
    "Portage Park": (-0.010, 0.004),
}
label_groups = ["Original 6 (ShotSpotter)", "Control (highest-hardship non-SS)"]
for _, row in gdf_plot[gdf_plot["group"].isin(label_groups)].iterrows():
    c = row.geometry.centroid
    dx, dy = LABEL_OFFSETS.get(row["ca_name"], (0.0, 0.0))
    ax.text(c.x + dx, c.y + dy, row["ca_name"], fontsize=6.5, ha="center", va="center",
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.18", fc="white", alpha=0.85, ec="none"))

ax.axis("off")
ax.set_title("Chicago Community Areas:\nShotSpotter Coverage vs. Control Areas",
             fontsize=15, fontweight="bold", pad=12)

legend_elements = [
    mpatches.Patch(color=C_BLUE, label=f"Original 6 study areas"),
    mpatches.Patch(color="#A6C8E0", label=f"Other ShotSpotter ({len(ALL_SS)-6} areas)"),
    mpatches.Patch(color=C_RED, label=f"Control: non-SS ({len(CONTROLS)} highest hardship)"),
    mpatches.Patch(color="#E8E8E8", label=f"Other non-SS ({len(NON_SS)-len(CONTROLS)} areas)"),
]
ax.legend(handles=legend_elements, loc="lower left", frameon=True,
          framealpha=0.95, edgecolor=COLORS["grid"], fontsize=10)
add_source_note(fig, "Source: City of Chicago Community Areas (geojson, 2024).")
fig.tight_layout()
savefig(fig, "did_map_treatment_control.png")

# ---- PLOT 6 : Robustness Summary ----
print("\nPlot 6: Robustness Summary")

all_regs = [reg_6, reg_all, reg_6nc, reg_allnc,
            reg_6_no16, reg_all_no16, reg_6_both, reg_all_both,
            nfs_reg_6, nfs_reg_all]

fig, ax = plt.subplots(figsize=(12, 7.5))
y_pos = list(range(len(all_regs)))

# Compute right margin to fit value labels (largest CI hi + headroom)
max_hi = max(r["ci_hi"] for r in all_regs)
xlim_right = max_hi + 8

# Light gray null band for visual emphasis
ax.axvspan(-0.5, 0.5, color=COLORS["grid"], alpha=0.5, zorder=0)

for i, r in enumerate(all_regs):
    is_nfs = "Non-fatal" in r["label"]
    color = C_GREEN if r["ci_lo"] < 0 and r["ci_hi"] > 0 else (
            C_BLUE if r["b3"] < 0 else C_RED)
    if is_nfs:
        color = COLORS["purple"]

    ax.plot([r["ci_lo"], r["ci_hi"]], [i, i], color=color, lw=2.5, zorder=2)
    ax.plot(r["b3"], i, "D" if is_nfs else "o", color=color, markersize=8,
            zorder=3, markeredgecolor="white", markeredgewidth=0.6)
    sig = sig_stars(r["p"])
    ax.text(r["ci_hi"] + 0.4, i, f"{r['b3']:+.2f}{sig}",
            va="center", fontsize=9, color=color, fontweight="bold")

ax.axvline(0, color="#333333", ls="-", lw=1.2, zorder=1)
ax.set_yticks(y_pos)
ax.set_yticklabels([r["label"] for r in all_regs], fontsize=10)
ax.invert_yaxis()
ax.set_xlim(min(r["ci_lo"] for r in all_regs) - 2, xlim_right)
ax.set_xlabel("DiD estimate (β₃: Treat × Post coefficient)\n"
              "Negative = fewer homicides in ShotSpotter areas relative to control")
ax.set_title("Robustness of DiD Estimates Across Specifications",
             fontweight="bold", pad=12)
ax.xaxis.grid(True)
ax.spines["left"].set_visible(False)
ax.tick_params(axis="y", length=0)

add_source_note(fig, SOURCE_DEFAULT,
                sample_n="Significance: *** p<0.01, ** p<0.05, * p<0.1")
fig.subplots_adjust(top=0.92, bottom=0.16, left=0.28, right=0.97)
savefig(fig, "robustness_summary.png")

# ---- PLOT 7 : Non-Fatal Shootings DiD Trends ----
print("\nPlot 7: Non-Fatal Shootings DiD Trends")
t_nfs, c_nfs = annual_group_nfs("treat_all")

fig, ax = plt.subplots(figsize=(12, 7))

base_t_nfs = t_nfs.loc[2016] if t_nfs.loc[2016] > 0 else 1
base_c_nfs = c_nfs.loc[2016] if c_nfs.loc[2016] > 0 else 1
idx_t_nfs = t_nfs / base_t_nfs * 100
idx_c_nfs = c_nfs / base_c_nfs * 100

ax.axvspan(2020, 2021, alpha=0.5, color=COLORS["covid_shade"],
           label="COVID (2020–21)", zorder=0)
ax.axhline(100, color="#cccccc", ls=":", lw=0.8, zorder=0)

ax.plot(idx_t_nfs.index, idx_t_nfs.values, "o-", color=C_BLUE, lw=2.2,
        markersize=6, label=f"ShotSpotter neighborhoods (n={len(ALL_SS)})", zorder=3)
ax.plot(idx_c_nfs.index, idx_c_nfs.values, "s-", color=C_RED, lw=2.2,
        markersize=6, label=f"Rest of Chicago (n={77 - len(ALL_SS)})", zorder=3)

ax.axvline(2017, color="#333", ls="--", lw=1.2, zorder=1)
y_top = max(idx_t_nfs.max(), idx_c_nfs.max())
ax.text(2017.12, y_top * 0.98, "ShotSpotter\ninstalled (2017)",
        fontsize=10, color="#333", va="top", fontweight="bold")

ax.set_xlabel("Year")
ax.set_ylabel("Non-fatal shootings (indexed: 2016 = 100)")
ax.set_title("Difference-in-Differences: Non-Fatal Gunshot Injuries\n"
             "(ShotSpotter vs. Rest of Chicago, Indexed to 2016 = 100)",
             fontweight="bold")
ax.set_xticks(range(2009, 2024))
ax.tick_params(axis="x", rotation=45)
ax.legend(loc="upper left", frameon=False, fontsize=10)
ax.yaxis.grid(True)

add_source_note(fig, SOURCE_DEFAULT)
fig.subplots_adjust(top=0.90, bottom=0.16, left=0.08, right=0.97)
savefig(fig, "nonfatal_did_trends.png")

# ==================================================================
# SECTION 5 : TASK 5  --  SUMMARY TABLES
# ==================================================================
print("\n" + "=" * 72)
print("TASK 5: SUMMARY TABLES")
print("=" * 72)

lines = []
lines.append("# Results Summary")
lines.append("")
lines.append("**ACE 592 Final Project: ShotSpotter & Gun Homicide DiD Analysis**")
lines.append("")

# 5A: Master Results Table
lines.append("## DiD Regression Results")
lines.append("")
lines.append("| Specification | DiD (beta3) | Std Error | t-stat | p-value | 95% CI | N | R-sq |")
lines.append("|---|---|---|---|---|---|---|---|")
for r in [reg_6, reg_all, reg_6nc, reg_allnc, reg_6_no16, reg_all_no16,
          reg_6_both, reg_all_both, nfs_reg_6, nfs_reg_all]:
    sig = "***" if r["p"] < 0.01 else "**" if r["p"] < 0.05 else "*" if r["p"] < 0.1 else ""
    r2_str = f'{r.get("r2", 0):.4f}' if "r2" in r else "---"
    lines.append(
        f'| {r["label"]} | {r["b3"]:+.3f} | {r["se"]:.3f} | '
        f'{r.get("t", r["b3"]/r["se"]):.3f} | {r["p"]:.4f}{sig} | '
        f'[{r["ci_lo"]:.3f}, {r["ci_hi"]:.3f}] | {r["n"]:,} | {r2_str} |')
lines.append("")

# Simple DiD tables
lines.append("## Simple DiD Estimates (Annualized)")
lines.append("")
lines.append("| Specification | Treatment Pre | Treatment Post | Control Pre | Control Post | DiD |")
lines.append("|---|---|---|---|---|---|")
for r in [r2a, r2b, r2c_6, r2c_all, r2_no16_6, r2_no16_all, r2_both_6, r2_both_all]:
    lines.append(
        f'| {r["label"]} | {r["t_pre"]:.1f} | {r["t_post"]:.1f} | '
        f'{r["c_pre"]:.1f} | {r["c_post"]:.1f} | {r["did"]:+.1f} |')
lines.append("")

# Interpretation
lines.append("### Interpretation")
lines.append("")
lines.append(f"The DiD estimate using all {len(ALL_SS)} ShotSpotter neighborhoods "
             f"vs. the remaining {77 - len(ALL_SS)} community areas is "
             f"**{reg_all['b3']:+.3f}** gun homicides per community area per year "
             f"(p={reg_all['p']:.4f}).")
lines.append("")
if reg_all["b3"] < 0:
    lines.append("A **negative** coefficient means ShotSpotter neighborhoods experienced "
                 "a relative **decrease** in gun homicides compared to the rest of Chicago. "
                 "While both groups saw increases after 2017, the increase was smaller "
                 "in ShotSpotter areas.")
else:
    lines.append("A **positive** coefficient means ShotSpotter neighborhoods experienced "
                 "a relative **increase** in gun homicides compared to the rest of Chicago.")
lines.append("")

# 5B: T-Test Table
lines.append("## T-Test Results: Before vs. After ShotSpotter (Monthly Rates)")
lines.append("")
lines.append("| Neighborhood | Before (mean/mo) | After (mean/mo) | Change | t-stat | p-value | Sig |")
lines.append("|---|---|---|---|---|---|---|")
for r in ttest_results:
    lines.append(
        f'| {r["hood"]} | {r["before_mean"]:.3f} | {r["after_mean"]:.3f} | '
        f'{r["pct_full"]:+.1f}% | {r["t"]:.3f} | {r["p"]:.4f} | {r["sig"]} |')
lines.append(
    f'| **All 6 Combined** | {bef.mean():.3f} | {aft.mean():.3f} | '
    f'{pct_c:+.1f}% | {t_c:.3f} | {p_c:.4f} | {sig_c} |')
lines.append("")
lines.append("Note: Monthly counts include zero-homicide months. "
             "Welch's t-test (unequal variances). "
             "Before = Jan 2009 - Dec 2016; After = Feb 2017 - Dec 2023.")
lines.append("")

# ShotSpotter neighborhoods table
lines.append("## ShotSpotter Coverage: All Identified Neighborhoods")
lines.append("")
lines.append("| Community Area | Total Alerts | Coverage (months) | Earliest Alert | In Original 6? |")
lines.append("|---|---|---|---|---|")
for name in ALL_SS:
    row = ss_profile.loc[name]
    orig = "Yes" if name in ORIG_6 else "No"
    lines.append(
        f'| {name} | {row.total_alerts:,} | {row.months} | '
        f'{row.earliest.strftime("%Y-%m-%d")} | {orig} |')
lines.append("")

# Event study coefficients
lines.append("## Event Study Coefficients (Base Year = 2016)")
lines.append("")
lines.append("| Year | Relative | Coefficient | 95% CI Low | 95% CI High |")
lines.append("|---|---|---|---|---|")
for yr_r, coef, lo, hi in zip(es_rel, es_coefs, es_lo, es_hi):
    actual_yr = yr_r + 2017
    base = " (base)" if yr_r == -1 else ""
    lines.append(f'| {actual_yr}{base} | {yr_r:+d} | {coef:+.3f} | {lo:.3f} | {hi:.3f} |')
lines.append("")

# Write to file
report_path = os.path.join(P, "docs", "results_summary_original.md")
with open(report_path, "w", encoding="utf-8") as f:
    f.write("\n".join(lines))
print(f"\nWrote: docs/results_summary_original.md ({os.path.getsize(report_path)/1024:.0f} KB)")

print("\n" + "=" * 72)
print("ALL TASKS COMPLETE")
print("=" * 72)
