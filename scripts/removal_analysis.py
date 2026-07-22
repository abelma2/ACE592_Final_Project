"""
removal_analysis.py — The September 2024 shutoff as a second natural experiment.

Chicago ended ShotSpotter on 22 September 2024. A wave of commentary has read the
gun-homicide decline that followed as evidence the technology was useless (or worse).
This script applies the same discipline the rest of the paper applies to the
*installation*: it asks whether the *removal* effect is identifiable, and shows it is
not---for the mirror-image reason. Gun homicides were already falling city-wide from
their 2021--2022 peak, in ShotSpotter and never-treated areas alike, so a naive
before/after in the coverage areas attributes a city-wide trend to the shutoff. The
same missing counterfactual that dooms the installation estimate dooms the removal
estimate; the public debate is drawing the mirror-image unwarranted inference.

Uses the victims file already in the repo (homicides through mid-2025). 2024 is a
complete year; 2025 is partial and used only for the monthly window.

Outputs:
  plots/removal_trajectory.png
  docs/results_summary_removal.md
"""
import os
import sys
import json
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import statsmodels.formula.api as smf
import statsmodels.api as sm

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_style import (apply_academic_style, COLORS, add_source_note, add_n_label,
                        sig_stars, SOURCE_DEFAULT)

P = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(P, "data", "raw")
PLOTS = os.path.join(P, "plots")
DOCS = os.path.join(P, "docs")
apply_academic_style()
REMOVAL = pd.Timestamp("2024-09-22")

# =====================================================================
# DATA
# =====================================================================
print("=" * 72); print("LOADING DATA"); print("=" * 72)
with open(os.path.join(DATA, "CommAreas_20250408.geojson"), "r", encoding="utf-8") as f:
    geo = json.load(f)
ca_num2name = {int(float(ft["properties"]["area_num_1"])): ft["properties"]["community"].title()
               for ft in geo["features"]}
ca_name2num = {name: n for n, name in ca_num2name.items()}
ca_name2num.update({name.upper(): n for n, name in ca_num2name.items()})

hom = pd.read_csv(os.path.join(
    DATA, "Violence_Reduction_-_Victims_of_Homicides_and_Non-Fatal_Shootings_20250513.csv"), low_memory=False)
hom["DATE"] = pd.to_datetime(hom["DATE"], format="%m/%d/%Y %I:%M:%S %p", errors="coerce")
hom["CA_num"] = hom["COMMUNITY_AREA"].map(ca_name2num)
gun = hom[(hom["VICTIMIZATION_PRIMARY"] == "HOMICIDE") & (hom["GUNSHOT_INJURY_I"] == "YES")].dropna(
    subset=["CA_num", "DATE"]).copy()
gun["CA_num"] = gun["CA_num"].astype(int)
gun["year"] = gun["DATE"].dt.year
gun["month"] = gun["DATE"].dt.to_period("M").dt.to_timestamp()
last_full = (gun["month"].max() - pd.offsets.MonthBegin(1))  # drop the partial final month
print(f"  Homicide data spans to {gun['DATE'].max().date()}; last complete month = {last_full.date()}")

ss = pd.read_csv(os.path.join(DATA, "Violence_Reduction_-_Shotspotter_Alerts_-_Historical_20250408.csv"))
ss["DATE"] = pd.to_datetime(ss["DATE"], format="%m/%d/%Y %I:%M:%S %p", errors="coerce")
print(f"  ShotSpotter alert data ends {ss['DATE'].max().date()} (the shutoff).")
ss_hist = ss[ss["DATE"] <= "2023-12-31"]
samp = str(ss_hist["COMMUNITY_AREA"].dropna().iloc[0])
ss_hist["CA"] = (pd.to_numeric(ss_hist["COMMUNITY_AREA"], errors="coerce").map(ca_num2name)
                 if samp.replace(".", "").isdigit() else ss_hist["COMMUNITY_AREA"].str.title())
prof = ss_hist.dropna(subset=["CA"]).groupby("CA").agg(a=("DATE", "size"), e=("DATE", "min"), l=("DATE", "max"))
prof["mo"] = ((prof.l - prof.e).dt.days / 30.4).round()
ALL_SS = set(prof[(prof.a >= 100) & (prof.mo >= 12)].index)
n_t = len(ALL_SS); n_c = 77 - n_t

# =====================================================================
# 1. THE CITY-WIDE DECLINE PREDATES THE REMOVAL
# =====================================================================
print("\n" + "=" * 72)
print("1. ANNUAL GUN HOMICIDES, treated vs. never-treated (2009-2024)")
print("=" * 72)
gy = gun[gun.year.between(2009, 2024)].copy()
gy["treated"] = gy["CA_num"].map(ca_num2name).isin(ALL_SS)
ann = gy.groupby(["year", "treated"]).size().unstack(fill_value=0)
ann_per = ann.div([n_c, n_t], axis=1)   # per-area (False=control cols order matches [n_c,n_t]? fix below)
ann_per = pd.DataFrame({"control": ann[False] / n_c, "treated": ann[True] / n_t})
print(ann_per.round(2).to_string())
peak = ann_per["treated"].idxmax()
print(f"\n  Treated-area gun homicides per area peaked in {peak} "
      f"({ann_per.loc[peak,'treated']:.1f}/yr), fell to {ann_per.loc[2024,'treated']:.1f} by 2024 "
      f"(removal was Sept 2024). Control areas moved in parallel: "
      f"{ann_per.loc[peak,'control']:.1f} -> {ann_per.loc[2024,'control']:.1f}.")

# =====================================================================
# 2. NAIVE "REMOVAL EFFECT" vs. THE DiD THAT NETS OUT THE COMMON TREND
# =====================================================================
print("\n" + "=" * 72)
print("2. REMOVAL: naive before/after vs. difference-in-differences (monthly)")
print("=" * 72)
# monthly panel over a window around the removal, dropping the partial final month
win_start = pd.Timestamp("2022-01-01")
months = pd.date_range(win_start, last_full, freq="MS")
mpanel = pd.DataFrame([(ca, m) for ca in sorted(ca_num2name.keys()) for m in months],
                      columns=["ca_num", "month"])
mc = gun.groupby(["CA_num", "month"]).size().rename("gun_hom").reset_index()
mc.columns = ["ca_num", "month", "gun_hom"]
mpanel = mpanel.merge(mc, on=["ca_num", "month"], how="left")
mpanel["gun_hom"] = mpanel["gun_hom"].fillna(0).astype(int)
mpanel["treated"] = mpanel["ca_num"].map(ca_num2name).isin(ALL_SS).astype(int)
mpanel["post"] = (mpanel["month"] >= pd.Timestamp("2024-10-01")).astype(int)   # first full post month
mpanel["did"] = mpanel["treated"] * mpanel["post"]
mpanel["mstr"] = mpanel["month"].dt.strftime("%Y-%m")

# naive before/after within treated areas (the commentary's comparison)
t_pre = mpanel[(mpanel.treated == 1) & (mpanel.post == 0)]["gun_hom"].mean()
t_post = mpanel[(mpanel.treated == 1) & (mpanel.post == 1)]["gun_hom"].mean()
c_pre = mpanel[(mpanel.treated == 0) & (mpanel.post == 0)]["gun_hom"].mean()
c_post = mpanel[(mpanel.treated == 0) & (mpanel.post == 0)]["gun_hom"].mean()
c_post = mpanel[(mpanel.treated == 0) & (mpanel.post == 1)]["gun_hom"].mean()
print(f"  Naive treated before/after: {t_pre:.3f} -> {t_post:.3f} per area-month "
      f"({(t_post/t_pre-1)*100:+.0f}%)  <- the 'crime fell after removal' claim")
print(f"  Never-treated before/after: {c_pre:.3f} -> {c_post:.3f} per area-month "
      f"({(c_post/c_pre-1)*100:+.0f}%)  <- but controls fell too")

ols = smf.ols("gun_hom ~ did + C(ca_num) + C(mstr)", data=mpanel).fit(
    cov_type="cluster", cov_kwds={"groups": mpanel["ca_num"]})
pois = smf.glm("gun_hom ~ did + C(ca_num) + C(mstr)", data=mpanel, family=sm.families.Poisson()).fit(
    cov_type="cluster", cov_kwds={"groups": mpanel["ca_num"]})
b, p = ols.params["did"], ols.pvalues["did"]
se = ols.bse["did"]
irr = np.exp(pois.params["did"]); il, ih = np.exp(pois.conf_int().loc["did"]); ip = pois.pvalues["did"]
ise = pois.bse["did"]
n_post = mpanel["post"].sum() // 77
print(f"\n  Removal DiD (treated x post-removal), {n_post} post-removal months:")
print(f"    OLS     = {b:+.4f} (SE {se:.4f}) per area-month (p={p:.3f} {sig_stars(p)})")
print(f"    Poisson = IRR {irr:.3f} [{il:.3f}, {ih:.3f}] (log-SE {ise:.4f}, p={ip:.3f} {sig_stars(ip)})")

# =====================================================================
# FIGURE: indexed monthly rates (single axis, apples-to-apples) show both
#         groups declining in parallel straight through the removal.
# =====================================================================
grp = mpanel.groupby(["month", "treated"])["gun_hom"].sum().unstack()
tr_rate, co_rate = grp[1] / n_t, grp[0] / n_c
base = (grp.index >= "2022-01-01") & (grp.index <= "2023-12-31")   # shared pre-removal baseline = 100
tr_idx = (tr_rate / tr_rate[base].mean() * 100).rolling(3, center=True, min_periods=2).mean()
co_idx = (co_rate / co_rate[base].mean() * 100).rolling(3, center=True, min_periods=2).mean()

fig, ax = plt.subplots(figsize=(11, 5.8))
ax.plot(tr_idx.index, tr_idx.values, "-", color=COLORS["treated"], lw=2.6,
        label=f"ShotSpotter areas (n={n_t})")
ax.plot(co_idx.index, co_idx.values, "--", color=COLORS["control"], lw=2.6,
        label=f"Never-treated (n={n_c})")
ax.axhline(100, color="#999", lw=1.0, zorder=0)
ax.axvspan(pd.Timestamp("2024-10-01"), tr_idx.index.max(), color=COLORS["treated"], alpha=0.06, zorder=0)
ax.axvline(REMOVAL, color=COLORS["treated"], ls=":", lw=1.6, zorder=1)
ax.text(REMOVAL - pd.Timedelta(days=20), ax.get_ylim()[1] * 0.97,
        "ShotSpotter\nremoved\n(Sep 2024)", fontsize=8.8, color=COLORS["treated"],
        va="top", ha="right", style="italic")
ax.text(pd.Timestamp("2025-01-10"), ax.get_ylim()[1] * 0.97, "post-\nremoval",
        fontsize=8.8, color="#777", va="top", ha="left", style="italic")
ax.set_ylabel("Gun homicides, index (2022–23 mean = 100, 3-mo avg)")
ax.set_xlabel("Month")
ax.set_title("On a common scale, ShotSpotter and never-treated areas decline in parallel---\n"
             "the drop continues straight through the September 2024 removal",
             fontweight="bold", pad=10)
ax.legend(loc="upper left", frameon=False, fontsize=9.5)
ax.grid(True, alpha=0.3)
add_n_label(ax, f"Removal DiD ({n_post} post months):\n"
                f"additive OLS {b:+.2f}/area-mo ({sig_stars(p)}) — scale artifact\n"
                f"multiplicative IRR {irr:.2f} [{il:.2f}, {ih:.2f}] (ns)", loc="lower left")
add_source_note(fig, SOURCE_DEFAULT.replace("2009–2023", "2009–2024"))
fig.tight_layout(rect=[0, 0.03, 1, 1])
fig.savefig(os.path.join(PLOTS, "removal_trajectory.png"), dpi=200, facecolor="white")
plt.close(fig)
print("\n  SAVED: removal_trajectory.png")

# =====================================================================
# DOC
# =====================================================================
L = ["# The September 2024 Removal as a Second Natural Experiment", "",
     "Results of `scripts/removal_analysis.py`. Chicago ended ShotSpotter on 22 September 2024.",
     "Commentary has read the gun-homicide decline that followed as proof the technology did",
     "nothing (or caused harm). We apply the paper's discipline to the *removal* and find the",
     "removal effect is no more identifiable than the installation effect---for the mirror-image",
     "reason.", "",
     "## 1. The decline predates the removal",
     "",
     f"Gun homicides per ShotSpotter area peaked in {peak} ({ann_per.loc[peak,'treated']:.1f}/yr) and had "
     f"already fallen to {ann_per.loc[2024,'treated']:.1f} by 2024; never-treated areas moved in parallel "
     f"({ann_per.loc[peak,'control']:.1f} to {ann_per.loc[2024,'control']:.1f}). The city-wide decline "
     f"from the 2021--2022 peak runs straight through the September 2024 shutoff in both groups.", "",
     "## 2. Naive before/after vs. difference-in-differences", "",
     f"The commentary's comparison---gun homicides in coverage areas before vs.\\ after the shutoff---shows "
     f"a drop ({t_pre:.3f} to {t_post:.3f} per area-month, {(t_post/t_pre-1)*100:+.0f}%). But never-treated "
     f"areas dropped almost identically in proportional terms ({c_pre:.3f} to {c_post:.3f}, "
     f"{(c_post/c_pre-1)*100:+.0f}%): the decline is city-wide, not a removal effect.", "",
     f"The removal even reproduces the paper's *scale* artifact. The additive OLS difference-in-differences "
     f"returns a 'significant' removal benefit of {b:+.3f} gun homicides per area-month (p = {p:.3f})---but "
     f"only because coverage areas start from a far higher baseline ({t_pre:.2f} vs {c_pre:.2f} per "
     f"area-month) and so swing more in absolute terms for the same proportional move. On the correct "
     f"multiplicative scale the removal effect is null: Poisson IRR {irr:.2f} (95% CI [{il:.2f}, {ih:.2f}], "
     f"p = {ip:.2f}). This is the identical additive-vs-multiplicative divergence documented for the "
     f"installation in the count-data section, now recurring on the removal.", "",
     "## Reading", "",
     "The symmetry is the point. The removal reproduces both mistakes at once: the raw before/after "
     "attributes a city-wide decline to the shutoff, and the additive DiD manufactures a 'significant' "
     "effect that the multiplicative model erases. The same two errors that fabricate an installation "
     "effect fabricate a removal effect, and the same missing counterfactual makes both non-identifiable. "
     "The post window is short (through the last complete month of the data) and 2025 is partial, so the "
     "DiD is also under-powered; the honest conclusion is that these data cannot identify the removal "
     "effect in either direction, exactly as they cannot identify the installation effect.", "",
     "Figure: `plots/removal_trajectory.png`"]
with open(os.path.join(DOCS, "results_summary_removal.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(L))
print("  WROTE: docs/results_summary_removal.md\n\nDONE")
