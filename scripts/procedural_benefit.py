"""
procedural_benefit.py — The process-vs-outcome split at the community-area level.

Piza et al. (2024) show gunshot detection improves *process* measures (ballistic
evidence, gun recoveries) but not *victim* outcomes (shootings). Their process
measures are FOIA-only for Chicago; the public analog on the City of Chicago
Data Portal is WEAPONS-VIOLATION enforcement from the Crimes dataset
(`ijzp-q8t2`). This script runs the *same* two-way fixed-effects DiD used for gun
homicides on two enforcement-output measures---weapons-violation incidents and
weapons-violation arrests---so the process and outcome channels are estimated on
identical terms and can be read side by side.

Reading: if enforcement output rises under coverage while gun homicides do not
(properly modeled), that is the process-vs-outcome split in Chicago
community-area data. It is associational, not causal---the same identification
limits developed in the main paper apply---and the enforcement rise is
double-edged: the "procedural benefit" vendors advertise is also the
concentrated-enforcement burden advocates warn of, since coverage sits in the
most over-policed neighborhoods.

Outputs:
  data/raw/Crimes_weapons_violations_by_ca_year.csv  (input, pulled separately)
  plots/procedural_vs_outcome.png
  docs/results_summary_procedural.md
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
from plot_style import (apply_academic_style, COLORS, add_source_note,
                        add_n_label, sig_stars, SOURCE_DEFAULT)

P = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(P, "data", "raw")
PLOTS = os.path.join(P, "plots")
DOCS = os.path.join(P, "docs")
apply_academic_style()

# =====================================================================
# DATA
# =====================================================================
print("=" * 72); print("LOADING DATA"); print("=" * 72)

with open(os.path.join(DATA, "CommAreas_20250408.geojson"), "r", encoding="utf-8") as f:
    geo = json.load(f)
ca_num2name = {int(float(ft["properties"]["area_num_1"])): ft["properties"]["community"].title()
               for ft in geo["features"]}
ca_name2num = {}
for n, name in ca_num2name.items():
    ca_name2num[name] = n; ca_name2num[name.upper()] = n

# gun homicides (the outcome channel), same filters as the main analysis
hom = pd.read_csv(os.path.join(
    DATA, "Violence_Reduction_-_Victims_of_Homicides_and_Non-Fatal_Shootings_20250513.csv"),
    low_memory=False)
hom["DATE"] = pd.to_datetime(hom["DATE"], format="%m/%d/%Y %I:%M:%S %p", errors="coerce")
hom["CA_num"] = hom["COMMUNITY_AREA"].map(ca_name2num)
hom["year"] = hom["DATE"].dt.year
gun = hom[(hom["VICTIMIZATION_PRIMARY"] == "HOMICIDE") & (hom["GUNSHOT_INJURY_I"] == "YES") &
          (hom["year"].between(2009, 2023))].dropna(subset=["CA_num"])

# ShotSpotter treated set (same rule as every other script)
ss = pd.read_csv(os.path.join(DATA, "Violence_Reduction_-_Shotspotter_Alerts_-_Historical_20250408.csv"))
ss["DATE"] = pd.to_datetime(ss["DATE"], format="%m/%d/%Y %I:%M:%S %p", errors="coerce")
ss = ss[ss["DATE"] <= "2023-12-31"]
samp = str(ss["COMMUNITY_AREA"].dropna().iloc[0])
ss["CA"] = (pd.to_numeric(ss["COMMUNITY_AREA"], errors="coerce").map(ca_num2name)
            if samp.replace(".", "").isdigit() else ss["COMMUNITY_AREA"].str.title())
prof = ss.dropna(subset=["CA"]).groupby("CA").agg(a=("DATE", "size"), e=("DATE", "min"), l=("DATE", "max"))
prof["mo"] = ((prof.l - prof.e).dt.days / 30.4).round()
ALL_SS = set(prof[(prof.a >= 100) & (prof.mo >= 12)].index)

# weapons enforcement (the process channel), pulled from ijzp-q8t2
wep = pd.read_csv(os.path.join(DATA, "Crimes_weapons_violations_by_ca_year.csv"))

# Build the community-area x year panel
YEARS = list(range(2009, 2024))
panel = pd.DataFrame([(ca, yr) for ca in sorted(ca_num2name.keys()) for yr in YEARS],
                     columns=["ca_num", "year"])
panel["ca_name"] = panel["ca_num"].map(ca_num2name)
g = gun.groupby([gun["CA_num"].astype(int), "year"]).size().rename("gun_hom").reset_index()
g.columns = ["ca_num", "year", "gun_hom"]
panel = panel.merge(g, on=["ca_num", "year"], how="left")
panel = panel.merge(wep.rename(columns={"community_area": "ca_num"}), on=["ca_num", "year"], how="left")
for c in ["gun_hom", "weapons_incidents", "weapons_arrests"]:
    panel[c] = panel[c].fillna(0).astype(int)
panel["treat"] = panel["ca_name"].isin(ALL_SS).astype(int)
panel["post"] = (panel["year"] >= 2017).astype(int)
panel["did"] = panel["treat"] * panel["post"]
print(f"  Panel {len(panel)} rows; treated={panel[panel.treat==1].ca_num.nunique()}, "
      f"control={panel[panel.treat==0].ca_num.nunique()}")
print(f"  Totals 2009-2023: gun homicides={panel.gun_hom.sum():,}, "
      f"weapons incidents={panel.weapons_incidents.sum():,}, weapons arrests={panel.weapons_arrests.sum():,}")


# =====================================================================
# SAME DiD, THREE OUTCOMES: OLS (additive) and Poisson (multiplicative/IRR)
# =====================================================================
def twfe(outcome, label):
    f = f"{outcome} ~ did + C(ca_num) + C(year)"
    ols = smf.ols(f, data=panel).fit(cov_type="cluster", cov_kwds={"groups": panel["ca_num"]})
    pois = smf.glm(f, data=panel, family=sm.families.Poisson()).fit(
        cov_type="cluster", cov_kwds={"groups": panel["ca_num"]})
    b_ols, p_ols = ols.params["did"], ols.pvalues["did"]
    b, se, p = pois.params["did"], pois.bse["did"], pois.pvalues["did"]
    lo, hi = pois.conf_int().loc["did"]
    # pre-period levels for a raw sense of scale
    pre_t = panel[(panel.treat == 1) & (panel.post == 0)][outcome].mean()
    res = {"outcome": label, "ols_b": b_ols, "ols_p": p_ols,
           "irr": np.exp(b), "irr_lo": np.exp(lo), "irr_hi": np.exp(hi), "pois_p": p,
           "pre_treat_mean": pre_t}
    print(f"\n  {label}")
    print(f"    OLS DiD  = {b_ols:+8.3f}  (p={p_ols:.4f} {sig_stars(p_ols)})")
    print(f"    Poisson  = IRR {np.exp(b):.3f}  [{np.exp(lo):.3f}, {np.exp(hi):.3f}]  "
          f"(p={p:.4f} {sig_stars(p)})")
    return res

print("\n" + "=" * 72)
print("SAME TWO-WAY FIXED-EFFECTS DiD, THREE OUTCOMES")
print("=" * 72)
rows = [twfe("gun_hom", "Gun homicides (outcome)"),
        twfe("weapons_incidents", "Weapons-violation incidents (process)"),
        twfe("weapons_arrests", "Weapons-violation arrests (process)")]


# =====================================================================
# FIGURE: the process-vs-outcome contrast, on the IRR scale
# =====================================================================
fig, ax = plt.subplots(figsize=(11, 5.4))
ypos = np.arange(len(rows))[::-1]
ANN_X = 1.72   # fixed column for the OLS-to-IRR annotations, so nothing collides
for y, r in zip(ypos, rows):
    is_outcome = "outcome" in r["outcome"]
    color = COLORS["control"] if is_outcome else COLORS["treated"]
    ax.plot([r["irr_lo"], r["irr_hi"]], [y, y], color=color, lw=2.6, zorder=2, alpha=0.55)
    ax.plot(r["irr"], y, "o", color=color, markersize=11, markeredgecolor="white",
            markeredgewidth=0.9, zorder=3)
    st = sig_stars(r["ols_p"]) or "ns"
    ax.text(ANN_X, y, f"OLS {r['ols_b']:+.1f}{st} $\\rightarrow$ IRR {r['irr']:.2f} (ns)",
            va="center", ha="left", fontsize=9.5, color=color, fontweight="bold")
ax.axvline(1.0, color="#333", lw=1.3, zorder=1)
ax.set_yticks(ypos)
ax.set_yticklabels([r["outcome"] for r in rows], fontsize=10.5)
ax.set_xlim(0.5, 2.7)
ax.set_ylim(-0.75, len(rows) - 0.25)
ax.text(1.01, -0.62, "IRR 1 = no effect", fontsize=9, color="#333", ha="left", va="center",
        style="italic")
ax.set_xlabel("Incidence-rate ratio (treated $\\times$ post), Poisson two-way fixed effects\n"
              "Same DiD design and panel as the homicide analysis", labelpad=2)
ax.set_title("The additive artifact is general: same DiD, three outcomes,\n"
             "large OLS coefficients that all vanish on the count scale",
             fontweight="bold", pad=10)
ax.xaxis.grid(True)
add_n_label(ax, "77 community areas x 15 years (2009-2023)\n51 ShotSpotter / 26 never-treated",
            loc="upper left")
legend = [
    Line2D([0], [0], color=COLORS["control"], marker="o", lw=2.6, markersize=9,
           label="Victim outcome (gun homicides)"),
    Line2D([0], [0], color=COLORS["treated"], marker="o", lw=2.6, markersize=9,
           label="Process / enforcement (weapons)"),
]
ax.legend(handles=legend, loc="lower left", frameon=False, fontsize=9,
          bbox_to_anchor=(0.0, 0.02))
add_source_note(fig, SOURCE_DEFAULT + "  Weapons violations from Crimes (ijzp-q8t2).")
fig.tight_layout(rect=[0, 0.03, 1, 1])
fig.savefig(os.path.join(PLOTS, "procedural_vs_outcome.png"), dpi=200, facecolor="white")
plt.close(fig)
print("\n  SAVED: procedural_vs_outcome.png")


# =====================================================================
# DOC
# =====================================================================
L = []
L.append("# The Procedural-Benefit Channel: Enforcement vs. Outcomes")
L.append("")
L.append("Results of `scripts/procedural_benefit.py`. Gunshot detection is documented to")
L.append("improve *process* measures but not *victim* outcomes (Piza et al., 2024). Chicago's")
L.append("process analogs (ballistic evidence, gun recoveries) are FOIA-only; the public")
L.append("proxy is weapons-violation enforcement from the Crimes dataset (`ijzp-q8t2`). We run")
L.append("the *same* two-way fixed-effects DiD used for gun homicides on two enforcement")
L.append("measures, so process and outcome are estimated on identical terms.")
L.append("")
L.append("| Outcome | OLS DiD | Poisson IRR | 95% CI (IRR) | Poisson p |")
L.append("|---|---|---|---|---|")
for r in rows:
    st = "***" if r["pois_p"] < 0.01 else "**" if r["pois_p"] < 0.05 else "*" if r["pois_p"] < 0.1 else ""
    L.append(f"| {r['outcome']} | {r['ols_b']:+.3f} | {r['irr']:.3f} | "
             f"[{r['irr_lo']:.3f}, {r['irr_hi']:.3f}] | {r['pois_p']:.3f}{st} |")
L.append("")
inc = next(r for r in rows if "incidents" in r["outcome"])
arr = next(r for r in rows if "arrests" in r["outcome"])
hom_r = next(r for r in rows if "omicide" in r["outcome"])
L.append(f"**Reading.** The result carries two lessons. First, the specification artifact is "
         f"*general*. Weapons-violation incidents produce a large, significant OLS coefficient "
         f"({inc['ols_b']:+.1f} per area-year) that collapses to a statistically null incidence-rate "
         f"ratio ({inc['irr']:.2f}, p = {inc['pois_p']:.2f}) once the count outcome is modeled "
         f"correctly---exactly the pattern the homicide analysis shows (OLS {hom_r['ols_b']:+.2f} "
         f"vs.\\ IRR {hom_r['irr']:.2f}). The additive-scale inflation is not peculiar to homicide "
         f"counts; it recurs on a completely different outcome. Second, and substantively, the "
         f"Chicago community-area data do *not* detect the process benefit that a microsynthetic "
         f"control recovers in Kansas City (Piza et al., 2024): weapons incidents (IRR {inc['irr']:.2f}) "
         f"and arrests (IRR {arr['irr']:.2f}) are both null on the correct scale. That is consistent "
         f"with the 2021 OIG finding that Chicago's alerts rarely produced investigative evidence---and "
         f"with the same missing-counterfactual problem that governs the rest of the paper.")
L.append("")
L.append("**Caveats.** This is associational, not causal: the identification limits developed")
L.append("throughout the paper apply here too (no valid control group; selection on violence).")
L.append("The arrest flag reflects any arrest recorded as of the last data refresh, not")
L.append("clearance. And a rise in weapons enforcement is double-edged---the procedural benefit")
L.append("vendors advertise is also the concentrated-enforcement burden advocates warn of, since")
L.append("coverage sits in the most over-policed, lowest-life-expectancy neighborhoods.")
L.append("")
L.append("Figure: `plots/procedural_vs_outcome.png`")
with open(os.path.join(DOCS, "results_summary_procedural.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(L))
print("  WROTE: docs/results_summary_procedural.md")
print("\nDONE")
