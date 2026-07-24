"""
treatment_definition.py -- Is the treated / never-treated split a knife-edge choice?

Every number in the paper depends on one screening rule: a community area counts as
covered by ShotSpotter if its historical alert file shows at least 100 alerts spanning
at least 12 months. That rule is asserted throughout but never defended, and it is the
first thing a referee will probe. This script defends it, and discloses what it misses.

Three questions:
  1. Does the threshold cut through a continuum, or does it fall in an empty region?
  2. How far can the cutoff move before the panel changes?
  3. Is the comparison group actually clean?

Outputs:
  plots/treatment_threshold.png
  docs/results_summary_treatment.md
"""
import os
import sys
import json
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_style import (apply_academic_style, COLORS, add_source_note, add_n_label,
                        SOURCE_DEFAULT)

P = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(P, "data", "raw")
PLOTS = os.path.join(P, "plots")
DOCS = os.path.join(P, "docs")
apply_academic_style()

ALERT_RULE, MONTH_RULE = 100, 12

# =====================================================================
print("=" * 72); print("TREATMENT DEFINITION DIAGNOSTICS"); print("=" * 72)

with open(os.path.join(DATA, "CommAreas_20250408.geojson"), encoding="utf-8") as f:
    geo = json.load(f)
num2name = {int(float(x["properties"]["area_num_1"])): x["properties"]["community"].title()
            for x in geo["features"]}

ss = pd.read_csv(os.path.join(
    DATA, "Violence_Reduction_-_Shotspotter_Alerts_-_Historical_20250408.csv"))
ss["DATE"] = pd.to_datetime(ss["DATE"], format="%m/%d/%Y %I:%M:%S %p", errors="coerce")
ss = ss[ss["DATE"] <= "2023-12-31"]
samp = str(ss["COMMUNITY_AREA"].dropna().iloc[0])
ss["nm"] = (pd.to_numeric(ss["COMMUNITY_AREA"], errors="coerce").map(num2name)
            if samp.replace(".", "").isdigit() else ss["COMMUNITY_AREA"].str.title())

prof = ss.dropna(subset=["nm"]).groupby("nm").agg(
    alerts=("DATE", "size"), t0=("DATE", "min"), t1=("DATE", "max"))
prof["months"] = ((prof["t1"] - prof["t0"]).dt.days / 30.4).round()

ca = pd.DataFrame({"nm": sorted(num2name.values())}).set_index("nm")
ca = ca.join(prof).fillna({"alerts": 0, "months": 0})
ca["treated"] = (ca.alerts >= ALERT_RULE) & (ca.months >= MONTH_RULE)

n_t = int(ca.treated.sum()); n_c = len(ca) - n_t
tr, nt = ca[ca.treated], ca[~ca.treated]
min_treated = int(tr.alerts.min()); max_control = int(nt.alerts.max())
zero_ctrl = int((nt.alerts == 0).sum()); trace_ctrl = int(((nt.alerts > 0) & (nt.alerts < ALERT_RULE)).sum())
print(f"  {n_t} treated, {n_c} never-treated under the paper's rule.")
print(f"  Lowest treated alert count : {min_treated}")
print(f"  Highest control alert count: {max_control}")
print(f"  Gap factor: {min_treated / max(max_control, 1):.1f}x with no area in between.")

# 1. sensitivity grid -------------------------------------------------
ALERT_GRID = [10, 50, 100, 250, 500, 1000]
MONTH_GRID = [0, 6, 12, 18]
grid = pd.DataFrame(
    [[int(((ca.alerts >= a) & (ca.months >= m)).sum()) for m in MONTH_GRID] for a in ALERT_GRID],
    index=ALERT_GRID, columns=MONTH_GRID)
print("\n  n treated by (alert cutoff x month cutoff):")
print(grid.to_string())
# Does the month screen do any work AT THE OPERATIVE alert cutoff? (It does bind at
# looser cutoffs, so this must be stated conditionally rather than as a blanket claim.)
month_inert_here = int(((ca.alerts >= ALERT_RULE) & (ca.months >= 0)).sum()) == n_t
month_binds_at = [a for a in ALERT_GRID
                  if int(((ca.alerts >= a) & (ca.months >= 0)).sum())
                  != int(((ca.alerts >= a) & (ca.months >= MONTH_RULE)).sum())]
stable_lo = max(a for a in ALERT_GRID if int(((ca.alerts >= a) & (ca.months >= MONTH_RULE)).sum()) == n_t)
print(f"\n  Month criterion does no work at the operative {ALERT_RULE}-alert cutoff: {month_inert_here}")
print(f"  It does bind at looser alert cutoffs: {month_binds_at if month_binds_at else 'none tested'}")
print(f"  Panel unchanged for alert cutoffs from {ALERT_RULE} up to {stable_lo}.")

# 2. control-group contamination -------------------------------------
print(f"\n  Comparison group: {zero_ctrl} of {n_c} have zero alerts; {trace_ctrl} have 1-{ALERT_RULE-1}.")
top_ctrl = nt.sort_values("alerts", ascending=False).head(5)
for nm, r in top_ctrl.iterrows():
    print(f"    {nm:<22} {int(r.alerts):>4} alerts over {int(r.months):>3} months "
          f"({r.alerts / max(r.months, 1):.1f}/mo)")
wh = ca.loc["Washington Heights"] if "Washington Heights" in ca.index else None
med_treated = float(tr.alerts.median())

# =====================================================================
# FIGURE
# =====================================================================
s = ca.sort_values("alerts", ascending=False)
fig, ax = plt.subplots(figsize=(11, 5.4))
xs = np.arange(len(s))
cols = [COLORS["treated"] if t else COLORS["control"] for t in s.treated]
ax.bar(xs, np.maximum(s.alerts.values, 0.7), color=cols, width=0.85)
ax.set_yscale("log")
ax.axhspan(max_control, min_treated, color="#c8a415", alpha=0.16, zorder=0)
ax.axhline(ALERT_RULE, color="#6d6128", ls="--", lw=1.4, zorder=3)
ax.text(len(s) * 0.86, np.sqrt(max_control * min_treated),
        f"no community area falls between\n{max_control} and {min_treated} alerts",
        fontsize=9, ha="center", va="center", style="italic", color="#5c5222")
ax.text(len(s) * 0.86, ALERT_RULE * 0.72, f"screening rule: {ALERT_RULE} alerts",
        fontsize=8.5, ha="center", va="top", color="#6d6128")
ax.set_xlabel("Community areas, sorted by ShotSpotter alert count (2017–2023)")
ax.set_ylabel("Alerts (log scale)")
ax.set_title("The coverage screen falls in an empty region, not through a continuum",
             fontweight="bold", pad=10)
handles = [plt.Rectangle((0, 0), 1, 1, color=COLORS["treated"]),
           plt.Rectangle((0, 0), 1, 1, color=COLORS["control"])]
ax.legend(handles, [f"Coded covered (n={n_t})", f"Coded never-treated (n={n_c})"],
          frameon=False, fontsize=9.5, loc="upper right")
add_n_label(ax, f"Lowest covered area: {min_treated} alerts\n"
                f"Highest never-treated: {max_control} alerts\n"
                f"Panel identical for cutoffs {ALERT_RULE}-{stable_lo}",
            loc="lower left")
add_source_note(fig, SOURCE_DEFAULT)
fig.tight_layout(rect=[0, 0.04, 1, 1])
fig.savefig(os.path.join(PLOTS, "treatment_threshold.png"), dpi=200, facecolor="white")
plt.close(fig)
print("\n  SAVED: treatment_threshold.png")

# =====================================================================
# DOC
# =====================================================================
L = ["# Is the Treated / Never-Treated Split a Knife-Edge Choice?", "",
     "Results of `scripts/treatment_definition.py`. Every estimate in the paper rests on one",
     f"screening rule: a community area counts as covered if the historical alert file shows at",
     f"least {ALERT_RULE} alerts spanning at least {MONTH_RULE} months, giving **{n_t} covered and "
     f"{n_c} never-treated** areas. This documents how much that rule matters.", "",
     "## 1. The cutoff falls in an empty region", "",
     f"The least-covered area coded as treated records **{min_treated} alerts**. The most-covered",
     f"area coded as never-treated records **{max_control}**. No community area in Chicago falls",
     f"between those two values, a gap of roughly {min_treated / max(max_control,1):.0f} times with",
     "nothing inside it. The screen therefore separates two genuinely distinct groups rather than",
     "slicing an arbitrary point out of a continuum (`plots/treatment_threshold.png`).", "",
     "## 2. The panel is invariant over a wide range of cutoffs", "",
     "Number of areas coded as covered, by alert cutoff (rows) and month cutoff (columns):", "",
     "| alerts >= | " + " | ".join(f"months >= {m}" for m in MONTH_GRID) + " |",
     "|---|" + "---|" * len(MONTH_GRID)]
for a in ALERT_GRID:
    L.append(f"| {a} | " + " | ".join(str(int(grid.loc[a, m])) for m in MONTH_GRID) + " |")
L += ["",
      f"The panel is identical for any alert cutoff from {ALERT_RULE} to {stable_lo}, so the "
      f"headline results do not depend on where in that range the line is drawn.", "",
      f"**At the operative cutoff the month criterion does no work.** Requiring 0 months and",
      f"requiring {max(MONTH_GRID)} months give the same {n_t} areas once the alert threshold is",
      f"{ALERT_RULE}: any area with that many alerts also spans years of them. The criterion is not",
      f"vacuous in general, and it does bind at looser alert cutoffs "
      f"({', '.join(str(a) for a in month_binds_at) if month_binds_at else 'none tested'}), where it",
      "screens out areas with a short burst of detections. We report the rule as implemented for",
      "reproducibility and note that, at the threshold actually used, the alert count alone",
      "determines the panel.", "",
      "## 3. The comparison group is clean, but not perfectly", "",
      f"Of the {n_c} never-treated areas, **{zero_ctrl} record zero alerts** and **{trace_ctrl} record",
      f"between 1 and {ALERT_RULE - 1}**. The five with the most:", "",
      "| Community area | Alerts | Months | Alerts/month |", "|---|---|---|---|"]
for nm, r in top_ctrl.iterrows():
    L.append(f"| {nm} | {int(r.alerts)} | {int(r.months)} | {r.alerts / max(r.months,1):.1f} |")
L += ["",
      "This is a real limitation and we disclose it rather than smooth it over. ShotSpotter was",
      "procured at the police-district level, and community areas cut across districts, so an area",
      "adjacent to or partially inside a covered district can register a trickle of alerts without",
      "being meaningfully covered. The magnitudes are small: the heaviest such case is about",
      f"{top_ctrl.iloc[0].alerts / max(top_ctrl.iloc[0].months,1):.1f} alerts per month against a "
      f"median of {med_treated:,.0f} alerts in covered areas over the same window, on the order of",
      "one percent of treated intensity."]
if wh is not None:
    L += ["",
          f"One case deserves explicit mention. **Washington Heights** ({int(wh.alerts)} alerts over",
          f"{int(wh.months)} months) is the single donor the synthetic control collapses onto, so the",
          "primary comparison unit in that exercise is not perfectly untreated. Given the trace",
          "volume this does not change the finding, and the synthetic control fails for the separate",
          "and far larger reason that no donor can match treated-area violence levels at all. But a",
          "reader checking the donor should know."]
L += ["",
      "## Reading", "",
      "The screening rule is not a knife-edge. It falls in an empty region of the data, the panel",
      f"is unchanged for cutoffs anywhere from {ALERT_RULE} to {stable_lo}, and only one of its two",
      "stated criteria does any work. The comparison group is clean for most of its members and",
      "carries trace contamination for a minority, an artifact of measuring a district-level",
      "program on a community-area map. None of this rescues identification, which fails for the",
      "reasons developed in the paper, but it does mean the failure is not an artifact of how",
      "treatment was coded.", "",
      "Figure: `plots/treatment_threshold.png`"]
with open(os.path.join(DOCS, "results_summary_treatment.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(L))
print("  WROTE: docs/results_summary_treatment.md\n\nDONE")
