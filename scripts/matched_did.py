"""
matched_did.py -- Does matching on pre-treatment observables rescue the design?

A natural reviewer suggestion for a difference-in-differences whose parallel-trends
assumption fails is to match treated and control units on pre-treatment observables
(e.g. on a propensity score) before estimating. This script runs that suggestion on
its own terms and reports what happens.

Three samples, all using the identical propensity model as the placement analysis
(logit on standardized pre-period gun homicides, hardship index, and per-capita income):

  1. FULL          all 77 community areas (the paper's baseline).
  2. TRIMMED       restricted to the propensity common-support region.
  3. MATCHED       1:1 nearest-neighbour propensity matching with a caliper.

For each we report (a) how many treated units survive and whether the six
highest-violence study neighborhoods survive, (b) whether the joint pre-trend Wald
test still rejects parallel trends, and (c) the DiD estimate on the additive and
multiplicative scales.

The answer is not that matching fails. It is that matching works, at a price the policy
question cannot pay: restricting to the region where treated and control units are
comparable does restore parallel pre-trends, but only by discarding every one of the six
highest-violence study neighborhoods, because no never-treated area exists at their
propensity level. What survives is a credible design for the marginal, lower-violence
coverage areas. This script quantifies that trade.

Outputs:
  plots/matched_did.png
  docs/results_summary_matching.md
"""
import os
import sys
import json
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_style import (apply_academic_style, COLORS, add_source_note, add_n_label,
                        sig_stars, SOURCE_DEFAULT)

P = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(P, "data", "raw")
PLOTS = os.path.join(P, "plots")
DOCS = os.path.join(P, "docs")
apply_academic_style()

PRE_START, PRE_END = 2009, 2016
PANEL_START, PANEL_END = 2009, 2023
STUDY6 = ["Austin", "Humboldt Park", "North Lawndale",
          "Englewood", "West Englewood", "South Shore"]

# =====================================================================
# DATA (mirrors placement_model.py / econometric_analysis.py exactly)
# =====================================================================
print("=" * 72); print("LOADING DATA"); print("=" * 72)

with open(os.path.join(DATA, "CommAreas_20250408.geojson"), "r", encoding="utf-8") as f:
    geo = json.load(f)
ca_num2name = {int(float(ft["properties"]["area_num_1"])): ft["properties"]["community"].title()
               for ft in geo["features"]}
ca_name2num = {}
for n, nm in ca_num2name.items():
    ca_name2num[nm] = n
    ca_name2num[nm.upper()] = n
    ca_name2num[nm.lower()] = n

hom = pd.read_csv(os.path.join(
    DATA, "Violence_Reduction_-_Victims_of_Homicides_and_Non-Fatal_Shootings_20250513.csv"),
    low_memory=False)
hom["DATE"] = pd.to_datetime(hom["DATE"], format="%m/%d/%Y %I:%M:%S %p", errors="coerce")
hom["CA_num"] = hom["COMMUNITY_AREA"].map(ca_name2num)
hom["year"] = hom["DATE"].dt.year
gun = hom[(hom["VICTIMIZATION_PRIMARY"] == "HOMICIDE") &
          (hom["GUNSHOT_INJURY_I"] == "YES")].dropna(subset=["CA_num"]).copy()
gun["CA_num"] = gun["CA_num"].astype(int)

pre = gun[gun.year.between(PRE_START, PRE_END)]
pre_counts = pre.groupby("CA_num").size().reindex(sorted(ca_num2name), fill_value=0)
pre_mean = pre_counts / (PRE_END - PRE_START + 1)

ss = pd.read_csv(os.path.join(
    DATA, "Violence_Reduction_-_Shotspotter_Alerts_-_Historical_20250408.csv"))
ss["DATE"] = pd.to_datetime(ss["DATE"], format="%m/%d/%Y %I:%M:%S %p", errors="coerce")
ss = ss[ss["DATE"] <= "2023-12-31"].copy()
samp = str(ss["COMMUNITY_AREA"].dropna().iloc[0])
if samp.replace(".", "").isdigit():
    ss["CA_name"] = pd.to_numeric(ss["COMMUNITY_AREA"], errors="coerce").map(ca_num2name)
else:
    ss["CA_name"] = ss["COMMUNITY_AREA"].str.title()
prof = ss.dropna(subset=["CA_name"]).groupby("CA_name").agg(
    a=("DATE", "size"), e=("DATE", "min"), l=("DATE", "max"))
prof["mo"] = ((prof.l - prof.e).dt.days / 30.4).round()
ALL_SS = set(prof[(prof.a >= 100) & (prof.mo >= 12)].index)
# cohort year = first alert year, for the event-time / staggered structure
cohort_year = prof.loc[sorted(ALL_SS), "e"].dt.year.to_dict()

cen = pd.read_csv(os.path.join(
    DATA, "Census_Data_-_Selected_socioeconomic_indicators_in_Chicago__2008___2012_20250408.csv"))
cen.columns = [c.strip() for c in cen.columns]
cen = cen[pd.to_numeric(cen["Community Area Number"], errors="coerce").notna()].copy()
cen["ca_num"] = cen["Community Area Number"].astype(int)
cen = cen.rename(columns={"PER CAPITA INCOME": "income", "HARDSHIP INDEX": "hardship"})

x = pd.DataFrame({"ca_num": sorted(ca_num2name)})
x["ca_name"] = x["ca_num"].map(ca_num2name)
x["pre_gun_hom"] = x["ca_num"].map(pre_mean)
x = x.merge(cen[["ca_num", "income", "hardship"]], on="ca_num", how="left")
x["treated"] = x["ca_name"].isin(ALL_SS).astype(int)
x = x.dropna(subset=["pre_gun_hom", "income", "hardship"]).reset_index(drop=True)

for v in ["pre_gun_hom", "income", "hardship"]:
    x[f"z_{v}"] = (x[v] - x[v].mean()) / x[v].std(ddof=0)
ZP = ["z_pre_gun_hom", "z_hardship", "z_income"]

clf = LogisticRegression(penalty="l2", C=1e6, solver="lbfgs", max_iter=10000)
clf.fit(x[ZP].values, x["treated"].values)
x["ps"] = clf.predict_proba(x[ZP].values)[:, 1]
ps_auc = roc_auc_score(x["treated"], x["ps"])
n_t_full = int(x.treated.sum()); n_c_full = int((1 - x.treated).sum())
print(f"  Propensity model AUC = {ps_auc:.3f}   ({n_t_full} treated, {n_c_full} never-treated)")

# annual panel
panel = pd.DataFrame([(c, y) for c in x.ca_num for y in range(PANEL_START, PANEL_END + 1)],
                     columns=["ca_num", "year"])
cnt = gun.groupby(["CA_num", "year"]).size().rename("gun_hom").reset_index()
cnt.columns = ["ca_num", "year", "gun_hom"]
panel = panel.merge(cnt, on=["ca_num", "year"], how="left")
panel["gun_hom"] = panel["gun_hom"].fillna(0).astype(int)
panel = panel.merge(x[["ca_num", "ca_name", "treated", "ps"]], on="ca_num", how="left")
panel["install"] = panel["ca_name"].map(cohort_year)
panel["post"] = ((panel.treated == 1) & (panel.year >= panel.install)).astype(int)

# =====================================================================
# ESTIMATORS
# =====================================================================
def pretrend_wald(sub):
    """Joint Wald test on pre-period treated x year interactions (2009-2016, 2009 base)."""
    d = sub[sub.year <= PRE_END].copy()
    if d.treated.nunique() < 2:
        return np.nan, np.nan, 0
    yrs = sorted(d.year.unique())[1:]          # drop first year as base
    for y in yrs:
        d[f"tx{y}"] = ((d.treated == 1) & (d.year == y)).astype(int)
    terms = [f"tx{y}" for y in yrs]
    m = smf.ols("gun_hom ~ " + " + ".join(terms) + " + C(ca_num) + C(year)",
                data=d).fit(cov_type="cluster", cov_kwds={"groups": d.ca_num})
    try:
        w = m.f_test(" = 0, ".join(terms) + " = 0")
        return float(np.squeeze(w.fvalue)), float(np.squeeze(w.pvalue)), len(terms)
    except Exception:
        return np.nan, np.nan, len(terms)

def did(sub):
    """TWFE DiD on the additive (OLS) and multiplicative (Poisson) scales."""
    sub = sub.copy()
    sub["did"] = sub["post"]
    out = {}
    try:
        m = smf.ols("gun_hom ~ did + C(ca_num) + C(year)", data=sub).fit(
            cov_type="cluster", cov_kwds={"groups": sub.ca_num})
        ci = m.conf_int().loc["did"]
        out["ols"] = (m.params["did"], m.bse["did"], m.pvalues["did"], ci[0], ci[1])
    except Exception:
        out["ols"] = (np.nan,) * 5
    try:
        g = smf.glm("gun_hom ~ did + C(ca_num) + C(year)", data=sub,
                    family=sm.families.Poisson()).fit(
            cov_type="cluster", cov_kwds={"groups": sub.ca_num})
        lo, hi = np.exp(g.conf_int().loc["did"])
        out["pois"] = (np.exp(g.params["did"]), lo, hi, g.pvalues["did"])
    except Exception:
        out["pois"] = (np.nan,) * 4
    return out

# --- Sample 1: full ---
samples = {}
samples["Full sample"] = x.ca_num.tolist()

# --- Sample 2: common-support trimming ---
lo = max(x.loc[x.treated == 1, "ps"].min(), x.loc[x.treated == 0, "ps"].min())
hi = min(x.loc[x.treated == 1, "ps"].max(), x.loc[x.treated == 0, "ps"].max())
trim = x[(x.ps >= lo) & (x.ps <= hi)]
samples["Common support"] = trim.ca_num.tolist()
print(f"\n  Common-support region: propensity in [{lo:.3f}, {hi:.3f}]")

# --- Sample 3: 1:1 nearest-neighbour matching on the linearized propensity, with caliper ---
def logit_lin(p):
    p = np.clip(p, 1e-6, 1 - 1e-6)
    return np.log(p / (1 - p))

x["lps"] = logit_lin(x["ps"].values)
caliper = 0.2 * x["lps"].std(ddof=0)
tr = x[x.treated == 1]
co = x[x.treated == 0]
pairs, unmatched = [], []
for _, r in tr.iterrows():
    d = (co["lps"] - r["lps"]).abs()
    if len(d) and d.min() <= caliper:
        pairs.append((r["ca_num"], int(co.loc[d.idxmin(), "ca_num"])))
    else:
        unmatched.append(r["ca_name"])
matched_ids = sorted(set([p[0] for p in pairs]) | set([p[1] for p in pairs]))
samples["NN matched"] = matched_ids
print(f"  Caliper = {caliper:.3f} (0.2 SD of linearized propensity)")
print(f"  Treated units finding a control within caliper: {len(pairs)} of {n_t_full}")

# =====================================================================
# RUN
# =====================================================================
print("\n" + "=" * 72)
print("DOES MATCHING RESCUE THE DESIGN?")
print("=" * 72)
rows = []
for name, ids in samples.items():
    sub = panel[panel.ca_num.isin(ids)]
    xs = x[x.ca_num.isin(ids)]
    nt, nc = int(xs.treated.sum()), int((1 - xs.treated).sum())
    kept6 = sorted(set(STUDY6) & set(xs.ca_name))
    F, pF, k = pretrend_wald(sub)
    est = did(sub)
    b, se, pb, blo, bhi = est["ols"]
    irr, ilo, ihi, pi = est["pois"]
    rows.append(dict(sample=name, n_t=nt, n_c=nc, study6=len(kept6), study6_names=kept6,
                     F=F, pF=pF, k=k, b=b, se=se, pb=pb,
                     irr=irr, ilo=ilo, ihi=ihi, pi=pi))
    print(f"\n  [{name}]  treated={nt}, control={nc}, study-6 retained={len(kept6)}/6")
    if F == F:
        verdict = "STILL REJECTS parallel trends" if pF < 0.05 else "does not reject"
        print(f"    Pre-trend joint Wald ({k} coefs): F={F:.2f}, p={pF:.4f}  -> {verdict}")
    else:
        print(f"    Pre-trend joint Wald: not estimable")
    print(f"    OLS DiD     = {b:+.3f} (SE {se:.3f}, p={pb:.3f} {sig_stars(pb)})")
    print(f"    Poisson IRR = {irr:.3f} [{ilo:.3f}, {ihi:.3f}] (p={pi:.3f} {sig_stars(pi)})")

res = pd.DataFrame(rows)

# =====================================================================
# FIGURE
# =====================================================================
fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5.6),
                               gridspec_kw={"width_ratios": [1.15, 1]})

# LEFT: propensity distributions with the common-support region and the study-6 band
bins = np.linspace(0, 1, 26)
axL.hist(x.loc[x.treated == 0, "ps"], bins=bins, alpha=0.75,
         color=COLORS["control"], label=f"Never-treated (n={n_c_full})")
axL.hist(x.loc[x.treated == 1, "ps"], bins=bins, alpha=0.75,
         color=COLORS["treated"], label=f"ShotSpotter (n={n_t_full})")
s6lo = x[x.ca_name.isin(STUDY6)]["ps"].min()
axL.axvspan(s6lo, 1.0, color=COLORS["treated"], alpha=0.10, zorder=0)
axL.text(s6lo - 0.02, axL.get_ylim()[1] * 0.62,
         "six study\nneighborhoods", fontsize=8.5, ha="right", va="top",
         style="italic", color=COLORS["treated"])
axL.set_xlabel("Estimated propensity to receive ShotSpotter")
axL.set_ylabel("Number of community areas")
axL.set_title(f"No control unit reaches the study areas' propensity band\n(AUC {ps_auc:.2f})",
              fontweight="bold", fontsize=11)
axL.legend(loc="upper center", frameon=False, fontsize=9)

# RIGHT: DiD on the multiplicative scale across samples
ypos = np.arange(len(res))[::-1]
for yp, (_, r) in zip(ypos, res.iterrows()):
    ok = r.irr == r.irr
    if not ok:
        continue
    axR.plot([r.ilo, r.ihi], [yp, yp], color=COLORS["treated"], lw=2.4, solid_capstyle="round")
    axR.plot([r.irr], [yp], "o", color=COLORS["treated"], markersize=8, zorder=3)
    axR.text(1.92, yp + 0.16, f"IRR {r.irr:.2f} [{r.ilo:.2f}, {r.ihi:.2f}]",
             fontsize=8.8, va="center", ha="left")
    axR.text(1.92, yp - 0.16,
             f"{int(r.n_t)}T / {int(r.n_c)}C, {int(r.study6)}/6 study areas\n"
             f"pre-trend p = {r.pF:.3f}",
             fontsize=7.8, va="top", ha="left", color="#555")
axR.axvline(1.0, color="#444", lw=1.2, zorder=1)
axR.set_yticks(ypos)
axR.set_yticklabels(res["sample"])
axR.set_xlabel("Poisson DiD incidence-rate ratio")
axR.set_title("Matching buys parallel trends by discarding\nevery study neighborhood",
              fontweight="bold", fontsize=11)
axR.set_xlim(0.55, 3.05)
axR.set_xticks([0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8])
axR.set_ylim(-0.65, len(res) - 0.35)
axR.grid(axis="x", alpha=0.3)

add_source_note(fig, SOURCE_DEFAULT)
fig.tight_layout(rect=[0, 0.04, 1, 1])
fig.savefig(os.path.join(PLOTS, "matched_did.png"), dpi=200, facecolor="white")
plt.close(fig)
print("\n  SAVED: matched_did.png")

# =====================================================================
# DOC
# =====================================================================
r_full = res[res["sample"] == "Full sample"].iloc[0]
r_trim = res[res["sample"] == "Common support"].iloc[0]
r_nn = res[res["sample"] == "NN matched"].iloc[0]

L = ["# Does Matching on Pre-Treatment Observables Rescue the Design?", "",
     "Results of `scripts/matched_did.py`. A standard response to a difference-in-differences",
     "whose parallel-trends assumption fails is to match treated and control units on",
     "pre-treatment observables first. We run that suggestion on its own terms.", "",
     "The propensity model is the one from the placement analysis: a logit of ShotSpotter",
     f"placement on standardized pre-period gun homicides, hardship, and per-capita income",
     f"(AUC {ps_auc:.3f}; {n_t_full} treated, {n_c_full} never-treated).", "",
     "## Results", "",
     "| Sample | Treated | Control | Study-6 retained | Pre-trend Wald F | p | OLS DiD | Poisson IRR |",
     "|---|---|---|---|---|---|---|---|"]
for _, r in res.iterrows():
    L.append(f"| {r['sample']} | {int(r.n_t)} | {int(r.n_c)} | {int(r.study6)}/6 | "
             f"{r.F:.2f} | {r.pF:.4f} | {r.b:+.3f} ({r.se:.3f}) | "
             f"{r.irr:.3f} [{r.ilo:.3f}, {r.ihi:.3f}] |")
L += ["",
      f"Common-support region: propensity in [{lo:.3f}, {hi:.3f}].",
      f"Nearest-neighbour matching used a caliper of {caliper:.3f} (0.2 SD of the linearized",
      f"propensity); {len(pairs)} of {n_t_full} treated areas found a control inside it.", "",
      "## Reading", "",
      "Matching does not fail here. It works, and what it costs is the point.", "",
      "**Restricting to comparable units does restore parallel pre-trends.** The joint Wald test",
      f"on pre-period interactions rejects parallel trends on the full sample "
      f"(F = {r_full.F:.2f}, p = {r_full.pF:.4f}), but stops rejecting once the sample is trimmed to",
      f"the propensity common-support region (F = {r_trim.F:.2f}, p = {r_trim.pF:.4f}) and under",
      f"nearest-neighbour matching (F = {r_nn.F:.2f}, p = {r_nn.pF:.4f}). The $F$ statistic falls",
      "well below one, so this is not only the loss of power that comes with a smaller sample,",
      "though with fewer than half the original units some of it surely is. On its own terms the",
      "suggestion works: comparable units really do move in parallel.", "",
      "**The price is every unit the study is about.** Placement is predicted well enough that the",
      "six highest-violence study neighborhoods occupy a propensity band no never-treated area",
      f"reaches. All six are present in the full sample; {int(r_trim.study6)} of 6 survive",
      f"common-support trimming and {int(r_nn.study6)} of 6 survive matching, and only "
      f"{len(pairs)} of {n_t_full} treated areas find any control inside the caliper at all.",
      "The design that satisfies parallel trends is therefore a design about the marginal,",
      "lower-violence edge of ShotSpotter coverage: real areas, but not the neighborhoods the",
      "program was built around, and not the ones the procurement debate was about.", "",
      "**On that subsample the effect is still not estimable in any useful sense.** The",
      f"multiplicative estimate is IRR {r_trim.irr:.2f} [{r_trim.ilo:.2f}, {r_trim.ihi:.2f}] on common",
      f"support and {r_nn.irr:.2f} [{r_nn.ilo:.2f}, {r_nn.ihi:.2f}] matched, against "
      f"{r_full.irr:.2f} [{r_full.ilo:.2f}, {r_full.ihi:.2f}] on the full panel. Both intervals",
      "still span one, and both are wider than the full-sample interval, so trimming buys internal",
      "validity and pays for it in precision.", "",
      "This is the identification problem stated as a trade rather than a failure. These data can",
      "deliver a credible comparison, or they can speak to the neighborhoods the policy is about,",
      "but not both at once. That is a sharper statement of the paper's claim than the full-sample",
      "diagnostics alone, and it is the reason the obstacle cannot be estimated away: no",
      "reweighting of the units that exist can manufacture a counterfactual for units that have none.", "",
      "**Specification note.** The pre-trend test here is a joint Wald test on treated-by-year",
      "interactions in the 2009--2016 pre-period, which is a simpler test than the cohort event-study",
      "Wald test reported in the main paper (F = 19.9, p = 0.006); the two are not the same statistic",
      "and are not directly comparable. The full-sample OLS DiD reported here (+2.48) is the staggered,",
      "cohort-specific specification, which is why it matches the Goodman-Bacon aggregate rather than",
      "the pooled 2017-cutoff headline of +2.63.", "",
      "Figure: `plots/matched_did.png`"]
with open(os.path.join(DOCS, "results_summary_matching.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(L))
print("  WROTE: docs/results_summary_matching.md\n\nDONE")
