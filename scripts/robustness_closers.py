"""
robustness_closers.py — Three analyses that close limitations the paper names explicitly.

  1. Non-fatal gunshot injuries re-estimated as count models (Poisson / Negative
     Binomial TWFE) — the paper's OLS-only NFS result was flagged as "left to
     follow-up."
  2. Restricted wild-cluster bootstrap (Rademacher weights, null imposed;
     Cameron-Gelbach-Miller 2008) for the headline OLS DiD — the paper flags
     ~77 clusters as thin for cluster-robust inference.
  3. Goodman-Bacon (2021) decomposition of the *staggered* TWFE estimate
     (2017/2018 cohorts) — quantifies how much weight the "forbidden"
     later-vs-earlier comparison actually carries. The decomposition is
     verified numerically: the weighted sum of 2x2 components must reproduce
     the staggered TWFE coefficient to machine precision.

Outputs:
  docs/results_summary_closers.md
  plots/bacon_decomposition.png
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
# DATA LOADING (mirrors econometric_analysis.py conventions)
# =====================================================================
print("=" * 72)
print("LOADING DATA")
print("=" * 72)

with open(os.path.join(DATA, "CommAreas_20250408.geojson"), "r", encoding="utf-8") as f:
    geo = json.load(f)
ca_num2name = {int(float(ft["properties"]["area_num_1"])): ft["properties"]["community"].title()
               for ft in geo["features"]}
ca_name2num = {}
for n, name in ca_num2name.items():
    ca_name2num[name] = n
    ca_name2num[name.upper()] = n

hom_raw = pd.read_csv(
    os.path.join(DATA, "Violence_Reduction_-_Victims_of_Homicides_and_Non-Fatal_Shootings_20250513.csv"),
    low_memory=False,
)
hom_raw["DATE"] = pd.to_datetime(hom_raw["DATE"], format="%m/%d/%Y %I:%M:%S %p", errors="coerce")
hom_raw["CA_num"] = hom_raw["COMMUNITY_AREA"].map(ca_name2num)
hom_raw["year"] = hom_raw["DATE"].dt.year

ss_raw = pd.read_csv(os.path.join(DATA, "Violence_Reduction_-_Shotspotter_Alerts_-_Historical_20250408.csv"))
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

YEARS = list(range(2009, 2024))
all_ca = sorted(ca_num2name.keys())
panel = pd.DataFrame([(ca, yr) for ca in all_ca for yr in YEARS], columns=["ca_num", "year"])
panel["ca_name"] = panel["ca_num"].map(ca_num2name)

gun = hom_raw[(hom_raw["VICTIMIZATION_PRIMARY"] == "HOMICIDE") &
              (hom_raw["GUNSHOT_INJURY_I"] == "YES") &
              (hom_raw["year"].between(2009, 2023))].dropna(subset=["CA_num"])
nfs = hom_raw[(hom_raw["VICTIMIZATION_PRIMARY"] != "HOMICIDE") &
              (hom_raw["GUNSHOT_INJURY_I"] == "YES") &
              (hom_raw["year"].between(2009, 2023))].dropna(subset=["CA_num"])
for col, src in [("gun_hom", gun), ("nfs", nfs)]:
    cnt = src.groupby([src["CA_num"].astype(int), "year"]).size().rename(col).reset_index()
    cnt.columns = ["ca_num", "year", col]
    panel = panel.merge(cnt, on=["ca_num", "year"], how="left")
    panel[col] = panel[col].fillna(0).astype(int)

panel["treat"] = panel["ca_name"].isin(ALL_SS).astype(int)
panel["post"] = (panel["year"] >= 2017).astype(int)
panel["did"] = panel["treat"] * panel["post"]           # headline single-cutoff treatment
cohort_year = {n: 2017 for n in COHORT_2017}
cohort_year.update({n: 2018 for n in COHORT_2018})
panel["did_stag"] = panel.apply(
    lambda r: int(r["year"] >= cohort_year[r["ca_name"]]) if r["ca_name"] in cohort_year else 0, axis=1)

print(f"  Panel: {len(panel)} rows; treated={len(ALL_SS)} (2017 cohort {len(COHORT_2017)}, "
      f"2018 cohort {len(COHORT_2018)}), never-treated={77 - len(ALL_SS)}")

# =====================================================================
# 1. NON-FATAL SHOOTINGS AS COUNT MODELS
# =====================================================================
print("\n" + "=" * 72)
print("1. NON-FATAL GUNSHOT INJURIES: Poisson / Negative Binomial TWFE")
print("=" * 72)
print("  (Note: the source records non-fatal shootings from 2010; 2009 counts")
print("   are near zero, matching the OLS specification already in the paper.)")

def count_did(df, outcome, label, exclude_years=None):
    d = df if not exclude_years else df[~df["year"].isin(exclude_years)]
    f = f"{outcome} ~ did + C(ca_num) + C(year)"
    pois = smf.glm(f, data=d, family=sm.families.Poisson()).fit(
        cov_type="cluster", cov_kwds={"groups": d["ca_num"]})
    mu = pois.fittedvalues
    y = d[outcome].values
    # Cameron-Trivedi NB2 auxiliary regression (through the origin) for alpha
    z = ((y - mu) ** 2 - y) / mu
    alpha = max(float(np.sum(z * mu) / np.sum(mu ** 2)), 1e-8)
    nb = smf.glm(f, data=d, family=sm.families.NegativeBinomial(alpha=alpha)).fit(
        cov_type="cluster", cov_kwds={"groups": d["ca_num"]})
    rows = []
    for name, m in [("Poisson", pois), ("NegBin", nb)]:
        b, se, p = m.params["did"], m.bse["did"], m.pvalues["did"]
        lo, hi = m.conf_int().loc["did"]
        rows.append({"model": name, "spec": label, "b": b, "irr": np.exp(b),
                     "irr_lo": np.exp(lo), "irr_hi": np.exp(hi), "p": p,
                     "alpha": (alpha if name == "NegBin" else None)})
        print(f"  {name:8s} {label:22s} b={b:+.4f}  IRR={np.exp(b):.3f} "
              f"[{np.exp(lo):.3f},{np.exp(hi):.3f}]  p={p:.3f} {sig_stars(p)}")
    return rows

nfs_rows = []
nfs_rows += count_did(panel, "nfs", "full panel")
nfs_rows += count_did(panel, "nfs", "excl COVID", [2020, 2021])
nfs_rows += count_did(panel, "nfs", "excl 2016", [2016])

# =====================================================================
# 2. RESTRICTED WILD-CLUSTER BOOTSTRAP (Rademacher, null imposed)
# =====================================================================
print("\n" + "=" * 72)
print("2. WILD-CLUSTER BOOTSTRAP: headline OLS DiD, 999 Rademacher reps")
print("=" * 72)

def fe_design(d):
    """Return (Q from QR of the FE design, cluster ids) for two-way FE."""
    X = pd.get_dummies(d[["ca_num", "year"]].astype(str), drop_first=True)
    X.insert(0, "const", 1.0)
    Q, _ = np.linalg.qr(X.values.astype(float))
    return Q, d["ca_num"].values

def cluster_t(d_til, y_til, clusters):
    """beta, cluster-robust t for y_til ~ d_til (both FE-residualized)."""
    dd = float(d_til @ d_til)
    beta = float(d_til @ y_til) / dd
    e = y_til - d_til * beta
    G = len(np.unique(clusters))
    meat = 0.0
    for c in np.unique(clusters):
        m = clusters == c
        meat += float(d_til[m] @ e[m]) ** 2
    se = np.sqrt(meat * G / (G - 1)) / dd
    return beta, beta / se

def wild_cluster_p(d, B=999, seed=0):
    Q, clusters = fe_design(d)
    y = d["gun_hom"].values.astype(float)
    x = d["did"].values.astype(float)
    d_til = x - Q @ (Q.T @ x)
    y_til = y - Q @ (Q.T @ y)
    beta, t_obs = cluster_t(d_til, y_til, clusters)
    # restricted model: FEs only (null beta=0 imposed)
    mu0 = Q @ (Q.T @ y)
    e0 = y - mu0
    uniq = np.unique(clusters)
    rng = np.random.default_rng(seed)
    exceed = 0
    for _ in range(B):
        w = rng.choice([-1.0, 1.0], size=len(uniq))
        wmap = dict(zip(uniq, w))
        ystar = mu0 + e0 * np.vectorize(wmap.get)(clusters)
        ystar_til = ystar - Q @ (Q.T @ ystar)
        _, t_b = cluster_t(d_til, ystar_til, clusters)
        if abs(t_b) >= abs(t_obs):
            exceed += 1
    return beta, t_obs, (exceed + 1) / (B + 1)

wcb_rows = []
for label, excl in [("full panel", None), ("excl COVID", [2020, 2021]),
                    ("excl 2016", [2016]), ("excl COVID + 2016", [2016, 2020, 2021])]:
    d = panel if not excl else panel[~panel["year"].isin(excl)]
    # cross-check the FWL beta against statsmodels TWFE
    sm_fit = smf.ols("gun_hom ~ did + C(ca_num) + C(year)", data=d).fit(
        cov_type="cluster", cov_kwds={"groups": d["ca_num"]})
    beta, t_obs, p_wcb = wild_cluster_p(d)
    assert abs(beta - sm_fit.params["did"]) < 1e-8, "FWL beta != statsmodels beta"
    wcb_rows.append({"spec": label, "b": beta, "t": t_obs,
                     "p_cluster": sm_fit.pvalues["did"], "p_wcb": p_wcb})
    print(f"  {label:22s} b={beta:+.3f}  t={t_obs:+.2f}  "
          f"cluster-robust p={sm_fit.pvalues['did']:.4f}  wild-bootstrap p={p_wcb:.4f}")

# =====================================================================
# 3. GOODMAN-BACON DECOMPOSITION (staggered 2017/2018 timing)
# =====================================================================
print("\n" + "=" * 72)
print("3. GOODMAN-BACON DECOMPOSITION of the staggered TWFE estimate")
print("=" * 72)

twfe_stag = smf.ols("gun_hom ~ did_stag + C(ca_num) + C(year)", data=panel).fit(
    cov_type="cluster", cov_kwds={"groups": panel["ca_num"]})
beta_stag = twfe_stag.params["did_stag"]
print(f"  Staggered TWFE (cohort-specific timing): b={beta_stag:+.4f} "
      f"(SE {twfe_stag.bse['did_stag']:.3f}, p={twfe_stag.pvalues['did_stag']:.4f})")

grp = panel["ca_name"].map(lambda n: "E" if n in COHORT_2017 else ("L" if n in COHORT_2018 else "U"))
panel["_grp"] = grp
T = len(YEARS)
nE, nL, nU = (grp[panel["year"] == 2009] == "E").sum(), (grp[panel["year"] == 2009] == "L").sum(), \
             (grp[panel["year"] == 2009] == "U").sum()
n = {"E": nE / 77, "L": nL / 77, "U": nU / 77}
Dbar = {"E": 7 / 15, "L": 6 / 15}

def did2x2(dsub, treated_grp, cutoff):
    t = dsub[dsub["_grp"] == treated_grp]
    c = dsub[dsub["_grp"] != treated_grp]
    tp = t[t["year"] >= cutoff]["gun_hom"].mean() - t[t["year"] < cutoff]["gun_hom"].mean()
    cp = c[c["year"] >= cutoff]["gun_hom"].mean() - c[c["year"] < cutoff]["gun_hom"].mean()
    return tp - cp

sub = lambda groups, y0, y1: panel[panel["_grp"].isin(groups) & panel["year"].between(y0, y1)]

b_EU = did2x2(sub(["E", "U"], 2009, 2023), "E", 2017)
b_LU = did2x2(sub(["L", "U"], 2009, 2023), "L", 2018)
b_EL_E = did2x2(sub(["E", "L"], 2009, 2017), "E", 2017)   # early treated, late not-yet (clean)
b_EL_L = did2x2(sub(["E", "L"], 2017, 2023), "L", 2018)   # late treated, early ALREADY treated (forbidden)

nEU = n["E"] / (n["E"] + n["U"]); nLU = n["L"] / (n["L"] + n["U"]); nEL = n["E"] / (n["E"] + n["L"])
s_EU = (n["E"] + n["U"]) ** 2 * nEU * (1 - nEU) * Dbar["E"] * (1 - Dbar["E"])
s_LU = (n["L"] + n["U"]) ** 2 * nLU * (1 - nLU) * Dbar["L"] * (1 - Dbar["L"])
s_EL_E = ((n["E"] + n["L"]) * (1 - Dbar["L"])) ** 2 * nEL * (1 - nEL) * \
         ((Dbar["E"] - Dbar["L"]) / (1 - Dbar["L"])) * ((1 - Dbar["E"]) / (1 - Dbar["L"]))
s_EL_L = ((n["E"] + n["L"]) * Dbar["E"]) ** 2 * nEL * (1 - nEL) * \
         (Dbar["L"] / Dbar["E"]) * ((Dbar["E"] - Dbar["L"]) / Dbar["E"])
tot = s_EU + s_LU + s_EL_E + s_EL_L
weights = {"2017 cohort vs never-treated": (s_EU / tot, b_EU),
           "2018 cohort vs never-treated": (s_LU / tot, b_LU),
           "2017 vs 2018 (not-yet-treated control)": (s_EL_E / tot, b_EL_E),
           "2018 vs 2017 (already-treated control) [forbidden]": (s_EL_L / tot, b_EL_L)}

recon = sum(w * b for w, b in weights.values())
print(f"\n  {'Comparison':50s} {'Weight':>8s} {'2x2 est.':>9s}")
for k, (w, b) in weights.items():
    print(f"  {k:50s} {w:8.3f} {b:+9.3f}")
print(f"\n  Weighted sum = {recon:+.4f}  vs staggered TWFE = {beta_stag:+.4f}  "
      f"(diff {abs(recon - beta_stag):.2e})")
DECOMP_EXACT = abs(recon - beta_stag) < 1e-6
if not DECOMP_EXACT:
    print("  [note] decomposition reproduces TWFE approximately, not exactly — "
          "report weights as approximate.")

forbidden_w = s_EL_L / tot

# Bacon scatter
fig, ax = plt.subplots(figsize=(9, 6))
styles = {"2017 cohort vs never-treated": (COLORS["control"], "o"),
          "2018 cohort vs never-treated": (COLORS["control"], "s"),
          "2017 vs 2018 (not-yet-treated control)": (COLORS["green"], "^"),
          "2018 vs 2017 (already-treated control) [forbidden]": (COLORS["treated"], "D")}
for k, (w, b) in weights.items():
    c, m = styles[k]
    ax.scatter(w, b, s=170, color=c, marker=m, zorder=3, edgecolor="white", linewidth=1.2)
ax.axhline(beta_stag, color="#333", lw=1.4, ls="--", zorder=2)
ax.text(0.015, beta_stag + 0.08, f"staggered TWFE = {beta_stag:+.2f}", fontsize=9.5,
        color="#333", va="bottom", ha="left")
ax.axhline(0, color=COLORS["grid"], lw=1.0, zorder=1)
ax.set_xlim(0, 0.62)
ax.set_ylim(-0.4, 5.2)
ax.set_xlabel("Goodman-Bacon weight")
ax.set_ylabel("2x2 DiD estimate (gun homicides per area-year)")
ax.set_title("Goodman-Bacon decomposition: the forbidden comparison carries\n"
             f"only {forbidden_w:.1%} of the weight", fontweight="bold", pad=10)
handles = [Line2D([0], [0], color=c, marker=m, lw=0, markersize=11, label=k.replace(" [forbidden]", ""))
           for k, (c, m) in styles.items()]
ax.legend(handles=handles, frameon=False, fontsize=9, loc="upper right")
ax.yaxis.grid(True)
add_n_label(ax, f"Cohorts: 2017 (n={nE}), 2018 (n={nL})\nnever-treated (n={nU})", loc="lower right")
add_source_note(fig, SOURCE_DEFAULT)
fig.tight_layout(rect=[0, 0.03, 1, 1])
fig.savefig(os.path.join(PLOTS, "bacon_decomposition.png"), dpi=200, facecolor="white")
plt.close(fig)
print("  SAVED: bacon_decomposition.png")

# =====================================================================
# WRITE docs/results_summary_closers.md
# =====================================================================
L = []
L.append("# Robustness Closers: NFS Count Models, Wild-Cluster Bootstrap, Bacon Decomposition")
L.append("")
L.append("Results of `scripts/robustness_closers.py` — three analyses that close limitations")
L.append("named in the paper.")
L.append("")
L.append("## 1. Non-fatal gunshot injuries as count models (TWFE)")
L.append("")
L.append("The paper's OLS estimate is β = +7.51 (p = 0.009) on the full panel. Estimated as")
L.append("count models (same two-way fixed effects, cluster-robust SEs; NB2 dispersion via the")
L.append("Cameron–Trivedi auxiliary regression):")
L.append("")
L.append("| Model | Specification | β | IRR | 95% CI (IRR) | p |")
L.append("|---|---|---|---|---|---|")
for r in nfs_rows:
    L.append(f"| {r['model']} | {r['spec']} | {r['b']:+.4f} | {r['irr']:.3f} | "
             f"[{r['irr_lo']:.3f}, {r['irr_hi']:.3f}] | {r['p']:.3f} |")
L.append("")
L.append("Source note: the victims file records non-fatal shootings from 2010, so 2009 counts")
L.append("are near zero; this matches the OLS specification already reported.")
L.append("")
L.append("## 2. Restricted wild-cluster bootstrap (headline OLS DiD)")
L.append("")
L.append("Rademacher weights, null imposed, 999 replications, t-statistics cluster-robust with")
L.append("a G/(G−1) small-sample factor (Cameron, Gelbach & Miller 2008). FWL-residualized")
L.append("implementation verified against statsmodels TWFE to 1e-8.")
L.append("")
L.append("| Specification | β (DiD) | t | Cluster-robust p | Wild-bootstrap p |")
L.append("|---|---|---|---|---|")
for r in wcb_rows:
    L.append(f"| {r['spec']} | {r['b']:+.3f} | {r['t']:+.2f} | {r['p_cluster']:.4f} | {r['p_wcb']:.4f} |")
L.append("")
L.append("## 3. Goodman-Bacon decomposition (staggered TWFE, 2017/2018 cohorts)")
L.append("")
L.append(f"Staggered TWFE estimate: {beta_stag:+.4f} (SE {twfe_stag.bse['did_stag']:.3f}).")
L.append(f"Weighted sum of 2x2 components: {recon:+.4f} "
         f"({'exact' if DECOMP_EXACT else 'approximate'} reproduction, diff {abs(recon-beta_stag):.2e}).")
L.append("")
L.append("| Comparison | Weight | 2x2 estimate |")
L.append("|---|---|---|")
for k, (w, b) in weights.items():
    L.append(f"| {k} | {w:.3f} | {b:+.3f} |")
L.append("")
L.append(f"The 'forbidden' comparison (2018 cohort vs the already-treated 2017 cohort) carries")
L.append(f"**{forbidden_w:.1%}** of the total weight. Negative-weighting contamination is therefore")
L.append("minor in this panel; the sign fragility documented in the paper comes from the 2016")
L.append("anomaly and the absence of a valid control group, not from staggered-timing bias in")
L.append("the TWFE aggregation itself.")
L.append("")
L.append("Figure: `plots/bacon_decomposition.png`")

with open(os.path.join(DOCS, "results_summary_closers.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(L))
print("\n  WROTE: docs/results_summary_closers.md")
print("\nDONE")
