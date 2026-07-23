# Detection Without Identification: ShotSpotter & Gun Homicide in Chicago

![Python](https://img.shields.io/badge/Python-3.12-blue.svg)
![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)
![Methods](https://img.shields.io/badge/methods-DiD%20%7C%20event%20study%20%7C%20synthetic%20control%20%7C%20placement%20model-orange.svg)

**The question.** Chicago's ShotSpotter gunshot-detection system was deployed across
51 of 77 community areas in 2017–2018 and cancelled in 2023. That gunshot detection
*doesn't* reduce gun violence is, by now, well established, including for Chicago at the
police-district level (Connealy et al. 2024, *J. Experimental Criminology*) and nationally
(Doucette et al. 2021, *J. Urban Health*), with a 2026 meta-analysis of 44 estimates
pooling to a null (Huff, Dunlap & Pearson 2026, *Am. J. Criminal Justice*). So this project asks a sharper, methodological
question: **can the community-area data that most local policy argument relies on even
*identify* such an effect, and what happens when a standard evaluation is run on them
without asking?**

**The finding.** It cannot, and the way a routine evaluation fails is the point. A naive
difference-in-differences says ShotSpotter areas saw **more** gun homicides (+2.6 per
area-year, p < 0.001). But that headline is a statistical artifact: it collapses to a null
once homicide *counts* are modeled correctly, it fails a parallel-trends test, it
reproduces equally "significant" effects when treatment is faked in pre-rollout years, it
has no valid synthetic-control comparison, and even the modern Callaway–Sant'Anna estimator
only turns it negative because of a single outlier year. With no clean control group and a
pre-period dominated by one anomalous year, these data cannot identify a causal effect in
either direction. This is a worked cautionary case in the LaLonde (1986) tradition: a
coefficient can survive fixed effects and robustness restrictions and still fail every
condition a causal reading requires.

**What the data *do* show.** A placement model pins down *why* no counterfactual exists:
ShotSpotter was sited where violence, structural disadvantage, and (even conditional on
those) the **Black and Hispanic share** of a neighborhood were highest, and so tightly
that no comparable untreated neighborhood remains. The design rules out a large *deterrence*
effect, but cannot exclude the survival-channel benefit a boundary regression discontinuity
attributes to faster emergency response, a magnitude that sits *inside* the confidence
interval this panel can resolve.

**How we get there** (and how you can reproduce it from the data and code below) is a
sequence of increasingly demanding causal-inference designs on a community-area × year
panel of Chicago gun homicides, 2009–2023, plus a cross-sectional model of *treatment
assignment* itself. The figures below walk the evidence in order; the formatted paper is
[paper/final_paper.pdf](paper/final_paper.pdf) and [NARRATIVE.md](NARRATIVE.md) is the
long-form write-up.

![OLS headline vs. count-data specifications](plots/specification_comparison.png)

*The whole story in one figure. The OLS effect is the additive +2.63 coefficient shown
on the multiplicative scale for comparison (implied IRR ≈ 1.35, p < 0.001); the Poisson
and Negative Binomial points are genuinely-estimated incidence-rate ratios. Both
count-data models collapse the headline to a statistically null IRR ≈ 1.07.*

## Headline result at a glance

| Specification | Estimate | 95% CI | p |
|---|---|---|---|
| **OLS DiD** | +2.63 gun homicides / area-year | [+1.38, +3.88] | **<0.001*** |
| **Poisson DiD** (two-way fixed effects) | IRR 1.07 | [0.87, 1.33] | 0.51 (ns) |
| **Negative Binomial DiD** (dispersion estimated) | IRR 1.07 | [0.86, 1.33] | 0.55 (ns) |
| **ETWFE-Poisson** (Wooldridge 2023, heterogeneity-robust) | IRR 1.05 | [0.79, 1.33] | ns |

The +2.63 OLS coefficient is *robust* to adding community-area and year fixed
effects; it does **not** come from a missing-controls problem. It comes from
modeling integer homicide counts as a continuous, homoscedastic outcome. On the
natural multiplicative (count-data) scale, the effect is statistically null.

The +2.63 collapses to the same null under the heterogeneity-robust count estimator the
staggered-rollout literature calls for (a Wooldridge 2023 ETWFE-Poisson), so it is not an
artifact of the plain Poisson's fixed-effects weighting.

The full write-up is in **[NARRATIVE.md](NARRATIVE.md)**; the formatted paper is
**[paper/final_paper.pdf](paper/final_paper.pdf)**. Regenerable result tables live in
[docs/results_summary_original.md](docs/results_summary_original.md) (original pipeline),
[docs/results_summary_reanalysis.md](docs/results_summary_reanalysis.md)
(re-analysis + robustness),
[docs/results_summary_placement.md](docs/results_summary_placement.md)
(treatment-assignment / equity model),
[docs/results_summary_etwfe.md](docs/results_summary_etwfe.md) (ETWFE robustness), and
[docs/results_summary_removal.md](docs/results_summary_removal.md) (the September 2024
removal experiment).

## Methods & skills demonstrated

- **Causal inference:** difference-in-differences (pooled and two-way fixed
  effects), Callaway–Sant'Anna heterogeneity-robust staggered DiD, cohort-specific
  event study, synthetic control with in-space permutation inference, pre-period
  placebo / falsification tests, and a second natural experiment exploiting the
  September 2024 ShotSpotter shutoff.
- **Selection / treatment-assignment modeling:** cross-sectional linear-probability
  and logistic models of *which* neighborhoods got ShotSpotter, standardized-effect
  balance tables, propensity-score overlap / common-support diagnostics, and an equity
  test of whether racial composition predicts placement conditional on violence and
  disadvantage.
- **Count-data econometrics:** Poisson and Negative Binomial regression (with the
  dispersion estimated), incidence-rate ratios, clustered standard errors, and a
  Wooldridge (2023) extended two-way fixed-effects Poisson (the heterogeneity-robust
  count analog of Callaway–Sant'Anna).
- **Robustness & diagnostics:** joint parallel-trends Wald test, COVID- and
  2016-sensitivity checks, multiple-specification forest plots.
- **Geospatial analysis:** GeoPandas, choropleth and kernel-density hotspot maps,
  CRS-correct metric grids (EPSG:26971), `contextily` basemaps.
- **Python data stack:** pandas, NumPy, statsmodels, SciPy, matplotlib, seaborn.
- **Reproducible research:** a relative-path pipeline that runs from any clone,
  scripted figure generation, and documented data provenance.

## How the analysis works

### 1. The setting: ShotSpotter was deployed where violence already was

![Treatment vs. control community areas](plots/did_map_treatment_control.png)

ShotSpotter coverage (blue) blankets the high-violence South and West sides, while the
never-treated areas sit elsewhere in the city. Before deployment, the 51 covered areas
averaged roughly **4.6 times** the gun homicides *per area* of the 26 uncovered ones
(about nine times in aggregate, since there are nearly twice as many covered areas).
Treatment was assigned *because of* violence, not at random, which is the central
problem every design here has to confront, and the reason the comparison leans
entirely on the parallel-trends assumption.

### 2. As a detector, the technology works

![Hourly pattern of alerts and homicides](plots/hourly_pattern.png)

ShotSpotter alerts and gun homicides rise and fall together over the day (both
bottoming out mid-morning and peaking late at night) and they peak in the same summer
months. The system detects the gunfire it is designed to detect. The question is
whether *detection* translated into *fewer homicides*.

### 3. The headline says "more homicides", but it is a modeling artifact

The two-way fixed-effects OLS DiD returns **+2.63 gun homicides per area-year
(p < 0.001)**. That number is robust to fixed effects, COVID exclusions, and dropping
2016, so it looks bulletproof. It is not. Annual homicide counts are small
non-negative integers whose variance grows with their mean; OLS treats them as
continuous and equal-variance. Estimated with the right model for counts (Poisson and
Negative Binomial with the same fixed effects) the effect is **statistically null**: a
small, non-significant *positive* point estimate (IRR ≈ 1.07, p ≈ 0.5; see the headline
figure above), not the dramatic +2.6 the additive scale implies. The count models agree
with OLS in *direction* but shrink the effect to noise on the natural multiplicative
scale. And the null is not an artifact of the plain Poisson's fixed-effects weighting under
the staggered rollout: the heterogeneity-robust count estimator the literature prescribes
(a Wooldridge (2023) extended two-way fixed-effects Poisson) returns the same null (overall
ATT IRR 1.05, 95% CI [0.79, 1.33]), with cohort-time effects at one at every horizon.

### 4. Parallel trends fail

![Cohort-aware event study](plots/event_study_cohort_aware.png)

In this cohort event study every pre-treatment coefficient is large and negative, an
artifact of 2016 (a Chicago-wide homicide-record year) sitting next to the base period.
The robust evidence is a joint Wald test on the pre-treatment coefficients, which is
invariant to the base year and rejects parallel trends (F = 19.9, p = 0.006). Once the
design fails this test, the DiD contrast cannot be read as causal.

### 5. Fake treatment dates produce the same "effect"

![Falsification / placebo test](plots/falsification_plot.png)

If +2.6 were a real treatment effect, assigning *fake* treatment dates years before
ShotSpotter existed should yield nothing. Instead, the 1999 and 2003 placebos come back
statistically significant (p < 0.01) and of comparable magnitude (β = −1.34 and −1.93).
They are opposite-signed, but that is the point: a design that manufactures "significant"
effects in years before the program existed is capturing pre-existing trend differences
across neighborhoods, not the intervention.

### 6. Placement was predictable, and it tracked race, not just violence

![What distinguished ShotSpotter areas](plots/placement_coefficients.png)

Why is there no counterfactual? Because treatment was *nearly deterministic*. A
cross-sectional model of ShotSpotter placement across all 77 community areas shows that
structural disadvantage, the hardship index (AUC 0.93) and low per-capita income (0.92),
separates covered from uncovered areas even more sharply than the pre-period gun-homicide
*count* itself (0.79). And placement mapped onto Chicago's racial geography: ShotSpotter
areas are on average **49% Black and 32% Hispanic, versus 14% and 17%** in never-treated
areas, and (conditional on violence, hardship, and income) a one-SD increase in the Black
share still raises the probability of placement by +33 points (p = 0.002) and the Hispanic
share by +28 (p = 0.003). Race, hardship, and violence are deeply entangled in a segregated
city, so this is a descriptive pattern, not proof that race drove siting *independently* of
disadvantage, but it is the concentrated-surveillance pattern the MacArthur Justice Center
and the city's Inspector General already documented; our marginal addition is the
*conditional* result (race predicts placement even after violence and hardship are held fixed).

![Propensity-score overlap](plots/placement_propensity_overlap.png)

The consequence for identification is stark. A propensity model separates the two groups
with an AUC of 0.94, and the six highest-violence study neighborhoods all sit at a placement
propensity of **≥ 0.97, a band no never-treated area reaches.** This is the selection-side
view of the next step: the city installed ShotSpotter in essentially every highest-violence,
highest-hardship neighborhood, leaving no comparable untreated unit from which to build a
counterfactual. (Racial composition is measured from the 2019–2023 ACS, a proxy for the
persistent pre-rollout composition.)

### 7. There is no valid comparison group (synthetic control)

![Synthetic control panels](plots/synthetic_control_panels.png)

For every study neighborhood, the synthetic counterfactual (red) sits far below the
actual series (blue); even in the pre-treatment window it was explicitly fit on. No
weighted combination of the 26 never-treated, lower-violence areas can reproduce a
high-violence treated unit, so the optimizer collapses onto a single donor and the
method is **infeasible**. The city left no comparable untreated neighborhood; the
absence of a counterfactual is itself the finding.

### 8. The "modern estimator" sign rides on a single year

![Callaway–Sant'Anna event study](plots/callaway_santanna_event_study.png)

The gold-standard estimator for staggered rollouts, Callaway & Sant'Anna (2021),
returns an overall ATT of **−2.45 (p < 0.01)** on the full panel, the *opposite* sign of
OLS. But that negative is not a property of the method: it is driven almost entirely by
the 2016 record-homicide year sitting at the *k* = −1 baseline. **Exclude 2016 and the
same estimator returns +1.90 (SE 0.72), back to the sign of OLS.** The event study
makes this visible: only the *k* = −1 (2016) pre-coefficient is significant; the rest
hover at zero.

So read honestly, the estimators do not show a mysterious "method-determined sign":

| Estimator | Full panel | What actually moves it |
|---|---|---|
| Naive OLS DiD | **+2.63** (p < 0.001) | additive scale, inflated by high-count areas |
| Poisson / Neg. Binomial | IRR ≈ **1.07** (ns) | positive in sign, but null on the multiplicative scale |
| Callaway–Sant'Anna | **−2.45** (p < 0.01) | flips to **+1.90** once the 2016 outlier year is dropped |

OLS and the count models actually agree in *direction* (both positive); only Callaway–
Sant'Anna goes negative, and only because of one anomalous year. The real lesson is not
"any answer is possible," but that every estimate is dominated by a single outlier year
on top of a design with **no valid control group** (Step 7) and a **failed pre-trend
test** (Step 4). With no clean counterfactual, the causal effect is simply **not
identifiable from these data**, in either direction.

### 9. A second natural experiment: the September 2024 removal

![The September 2024 removal](plots/removal_trajectory.png)

Chicago switched ShotSpotter off on 22 September 2024, and the gun-homicide decline that
followed has been read in public commentary as proof the technology never worked. Applying
the same discipline to the *removal* shows it is no more identifiable than the installation,
and it fails in the same two ways, which is the cleanest confirmation of the whole argument.
Gun homicides were already falling city-wide from the 2021–22 peak, in ShotSpotter and
never-treated areas alike (coverage areas **−37%**, never-treated **−31%** over the same
post-removal months), so the raw "crime fell after ShotSpotter left" is a city-wide trend,
not a removal effect. And the additive removal DiD even manufactures a "significant" benefit
(**−0.25 homicides/area-month, p < 0.001**) that vanishes on the multiplicative scale
(**Poisson IRR 0.91, 95% CI [0.54, 1.52]**), the identical scale artifact as the
installation. Both the raw drop and the additive DiD are mirages; the public debate has
confidently drawn both unwarranted inferences.

## Bottom line

ShotSpotter did what it was built to do: *detect* gunfire that 911 calls miss. But
across every design that respects the data, there is **no credible evidence of an effect
on gun homicides in either direction**. The honest reading is neither "ShotSpotter
increased violence" (the +2.6 is an artifact of the additive scale) nor "ShotSpotter cut
violence" (that reads a 2016-driven CS estimate too literally); it is that, with no
valid control group and a pre-period dominated by one outlier year, **these data simply
cannot identify a causal effect**. They rule out a large *deterrence* effect (a broad,
sustained drop in shootings of the kind a $3M detection contract was meant to deliver) but
not the smaller, survival-channel mortality benefit a boundary regression discontinuity
attributes to faster emergency response, a magnitude that falls *inside* the interval this
panel can resolve. That non-result is itself the policy finding, and it is consistent with
the city's 2023 decision to end the contract. The 2024 removal, analyzed the same way
(Step 9), reproduces both failures, so the non-identification is a property of the setting,
not of how the installation happens to be timed. Full argument and limitations:
[NARRATIVE.md](NARRATIVE.md).

## Limitations & next steps

We'd rather name the soft spots than have a reviewer find them:

- **The identification ceiling is the data, not just the method.** Chicago deployed
  ShotSpotter in *every* high-violence neighborhood, so there is no untreated unit at
  the treated units' violence level. No estimator can manufacture a counterfactual the
  data don't contain, which is why the synthetic control is reported as *infeasible*
  rather than as a result.
- **The placement model is descriptive.** It documents *that* siting tracked disadvantage,
  violence, and racial composition, and quantifies the resulting lack of common support,
  but with 77 heavily collinear units it does not identify a causal assignment rule, and its
  racial-composition measure is the 2019–2023 ACS used as a proxy for the (highly
  persistent) pre-rollout composition.
- **Two event-study implementations.** The cohort event study in
  `econometric_analysis.py` is a transparent hand-rolled difference-in-means with a
  bootstrap; `did_callaway_santanna.py` is the packaged, heterogeneity-robust
  Callaway–Sant'Anna estimator and should be treated as the authoritative
  staggered-DiD result. The former is kept for intuition, not as the final word.
- **Inference.** Standard errors cluster on community area (~77 clusters); a restricted
  wild-cluster bootstrap (999 Rademacher replications) leaves the OLS DiD significant
  (p ≈ 0.001) in every specification; the headline artifact is a *specification* problem,
  not a standard-error problem. A Goodman–Bacon decomposition likewise shows the
  "forbidden" already-treated comparison carries only 4.6% of the staggered TWFE weight,
  so the fragility is the 2016 anomaly and the missing control group, not timing bias.
- **Non-fatal shootings**, re-estimated as count models, tell the same story even more
  starkly: the OLS +7.51 (p = 0.009) flips *below* one and goes statistically null
  (IRR 0.84–0.92, p ≥ 0.23).
- **Scope.** The panel is annual (masking within-year dynamics), and the analysis
  observes neither police response times, dispatch decisions, nor arrests: the actual
  channel through which detection could plausibly reduce violence.

## Repository layout

```
.
├── README.md
├── NARRATIVE.md                       Full write-up of findings and caveats
├── LICENSE
├── requirements.txt
│
├── data/
│   ├── raw/                           City of Chicago portal inputs (large files
│   │                                  git-ignored, see "Data" below)
│   └── processed/                     Cleaned subsets produced by notebooks/clean_*.ipynb
│
├── notebooks/                         Jupyter notebooks (cleaning + spatial work)
│   ├── clean_homicides_shotspotter.ipynb
│   ├── clean_poverty.ipynb
│   ├── heatmap.ipynb
│   └── murder_map.ipynb
│
├── scripts/                           Standalone Python analysis pipeline
│   ├── plot_style.py                  Shared academic-paper style (palette, helpers)
│   ├── did_analysis.py                OLS DiD, event study, t-tests
│   ├── econometric_analysis.py        TWFE DiD, cohort event study, synthetic control
│   ├── econometric_robustness.py      Poisson/NegBin DiD, pre-trend Wald, permutation SC
│   ├── did_callaway_santanna.py       Callaway–Sant'Anna heterogeneity-robust DiD
│   ├── placement_model.py             Treatment-assignment / equity model + propensity overlap
│   ├── robustness_closers.py          NFS count models, wild-cluster bootstrap, Bacon decomposition
│   ├── procedural_benefit.py          Process-vs-outcome test: weapons enforcement DiD (Crimes dataset)
│   ├── etwfe_count.py                 Wooldridge (2023) ETWFE-Poisson: heterogeneity-robust count ATT
│   ├── removal_analysis.py            Sept 2024 shutoff as a second natural experiment (non-identifiable)
│   ├── make_narrative_plots.py        Specification-comparison + falsification figures
│   ├── make_pretty_plots.py           Descriptive figures (monthly, hourly, hotspots)
│   └── plot_covid_sensitivity.py
│
├── plots/                             All generated figures (PNG, 200 DPI)
├── paper/
│   ├── final_paper.tex / .pdf         Formatted paper (LaTeX source + PDF)
│   └── methods_memo.tex / .pdf        ~8-page methods-memo version
└── docs/
    ├── results_summary_original.md    Original DiD tables (from did_analysis.py)
    ├── results_summary_reanalysis.md  Re-analysis + robustness tables (econometric_*.py)
    ├── results_summary_placement.md   Placement / equity model tables (placement_model.py)
    ├── results_summary_closers.md     NFS count models, wild bootstrap, Bacon (robustness_closers.py)
    ├── results_summary_procedural.md   Process-vs-outcome / weapons-enforcement DiD (procedural_benefit.py)
    ├── results_summary_etwfe.md        Wooldridge ETWFE-Poisson robustness (etwfe_count.py)
    ├── results_summary_removal.md       Sept 2024 removal second natural experiment (removal_analysis.py)
    ├── Pov_data_table.docx            Poverty summary table (generated by clean_poverty.ipynb)
    └── shotspotter_analysis.docx      Original course write-up
```

## Data

All inputs are public products of the [Chicago Data Portal](https://data.cityofchicago.org/):

| File (in `data/raw/`) | Source | In repo? |
|---|---|---|
| `Violence_Reduction_-_Victims_of_Homicides_and_Non-Fatal_Shootings_*.csv` | [Victims of Homicides & Non-Fatal Shootings](https://data.cityofchicago.org/Public-Safety/Violence-Reduction-Victims-of-Homicides-and-Non-Fa/gumc-mgzr/about_data) | No (download) |
| `Violence_Reduction_-_Shotspotter_Alerts_-_Historical_*.csv` | [ShotSpotter Alerts (Historical)](https://data.cityofchicago.org/Public-Safety/Violence-Reduction-Shotspotter-Alerts-Historical/3h7q-7mdb/about_data) | No (download) |
| `transportation_*.csv` | [Transportation network](https://data.cityofchicago.org/Transportation/transportation/7ez8-272k/about_data) (street overlays only) | No (download) |
| `CommAreas_*.geojson` | [Community Areas](https://data.cityofchicago.org/Facilities-Geographic-Boundaries/Boundaries-Community-Areas-current-/cauq-8yn6) | Yes |
| `Census_Data_-_Selected_socioeconomic_indicators_*.csv` | [Census: socioeconomic indicators](https://data.cityofchicago.org/Health-Human-Services/Census-Data-Selected-socioeconomic-indicators-in-C/kn9c-c2s2/about_data) | Yes |
| `ACS_5yr_race_by_community_area_2023.csv` | [ACS 5-Year Data by Community Area](https://data.cityofchicago.org/d/t68z-cikk) (racial composition for the placement model) | Yes |
| `Public_Health_Statistics_Life_Expectancy_By_Community_Area.csv` | [Life Expectancy by Community Area](https://data.cityofchicago.org/d/qjr3-bm53) (2010, equity measure for the placement model) | Yes |
| `Crimes_weapons_violations_by_ca_year.csv` | [Crimes – 2001 to Present](https://data.cityofchicago.org/d/ijzp-q8t2) (weapons-violation incidents/arrests per CA-year, for the procedural-channel analysis) | Yes |

The three large raw files (~110 MB combined) are **git-ignored** to keep the
repository lightweight and to avoid re-hosting victim-level data. They are freely
re-downloadable from the links above; drop them into `data/raw/`. The analysis never
uses the victim-name columns; the cleaned files in `data/processed/` contain only
aggregated, de-identified fields.

## Reproducing the analysis

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Download the three git-ignored raw files into data/raw/ (see "Data" above).

# 3. (Optional) regenerate the cleaned subsets in data/processed/ by running the
#    notebooks under notebooks/. The committed copies are already sufficient for
#    the descriptive scripts; the regression scripts read directly from data/raw/.

# 4. Original DiD pipeline (writes docs/results_summary_original.md + figures)
python scripts/did_analysis.py

# 5. Re-analysis + robustness (TWFE, count-data, synthetic control,
#    pre-period falsification, permutation inference, Callaway–Sant'Anna)
python scripts/econometric_analysis.py
python scripts/econometric_robustness.py
python scripts/did_callaway_santanna.py    # needs the `differences` package

# 6. Placement / equity model (treatment-assignment LPM + logit, propensity
#    overlap, racial-composition equity test)
python scripts/placement_model.py
python scripts/robustness_closers.py   # NFS count models, wild-cluster bootstrap, Bacon decomposition
python scripts/procedural_benefit.py   # process-vs-outcome: weapons-enforcement DiD (Crimes dataset)
python scripts/etwfe_count.py          # Wooldridge ETWFE-Poisson (heterogeneity-robust count ATT)
python scripts/removal_analysis.py     # Sept 2024 shutoff as a second natural experiment

# 7. Narrative + descriptive figures
python scripts/make_narrative_plots.py
python scripts/make_pretty_plots.py
python scripts/plot_covid_sensitivity.py
```

Every script resolves paths relative to its own location, so the pipeline runs from any
clone with no path edits. (The hotspot maps fetch basemap tiles over the network via
`contextily`; they degrade gracefully to no-basemap when offline.)

## Contact

Austin Belman, abelma2@illinois.edu (repository maintainer)

## License

MIT, see [LICENSE](LICENSE).
