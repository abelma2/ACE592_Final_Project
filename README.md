# ShotSpotter & Gun Homicide in Chicago — A Difference-in-Differences Study

![Python](https://img.shields.io/badge/Python-3.12-blue.svg)
![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)
![Methods](https://img.shields.io/badge/methods-DiD%20%7C%20event%20study%20%7C%20synthetic%20control-orange.svg)

**The question.** Did Chicago's ShotSpotter gunshot-detection system — deployed
across 51 of 77 community areas in 2017–2018 and cancelled in 2023 — actually
reduce gun homicides?

**The finding.** No — and, more importantly, *the data cannot credibly show that it
did or didn't.* A naive difference-in-differences says ShotSpotter areas saw **more**
gun homicides (+2.6 per area-year, p < 0.001). But that headline is a statistical
artifact: it evaporates the moment homicide *counts* are modeled correctly, it fails
a parallel-trends test, it reproduces equally "significant" effects when treatment is
faked in pre-rollout years, and it has no valid synthetic-control comparison.
ShotSpotter worked as a *detector* of gunfire; the claim that it *reduced* gun
homicides does not survive scrutiny.

**How I get there** — and how you can reproduce it from the data and code below — is
a sequence of five increasingly demanding causal-inference designs on a
community-area × year panel of Chicago gun homicides, 2009–2023. The rest of this
README walks the evidence; [NARRATIVE.md](NARRATIVE.md) is the full write-up.

![OLS headline vs. count-data specifications](plots/specification_comparison.png)

*The OLS headline (+7.4% implied, p < 0.001) collapses to a statistically null
incidence-rate ratio under both proper count-data models.*

## Headline result at a glance

| Specification | Estimate | 95% CI | p |
|---|---|---|---|
| **OLS DiD** | +2.63 gun homicides / area-year | [+1.38, +3.88] | **<0.001*** |
| **Poisson DiD** (two-way fixed effects) | IRR 1.07 | [0.87, 1.33] | 0.51 (ns) |
| **Negative Binomial DiD** | IRR 0.97 | [0.72, 1.31] | 0.84 (ns) |

The +2.63 OLS coefficient is *robust* to adding community-area and year fixed
effects — it does **not** come from a missing-controls problem. It comes from
modeling integer homicide counts as a continuous, homoscedastic outcome. On the
natural multiplicative (count-data) scale, the effect is statistically null.

The full write-up is in **[NARRATIVE.md](NARRATIVE.md)**; the formatted paper is
**[paper/final_paper.pdf](paper/final_paper.pdf)**. Regenerable result tables live
in [docs/results_summary_original.md](docs/results_summary_original.md) (original
pipeline) and [docs/results_summary_reanalysis.md](docs/results_summary_reanalysis.md)
(re-analysis + robustness).

## Methods & skills demonstrated

- **Causal inference:** difference-in-differences (pooled and two-way fixed
  effects), cohort-specific staggered-rollout event study, synthetic control
  with in-space permutation inference, pre-period placebo / falsification tests.
- **Count-data econometrics:** Poisson and Negative Binomial regression,
  incidence-rate ratios, clustered standard errors.
- **Robustness & diagnostics:** joint parallel-trends Wald test, COVID- and
  2016-sensitivity checks, multiple-specification forest plots.
- **Geospatial analysis:** GeoPandas, choropleth and kernel-density hotspot maps,
  CRS-correct metric grids (EPSG:26971), `contextily` basemaps.
- **Python data stack:** pandas, NumPy, statsmodels, SciPy, matplotlib, seaborn.
- **Reproducible research:** a relative-path pipeline that runs from any clone,
  scripted figure generation, and documented data provenance.

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
│   │                                  git-ignored — see "Data" below)
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
│   ├── make_narrative_plots.py        Specification-comparison + falsification figures
│   ├── make_pretty_plots.py           Descriptive figures (monthly, hourly, hotspots)
│   └── plot_covid_sensitivity.py
│
├── plots/                             All generated figures (PNG, 200 DPI)
├── paper/
│   └── final_paper.pdf                Formatted paper
└── docs/
    ├── results_summary_original.md    Original DiD tables (from did_analysis.py)
    ├── results_summary_reanalysis.md  Re-analysis + robustness tables (econometric_*.py)
    ├── Pov_data_table.docx            Poverty summary table (generated by clean_poverty.ipynb)
    └── shotspotter_analysis.docx      Original course write-up
```

## Data

All inputs are public products of the [Chicago Data Portal](https://data.cityofchicago.org/):

| File (in `data/raw/`) | Source | In repo? |
|---|---|---|
| `Violence_Reduction_-_Victims_of_Homicides_and_Non-Fatal_Shootings_*.csv` | [Victims of Homicides & Non-Fatal Shootings](https://data.cityofchicago.org/Public-Safety/Violence-Reduction-Victims-of-Homicides-and-Non-Fa/gumc-mgzr/about_data) | No — download |
| `Violence_Reduction_-_Shotspotter_Alerts_-_Historical_*.csv` | [ShotSpotter Alerts (Historical)](https://data.cityofchicago.org/Public-Safety/Violence-Reduction-Shotspotter-Alerts-Historical/3h7q-7mdb/about_data) | No — download |
| `transportation_*.csv` | [Transportation network](https://data.cityofchicago.org/Transportation/transportation/7ez8-272k/about_data) (street overlays only) | No — download |
| `CommAreas_*.geojson` | [Community Areas](https://data.cityofchicago.org/Facilities-Geographic-Boundaries/Boundaries-Community-Areas-current-/cauq-8yn6) | Yes |
| `Census_Data_-_Selected_socioeconomic_indicators_*.csv` | [Census — socioeconomic indicators](https://data.cityofchicago.org/Health-Human-Services/Census-Data-Selected-socioeconomic-indicators-in-C/kn9c-c2s2/about_data) | Yes |

The three large raw files (~110 MB combined) are **git-ignored** to keep the
repository lightweight and to avoid re-hosting victim-level data. They are freely
re-downloadable from the links above — drop them into `data/raw/` (the cleaned
subsets the analysis filters from are derived from these). The analysis itself
never uses the victim-name columns; the cleaned files in `data/processed/` contain
only aggregated, de-identified fields.

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
#    pre-period falsification, permutation inference)
python scripts/econometric_analysis.py
python scripts/econometric_robustness.py

# 6. Narrative + descriptive figures
python scripts/make_narrative_plots.py
python scripts/make_pretty_plots.py
python scripts/plot_covid_sensitivity.py
```

Every script resolves paths relative to its own location, so the pipeline runs
from any clone with no path edits. (The hotspot maps fetch basemap tiles over the
network via `contextily`; they degrade gracefully to no-basemap when offline.)

## Why the headline result doesn't hold

![Cohort-aware event study](plots/event_study_cohort_aware.png)

*Every pre-treatment coefficient is large and negative — a clear parallel-trends
violation driven by 2016 (a Chicago-wide homicide-record year) sitting next to the
event-study base period. A joint Wald test rejects parallel trends (F = 19.9,
p = 0.006).*

The +2.6 OLS coefficient survives unit and year fixed effects, so the problem is
not omitted controls. The problem is everything else — the model family, the
staggered rollout collapsed into a single 2017 cutoff, and an anomalous base year:

- **Count-data models are null** (Poisson IRR 1.07, p = 0.51; NegBin IRR 0.97, p = 0.84).
- **Placebos are "significant"** — fake 1999/2003 treatment dates produce effects of
  similar magnitude, so the design is picking up pre-existing trends, not treatment.
- **Parallel trends fail** — the joint Wald test rejects (F = 19.9, p = 0.006).
- **Synthetic control is infeasible** — the 26 never-treated areas are all
  lower-violence than any treated neighborhood, so there is no valid donor pool.

## Bottom line

ShotSpotter did what it was built to do — *detect* gunfire that 911 calls miss. But
across every design that respects the data, there is **no evidence it reduced gun
homicides**, in either direction. The honest reading is neither "ShotSpotter
increased violence" (the +2.6 is an artifact of the wrong model) nor "ShotSpotter cut
violence" (that ignores an anomalous 2016 base year) — it is that **these data simply
cannot identify a causal effect**, while ruling out the large reduction a $3M
detection contract was meant to deliver. That null is itself the policy finding, and
it is consistent with the city's 2023 decision to end the contract. Full argument and
limitations: [NARRATIVE.md](NARRATIVE.md).

## Contact

Austin Belman — abelma2@illinois.edu

## License

MIT — see [LICENSE](LICENSE).
