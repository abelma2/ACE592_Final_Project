# ShotSpotter Placement Model: What Predicted Treatment?

Results of `scripts/placement_model.py`. This is the *selection* / first-stage
analysis: a cross-section of all 77 Chicago community areas in which ShotSpotter
placement is modeled on pre-period (2009-2016) gun-homicide levels and the
2008-2012 Census socioeconomic indicators. It converts the paper's asserted
selection-on-violence story into a directly estimated fact, and shows the
selection-side reason no valid counterfactual exists.

Sample: 77 community areas (51 ShotSpotter, 26 never-treated). Treatment defined as in the main analysis (≥100 alerts and ≥12 months of coverage).

## 1. Balance table (treated vs. never-treated)

`Std. diff` is the standardized mean difference (treated − control, in pooled SDs). `AUC` is the direction-free area under the ROC curve for that single characteristic as a classifier of placement (0.5 = no discriminating power, 1.0 = perfect).

| Characteristic | Treated mean | Control mean | Std. diff | Univariate AUC |
|---|---|---|---|---|
| Hardship index | 63.90 | 21.27 | +2.09 | 0.931 |
| Per-capita income ($000s) | 18.65 | 39.11 | -1.72 | 0.919 |
| Unemployed, 16+ (%) | 18.70 | 8.84 | +1.66 | 0.919 |
| Households below poverty (%) | 26.44 | 12.60 | +1.45 | 0.885 |
| No high-school diploma, 25+ (%) | 24.90 | 11.40 | +1.35 | 0.860 |
| Life expectancy, 2010 (years) | 76.16 | 80.43 | -1.20 | 0.806 |
| Pre-period gun homicides (mean/yr, 2009-16) | 7.48 | 1.62 | +0.89 | 0.793 |
| Dependency: under 18 or over 64 (%) | 38.42 | 30.51 | +1.25 | 0.787 |
| Crowded housing (%) | 5.84 | 3.12 | +0.78 | 0.745 |

The health penalty is the most human summary of the gap: ShotSpotter areas averaged **76.2 years** of life expectancy in 2010 versus **80.4** in never-treated areas---a **-4.3-year** difference. Life expectancy is a descriptive disadvantage summary here, not a siting predictor (the city did not site on it); it is kept out of the placement regressions below and reported only to characterize who bears the coverage. Source: Chicago Public Health Statistics, Life Expectancy by Community Area (`qjr3-bm53`).

## 2. Placement model and the violence horse-race

Linear-probability models (HC1 robust SEs) with standardized predictors, so each
coefficient is the change in the probability of placement per one-SD increase.
The logit AUC and pseudo-R² come from a near-unpenalized logistic fit on the same
predictors. The question: does anything predict placement *on top of* violence?

| Model | LPM R² | Logit AUC | Pseudo-R² |
|---|---|---|---|
| Violence only | 0.155 | 0.793 | 0.197 |
| Disadvantage only (hardship+income) | 0.505 | 0.934 | 0.512 |
| Violence + hardship + income | 0.506 | 0.935 | 0.519 |
| Full (violence + all SES) | 0.541 | 0.959 | 0.624 |

Standardized coefficients, parsimonious model (`Violence + hardship + income`):

| Predictor | b (Δ prob per SD) | SE | p |
|---|---|---|---|
| z_pre_gun_hom | +0.0157 | 0.0263 | 0.5495 |
| z_hardship | +0.2713 | 0.0675 | 0.0001 |
| z_income | -0.0646 | 0.0527 | 0.2204 |

**Headline.** The dominant predictor of placement is structural disadvantage, not the gun-homicide count itself. The hardship index alone classifies placement with AUC 0.93 (LPM R² 0.51), versus AUC 0.79 (R² 0.15) for pre-period violence alone. Adding violence on top of hardship barely moves fit (R² 0.51 → 0.51), and the violence coefficient is not significant once hardship is held fixed (+1.6 pp, p = 0.55). Read carefully: violence and hardship are entangled (the same neighborhoods score high on both), so this does not say violence was irrelevant — it says the city's siting tracked the broad map of disinvestment at least as tightly as it tracked the specific gun-homicide burden. This refines, rather than overturns, the paper's selection-on-violence premise: treatment was assigned on characteristics (hardship, violence) that are themselves strong predictors of the outcome, which is what makes a clean counterfactual impossible.

**Caveat on the violence measure.** Pre-period violence enters as the gun-homicide *count* (mean/yr, 2009-2016), consistent with the rest of the paper and with the plausible premise that the city responded to the absolute burden of gunfire. A per-capita *rate* would be the natural robustness check but requires community-area population, which the 2008-2012 indicators file does not carry. Adding it (one ACS pull) is a listed extension.

## 2b. Equity: did placement track racial composition?

ShotSpotter areas are on average 49% Black and 32% Hispanic, against 14% and 17% in never-treated areas; the White non-Hispanic share is the single sharpest separator in the entire analysis (treated 13% vs control 56%, AUC 0.95).

| Composition | Treated mean | Control mean | Std. diff | Univariate AUC |
|---|---|---|---|---|
| Black share of population (%) | 48.9 | 13.5 | +1.03 | 0.733 |
| Hispanic/Latino share (%) | 32.0 | 17.4 | +0.55 | 0.537 |
| White, non-Hispanic share (%) | 12.6 | 55.7 | -2.80 | 0.948 |

Adding Black and Hispanic share to the placement model lifts logit AUC from 0.935 to 0.956, and the association survives conditioning: in an LPM on pre-period violence, hardship, income, and race, a one-SD increase in the Black share is associated with a +33.3 pp change in placement probability (p = 0.002) and the Hispanic share +28.2 pp (p = 0.003), while violence and hardship themselves turn insignificant. Read carefully: race, hardship, and violence are deeply entangled in a segregated city, and with 77 areas the model cannot cleanly separate them, so this is *not* evidence that race drove siting independently of disadvantage. What it shows, descriptively, is that ShotSpotter placement tracked Chicago's racial geography at least as tightly as it tracked violence or hardship---the concentrated-surveillance pattern advocacy groups described.

*Vintage caveat:* the community-area race figures are the 2019--2023 ACS 5-year (the only community-area race table the Chicago open-data portal exposes), which post-dates the 2017/2018 rollout. Chicago's community-area racial composition is highly persistent across decades, so we use it as a proxy for the pre-rollout composition; a 2008--2012-vintage race table (via the Census tract API) would remove the mismatch.

## 3. Propensity-score overlap

Propensity model (violence + hardship + income) AUC = 0.935. Median propensity is 0.97 for ShotSpotter areas versus 0.19 for never-treated areas. Separation is strong but not complete — the two groups mix in the middle of the scale, where a handful of lower-hardship areas received ShotSpotter and a handful of higher-hardship areas did not.

The decisive fact is at the top of the scale, and it needs no hand-picked subset. The highest propensity among all 26 never-treated areas is 0.873, and **36 of the 51 treated areas lie above every single never-treated unit**. At a threshold of 0.90, 35 of 51 treated areas qualify against 0 of 26 never-treated; at 0.95 it is 32 versus 0, and at 0.97 27 versus 0. Those 36 areas span 0.6 to 37.4 gun homicides per year, so the support failure is not confined to the violence tail: hardship drives it too. The placement model therefore states the identification problem from the selection side: the city installed ShotSpotter across essentially the whole high-hardship, high-violence region of the city, leaving no comparably situated untreated unit from which to build a counterfactual for those areas.

The six neighborhoods used as illustrative cases elsewhere in the paper (Austin, Humboldt Park, North Lawndale, Englewood, West Englewood, South Shore) sit inside that region, all with propensity at least 0.97, which 0 of 26 never-treated areas reach. They are carried over from the original analysis and are ranks 1, 2, 3, 4, 8 and 9 of 51 by pre-period violence, not the literal top six; nothing in the identification argument depends on which of the high-propensity areas are displayed.

## Figures

- `plots/placement_coefficients.png` — standardized differences (treated vs. never-treated), by characteristic
- `plots/placement_propensity_overlap.png` — propensity distributions, treated vs. control

## Data note

Race data: `data/raw/ACS_5yr_race_by_community_area_2023.csv` (Chicago Data Portal dataset t68z-cikk, ACS 2019--2023 5-year by community area). All other inputs as in the main pipeline.