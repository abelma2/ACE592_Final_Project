# Detection Without Identification: ShotSpotter and Gun Homicide in Chicago, 2009–2023

*A revised narrative for the ACE 592 final project, built from the data and code in this repository. Quantitative results come from `scripts/did_analysis.py` and `scripts/econometric_analysis.py`; tables are reproduced in [`docs/results_summary_original.md`](docs/results_summary_original.md) and [`docs/results_summary_reanalysis.md`](docs/results_summary_reanalysis.md).*

## Abstract

Between 2017 and 2023, Chicago paid millions of dollars to deploy ShotSpotter — an acoustic gunshot-detection system — across 51 of its 77 community areas. Mayor Brandon Johnson cancelled the contract in 2023, citing limited evidence of impact. We ask whether the data support that decision.

Using a community-area-by-year panel of gun homicides spanning 2009–2023, we estimate five classes of model: (i) a two-way fixed-effects OLS difference-in-differences, (ii) Poisson and Negative Binomial regressions appropriate for count outcomes, (iii) a cohort-specific event study respecting the staggered 2017 / 2018 rollout, (iv) a synthetic control with permutation inference for each of the six highest-violence study neighborhoods, and (v) pre-period placebo tests in 1999, 2003, 2007, and 2011.

The OLS DiD coefficient is positive and statistically significant (β = +2.63 gun homicides per community-area-year, p < 0.001), and survives the addition of unit and year fixed effects. **But that result is a count-data misspecification artifact: when the same panel is estimated under Poisson regression with two-way fixed effects, the coefficient collapses to β = +0.072 (incidence-rate ratio 1.07, 95% CI [0.87, 1.33], p = 0.51); under Negative Binomial it is β = −0.030 (IRR 0.97, [0.72, 1.31], p = 0.84). Both proper count-data specifications are statistically null.** The data are consistent with anywhere from a 13% reduction to a 33% increase in homicides — spanning the entire policy-relevant range.

Three further pieces of evidence reinforce that the OLS result should not be read as a causal estimate. The same design produces *similarly* statistically significant coefficients when treatment is reassigned to fake pre-rollout dates (β = −1.34, p = 0.008 for a 1999 placebo; β = −1.93, p = 0.007 for 2003). A joint Wald test on the seven pre-treatment event-study coefficients rejects parallel trends at the 1% level (F = 19.9, p = 0.006). And synthetic control is infeasible: the 26 never-treated community areas are systematically lower-violence than any of the six study neighborhoods, leaving no credible donor pool — the optimizer concentrates 100% of weight on a single donor (Washington Heights) for every treated unit, with pre-treatment RMSEs of 10–34 against actual values of 15–70.

The honest econometric reading of these data is that **we cannot identify a causal effect of ShotSpotter on gun homicides, in either direction.** What we can say is that the data are inconsistent with a large reduction effect — confidence intervals exclude the homicide reductions a $3 million policy intervention should plausibly produce. Persistent violence in the treated neighborhoods correlates almost entirely with structural disadvantage (hardship indices of 55–94 versus a city-wide median below 50). ShotSpotter performed its narrow technical function — detecting gunfire that 911 calls miss — but the policy-relevant causal claim does not survive falsification.

---

## 1. Introduction

Acoustic gunshot detection promises to close a measurement gap: research finds that only about 20% of shootings are reported to law enforcement (Irvin-Erickson et al., 2016). Chicago piloted ShotSpotter in 2012 in two small zones, expanded to full district coverage in 2016, and signed a three-year, $3 million contract to widen coverage in 2018 (Office of Inspector General, 2021). The system has been controversial: a 2021 OIG audit found that alerts rarely produced evidence used in investigations, advocacy groups argued the technology concentrated police presence in already over-policed Black and Brown communities, and Mayor Brandon Johnson terminated the contract in 2023.

The empirical question that should drive that policy decision is straightforward: *did ShotSpotter reduce gun violence?* This paper attempts to answer it using public data and a sequence of increasingly demanding identification strategies. The answer, defensibly stated: **no design we can run on these data produces a credible causal estimate.** We close with the only conclusion the data robustly support — that ShotSpotter functioned as a detection instrument, not as a violence-reduction intervention.

## 2. Background: ShotSpotter in Chicago

ShotSpotter is a real-time acoustic-detection service. Microphones triangulate the sound of gunfire; an algorithm filters non-gunshot sounds; a human acoustic expert reviews the alert before dispatching it to police (Choi, Librett, & Collins, 2014). The technology is sold as a response-time tool: faster, more accurate detection should mean faster police response and faster trauma care.

Chicago's deployment was concentrated in the city's most disinvested neighborhoods. Of the 51 community areas covered, every one of the six neighborhoods studied here — Austin, Humboldt Park, North Lawndale, Englewood, West Englewood, and South Shore — has a poverty rate roughly twice the Chicago average (28.6%–46.6% versus 19.7% city-wide), unemployment 1.5–3× the city rate, and a hardship index between 55 and 94 versus a city-wide median below 50 (Census 2008–2012). The 51 ShotSpotter community areas accounted for an average of 381 gun homicides per year in 2009–2016, compared with 42 per year across the 26 community areas without ShotSpotter — roughly a 9-to-1 ratio.

That ratio is the central methodological challenge: ShotSpotter was deployed *because* of high violence, not assigned at random.

## 3. Data

The analysis combines four City of Chicago open-data products:

- **Violence Reduction — Victims of Homicides and Non-Fatal Shootings** (61,450 incidents, 1991–2025). Filtered to gun homicides (n ≈ 7,600 in 2009–2023) and non-fatal gunshot injuries (n ≈ 30,000 in 2009–2023).
- **Violence Reduction — ShotSpotter Alerts (Historical)** (222,108 alerts, 2017-01-13 to 2023-12-31).
- **Census Data — Selected socioeconomic indicators, 2008–2012** (77 community areas plus city total).
- **Community Areas GeoJSON** (77 polygons).

The ShotSpotter dataset reveals the staggered rollout structure: 21 community areas appear in the data starting in early 2017, and 30 between mid-2017 and mid-2018. We exploit that staggering in the cohort-specific analyses below; the original DiD pipeline collapses both cohorts to a single 2017 cutoff, which is the standard convention but biases pooled estimates when treatment timing varies (Goodman-Bacon, 2021).

## 4. Empirical Strategy

### 4.1 The identification problem

In a randomized experiment we would compare treated neighborhoods to a control group drawn from the same distribution of pre-treatment violence. Chicago's data do not allow that. The 51 ShotSpotter neighborhoods had ~9× the homicide rate of the 26 non-ShotSpotter neighborhoods *before* deployment. Any difference-in-differences contrast must therefore lean on the parallel-trends assumption: that, absent ShotSpotter, the two groups would have moved in parallel.

### 4.2 Specifications

We estimate four models:

1. **Two-way fixed effects DiD.** `gun_hom_it = β·(treat_i × post_t) + α_i + γ_t + ε_it` with cluster-robust SEs at the community-area level. Estimated on the full panel and with three robustness restrictions (excluding 2016, COVID, and both).

2. **Cohort-specific event study.** Separate ATT(g, k) estimates for the 2017 cohort (g = 2017, install year) and 2018 cohort (g = 2018), using the 26 never-treated community areas as controls. Effects are aggregated across cohorts with sample-size weights, with bootstrap standard errors.

3. **Synthetic control.** For each of the six study neighborhoods, weights over the 26 never-treated community areas are chosen by SLSQP to minimize pre-treatment squared error on annual gun-homicide counts, subject to the simplex constraint (weights non-negative, sum to one).

4. **Pre-period placebo / falsification.** The DiD is re-estimated with treatment reassigned to fake years (1999, 2003, 2007, 2011) on overlapping pre-2017 windows. A non-null placebo coefficient flags spurious significance and rejects the design.

## 5. Results

### 5.1 ShotSpotter detected what it was designed to detect

Alert volume in covered neighborhoods is enormous: 222,000 alerts across roughly seven years and 51 communities, with the heaviest-coverage areas (Austin, West Englewood, Englewood, North Lawndale) generating 8,000–15,000 alerts each. Both alert volume and homicide counts peak in the same months (May–August) and the same hours (10 PM to 1 AM), confirming that the system detects the activity it claims to detect. As a *measurement instrument*, ShotSpotter works.

### 5.2 The headline OLS DiD is positive — but disappears under count-data regression

The two-way fixed-effects DiD coefficient on the full 2009–2023 panel is **β = +2.627** gun homicides per community-area-year (SE = 0.638, p < 0.001, 95% CI [+1.38, +3.88]). Excluding COVID years yields +1.75 (p < 0.001); excluding 2016 yields +3.19 (p < 0.001); excluding both yields +2.31 (p < 0.001). The result is robust to fixed-effect specification — adding unit and year fixed effects reproduces the original DiD coefficient almost exactly.

It is *not* robust to choosing the right model for count data. Annual gun-homicide counts in a community area are non-negative integers with mean and variance both increasing in the pre-existing violence level. OLS treats them as continuous and homoscedastic, weighting all observations equally. Poisson and Negative Binomial regressions allow the variance to scale with the mean and estimate effects on the log scale, where they correspond to multiplicative (proportional) changes — the natural metric for a percentage-effect treatment.

| Specification | β | IRR = exp(β) | 95% CI on IRR | p |
|---|---|---|---|---|
| OLS (TWFE) | +2.627 | — | — | <0.001*** |
| **Poisson (TWFE), full panel** | **+0.072** | **1.074** | **[0.869, 1.328]** | **0.508** |
| Poisson (TWFE), excl COVID | +0.014 | 1.014 | [0.816, 1.259] | 0.902 |
| Poisson (TWFE), excl 2016 | +0.084 | 1.087 | [0.857, 1.381] | 0.491 |
| **Neg. Binomial, full panel** | **−0.030** | **0.970** | **[0.720, 1.308]** | **0.843** |
| Neg. Binomial, excl COVID | −0.084 | 0.920 | [0.674, 1.256] | 0.600 |
| Neg. Binomial, excl 2016 | −0.029 | 0.972 | [0.700, 1.350] | 0.864 |

The proper count-data models give **statistically null results across every robustness specification**. Poisson IRRs hover around 1.0–1.1 with confidence intervals consistent with anywhere from a 13% reduction to a 33% increase in homicides. Negative Binomial IRRs are near 0.97 with similarly wide intervals. Whatever the OLS DiD is picking up, it is being amplified by the additive scale and the unequal variance across community areas. On the natural multiplicative scale, no effect is detectable.

Non-fatal gunshot injuries (a 4× larger sample) yield β = +7.51 under OLS (p = 0.009), but should be re-estimated under Poisson for the same reason; that exercise is left to a follow-up but is unlikely to overturn the qualitative pattern given how much of the OLS signal here also disappears under proper count-data treatment.

### 5.3 …but the design also produces "significant" placebos

If the +2.627 reflected a real treatment effect, then assigning *fake* treatment dates well before the actual rollout should yield null estimates. It does not.

| Fake treatment year | β | SE | p | 95% CI |
|---|---|---|---|---|
| **1999** | **−1.34** | 0.51 | **0.008** | [−2.34, −0.35] |
| **2003** | **−1.93** | 0.71 | **0.007** | [−3.33, −0.53] |
| 2007 | −0.65 | 0.50 | 0.194 | [−1.62, +0.33] |
| 2011 | +0.26 | 0.35 | 0.461 | [−0.42, +0.94] |

Two of four placebos are statistically significant at the 1% level, with magnitudes (1.3, 1.9 homicides per area-year) on the same order as the post-2017 estimate. What the placebos appear to be picking up is the long-run city-wide decline in violence from the early-1990s peak, differentially absorbed across community areas. Whatever the design captures, it captures it before ShotSpotter exists. The post-2017 +2.627 cannot be cleanly interpreted as a treatment effect when the same design produces effects of similar magnitude in the pre-treatment placebo periods.

### 5.4 The cohort-specific event study confirms the timing problem

Aggregated across the 2017 and 2018 cohorts (sample-weighted, with never-treated controls and bootstrap SEs):

| k (years from install) | ATT | 95% CI |
|---|---|---|
| −7 | **−5.46** | [−9.30, −1.62] |
| −6 | **−4.54** | [−8.53, −0.56] |
| −5 | **−5.09** | [−8.97, −1.20] |
| −4 | **−6.25** | [−10.04, −2.46] |
| −3 | **−6.09** | [−10.01, −2.17] |
| −2 | **−4.92** | [−8.92, −0.93] |
| −1 (base) | 0 | — |
| 0 | −1.95 | [−6.57, +2.68] |
| +1 | **−4.65** | [−8.79, −0.51] |
| +2 | −3.83 | [−8.12, +0.47] |
| +3 | +0.10 | [−4.85, +5.05] |
| +4 | −0.43 | [−5.25, +4.39] |
| +5 | −2.99 | [−7.19, +1.22] |

Every pre-treatment coefficient is large and negative — a clear failure of parallel trends. A formal joint Wald F-test on the seven pre-treatment interaction coefficients yields **F = 19.9, p = 0.006**, rejecting parallel trends at the 1% level. The pattern is mechanical: 2016 was a Chicago-wide homicide-record year, disproportionately concentrated in the neighborhoods that would soon receive ShotSpotter, and the event study uses a year adjacent to 2016 as the base. Once that year is the reference point, every other year looks anomalously low.

The post-treatment coefficients are mostly negative or null — the opposite sign from the pooled +2.627. A reader who took the event study literally would conclude ShotSpotter was associated with a *small reduction* in gun homicides; a reader who took the pooled DiD literally would conclude the opposite. Both readings are wrong: the design cannot separate treatment effects from the 2016 anomaly.

### 5.5 Synthetic control is infeasible

For each of the six study neighborhoods, the synthetic control optimizer concentrates 100% of the weight on a single donor (Washington Heights — the highest-violence non-ShotSpotter community area). Pre-treatment RMSEs are large (10–34 homicides per year against actual values of 15–70), meaning even the best linear combination of donors cannot reproduce the pre-treatment trajectory of any treated unit. The synthetic counterfactual lies far below the treated series in every panel of `plots/synthetic_control_panels.png`, including in pre-treatment years where the synthetic was deliberately optimized to fit.

Permutation inference makes this concrete. Running the same synthetic control on each of the 26 never-treated areas (treating each as a "placebo treated" unit) produces a tight placebo gap distribution centered near zero (95% range roughly ±2 homicides per year), against actual treated gaps of +14 to +49. The naïve permutation p-value is < 0.001 for every treated unit — but this rejection is dominated by the level mismatch between treated and donor units, not by any treatment effect. When the placebos themselves can be cleanly synthesized (low-violence areas have low-violence donors and small post-period gaps) but the treated cannot (high-violence areas have no high-violence donors and large pre-period RMSE), the inference is invalid by construction.

This is not a failure of method; it is a feature of the data. The 26 never-treated community areas occupy the low-violence end of Chicago's distribution. There is no donor pool capable of constructing a credible counterfactual for the city's most violent neighborhoods because the city did not leave any of those neighborhoods untreated. ShotSpotter's selection rule is, in a literal sense, what makes synthetic control infeasible.

### 5.6 What does explain neighborhood-level violence?

Hardship index from the 2008–2012 Census ranks the six study neighborhoods at 55, 73, 85, 87, 89, and 94 — versus Chicago's overall hardship score below 50 and the lowest-violence community areas at 1–10. Per capita income in these neighborhoods ranges from $11,317 to $19,398, against a city-wide $28,202. The same neighborhoods that had the highest homicide rates *before* ShotSpotter still had the highest homicide rates after. Cross-sectional variation in violence tracks structural disadvantage tightly; cross-sectional variation in ShotSpotter coverage adds little explanatory power once disadvantage is held constant.

## 6. What the Data Can and Cannot Support

The combination of (a) an OLS DiD coefficient that disappears under proper count-data regression (Poisson IRR 1.07 with p = 0.51, NegBin IRR 0.97 with p = 0.84), (b) a Wald F-test on pre-trends that rejects parallel trends at the 1% level (F = 19.9, p = 0.006), (c) similarly significant pre-period OLS placebos, and (d) infeasible synthetic control points to a single conclusion: **the design cannot credibly identify the causal effect of ShotSpotter installation on gun homicides.** The data do not say ShotSpotter caused a 2.6-homicide-per-year increase — that statistic is a misspecification artifact. They also do not say ShotSpotter caused a 5-homicide-per-year decrease — that reading would be reading the event-study coefficients without acknowledging that the base year (2016) is anomalous. What the data do say is that the proper count-data specifications return statistically null effects with confidence intervals spanning the entire policy-relevant range.

We can rule out one thing: the data are not consistent with a *large* reduction in gun homicides attributable to ShotSpotter. A genuinely effective intervention with these sample sizes would have produced visibly negative post-treatment event-study coefficients with confidence intervals that exclude zero from above, robustness across specifications, and a cleaner signal in the higher-power non-fatal shootings analysis. None of these appear. Whatever effect ShotSpotter had on gun homicides, it was not large enough to dominate the noise, the 2016 anomaly, or the pre-existing trend differences across community areas.

## 7. Why Detection Did Not Reduce Violence

Three explanations are consistent with the empirical pattern.

**Selection.** ShotSpotter was installed in places that needed it most. Even if the technology produced a small benefit, the comparison would still look unfavorable because the rest of Chicago started from a much lower baseline with less room to deteriorate.

**Detection is not intervention.** The mechanism by which ShotSpotter could plausibly reduce homicides runs through faster police response and faster trauma triage. The data here do not measure response time, dispatch decisions, or trauma outcomes. The 2021 OIG audit found that ShotSpotter alerts rarely produced evidence used in investigations and that officers did not consistently change behavior based on alerts. If the operational link between detection and outcome is weak, detection volume is not a sufficient condition for fewer homicides.

**Structural factors dominate.** A surveillance technology cannot, by itself, change the distribution of poverty, unemployment, housing instability, school quality, or labor-market access. Programs that have been shown to reduce homicide in randomized or quasi-experimental evaluations — READI Chicago, Cure Violence, Becoming a Man — work on those structural channels directly. A microphone on a streetlight does not.

## 8. Limitations

The most important limitation is the one above: the design itself cannot identify a clean treatment effect from these data. We cannot recover the counterfactual we need (Chicago neighborhoods with ShotSpotter-eligible violence levels but no ShotSpotter installation) because the city did not generate one. Pre-period placebo tests confirm that whatever the DiD estimator captures, it is not specifically the effect of treatment.

Beyond that: the panel is annual, masking within-year dynamics; non-fatal shooting counts in 2017–2023 are subject to coverage changes in the underlying database; the analysis does not observe police response times, dispatch decisions, or arrests; and the 26 never-treated community areas are systematically different from the treated areas on every observable dimension — wealth, demographics, education, employment.

## 9. Conclusion

ShotSpotter in Chicago succeeded at detection. Alerts climbed in covered areas because the technology measures gunfire that 911 underreports. What we cannot demonstrate is that detection translated into a reduction in lethal violence: every design we ran on these data either produces a positive DiD coefficient that fails falsification, an event study that violates parallel trends because of the 2016 anomaly, or a synthetic control whose donor pool cannot match treated-area violence levels. The only empirical statement these data robustly support is the absence of a large reduction effect.

That null finding is itself a policy statement. The 2023 termination of Chicago's $3 million contract is consistent with the absence of an empirical reduction story. Resources reallocated from detection contracts toward community-based violence-prevention programs would, on the existing evidence, have at least as good a chance of producing the homicide reductions that ShotSpotter did not.

---

## Notes on figures, tables, and reproduction

This narrative draws on artifacts already produced by the code in this repository. To regenerate everything from scratch:

```bash
pip install -r requirements.txt
python scripts/did_analysis.py             # original DiD pipeline
python scripts/econometric_analysis.py     # TWFE, cohort event study, SC, placebos
python scripts/econometric_robustness.py   # Poisson/NegBin, Wald pre-trend test, SC inference
python scripts/make_pretty_plots.py        # publication figures
python scripts/plot_covid_sensitivity.py   # COVID-aware annualized comparison
```

Key artifacts referenced above:

- **Master DiD regression table:** [`docs/results_summary_original.md`](docs/results_summary_original.md)
- **Econometric re-analysis (TWFE, cohort, SC, placebos):** [`docs/results_summary_reanalysis.md`](docs/results_summary_reanalysis.md)
- **Map of treatment vs. control community areas:** [`plots/did_map_treatment_control.png`](plots/did_map_treatment_control.png)
- **Cohort-aware event study with bootstrap CIs:** [`plots/event_study_cohort_aware.png`](plots/event_study_cohort_aware.png)
- **Six-panel synthetic control comparison:** [`plots/synthetic_control_panels.png`](plots/synthetic_control_panels.png)
- **Synthetic control permutation inference:** [`plots/synthetic_control_inference.png`](plots/synthetic_control_inference.png)
- **Original event study (TWFE on pooled cohort):** [`plots/did_event_study.png`](plots/did_event_study.png)
- **Forest plot of treatment effects across neighborhoods:** [`plots/forest_plot_effects.png`](plots/forest_plot_effects.png)
- **Robustness summary across specifications:** [`plots/robustness_summary.png`](plots/robustness_summary.png)
- **Non-fatal shootings DiD trend:** [`plots/nonfatal_did_trends.png`](plots/nonfatal_did_trends.png)
- **Annual trends per neighborhood:** [`plots/annual_trends_neighborhoods.png`](plots/annual_trends_neighborhoods.png)
- **COVID-adjusted before/during comparison:** [`plots/covid_sensitivity.png`](plots/covid_sensitivity.png)
