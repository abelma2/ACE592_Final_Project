# Detection Without Identification: ShotSpotter and Gun Homicide in Chicago, 2009–2023

*A revised narrative of this study, built from the data and code in this repository. Quantitative results come from `scripts/did_analysis.py`, `scripts/econometric_analysis.py`, and `scripts/placement_model.py`; tables are reproduced in [`docs/results_summary_original.md`](docs/results_summary_original.md), [`docs/results_summary_reanalysis.md`](docs/results_summary_reanalysis.md), and [`docs/results_summary_placement.md`](docs/results_summary_placement.md).*

## Abstract

Between 2017 and 2023, Chicago paid millions of dollars to deploy ShotSpotter (an acoustic gunshot-detection system) across 51 of its 77 community areas. Mayor Brandon Johnson cancelled the contract in 2023, citing limited evidence of impact. That gunshot detection *doesn't* reduce gun violence is, by now, well established, for Chicago at the police-district level (Connealy et al., 2024) and nationally (Doucette et al., 2021). The question we take up is not whether that null replicates but whether the public, community-area-level data most local policy argument relies on can *identify* such an effect at all, and what a routine evaluation on them looks like when it cannot.

Using a community-area-by-year panel of gun homicides spanning 2009–2023, we estimate five classes of model: (i) a two-way fixed-effects OLS difference-in-differences, (ii) Poisson and Negative Binomial regressions appropriate for count outcomes, (iii) cohort-specific and Callaway–Sant'Anna (2021) event studies respecting the staggered 2017 / 2018 rollout, (iv) a synthetic control with permutation inference for each of the six highest-violence study neighborhoods, and (v) pre-period placebo tests in 1999, 2003, 2007, and 2011.

The OLS DiD coefficient is positive and statistically significant (β = +2.63 gun homicides per community-area-year, p < 0.001), and survives the addition of unit and year fixed effects. **But that result is a count-data misspecification artifact: when the same panel is estimated under Poisson regression with two-way fixed effects, the coefficient collapses to β = +0.072 (incidence-rate ratio 1.07, 95% CI [0.87, 1.33], p = 0.51); under Negative Binomial (with the dispersion estimated, not fixed) it is β = +0.067 (IRR 1.07, [0.86, 1.33], p = 0.55). Both proper count-data specifications are statistically null.** The data are consistent with anywhere from a 13% reduction to a 33% increase in homicides, spanning the entire policy-relevant range. A heterogeneity-robust Poisson extended two-way fixed-effects estimator (Wooldridge, 2023) returns the same null (IRR 1.05, 95% CI [0.79, 1.33]), so the result is not a fixed-effects weighting artifact. The modern Callaway–Sant'Anna estimator returns a *negative* overall ATT (−2.45 per area-year) on the full panel, but that sign rides on a single year: dropping the 2016 homicide-record year flips it back to +1.90. The instability is not estimator disagreement about a well-identified parameter; it is the symptom of a design with no valid control group and a pre-period contaminated by one outlier year, which cannot identify a causal effect in either direction.

Three further pieces of evidence reinforce that the OLS result should not be read as a causal estimate. The same design produces *similarly* statistically significant coefficients when treatment is reassigned to fake pre-rollout dates (β = −1.34, p = 0.008 for a 1999 placebo; β = −1.93, p = 0.007 for 2003). A joint Wald test on the seven pre-treatment event-study coefficients rejects parallel trends at the 1% level (F = 19.9, p = 0.006). And synthetic control is infeasible: the 26 never-treated community areas are systematically lower-violence than any of the six study neighborhoods, leaving no credible donor pool: the optimizer concentrates 100% of weight on a single donor (Washington Heights) for every treated unit, with pre-treatment RMSEs of 10–34 against actual values of 15–70.

A placement model makes the selection explicit: ShotSpotter siting is predicted by structural disadvantage and pre-period violence with near-perfect separation (propensity AUC 0.94), and (even conditional on those) by the Black and Hispanic share of a neighborhood, so tightly that the six highest-violence study neighborhoods have no comparably situated untreated counterpart. The honest econometric reading is that **we cannot identify a causal effect of ShotSpotter on gun homicides, in either direction.** The same non-identification recurs on an independent event: Chicago's September 2024 shutoff, where gun homicides were already falling city-wide in coverage and comparison areas alike (−37% versus −31% over the same months), so the widely cited post-removal decline reflects a city-wide trend, not a removal effect. The data rule out a large *deterrence* effect, but cannot exclude the survival-channel mortality benefit a boundary regression discontinuity attributes to faster emergency response, a magnitude that lies inside our own confidence interval.

The generalizable point, and the paper's contribution, is the failure class rather than the program. We call it **saturation targeting**: assignment is predicted by the same structural variables that predict the outcome, *and* coverage among units above the operative risk threshold approaches one. Both are ordinary features of competent, need-based administration, and together they destroy the counterfactual before any estimator runs. The failure leaves three signatures that are testable on public data before an effect is estimated (no overlap; an additive estimator reporting the baseline gap; instability that a Goodman–Bacon decomposition shows is *not* a staggered-timing problem), and the modern difference-in-differences toolkit repairs none of them, because those estimators correct weighting under heterogeneous timing rather than missing overlap. Matching makes the trade explicit: trimming to common support does restore parallel trends, but only by discarding all six study neighborhoods, so a credible comparison and a policy-relevant one cannot be had at once. The same arc runs through the COPS literature (Evans & Owens 2007 → Worrall & Kovandzic 2007 → Mello 2019): when treatment is awarded because of the outcome, the binding constraint is the assignment mechanism, and progress comes from finding variation the mechanism did not determine, not from refining the estimator.

---

## 1. Introduction

Acoustic gunshot detection promises to close a measurement gap: research finds that only about 20% of shootings are reported to law enforcement (Irvin-Erickson et al., 2016). Chicago piloted ShotSpotter in 2012 in two small zones, expanded to full district coverage in 2016, and signed a three-year, $3 million contract to widen coverage in 2018 (Office of Inspector General, 2021). The system has been controversial: a 2021 OIG audit found that alerts rarely produced evidence used in investigations, advocacy groups argued the technology concentrated police presence in already over-policed Black and Brown communities, and Mayor Brandon Johnson terminated the contract in 2023.

For Chicago, the empirical question (*did ShotSpotter reduce gun violence?*) is largely answered: district-level (Connealy et al., 2024) and national (Doucette et al., 2021) evaluations find no reduction in shootings, homicides, or clearance, and a district-level study finds detection even slowed dispatch and lowered arrest probability (Topper & Ferrazares, 2024). This narrative asks a different, methodological question: can the community-area data most local policy argument relies on *identify* such an effect at all, and what happens when a standard evaluation is run on them without asking? The answer, defensibly stated: **no design we can run on these data produces a credible causal estimate**, and the way a routine evaluation fails is itself the finding. We close with the only conclusion the data robustly support: ShotSpotter functioned as a detection instrument, not as a demonstrated violence-reduction intervention.

## 2. Background: ShotSpotter in Chicago

ShotSpotter is a real-time acoustic-detection service. Microphones triangulate the sound of gunfire; an algorithm filters non-gunshot sounds; a human acoustic expert reviews the alert before dispatching it to police (Choi, Librett, & Collins, 2014). The technology is sold as a response-time tool: faster, more accurate detection should mean faster police response and faster trauma care.

Chicago's deployment was concentrated in the city's most disinvested neighborhoods. Of the 51 community areas covered, every one of the six neighborhoods studied here (Austin, Humboldt Park, North Lawndale, Englewood, West Englewood, and South Shore) has a poverty rate roughly twice the Chicago average (28.6%–46.6% versus 19.7% city-wide), unemployment 1.5–3× the city rate, and a hardship index between 55 and 94 versus a city-wide median below 50 (Census 2008–2012). The 51 ShotSpotter community areas accounted for an average of 381 gun homicides per year in 2009–2016, compared with 42 per year across the 26 community areas without ShotSpotter, roughly a 9-to-1 ratio.

That ratio is the central methodological challenge: ShotSpotter was deployed *because* of concentrated violence and disadvantage, not assigned at random, a selection rule we estimate directly in §5.6.

The evaluation record is, by now, fairly settled. No rigorous study finds gunshot detection reduces shootings or homicides, in Chicago (Connealy et al., 2024) or nationally (Doucette et al., 2021), and a district-level study finds it slowed dispatch and lowered arrest probability (Topper & Ferrazares, 2024, forthcoming *Journal of Human Resources*); a 2026 meta-analysis of 44 estimates across eight studies pools to a null relative incidence-rate ratio of 1.02 (Huff, Dunlap, & Pearson, 2026). The one channel with a credible claim to benefit is medical: a boundary regression discontinuity attributes a lower shooting-fatality rate inside coverage to faster emergency response (University of Chicago Crime Lab, 2024, not peer-reviewed). Against that backdrop, this project's contribution is not the null itself but a demonstration of *why* the community-area data cannot identify it, plus a model of the selection that makes identification impossible.

## 3. Data

The analysis combines four City of Chicago open-data products:

- **Violence Reduction: Victims of Homicides and Non-Fatal Shootings** (61,450 incidents, 1991–2025). Filtered to gun homicides (n ≈ 7,600 in 2009–2023) and non-fatal gunshot injuries (n ≈ 30,000 in 2009–2023).
- **Violence Reduction: ShotSpotter Alerts (Historical)** (222,108 alerts, 2017-01-13 to 2024-09-22, when the city ended ShotSpotter; the 191,824 through 2023 fall in the analysis window). [Chicago Data Portal 3h7q-7mdb]
- **Census Data: Selected socioeconomic indicators, 2008–2012** (77 community areas plus city total).
- **Community Areas GeoJSON** (77 polygons).

The ShotSpotter dataset reveals the staggered rollout structure: 21 community areas appear in the data starting in early 2017, and 30 between mid-2017 and mid-2018. We exploit that staggering in the cohort-specific analyses below; the original DiD pipeline collapses both cohorts to a single 2017 cutoff, which is the standard convention but biases pooled estimates when treatment timing varies (Goodman-Bacon, 2021).

## 4. Empirical Strategy

### 4.1 The identification problem

In a randomized experiment we would compare treated neighborhoods to a control group drawn from the same distribution of pre-treatment violence. Chicago's data do not allow that. The 51 ShotSpotter neighborhoods had ~9× the homicide rate of the 26 non-ShotSpotter neighborhoods *before* deployment. Any difference-in-differences contrast must therefore lean on the parallel-trends assumption: that, absent ShotSpotter, the two groups would have moved in parallel.

### 4.2 Specifications

We estimate five classes of model:

1. **Two-way fixed effects DiD.** `gun_hom_it = β·(treat_i × post_t) + α_i + γ_t + ε_it` with cluster-robust SEs at the community-area level. Estimated on the full panel and with three robustness restrictions (excluding 2016, COVID, and both).

2. **Count-data regression.** Poisson and Negative Binomial DiD with the same two-way fixed effects, appropriate for the integer homicide counts; the NB dispersion is estimated (Cameron–Trivedi NB2 auxiliary regression) rather than fixed.

3. **Event studies.** A cohort-specific ATT(g, k) hand-rolled estimator (2017 and 2018 cohorts vs. the 26 never-treated areas, sample-weighted, bootstrap SEs) and the modern **Callaway & Sant'Anna (2021)** group-time ATT (via the `differences` package), which is the authoritative staggered-rollout estimate.

4. **Synthetic control.** For each of the six study neighborhoods, weights over the 26 never-treated community areas are chosen by SLSQP to minimize pre-treatment squared error on annual gun-homicide counts, subject to the simplex constraint (weights non-negative, sum to one).

5. **Pre-period placebo / falsification.** The DiD is re-estimated with treatment reassigned to fake years (1999, 2003, 2007, 2011) on overlapping pre-2017 windows. A non-null placebo coefficient flags spurious significance and rejects the design.

Separately from these outcome models, we estimate a cross-sectional model of *treatment assignment* itself (a linear-probability and logistic model of ShotSpotter placement across all 77 community areas on pre-period (2009–2016) gun-homicide levels, the 2008–2012 Census indicators, and (2019–2023 ACS) racial composition) to characterize the selection mechanism and the resulting lack of common support (§5.6).

## 5. Results

### 5.1 ShotSpotter detected what it was designed to detect

Alert volume in covered neighborhoods is enormous: 222,000 alerts across roughly seven years and 51 communities, with the heaviest-coverage areas (Austin, West Englewood, Englewood, North Lawndale) generating 8,000–15,000 alerts each. Both alert volume and homicide counts peak in the same months (May–August) and the same hours (10 PM to 1 AM), an expected pattern, since alerts and homicides both track the diurnal cycle of gunfire, and read here only as validation that the system detects the activity it claims to detect. As a *measurement instrument*, ShotSpotter works.

### 5.2 The headline OLS DiD is positive, but disappears under count-data regression

The two-way fixed-effects DiD coefficient on the full 2009–2023 panel is **β = +2.627** gun homicides per community-area-year (SE = 0.638, p < 0.001, 95% CI [+1.38, +3.88]). Excluding COVID years yields +1.75 (p < 0.001); excluding 2016 yields +3.19 (p < 0.001); excluding both yields +2.31 (p < 0.001). The result is robust to fixed-effect specification: adding unit and year fixed effects reproduces the original DiD coefficient almost exactly.

It is *not* robust to choosing the right model for count data. Annual gun-homicide counts in a community area are non-negative integers with mean and variance both increasing in the pre-existing violence level. OLS treats them as continuous and homoscedastic, weighting all observations equally. Poisson and Negative Binomial regressions allow the variance to scale with the mean and estimate effects on the log scale, where they correspond to multiplicative (proportional) changes, the natural metric for a percentage-effect treatment.

| Specification | β | IRR = exp(β) | 95% CI on IRR | p |
|---|---|---|---|---|
| OLS (TWFE) | +2.627 | n/a | n/a | <0.001*** |
| **Poisson (TWFE), full panel** | **+0.072** | **1.074** | **[0.869, 1.328]** | **0.508** |
| Poisson (TWFE), excl COVID | +0.014 | 1.014 | [0.816, 1.259] | 0.902 |
| Poisson (TWFE), excl 2016 | +0.084 | 1.087 | [0.857, 1.381] | 0.491 |
| **Neg. Binomial, full panel** | **+0.067** | **1.069** | **[0.860, 1.329]** | **0.548** |
| Neg. Binomial, excl COVID | +0.014 | 1.014 | [0.810, 1.269] | 0.902 |
| Neg. Binomial, excl 2016 | +0.074 | 1.077 | [0.844, 1.373] | 0.551 |

The proper count-data models give **statistically null results across every robustness specification**. Poisson IRRs hover around 1.0–1.1 with confidence intervals consistent with anywhere from a 13% reduction to a 33% increase in homicides. The Negative Binomial dispersion, estimated rather than fixed (α ≈ 0.02 via the Cameron–Trivedi NB2 auxiliary regression), is small (the community-area-by-year counts are only mildly overdispersed) so the Negative Binomial IRRs land essentially on top of the Poisson estimates, near 1.07 with similarly wide intervals. They also sit on top of the meta-analytic pooled effect for gunshot detection (relative IRR 1.02, 95% CI [0.90, 1.16]; Huff, Dunlap, & Pearson, 2026), estimated with the same incidence-rate-ratio metric, evidence the multiplicative-scale null is a property of the intervention, not of this panel. Whatever the OLS DiD is picking up, it is being amplified by the additive scale and the unequal variance across community areas. On the natural multiplicative scale, no effect is detectable.

Non-fatal gunshot injuries (a 4× larger sample) yield β = +7.51 under OLS (p = 0.009), and the count-data correction is even starker there: estimated as Poisson or Negative Binomial with the same fixed effects, the point estimate flips *below* one and is statistically null in every specification (IRRs 0.84–0.92, p ≥ 0.23). The higher-powered outcome delivers the same verdict: the additive OLS coefficient is a scale artifact.

The pattern is not even specific to gun-violence outcomes. Applying the identical DiD to *enforcement* output (weapons-violation incidents and arrests from the Crimes dataset, the public analog of the ballistic-evidence and gun-recovery "process" measures that Piza et al. (2024) found rose under gunshot detection in Kansas City) reproduces it exactly: weapons incidents show an OLS DiD of +51.8 per area-year (p < 0.001) that collapses to a null incidence-rate ratio of 1.13 (p = 0.48); arrests likewise (OLS +22.0; IRR 0.95). Two lessons: the additive-scale artifact is *general* (it recurs, unchanged, on a completely different outcome), and (unlike localized Kansas City) Chicago's community-area data detect no enforcement footprint of coverage on the correct scale, consistent with the 2021 OIG finding that alerts rarely produced investigative evidence. Whatever enforcement coverage concentrated, it fell on the lowest-life-expectancy, most heavily policed neighborhoods: the vendor's "procedural benefit" and the advocates' surveillance burden are the same coin.

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
| −1 (base) | 0 | n/a |
| 0 | −1.95 | [−6.57, +2.68] |
| +1 | **−4.65** | [−8.79, −0.51] |
| +2 | −3.83 | [−8.12, +0.47] |
| +3 | +0.10 | [−4.85, +5.05] |
| +4 | −0.43 | [−5.25, +4.39] |
| +5 | −2.99 | [−7.19, +1.22] |

Every pre-treatment coefficient is large and negative, a clear failure of parallel trends. A formal joint Wald F-test on the seven pre-treatment interaction coefficients yields **F = 19.9, p = 0.006**, rejecting parallel trends at the 1% level. The pattern is mechanical: 2016 was a Chicago-wide homicide-record year, disproportionately concentrated in the neighborhoods that would soon receive ShotSpotter, and the event study uses a year adjacent to 2016 as the base. Once that year is the reference point, every other year looks anomalously low.

The post-treatment coefficients are mostly negative or null, the opposite sign from the pooled +2.627. A reader who took the event study literally would conclude ShotSpotter was associated with a *small reduction* in gun homicides; a reader who took the pooled DiD literally would conclude the opposite. Both readings are wrong: the design cannot separate treatment effects from the 2016 anomaly.

The hand-rolled cohort event study above is intuitive but not the modern estimator the staggered-rollout literature prescribes. As a final check we re-estimate the effect with the **Callaway & Sant'Anna (2021)** group-time ATT (the `differences` package, never-treated controls, influence-function inference). On the full panel the overall ATT is **−2.45 per community-area-year (SE 0.89)**, negative, the opposite sign of the pooled OLS, and stable to using not-yet-treated rather than never-treated controls (−2.51). But this negative sign is *not* a property of the estimator. It is driven almost entirely by the 2016 record-homicide year, which sits at the *k* = −1 baseline for the 2017 cohort: the event study shows only the *k* = −1 (2016) pre-coefficient is significant while the rest hover at zero, and **excluding 2016 and re-estimating returns +1.90 (SE 0.72), the same sign as the OLS DiD**. The apparent OLS-vs-Callaway–Sant'Anna sign reversal is therefore a single-year sensitivity, not a sign that the estimators fundamentally disagree about a well-identified parameter. What it confirms is the opposite of identification: with no clean counterfactual and a pre-period dominated by one anomalous year, the data do not pin down even the sign of the effect (`plots/callaway_santanna_event_study.png`). A Goodman–Bacon decomposition completes the diagnosis: the "forbidden" already-treated comparison carries only 4.6% of the staggered TWFE weight (and the weighted 2×2 components reproduce the aggregate exactly), so the instability is the 2016 anomaly and the missing control group, not staggered-timing contamination (`plots/bacon_decomposition.png`). One loose end remains: the headline null is itself a two-way-FE Poisson, which under staggered timing inherits the same weighting concern as linear TWFE, and Callaway–Sant'Anna is linear and does not transfer to the multiplicative outcome the count argument calls for. So we estimate the count analog directly: a Wooldridge (2023) extended two-way fixed-effects Poisson that saturates the treatment in cohort and calendar time (`scripts/etwfe_count.py`). It returns the same null: an overall ATT incidence-rate ratio of **1.05 (95% CI [0.79, 1.33]**, community-area cluster bootstrap), with cohort-time effects at one at every horizon (`plots/etwfe_event_study.png`). The heterogeneity-robust count estimator removes the last escape hatch: the count-model null is not an artifact of the plain PPML's weighting.

### 5.5 Synthetic control is infeasible

For each of the six study neighborhoods, the synthetic control optimizer concentrates 100% of the weight on a single donor (Washington Heights, the highest-violence non-ShotSpotter community area). Pre-treatment RMSEs are large (10–34 homicides per year against actual values of 15–70), meaning even the best linear combination of donors cannot reproduce the pre-treatment trajectory of any treated unit. The synthetic counterfactual lies far below the treated series in every panel of `plots/synthetic_control_panels.png`, including in pre-treatment years where the synthetic was deliberately optimized to fit.

Permutation inference makes this concrete. Running the same synthetic control on each of the 26 never-treated areas (treating each as a "placebo treated" unit) produces a tight placebo gap distribution centered near zero (95% range roughly ±2 homicides per year), against actual treated gaps of +14 to +49. The naïve permutation p-value is < 0.001 for every treated unit, but this rejection is dominated by the level mismatch between treated and donor units, not by any treatment effect. When the placebos themselves can be cleanly synthesized (low-violence areas have low-violence donors and small post-period gaps) but the treated cannot (high-violence areas have no high-violence donors and large pre-period RMSE), the inference is invalid by construction.

This is not a failure of method; it is a feature of the data. The 26 never-treated community areas occupy the low-violence end of Chicago's distribution. There is no donor pool capable of constructing a credible counterfactual for the city's most violent neighborhoods because the city did not leave any of those neighborhoods untreated. ShotSpotter's selection rule is, in a literal sense, what makes synthetic control infeasible.

### 5.6 What explains neighborhood violence, and ShotSpotter placement

Hardship index from the 2008–2012 Census ranks the six study neighborhoods at 55, 73, 85, 87, 89, and 94, versus Chicago's overall hardship score below 50 and the lowest-violence community areas at 1–10. Per capita income in these neighborhoods ranges from $11,317 to $19,398, against a city-wide $28,202. The same neighborhoods that had the highest homicide rates *before* ShotSpotter still had the highest homicide rates after. Cross-sectional variation in violence tracks structural disadvantage tightly; cross-sectional variation in ShotSpotter coverage adds little explanatory power once disadvantage is held constant.

**Placement was nearly deterministic.** We model ShotSpotter placement across all 77 community areas on pre-period gun-homicide levels and the 2008–2012 Census indicators ([`plots/placement_coefficients.png`](plots/placement_coefficients.png)). Structural disadvantage predicts siting even more tightly than the homicide count itself: the hardship index separates treated from never-treated areas with an AUC of 0.93 and per-capita income with 0.92, against 0.79 for the pre-period gun-homicide count, and once hardship is held fixed the violence count adds essentially nothing. The city sited ShotSpotter along the map of disinvestment, of which concentrated gun violence is one facet.

The consequence for identification is the selection-side statement of §5.5. A propensity model (violence, hardship, income) separates the two groups with an AUC of 0.94; the median propensity is 0.97 for ShotSpotter areas against 0.19 for never-treated areas; and the six highest-violence study neighborhoods all carry a placement propensity of at least 0.97, a band *no* never-treated area reaches ([`plots/placement_propensity_overlap.png`](plots/placement_propensity_overlap.png)). The city installed ShotSpotter in essentially every highest-violence, highest-hardship neighborhood, leaving no comparable untreated unit from which to build a counterfactual, which is exactly why the synthetic control in §5.5 is infeasible.

**Placement tracked race.** ShotSpotter areas are on average 49% Black and 32% Hispanic, against 14% and 17% in never-treated areas, and the white non-Hispanic share is the single sharpest separator in the whole analysis (13% versus 56%, AUC 0.95). The pattern survives conditioning: adding Black and Hispanic share to the placement model raises its AUC from 0.94 to 0.96, and in a linear-probability model with pre-period violence, hardship, and income, a one-SD rise in the Black share is associated with a +33 percentage-point change in placement probability (p = 0.002) and the Hispanic share +28 points (p = 0.003), while violence and hardship themselves turn insignificant. Race, hardship, and violence are deeply entangled in a segregated city, and with 77 areas the design cannot separate them cleanly; this is *not* evidence that race drove siting independently of disadvantage. What it shows, descriptively, is that ShotSpotter placement tracked Chicago's racial geography at least as tightly as it tracked violence or hardship. That coverage blankets the city's highest-Black-and-Latino police districts is already documented (MacArthur Justice Center, 2021; Office of Inspector General, 2021); our marginal contribution is the *conditional* result: that racial composition predicts placement even once pre-period violence and hardship are held fixed. (Community-area race is measured from the 2019–2023 ACS, which post-dates the rollout; because Chicago's neighborhood racial composition is highly persistent, we read it as a proxy for the pre-rollout map.)

### 5.7 A second natural experiment: the September 2024 removal

Chicago switched ShotSpotter off on 22 September 2024, and the gun-homicide decline that followed has been read in public commentary as proof the technology never worked. The shutoff is a genuine second experiment, and applying the paper's discipline to it shows the removal effect is no more identifiable than the installation effect, and fails in precisely the same two ways (`scripts/removal_analysis.py`, `plots/removal_trajectory.png`).

**The common trend.** Gun homicides had been falling city-wide from their 2021–2022 peak for two full years before the shutoff, in ShotSpotter and never-treated areas alike. On a common indexed scale the two groups decline in near-lockstep straight through the removal date: coverage-area homicides fell 37% from the pre- to the post-removal window, but never-treated areas fell 31% over the very same months. The naive "crime fell after ShotSpotter left" comparison attributes a city-wide trend to the shutoff: the mirror image of a naive "crime rose after ShotSpotter arrived" comparison.

**The scale artifact.** Netting out the common trend with a difference-in-differences does not rescue identification, because the additive DiD reproduces the artifact of §5.2: it returns a "significant" removal *benefit* of −0.25 gun homicides per area-month (p < 0.001), but only because coverage areas start from a far higher baseline (0.85 versus 0.20 per area-month) and so move more in absolute terms for the same proportional change. On the correct multiplicative scale the removal effect is null (Poisson IRR 0.91, 95% CI [0.54, 1.52]). The post window is also short: seven complete months, with 2025 only partially observed, so the estimate is under-powered on top of being non-identified. Both the raw drop and the additive DiD are mirages: these data can no more identify the removal effect than the installation effect, in either direction, and the public debate has confidently drawn both unwarranted inferences.

## 6. What the Data Can and Cannot Support

The combination of (a) an OLS DiD coefficient that disappears under proper count-data regression (Poisson IRR 1.07 with p = 0.51, NegBin IRR 1.07 with p = 0.55), (b) a Wald F-test on pre-trends that rejects parallel trends at the 1% level (F = 19.9, p = 0.006), (c) similarly significant pre-period OLS placebos, (d) infeasible synthetic control, and (e) a modern Callaway–Sant'Anna estimate whose sign reverses (−2.45 to +1.90) when a single anomalous year is dropped, points to a single conclusion: **the design cannot credibly identify the causal effect of ShotSpotter installation on gun homicides.** The data do not say ShotSpotter caused a 2.6-homicide-per-year increase; that statistic is a misspecification artifact. They also do not say ShotSpotter caused a 5-homicide-per-year decrease; that reading would be reading the event-study coefficients without acknowledging that the base year (2016) is anomalous. What the data do say is that the proper count-data specifications return statistically null effects with confidence intervals spanning the entire policy-relevant range.

We can bound the answer in one direction: the data are inconsistent with a large *deterrence* effect: a broad, sustained fall in the number of shootings of the kind detection volume alone was meant to produce. A treatment that large, with these sample sizes, would have produced visibly negative post-treatment event-study coefficients with confidence intervals that exclude zero from above, robustness across specifications, and a cleaner signal in the higher-power non-fatal shootings analysis. None of these appear.

What the design *cannot* rule out is a benefit operating through a different channel. A boundary regression-discontinuity analysis attributes roughly 85 fewer deaths per year to faster emergency response inside coverage, working through victim survival rather than deterrence (University of Chicago Crime Lab, 2024; a working analysis, not yet peer-reviewed). A mortality reduction of that magnitude falls squarely *inside* our own count-model confidence interval, and a community-area-by-year homicide panel (which observes neither response times nor the conversion of shootings into deaths) is simply not built to detect it. Ruling out a large deterrence effect is not the same as ruling out every benefit, and we do not claim the stronger statement. Whatever effect ShotSpotter had on gun homicides through deterrence was not large enough to dominate the noise, the 2016 anomaly, or the pre-existing trend differences across community areas.

## 7. Why Detection Did Not Reduce Violence

Three explanations are consistent with the empirical pattern.

**Selection.** ShotSpotter was installed in places that needed it most. Even if the technology produced a small benefit, the comparison would still look unfavorable because the rest of Chicago started from a much lower baseline with less room to deteriorate.

**Detection is not intervention.** The mechanism by which ShotSpotter could plausibly reduce homicides runs through faster police response and faster trauma triage. The data here do not measure response time, dispatch decisions, or trauma outcomes. The 2021 OIG audit found that ShotSpotter alerts rarely produced evidence used in investigations and that officers did not consistently change behavior based on alerts. If the operational link between detection and outcome is weak, detection volume is not a sufficient condition for fewer homicides.

**Structural factors dominate.** A surveillance technology cannot, by itself, change the distribution of poverty, unemployment, housing instability, school quality, or labor-market access. Programs that have been shown to reduce homicide in randomized or quasi-experimental evaluations (READI Chicago, Cure Violence, Becoming a Man) work on those structural channels directly. A microphone on a streetlight does not.

## 8. Limitations

The most important limitation is the one above: the design itself cannot identify a clean treatment effect from these data. We cannot recover the counterfactual we need (Chicago neighborhoods with ShotSpotter-eligible violence levels but no ShotSpotter installation) because the city did not generate one. Pre-period placebo tests confirm that whatever the DiD estimator captures, it is not specifically the effect of treatment.

Beyond that: the panel is annual, masking within-year dynamics; a restricted wild-cluster bootstrap leaves the OLS coefficient significant (p ≈ 0.001) in every specification, confirming the problem is the specification and the design rather than thin-cluster inference; non-fatal shooting counts in 2017–2023 are subject to coverage changes in the underlying database; the count-model null itself rests on a two-way fixed-effects Poisson/NegBin, which under a staggered rollout inherits the same heterogeneity-weighting concern as linear TWFE and would ideally be re-estimated with a nonlinear, heterogeneity-robust count DiD (Wooldridge, 2023; Moreau-Kastler, 2025), a caveat that, if anything, only reinforces the non-identification conclusion; the analysis does not observe police response times, dispatch decisions, or arrests; and the 26 never-treated community areas are systematically different from the treated areas on every observable dimension: wealth, demographics, education, employment.

## 9. Conclusion

ShotSpotter in Chicago succeeded at detection. Alerts climbed in covered areas because the technology measures gunfire that 911 underreports. What we cannot demonstrate is that detection translated into a reduction in lethal violence: every design we ran on these data either produces a positive DiD coefficient that fails falsification, an event study that violates parallel trends because of the 2016 anomaly, or a synthetic control whose donor pool cannot match treated-area violence levels. With no clean counterfactual and a pre-period dominated by one outlier year, these data cannot identify a causal effect in either direction. They are inconsistent with a large *deterrence* effect, though not with the smaller, survival-channel mortality benefit a boundary discontinuity design attributes to faster response, an effect this panel is not built to see.

That non-result is itself a policy statement. The 2023 termination of Chicago's $3 million contract is consistent with the absence, in the aggregate data, of the broad violence reduction the technology was marketed to deliver. Resources reallocated from detection contracts toward community-based violence-prevention programs (READI Chicago, Cure Violence, Becoming a Man) would, on the existing quasi-experimental evidence, have at least as good a chance of producing the homicide reductions ShotSpotter has not been shown to deliver at scale.

---

## Notes on figures, tables, and reproduction

This narrative draws on artifacts already produced by the code in this repository. To regenerate everything from scratch:

```bash
pip install -r requirements.txt
python scripts/did_analysis.py             # original DiD pipeline
python scripts/econometric_analysis.py     # TWFE, cohort event study, SC, placebos
python scripts/econometric_robustness.py   # Poisson/NegBin, Wald pre-trend test, SC inference
python scripts/did_callaway_santanna.py    # Callaway-Sant'Anna heterogeneity-robust DiD
python scripts/placement_model.py          # treatment-assignment / equity model (+ ACS race)
python scripts/robustness_closers.py       # NFS count models, wild-cluster bootstrap, Bacon decomposition
python scripts/make_pretty_plots.py        # publication figures
python scripts/plot_covid_sensitivity.py   # COVID-aware annualized comparison
```

Key artifacts referenced above:

- **Master DiD regression table:** [`docs/results_summary_original.md`](docs/results_summary_original.md)
- **Econometric re-analysis (TWFE, cohort, SC, placebos):** [`docs/results_summary_reanalysis.md`](docs/results_summary_reanalysis.md)
- **Placement / equity model (LPM, logit, propensity overlap, race):** [`docs/results_summary_placement.md`](docs/results_summary_placement.md)
- **Robustness closers (NFS count models, wild bootstrap, Bacon):** [`docs/results_summary_closers.md`](docs/results_summary_closers.md)
- **Goodman–Bacon decomposition:** [`plots/bacon_decomposition.png`](plots/bacon_decomposition.png)
- **Map of treatment vs. control community areas:** [`plots/did_map_treatment_control.png`](plots/did_map_treatment_control.png)
- **What distinguished ShotSpotter areas (placement):** [`plots/placement_coefficients.png`](plots/placement_coefficients.png)
- **Propensity-score overlap (selection-side view of no counterfactual):** [`plots/placement_propensity_overlap.png`](plots/placement_propensity_overlap.png)
- **Cohort-aware event study with bootstrap CIs:** [`plots/event_study_cohort_aware.png`](plots/event_study_cohort_aware.png)
- **Callaway–Sant'Anna heterogeneity-robust event study:** [`plots/callaway_santanna_event_study.png`](plots/callaway_santanna_event_study.png)
- **Six-panel synthetic control comparison:** [`plots/synthetic_control_panels.png`](plots/synthetic_control_panels.png)
- **Synthetic control permutation inference:** [`plots/synthetic_control_inference.png`](plots/synthetic_control_inference.png)
- **Original event study (TWFE on pooled cohort):** [`plots/did_event_study.png`](plots/did_event_study.png)
- **Forest plot of treatment effects across neighborhoods:** [`plots/forest_plot_effects.png`](plots/forest_plot_effects.png)
- **Robustness summary across specifications:** [`plots/robustness_summary.png`](plots/robustness_summary.png)
- **Non-fatal shootings DiD trend:** [`plots/nonfatal_did_trends.png`](plots/nonfatal_did_trends.png)
- **Annual trends per neighborhood:** [`plots/annual_trends_neighborhoods.png`](plots/annual_trends_neighborhoods.png)
- **COVID-adjusted before/during comparison:** [`plots/covid_sensitivity.png`](plots/covid_sensitivity.png)
