# Econometric Re-analysis: ShotSpotter & Gun Homicide DiD

This document presents the results of `scripts/econometric_analysis.py`,
which re-estimates the project's main causal claims using methods that
address the limitations of the original specification:

- **Two-way fixed effects** (community-area + year) instead of pooled DiD
- **Cohort-specific event study** that respects the staggered 2017/2018 rollout
- **Synthetic control** for each of the 6 study neighborhoods
- **Pre-period placebo tests** in 1999, 2003, 2007, and 2011

## 1. Two-Way Fixed Effects DiD

Specification: `gun_hom_it = β·(treat_i × post_t) + α_i + γ_t + ε_it`
with cluster-robust SEs at the community-area level.

### Gun homicides

| Specification | β (DiD) | SE | p | 95% CI | N |
|---|---|---|---|---|---|
| All 51 SS areas, 2009-2023 | +2.627 | 0.638 | 0.0000*** | [+1.376, +3.878] | 1,155 |
| All 51 SS areas, excl COVID (2020-21) | +1.749 | 0.502 | 0.0005*** | [+0.766, +2.732] | 1,001 |
| All 51 SS areas, excl 2016 | +3.189 | 0.767 | 0.0000*** | [+1.687, +4.692] | 1,078 |
| All 51 SS areas, excl COVID + 2016 | +2.311 | 0.624 | 0.0002*** | [+1.087, +3.535] | 924 |

### Non-fatal gunshot injuries

| Specification | β (DiD) | SE | p | 95% CI | N |
|---|---|---|---|---|---|
| All 51 SS areas, NFS, 2009-2023 | +7.512 | 2.870 | 0.0089*** | [+1.886, +13.138] | 1,155 |
| All 51 SS areas, NFS, excl COVID | +3.705 | 2.420 | 0.1258 | [-1.039, +8.449] | 1,001 |
| All 51 SS areas, NFS, excl 2016 | +10.496 | 3.578 | 0.0034*** | [+3.483, +17.508] | 1,078 |

## 2. Cohort-Specific Event Study

Cohorts: 21 community areas treated in 2017, 30 in 2018, 26 never-treated controls.

### Aggregated event-time ATT (sample-weighted across cohorts)

| Event time k | ATT | SE | 95% CI |
|---|---|---|---|
| -7 | -5.456 | 1.959 | [-9.296, -1.615] |
| -6 | -4.544 | 2.032 | [-8.527, -0.561] |
| -5 | -5.086 | 1.983 | [-8.973, -1.199] |
| -4 | -6.251 | 1.934 | [-10.041, -2.461] |
| -3 | -6.090 | 1.999 | [-10.008, -2.173] |
| -2 | -4.922 | 2.039 | [-8.919, -0.925] |
| -1 | +0.000 | 2.480 | [-4.860, +4.860] |
| +0 | -1.948 | 2.359 | [-6.572, +2.676] |
| +1 | -4.652 | 2.113 | [-8.794, -0.510] |
| +2 | -3.827 | 2.190 | [-8.119, +0.466] |
| +3 | +0.101 | 2.526 | [-4.850, +5.052] |
| +4 | -0.428 | 2.459 | [-5.248, +4.393] |
| +5 | -2.985 | 2.145 | [-7.189, +1.219] |

## 3. Synthetic Control

For each of the 6 study neighborhoods, weights over the 26 never-treated
community areas are chosen to minimize the pre-treatment squared error on
annual gun homicide counts. A positive 'gap' means the treated unit had
more homicides post-treatment than its synthetic counterfactual.

| Neighborhood | Cohort | Pre-tx RMSE | Post-tx gap (per year) | Top donors |
|---|---|---|---|---|
| Austin | 2017 | 33.66 | +48.57 | Washington Heights=1.00 |
| Humboldt Park | 2017 | 13.54 | +16.57 | Washington Heights=1.00 |
| North Lawndale | 2017 | 9.95 | +21.29 | Washington Heights=1.00 |
| Englewood | 2017 | 16.58 | +15.43 | Washington Heights=1.00 |
| West Englewood | 2017 | 16.48 | +14.29 | Washington Heights=1.00 |
| South Shore | 2018 | 9.74 | +20.33 | Washington Heights=1.00 |

> **Synthetic control is infeasible here — read these gaps with care.** For every
> treated unit the optimizer concentrates 100% of the weight on a single donor
> (Washington Heights, the highest-violence never-treated area) and the SLSQP solve
> does not converge to an interior optimum. The pre-treatment RMSE is large relative
> to the outcome level (e.g. Austin's pre-RMSE of ~34 against an actual level of
> 15-70/yr), so the synthetic counterfactual cannot reproduce the treated series even
> in the pre-period it was fit on. The 26 never-treated community areas are
> systematically lower-violence than any treated neighborhood, so no convex
> combination can match them. The positive 'gaps' below therefore reflect this
> level mismatch, **not** a treatment effect. See NARRATIVE.md §5.5.

## 4. Pre-Period Placebo Tests

Treatment is reassigned to fake years before the actual ShotSpotter rollout.
A non-null placebo coefficient suggests that the post-2017 estimate
reflects pre-existing trend differences rather than a causal effect.

| Fake treatment year | β | SE | p | 95% CI |
|---|---|---|---|---|
| 1999 | -1.344 | 0.507 | 0.0081*** | [-2.338, -0.349] |
| 2003 | -1.929 | 0.714 | 0.0069*** | [-3.328, -0.531] |
| 2007 | -0.645 | 0.497 | 0.1942 | [-1.618, +0.329] |
| 2011 | +0.256 | 0.347 | 0.4611 | [-0.425, +0.937] |

## Generated figures

- `plots/event_study_cohort_aware.png` — staggered event-study with bootstrap CIs
- `plots/synthetic_control_panels.png` — 6-panel actual vs. synthetic comparison

---
# Robustness Layer (econometric_robustness.py)
Three additions on top of the prior re-analysis:
- Count-data models (Poisson, Negative Binomial)
- Joint Wald test on pre-trends
- Permutation inference for synthetic control

## 5. Count-Data Regression
Coefficients are on the log-mean scale; IRR = exp(beta).

### Poisson with two-way FE

| Specification | beta | SE | p | IRR | IRR 95% CI |
|---|---|---|---|---|---|
| Full panel 2009-2023 | +0.0715 | 0.1082 | 0.5083 | 1.0742 | [0.8690, 1.3278] |
| Excl COVID (2020-21) | +0.0136 | 0.1105 | 0.9022 | 1.0137 | [0.8162, 1.2589] |
| Excl 2016 | +0.0838 | 0.1218 | 0.4912 | 1.0874 | [0.8565, 1.3806] |

### Negative Binomial (alpha=1) with two-way FE

| Specification | beta | SE | p | IRR | IRR 95% CI |
|---|---|---|---|---|---|
| Full panel 2009-2023 | -0.0302 | 0.1524 | 0.8427 | 0.9702 | [0.7197, 1.3079] |
| Excl COVID (2020-21) | -0.0835 | 0.1590 | 0.5997 | 0.9199 | [0.6736, 1.2564] |
| Excl 2016 | -0.0286 | 0.1676 | 0.8644 | 0.9718 | [0.6997, 1.3497] |

## 6. Joint Wald Test for Pre-Trends

Full TWFE event-study `gun_hom ~ C(ca_num) + C(year) + treat:1[year=k]`
for k in 2009..2023 (omitting 2016 as base).

H0: all 7 pre-treatment interaction coefficients (2009-2015) = 0

- F = 19.903
- p = 0.005782
- **Verdict:** REJECT parallel trends at the 1% level.

## 7. Permutation Inference for Synthetic Control

For each treated unit, the post-treatment gap (actual - synthetic) is
compared to the distribution of analogous gaps obtained by running
synthetic control on each of the 26 never-treated areas. Permutation
p-value = fraction of placebos with absolute gap >= absolute treated gap.

| Neighborhood | Cohort | Actual gap | Pre-tx RMSE | Placebo n | Perm. p |
|---|---|---|---|---|---|
| Austin | 2017 | +48.57 | 33.66 | 26 | 0.000 |
| Humboldt Park | 2017 | +16.57 | 13.54 | 26 | 0.000 |
| North Lawndale | 2017 | +21.29 | 9.95 | 26 | 0.000 |
| Englewood | 2017 | +15.43 | 16.58 | 26 | 0.000 |
| West Englewood | 2017 | +14.29 | 16.48 | 26 | 0.000 |
| South Shore | 2018 | +20.33 | 9.74 | 26 | 0.000 |

> **These p-values are invalid by construction — do not read them as treatment effects.**
> The permutation p of 0.000 for every unit is an artifact of the synthetic-control
> infeasibility documented in §3: the treated units' pre-treatment RMSE (10-34) is far
> larger than the placebo never-treated units' (which fit themselves well), so the treated
> 'gaps' land far outside a tight placebo distribution **because of pre-treatment level
> mismatch, not because of any post-treatment effect.** The standard RMSE-ratio guard
> cannot rescue this when no donor pool can match the treated units. The honest reading is
> that synthetic control is infeasible for these data, not that the effect is significant.
> See NARRATIVE.md §5.5.

Figure: `plots/synthetic_control_inference.png`
