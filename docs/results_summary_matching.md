# Does Matching on Pre-Treatment Observables Rescue the Design?

Results of `scripts/matched_did.py`. A standard response to a difference-in-differences
whose parallel-trends assumption fails is to match treated and control units on
pre-treatment observables first. We run that suggestion on its own terms.

The propensity model is the one from the placement analysis: a logit of ShotSpotter
placement on standardized pre-period gun homicides, hardship, and per-capita income
(AUC 0.935; 51 treated, 26 never-treated).

## Results

| Sample | Treated | Control | Study-6 retained | Pre-trend Wald F | p | OLS DiD | Poisson IRR |
|---|---|---|---|---|---|---|---|
| Full sample | 51 | 26 | 6/6 | 2.67 | 0.0160 | +2.481 (0.725) | 1.016 [0.843, 1.224] |
| Common support | 15 | 17 | 0/6 | 0.66 | 0.7044 | +0.812 (0.459) | 1.306 [0.967, 1.762] |
| NN matched | 17 | 8 | 0/6 | 0.87 | 0.5453 | +0.397 (0.479) | 1.201 [0.841, 1.715] |

Common-support region: propensity in [0.124, 0.873].
Nearest-neighbour matching used a caliper of 0.634 (0.2 SD of the linearized
propensity); 17 of 51 treated areas found a control inside it.

## Reading

Matching does not fail here. It works, and what it costs is the point.

**Restricting to comparable units does restore parallel pre-trends.** The joint Wald test
on pre-period interactions rejects parallel trends on the full sample (F = 2.67, p = 0.0160), but stops rejecting once the sample is trimmed to
the propensity common-support region (F = 0.66, p = 0.7044) and under
nearest-neighbour matching (F = 0.87, p = 0.5453). The $F$ statistic falls
well below one, so this is not only the loss of power that comes with a smaller sample,
though with fewer than half the original units some of it surely is. On its own terms the
suggestion works: comparable units really do move in parallel.

**The price is most of the coverage the study is about.** Placement is predicted well
enough that no never-treated area exceeds a propensity of 0.873, while 36 of the
51 ShotSpotter areas sit above that ceiling and 35 sit at or above 0.90,
against 0 of the 26 never-treated areas. Those 36 areas span 0.6 to 37.4 gun homicides a
year, so the failure is not confined to the violence tail. The six study neighborhoods,
carried over from the original analysis and used here because synthetic control needs a
readable number of panels, all sit inside that region: all six are present in the full sample; 0 of 6 survive common-support trimming and 0 of 6 survive matching, and only 17 of 51 treated areas find any control inside the caliper at all.
The design that satisfies parallel trends is therefore a design about the marginal,
low-propensity edge of ShotSpotter coverage: real areas, but not the neighborhoods the
program was built around, and not the ones the procurement debate was about.

**On that subsample the effect is still not estimable in any useful sense.** The
multiplicative estimate is IRR 1.31 [0.97, 1.76] on common
support and 1.20 [0.84, 1.72] matched, against 1.02 [0.84, 1.22] on the full panel. Both intervals
still span one, and both are wider than the full-sample interval, so trimming buys internal
validity and pays for it in precision.

This is the identification problem stated as a trade rather than a failure. These data can
deliver a credible comparison, or they can speak to the neighborhoods the policy is about,
but not both at once. That is a sharper statement of the paper's claim than the full-sample
diagnostics alone, and it is the reason the obstacle cannot be estimated away: no
reweighting of the units that exist can manufacture a counterfactual for units that have none.

**Specification note.** The pre-trend test here is a joint Wald test on treated-by-year
interactions in the 2009--2016 pre-period, which is a simpler test than the cohort event-study
Wald test reported in the main paper (F = 19.9, p = 0.006); the two are not the same statistic
and are not directly comparable. The full-sample OLS DiD reported here (+2.48) is the staggered,
cohort-specific specification, which is why it matches the Goodman-Bacon aggregate rather than
the pooled 2017-cutoff headline of +2.63.

Figure: `plots/matched_did.png`