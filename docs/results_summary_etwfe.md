# Heterogeneity-Robust Count Estimator: Wooldridge ETWFE-Poisson

Results of `scripts/etwfe_count.py`. The paper's headline Poisson is a plain two-way-FE
PPML, subject under staggered adoption to the same negative-weighting concern as linear
TWFE. This closes that gap with the estimator the count argument actually calls for: a
Poisson ETWFE (Wooldridge, 2023) that saturates the treatment in cohort and calendar
time, so each cohort-by-year ATT uses clean comparisons, then aggregates.

**Overall ATT: incidence-rate ratio 1.048, 95% CI [0.791, 1.330]** (community-area cluster bootstrap, 400 draws).

| Years since install (k) | ATT (IRR) | 95% CI |
|---|---|---|
| 0 | 1.123 | [0.842, 1.475] |
| 1 | 0.938 | [0.650, 1.340] |
| 2 | 0.967 | [0.663, 1.362] |
| 3 | 1.279 | [0.951, 1.849] |
| 4 | 0.975 | [0.738, 1.301] |
| 5 | 0.953 | [0.668, 1.296] |
| 6 | 1.299 | [0.884, 2.067] |

**Reading.** The correct heterogeneity-robust *count* estimator returns the same null as
the plain Poisson (IRR near 1 at every horizon, overall CI spanning 1). Closing the
acknowledged soft spot does not change the conclusion---if anything it removes the last
escape hatch: the null is not an artifact of the plain two-way-FE PPML's weighting.

Figure: `plots/etwfe_event_study.png`