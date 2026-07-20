# Robustness Closers: NFS Count Models, Wild-Cluster Bootstrap, Bacon Decomposition

Results of `scripts/robustness_closers.py` — three analyses that close limitations
named in the paper.

## 1. Non-fatal gunshot injuries as count models (TWFE)

The paper's OLS estimate is β = +7.51 (p = 0.009) on the full panel. Estimated as
count models (same two-way fixed effects, cluster-robust SEs; NB2 dispersion via the
Cameron–Trivedi auxiliary regression):

| Model | Specification | β | IRR | 95% CI (IRR) | p |
|---|---|---|---|---|---|
| Poisson | full panel | -0.1017 | 0.903 | [0.691, 1.181] | 0.457 |
| NegBin | full panel | -0.1684 | 0.845 | [0.640, 1.115] | 0.234 |
| Poisson | excl COVID | -0.0926 | 0.912 | [0.710, 1.171] | 0.469 |
| NegBin | excl COVID | -0.1573 | 0.854 | [0.653, 1.118] | 0.251 |
| Poisson | excl 2016 | -0.0868 | 0.917 | [0.686, 1.226] | 0.558 |
| NegBin | excl 2016 | -0.1748 | 0.840 | [0.623, 1.131] | 0.251 |

Source note: the victims file records non-fatal shootings from 2010, so 2009 counts
are near zero; this matches the OLS specification already reported.

## 2. Restricted wild-cluster bootstrap (headline OLS DiD)

Rademacher weights, null imposed, 999 replications, t-statistics cluster-robust with
a G/(G−1) small-sample factor (Cameron, Gelbach & Miller 2008). FWL-residualized
implementation verified against statsmodels TWFE to 1e-8.

| Specification | β (DiD) | t | Cluster-robust p | Wild-bootstrap p |
|---|---|---|---|---|
| full panel | +2.627 | +4.29 | 0.0000 | 0.0010 |
| excl COVID | +1.749 | +3.65 | 0.0005 | 0.0010 |
| excl 2016 | +3.189 | +4.35 | 0.0000 | 0.0010 |
| excl COVID + 2016 | +2.311 | +3.89 | 0.0002 | 0.0010 |

## 3. Goodman-Bacon decomposition (staggered TWFE, 2017/2018 cohorts)

Staggered TWFE estimate: +2.4813 (SE 0.725).
Weighted sum of 2x2 components: +2.4813 (exact reproduction, diff 4.44e-16).

| Comparison | Weight | 2x2 estimate |
|---|---|---|
| 2017 cohort vs never-treated | 0.375 | +4.642 |
| 2018 cohort vs never-treated | 0.517 | +0.946 |
| 2017 vs 2018 (not-yet-treated control) | 0.062 | +3.765 |
| 2018 vs 2017 (already-treated control) [forbidden] | 0.046 | +0.396 |

The 'forbidden' comparison (2018 cohort vs the already-treated 2017 cohort) carries
**4.6%** of the total weight. Negative-weighting contamination is therefore
minor in this panel; the sign fragility documented in the paper comes from the 2016
anomaly and the absence of a valid control group, not from staggered-timing bias in
the TWFE aggregation itself.

Figure: `plots/bacon_decomposition.png`