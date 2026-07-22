# The Procedural-Benefit Channel: Enforcement vs. Outcomes

Results of `scripts/procedural_benefit.py`. Gunshot detection is documented to
improve *process* measures but not *victim* outcomes (Piza et al., 2024). Chicago's
process analogs (ballistic evidence, gun recoveries) are FOIA-only; the public
proxy is weapons-violation enforcement from the Crimes dataset (`ijzp-q8t2`). We run
the *same* two-way fixed-effects DiD used for gun homicides on two enforcement
measures, so process and outcome are estimated on identical terms.

| Outcome | OLS DiD | Poisson IRR | 95% CI (IRR) | Poisson p |
|---|---|---|---|---|
| Gun homicides (outcome) | +2.627 | 1.074 | [0.869, 1.328] | 0.508 |
| Weapons-violation incidents (process) | +51.782 | 1.125 | [0.810, 1.563] | 0.481 |
| Weapons-violation arrests (process) | +21.961 | 0.952 | [0.620, 1.464] | 0.824 |

**Reading.** The result carries two lessons. First, the specification artifact is *general*. Weapons-violation incidents produce a large, significant OLS coefficient (+51.8 per area-year) that collapses to a statistically null incidence-rate ratio (1.13, p = 0.48) once the count outcome is modeled correctly---exactly the pattern the homicide analysis shows (OLS +2.63 vs.\ IRR 1.07). The additive-scale inflation is not peculiar to homicide counts; it recurs on a completely different outcome. Second, and substantively, the Chicago community-area data do *not* detect the process benefit that a microsynthetic control recovers in Kansas City (Piza et al., 2024): weapons incidents (IRR 1.13) and arrests (IRR 0.95) are both null on the correct scale. That is consistent with the 2021 OIG finding that Chicago's alerts rarely produced investigative evidence---and with the same missing-counterfactual problem that governs the rest of the paper.

**Caveats.** This is associational, not causal: the identification limits developed
throughout the paper apply here too (no valid control group; selection on violence).
The arrest flag reflects any arrest recorded as of the last data refresh, not
clearance. And a rise in weapons enforcement is double-edged---the procedural benefit
vendors advertise is also the concentrated-enforcement burden advocates warn of, since
coverage sits in the most over-policed, lowest-life-expectancy neighborhoods.

Figure: `plots/procedural_vs_outcome.png`