# The September 2024 Removal as a Second Natural Experiment

Results of `scripts/removal_analysis.py`. Chicago ended ShotSpotter on 22 September 2024.
Commentary has read the gun-homicide decline that followed as proof the technology did
nothing (or caused harm). We apply the paper's discipline to the *removal* and find the
removal effect is no more identifiable than the installation effect---for the mirror-image
reason.

## 1. The decline predates the removal

Gun homicides per ShotSpotter area peaked in 2021 (13.7/yr) and had already fallen to 9.1 by 2024; never-treated areas moved in parallel (2.0 to 1.9). The city-wide decline from the 2021--2022 peak runs straight through the September 2024 shutoff in both groups.

## 2. Naive before/after vs. difference-in-differences

The commentary's comparison---gun homicides in coverage areas before vs.\ after the shutoff---shows a drop (0.849 to 0.575 per area-month, -32%). But never-treated areas dropped almost identically in proportional terms (0.199 to 0.154, -23%): the decline is city-wide, not a removal effect.

The removal even reproduces the paper's *scale* artifact. The additive OLS difference-in-differences returns a 'significant' removal benefit of -0.228 gun homicides per area-month (p = 0.001)---but only because coverage areas start from a far higher baseline (0.85 vs 0.20 per area-month) and so swing more in absolute terms for the same proportional move. On the correct multiplicative scale the removal effect is null: Poisson IRR 0.88 (95% CI [0.53, 1.46], p = 0.62). This is the identical additive-vs-multiplicative divergence documented for the installation in the count-data section, now recurring on the removal.

## Reading

The symmetry is the point. The removal reproduces both mistakes at once: the raw before/after attributes a city-wide decline to the shutoff, and the additive DiD manufactures a 'significant' effect that the multiplicative model erases. The same two errors that fabricate an installation effect fabricate a removal effect, and the same missing counterfactual makes both non-identifiable. The post window is short (through the last complete month of the data) and 2025 is partial, so the DiD is also under-powered; the honest conclusion is that these data cannot identify the removal effect in either direction, exactly as they cannot identify the installation effect.

Figure: `plots/removal_trajectory.png`