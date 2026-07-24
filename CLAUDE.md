# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

A reproducible econometrics study, **not a software application**. The deliverable is an academic paper — `paper/final_paper.tex` (→ `final_paper.pdf`) — arguing that the causal effect of Chicago's ShotSpotter gunshot-detection system on gun homicides **cannot be identified** from public community-area data (the positive OLS difference-in-differences coefficient is a count-data scale artifact, and there is no valid control group). `scripts/` is the evidence pipeline behind that argument.

**Framing constraint (load-bearing).** The paper's *own* contribution is **non-identification — "cannot identify in either direction."** It is NOT "ShotSpotter does not reduce violence" (that null is prior work by others — Connealy 2024, Doucette 2021, Huff 2026 — which the paper cites, not claims). Do not let edits drift the paper's own claim into an effectiveness verdict.

## Commands

```bash
pip install -r requirements.txt          # pandas / statsmodels / geopandas stack, Python 3
```

Each analysis script is **standalone** — it reads `data/raw/`, and writes result tables to `docs/results_summary_*.md` and figures to `plots/*.png`. There is no top-level driver; run them individually:

```bash
python scripts/econometric_analysis.py    # TWFE DiD, cohort event study, synthetic control
python scripts/econometric_robustness.py  # Poisson/NegBin DiD, pre-trend Wald test, permutation SC
python scripts/did_callaway_santanna.py   # Callaway–Sant'Anna — needs the `differences` package (pinned 0.3.x)
python scripts/placement_model.py         # treatment-assignment / equity model + propensity overlap
python scripts/etwfe_count.py             # Wooldridge (2023) ETWFE-Poisson robustness
python scripts/removal_analysis.py        # Sept-2024 shutoff as a second natural experiment
# ...one script per analysis; README "Reproducing the analysis" lists the canonical run order
```

Build the paper — run `pdflatex` **twice** so cross-references resolve (keep the `.aux` between passes):

```bash
cd paper && pdflatex -interaction=nonstopmode -halt-on-error final_paper.tex && pdflatex -interaction=nonstopmode -halt-on-error final_paper.tex
```

Same for `methods_memo.tex`. A clean build has no `undefined` in the log and no LaTeX warnings; `final_paper.pdf` is ~36 pages.

**There is a consistency check, and it is the closest thing here to a test suite.** Run it after changing any analysis and before committing:

```bash
python scripts/check_consistency.py     # exits non-zero on failure; stdlib only
```

It re-reads the source-of-truth `docs/results_summary_*.md` and asserts the paper states those numbers, verifies every occurrence of a confidence interval carries the same point estimate (the characteristic failure here is updating one instance of a number and missing another), and checks reference integrity, figure availability, panel-count arithmetic, the no-em-dash rule, and the framing constraint below. It has been validated by deliberately re-introducing three past regressions and confirming each is caught.

Beyond that, verifying a change still means: re-run the affected script, inspect its `docs/results_summary_*.md`, and recompile the paper.

## Architecture & the consistency hazard

```
data/raw/*.csv  →  scripts/*.py  →  docs/results_summary_*.md   (the numbers)
                                 →  plots/*.png                  (the figures, 200 DPI)
                                 →  paper/*.tex                  (\includegraphics + hand-transcribed numbers)
```

**The #1 hazard: numbers in the paper are hand-transcribed from the scripts' `docs/results_summary_*.md` outputs — nothing binds them automatically.** When you change an analysis you must (a) re-run its script to refresh the doc + figure, (b) update by hand every place the number appears in `paper/final_paper.tex`, and (c) keep the companion write-ups in sync: `NARRATIVE.md`, `paper/methods_memo.tex`, and `README.md` all restate the same headline numbers and story. A single result (e.g. an incidence-rate ratio) typically appears in ~5 places across these files.

**Treated-group definition (recurs across scripts, must stay consistent).** A community area is "treated" (ShotSpotter) if the historical alerts file shows **≥100 alerts spanning ≥12 months** → **51 treated, 26 never-treated (77 total)**, in two staggered cohorts: **2017 (n=21) and 2018 (n=30)**. Scripts re-derive this independently via the `ALL_SS = set(prof[(prof.a >= 100) & (prof.mo >= 12)].index)` pattern. Changing the threshold changes the entire panel and every downstream number.

**Shared plotting contract.** Every figure script imports `scripts/plot_style.py` and calls `apply_academic_style()`, draws from the `COLORS` palette (`treated` / `control` keys), and stamps provenance with `add_source_note()` / `add_n_label()` / `SOURCE_DEFAULT`. New figures should reuse these helpers so the paper's plots stay visually uniform.

**Appendix structure.** The procedural/weapons-enforcement test is Appendix A (`app:procedural`) and the granular event-study coefficient table is Appendix B (`app:eventstudy`); the main body carries a one-line pointer to each. `\appendix` sits before the (unnumbered) Data Availability and References sections — keep it there so those back-matter sections stay unnumbered.

## Data & git gotchas

- Three large raw files (victims, ShotSpotter alerts, transportation; ~110 MB combined) are **git-ignored** — re-download from the Chicago Data Portal links in `README.md` into `data/raw/`. The committed CSVs/GeoJSON already cover the descriptive and placement scripts.
- Stage commits with **explicit paths, not `git add -A`** — raw-data and one-off PDF files here have been swept into commits that way before (`support_text/` is git-ignored for the same reason).
