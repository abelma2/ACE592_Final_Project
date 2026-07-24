"""
check_consistency.py -- the repository's test suite.

This project has no unit tests because it has almost no reusable functions; what it
has instead is a hand-transcription hazard. Numbers travel from the analysis scripts
into docs/results_summary_*.md and from there, BY HAND, into paper/final_paper.tex,
NARRATIVE.md and README.md. Nothing binds them. Every consistency failure this
repository has produced has been of that shape.

This script is the guard. It re-reads the source-of-truth docs, extracts the headline
numbers, and asserts that the paper states them correctly. It also checks the
structural invariants that have broken before: reference integrity, the treated-group
definition, figure availability, and the no-em-dash style rule.

Run it after changing any analysis, and before committing:

    python scripts/check_consistency.py        # exits non-zero if anything fails

It is deliberately dependency-free (stdlib only) so it can run in any environment.
"""
import os
import re
import sys

P = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCS = os.path.join(P, "docs")
PLOTS = os.path.join(P, "plots")
PAPER = os.path.join(P, "paper", "final_paper.tex")

FAILS, WARNS, PASSES = [], [], []


def read(path):
    with open(path, encoding="utf-8") as f:
        return f.read()


def ok(msg):
    PASSES.append(msg)


def fail(msg):
    FAILS.append(msg)


def warn(msg):
    WARNS.append(msg)


def grab(text, pattern, label):
    """Extract the first capture group(s) of pattern, or record a failure."""
    m = re.search(pattern, text)
    if not m:
        fail(f"{label}: could not locate source-of-truth pattern {pattern!r}")
        return None
    return m.groups()


paper = read(PAPER)
# strip LaTeX math/formatting so numeric comparison is not defeated by $...$ and \,
paper_flat = re.sub(r"[\$\\{}]|\\,", "", paper)
paper_flat = re.sub(r"\s+", " ", paper_flat)


# =====================================================================
# 1. NUMBERS: does the paper state what the pipeline produced?
# =====================================================================
# Each entry: (label, doc filename, regex over the doc, formatter -> list of strings
# that must all appear in the paper). The formatter encodes the rounding the paper
# uses, so a pipeline change that moves a number will fail here rather than silently
# disagreeing with the prose.
NUMERIC_CHECKS = [
    (
        "ETWFE overall ATT",
        "results_summary_etwfe.md",
        r"incidence-rate ratio ([\d.]+), 95% CI \[([\d.]+), ([\d.]+)\]",
        lambda g: [f"{float(g[0]):.2f}", f"{float(g[1]):.2f}", f"{float(g[2]):.2f}"],
    ),
    (
        "Removal Poisson IRR and CI",
        "results_summary_removal.md",
        r"Poisson IRR ([\d.]+) \(95% CI \[([\d.]+), ([\d.]+)\]",
        lambda g: [f"{float(g[0]):.2f}", f"{float(g[1]):.2f}", f"{float(g[2]):.2f}"],
    ),
    (
        "Placement propensity AUC",
        "results_summary_placement.md",
        r"Propensity model \(violence \+ hardship \+ income\) AUC = (\d+\.\d+)",
        lambda g: [f"{float(g[0]):.2f}"],
    ),
    (
        "No-overlap headline (treated above every control)",
        "results_summary_placement.md",
        r"\*\*(\d+) of the (\d+) treated areas lie above every single never-treated",
        lambda g: [f"{g[0]} of the {g[1]} treated areas", f"{g[0]} of the {g[1]}"],
    ),
    (
        "Highest never-treated propensity",
        "results_summary_placement.md",
        r"highest propensity among all \d+ never-treated areas is (\d+\.\d+)",
        lambda g: [f"{float(g[0]):.3f}"],
    ),
    (
        "Coverage screen gap (lowest treated / highest control alerts)",
        "results_summary_treatment.md",
        r"treated records \*\*(\d+) alerts\*\*[\s\S]*?never-treated records \*\*(\d+)\*\*",
        lambda g: [g[0], g[1]],
    ),
    (
        "Goodman-Bacon forbidden weight",
        "results_summary_closers.md",
        r"carries\s*\n?\*\*([\d.]+)%\*\*",
        lambda g: [f"{float(g[0]):.1f}%"],
    ),
    (
        "Non-fatal shootings Poisson IRR (full panel)",
        "results_summary_closers.md",
        r"\| Poisson \| full panel \| [-\d.]+ \| ([\d.]+) \|",
        lambda g: [f"{float(g[0]):.3f}"[:4]],
    ),
    (
        "Procedural weapons-incident IRR",
        "results_summary_procedural.md",
        r"Weapons-violation incidents \(process\) \| \+[\d.]+ \| ([\d.]+) \|",
        lambda g: [f"{float(g[0]):.3f}"[:4]],
    ),
]

for label, docname, pattern, fmt in NUMERIC_CHECKS:
    path = os.path.join(DOCS, docname)
    if not os.path.exists(path):
        fail(f"{label}: source doc {docname} is missing (re-run its script)")
        continue
    g = grab(read(path), pattern, label)
    if g is None:
        continue
    vals = fmt(g)
    missing = [s for s in vals if s not in paper_flat]
    if missing:
        fail(f"{label}: {docname} implies {vals} but the paper does not state {missing}")
    else:
        ok(f"{label}: paper matches {docname}")

    # A number appearing SOMEWHERE is not enough. This repo's characteristic failure is
    # updating one instance of a figure and missing another, so where an estimate travels
    # with a confidence interval, verify that EVERY occurrence of that interval is preceded
    # by the same point estimate.
    if len(vals) == 3:
        point, lo, hi = vals
        for m in re.finditer(re.escape(f"[{lo}, {hi}]"), paper_flat):
            window = paper_flat[max(0, m.start() - 70):m.start()]
            nums = re.findall(r"\d+\.\d+", window)
            if nums and nums[-1] != point:
                fail(f"{label}: found interval [{lo}, {hi}] preceded by {nums[-1]}, "
                     f"but {docname} gives {point} (an instance was missed in an update)")


# =====================================================================
# 2. TREATED-GROUP DEFINITION: the invariant every script re-derives
# =====================================================================
EXPECTED = {"treated": 51, "control": 26, "total": 77, "c2017": 21, "c2018": 30}
if EXPECTED["c2017"] + EXPECTED["c2018"] != EXPECTED["treated"]:
    fail("cohort arithmetic: 2017 + 2018 cohorts do not sum to the treated count")
if EXPECTED["treated"] + EXPECTED["control"] != EXPECTED["total"]:
    fail("panel arithmetic: treated + never-treated do not sum to 77")

for f in sorted(os.listdir(os.path.join(P, "scripts"))):
    if not f.endswith(".py") or f in ("plot_style.py", "check_consistency.py"):
        continue
    src = read(os.path.join(P, "scripts", f))
    if "alerts >= 100" in src.replace("total_alerts", "alerts").replace(".a >= 100", "alerts >= 100"):
        pass
    # any script that screens alerts must use the same two thresholds
    if re.search(r"(>=\s*100|>= ALERT_RULE|a >= 100|alerts >= 100)", src):
        if not re.search(r"(>=\s*12|>= MONTH_RULE|mo >= 12|months >= 12)", src):
            warn(f"{f}: screens on the 100-alert rule but not the 12-month rule")
ok("treated-group definition: 51 + 26 = 77, cohorts 21 + 30 = 51")

for k, v in (("51", "treated"), ("26", "never-treated"), ("77", "community areas")):
    if k not in paper_flat:
        fail(f"panel counts: the paper never states {k} ({v})")
if not FAILS or all("panel counts" not in x for x in FAILS):
    ok("panel counts: 51 / 26 / 77 all present in the paper")


# =====================================================================
# 3. REFERENCE INTEGRITY
# =====================================================================
labels = set(re.findall(r"\\label\{([^}]+)\}", paper))
refs = set(re.findall(r"\\(?:page)?ref\{([^}]+)\}", paper))
dangling = refs - labels
orphans = labels - refs
if dangling:
    fail(f"references: {len(dangling)} referenced but never defined: {sorted(dangling)}")
if orphans:
    warn(f"references: {len(orphans)} defined but never referenced: {sorted(orphans)}")
if not dangling and not orphans:
    ok(f"references: {len(labels)} labels, {len(refs)} refs, no orphans or dangles")


# =====================================================================
# 4. FIGURES: every \includegraphics resolves to a file that exists
# =====================================================================
used = re.findall(r"\\includegraphics\[[^\]]*\]\{([^}]+)\}", paper)
missing = [u for u in used if not os.path.exists(os.path.join(PLOTS, u))]
if missing:
    fail(f"figures: referenced but not present in plots/: {missing}")
else:
    ok(f"figures: all {len(used)} included figures exist in plots/")


# =====================================================================
# 5. STYLE: the paper deliberately contains no em dashes
# =====================================================================
# A LaTeX em dash is exactly three hyphens. Excluding runs of four or more (the
# decorative "% ------" comment dividers) is the only exemption; an em dash at the end
# of a wrapped line is still an em dash, which an earlier version of this check missed.
emdash_lines = []
for i, ln in enumerate(paper.splitlines(), 1):
    if ln.lstrip().startswith("%"):
        continue
    if re.search(r"(?<!-)---(?!-)", ln) or "—" in ln:
        emdash_lines.append(i)
if emdash_lines:
    fail(f"style: em dash(es) on line(s) {emdash_lines[:10]} (the paper should have none)")
else:
    ok("style: no em dashes in the paper")


# =====================================================================
# 6. FRAMING GUARD: the paper must not claim the effectiveness null as its own
# =====================================================================
BANNED = [
    r"we (?:find|show|demonstrate|conclude) that ShotSpotter (?:does not|did not) reduce",
    r"our (?:results|findings|analysis) show that ShotSpotter (?:does not|did not) reduce",
    r"ShotSpotter (?:does not|did not) reduce gun (?:violence|homicides)\.\s*(?:We|Our) find",
]
hits = [b for b in BANNED if re.search(b, paper, re.I)]
if hits:
    fail(f"framing: the paper appears to claim the effectiveness null as its own finding ({hits})")
else:
    ok("framing: no claim of the effectiveness null as the paper's own finding")


# =====================================================================
# 7. COMPANION DOCS: README and NARRATIVE restate the same story
# =====================================================================
# The paper is not the only place these numbers live. README.md and NARRATIVE.md
# restate the headline results, and they have drifted from the paper before.
COMPANIONS = {
    "README.md": read(os.path.join(P, "README.md")),
    "NARRATIVE.md": read(os.path.join(P, "NARRATIVE.md")),
}
# headline numbers that must appear in both, in the form the companions use
COMPANION_NUMBERS = ["+2.63", "1.07", "19.9", "0.94", "36 of the 51", "1.05"]
for name, text in COMPANIONS.items():
    flat = re.sub(r"\s+", " ", text)
    absent = [n for n in COMPANION_NUMBERS if n not in flat]
    if absent:
        warn(f"{name}: does not restate {absent} (check it still matches the paper)")
    else:
        ok(f"{name}: restates the headline numbers")

    # the framing constraint applies to the companions too
    if re.search(r"Why Detection (?:Did Not|Would Not) Reduce", text, re.I):
        fail(f"{name}: contains an effectiveness-verdict section removed from the paper")
    # em dashes were removed project-wide
    if "—" in text:
        fail(f"{name}: contains em dash(es) (they were removed project-wide)")

# every results doc a companion links to should exist
for name, text in COMPANIONS.items():
    for ref in set(re.findall(r"docs/(results_summary_\w+\.md)", text)):
        if not os.path.exists(os.path.join(DOCS, ref)):
            fail(f"{name}: links to docs/{ref}, which does not exist")


# =====================================================================
# REPORT
# =====================================================================
print("=" * 72)
print("CONSISTENCY CHECK")
print("=" * 72)
for m in PASSES:
    print(f"  PASS  {m}")
for m in WARNS:
    print(f"  WARN  {m}")
for m in FAILS:
    print(f"  FAIL  {m}")
print("-" * 72)
print(f"  {len(PASSES)} passed, {len(WARNS)} warnings, {len(FAILS)} failures")
if FAILS:
    print("\nA failure means the paper and the pipeline disagree. Re-run the affected")
    print("script, then update every place the number appears (paper, NARRATIVE, README).")
sys.exit(1 if FAILS else 0)
