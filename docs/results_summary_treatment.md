# Is the Treated / Never-Treated Split a Knife-Edge Choice?

Results of `scripts/treatment_definition.py`. Every estimate in the paper rests on one
screening rule: a community area counts as covered if the historical alert file shows at
least 100 alerts spanning at least 12 months, giving **51 covered and 26 never-treated** areas. This documents how much that rule matters.

## 1. The cutoff falls in an empty region

The least-covered area coded as treated records **380 alerts**. The most-covered
area coded as never-treated records **99**. No community area in Chicago falls
between those two values, a gap of roughly 4 times with
nothing inside it. The screen therefore separates two genuinely distinct groups rather than
slicing an arbitrary point out of a continuum (`plots/treatment_threshold.png`).

## 2. The panel is invariant over a wide range of cutoffs

Number of areas coded as covered, by alert cutoff (rows) and month cutoff (columns):

| alerts >= | months >= 0 | months >= 6 | months >= 12 | months >= 18 |
|---|---|---|---|---|
| 10 | 57 | 57 | 56 | 56 |
| 50 | 54 | 54 | 54 | 54 |
| 100 | 51 | 51 | 51 | 51 |
| 250 | 51 | 51 | 51 | 51 |
| 500 | 47 | 47 | 47 | 47 |
| 1000 | 36 | 36 | 36 | 36 |

The panel is identical for any alert cutoff from 100 to 250, so the headline results do not depend on where in that range the line is drawn.

**At the operative cutoff the month criterion does no work.** Requiring 0 months and
requiring 18 months give the same 51 areas once the alert threshold is
100: any area with that many alerts also spans years of them. The criterion is not
vacuous in general, and it does bind at looser alert cutoffs (10), where it
screens out areas with a short burst of detections. We report the rule as implemented for
reproducibility and note that, at the threshold actually used, the alert count alone
determines the panel.

## 3. The comparison group is clean, but not perfectly

Of the 26 never-treated areas, **15 record zero alerts** and **11 record
between 1 and 99**. The five with the most:

| Community area | Alerts | Months | Alerts/month |
|---|---|---|---|
| Morgan Park | 99 | 69 | 1.4 |
| Washington Heights | 71 | 68 | 1.0 |
| Dunning | 63 | 68 | 0.9 |
| Portage Park | 40 | 63 | 0.6 |
| Albany Park | 17 | 8 | 2.1 |

This is a real limitation and we disclose it rather than smooth it over. ShotSpotter was
procured at the police-district level, and community areas cut across districts, so an area
adjacent to or partially inside a covered district can register a trickle of alerts without
being meaningfully covered. The magnitudes are small: the heaviest such case is about
1.4 alerts per month against a median of 2,536 alerts in covered areas over the same window, on the order of
one percent of treated intensity.

One case deserves explicit mention. **Washington Heights** (71 alerts over
68 months) is the single donor the synthetic control collapses onto, so the
primary comparison unit in that exercise is not perfectly untreated. Given the trace
volume this does not change the finding, and the synthetic control fails for the separate
and far larger reason that no donor can match treated-area violence levels at all. But a
reader checking the donor should know.

## Reading

The screening rule is not a knife-edge. It falls in an empty region of the data, the panel
is unchanged for cutoffs anywhere from 100 to 250, and only one of its two
stated criteria does any work. The comparison group is clean for most of its members and
carries trace contamination for a minority, an artifact of measuring a district-level
program on a community-area map. None of this rescues identification, which fails for the
reasons developed in the paper, but it does mean the failure is not an artifact of how
treatment was coded.

Figure: `plots/treatment_threshold.png`