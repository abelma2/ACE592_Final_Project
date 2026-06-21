# Results Summary

**ACE 592 Final Project: ShotSpotter & Gun Homicide DiD Analysis**

## DiD Regression Results

| Specification | DiD (beta3) | Std Error | t-stat | p-value | 95% CI | N | R-sq |
|---|---|---|---|---|---|---|---|
| 6 Neighborhoods, full period | +7.388 | 2.815 | 2.625 | 0.0087*** | [1.872, 12.904] | 1,155 | 0.4068 |
| All SS Neighborhoods, full period | +2.627 | 0.614 | 4.282 | 0.0000*** | [1.425, 3.830] | 1,155 | 0.1449 |
| 6 Neighborhoods, excl COVID | +4.455 | 2.415 | 1.844 | 0.0651* | [-0.279, 9.189] | 1,001 | 0.3968 |
| All SS Neighborhoods, excl COVID | +1.749 | 0.480 | 3.647 | 0.0003*** | [0.809, 2.689] | 1,001 | 0.1379 |
| 6 Neighborhoods, excl 2016 | +9.627 | 3.247 | 2.965 | 0.0030*** | [3.263, 15.990] | 1,078 | 0.4245 |
| All SS Neighborhoods, excl 2016 | +3.189 | 0.735 | 4.340 | 0.0000*** | [1.749, 4.630] | 1,078 | 0.1597 |
| 6 Neighborhoods, excl COVID+2016 | +6.694 | 2.823 | 2.371 | 0.0177** | [1.160, 12.227] | 924 | 0.4183 |
| All SS Neighborhoods, excl COVID+2016 | +2.311 | 0.595 | 3.886 | 0.0001*** | [1.145, 3.477] | 924 | 0.1551 |
| Non-fatal: 6 hoods, full | +21.194 | 12.251 | 1.730 | 0.0836* | [-2.817, 45.205] | 1,155 | --- |
| Non-fatal: All SS, full | +7.512 | 2.759 | 2.723 | 0.0065*** | [2.105, 12.919] | 1,155 | --- |

## Simple DiD Estimates (Annualized)

| Specification | Treatment Pre | Treatment Post | Control Pre | Control Post | DiD |
|---|---|---|---|---|---|
| 2A: Original 6 Neighborhoods | 134.1 | 188.6 | 289.1 | 408.9 | -65.3 |
| 2B: All ShotSpotter Neighborhoods | 381.2 | 541.9 | 42.0 | 55.6 | +147.0 |
| 2C-6:  6 Neighborhoods, excl COVID | 134.1 | 168.4 | 289.1 | 378.4 | -55.0 |
| 2C-all: All SS, excl COVID | 381.2 | 493.2 | 42.0 | 53.6 | +100.3 |
| Extra: 6 Hoods, excl 2016 | 118.7 | 188.6 | 265.7 | 408.9 | -73.3 |
| Extra: All SS, excl 2016 | 345.9 | 541.9 | 38.6 | 55.6 | +179.0 |
| Extra: 6 Hoods, excl COVID+2016 | 118.7 | 168.4 | 265.7 | 378.4 | -63.0 |
| Extra: All SS, excl COVID+2016 | 345.9 | 493.2 | 38.6 | 53.6 | +132.3 |

### Interpretation

The DiD estimate using all 51 ShotSpotter neighborhoods vs. the remaining 26 community areas is **+2.627** gun homicides per community area per year (p=0.0000).

A **positive** coefficient means ShotSpotter neighborhoods experienced a relative **increase** in gun homicides compared to the rest of Chicago.

## T-Test Results: Before vs. After ShotSpotter (Monthly Rates)

| Neighborhood | Before (mean/mo) | After (mean/mo) | Change | t-stat | p-value | Sig |
|---|---|---|---|---|---|---|
| Austin | 3.115 | 4.735 | +52.0% | -4.549 | 0.0000 | *** |
| Humboldt Park | 1.646 | 2.120 | +28.8% | -1.909 | 0.0581 | * |
| North Lawndale | 1.375 | 2.530 | +84.0% | -4.053 | 0.0001 | *** |
| Englewood | 1.812 | 2.012 | +11.0% | -0.703 | 0.4831 | ns |
| West Englewood | 1.885 | 1.952 | +3.5% | -0.265 | 0.7910 | ns |
| South Shore | 1.344 | 2.373 | +76.6% | -4.313 | 0.0000 | *** |
| **All 6 Combined** | 11.177 | 15.723 | +40.7% | -4.724 | 0.0000 | *** |

Note: Monthly counts include zero-homicide months. Welch's t-test (unequal variances). Before = Jan 2009 - Dec 2016; After = Feb 2017 - Dec 2023.

## ShotSpotter Coverage: All Identified Neighborhoods

| Community Area | Total Alerts | Coverage (months) | Earliest Alert | In Original 6? |
|---|---|---|---|---|
| Archer Heights | 1,076 | 69 | 2018-03-30 | No |
| Armour Square | 586 | 71 | 2018-01-22 | No |
| Ashburn | 3,001 | 71 | 2018-02-07 | No |
| Auburn Gresham | 7,290 | 75 | 2017-09-23 | No |
| Austin | 14,813 | 82 | 2017-03-16 | Yes |
| Avalon Park | 1,994 | 71 | 2018-01-30 | No |
| Avondale | 431 | 68 | 2018-05-07 | No |
| Belmont Cragin | 4,284 | 68 | 2018-04-24 | No |
| Bridgeport | 1,031 | 72 | 2017-12-22 | No |
| Brighton Park | 2,764 | 78 | 2017-06-24 | No |
| Burnside | 935 | 74 | 2017-10-19 | No |
| Calumet Heights | 2,452 | 71 | 2018-01-30 | No |
| Chatham | 5,895 | 75 | 2017-10-06 | No |
| Chicago Lawn | 6,539 | 83 | 2017-01-21 | No |
| Clearing | 762 | 69 | 2018-04-08 | No |
| Douglas | 583 | 67 | 2018-05-25 | No |
| East Garfield Park | 3,496 | 83 | 2017-02-15 | No |
| East Side | 2,536 | 71 | 2018-01-31 | No |
| Englewood | 8,479 | 84 | 2017-01-13 | Yes |
| Fuller Park | 672 | 78 | 2017-06-27 | No |
| Gage Park | 3,261 | 78 | 2017-06-20 | No |
| Garfield Ridge | 1,424 | 68 | 2018-04-13 | No |
| Grand Boulevard | 2,680 | 68 | 2018-05-16 | No |
| Greater Grand Crossing | 8,120 | 83 | 2017-01-23 | No |
| Hegewisch | 884 | 66 | 2018-07-05 | No |
| Hermosa | 1,407 | 68 | 2018-04-25 | No |
| Humboldt Park | 7,263 | 83 | 2017-01-19 | Yes |
| Hyde Park | 380 | 67 | 2018-05-28 | No |
| Kenwood | 817 | 68 | 2018-05-16 | No |
| Logan Square | 606 | 68 | 2018-04-27 | No |
| Lower West Side | 769 | 74 | 2017-11-01 | No |
| Mckinley Park | 1,087 | 71 | 2018-01-20 | No |
| Montclare | 786 | 68 | 2018-04-27 | No |
| Near West Side | 448 | 82 | 2017-03-09 | No |
| New City | 9,329 | 79 | 2017-06-16 | No |
| North Lawndale | 8,829 | 82 | 2017-02-26 | Yes |
| Oakland | 461 | 67 | 2018-05-18 | No |
| Pullman | 978 | 69 | 2018-03-13 | No |
| Riverdale | 1,922 | 70 | 2018-03-13 | No |
| Roseland | 8,702 | 75 | 2017-10-09 | No |
| South Chicago | 8,359 | 71 | 2018-01-29 | No |
| South Deering | 3,328 | 71 | 2018-01-30 | No |
| South Lawndale | 6,956 | 74 | 2017-10-28 | No |
| South Shore | 9,510 | 71 | 2018-01-30 | Yes |
| Washington Park | 2,785 | 71 | 2018-02-01 | No |
| West Elsdon | 1,077 | 69 | 2018-03-28 | No |
| West Englewood | 9,775 | 84 | 2017-01-14 | Yes |
| West Garfield Park | 5,550 | 82 | 2017-02-16 | No |
| West Lawn | 1,706 | 69 | 2018-03-30 | No |
| West Pullman | 7,962 | 70 | 2018-03-13 | No |
| Woodlawn | 3,703 | 71 | 2018-02-01 | No |

## Event Study Coefficients (Base Year = 2016)

| Year | Relative | Coefficient | 95% CI Low | 95% CI High |
|---|---|---|---|---|
| 2009 | -8 | -4.784 | -7.667 | -1.900 |
| 2010 | -7 | -4.480 | -7.690 | -1.270 |
| 2011 | -6 | -4.108 | -6.937 | -1.279 |
| 2012 | -5 | -3.606 | -6.630 | -0.583 |
| 2013 | -4 | -5.138 | -8.165 | -2.111 |
| 2014 | -3 | -4.845 | -7.789 | -1.900 |
| 2015 | -2 | -4.518 | -6.918 | -2.118 |
| 2016 (base) | -1 | +0.000 | 0.000 | 0.000 |
| 2017 | +0 | -0.681 | -2.647 | 1.285 |
| 2018 | +1 | -2.977 | -5.269 | -0.684 |
| 2019 | +2 | -3.489 | -5.977 | -1.000 |
| 2020 | +3 | -0.173 | -2.260 | 1.913 |
| 2021 | +4 | +1.949 | -0.337 | 4.236 |
| 2022 | +5 | -2.510 | -5.123 | 0.104 |
| 2023 | +6 | -1.273 | -3.683 | 1.137 |
