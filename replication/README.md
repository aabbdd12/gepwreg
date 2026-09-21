# Replication — Araar (2026), version 3

Percentile weights regression: what the outcome-ranked estimator measures, a
two-step estimator of the effect at a quantile, and inference under complex
survey design.

## What you need

Stata 16 or later, and `gepwreg` **1.4** — the version the paper was produced
with, which is this branch:

```stata
net install gepwreg, from("https://raw.githubusercontent.com/aabbdd12/gepwreg/v14") replace
net get     gepwreg, from("https://raw.githubusercontent.com/aabbdd12/gepwreg/v14") replace
```

The second line also brings `bkf98I.dta`, the Burkina Faso 1998 household
survey extract the application uses. Put it in the folder you run from.
`p3_replications.do`, `p4_twostep_vs_rif.do` and `p7_burkina.do` compare
against `rifhdreg` of Rios-Avila (2020), which is not part of this package:

```stata
ssc install rif
```

The command comes from the SSC package `rif`, not from a package of its own
name, and **the version matters**: SSC installs `rifhdreg 2.55` (August 2021),
which is the one these results were produced with. The Stata Journal archive
of the same article still ships version 2.5 (July 2019), whose header records
a bug in the definition of the estimation sample — the very thing a comparison
of two estimators on one sample depends on. Each of the three files begins
with `which rifhdreg`, which prints the version it found.

Each programme begins with `which gepwreg` and it must report **v1.4**. With
any other version the numbers will not be the ones in the paper — 1.3
estimates a different object, and later versions may change the bandwidth or
the variance.

Five of the programmes are in Python and need `numpy`; `pwr_exact_truth.py`
and `p10_rif_gap.py` also need `scipy`. `make_tables.py` needs neither.

## The tables of the paper

Run each file from the folder that contains it. The results are written to
`results/` beside it, which is created if it does not exist.

| file | section | writes | time |
|---|---|---|---|
| `p1_onestep_attenuation.do` | 3 | `p1_onestep.csv`, `p1_sweep.csv`, `p1_goldberger.csv` | minutes |
| `p2_closed_doors.do` | 3 | `p2_closed_doors.csv`, `p2_band.csv` | minutes |
| `p3_replications.do` | 3 | `p3_note2016.csv`, `p3_mc2023.csv` | minutes |
| `p4_twostep_vs_rif.do` | 6 | `p4_twostep.csv` | minutes |
| `p5_coverage_mc.do` | 6 | `p5_coverage.csv` | long |
| `p6_svy_mc.do` | 6 | `p6_svy.csv` | long |
| `p7_burkina.do` | 7 | `p7_tableA.csv` … `p7_tableD.csv` | 20–30 min |
| `p7d_tableD.do` | 7 | `p7_tableD.csv` | about 10 min |
| `p10_rif_gap.py` | 4 | `p10_rif_gap.csv` | about 1 min |

`p10_rif_gap.py` is Python, not Stata, and has to run before
`make_tables.py` like the rest. `p5` and `p6` are Monte Carlo experiments and
are the slow ones. `p7d` rebuilds Table D alone, without rerunning the whole application: it runs the
same code as the Table D section of `p7_burkina.do` and writes the same
`p7_tableD.csv`, byte for byte, which is a check in itself.

Then, with Python 3:

```
python make_tables.py
```

which reads `results/` and writes the LaTeX tables of the paper into
`tables/` beside it. It needs nothing beyond the standard library. If a
result is missing it stops and names the do-file that writes it, so the
tables are never built from a partial run.

## The checks and the derivations

These print their results; `p8_merr.do` and `p9_merr_curvature.do` also write
a CSV each. They are what the paper rests on where it does not report a
table.

| file | section | what it establishes | time |
|---|---|---|---|
| `gepwreg_gmm_check.do` | 4 | the two-step as a stacked GMM, and the four standard errors side by side | about 2 min |
| `gepwreg_test15.do` | — | the command against closed forms and known effects | minutes |
| `pwr_bandwidth_check.py` | 4 | the plug-in bandwidth against Silverman's rule and against the oracle | 2–3 min |
| `pwr_bandwidth_variance.py` | 5 | what choosing the bandwidth on the sample costs in variance | 1–2 min |
| `pwr_exact_truth.py` | 6 | the exact θ(τ) of the Monte Carlo designs, by numerical integration | about 5 min |
| `p8_merr.do` | 5 | measurement error in a heterogeneity variable, against a known truth; writes `results/p8_merr.csv` | about 1 min |
| `p9_merr_curvature.do` | 5 | the limit: the correction fires on a variable with no error but a curved effect; writes `results/p9_merr_curvature.csv` | about 1 min |
| `p8_truth.py` | 5 | the values of θ(0.75) the two files above compare against, with their convergence tables | about 2 min |

`pwr_bandwidth_variance.py` imports from `pwr_bandwidth_check.py`, so those
two must stay in the same folder. The times above were
measured on one machine and are indicative.

## On the numbers

The Monte Carlo files set their own seed. `p7_burkina.do` and
`p7d_tableD.do` set none and do not need one: `gepwreg` seeds its own
bootstrap, `rseed()` being 12345 unless you change it. Across major Stata
versions the random-number stream can differ; the conclusions do not depend
on it.

What reproduces, and how exactly. Every coefficient in this package returns
to the last bit, and so does every bootstrap standard error: Tables B and D
of the application, which carry ten bootstrap standard errors at B = 200 and
B = 50, come back identical byte for byte. The analytical standard errors can
move a little — by a few parts in a million for the Taylor standard errors of
Table A, and in the third decimal for the profile along household size of
Table C, where eight of the forty entries round differently. No coefficient
and no conclusion of the paper depends on that.

One number deserves its own warning. The standard error of a bootstrap
standard error is about 1/√(2B) — five percent at B = 200. Where the paper
compares a bootstrap standard error with an analytical one at the percent
level, as in `gepwreg_gmm_check.do`, B is in the thousands for that reason,
and lowering it will appear to produce a disagreement that is not there.

## What is not here

The paper's LaTeX source and its figures. The estimates are the object of
this folder; the typesetting is not.
