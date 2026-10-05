# Replication — Araar (2026), version 5

Percentile weights regression: what the outcome-ranked estimator measures, a
two-step estimator of the effect at a quantile, and inference under complex
survey design.

## What you need

Stata 16 or later, and `gepwreg` **1.5.0**, which is this branch:

```stata
net install gepwreg, from("https://raw.githubusercontent.com/aabbdd12/gepwreg/v15") replace
net get     gepwreg, from("https://raw.githubusercontent.com/aabbdd12/gepwreg/v15") replace
```

The second line also brings `bkf98I.dta`, the Burkina Faso 1998 household
survey extract the application uses. Put it in the folder you run from.

The programmes `p1` to `p10` and the checks come unchanged from the
replication of the previous version of the paper (branch `v14`), which they
name in their headers: they were run with `gepwreg` 1.4. Version 1.5.0
changes no result of version 1.4 except for a binary regressor that is also a
heterogeneity variable, which none of them uses, so they give the same
numbers with either. The programmes `p11` to `p14` are new in version 5 of
the paper (Section 7 and Section 9.6) and need version 1.5.0.

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

Each programme begins with `which gepwreg`: `p1` to `p10` name **v1.4**,
under which they were run, and give the same numbers with **v1.5.0**; `p11`,
`p12` and `p14` need **v1.5.0**; `p13` does not call `gepwreg`. Version 1.3
estimates a different object.

Five of the programmes are in Python and need `numpy`; `pwr_exact_truth.py`
and `p10_rif_gap.py` also need `scipy`. `make_tables.py` needs neither.

## The tables of the paper

The sections are those of version 5 of the paper.

Run each file from the folder that contains it. The results are written to
`results/` beside it, which is created if it does not exist.

| file | section | writes | time |
|---|---|---|---|
| `p1_onestep_attenuation.do` | 3 | `p1_onestep.csv`, `p1_sweep.csv`, `p1_goldberger.csv` | minutes |
| `p2_closed_doors.do` | 3 | `p2_closed_doors.csv`, `p2_band.csv` | minutes |
| `p3_replications.do` | 3 | `p3_note2016.csv`, `p3_mc2023.csv` | minutes |
| `p4_twostep_vs_rif.do` | 8 | `p4_twostep.csv` | minutes |
| `p5_coverage_mc.do` | 8 | `p5_coverage.csv` | long |
| `p6_svy_mc.do` | 8 | `p6_svy.csv` | long |
| `p7_burkina.do` | 9 | `p7_tableA.csv` … `p7_tableD.csv` | 20–30 min |
| `p7d_tableD.do` | 9 | `p7_tableD.csv` | about 10 min |
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

## The profile across quantiles (version 5: Section 7 and Section 9.6)

| file | section | what it gives | writes | time |
|---|---|---|---|---|
| `p11_profile_burkina.do` | 7.4 | Table 8: the group, the linear correction and the exact effect of the mean household under `het(qr)`, Burkina | `p11_profile_burkina.csv` | about 2 min |
| `p12_profile_dgp.do` | 7.4 | Table 9: the linear correction against the exact estimator where the truth is known | `p12_profile_dgp.csv` | about 10 min |
| `p13_profile_delta.do` | 7.5 | the delta-method standard error of the effect for a given household against the bootstrap, with the Powell and Hendricks–Koenker densities | `p13_profile_delta.log` and the draws | about 5 min |
| `p14_profile_figures.do` | 9.6 | Figures 1 and 2 and the numbers quoted with them | `fig_profile_hetz.png`, `fig_profile_hetqr.png`, `p14_profile_*.csv` | about 2 min |

`p12_profile_dgp.do` runs `p12_core_mata.do`, copies of the Mata functions
of `gepwreg` 1.5.0 that it calls; `p13_profile_delta.do` runs
`p13a_E_draws.do` and `p13b_E_densities.do`, which program the estimator
with official quantile regressions. All of them must stay in this folder.
Tables 8 and 9 are printed by `p11` and `p12` and were typed into the paper;
`make_tables.py` builds the tables of the other sections.

## The checks and the derivations

These print their results; `p8_merr.do` and `p9_merr_curvature.do` also write
a CSV each. They are what the paper rests on where it does not report a
table.

| file | section | what it establishes | time |
|---|---|---|---|
| `gepwreg_gmm_check.do` | 4 | the two-step as a stacked GMM, and the four standard errors side by side | about 2 min |
| `gepwreg_test15.do` | — | the command against closed forms and known effects | minutes |
| `pwr_bandwidth_check.py` | 4 | the plug-in bandwidth against Silverman's rule and against the oracle | 2–3 min |
| `pwr_bandwidth_variance.py` | 6 | what choosing the bandwidth on the sample costs in variance | 1–2 min |
| `pwr_exact_truth.py` | 8 | the exact θ(τ) of the Monte Carlo designs, by numerical integration | about 5 min |
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

The paper's LaTeX source. The estimates are the object of this folder, with
the two figures of Section 9.6, which `p14` writes; the typesetting is not.
