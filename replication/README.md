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
`p4_twostep_vs_rif.do` compares against `rifhdreg`, which is not part of this
package (`ssc install rifhdreg`); that file is the only one that needs it.

Each do-file begins with `which gepwreg` and it must report **v1.4**. With any
other version the numbers will not be the ones in the paper — 1.3 estimates a
different object, and later versions may change the bandwidth or the variance.

## Running

Run each file from the folder that contains it. The results are written to
`results/` beside it, which is created if it does not exist.

| file | section | output | time |
|---|---|---|---|
| `p1_onestep_attenuation.do` | 3 | `p1_onestep.csv`, `p1_sweep.csv`, `p1_goldberger.csv` | minutes |
| `p2_closed_doors.do` | 3 | `p2_closed_doors.csv`, `p2_band.csv` | minutes |
| `p3_replications.do` | 3 | `p3_note2016.csv`, `p3_mc2023.csv` | minutes |
| `p4_twostep_vs_rif.do` | 6 | `p4_twostep.csv` | minutes |
| `p5_coverage_mc.do` | 4 | `p5_coverage.csv` | long |
| `p6_svy_mc.do` | 4 | `p6_svy.csv` | long |
| `p7_burkina.do` | 5 | `p7_tableA.csv` … `p7_tableD.csv` | 20–30 min |
| `p7d_tableD.do` | 5 | `p7_tableD.csv` | about 10 min |

`p5` and `p6` are Monte Carlo experiments and are the slow ones.

Then, with Python 3:

```
python make_tables.py
```

which reads `results/` and writes the LaTeX tables of the paper into
`tables/` beside it. It needs nothing beyond the standard library. If a
result is missing it stops and names the do-file that writes it, so the
tables are never built from a partial run.

## What is not here

The paper's LaTeX source and its figures. The estimates are the object of this
folder; the typesetting is not.

The Monte Carlo files set their own seed. `p7_burkina.do` and
`p7d_tableD.do` set none and do not need one: `gepwreg` seeds its own
bootstrap, `rseed()` being 12345 unless you change it. Either way the results
reproduce exactly on the same Stata version. Across major Stata versions the
random-number stream can differ; the conclusions do not depend on it.
