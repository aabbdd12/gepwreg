# gepwreg: Percentile Weights Regression for Stata

`gepwreg` estimates the **Percentile Weights Regression (PWR)** of Araar
(2016, 2023): the local marginal effect of the regressors on an outcome, for
the observations whose rank in the distribution of the outcome — or of any
other ranking variable — lies in the neighbourhood of a target percentile.
Unlike conditional quantile regression, the ranking is on the outcome itself
rather than on its predicted part; unlike unconditional quantile regression,
the coefficient is a local slope rather than the marginal shift of a
quantile.

Version 1.3 implements the results of Araar (2026):

* **analytical standard errors** in closed form, by linearisation (the
  functional delta method), which account for the sampling variability of
  the estimated kernel weights — the naive weighted-least-squares standard
  errors understate the uncertainty by 30 to 45 percent;
* an **MSE-optimal bandwidth**, chosen by a two-step plug-in rule
  (the Silverman rule remains available with `silverman`);
* **complex survey design**: Taylor linearisation under strata and primary
  sampling units, used automatically once the data are `svyset`;
* a **pairs bootstrap** that re-estimates the kernel weights at every
  replication, for validation of the analytical standard errors;
* **generalised ranking** on any variable other than the outcome
  (`rankvar()`).

A dialog box (`db gepwreg`) covers every option.

## Installation

```stata
net install gepwreg, from("https://raw.githubusercontent.com/aabbdd12/gepwreg/main") replace
net get     gepwreg, from("https://raw.githubusercontent.com/aabbdd12/gepwreg/main") replace
```

The second line retrieves `bkf98I.dta`, the Burkina Faso 1998 household
survey extract used in the paper and in the help file. Stata 16 or later; no
dependencies.

## Quick start

```stata
use bkf98I, clear
generate lexp  = ln(exppc)
generate male  = (sex == 1)
generate urban = (zone == 2)

gepwreg lexp size male urban i.gse [pw=weight], per(0.5)      // IF-corrected s.e.

svyset psu [pw=weight], strata(strata)
gepwreg lexp size male urban i.gse, per(0.5)                  // Taylor s.e., automatically

gepwreg lexp size male urban i.gse [pw=weight], per(0.5) boot(500)
gepwreg_setable                                               // naive vs IF vs bootstrap s.e.
```

`help gepwreg` documents the options, the three standard-error estimators,
the bandwidth rules and the stored results.

## Documentation

The methods — the linearisation variance, the MSE-optimal bandwidth, the
Taylor variance under stratified cluster sampling and the Monte Carlo
evidence — are derived in the paper archived on Zenodo, which also documents
the command in an appendix:

Araar, A. (2026). *Exploring Heterogeneous Effects: Quantile Models and
Percentile Weights Regression — Analytical Standard Errors, MSE-Optimal
Bandwidth, and Inference under Complex Survey Design.* Zenodo.
[10.5281/zenodo.20315684](https://doi.org/10.5281/zenodo.20315684)

The source and PDF of that paper (`pwr_paper_rev.tex`, `pwr_paper_rev.pdf`)
are kept in this repository. The original notes are Araar (2016), *Percentile
weights regression*, PEP technical note, and Araar (2023), PEP Working Paper
2023-15.

## Citing

```
Araar, A. (2026). Exploring Heterogeneous Effects: Quantile Models and
Percentile Weights Regression. Zenodo. https://doi.org/10.5281/zenodo.20315684
```

A `CITATION.cff` file is included for reference managers. The package itself
is version 1.3.1 at https://github.com/aabbdd12/gepwreg.

## License

MIT — see `LICENSE`. Stata is a registered trademark of StataCorp LLC.

## Author

Abdelkrim Araar, Université Laval and Partnership for Economic Policy (PEP)
— aabd@ecn.ulaval.ca
