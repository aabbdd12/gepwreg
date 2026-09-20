# gepwreg: Percentile Weights Regression for Stata

`gepwreg` estimates how the effect of the covariates on an outcome varies
along a distribution, with the **percentile weights** of Araar (2016):
Gaussian kernel weights on the percentile scale, centred at a target
percentile. Version 1.4 estimates three objects:

* the **effect at the τ-quantile of the outcome** — the mean effect of the
  covariates among the households whose outcome is at that quantile (the
  unconditional quantile partial effect of Firpo, Fortin and Lemieux 2009);
* the **effect at the initial quantile** — among the households whose
  outcome *before* the contribution of the covariates is at that quantile
  (who benefits among the initially poor);
* the **profile of the effect along a covariate** — among the households at
  the τ-quantile of an observed variable (`rankvar()`).

The first two are estimated in two steps: a unit-level effect is obtained
for every household from a heterogeneity model fitted on the whole sample
(conditional quantile regression, `het(qr)`, the default; or interactions
with observed variables, `het(varlist)`), and the effects are averaged with
the percentile weights on the rank of the outcome. With `het(varlist)` the
estimator has a closed form — the heterogeneity of the effect times the
composition of the group — and an analytical variance by linearisation
that accounts for step 1, for the composition and for the estimated ranks;
Taylor linearisation applies under a complex survey design declared by
`svyset`; a bootstrap re-runs the whole procedure.

## What changed in 1.4

One release since 1.3. The version numbers in between were development states
and were never distributed.

### The estimator

Up to version 1.3 the effect at a quantile of the outcome was estimated by the
kernel-weighted regression of the outcome on the covariates with the weights on
the rank of the outcome itself. That regression selects the right households but
cannot measure the effect among them: inside a window on the outcome the outcome
does not vary, so the weighted slope converges to a fraction of the effect that
shrinks with the bandwidth (Goldberger 1981, selection on the dependent
variable). It equals the effect only when the outcome is a deterministic
monotone function of a single regressor, the case of the examples of the 2016
and 2023 notes; on survey data it is attenuated by one to two orders of
magnitude. It is kept as `rankdep`, labelled as a descriptive slope. The
derivation, the simulations and the application are in version 3 of the Zenodo
paper. Results obtained with version 1.3 on the outcome rank should be
re-estimated; results obtained with `rankvar()` are unaffected.

In its place, a two-step estimator: unit-level effects from a heterogeneity
model fitted on the whole sample — a conditional quantile regression, `het(qr)`,
or interactions with observed variables, `het(varlist)` — averaged with the
percentile weights on the rank of the outcome. New options `initial()` and
`xref()` give the effect at the initial quantile, the outcome before the
contribution of the listed regressors.

### Bandwidth

The bandwidth of the second step is chosen by the MSE-optimal plug-in rather
than by Silverman's rule, as it already was for `rankvar()` and `rankdep`. Step
2 averages the unit-level effects along the rank of the outcome, which is
uniform by construction, so the design-density term of the local-constant bias
vanishes and only the curvature of the profile survives; the coefficients do not
share units, so the rule minimises the average *relative* mean squared error and
one bandwidth serves the whole call. The unit-level effects do not depend on the
bandwidth, so the plug-in costs three weighted means and no refitting. On
simulations at n = 2,000 to 20,000 it sits on average 13 percent above the root
mean squared error of the best possible bandwidth, against 33 percent for
Silverman's rule, with a worst case of 52 percent against 216. `silverman`
restores the older behaviour and `band()` fixes the bandwidth by hand.

### Measurement error

`merr` corrects `het(varlist)` for classical measurement error in the
heterogeneity variables. A noisy heterogeneity variable does not simply flatten
the profile: the composition term is untouched, so the bias is a level shift
plus a shrinkage of the slope — and the level shift carries the **main** effect
of the heterogeneity variable, not its interaction, so the reported effect can
be over-stated several-fold while the interaction coefficient still looks
right.

One error variance per heterogeneity variable is estimated from the third
moments of the step-1 residual, identification coming from the skewness of the
variable or from its correlation with the regressors (Ben-Moshe,
D'Haultfœuille and Lewbel, 2017, specialised to the step-1 design). The
correction activates per variable only when `t = sigma2/se` exceeds `tcrit()`,
default 2; otherwise the command falls back to the uncorrected fit.
`e(merr_s2)`, `e(merr_t)` and `e(merr_keep)` report the decision, and
`e(merr_diag)` the implied reliability, the skewness and the R-squared on the
regressors, so a variable that is unidentifiable can be told from one this
sample merely fails to resolve.

The influence function accounts for the estimated variance, so the analytical
and survey-design standard errors remain valid where the correction is active;
`boot()` redoes the correction and its threshold in every resample. At the
threshold itself the estimator is a pre-test estimator and no standard error is
valid uniformly; report both profiles there.

The command refuses what it cannot identify: a heterogeneity variable taking
two values — if z takes only a and b then z² is affine in z, so the
third-moment conditions reduce to moments the first step has already set to
zero — a singular step-1 design, and a corrected moment matrix that is not
positive definite. Error in a categorical variable is misclassification, not
classical, and is not covered; nor is measurement error in the regressors
themselves. And the third moments cannot separate a missing functional form
from measurement error: fit the enriched model without `merr` to see which it
is.

### Standard errors

Taylor linearisation by default under `svyset` (`vce(if)` keeps the
influence-function ones), `e(V_IF)` always stored, and a bootstrap that re-runs
the whole procedure. The scores of `rankvar` and `rankdep` are centred — the
constant of the indirect term was omitted up to 1.3, a difference of second
order.

Standard errors no longer depend on the order of the observations. The rank
term of the influence function is summed over the set of units above *i* rather
than over the rows below it in the sorted data; on a ranking variable with seven
distinct values over 8,478 observations the reported standard errors moved by
8e-04 when the rows were permuted. The estimates were never affected.

### Weights

One rule, stated in the header of the output. A weight declared by `svyset` is
captured before anything is computed and used for everything — coefficients,
bandwidth, correction, ranks, variance — and a weight on the command line is
then ignored, with a note. The command weight is used only when `svyset`
declares none. `e(wtype)`, `e(wexp)` and `e(wsrc)` record the source.

### Also

Factor variables: every non-base level gets its own indicator, where up to 1.3
the second level was pooled with the base. A PSU count per stratum.
`gepwreg_setable`, a post-estimation comparison of the standard errors, in its
own file and available as the `setable` option. A rewritten dialog box covering
every option. Progress dots on every bootstrap path, with `nodots`. Execution
time in the header and in `e(etime)`.

## Installation

### Version 1.4 (September 2026) — current

```stata
net install gepwreg, from("https://raw.githubusercontent.com/aabbdd12/gepwreg/v14") replace
net get     gepwreg, from("https://raw.githubusercontent.com/aabbdd12/gepwreg/v14") replace
```

### Version 1.3 (June 2026) — to reproduce earlier results

```stata
net install gepwreg, from("https://raw.githubusercontent.com/aabbdd12/gepwreg/main") replace
net get     gepwreg, from("https://raw.githubusercontent.com/aabbdd12/gepwreg/main") replace
```

In both cases the second line retrieves `bkf98I.dta`, the Burkina Faso 1998
household survey extract used in the paper and in the help file. Stata 16 or
later; no dependencies.

`replace` overwrites whichever version is installed, so those are also the two
lines to switch between them. `which gepwreg` says which one you have: the
header of 1.4 is dated September 2026, that of 1.3 June 2026.

Results obtained with 1.3 on the rank of the outcome should be re-estimated —
see *What changed in 1.4* above. Results obtained with `rankvar()` are
unaffected.

The `replication/` folder of the `v14` branch holds the programmes that
produce every number of the paper, with their own README.

## Quick start

```stata
use bkf98I, clear
generate lexp  = ln(exppc)
generate male  = (sex == 1)
generate urban = (zone == 2)

gepwreg lexp size male i.gse [pw=weight], per(0.25) het(urban size)   // effect at the first quartile, IF s.e.
gepwreg lexp size male i.gse [pw=weight], per(0.25) boot(200)         // same, het(qr), bootstrap s.e.

svyset psu [pw=weight], strata(strata)
gepwreg lexp size male i.gse, per(0.25) het(urban size)               // Taylor s.e., automatically
gepwreg lexp urban size male i.gse, per(0.25) het(size) initial(urban)             // effect for the initially poor
gepwreg lexp male urban i.gse, per(0.25) rankvar(size)                // profile along household size
gepwreg_setable                                                       // the standard errors side by side
```

`help gepwreg` documents the objects, the estimators, the options, the
standard errors and the stored results.

## Documentation

One paper per version, and they are different papers.

**Version 1.4**, the `v14` branch:

Araar, A. (2026). *Effects at a quantile of the outcome: a two-step
percentile-weights estimator, with analytical and survey-design inference*.
Zenodo. [10.5281/zenodo.22845707](https://doi.org/10.5281/zenodo.22845707)

**Version 1.3**, the `main` branch:

Araar, A. (2026). *Exploring Heterogeneous Effects: Quantile Models and
Percentile Weights Regression*. Zenodo.
[10.5281/zenodo.20315684](https://doi.org/10.5281/zenodo.20315684)

Those records are what to read: the repository carries the software, not the
writing. The first DOI stands for all versions of its record and resolves to
the latest; the second points at the fixed version that documents 1.3. The
earlier statements of the method are Araar (2016), *Percentile weights
regression*, PEP technical note, and Araar (2023), PEP Working Paper 2023-15.

## Citing

For version 1.4:

```
Araar, A. (2026). Effects at a quantile of the outcome: a two-step
percentile-weights estimator, with analytical and survey-design inference.
Zenodo. https://doi.org/10.5281/zenodo.22845707
```

For results obtained with version 1.3:

```
Araar, A. (2026). Exploring Heterogeneous Effects: Quantile Models and
Percentile Weights Regression. Zenodo. https://doi.org/10.5281/zenodo.20315684
```

A `CITATION.cff` file is included for reference managers; on each branch it
carries that branch's version and its paper. The package itself is at
https://github.com/aabbdd12/gepwreg.

## License

MIT — see `LICENSE`. Stata is a registered trademark of StataCorp LLC.

## Author

Abdelkrim Araar, Université Laval and Partnership for Economic Policy (PEP)
— aabd@ecn.ulaval.ca
