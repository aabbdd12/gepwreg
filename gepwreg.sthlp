{smcl}
{* *! gepwreg.sthlp  v1.0  Araar A.  2024}{...}
{vieweralsosee "qreg" "help qreg"}{...}
{vieweralsosee "rifreg" "help rifreg"}{...}
{vieweralsosee "gepwe" "help gepwe"}{...}
{viewerjumpto "Syntax" "gepwreg##syntax"}{...}
{viewerjumpto "Description" "gepwreg##description"}{...}
{viewerjumpto "Options" "gepwreg##options"}{...}
{viewerjumpto "Stored results" "gepwreg##results"}{...}
{viewerjumpto "Standard errors" "gepwreg##se"}{...}
{viewerjumpto "Bandwidth" "gepwreg##bandwidth"}{...}
{viewerjumpto "Examples" "gepwreg##examples"}{...}
{viewerjumpto "References" "gepwreg##references"}{...}
{hline}
{title:Title}

{phang}
{bf:gepwreg} {hline 2} Percentile Weights Regression with analytical
standard errors

{hline}
{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:gepwreg}
{depvar}
[{indepvars}]
{ifin}
{weight}
{cmd:,}
[{it:options}]

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt per:centile(#)}}target percentile {it:tau}; default {cmd:percentile(0.5)}{p_end}
{synopt:{opt cband(#)}}Silverman bandwidth constant; default {cmd:cband(0.9)}{p_end}
{synopt:{opt band(#)}}manual bandwidth override (0 = use Silverman){p_end}
{synopt:{opt rankvar(varname)}}ranking variable; default = {it:depvar} (the outcome){p_end}
{syntab:Standard errors}
{synopt:{it:(automatic)}}Taylor SE if {cmd:svyset} declared; IF-corrected otherwise{p_end}
{synopt:{opt vce(svy)}}force Taylor SE explicitly (optional){p_end}
{synopt:{opt boot(#)}}pairs bootstrap — command line only; default {cmd:boot(0)}{p_end}
{syntab:Display}
{synopt:{opt l:evel(#)}}confidence level; default {cmd:level(95)}{p_end}
{synopt:{opt noc:onstant}}suppress constant term{p_end}
{synoptline}

{phang}
{cmd:fweight}s, {cmd:aweight}s, and {cmd:pweight}s are allowed; see
{help weight}.{p_end}

{phang}
After estimation, {cmd:gepwreg_setable} displays a comparison of naive,
IF-corrected, and bootstrap SE. The {cmd:boot(B)} option (command line only)
is required for bootstrap SE.

{hline}
{marker description}{...}
{title:Description}

{pstd}
{cmd:gepwreg} estimates the Percentile Weights Regression (PWR) model of
{browse "http://dasp.ecn.ulaval.ca/stata_adds/gepwe/PEP_Notes_Araar_01.pdf":Araar (2016)}.
It estimates the local marginal effect of the regressors on the outcome
for observations whose {it:outcome rank} is in the neighbourhood of a
target percentile {it:tau}:

{pmore}
{it:beta(tau)} = d{it:y}/d{it:x} |_{{it:p} = {it:tau}}

{pstd}
Unlike conditional quantile regression ({help qreg}), the percentile
ranking is based on the {it:dependent variable} rather than the predicted
component. Unlike unconditional quantile regression
({help rifreg}), the coefficient measures a local slope rather than a
marginal location shift.

{pstd}
The three-step procedure is:

{phang2}
(1) Compute the weighted empirical CDF rank {it:p-hat_i} = F-hat_n({it:y_i}).{p_end}

{phang2}
(2) Assign Gaussian kernel weights centred at {it:tau}:{p_end}

{pmore2}
{it:w_i} = exp(-0.25*((p-hat_i - tau)/h)^2) / (h * sqrt(2*pi) * n)

{phang2}
(3) Estimate beta by weighted least squares (WLS) using {it:w_i}.{p_end}

{pstd}
{bf:Standard errors} are selected automatically: if {cmd:svyset} has been
declared, Taylor linearisation SE is used; otherwise IF-corrected analytical
SE (Deville 1999) is applied. Naive WLS SE are {it:inconsistent} and stored
in {cmd:e(V_naive)} for reference only. See {help gepwreg##se:Standard errors} below.

{hline}
{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt percentile(#)} specifies the target percentile {it:tau} in the open
interval (0, 1). The default is {cmd:percentile(0.5)} (the median).

{phang}
{opt cband(#)} sets the Silverman (1986) bandwidth constant. The
bandwidth is computed as {it:h} = {it:cband} * {it:sigma_p} * {it:n}^(-1/5),
where {it:sigma_p} is a robust scale estimate of the percentile ranks.
The default is {cmd:cband(0.9)}.

{phang}
{opt band(#)} overrides the automatic bandwidth with a fixed value.
Useful for sensitivity analysis or when comparing results across
different percentiles with a common bandwidth. Setting {cmd:band(0)}
(the default) restores the Silverman rule.

{phang}
{opt rankvar(varname)} specifies an alternative variable {it:z}
on which to compute percentile ranks. {bf:By default, the ranking variable
equals the outcome variable} {it:depvar} --- this is the standard PWR of
Araar (2016). When {opt rankvar(z)} is specified, the rank of observation
{it:i} is computed as {it:p_i = F_n(z_i)}, so that {cmd:gepwreg} estimates
the marginal effect of {it:indepvars} on {it:depvar} for observations at
the {it:tau}-th percentile of {it:z}.

{pmore}
This generalisation is useful when interest lies in distributional effects
along a dimension different from the outcome. For example,
{cmd:rankvar(income)} with {cmd:depvar(lexp)} identifies the effect of
household characteristics on consumption for households at each percentile
of the {it:income} distribution.

{pmore}
The bandwidth formula, SE derivation (including the indirect IF term),
and all asymptotic results remain valid under {it:F_Z}. In the IF score,
the indirect term sums over observations with {it:z_i >= z_j} instead of
{it:y_i >= y_j}. All three SE types (SE_IF, SE_SVY, bootstrap) respect
the z-based ordering.

{phang}
{opt vce(svy)} explicitly requests Taylor linearisation SE under the
survey design declared by {helpb svyset}. This option is {bf:optional} —
{cmd:gepwreg} automatically uses Taylor SE when {cmd:svyset} has been
declared, without requiring {opt vce(svy)}.
The PSU, strata, and probability weights are read directly from {cmd:svyset}.
At least 2 PSUs per stratum are required; a warning is issued when
any stratum has fewer than 5 PSUs.

{phang}
{opt boot(#)} requests pairs bootstrap SE with {it:#} replications.
At each replication, kernel weights are {it:fully reconstructed} from
the resampled data. This option is available {bf:via command line only}
(not in the dialog box). Setting {cmd:boot(0)} (the default) skips the
bootstrap. Note: the pairs bootstrap ignores PSU and strata; for
survey-consistent bootstrap SE use the external
{cmd:bootstrap, strata() cluster() :} prefix. We recommend {it:#} >= 1000.

{phang}
{opt level(#)} sets the confidence level for the displayed confidence
intervals. The default is controlled by {help set level}.

{phang}
{opt noconstant} suppresses the constant term in the regression.

{hline}
{marker se}{...}
{title:Standard errors}

{pstd}
{cmd:gepwreg} selects the SE type {bf:automatically}:

{phang2}
{bf:svyset declared} (PSU and/or strata): Taylor linearisation SE.{p_end}
{phang2}
{bf:No svyset}: IF-corrected analytical SE (Deville 1999).{p_end}

{pstd}
Three SE estimators are always computed and stored:

{p2colset 5 28 30 2}{...}
{p2col:{bf:e(V)}}main SE: Taylor if svyset declared, IF otherwise{p_end}
{p2col:{bf:e(V_naive)}}naive WLS SE (inconsistent — for reference only){p_end}
{p2col:{bf:e(V_boot)}}bootstrap SE (stored if {cmd:boot()}>0){p_end}

{pstd}
{ul:Naive SE (inconsistent).}
The standard WLS formula treats the kernel weights {it:w_i} as fixed
constants. Because {it:w_i} depends on the empirical CDF of {it:y}
(which is itself a statistic), ignoring this dependence yields standard
errors that underestimate the true uncertainty. In practice the
underestimation is typically 30-45%.

{pstd}
{ul:Analytical IF-corrected SE (consistent).}
{cmd:gepwreg} accounts for the two-step structure (CDF estimation + WLS)
using the influence-function variance formula of Deville (1999):

{pmore}
V-hat = A^{-1} * (1/n^2 * sum_j psi_j*psi_j') * A^{-1}

{pstd}
where the total influence of observation {it:j} is:

{pmore}
psi_j = w_j * x_j * e_j  +  (1/n) * sum_{{i: y_i >= y_j}} dw_i * x_i * e_i

{pstd}
The first term is the direct score contribution; the second captures
the effect of observation {it:j} on the CDF rank (and hence the kernel
weight) of all observations above it. The derivative of the kernel
weight is:

{pmore}
dw_i = -0.5 * (p-hat_i - tau) / h^2 * w_i

{pstd}
This formula is computed in O(n*log(n)) time using a reverse cumulative
sum after sorting by {it:y}.

{pstd}
{ul:Bootstrap SE.}
When {cmd:boot(#)} > 0, a pairs bootstrap is performed. At each
replication, the kernel weights are {it:fully reconstructed} from the
resampled outcome and survey weights — not held fixed. This is the
numerically correct bootstrap for the PWR estimator and serves as the
gold standard for validating the analytical SE.

{pstd}
{ul:Recommendation.}
Declare {cmd:svyset} before {cmd:gepwreg} when working with complex
survey data — Taylor SE is applied automatically.
For data without survey design, IF-corrected SE is applied automatically.
For bootstrap validation, use {cmd:boot(B)} via command line and
{cmd:gepwreg_setable} to compare SE estimates.

{hline}
{marker bandwidth}{...}
{title:Bandwidth selection}

{pstd}
The default bandwidth follows Silverman (1986):

{pmore}
h = cband * sigma_p * n^(-1/5)

{pstd}
where sigma_p = min(sd(p-hat), IQR(p-hat)/1.34) is a robust scale
of the percentile ranks. The constant {it:cband} = 0.9 is the standard
Gaussian kernel multiplier.

{pstd}
This rule was designed to minimise the integrated squared error of a
{it:density estimate}, not the MSE of regression coefficients. It is a
reasonable default but may not be optimal. In particular:

{phang2}
(i) When the local fit is strong (high local R-squared near {it:tau}),
a smaller bandwidth may be preferred to sharpen the local estimate.{p_end}

{phang2}
(ii) Near the boundary (tau close to 0 or 1), the Silverman rule may
produce too few effective observations.{p_end}

{pstd}
For sensitivity analysis, re-run {cmd:gepwreg} with different values of
{cmd:band()} or {cmd:cband()} and compare results. The effective sample
size {cmd:e(N_eff)} reported after estimation provides a useful
diagnostic: values below 50 suggest the bandwidth may be too narrow.

{hline}
{marker results}{...}
{title:Stored results}

{pstd}
{cmd:gepwreg} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(tau)}}target percentile{p_end}
{synopt:{cmd:e(h)}}bandwidth used{p_end}
{synopt:{cmd:e(N_eff)}}effective sample size{p_end}
{synopt:{cmd:e(boot)}}number of bootstrap replications{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}main VCV: Taylor if svyset declared, IF otherwise{p_end}
{synopt:{cmd:e(V_IF)}}IF-based analytical VCV (always stored){p_end}
{synopt:{cmd:e(V_naive)}}naive WLS variance matrix (inconsistent, reference only){p_end}
{synopt:{cmd:e(V_boot)}}bootstrap variance matrix (if {cmd:boot()}>0){p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:gepwreg}{p_end}
{synopt:{cmd:e(rankvar)}}ranking variable (= depvar if default){p_end}
{synopt:{cmd:e(ranklabel)}}label: "{it:varname} (instead of depvar)" or "y (outcome)"{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(title)}}{cmd:Percentile Weights Regression}{p_end}
{synopt:{cmd:e(SE_type)}}SE type applied: {cmd:IF-corrected} or {cmd:Taylor linearisation}{p_end}
{synopt:{cmd:e(wgt)}}weight variable used (e.g. {cmd:PW[weight]}){p_end}

{hline}
{marker examples}{...}
{title:Examples}

{pstd}{ul:Basic usage}{p_end}

{phang2}{cmd:. use data.dta}{p_end}
{phang2}{cmd:. gen fw = hhsize * wta_hh}{p_end}
{phang2}{cmd:. gen lexp = log(pcexp)}{p_end}

{phang2}{cmd:. gepwreg lexp i.hhedlevel hhagey hhsize i.hhsex i.rururb [fw=fw], per(0.5)}{p_end}

{pstd}{ul:Compare with OLS and quantile regression}{p_end}

{phang2}{cmd:. gepwreg  lexp i.hhedlevel hhagey hhsize [fw=fw], per(0.5)}{p_end}
{phang2}{cmd:. eststo pwr50}{p_end}
{phang2}{cmd:. qreg   lexp i.hhedlevel hhagey hhsize [fw=fw], quantile(0.5)}{p_end}
{phang2}{cmd:. eststo qr50}{p_end}
{phang2}{cmd:. esttab pwr50 qr50}{p_end}

{pstd}{ul:Multiple percentiles}{p_end}

{phang2}{cmd:. foreach tau in 0.10 0.25 0.50 0.75 0.90 {c -(}}{p_end}
{phang2}{cmd:.     gepwreg lexp hhagey hhsize [fw=fw], per(`tau')}{p_end}
{phang2}{cmd:.     eststo pwr`=100*`tau''}{p_end}
{phang2}{cmd:. {c )-}}{p_end}
{phang2}{cmd:. esttab pwr10 pwr25 pwr50 pwr75 pwr90}{p_end}

{pstd}{ul:With bootstrap validation}{p_end}

{phang2}{cmd:. gepwreg lexp i.hhedlevel hhagey hhsize [fw=fw], per(0.5) boot(1000)}{p_end}
{phang2}{cmd:. gepwreg_setable}{p_end}

{pstd}{ul:Manual bandwidth}{p_end}

{phang2}{cmd:. gepwreg lexp hhagey hhsize [fw=fw], per(0.5) band(0.04)}{p_end}
{phang2}{cmd:. di "h = " e(h) "  N_eff = " e(N_eff)}{p_end}

{pstd}{ul:SE comparison table}{p_end}

{phang2}{cmd:. gepwreg lexp i.hhedlevel hhagey hhsize [fw=fw], per(0.5) boot(500)}{p_end}
{phang2}{cmd:. gepwreg_setable}{p_end}

{hline}
{marker references}{...}
{title:References}

{phang}
Araar, A. (2016). Percentile weights regression. {it:PEP Technical Note},
Université Laval.
{browse "http://dasp.ecn.ulaval.ca/stata_adds/gepwe/PEP_Notes_Araar_01.pdf"}

{phang}
Araar, A. (2023). Exploring heterogeneous effects: Quantile models and
percentile weights regression. {it:PEP Working Paper Series}, 2023-15.

{phang}
Deville, J.-C. (1999). Variance estimation for complex statistics and
estimators: Linearization and residual techniques.
{it:Survey Methodology}, 25(2):193-203.

{phang}
Firpo, S., Fortin, N. M., and Lemieux, T. (2009). Unconditional quantile
regressions. {it:Econometrica}, 77(3):953-973.

{phang}
Koenker, R. and Bassett, G. (1978). Regression quantiles.
{it:Econometrica}, 46(1):33-50.

{phang}
Newey, W. K. and McFadden, D. (1994). Large sample estimation and
hypothesis testing. In Engle, R. F. and McFadden, D. (Eds.),
{it:Handbook of Econometrics}, Vol. 4, pp. 2111-2245. Elsevier.

{phang}
Rios-Avila, F. (2020). Recentered influence functions (RIFs) in Stata:
RIF regression and RIF decomposition.
{it:The Stata Journal}, 20(1):51-94.

{phang}
Silverman, B. W. (1986).
{it:Density Estimation for Statistics and Data Analysis}.
Chapman & Hall, London.

{hline}
{title:Author}

{pstd}
Abdelkrim Araar{break}
Université Laval, Québec, Canada{break}
{browse "mailto:aabd@ecn.ulaval.ca":aabd@ecn.ulaval.ca}

{hline}
{title:Also see}

{psee}
{helpb qreg}: Quantile regression{p_end}
{psee}
{helpb rifreg}: RIF regression (UQR){p_end}
{psee}
{helpb gepwe}: Generate percentile kernel weights{p_end}
{psee}
{helpb regress}: OLS regression{p_end}
{hline}
