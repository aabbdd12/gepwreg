{smcl}
{* *! gepwreg.sthlp  v1.4  19sep2026  Araar A.}{...}
{vieweralsosee "qreg" "help qreg"}{...}
{vieweralsosee "rifhdreg" "help rifhdreg"}{...}
{vieweralsosee "gepwe" "help gepwe"}{...}
{viewerjumpto "Syntax" "gepwreg##syntax"}{...}
{viewerjumpto "Description" "gepwreg##description"}{...}
{viewerjumpto "Objects and estimators" "gepwreg##methods"}{...}
{viewerjumpto "Options" "gepwreg##options"}{...}
{viewerjumpto "Standard errors" "gepwreg##se"}{...}
{viewerjumpto "Bandwidth" "gepwreg##bandwidth"}{...}
{viewerjumpto "Stored results" "gepwreg##results"}{...}
{viewerjumpto "Examples" "gepwreg##examples"}{...}
{viewerjumpto "References" "gepwreg##references"}{...}
{hline}
{title:Title}

{phang}
{bf:gepwreg} {hline 2} Percentile weights regression: the effect of the
covariates for the units at a quantile of the outcome, along the rank of a
covariate, or at the initial quantile, with analytical standard errors and
survey-design inference

{hline}
{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:gepwreg}
{depvar}
{indepvars}
{ifin}
{weight}
{cmd:,}
[{it:options}]

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt per:centile(#)}}target percentile {it:tau}; default {cmd:percentile(0.5)}{p_end}
{synopt:{opt het(varlist)}}heterogeneity model: interactions of {it:indepvars} with {it:varlist}{p_end}
{synopt:{opt het(qr)}}heterogeneity model: conditional quantile regression (the default){p_end}
{synopt:{opt init:ial(varlist|_all)}}effect at the {it:tau}-quantile of the outcome before the contribution of the listed regressors{p_end}
{synopt:{opt xref(numlist)}}reference values of {it:indepvars} for {opt initial()}; default zeros{p_end}
{synopt:{opt qgrid(numlist)}}quantiles of the step-1 quantile regressions; default 0.05(0.05)0.95{p_end}
{synopt:{opt rankvar(varname)}}profile of the effect along the rank of {it:varname}; one ordered variable, not a factor-variable term{p_end}
{synopt:{opt rankdep}}the one-step weighted regression on outcome-ranked units (version 1.3){p_end}

{syntab:Measurement error in the heterogeneity variables}
{synopt:{opt merr}}correct {opt het()} for classical measurement error in the heterogeneity variables{p_end}
{synopt:{opt tcrit(#)}}activation threshold of the correction; default {cmd:tcrit(2)}{p_end}
{syntab:Bandwidth}
{synopt:{it:(default)}}MSE-optimal plug-in bandwidth, every estimator{p_end}
{synopt:{opt optbw}}the plug-in asked for explicitly (same as the default){p_end}
{synopt:{opt sil:verman}}Silverman rule instead of the plug-in{p_end}
{synopt:{opt cband(#)}}Silverman constant; default {cmd:cband(0.9)}{p_end}
{synopt:{opt band(#)}}fixed bandwidth on the percentile scale{p_end}
{syntab:Standard errors}
{synopt:{it:(automatic)}}Taylor SE if {cmd:svyset} declares a PSU or strata; influence-function SE otherwise{p_end}
{synopt:{opt vce(svy)}}Taylor SE explicitly{p_end}
{synopt:{opt vce(if)}}influence-function SE even on {cmd:svyset} data{p_end}
{synopt:{opt boot(#)}}pairs bootstrap (PSUs within strata under a design); default {cmd:boot(0)}{p_end}
{synopt:{opt seed(#)}}seed of the bootstrap; default {cmd:seed(12345)}{p_end}
{synopt:{opt nodots}}suppress the progress dots of the bootstrap{p_end}
{synopt:{opt set:able}}display the standard errors side by side after the table (same as {cmd:gepwreg_setable}){p_end}
{syntab:Display}
{synopt:{opt l:evel(#)}}confidence level; default {cmd:level(95)}{p_end}
{synopt:{opt noc:onstant}}suppress the constant ({opt rankvar}/{opt rankdep} only){p_end}
{synoptline}

{phang}
{cmd:fweight}s, {cmd:aweight}s, and {cmd:pweight}s are allowed; see
{help weight}.  One rule decides which weights are used, and the header of
the output names the source in every run.  If {cmd:svyset} declares a weight
it is captured before anything is computed and used for everything -- the
step-1 coefficients, the bandwidth, the measurement-error correction, the
ranks and the variance -- and a weight expression on the command line is then
ignored, with a note.  Two declarations that can disagree would otherwise
estimate with one weight and build the design variance for the plan of
another.  The weight on the command line is used only when {cmd:svyset}
declares none; with neither, the estimation is unweighted.  {cmd:e(wtype)},
{cmd:e(wexp)} and {cmd:e(wsrc)} record what was used.{p_end}

{phang}
{it:indepvars} and the {opt het()} variables may contain factor variables;
see {help fvvarlist}.  Base levels are carried in {cmd:e(b)} as zero entries
and hidden in the table.  {opt rankvar()} is the exception: it names the
variable whose rank carries the profile, which needs a single ordered scale,
so a factor-variable term is refused there.{p_end}

{phang}
The post-estimation command

{p 8 17 2}
{cmd:gepwreg_setable} [{cmd:,} {opt naive}]

{phang}
displays the consistent standard errors available after the last
estimation side by side (influence function, Taylor under a design,
bootstrap) with the ratio of the two most informative ones.  The main
table already reports the appropriate one; the comparison is a check of
the asymptotic approximation, to be run when {cmd:e(N_eff)} is small, when
few PSUs fall in a stratum, or when a category is rare in the group: a
ratio far from 1 says to read the bootstrap.  {opt naive} adds the WLS
variance with the kernel weights taken as fixed, {cmd:e(V_naive)}
({opt rankvar}, {opt rankdep}), which is inconsistent and shown only for
the replication of the comparisons of the paper.

{hline}
{marker description}{...}
{title:Description}

{pstd}
{cmd:gepwreg} estimates how the effect of the covariates on an outcome varies
along a distribution, with the percentile weights of Araar (2016): Gaussian
kernel weights on the percentile scale, centred at a target percentile
{it:tau}.  Version 1.4 estimates three objects:

{phang2}
(a) the {bf:effect at the {it:tau}-quantile of the outcome}: the mean effect
of {it:indepvars} among the units whose outcome is at the {it:tau}-quantile
of {it:depvar} (the default, and {opt het()});{p_end}

{phang2}
(b) the {bf:effect at the initial quantile}: the mean effect among the units
whose outcome {it:before} the contribution of one or more regressors,
{it:y0 = y - b_i (x - xref)}, is at the {it:tau}-quantile of {it:y0}
({opt initial()});{p_end}

{phang2}
(c) the {bf:profile of the effect along a covariate}: the effect of
{it:indepvars} among the units at the {it:tau}-quantile of an observed
variable {it:z} ({opt rankvar(z)}).{p_end}

{pstd}
Objects (a) and (b) are estimated in two steps: a unit-level effect
{it:b_i} is estimated for every unit from a model of the heterogeneity of
the effect fitted on the whole sample, and the effects are averaged with
the percentile weights on the rank of {it:y} (or of {it:y0}).  Object (c)
is the kernel-weighted regression of {it:y} on {it:x} with the weights on
the rank of {it:z}.

{pstd}
{bf:What changed in 1.4.}  Up to version 1.3, object (a) was estimated
by the kernel-weighted regression of {it:y} on {it:x} with the weights on
the rank of {it:y} itself.  That regression selects the right units but
cannot measure the effect among them: inside a window on {it:y} the
outcome does not vary, so the weighted slope converges to a fraction of
the effect that shrinks with the bandwidth and vanishes as {it:h} goes to
zero (Goldberger 1981: selection on the dependent variable).  It equals
the effect only when {it:y} is a deterministic monotone function of a
single regressor, the case of the examples of Araar (2016, 2023).  On
survey data it is attenuated by one to two orders of magnitude, and its
sign can follow the heteroskedasticity of the error rather than the
effect.  It remains available with {opt rankdep}, labelled as a descriptive
slope.  The derivation, the simulations and the Burkina Faso application
are in Araar (2026, version 3).

{hline}
{marker methods}{...}
{title:Objects and estimators}

{pstd}
{ul:Notation.}  {it:p_i = F_n(y_i)} is the (weighted) percentile rank of
unit {it:i}; the percentile weight is

{pmore}
{it:w_i} = exp(-0.25*((p_i - tau)/h)^2) / (h * sqrt(2*pi) * n)

{pstd}
and {it:fw_i} is the sampling weight.  The group at {it:tau} is the set of
units with a large {it:w_i}; its effective size is
{it:N_eff} = (sum {it:fw_i w_i})^2 / sum ({it:fw_i w_i})^2.

{pstd}
{ul:Effect at the {it:tau}-quantile of {it:y} (default; het()).}  The
object is E[ dy/dx | y = q_tau ], the unconditional quantile partial
effect of Firpo, Fortin and Lemieux (2009), who show that it equals both
the mean derivative of the units at {it:q_tau} and the effect of a
marginal shift of {it:x} on the population quantile {it:q_tau}.  The
estimator is

{pmore}
theta(tau) = sum_i fw_i w_i b_i / sum_i fw_i w_i

{pstd}
where {it:b_i} is the effect of {it:x} for unit {it:i} from the
heterogeneity model of step 1.  The weights {it:w_i} on the rank of
{it:y} estimate the weight f(q_tau|x)/f(q_tau) of the FFL representation
without a conditional density.  Two heterogeneity models are offered:

{phang2}
{opt het(varlist)}: {it:y} is regressed on {it:x}, {it:z} and all the
interactions {it:x} x {it:z} on the whole sample, and
{it:b_ij} = {it:g_j} + sum_l {it:g_jl z_il}.  The estimator then has the
closed form theta_j = {it:g_j} + sum_l {it:g_jl} zbar_l(tau), where
zbar_l(tau) is the weighted mean of {it:z_l} in the group at {it:tau}:
the profile of the effect along {it:y} is the heterogeneity of the effect
times the composition of the group, and the table {cmd:e(decomp)} shows
that decomposition.  This model captures the heterogeneity that runs
through the observables {it:z}; heterogeneity ordered with {it:y} within
{it:z} is averaged out.{p_end}

{phang2}
{opt het(qr)} (the default): conditional quantile regressions of {it:y} on
{it:x} at the quantiles of {opt qgrid()}; the conditional rank {it:u_i} of
unit {it:i} is the midpoint of the two grid quantiles between which
{it:y_i} falls given {it:x_i}, and {it:b_i} is the quantile-regression
slope interpolated at {it:u_i}.  Under {it:y = h(x, e)}, {it:e}
independent of {it:x} and {it:h} increasing in {it:e}, {it:b_i} is the
unit's own derivative, and the estimator is the FFL representation of the
UQPE as the average of the conditional quantile partial effects.  This
model captures the heterogeneity ordered with {it:y} given {it:x},
observed or not; its price is the specification of the conditional
quantile function (linear in {it:x} at each quantile).{p_end}

{pstd}
The two models and the RIF regression of {helpb rifhdreg} target the same
object and can be compared: a difference between {opt het(z)} and
{opt het(qr)} signals heterogeneity ordered with {it:y} that {it:z} does
not carry.

{pstd}
{ul:Effect at the initial quantile (initial()).}  The object is
E[ b_i | y0 = q_tau(y0) ] with {it:y0_i = y_i - sum_j b_ij (x_ij - xref_j)},
the sum running over the regressors listed in {opt initial()} ({cmd:_all}
for all of them): the outcome net of their contribution at the reference
values {opt xref()}.  It answers "who benefits among the initially poor", the
question of programme evaluation (a binary {it:x} with {it:xref} = 0).
The unit-level effects of step 1 give {it:y0}, and the percentile weights
on the rank of {it:y0} give the estimator.  The object is identified by
the heterogeneity model (or by rank invariance, Heckman, Smith and Clements
1997).  {cmd:e(b_atq)} stores the effect at the quantile of {it:y} from the
same run.

{pstd}
{ul:Profile along a covariate (rankvar(z)).}  The object is
E[ b_i | rank(z) = tau ]: the effect of {it:x} for the units at the
{it:tau}-quantile of an observed {it:z} (age of the head, an asset index,
the predicted outcome {it:x'b}).  Because the weights do not involve
{it:y}, the kernel-weighted regression of {it:y} on {it:x} is a
varying-coefficient regression and estimates this object with the usual
smoothing bias; the influence-function, Taylor and bootstrap standard
errors of version 1.3 apply unchanged.

{pstd}
{ul:One-step weighted regression on outcome-ranked units (rankdep).}  The
weighted slope among units with {it:y} near {it:q_tau}.  Under joint
normality and a constant effect {it:b}, its limit is
{it:b} x [Var_w(y)/Var(y)] / [Var_w(x)/Var(x)] (Goldberger 1981), the
share of the variance of {it:y} that the kernel retains: about 1 percent
at the Silverman bandwidth.  In general its limit is of order {it:h}^2 and
carries the slope of the conditional density of {it:y} given {it:x} at
{it:q_tau}, not the effect.  No bandwidth, truncation correction,
variance correction or reweighting repairs it, because the information
about the effect is not inside a level set of {it:y}.

{hline}
{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt percentile(#)} specifies the target percentile {it:tau} in (0, 1).
The default is {cmd:percentile(0.5)}.

{phang}
{opt het(varlist)} selects the heterogeneity model on observables: the
effect of each regressor is allowed to vary linearly with each variable of
{it:varlist} (interactions).  {it:varlist} may include regressors of
{it:indepvars} and factor variables.  The unit-level effect is the
derivative of the step-1 model with respect to the regressor; when a
regressor is also a {opt het()} variable, that derivative includes its own
quadratic term (d/dx of {it:d x^2} = 2{it:d x}) and the interactions of
the other regressors with it, so that the effect of household size, say,
is the full marginal effect of size and not the coefficient of one
column.

{phang}
{opt het(qr)} selects the conditional quantile regression as the
heterogeneity model.  This is the default when neither {opt het()},
{opt rankvar()} nor {opt rankdep} is specified.  Standard errors require
{opt boot()}.

{phang}
{opt initial(varlist|_all)} requests the effect at the {it:tau}-quantile
of the outcome before the contribution of the listed regressors (a
treatment or programme variable, typically; {cmd:_all} nets out every
regressor), {it:y0}.  Standard errors by the pairs bootstrap;
{cmd:boot(200)} is assumed when {opt boot()} is not given, and {opt boot(#)}
sets the number of replications here as everywhere else -- there is one
option for the bootstrap, not one per estimator.  The replications show the
usual progress dots, with a count at the end of each row of fifty and an
{cmd:x} in place of a dot for a replication that failed; {opt nodots}
suppresses them.  The header reports the execution time in seconds, also in
{cmd:e(etime)}, which is the direct way to see what any option costs.  It is
measured with Stata's timer 99: a value kept there across a {cmd:gepwreg}
call is cleared.

{phang}
{opt xref(numlist)} gives the reference values of the regressors used to
define {it:y0}: one per term of the coefficient table in its order (the
values on base levels are ignored), or one per non-base regressor; the
default is zero for all.

{phang}
{opt qgrid(numlist)} sets the quantiles of the step-1 quantile regressions
of {opt het(qr)}; the default is 0.05, 0.10, ..., 0.95.

{phang}
{opt rankvar(varname)} estimates the profile of the effect along the rank
of {it:varname} by the kernel-weighted regression; {it:varname} must not be
the outcome.

{phang}
{opt rankdep} estimates the one-step weighted regression on outcome-ranked
units of version 1.3.  Its slope is a descriptive quantity, not the
effect; see {help gepwreg##methods:Objects and estimators}.

{dlgtab:Measurement error}

{phang}
{opt merr} corrects the step-1 fit of {opt het()} for classical measurement
error in the heterogeneity variables.  It is available only with
{opt het(}{it:varlist}{opt )}: the correction is a subtraction on the matrix of
second moments, which least squares provides and the check function of
{opt het(qr)} does not.

{pmore}
Why it matters.  A noisy heterogeneity variable does not simply flatten the
profile.  The composition term is untouched, because the measurement error is
independent of the rank, so the whole bias is a level shift plus a shrinkage of
the slope -- and the level shift carries the MAIN effect of the heterogeneity
variable, not its interaction.  The reported effect of {it:indepvars} can
therefore be over-stated by a factor of two to six while the interaction
coefficient looks almost right.

{pmore}
How it works.  The error variance is estimated from the sample itself, from the
third moments of the step-1 residual, one variance per heterogeneity variable.
Identification comes from the skewness of the heterogeneity variable, or from
its correlation with the regressors; a symmetric variable uncorrelated with
{it:indepvars} carries no information about the error, and the correction then
abstains.  This is the identification of Ben-Moshe, D'Haultfoeuille and Lewbel
(2017), specialised to the step-1 design.

{pmore}
When it activates.  For each heterogeneity variable the command forms
{it:t} = sigma2 / se(sigma2) and corrects only when {it:t} exceeds
{opt tcrit()}, whose default is 2.  Where the correction is not identified,
where the maintained assumptions fail, or where there is no measurement error
at all, the estimated variance goes to zero and the estimator falls back to the
uncorrected {opt het()} fit.  The output names the decision for every variable,
and {cmd:e(merr_s2)}, {cmd:e(merr_t)} and {cmd:e(merr_keep)} hold it.

{pmore}
What it costs.  The correction removes a bias and pays for it in variance.  In
the example below, on a heterogeneity variable whose reliability is 0.52, it
takes an over-statement of 16 per cent down to 0.7 per cent.  In repeated draws
of that design the uncorrected profile was wrong by 14 to 17 per cent every
time while the corrected one moved over -3.5 to +2 per cent; as the reliability
falls the correction stays centred on the truth and its dispersion grows -- at
a reliability of 0.25 the estimated error variance is still right, and the
corrected profile ranges over -11 to +12 per cent.  Report the corrected and
the uncorrected profile together, not the corrected one alone.

{pmore}
What it requires.  The measurement error must be classical: independent of the
regressors, of the true value and of the equation error, with mean zero.  No
distributional assumption is needed for the correction to be VALID, but the
precision of the estimated variance depends on the shape of the heterogeneity
variable, because the identifying moments are third moments.  A skewed variable
correlated with the regressors is the favourable case; a symmetric variable
independent of them carries almost no information about its own error, however
large that error is, and the command abstains rather than guess -- {it:z2} in
the example has a reliability of 0.57 and is left alone.

{pmore}
What the output reports.  Beside {it:sigma2} and its {it:t}, the table gives
for each heterogeneity variable the implied reliability,
1 - {it:sigma2}/Var({it:z~}), the skewness of {it:z~}, and the R-squared of
{it:z~} on the regressors: the two conditions above, one per column, with no
threshold applied to either.  They separate the two reasons a variable can be
left alone.  One that is skewed or correlated with the regressors but whose
{it:t} is small is identifiable and merely not resolved by this sample, and
more data would help.  One whose skewness and R-squared are both near zero
carries no information about its own error, and no sample size ever will.
The reliability says how much work the correction is doing, which is how much
variance it costs.  {cmd:e(merr_diag)} holds the three columns.

{pmore}
What it refuses, and why.  A heterogeneity variable taking two values, an
indicator of a factor variable included.  If {it:z} takes only {it:a} and
{it:b} then {it:z}^2 = ({it:a}+{it:b}) {it:z} - {it:ab}, so the third-moment
conditions are linear combinations of moments that step 1 has already set to
zero: the error variance is not identified at all, and a {it:t} computed from
it is a ratio of two numerical residuals.  Error in a categorical variable is
misclassification in any case, mechanically correlated with the true value,
which is not what is corrected here; the misclassification literature is the
place to look.  The command also refuses a singular step-1 design -- putting
{cmd:i.z} among the regressors while {it:z} is in {opt het()} makes their
interaction proportional to the indicators -- and a corrected moment matrix
that is not positive definite, which is what an error variance too large for
the design produces.  Measurement error in {it:indepvars} themselves is not
covered.  If a heterogeneity variable is also one of the {it:indepvars}, the
command refuses: the error would reach the regressor block as well, and that is
the case where the effect of interest is not defined.

{pmore}
The limit that no check can remove, and the comparison that reads it.  The
identifying moments cannot tell a missing functional form from measurement
error: if the effect is not linear in the heterogeneity variable, the curvature
lands in the step-1 residual and the third moments read it as an error
variance.  On a variable measured EXACTLY whose effect is quadratic, the
command returns an error variance with {it:t} = 11 and moves the profile part
of the way towards what the correctly specified model gives with no correction
at all.

{pmore}
Separating the two cases costs one command.  Fit the enriched model WITHOUT
{opt merr} -- a square, or {cmd:i.} for a categorical code -- and see where it
puts the profile.  If it lands where {opt merr} put it, {opt merr} was doing
the work of the missing term.  If it stays where the uncorrected linear model
was, the correction is doing something no functional form will do.  The two
cases are mirror images: in the quadratic design the corrected profile sits 4
per cent from the enriched one and the uncorrected profile 11 per cent from it;
in the design with real measurement error the distances are the other way
round, 12 and 3.  Do not enrich the model and keep {opt merr} at the same time:
the square of a noisy variable does not carry classical error -- if
{it:z~} = {it:z} + {it:u} then {it:z~}^2 = {it:z}^2 + 2{it:zu} + {it:u}^2,
whose error is correlated with the true value -- and the command will refuse
it.  Entering a categorical code as though it were continuous is the usual way
into this trap, which is why such a variable belongs in {cmd:i.} form, where
the correction is refused outright.

{pmore}
Inference.  The influence function accounts for the estimation of the error
variance, so the reported standard errors and the survey-design standard errors
remain valid in the regime where the correction is active.  At the activation
threshold itself the estimator is a pre-test estimator and no standard error is
valid uniformly there; where {it:t} is close to {opt tcrit()}, report the
corrected and uncorrected profiles side by side.  {opt boot(#)} re-runs the
whole procedure, the correction and its threshold included, in every resample.

{dlgtab:Bandwidth}

{phang}
{opt cband(#)}, {opt band(#)}, {opt silverman}, {opt optbw}: see
{help gepwreg##bandwidth:Bandwidth}.

{dlgtab:Standard errors}

{phang}
{opt vce(svy)} requests the Taylor linearisation under the {helpb svyset}
design; it is automatic when {cmd:svyset} declares a PSU or strata.
{opt vce(if)} keeps the influence-function standard errors on {cmd:svyset}
data.

{phang}
{opt boot(#)} requests the bootstrap with {it:#} replications: pairs of
observations, or PSUs within strata when a design is declared.  The whole
procedure (ranks, step 1, step 2) is re-run at each replication.

{phang}
{opt seed(#)} sets the seed of the bootstrap; the random-number state of
the session is restored afterwards.

{phang}
{opt setable} displays, after the coefficient table, the consistent
standard errors available for the estimation side by side, as the
post-estimation command {cmd:gepwreg_setable} does (a check, not a
choice); the command remains available on its own, for instance after
{cmd:estimates restore}.

{phang}
{opt level(#)} sets the confidence level of the table.

{phang}
{opt noconstant} suppresses the constant of the weighted regression
({opt rankvar}, {opt rankdep}); step 1 of the two-step keeps its constant.

{hline}
{marker se}{...}
{title:Standard errors}

{pstd}
{ul:Two-step, het(z).}  The influence function of theta_j is

{pmore}
psi_i(theta_j) = psi_i(g_j) + sum_l [ psi_i(g_jl) zbar_l + g_jl psi_i(zbar_l) ]

{pstd}
with psi(g) the OLS influence of step 1 and psi(zbar_l) the influence of
the weighted mean of {it:z_l} in the group, which has a direct term,
{it:n fw_i w_i (z_il - zbar_l) / S}, and a rank term from the estimated
percentile ranks inside the weights, computed by a reverse cumulative sum
over the units above {it:i} as in Deville (1999).  {cmd:e(V_IF)} is
psi'psi / n^2 with centred scores; under a survey design {cmd:e(V)} is the
Taylor variance of the same scores over PSUs within strata.  The three
components can be read from a bootstrap comparison: step-1 error alone
(delta method), plus the composition term, plus the rank term.

{pstd}
{ul:Ties in the ranking variable.}  Every observation of a tie group receives
the same estimated rank, so the kernel weights and the point estimates do not
depend on the order of the observations.  The rank term of the influence
function is summed over the SET of units whose rank is at least that of
{it:i}, not over the rows below it in the sorted data, so the standard errors
do not either.  Up to version 1.3 that sum was positional and the reported
standard errors moved when the data were permuted -- by 8e-04 on a ranking
variable with seven distinct values over 8,478 observations, by 7e-07 on a
continuous one.  The estimates themselves were never affected.

{pstd}
{ul:Two-step, het(qr) and initial().}  Pairs bootstrap ({opt boot()}),
PSUs within strata under a design.

{pstd}
{ul:rankvar and rankdep.}  The influence-function variance of versions
1.3 (Deville 1999), with the scores now centred; the naive WLS
variance is stored in {cmd:e(V_naive)} for reference.

{hline}
{marker bandwidth}{...}
{title:Bandwidth}

{pstd}
The kernel is on the percentile scale, and since 1.4 the MSE-optimal
plug-in is the default for every estimator.  {opt silverman} reverts to
{it:h} = {it:cband} x min(sd(p), IQR(p)/1.34) x {it:n}^(-1/5) and
{opt band(#)} fixes {it:h} by hand.

{pstd}
In the two-step, step 2 is a Nadaraya-Watson average of the unit-level
effects along the rank of {it:y}.  That rank is uniform by construction, so
the design-density term of the local constant bias vanishes and
MSE_j({it:h}) = {it:h}^4 theta_j''(tau)^2 + s_j^2 / (2 sqrt(2 pi) {it:n h}),
with s_j^2 = Var(b_ij | p_i = tau).  The coefficients do not share units, so
what is minimised is the average {it:relative} MSE and one bandwidth serves
the whole call:

{pmore}
{it:h*} = [ {it:k} / (8 sqrt(2 pi) {it:n} sum_j (theta_j'' / s_j)^2 ) ]^(1/5),

{pstd}
clipped to [0.5, 3] x {it:h_Silverman}.  The pilot is at
2 x {it:h_Silverman}, with a second difference of step
min(0.05, {it:tau}/2, (1-{it:tau})/2).  It is cheap: the unit-level effects
do not depend on {it:h}, so the pilot costs three weighted means and no
refitting.

{pstd}
For {opt rankvar} and {opt rankdep} the plug-in of Araar (2026) applies to
the local mean of the ranking variable: pilot 2 x {it:h_Silverman},
curvature at {it:tau} +/- 0.02, {it:h*} clipped to
[0.3, 5] x {it:h_Silverman}.

{pstd}
On simulations at {it:n} = 2,000 to 20,000 the two-step plug-in sits on
average 13 percent above the root mean squared error of the best possible
bandwidth, against 33 percent for the Silverman rule, with a worst case of
52 percent against 216.  It is least reliable at {it:tau} close to 0 or 1
when the profile of the effect is strongly curved there, because the pilot's
second difference then reads the profile where the kernel is cut by the
boundary of the rank; {opt silverman} is the fallback in that case.
{cmd:e(N_eff)} is the effective size of the group; values below 50 suggest a
bandwidth too narrow for the sample.

{pstd}
{bf:One warning about the standard errors.}  The influence-function variance
is derived with the bandwidth held fixed, and a bandwidth chosen on the
sample adds a term it cannot contain: the estimate inherits the variability
of the choice.  The term is negligible wherever the bandwidth does not matter
much, which is most of the range, and material where the MSE curve is steep
-- at {it:tau} near 0 or 1 with a curved profile.  On simulations at
{it:n} = 4,000 it is worth 23 percent of the standard deviation at
{it:tau} = 0.9, 4 percent at 0.5 and nothing at 0.1, and the coverage of the
nominal 95 percent interval falls to 0.82 in the first case.  {opt boot(#)}
re-chooses the bandwidth in every resample and prices it correctly (ratio to
the Monte Carlo dispersion 1.05, 0.99, 1.03 at the three quantiles), so use
it near the tails or whenever {cmd:e(h)} comes back near 3 x {it:h_Silverman}
or 0.5 x {it:h_Silverman}.  {opt silverman} and {opt band(#)} make the
bandwidth deterministic and restore the analytical variance exactly.

{hline}
{marker results}{...}
{title:Stored results}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(tau)}}target percentile{p_end}
{synopt:{cmd:e(h)}}bandwidth on the rank of the ranking variable{p_end}
{synopt:{cmd:e(h0)}}bandwidth on the rank of {it:y0} ({opt initial()}){p_end}
{synopt:{cmd:e(N_eff)}}effective size of the group at {it:tau}{p_end}
{synopt:{cmd:e(N_eff0)}}effective size on the rank of {it:y0} ({opt initial()}){p_end}
{synopt:{cmd:e(y0_regressors)}}indicator of the regressors netted out of {it:y0}{p_end}
{synopt:{cmd:e(boot)}}bootstrap replications{p_end}
{synopt:{cmd:e(do_svy)}}1 if the Taylor variance was used{p_end}
{synopt:{cmd:e(etime)}}execution time in seconds, as shown in the header{p_end}
{synopt:{cmd:e(merr)}}1 if {opt merr} was requested{p_end}
{synopt:{cmd:e(tcrit)}}activation threshold of the correction{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}the estimated effects{p_end}
{synopt:{cmd:e(V)}}their variance: Taylor under a design, influence function otherwise, bootstrap for {opt het(qr)} and {opt initial()}{p_end}
{synopt:{cmd:e(V_IF)}}influence-function variance ({opt het(z)}, {opt rankvar}, {opt rankdep}){p_end}
{synopt:{cmd:e(V_svy)}}Taylor variance (when used){p_end}
{synopt:{cmd:e(V_boot)}}bootstrap variance (when {opt boot()} > 0){p_end}
{synopt:{cmd:e(V_naive)}}naive WLS variance ({opt rankvar}, {opt rankdep}){p_end}
{synopt:{cmd:e(b_atq)}}effect at the quantile of {it:y} in an {opt initial()} run{p_end}
{synopt:{cmd:e(b_step1)}}step-1 coefficients of {opt het(z)}{p_end}
{synopt:{cmd:e(zbar_tau)}}means of the {opt het()} variables in the group at {it:tau} (at the {it:tau}-quantile of y0 with {opt initial()}){p_end}
{synopt:{cmd:e(zbar_pop)}}their population means{p_end}
{synopt:{cmd:e(decomp)}}composition of the group and contribution to each effect (group at the {it:tau}-quantile of y0 with {opt initial()}){p_end}
{synopt:{cmd:e(b_qgrid)}}step-1 quantile-regression coefficients of {opt het(qr)}{p_end}
{synopt:{cmd:e(xref)}}reference values of {it:indepvars}{p_end}
{synopt:{cmd:e(merr_s2)}}estimated error variance of each {opt het()} variable ({opt merr}){p_end}
{synopt:{cmd:e(merr_t)}}its {it:t} statistic, {it:sigma2}/se({it:sigma2}) ({opt merr}){p_end}
{synopt:{cmd:e(merr_keep)}}1 where the correction was applied, 0 where the command abstained ({opt merr}){p_end}
{synopt:{cmd:e(merr_diag)}}implied reliability, skewness of {it:z~} and R2 of {it:z~} on the regressors, one row per {opt het()} variable ({opt merr}){p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:gepwreg}{p_end}
{synopt:{cmd:e(method)}}{cmd:twostep}, {cmd:initial}, {cmd:rankvar} or {cmd:rankdep}{p_end}
{synopt:{cmd:e(het)}}{opt het()} variables, or {cmd:qr}{p_end}
{synopt:{cmd:e(qgrid)}}quantile grid of {opt het(qr)}{p_end}
{synopt:{cmd:e(rankvar)}}ranking variable{p_end}
{synopt:{cmd:e(bw_method)}}bandwidth rule{p_end}
{synopt:{cmd:e(SE_type)}}standard errors reported{p_end}
{synopt:{cmd:e(wgt)}}weights used, as printed in the header{p_end}
{synopt:{cmd:e(wtype)}}{cmd:pweight}, {cmd:aweight} or {cmd:fweight}{p_end}
{synopt:{cmd:e(wexp)}}the weight expression{p_end}
{synopt:{cmd:e(wsrc)}}{cmd:svy} if the weights came from {cmd:svyset}, {cmd:cmd} if from the command line, {cmd:none} if unweighted{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}

{hline}
{marker examples}{...}
{title:Examples}

{pstd}
The examples use {cmd:bkf98I.dta}, the Burkina Faso 1998 household survey
extract distributed with the package ({cmd:net get gepwreg}): 8,478
households, 10 strata, 425 primary sampling units, sampling weight
{cmd:weight}.  The file carries its {cmd:svyset} settings, so the weights and
the design are taken automatically and no weight expression is written on the
command line; writing one would only produce a note saying it is ignored.  On
data that are not {cmd:svyset}, write {cmd:[pw=}{it:w}{cmd:]} as usual.{p_end}

{pstd}{ul:Setup}{p_end}

{phang2}{cmd:. use bkf98I, clear}{p_end}
{phang2}{cmd:. generate lexp  = ln(exppc)}{p_end}
{phang2}{cmd:. generate male  = (sex == 1)}{p_end}
{phang2}{cmd:. generate urban = (zone == 2)}{p_end}

{pstd}{ul:Effect of the covariates for the households at the first quartile of expenditure}{p_end}

{phang2}{cmd:. gepwreg lexp size male i.gse, per(0.25) het(urban size)}{p_end}
{phang2}{cmd:. gepwreg lexp size male i.gse, per(0.25) boot(200)}{p_end}
{phang2}{cmd:. rifhdreg lexp size male i.gse [pw=weight], rif(q(25))}{p_end}

{pstd}{ul:Survey design: Taylor standard errors, automatic once svyset}{p_end}

{phang2}{cmd:. svyset psu [pw=weight], strata(strata)}{p_end}
{phang2}{cmd:. gepwreg lexp size male i.gse, per(0.25) het(urban size)}{p_end}
{phang2}{cmd:. gepwreg lexp size male i.gse, per(0.25) het(urban size) boot(200)}{p_end}
{phang2}{cmd:. gepwreg_setable}{p_end}

{pstd}{ul:Profile across percentiles}{p_end}

{phang2}{cmd:. foreach tau in 0.10 0.25 0.50 0.75 0.90 {c -(}}{p_end}
{phang2}{cmd:.     gepwreg lexp size male i.gse, per(`tau') het(urban size)}{p_end}
{phang2}{cmd:.     estimates store q`=100*`tau''}{p_end}
{phang2}{cmd:. {c )-}}{p_end}
{phang2}{cmd:. estimates table q10 q25 q50 q75 q90, se}{p_end}

{pstd}{ul:Effect at the initial quantile (households poor before the urban premium)}{p_end}

{phang2}{cmd:. gepwreg lexp urban size male i.gse, per(0.25) het(size) initial(urban)}{p_end}

{pstd}{ul:Profile of the effect along household size}{p_end}

{phang2}{cmd:. gepwreg lexp male urban i.gse, per(0.25) rankvar(size)}{p_end}

{pstd}{ul:The one-step regression of version 1.3, for comparison}{p_end}

{phang2}{cmd:. gepwreg lexp size male i.gse, per(0.25) rankdep}{p_end}

{pstd}{ul:Measurement error in a heterogeneity variable}{p_end}

{pstd}
A simulation, because a correction for bias can only be shown against a truth
that is known.  {it:z1} is skewed and correlated with {it:x}; {it:z2} is
symmetric and independent of it.  Both are observed with error.{p_end}

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set seed 20260919}{p_end}
{phang2}{cmd:. set obs 8000}{p_end}
{phang2}{cmd:. generate double g   = rnormal()}{p_end}
{phang2}{cmd:. generate double z1  = exp(0.45*g)}{p_end}
{phang2}{cmd:. generate double z2  = runiform()}{p_end}
{phang2}{cmd:. generate double x   = 1 + 0.45*g + sqrt(1-0.45^2)*rnormal()}{p_end}
{phang2}{cmd:. generate double b   = 0.5 + 1.5*z1 + 1.0*z2}{p_end}
{phang2}{cmd:. generate double y   = 5 + b*x + 3*z1 + 2*z2 + rnormal()}{p_end}
{phang2}{cmd:. generate double zt1 = z1 + 0.50*rnormal()}{p_end}
{phang2}{cmd:. generate double zt2 = z2 + 0.25*rnormal()}{p_end}

{phang2}{cmd:. gepwreg y x, het(z1 z2)   per(0.75)}{p_end}
{phang2}{cmd:. gepwreg y x, het(zt1 zt2) per(0.75)}{p_end}
{phang2}{cmd:. gepwreg y x, het(zt1 zt2) per(0.75) merr}{p_end}

{pstd}
The true error variances are 0.25 and 0.0625, and the true effect at
{it:tau} = 0.75 is 2.976.  The first run is the oracle, on the heterogeneity
variables themselves, and returns 2.963.  The second, on the proxies, returns
3.453 -- over-stated by 16 per cent.  The third estimates 0.2491 for the first
error variance, with {it:t} = 19, and abstains on the second, whose {it:t} is
1.3: {it:z2} is symmetric and independent of {it:x}, so its error is not
identified however large it is, and here its reliability is 0.57.  The
corrected effect is 2.956, within 0.7 per cent of the truth.{p_end}

{pstd}
Before reporting any corrected profile, run the comparison of the previous
section: {cmd:gepwreg y x, het(zt1 c.zt1#c.zt1 zt2) per(0.75)}, without
{opt merr}.  Here it returns 3.345, close to the uncorrected 3.453 and far from
the corrected 2.956 -- enriching the model does not do what the correction
does, which is the sign that the correction is answering measurement error and
not a missing term.{p_end}

{hline}
{marker references}{...}
{title:References}

{phang}
Araar, A. (2016). Percentile weights regression. {it:PEP Technical Note},
Universite Laval.
{browse "http://dasp.ecn.ulaval.ca/stata_adds/gepwe/PEP_Notes_Araar_01.pdf"}

{phang}
Araar, A. (2023). Exploring heterogeneous effects: Quantile models and
percentile weights regression. {it:PEP Working Paper Series}, 2023-15.

{phang}
Araar, A. (2026). Percentile weights regression: what the outcome-ranked
estimator measures, a two-step estimator of the effect at a quantile, and
inference under complex survey design (version 3).  Zenodo,
{browse "https://doi.org/10.5281/zenodo.20315684":10.5281/zenodo.20315684}.
The reference for the estimators, their standard errors and the survey
design variance implemented in this version, and the paper to cite for the
command.

{phang}
Deville, J.-C. (1999). Variance estimation for complex statistics and
estimators: Linearization and residual techniques.
{it:Survey Methodology}, 25(2):193-203.

{phang}
Firpo, S., Fortin, N. M., and Lemieux, T. (2009). Unconditional quantile
regressions. {it:Econometrica}, 77(3):953-973.

{phang}
Goldberger, A. S. (1981). Linear regression after selection.
{it:Journal of Econometrics}, 15(3):357-366.

{phang}
Heckman, J. J., Smith, J., and Clements, N. (1997). Making the most out
of programme evaluations and social experiments: Accounting for
heterogeneity in programme impacts. {it:Review of Economic Studies},
64(4):487-535.

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
{title:Version history}

{phang}
1.4 (September 2026).  One release since 1.3; the version numbers in between
were development states and were never distributed.

{pmore}
{ul:The estimator.}  The effect at the {it:tau}-quantile of the outcome is
estimated by a two-step estimator ({opt het(qr)} by default,
{opt het(varlist)}): unit-level effects from a heterogeneity model fitted on
the whole sample, averaged with the percentile weights on the rank of the
outcome.  The kernel-weighted regression on outcome-ranked units of 1.3
converges to a bandwidth-dependent fraction of the effect as soon as the
outcome has an error term; it is kept as {opt rankdep} and labelled
descriptive.  New: {opt initial()} and {opt xref()} for the effect at the
initial quantile, {opt qgrid()}.  {opt rankvar()} is unchanged and now says
why it refuses several variables or a factor-variable term.

{pmore}
{ul:Bandwidth.}  The MSE-optimal plug-in is the default for every estimator,
the second step of the two-step included; {opt silverman} and {opt band()}
override it.

{pmore}
{ul:Measurement error.}  {opt merr} corrects the step-1 fit of
{opt het(}{it:varlist}{opt )} for classical measurement error, one variance
per variable, identified from the third moments of the step-1 residual and
applied only where {it:t} exceeds {opt tcrit()}; the influence function
accounts for the estimated variance, and {opt boot()} redoes the correction
and its threshold in every resample.  The command refuses what it cannot
identify: a heterogeneity variable with two values, a singular step-1 design,
a corrected moment matrix that is not positive definite.  The table reports
the implied reliability, the skewness and the R-squared on the regressors per
variable, in {cmd:e(merr_diag)}.

{pmore}
{ul:Standard errors.}  Taylor linearisation by default under {cmd:svyset}
({opt vce(if)} keeps the IF standard errors), {cmd:e(V_IF)} always stored,
and a bootstrap that re-runs the whole procedure.  The scores of
{opt rankvar} and {opt rankdep} are centred -- the constant of the indirect
term was omitted up to 1.3, a difference of second order.  Standard errors no
longer depend on the order of the observations: the rank term is summed over
the set of units above {it:i} rather than over the rows below it, which the
estimates never needed and the variances did.

{pmore}
{ul:Weights.}  One rule, stated in the header: {cmd:svyset} first, the command
line only when {cmd:svyset} declares none; {cmd:e(wtype)}, {cmd:e(wexp)} and
{cmd:e(wsrc)} record the choice.

{pmore}
{ul:Also.}  Factor variables: every non-base level gets its own indicator, where
up to 1.3 the second level was pooled with the base.  A PSU count per stratum.
{cmd:gepwreg_setable} in its own file and as the {opt setable} option.  A
rewritten dialog box covering every option.  Progress dots on every bootstrap
path, with {opt nodots}.  Execution time in the header and in {cmd:e(etime)}.

{phang}
1.3 (June 2026).  MSE-optimal bandwidth as the default, Taylor linearisation
under {cmd:svyset}, ties in the ranking variable given the same rank,
{cmd:rankvar()}.

{hline}
{title:Author}

{pstd}
Abdelkrim Araar{break}
Universite Laval and Partnership for Economic Policy (PEP), Quebec, Canada{break}
{browse "mailto:aabd@ecn.ulaval.ca":aabd@ecn.ulaval.ca}{break}
Package: {browse "https://github.com/aabbdd12/gepwreg"}

{hline}
{title:Also see}

{psee}
{helpb qreg}: Quantile regression{p_end}
{psee}
{helpb rifhdreg}: RIF regression (UQR), Rios-Avila (2020){p_end}
{psee}
{helpb gepwe}: Generate percentile kernel weights{p_end}
{hline}
