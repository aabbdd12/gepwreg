* ---------------------------------------------------------------------
* Replication, Araar (2026), version 3.  Run this file from the folder
* that contains it; the results are written to the results/ folder
* beside it, which is created if it does not exist.
*
* It needs gepwreg 1.4, the version the paper was produced with:
*
*   net install gepwreg, from("https://raw.githubusercontent.com/aabbdd12/gepwreg/v14") replace
*
* The which gepwreg below must report v1.4.
* ---------------------------------------------------------------------
*
* Section 5.3: measurement error in the heterogeneity variables, on a design
* where the truth is known.  Two heterogeneity variables: z1 skewed and
* correlated with the regressor, z2 symmetric and independent of it.  Both
* are observed with error, and only the first is identified -- which is the
* point of the design, since it exercises the t rule in both directions.
*
* True error variances: 0.25 on z1, 0.0625 on z2.
* ---------------------------------------------------------------------
version 16
clear all
set more off
which gepwreg
capture mkdir "results"

set seed 20260919
set obs 8000
generate double g   = rnormal()
generate double z1  = exp(0.45*g)
generate double z2  = runiform()
generate double x   = 1 + 0.45*g + sqrt(1-0.45^2)*rnormal()
generate double b   = 0.5 + 1.5*z1 + 1.0*z2
generate double y   = 5 + b*x + 3*z1 + 2*z2 + rnormal()
generate double zt1 = z1 + 0.50*rnormal()
generate double zt2 = z2 + 0.25*rnormal()

local t = 0.75

* ---- the truth ------------------------------------------------------------
* theta(0.75) = E[b | y = q_0.75] is a population quantity and cannot be read
* off these 8,000 observations: a window of half-width 0.02 on the rank holds
* 320 units, whose mean carries a Monte Carlo standard error of about 0.06 --
* larger than every difference below.  The value is computed once, on a large
* draw, by p8_truth.py, which also shows that it has stopped moving both in
* the sample size and in the window width.
local truth = 2.9759
di as text _n "truth at tau = `t' (from p8_truth.py): " as result %8.4f `truth'

* the small-sample window is reported only to show how far it can sit from
* the truth, which is why it is not used as the reference
sort y
gen double pc = _n/_N
quietly summarize b if abs(pc - `t') <= 0.02
di as text "  (same quantity from these 8,000 units: " as result %8.4f r(mean) ///
   as text ", on " as result r(N) as text " units)"

* ---- the reliabilities, in closed form ------------------------------------
* Var(z1) is that of a lognormal with sigma = 0.45; Var(z2) is 1/12.  These
* are population quantities, so they are computed rather than estimated; the
* command reports its own reliability only where it corrects.
local vz1 = (exp(0.45^2)-1)*exp(0.45^2)
local vz2 = 1/12
local rel1 = `vz1'/(`vz1' + 0.25)
local rel2 = `vz2'/(`vz2' + 0.0625)
di as text "reliability of zt1 = " as result %5.3f `rel1' ///
   as text ",  of zt2 = " as result %5.3f `rel2' ///
   as text "   (z2 is the one left uncorrected)"

* ---- the oracle: the heterogeneity variables themselves -------------------
quietly gepwreg y x, het(z1 z2) per(`t')
local b_oracle = _b[x]

* ---- the naive fit: the proxies, no correction ----------------------------
quietly gepwreg y x, het(zt1 zt2) per(`t')
local b_naive = _b[x]

* ---- the corrected fit ----------------------------------------------------
gepwreg y x, het(zt1 zt2) per(`t') merr
local b_merr = _b[x]
matrix S2 = e(merr_s2)
matrix TM = e(merr_t)
matrix KP = e(merr_keep)
local s2_1 = S2[1,1]
local s2_2 = S2[1,2]
local t_1  = TM[1,1]
local t_2  = TM[1,2]
local k_1  = KP[1,1]
local k_2  = KP[1,2]

* ---- the comparison that separates error from a missing term --------------
* The enriched model WITHOUT the correction.  If it lands where merr put the
* profile, merr was doing the work of the missing term; if it stays where the
* uncorrected linear model was, the correction is answering measurement
* error.  See Section 5.5.
quietly gepwreg y x, het(zt1 c.zt1#c.zt1 zt2) per(`t')
local b_rich = _b[x]

* ---- out ------------------------------------------------------------------
tempname fh
file open `fh' using "results/p8_merr.csv", write replace
file write `fh' "quantity,value" _n
file write `fh' "tau,`t'" _n
file write `fh' "truth,`truth'" _n
file write `fh' "truth_source,p8_truth.py" _n
file write `fh' "oracle_on_z,`b_oracle'" _n
file write `fh' "naive_on_proxies,`b_naive'" _n
file write `fh' "corrected,`b_merr'" _n
file write `fh' "enriched_no_merr,`b_rich'" _n
file write `fh' "sigma2_z1,`s2_1'" _n
file write `fh' "sigma2_z2,`s2_2'" _n
file write `fh' "t_z1,`t_1'" _n
file write `fh' "t_z2,`t_2'" _n
file write `fh' "keep_z1,`k_1'" _n
file write `fh' "keep_z2,`k_2'" _n
file write `fh' "reliability_z1,`rel1'" _n
file write `fh' "reliability_z2,`rel2'" _n
file close `fh'

di as text _n "{hline 70}"
di as text "Section 5.3, at tau = `t'"
di as text "{hline 70}"
di as text "truth                          = " as result %8.4f `truth'
di as text "oracle, on z1 and z2           = " as result %8.4f `b_oracle'
di as text "naive, on the proxies          = " as result %8.4f `b_naive' ///
   as text "   (" as result %5.1f 100*(`b_naive'/`truth'-1) as text " per cent)"
di as text "corrected                      = " as result %8.4f `b_merr' ///
   as text "   (" as result %5.1f 100*(`b_merr'/`truth'-1) as text " per cent)"
di as text "enriched model, no correction  = " as result %8.4f `b_rich'
di as text _n "error variances (true: 0.25 and 0.0625)"
di as text "  z1: sigma2 = " as result %7.4f `s2_1' as text ", t = " ///
   as result %6.2f `t_1' as text ", corrected = " as result `k_1'
di as text "  z2: sigma2 = " as result %7.4f `s2_2' as text ", t = " ///
   as result %6.2f `t_2' as text ", corrected = " as result `k_2'
di as text "{hline 70}"
di as text "written: results/p8_merr.csv"
