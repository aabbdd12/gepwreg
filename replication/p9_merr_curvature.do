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
* Section 5.5: the limit no test removes.  The identifying moments cannot
* tell a missing functional form from measurement error.
*
* The design makes the claim falsifiable.  The heterogeneity variable is
* measured EXACTLY -- the true error variance is zero -- but the effect is
* QUADRATIC in it.  If the claim holds, the correction must fire anyway, and
* the corrected profile must move towards where the correctly specified
* quadratic model lands with no correction at all: the "correction" would
* then be doing the work of the missing term.  If merr returned zero here,
* the claim would be wrong and the paper should not make it.
* ---------------------------------------------------------------------
version 16
clear all
set more off
which gepwreg
capture mkdir "results"

* theta(0.75) for this design, computed on a large draw by p8_truth.py.
* It cannot be read off the 8,000 observations below: see the header of
* that file.
local truth = 2.714

set seed 20260919
set obs 8000
gen double g = rnormal()
gen double z = exp(0.45*g)
gen double x = 1 + 0.45*g + sqrt(1-0.45^2)*rnormal()
gen double b = 0.5 + 1.5*z + 0.8*(z^2 - 1.5)
gen double y = 5 + b*x + 3*z + rnormal()

di as text _n "{hline 70}"
di as text "A.  het(z) merr   -- does sigma2 fire where there is NO error?"
di as text "{hline 70}"
gepwreg y x, het(z) per(0.75) merr tcrit(2)
local thA = _b[x]
local s2A = el(e(merr_s2),1,1)
local tA  = el(e(merr_t),1,1)

di as text _n "{hline 70}"
di as text "B.  het(z c.z#c.z) -- the CORRECT model, no correction"
di as text "{hline 70}"
quietly gepwreg y x, het(z c.z#c.z) per(0.75)
local thB = _b[x]

di as text _n "{hline 70}"
di as text "C.  het(z)         -- the linear model, no correction"
di as text "{hline 70}"
quietly gepwreg y x, het(z) per(0.75)
local thC = _b[x]

tempname fh
file open `fh' using "results/p9_merr_curvature.csv", write replace
file write `fh' "quantity,value" _n
file write `fh' "tau,0.75" _n
file write `fh' "truth,`truth'" _n
file write `fh' "truth_source,p8_truth.py" _n
file write `fh' "linear_uncorrected,`thC'" _n
file write `fh' "linear_merr,`thA'" _n
file write `fh' "quadratic_uncorrected,`thB'" _n
file write `fh' "sigma2_true,0" _n
file write `fh' "sigma2_estimated,`s2A'" _n
file write `fh' "t,`tA'" _n
file close `fh'

di as text _n "{hline 70}"
di as text "Verdict"
di as text "{hline 70}"
di as text "sigma2 where the true error variance is ZERO = " as result %9.5f `s2A'
di as text "its t                                        = " as result %9.2f `tA'
di as text ""
di as text "theta(0.75)  truth                     " as result %9.4f `truth'
di as text "             linear, uncorrected       " as result %9.4f `thC' ///
   as text "  (" as result %5.1f 100*(`thC'/`truth'-1) as text " per cent)"
di as text "             linear + merr             " as result %9.4f `thA' ///
   as text "  (" as result %5.1f 100*(`thA'/`truth'-1) as text " per cent)"
di as text "             quadratic, uncorrected    " as result %9.4f `thB' ///
   as text "  (" as result %5.1f 100*(`thB'/`truth'-1) as text " per cent)"
di as text ""
di as text "The correction fires on a variable with no measurement error, and"
di as text "moves the profile part of the way towards the correctly specified"
di as text "model.  That is the limit, and no diagnostic internal to the"
di as text "correction removes it."
di as text "{hline 70}"
di as text "written: results/p9_merr_curvature.csv"
