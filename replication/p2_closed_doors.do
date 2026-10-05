* p2_closed_doors.do -- Paper (version 3), Section 3, "Why the outcome-ranked
* estimator cannot be repaired": four candidate corrections of the one-step
* weighted regression, each derived to fail, each shown to fail.
*   (i)   the Goldberger correction  b / {[Var_w(y)/Var(y)]/[Var_w(x)/Var(x)]}
*   (ii)  truncated regression on the band 0.35-0.45 of y (truncreg)
*   (iii) double weighting  w = K(rank y) / E[K | x]
*   (iv)  for contrast, selection on x instead of y (harmless)
* DGPs as in p1 (n = 100000; e ~ N(0,1)): A homoskedastic, B and C
* heteroskedastic, D heterogeneous.  Output: results/p2_closed_doors.csv,
* p2_band.csv.   Run: do p2_closed_doors

* ---------------------------------------------------------------------
* Replication, Araar (2026), version 3.  Run this file from the folder
* that contains it; the results are written to the results/ folder
* beside it, which is created if it does not exist.
*
* It needs gepwreg 1.4, the version the paper was produced with:
*
*   net install gepwreg, from("https://raw.githubusercontent.com/aabbdd12/gepwreg/v14") replace
*   net get     gepwreg, from("https://raw.githubusercontent.com/aabbdd12/gepwreg/v14") replace
*
* The second line also brings bkf98I.dta, which the application uses.
* The which gepwreg below must report v1.4; if it reports anything else,
* the numbers will not be the ones in the paper.
* ---------------------------------------------------------------------
version 16
clear all
set more off
which gepwreg
capture mkdir "results"

set seed 20260916
set obs 100000
gen double x = rnormal()
gen double e = rnormal()
gen double U = runiform()
gen double yA = 5 + 0.5*x + e
gen double yB = 5 + 0.5*x + exp( 0.5*x)*e
gen double yC = 5 + 0.5*x + exp(-0.5*x)*e
gen double yD = 5 + (1+U)*x + U + e
gen double dA = 0.5
gen double dB = 0.5 + 0.5*exp( 0.5*x)*e
gen double dC = 0.5 - 0.5*exp(-0.5*x)*e
gen double dD = 1 + U
quietly summarize x
local Vx = r(Var)

* ---- (i) Goldberger correction and (iii) double weighting, all DGPs ---------
tempname fh
file open `fh' using "results/p2_closed_doors.csv", write replace
file write `fh' "dgp,tau,true,onestep,ratio,corrected,double_weighted" _n
foreach c in A B C D {
    quietly summarize y`c'
    local Vy = r(Var)
    sort y`c'
    tempvar pc w K m w2
    gen double `pc' = _n/_N
    gen double `w'  = .
    gen double `w2' = .
    di as text _n "DGP `c'"
    di as text %6s "tau" %8s "true" %10s "one-step" %9s "ratio" %11s "corrected" %14s "double weight"
    foreach t in 0.1 0.25 0.5 0.75 0.9 {
        quietly summarize d`c' if abs(`pc' - `t') <= 0.05
        local tr = r(mean)
        quietly gepwreg y`c' x, per(`t') rankdep
        local b = _b[x]
        local h = e(h)
        quietly replace `w' = exp(-0.25*((`pc' - `t')/`h')^2)
        quietly summarize y`c' [aw=`w']
        local ry = r(Var)/`Vy'
        quietly summarize x [aw=`w']
        local rx = r(Var)/`Vx'
        local ratio = `ry'/`rx'
        * double weighting: divide the kernel weight by its mean within 100 bins of x
        tempvar bin Kbar
        quietly xtile `bin' = x, nq(100)
        quietly bysort `bin': egen double `Kbar' = mean(`w')
        quietly replace `w2' = `w'/`Kbar'
        quietly regress y`c' x [aw=`w2']
        local bdw = _b[x]
        drop `bin' `Kbar'
        di as text %6.2f `t' as result %8.3f `tr' %10.4f `b' %9.4f `ratio' %11.3f `=`b'/`ratio'' %14.4f `bdw'
        file write `fh' "`c',`t',`tr',`b',`ratio',`=`b'/`ratio'',`bdw'" _n
    }
}
file close `fh'

* ---- (ii) the hard band 0.35-0.45 of y, DGP A: OLS, truncreg, and selection on x
file open `fh' using "results/p2_band.csv", write replace
file write `fh' "estimator,slope,se" _n
quietly regress yA x
di as text _n "DGP A (beta = 0.5): OLS full sample = " as result %7.4f _b[x]
file write `fh' "OLS full sample,`=_b[x]',`=_se[x]'" _n
_pctile yA, p(35 45)
local ql = r(r1)
local qh = r(r2)
quietly regress yA x if yA >= `ql' & yA <= `qh'
di as text "OLS in the band 0.35-0.45 of y:      " as result %7.4f _b[x] as text "  (s.e. " %6.4f _se[x] ")"
file write `fh' "OLS in the band of y,`=_b[x]',`=_se[x]'" _n
quietly correlate x e if yA >= `ql' & yA <= `qh'
di as text "corr(x, e) in the band:              " as result %7.3f r(rho)
quietly correlate x e
di as text "corr(x, e) full sample:              " as result %7.3f r(rho)
capture noisily truncreg yA x if yA >= `ql' & yA <= `qh', ll(`ql') ul(`qh')
di as text "truncated regression on the band:    " as result %7.4f _b[x] as text "  (s.e. " %6.4f _se[x] ")"
file write `fh' "truncreg on the band of y,`=_b[x]',`=_se[x]'" _n
_pctile x, p(35 45)
local xl = r(r1)
local xh = r(r2)
quietly regress yA x if x >= `xl' & x <= `xh'
di as text "OLS in the band 0.35-0.45 of x:      " as result %7.4f _b[x] as text "  (s.e. " %6.4f _se[x] ")"
file write `fh' "OLS in the band 0.35-0.45 of x,`=_b[x]',`=_se[x]'" _n
_pctile x, p(25 75)
local xl = r(r1)
local xh = r(r2)
quietly regress yA x if x >= `xl' & x <= `xh'
di as text "OLS in the band 0.25-0.75 of x:      " as result %7.4f _b[x] as text "  (s.e. " %6.4f _se[x] ")"
file write `fh' "OLS in the band 0.25-0.75 of x,`=_b[x]',`=_se[x]'" _n
file close `fh'
di as text _n "written: results/p2_closed_doors.csv, p2_band.csv"
