* p1_onestep_attenuation.do -- Paper (version 3), Section 3, Table "What the
* outcome-ranked estimator returns": the one-step weighted regression on
* outcome-ranked units (gepwreg ... rankdep, the estimator of version 1.3)
* against the true effect of the units at the tau-quantile of y, on four DGPs,
* with the MSE-optimal bandwidth of the paper; the bandwidth sweep; and the
* retained-variance ratio of Goldberger (1981).
*
* DGPs (n = 100000, x ~ N(0,1), e ~ N(0,1), U ~ U(0,1)):
*   A: y = 5 + 0.5 x + e                        effect 0.5
*   B: y = 5 + 0.5 x + exp( 0.5 x) e            effect 0.5 + 0.5 exp( 0.5x) e  (scale effect)
*   C: y = 5 + 0.5 x + exp(-0.5 x) e            effect 0.5 - 0.5 exp(-0.5x) e
*   D: y = 5 + (1 + U) x + U + e                effect 1 + U
* "true" = mean effect of the units whose rank in y is within +-0.05 of tau.
* Output: results/p1_onestep.csv, p1_sweep.csv, p1_goldberger.csv
* Run from the project root: do p1_onestep_attenuation

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
quietly regress yA x
local olsA = _b[x]

* ---- Table 1: one-step against the truth, h_MSE ---------------------------
tempname fh
file open `fh' using "results/p1_onestep.csv", write replace
file write `fh' "dgp,tau,true,onestep,se_IF,h,N_eff,ols" _n
foreach c in A B C D {
    quietly regress y`c' x
    local ols = _b[x]
    sort y`c'
    tempvar pc
    gen double `pc' = _n/_N
    di as text _n "DGP `c'   OLS = " as result %7.4f `ols'
    di as text %6s "tau" %9s "true" %10s "one-step" %9s "s.e." %9s "h" %8s "N_eff"
    foreach t in 0.1 0.25 0.5 0.75 0.9 {
        quietly summarize d`c' if abs(`pc' - `t') <= 0.05
        local tr = r(mean)
        quietly gepwreg y`c' x, per(`t') rankdep
        di as text %6.2f `t' as result %9.3f `tr' %10.4f _b[x] %9.4f _se[x] %9.4f e(h) %8.0f e(N_eff)
        file write `fh' "`c',`t',`tr',`=_b[x]',`=_se[x]',`=e(h)',`=e(N_eff)',`ols'" _n
    }
}
file close `fh'

* ---- Table 2: bandwidth sweep at tau = 0.5, DGPs A and D --------------------
file open `fh' using "results/p1_sweep.csv", write replace
file write `fh' "dgp,h,onestep,N_eff,true,ols" _n
foreach c in A D {
    sort y`c'
    tempvar pc
    gen double `pc' = _n/_N
    quietly summarize d`c' if abs(`pc' - 0.5) <= 0.05
    local tr = r(mean)
    quietly regress y`c' x
    local ols = _b[x]
    di as text _n "DGP `c', tau = 0.5, true = " as result %6.3f `tr' as text ", OLS = " as result %6.3f `ols'
    foreach h in 0.005 0.01 0.02 0.05 0.1 0.2 0.5 2 {
        quietly gepwreg y`c' x, per(0.5) rankdep band(`h')
        di as text "  h = " as result %5.3f `h' as text "   one-step = " as result %7.4f _b[x] as text "   N_eff = " as result %7.0f e(N_eff)
        file write `fh' "`c',`h',`=_b[x]',`=e(N_eff)',`tr',`ols'" _n
    }
}
file close `fh'

* ---- Table 3: Goldberger's retained-variance ratio, DGP A -------------------
* b_w / beta = [Var_w(y)/Var(y)] / [Var_w(x)/Var(x)]  (joint normality, constant beta)
file open `fh' using "results/p1_goldberger.csv", write replace
file write `fh' "h,N_eff,Vw_y_over_V_y,Vw_x_over_V_x,ratio,b_over_beta" _n
sort yA
tempvar pc w
gen double `pc' = _n/_N
gen double `w' = .
quietly summarize yA
local Vy = r(Var)
quietly summarize x
local Vx = r(Var)
di as text _n "DGP A, tau = 0.5: the retained-variance identity"
di as text %7s "h" %8s "N_eff" %12s "Vw(y)/V(y)" %12s "Vw(x)/V(x)" %9s "ratio" %10s "b/beta"
foreach h in 0.02 0.05 0.1 0.2 0.3 0.5 1 3 {
    quietly gepwreg yA x, per(0.5) rankdep band(`h')
    local b = _b[x]
    local ne = e(N_eff)
    quietly replace `w' = exp(-0.25*((`pc' - 0.5)/`h')^2)
    quietly summarize yA [aw=`w']
    local ry = r(Var)/`Vy'
    quietly summarize x [aw=`w']
    local rx = r(Var)/`Vx'
    di as text %7.3f `h' as result %8.0f `ne' %12.4f `ry' %12.4f `rx' %9.4f `=`ry'/`rx'' %10.4f `=`b'/0.5'
    file write `fh' "`h',`ne',`ry',`rx',`=`ry'/`rx'',`=`b'/0.5'" _n
}
file close `fh'
di as text _n "written: results/p1_onestep.csv, p1_sweep.csv, p1_goldberger.csv"
