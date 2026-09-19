* p3_replications.do -- Paper (version 3), Section 3, "The examples of the
* earlier notes": (a) examples A1 and A3 of the 2016 note with the one-step
* estimator and the note's bandwidth (Silverman/3); (b) the Monte Carlo of the
* 2023 working paper (section 2.7, basic case) with the exact estimator, then
* with a larger error variance.
*
* (a) A1: x1 = _n-1, p = _n/1000, y = 1000 + p x1 + 1e-5 u.  y is a deterministic
*         monotone function of x1: the rank of y is the rank of x1, and the
*         one-step slope is the total derivative 2 tau (sd(x1 | y) = 0).
*     A3: x1 ~ U(0,1000), p = _n/1000 independent of x1, y = 2000 + p x1 + 60 u.
*         The note's coefficient is tau; ranking on y no longer ranks on x1.
* (b) N = 4000, x1 = 0.2 r^0.5 + N(1,1) with r the ascending rank of x1,
*     x2 ~ N(2,2), eps ~ N(0, sd), beta1(r) = (((r-1)^2 + (N-r)^2)/(10N))^0.5,
*     beta2 = 2, c = 10.  "true" = E[beta1 | rank y within +-0.05 of tau].
*     Estimators: one-step on y (Silverman, as in 2023), one-step ranked on x1
*     (rankvar), RIF-OLS (rifhdreg), OLS.  20 seeds, means and s.d.
* Output: results/p3_note2016.csv, p3_mc2023.csv
* Run: do p3_replications

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
which rifhdreg

* ---------------------------------------------------------------------------
* (a) the 2016 note
* ---------------------------------------------------------------------------
tempname fh
file open `fh' using "results/p3_note2016.csv", write replace
file write `fh' "example,tau,true,onestep_h3,onestep_hSil" _n

set seed 1234
set obs 1000
gen double x1 = _n - 1
gen double p  = _n/1000
gen double y  = 1000 + p*x1 + 0.00001*runiform()
sort y
gen double pc = _n/_N
quietly summarize pc, detail
local hs = 0.9*min(r(sd), (r(p75)-r(p25))/1.34)*_N^(-0.2)
di as text _n "2016 note, example A1 (deterministic): true 2 tau;  h_Sil = " as result %6.4f `hs'
foreach t in 0.05 0.15 0.25 0.35 0.45 0.55 0.65 0.75 0.85 0.95 {
    quietly gepwreg y x1, per(`t') rankdep band(`=`hs'/3')
    local b3 = _b[x1]
    quietly gepwreg y x1, per(`t') rankdep silverman
    di as text "tau = " as result %4.2f `t' as text "   true = " as result %5.2f 2*`t' ///
       as text "   one-step h/3 = " as result %6.3f `b3' as text "   h_Sil = " as result %6.3f _b[x1]
    file write `fh' "A1,`t',`=2*`t'',`b3',`=_b[x1]'" _n
}

clear
set seed 1234
set obs 1000
gen double x1 = runiform()*1000
gen double p  = _n/1000
gen double y  = 2000 + p*x1 + 60*runiform()
sort y
gen double pc = _n/_N
quietly summarize pc, detail
local hs = 0.9*min(r(sd), (r(p75)-r(p25))/1.34)*_N^(-0.2)
di as text _n "2016 note, example A3 (x1 independent of p): the note's coefficient is tau"
foreach t in 0.1 0.25 0.5 0.75 0.9 {
    quietly gepwreg y x1, per(`t') rankdep band(`=`hs'/3')
    local b3 = _b[x1]
    quietly gepwreg y x1, per(`t') rankdep silverman
    di as text "tau = " as result %4.2f `t' as text "   'true' = " as result %5.2f `t' ///
       as text "   one-step h/3 = " as result %6.3f `b3' as text "   h_Sil = " as result %6.3f _b[x1]
    file write `fh' "A3,`t',`t',`b3',`=_b[x1]'" _n
}
file close `fh'

* ---------------------------------------------------------------------------
* (b) the 2023 Monte Carlo, basic case, with sd(eps) = 1 (the paper), 10, 25
* ---------------------------------------------------------------------------
file open `fh' using "results/p3_mc2023.csv", write replace
file write `fh' "sd_eps,tau,R2,sdx1_given_y,true,onestep_mean,onestep_sd,rankvar_mean,rankvar_sd,rif_mean,rif_sd,ols_mean" _n
local N = 4000
local R = 20
foreach sig in 1 10 25 {
    di as text _n "2023 MC basic case, sd(eps) = `sig'   (`R' seeds; mean [sd])"
    di as text %5s "tau" %6s "R2" %10s "sd(x1|y)" %8s "true" %16s "one-step y" %16s "rankvar(x1)" %16s "RIF-OLS" %8s "OLS"
    foreach t in 0.05 0.15 0.25 0.35 0.45 0.55 0.65 0.75 0.85 0.95 {
        tempname M
        matrix `M' = J(`R', 7, .)
        forvalues s = 1/`R' {
            clear
            quietly set obs `N'
            set seed `=500 + `s''
            gen double z  = rnormal(1, 1)
            sort z
            gen double r  = _n
            gen double x1 = 0.2*sqrt(r) + z
            gen double x2 = rnormal(2, 2)
            gen double eps = rnormal(0, `sig')
            gen double b1 = sqrt(((r-1)^2 + (`N'-r)^2)/(10*`N'))
            gen double y  = 10 + b1*x1 + 2*x2 + eps
            quietly regress y x1 x2
            local r2  = e(r2)
            local ols = _b[x1]
            * sd(x1 | y) relative to sd(x1): residual s.d. within bins of 100 on y
            sort y
            gen double pc = _n/_N
            gen int ybin = ceil(_n/100)
            quietly bysort ybin: egen double x1m = mean(x1)
            quietly summarize x1
            local sdx = r(sd)
            gen double x1r = x1 - x1m
            quietly summarize x1r
            local sdx1y = r(sd)/`sdx'
            quietly summarize b1 if abs(pc - `t') <= 0.05
            local tr = r(mean)
            quietly gepwreg y x1 x2, per(`t') rankdep silverman
            local b_y = _b[x1]
            quietly gepwreg y x1 x2, per(`t') rankvar(x1) silverman
            local b_r = _b[x1]
            quietly rifhdreg y x1 x2, rif(q(`=100*`t''))
            local b_rif = _b[x1]
            matrix `M'[`s', 1] = `r2'
            matrix `M'[`s', 2] = `sdx1y'
            matrix `M'[`s', 3] = `tr'
            matrix `M'[`s', 4] = `b_y'
            matrix `M'[`s', 5] = `b_r'
            matrix `M'[`s', 6] = `b_rif'
            matrix `M'[`s', 7] = `ols'
        }
        mata: st_matrix("_p3m", mean(st_matrix("`M'")))
        mata: st_matrix("_p3s", sqrt(diagonal(variance(st_matrix("`M'")))'))
        di as text %5.2f `t' as result %6.2f _p3m[1,1] %10.2f _p3m[1,2] %8.2f _p3m[1,3] ///
           %9.2f _p3m[1,4] " [" %4.2f _p3s[1,4] "]" %9.2f _p3m[1,5] " [" %4.2f _p3s[1,5] "]" ///
           %9.2f _p3m[1,6] " [" %4.2f _p3s[1,6] "]" %8.2f _p3m[1,7]
        file write `fh' "`sig',`t',`=_p3m[1,1]',`=_p3m[1,2]',`=_p3m[1,3]',`=_p3m[1,4]',`=_p3s[1,4]',`=_p3m[1,5]',`=_p3s[1,5]',`=_p3m[1,6]',`=_p3s[1,6]',`=_p3m[1,7]'" _n
    }
}
file close `fh'
di as text _n "written: results/p3_note2016.csv, p3_mc2023.csv"
