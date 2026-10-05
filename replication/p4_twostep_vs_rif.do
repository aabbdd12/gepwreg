* p4_twostep_vs_rif.do -- Paper (version 3), Section 6, Table "The two-step
* estimator against RIF and against the truth", and the price of the
* heterogeneity model.
*
* DGPs (n = 6000, R = 25 seeds):
*   D3: y = 5 z + (1 + z) x,            z ~ U(0,1) observed, x ~ N(0,1)   effect 1 + z
*   H1: y = 5 + (1 + z) x + z + e,      z ~ U(0,1) observed, e ~ N(0,1)   effect 1 + z
*   D1: y = 100 z + (1 + z) x                                            effect 1 + z
*   U1: y = 5 + (1 + U) x + 3 U,        U unobserved, x ~ U(0,2),
*       z ~ N(0,1) observed and irrelevant                               effect 1 + U
*
* THE TRUTH IS EXACT.  Earlier versions used the mean effect of the units whose
* rank in y falls within +-0.05 of tau.  That is not the estimand: it is the
* estimand smoothed over a box of standard deviation 0.05/sqrt(3) = 0.029,
* against a kernel of standard deviation h*sqrt(2) = 0.065 at the Silverman
* rule with n = 6000.  The window therefore carries about a fifth of the
* estimator's own smoothing, with the same sign, and subtracting it removes a
* fifth of the bias at the Silverman bandwidth but only a twentieth of the bias
* at twice that bandwidth -- it tilts any bandwidth comparison towards the
* smaller one.  The tilt is at most 0.004 here, negligible against biases of
* order 0.1 but not against the bandwidth question.  Each DGP has a latent
* variable z (or U) carrying the effect and a closed-form law of y given it, so
*     theta(tau) = E_z[(1+z) f(q_tau|z)] / E_z[f(q_tau|z)]
* is a one-dimensional integral.  The values below come from
* python/pwr_exact_truth.py, which also verifies them against a 4,000,000-draw
* simulation with the window width extrapolated to zero, and against the closed
* form available for U1.
*
* Estimators: one-step on y (rankdep), two-step het(z) [het(r)] with the
* MSE-optimal bandwidth of version 1.4 and with the Silverman rule, two-step
* het(qr), RIF-OLS (rifhdreg).  Output: results/p4_twostep.csv
* Run: do p4_twostep_vs_rif    (about 45 minutes: qreg)

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

* This file compares against RIF regression, which needs rifhdreg.  It is not
* part of the gepwreg package: it belongs to the SSC module rif, and the line
* below installs it if it is missing.  The version these results were produced
* with is 2.55 of August 2021; the copy archived with the Stata Journal
* article is an earlier one, 2.5, whose header records a bug in the definition
* of the estimation sample.
capture which rifhdreg
if _rc {
    display as text "rifhdreg not found -- installing the SSC module rif"
    ssc install rif
}
which rifhdreg

local n = 6000
local R = 25
tempname fh
file open `fh' using "results/p4_twostep.csv", write replace
file write `fh' "dgp,tau,true,onestep_m,onestep_sd,hetz_m,hetz_sd,hetz_seIF,hetzsil_m,hetzsil_sd,hetqr_m,hetqr_sd,rif_m,rif_sd,h_ratio" _n

foreach dgp in D3 H1 D1 U1 {

    * exact theta(tau), python/pwr_exact_truth.py
    if "`dgp'" == "D3" {
        local TRUE_0p1  = 1.231125
        local TRUE_0p25 = 1.319275
        local TRUE_0p5  = 1.492782
        local TRUE_0p75 = 1.672133
        local TRUE_0p9  = 1.780590
    }
    if "`dgp'" == "H1" {
        local TRUE_0p1  = 1.464797
        local TRUE_0p25 = 1.447873
        local TRUE_0p5  = 1.462897
        local TRUE_0p75 = 1.513874
        local TRUE_0p9  = 1.588059
    }
    if "`dgp'" == "D1" {
        local TRUE_0p1  = 1.100110
        local TRUE_0p25 = 1.250125
        local TRUE_0p5  = 1.500150
        local TRUE_0p75 = 1.750175
        local TRUE_0p9  = 1.900190
    }
    if "`dgp'" == "U1" {
        local TRUE_0p1  = 1.182857
        local TRUE_0p25 = 1.289697
        local TRUE_0p5  = 1.532903
        local TRUE_0p75 = 1.683244
        local TRUE_0p9  = 1.799876
    }

    di as text _n "{hline 100}" _n "DGP `dgp'   (n = `n', `R' seeds; mean [sd], exact truth)" _n "{hline 100}"
    di as text %5s "tau" %8s "true" %15s "one-step y" %15s "het(z)" %8s "se_IF" ///
       %15s "het(z) Silv" %15s "het(qr)" %15s "RIF-OLS" %7s "h/hS"
    foreach t in 0.1 0.25 0.5 0.75 0.9 {
        local k = subinstr("`t'", ".", "p", .)
        local tr = `TRUE_`k''
        tempname M
        matrix `M' = J(`R', 7, .)
        forvalues s = 1/`R' {
            clear
            quietly set obs `n'
            set seed `=1000 + `s''
            if "`dgp'" == "D3" {
                gen double z = runiform()
                gen double x = rnormal()
                gen double y = 5*z + (1+z)*x
            }
            if "`dgp'" == "H1" {
                gen double z = runiform()
                gen double x = rnormal()
                gen double e = rnormal()
                gen double y = 5 + (1+z)*x + z + e
            }
            if "`dgp'" == "D1" {
                gen double z = runiform()
                gen double x = rnormal()
                gen double y = 100*z + (1+z)*x
            }
            if "`dgp'" == "U1" {
                gen double U = runiform()
                gen double x = 2*runiform()
                gen double z = rnormal()
                gen double y = 5 + (1+U)*x + 3*U
            }
            quietly gepwreg y x, per(`t') rankdep
            matrix `M'[`s', 1] = _b[x]
            quietly gepwreg y x, per(`t') het(z)
            matrix `M'[`s', 2] = _b[x]
            matrix `M'[`s', 3] = _se[x]
            local hopt = e(h)
            quietly gepwreg y x, per(`t') het(z) silverman
            matrix `M'[`s', 4] = _b[x]
            matrix `M'[`s', 7] = `hopt'/e(h)
            quietly gepwreg y x, per(`t') het(qr)
            matrix `M'[`s', 5] = _b[x]
            quietly rifhdreg y x, rif(q(`=100*`t''))
            matrix `M'[`s', 6] = _b[x]
        }
        mata: st_matrix("_p4m", mean(st_matrix("`M'")))
        mata: st_matrix("_p4s", sqrt(diagonal(variance(st_matrix("`M'")))'))
        di as text %5.2f `t' as result %8.4f `tr' ///
           %8.3f _p4m[1,1] " [" %5.3f _p4s[1,1] "]" ///
           %8.3f _p4m[1,2] " [" %5.3f _p4s[1,2] "]" %8.4f _p4m[1,3] ///
           %8.3f _p4m[1,4] " [" %5.3f _p4s[1,4] "]" ///
           %8.3f _p4m[1,5] " [" %5.3f _p4s[1,5] "]" ///
           %8.3f _p4m[1,6] " [" %5.3f _p4s[1,6] "]" %7.2f _p4m[1,7]
        file write `fh' "`dgp',`t',`tr',`=_p4m[1,1]',`=_p4s[1,1]',`=_p4m[1,2]',`=_p4s[1,2]',`=_p4m[1,3]',`=_p4m[1,4]',`=_p4s[1,4]',`=_p4m[1,5]',`=_p4s[1,5]',`=_p4m[1,6]',`=_p4s[1,6]',`=_p4m[1,7]'" _n
    }
}
file close `fh'
di as text _n "written: results/p4_twostep.csv"
