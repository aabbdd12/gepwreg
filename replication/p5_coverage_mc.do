* p5_coverage_mc.do -- Paper (version 3), Section 4, Table "Analytical
* standard errors: Monte Carlo": the influence-function variance of the
* two-step (het(z)) and of the profile along a covariate (rankvar(z)) against
* the Monte Carlo dispersion, with 95% coverage.
*
* DGP H1: y = 5 + (1 + z) x + z + e, z ~ U(0,1), x ~ N(0,1), e ~ N(0,1),
* n = 4000, R = 500 replications, tau = 0.1, 0.5, 0.9.
* Two coverages are reported.  (i) Against the Monte Carlo mean of the
* estimator: this isolates the variance estimator (the object of the table)
* from the kernel smoothing, whose bandwidth depends on n.  (ii) Against the
* estimator's own object.
*
* THE OBJECT IS EXACT.  For het(z) it is theta(tau) = E[1 + z | y = q_tau].
* Since y | z ~ N(5 + z, (1+z)^2 + 1), this is the one-dimensional integral
* E_z[(1+z) f(q|z)] / E_z[f(q|z)] with q solving E_z[F(q|z)] = tau; the values
* come from python/pwr_exact_truth.py.  Earlier versions used the mean effect
* inside a +-0.05 window of the rank, which is that integral smoothed over a
* box of standard deviation 0.029 and differs from it by up to 0.004 here.
* For rankvar(z) the object is E[1 + z | rank z = tau] = 1 + tau, already exact.
* The value of each estimator on a sample of 400,000 is written for the record:
* at that size the bandwidth has shrunk enough that it should sit on the exact
* value, which is the consistency check the window reference could not provide.
* Output: results/p5_coverage.csv
* Run: do p5_coverage_mc   (a few minutes)

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

local n = 4000
local R = 500

* ---- the exact objects (python/pwr_exact_truth.py) --------------------------
local true_z_0p1 = 1.464797
local true_z_0p5 = 1.462897
local true_z_0p9 = 1.588059
foreach t in 0.1 0.5 0.9 {
    local k = subinstr("`t'", ".", "p", .)
    local true_r_`k' = 1 + `t'
}

* ---- the estimators on a large sample, as a consistency check ---------------
clear
quietly set obs 400000
set seed 777
gen double z = runiform()
gen double x = rnormal()
gen double e = rnormal()
gen double y = 5 + (1+z)*x + z + e
foreach t in 0.1 0.5 0.9 {
    local k = subinstr("`t'", ".", "p", .)
    quietly gepwreg y x, per(`t') het(z)
    local plim_z_`k' = _b[x]
    local hbig = e(h)
    quietly gepwreg y x, per(`t') rankvar(z)
    local plim_r_`k' = _b[x]
    di as text "tau = `t'   theta(tau) exact = " as result %7.4f `true_z_`k'' ///
       as text "   het(z) on n = 400000: " as result %7.4f `plim_z_`k'' ///
       as text " (h = " as result %6.4f `hbig' as text ")" ///
       as text "   | 1 + tau = " as result %5.3f `true_r_`k'' ///
       as text "   rankvar(z) on n = 400000: " as result %7.4f `plim_r_`k''
}

* ---- Monte Carlo -------------------------------------------------------------
tempname fh
file open `fh' using "results/p5_coverage.csv", write replace
file write `fh' "estimator,tau,true_object,value_n400000,mean_est,mc_sd,mean_seIF,ratio_seIF_mcsd,coverage95_mcmean,coverage95_true" _n
foreach t in 0.1 0.5 0.9 {
    local k = subinstr("`t'", ".", "p", .)
    tempname M
    matrix `M' = J(`R', 4, .)
    forvalues s = 1/`R' {
        clear
        quietly set obs `n'
        set seed `=2000 + `s''
        gen double z = runiform()
        gen double x = rnormal()
        gen double e = rnormal()
        gen double y = 5 + (1+z)*x + z + e
        quietly gepwreg y x, per(`t') het(z)
        matrix `M'[`s', 1] = _b[x]
        matrix `M'[`s', 2] = _se[x]
        quietly gepwreg y x, per(`t') rankvar(z)
        matrix `M'[`s', 3] = _b[x]
        matrix `M'[`s', 4] = _se[x]
    }
    foreach est in z r {
        if "`est'" == "z" {
            local c1 = 1
            local name "het(z)"
        }
        else {
            local c1 = 3
            local name "rankvar(z)"
        }
        local c2 = `c1' + 1
        mata: M = st_matrix("`M'")
        mata: b = M[., `c1']
        mata: se = M[., `c2']
        mata: st_numscalar("_mb",  mean(b))
        mata: st_numscalar("_msd", sqrt(variance(b)))
        mata: st_numscalar("_mse", mean(se))
        mata: st_numscalar("_cov_m", mean(abs(b :- mean(b)) :<= 1.96 :* se))
        mata: st_numscalar("_cov_t", mean(abs(b :- `true_`est'_`k'') :<= 1.96 :* se))
        di as text "`name'  tau = `t':  own object = " as result %7.4f `true_`est'_`k'' ///
           as text "  mean = " as result %7.4f _mb as text "  MC sd = " as result %6.4f _msd ///
           as text "  mean s.e. IF = " as result %6.4f _mse as text "  ratio = " as result %5.3f _mse/_msd ///
           as text "  coverage (MC mean) = " as result %5.3f _cov_m as text "  (own object) = " as result %5.3f _cov_t
        file write `fh' "`name',`t',`true_`est'_`k'',`plim_`est'_`k'',`=_mb',`=_msd',`=_mse',`=_mse/_msd',`=_cov_m',`=_cov_t'" _n
    }
}
file close `fh'
di as text _n "written: results/p5_coverage.csv"
