* p6_svy_mc.do -- Paper (version 3), Section 4, Table "Inference under a
* complex survey design": the Taylor variance of the two-step (het(z)) under
* strata and PSUs against the Monte Carlo dispersion, with the IF variance
* that ignores the design for comparison, and one stratified cluster bootstrap.
*
* Superpopulation with clustering: 10 strata x 40 PSUs x 20 households
* (n = 8000).  Two independent PSU effects: a_c ~ N(0, 0.6) in the level of y
* and c_c ~ N(0, 0.6) in x (so that the intra-cluster correlation matters for
* the effect without making x endogenous), stratum effect s_h in the level,
* sampling weights w = 1 + 0.1 h (constant within a stratum).
*   y = 5 + s_h + (1 + z) x + z + a_c + e,   z ~ U(0,1), e ~ N(0,1)
*   x = c_c + N(0,1)
* R = 300 replications, tau = 0.25 and 0.75.  Output: results/p6_svy.csv
* Run: do p6_svy_mc   (a few minutes)
*
* THE TRUTH IS EXACT.  The cluster effects integrate out in closed form: given
* z and the stratum, Var[(1+z) c_c + a_c + e | z] = 0.36 (1+z)^2 + 0.36 + 1, so
*     y | z, h ~ N(5 + s_h + z, 1.36 [(1+z)^2 + 1])
* and theta(tau) = E[1 + z | y = q_tau] is a one-dimensional integral over z,
* mixed over the ten strata with the SAMPLING WEIGHTS w_h -- the estimator
* targets the weighted population, not the sample.  Values from
* python/pwr_exact_truth.py.  Earlier versions averaged the +-0.05 window
* truth over the replications, which is the estimand smoothed over a box.

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

local H = 10
local C = 40
local m = 20
local R = 300

capture program drop _p6draw
program define _p6draw
    args H C m
    clear
    quietly set obs `=`H'*`C''
    gen int strata = ceil(_n/`C')
    gen int psu    = _n
    gen double a   = rnormal(0, 0.6)
    gen double c   = rnormal(0, 0.6)
    gen double s   = 0.3*(strata - 5.5)/3
    expand `m'
    gen double z = runiform()
    gen double e = rnormal()
    gen double x = c + rnormal()
    gen double y = 5 + s + (1+z)*x + z + a + e
    gen double w = 1 + 0.1*strata
    gen double d = 1 + z
end

tempname fh
file open `fh' using "results/p6_svy.csv", write replace
file write `fh' "tau,true_exact,mean_est,mc_sd,mean_se_taylor,ratio_taylor,cov95_taylor,mean_se_IF_nodesign,ratio_IF,cov95_IF,boot_cluster_se_one_sample" _n

* exact theta(tau) for the weighted population (python/pwr_exact_truth.py)
local TRUE_0p25 = 1.452892
local TRUE_0p75 = 1.508806

foreach t in 0.25 0.75 {
    local k = subinstr("`t'", ".", "p", .)
    scalar _tr = `TRUE_`k''
    tempname M
    matrix `M' = J(`R', 3, .)
    forvalues r = 1/`R' {
        set seed `=3000 + `r''
        _p6draw `H' `C' `m'
        svyset psu [pw=w], strata(strata)
        quietly gepwreg y x, per(`t') het(z)
        matrix `M'[`r', 1] = _b[x]
        matrix `M'[`r', 2] = _se[x]
        matrix `M'[`r', 3] = sqrt(e(V_IF)[1,1])
    }
    mata: M = st_matrix("`M'")
    mata: st_numscalar("_mb",  mean(M[., 1]))
    mata: st_numscalar("_msd", sqrt(variance(M[., 1])))
    mata: st_numscalar("_mst", mean(M[., 2]))
    mata: st_numscalar("_msi", mean(M[., 3]))
    mata: st_numscalar("_cvt", mean(abs(M[., 1] :- mean(M[., 1])) :<= 1.96 :* M[., 2]))
    mata: st_numscalar("_cvi", mean(abs(M[., 1] :- mean(M[., 1])) :<= 1.96 :* M[., 3]))
    * one sample with the stratified cluster bootstrap for the record
    set seed 3000
    _p6draw `H' `C' `m'
    svyset psu [pw=w], strata(strata)
    quietly gepwreg y x, per(`t') het(z) boot(200)
    local sb = sqrt(e(V_boot)[1,1])
    di as text _n "tau = `t':  theta(tau) exact = " as result %7.4f _tr as text "   mean estimate = " as result %7.4f _mb ///
       as text "   MC sd = " as result %6.4f _msd
    di as text "   Taylor s.e. (mean) = " as result %6.4f _mst as text "  ratio = " as result %5.3f _mst/_msd ///
       as text "  coverage = " as result %5.3f _cvt
    di as text "   IF s.e. ignoring the design (mean) = " as result %6.4f _msi as text "  ratio = " as result %5.3f _msi/_msd ///
       as text "  coverage = " as result %5.3f _cvi
    di as text "   stratified cluster bootstrap, one sample (B = 200) = " as result %6.4f `sb'
    file write `fh' "`t',`=_tr',`=_mb',`=_msd',`=_mst',`=_mst/_msd',`=_cvt',`=_msi',`=_msi/_msd',`=_cvi',`sb'" _n
}
file close `fh'
di as text _n "written: results/p6_svy.csv"
