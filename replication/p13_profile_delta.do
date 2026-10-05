* p13_profile_delta.do -- Paper (version 5), Section 7.5: the delta-method
* standard error of the effect for a given household under het(qr) against
* the bootstrap of the whole procedure, with two densities for the
* quantile-regression influence function (Powell kernel, Hendricks-Koenker).
*   The perfect case (every coefficient varies strongly with u, n = 2000,
*   tau = 0.25): "within 4 to 13 percent ... and within 1 to 11 percent";
*   Burkina (bkf98I.dta, lexp on size male urban i.gse, unweighted,
*   tau = 0.25, the default reference): "by up to 40 percent (household
*   size: 0.90 ... with the Powell kernel, 0.59 with Hendricks-Koenker)".
* Runs p13a_E_draws.do (the point, its delta-method standard error and the
* B = 200 draws, saved in results/) and p13b_E_densities.do (the two
* standard errors against the draws).  The numbers are the columns kern/SD
* and HK/SD of p13b's two reports.
* Output: results/p13_profile_delta.log
* Run: do p13_profile_delta   (about five minutes)

* ---------------------------------------------------------------------
* Replication, Araar (2026), version 5.  Run this file from the folder
* that contains it, with bkf98I.dta there; the results are written to the
* results/ folder beside it, which is created if it does not exist.
* It does not call gepwreg: the estimator is programmed in p13a and p13b
* with official quantile regressions.
* ---------------------------------------------------------------------
version 16
clear all
set more off
capture mkdir "results"
capture log close _p13
log using "results/p13_profile_delta.log", text replace name(_p13)
do p13a_E_draws dgp
forvalues c = 1/3 {
    do p13a_E_draws bkf `c'
}
do p13a_E_draws bkfreport
do p13b_E_densities dgp
do p13b_E_densities bkf
log close _p13
