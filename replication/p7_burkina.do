* p7_burkina.do -- Paper (version 3), Section 5, the application: Burkina Faso
* 1998 household survey (bkf98I.dta, 8,478 households, 10 strata, 425 PSUs).
* Outcome: log per capita expenditure.  Regressors: household size, male head,
* socio-economic group of the head (i.gse); urban and size as heterogeneity
* variables.
*   Table A: effect at the quantiles of lexp, tau = 0.1 ... 0.9:
*            two-step het(urban size) with Taylor s.e.; two-step het(qr) with
*            the stratified cluster bootstrap (B = 50); RIF-OLS (rifhdreg);
*            the one-step of version 1.3 (rankdep) for the record.
*   Table B: the composition of the group at each tau (e(decomp)).
*   Table C: profile of the effects along household size (rankvar(size)).
*   Table D: effect of the urban location at the initial quantile
*            (initial(urban)) against the effect at the quantile of lexp.
* Output: results/p7_tableA.csv, p7_tableB.csv, p7_tableC.csv, p7_tableD.csv
* Run: do p7_burkina   (20-30 minutes: het(qr) bootstraps)

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

use bkf98I, clear
generate lexp  = ln(exppc)
generate male  = (sex == 1)
generate urban = (zone == 2)
svyset psu [pw=weight], strata(strata)
local xvars size male i.gse

* ---- Table A ---------------------------------------------------------------
tempname fa fb
file open `fa' using "results/p7_tableA.csv", write replace
file write `fa' "tau,estimator,variable,coef,se,N_eff" _n
file open `fb' using "results/p7_tableB.csv", write replace
file write `fb' "tau,zvar,mean_at_tau,mean_pop,contrib_size,contrib_male" _n
foreach t in 0.1 0.25 0.5 0.75 0.9 {
    di as text _n "{hline 76}" _n "tau = `t'" _n "{hline 76}"
    * two-step, observables, Taylor
    gepwreg lexp `xvars', per(`t') het(urban size)
    local ne = e(N_eff)
    foreach v in size male 2.gse 3.gse 4.gse 5.gse 6.gse 7.gse {
        file write `fa' "`t',het(urban size) Taylor,`v',`=_b[`v']',`=_se[`v']',`ne'" _n
    }
    tempname D
    matrix `D' = e(decomp)
    forvalues i = 1/2 {
        local zv : word `i' of urban size
        file write `fb' "`t',`zv',`=`D'[`i',1]',`=`D'[`i',2]',`=`D'[`i',3]',`=`D'[`i',4]'" _n
    }
    * two-step, conditional quantile regression, stratified cluster bootstrap
    gepwreg lexp `xvars', per(`t') het(qr) boot(50)
    foreach v in size male 2.gse 3.gse 4.gse 5.gse 6.gse 7.gse {
        file write `fa' "`t',het(qr) boot50,`v',`=_b[`v']',`=_se[`v']',`=e(N_eff)'" _n
    }
    * RIF-OLS
    rifhdreg lexp `xvars' [pw=weight], rif(q(`=100*`t''))
    foreach v in size male 2.gse 3.gse 4.gse 5.gse 6.gse 7.gse {
        file write `fa' "`t',RIF-OLS,`v',`=_b[`v']',`=_se[`v']',." _n
    }
    * the one-step of version 1.3
    gepwreg lexp `xvars', per(`t') rankdep
    foreach v in size male 2.gse 3.gse 4.gse 5.gse 6.gse 7.gse {
        file write `fa' "`t',one-step rankdep (<=1.4),`v',`=_b[`v']',`=_se[`v']',`=e(N_eff)'" _n
    }
}
file close `fa'
file close `fb'

* ---- Table C: profile along household size --------------------------------
tempname fc
file open `fc' using "results/p7_tableC.csv", write replace
file write `fc' "tau_of_size,variable,coef,se,N_eff" _n
foreach t in 0.1 0.25 0.5 0.75 0.9 {
    gepwreg lexp male urban i.gse, per(`t') rankvar(size)
    foreach v in male urban 2.gse 3.gse 4.gse 5.gse 6.gse 7.gse {
        file write `fc' "`t',`v',`=_b[`v']',`=_se[`v']',`=e(N_eff)'" _n
    }
}
file close `fc'

* ---- Table D: urban at the initial quantile --------------------------------
* Two heterogeneity models for the urban premium: on observables (het(size):
* the premium varies with household size only) and on the conditional rank
* (het(qr): the premium varies along the conditional distribution).
tempname fd
file open `fd' using "results/p7_tableD.csv", write replace
file write `fd' "tau,het,urban_at_initial_quantile,se_boot,urban_at_quantile_of_lexp,size_at_initial_quantile,size_at_quantile_of_lexp" _n
foreach t in 0.1 0.25 0.5 0.75 0.9 {
    gepwreg lexp urban size male i.gse, per(`t') het(size) initial(urban) boot(200)
    local b0 = _b[urban]
    local se0 = _se[urban]
    local b1 = e(b_atq)[1,1]
    local s0 = e(zbar_tau)[1,1]
    quietly gepwreg lexp urban size male i.gse, per(`t') het(size)
    local s1 = e(zbar_tau)[1,1]
    file write `fd' "`t',het(size),`b0',`se0',`b1',`s0',`s1'" _n
    gepwreg lexp urban size male i.gse, per(`t') het(qr) initial(urban) boot(50)
    file write `fd' "`t',het(qr),`=_b[urban]',`=_se[urban]',`=e(b_atq)[1,1]',.,." _n
}
file close `fd'
di as text _n "written: results/p7_tableA.csv ... p7_tableD.csv"
