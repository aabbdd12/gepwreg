* p7d_tableD.do -- Paper (version 3), Section 5, Table D only (the urban premium
* at the initial quantile), rerun after the fix of the decomposition display in
* initial mode and with the het(qr) column added.  Same code as the Table D
* section of p7_burkina.do.  Output: results/p7_tableD.csv
* Run: do p7d_tableD   (about 10 minutes)

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

use bkf98I, clear
generate lexp  = ln(exppc)
generate male  = (sex == 1)
generate urban = (zone == 2)
svyset psu [pw=weight], strata(strata)

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
di as text _n "written: results/p7_tableD.csv"
