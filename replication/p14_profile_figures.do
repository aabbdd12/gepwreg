* p14_profile_figures.do -- Paper (version 5), Section 9.6, Figures 1 and 2:
* the profile of four effects at 37 quantiles from 0.05 to 0.95 on the Burkina
* Faso 1998 survey (bkf98I.dta, svyset: the design is taken from the file),
* regressors size, male head, urban, socio-economic group, with the
* reference household, the simultaneous band and the test of constancy.
*   Figure 1: het(urban size), the mean household (ref(mean)), Taylor s.e.;
*   Figure 2: het(qr), the default real type (ref(default)), the stratified
*             cluster bootstrap, B = 50, one loop of draws for every quantile.
* These are examples 8 and 9 of the help of gepwreg.
* Output: results/fig_profile_hetz.png, results/fig_profile_hetqr.png (the
* figures); results/p14_profile_hetz.csv, results/p14_profile_hetqr.csv
* (e(profile): tau, then for each coefficient b se lb ub ulb uub rb rse rlb
* rub rulb ruub); results/p14_profile_summary.csv (the numbers quoted in
* Section 9.6: the estimates at per() and for the reference household, the
* p-values of the test of constancy, the composition of the group, the
* conditional rank of the reference household and where it is held).
* Run: do p14_profile_figures   (about two minutes)

* ---------------------------------------------------------------------
* Replication, Araar (2026), version 5.  Run this file from the folder
* that contains it; the results are written to the results/ folder
* beside it, which is created if it does not exist.
*
* It needs gepwreg 1.5.0:
*
*   net install gepwreg, from("https://raw.githubusercontent.com/aabbdd12/gepwreg/v15") replace
*   net get     gepwreg, from("https://raw.githubusercontent.com/aabbdd12/gepwreg/v15") replace
*
* The second line also brings bkf98I.dta.  The which gepwreg below must
* report v1.5.0.
* ---------------------------------------------------------------------
version 16
clear all
set more off
set varabbrev off
which gepwreg
capture mkdir "results"

use bkf98I, clear
generate lexp  = ln(exppc)
generate male  = (sex == 1)
generate urban = (zone == 2)
local coefs "size male 5.gse urban"

tempname fh B V BR VR F RF P Z ZP RU
file open `fh' using "results/p14_profile_summary.csv", write replace
file write `fh' "model,quantity,coefficient,value" _n

* ---- Figure 1: het(urban size), the mean household ----------------------------
gepwreg lexp size male urban i.gse, per(0.25) het(urban size) ref(mean) graph(`coefs', ref uniform)
graph export "results/fig_profile_hetz.png", name(gepwreg_profile) replace width(1400)
local m "het(urban size)"
matrix `B'  = e(b)
matrix `V'  = e(V)
matrix `BR' = e(b_ref)
matrix `VR' = e(V_ref)
matrix `F'  = e(prof_flat)
matrix `P'  = e(profile)
foreach c of local coefs {
    local j = colnumb(`B', "`c'")
    file write `fh' "`m',b at per(0.25),`c',`=`B'[1, `j']'" _n
    file write `fh' "`m',se at per(0.25),`c',`=sqrt(`V'[`j', `j'])'" _n
    local j = colnumb(`BR', "`c'")
    file write `fh' "`m',b of the reference,`c',`=`BR'[1, `j']'" _n
    file write `fh' "`m',se of the reference,`c',`=sqrt(`VR'[`j', `j'])'" _n
    local j = colnumb(`F', "`c'")
    file write `fh' "`m',p constant across tau,`c',`=`F'[2, `j']'" _n
}
* the composition of the group at per() against the population (e(decomp))
matrix `Z'  = e(zbar_tau)
matrix `ZP' = e(zbar_pop)
foreach z in urban size {
    local j = colnumb(`Z', "`z'")
    file write `fh' "`m',mean at per(0.25),`z',`=`Z'[1, `j']'" _n
    local j = colnumb(`ZP', "`z'")
    file write `fh' "`m',mean in the population,`z',`=`ZP'[1, `j']'" _n
}
preserve
clear
quietly svmat double `P', names(col)
export delimited using "results/p14_profile_hetz.csv", replace
restore

* ---- Figure 2: het(qr), a given household --------------------------------------
gepwreg lexp size male urban i.gse, per(0.25) ref(default) boot(50) graph(`coefs', ref uniform)
graph export "results/fig_profile_hetqr.png", name(gepwreg_profile) replace width(1400)
local m "het(qr)"
matrix `B'  = e(b)
matrix `V'  = e(V)
matrix `BR' = e(b_ref)
matrix `VR' = e(V_ref)
matrix `F'  = e(prof_flat)
matrix `RF' = e(prof_rflat)
matrix `P'  = e(profile)
matrix `RU' = e(prof_ru)
foreach c of local coefs {
    local j = colnumb(`B', "`c'")
    file write `fh' "`m',b at per(0.25),`c',`=`B'[1, `j']'" _n
    file write `fh' "`m',se at per(0.25),`c',`=sqrt(`V'[`j', `j'])'" _n
    local j = colnumb(`BR', "`c'")
    file write `fh' "`m',b of the reference,`c',`=`BR'[1, `j']'" _n
    file write `fh' "`m',se of the reference,`c',`=sqrt(`VR'[`j', `j'])'" _n
    local j = colnumb(`F', "`c'")
    file write `fh' "`m',p constant across tau,`c',`=`F'[2, `j']'" _n
    file write `fh' "`m',p constant across tau (reference),`c',`=`RF'[2, `j']'" _n
}
file write `fh' "`m',conditional rank of the reference at per(0.25),,`=e(u_ref)'" _n
* the first quantile of the figure from which the reference is held at the
* last quantile of the grid
local held = .
forvalues i = 1/`=rowsof(`P')' {
    if `RU'[`i', 2] == 1 & `held' == . local held = `P'[`i', 1]
}
file write `fh' "`m',reference held at the end of the grid from tau,,`held'" _n
preserve
clear
quietly svmat double `P', names(col)
export delimited using "results/p14_profile_hetqr.csv", replace
restore
file close `fh'
type "results/p14_profile_summary.csv"
