* ---------------------------------------------------------------------
* Replication, Araar (2026), version 3.  Run this file from the folder
* that contains it.
*
* It needs gepwreg 1.4, the version the paper was produced with:
*
*   net install gepwreg, from("https://raw.githubusercontent.com/aabbdd12/gepwreg/v14") replace
*   net get     gepwreg, from("https://raw.githubusercontent.com/aabbdd12/gepwreg/v14") replace
*
* The second line also brings bkf98I.dta, which Section 3 below uses.
* The which gepwreg below must report v1.4.
* ---------------------------------------------------------------------
version 16
clear all
set more off
set trace off
which gepwreg
which gepwreg_setable

* ---------------------------------------------------------------------------
* 1. H1: y = 5 + (1+z) x + z + e, z observed.  het(z) is exact for this DGP.
*    Checks: truth vs het(z); closed form by hand (regress + nlcom);
*            IF s.e. vs pairs bootstrap; het(qr) at the median.
* ---------------------------------------------------------------------------
set seed 20260916
set obs 20000
gen double x = rnormal()
gen double z = runiform()
gen double e = rnormal()
gen double y = 5 + (1+z)*x + z + e
gen double imp = 1 + z
gen double y0 = 5 + z + e
sort y
gen double pcy = _n/_N
sort y0
gen double pcy0 = _n/_N

di as text _n "{hline 76}" _n "TEST 1  H1: y = 5 + (1+z)x + z + e   het(z)" _n "{hline 76}"
quietly summarize pcy, detail
local hs = 0.9*min(r(sd), (r(p75)-r(p25))/1.34)*_N^(-0.2)
foreach t in 0.1 0.5 0.9 {
    quietly summarize imp if abs(pcy - `t') <= 0.05
    local tr = r(mean)
    * closed form by hand, before the command (regress would overwrite e())
    quietly gen double w = exp(-0.25*((pcy - `t')/`hs')^2)
    quietly summarize z [aw=w]
    local zb = r(mean)
    quietly regress y c.x##c.z
    local cf = _b[x] + `zb'*_b[c.x#c.z]
    drop w
    gepwreg y x, per(`t') het(z) boot(200)
    local th = _b[x]
    local se_if = sqrt(e(V_IF)[1,1])
    local se_bt = sqrt(e(V_boot)[1,1])
    di as text "tau = `t'   true = " as result %7.4f `tr' as text "   het(z) = " as result %7.4f `th' ///
       as text "   closed form by hand = " as result %7.4f `cf' ///
       as text "   s.e. IF = " as result %6.4f `se_if' as text "  boot = " as result %6.4f `se_bt'
}
gepwreg_setable

di as text _n "{hline 76}" _n "TEST 2  H1: het(qr) (default) at tau = 0.5, boot(20)" _n "{hline 76}"
quietly summarize imp if abs(pcy - 0.5) <= 0.05
local tr = r(mean)
gepwreg y x, per(0.5) boot(20)
di as text "true = " as result %7.4f `tr' as text "   het(qr) = " as result %7.4f _b[x] ///
   as text "   s.e. boot = " as result %6.4f sqrt(e(V)[1,1])
gepwreg y x, per(0.5)
di as text "without boot(): coefficient only, " as result %7.4f _b[x]

di as text _n "{hline 76}" _n "TEST 3  H1: initial -- effect at the tau-quantile of y0 = 5 + z + e" _n "{hline 76}"
foreach t in 0.1 0.5 0.9 {
    quietly summarize imp if abs(pcy0 - `t') <= 0.05
    local tr0 = r(mean)
    quietly summarize imp if abs(pcy - `t') <= 0.05
    local tr = r(mean)
    gepwreg y x, per(`t') het(z) initial boot(100)
    di as text "tau = `t'   true at y0 = " as result %7.4f `tr0' as text "   initial = " as result %7.4f _b[x] ///
       as text "   s.e. = " as result %6.4f sqrt(e(V)[1,1]) ///
       as text "   | true at y = " as result %7.4f `tr' as text "   e(b_atq) = " as result %7.4f e(b_atq)[1,1]
}

di as text _n "{hline 76}" _n "TEST 4  rankvar(z) and rankdep unchanged" _n "{hline 76}"
gepwreg y x, per(0.5) rankvar(z)
di as text "rankvar(z) at tau = 0.5: " as result %7.4f _b[x] as text "  (true 1 + 0.5 = 1.5)"
gepwreg y x, per(0.5) rankdep
di as text "rankdep at tau = 0.5: " as result %7.4f _b[x] as text "  (descriptive; far below 1.46)"
gepwreg_setable

* ---------------------------------------------------------------------------
* 4b. A regressor that is also a het() variable: y = 5 + (1+z)x - 0.3x^2 + z + e,
*     derivative 1 + z - 0.6x; het(z x) must return E[1 + z - 0.6x | rank y]
* ---------------------------------------------------------------------------
clear
set seed 20260918
set obs 20000
gen double x = rnormal()
gen double z = runiform()
gen double e = rnormal()
gen double y = 5 + (1+z)*x - 0.3*x^2 + z + e
gen double der = 1 + z - 0.6*x
sort y
gen double pcy = _n/_N
di as text _n "{hline 76}" _n "TEST 4b  y = 5 + (1+z)x - 0.3x^2 + z + e,  het(z x)" _n "{hline 76}"
foreach t in 0.1 0.5 0.9 {
    quietly summarize der if abs(pcy - `t') <= 0.05
    local tr = r(mean)
    quietly gepwreg y x, per(`t') het(z x) boot(100)
    di as text "tau = `t'   true = " as result %7.4f `tr' as text "   het(z x) = " as result %7.4f _b[x] ///
       as text "   s.e. IF = " as result %6.4f sqrt(e(V_IF)[1,1]) as text "  boot = " as result %6.4f sqrt(e(V_boot)[1,1])
}
gepwreg y x, per(0.5) het(z x)

* ---------------------------------------------------------------------------
* 5. Unobserved heterogeneity ordered with y: y = 5 + (1+U)x + 3U, x ~ U(0,2)
*    het(qr) recovers E[1+U | rank y]; het(z) with an irrelevant z is flat.
* ---------------------------------------------------------------------------
clear
set seed 20260917
set obs 20000
gen double U = runiform()
gen double x = 2*runiform()
gen double z = rnormal()
gen double y = 5 + (1+U)*x + 3*U
gen double imp = 1 + U
sort y
gen double pcy = _n/_N
di as text _n "{hline 76}" _n "TEST 5  y = 5 + (1+U)x + 3U, U unobserved" _n "{hline 76}"
foreach t in 0.1 0.5 0.9 {
    quietly summarize imp if abs(pcy - `t') <= 0.05
    local tr = r(mean)
    quietly gepwreg y x, per(`t')
    local bq = _b[x]
    quietly gepwreg y x, per(`t') het(z)
    di as text "tau = `t'   true = " as result %7.4f `tr' as text "   het(qr) = " as result %7.4f `bq' ///
       as text "   het(z irrelevant) = " as result %7.4f _b[x]
}

* ---------------------------------------------------------------------------
* 6. Burkina: factor variables, survey design, Taylor s.e., decomposition
* ---------------------------------------------------------------------------
di as text _n "{hline 76}" _n "TEST 6  bkf98I: het(urban size), svyset, factor variables" _n "{hline 76}"
use bkf98I, clear
generate lexp  = ln(exppc)
generate male  = (sex == 1)
generate urban = (zone == 2)
gepwreg lexp size male i.gse [pw=weight], per(0.25) het(urban size)
gepwreg_setable
svyset psu [pw=weight], strata(strata)
gepwreg lexp size male i.gse, per(0.25) het(urban size)
gepwreg lexp size male i.gse, per(0.25) het(urban size) vce(if) boot(100)
gepwreg_setable
gepwreg lexp size male i.gse, per(0.25) het(qr) boot(20)
gepwreg lexp size male i.gse, per(0.25) rankvar(size)
gepwreg lexp size male i.gse, per(0.25) rankdep
rifhdreg lexp size male i.gse [pw=weight], rif(q(25))
