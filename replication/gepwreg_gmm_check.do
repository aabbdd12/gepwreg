* ---------------------------------------------------------------------
* Replication, Araar (2026), version 3.  Run this file from the folder
* that contains it.  It generates its own data and needs no survey file.
*
* It needs gepwreg 1.4, the version the paper was produced with:
*
*   net install gepwreg, from("https://raw.githubusercontent.com/aabbdd12/gepwreg/v14") replace
*
* The which gepwreg below must report v1.4.
* ---------------------------------------------------------------------
*
* gepwreg_gmm_check.do -- the two-step (observables variant) as a stacked GMM,
* and its closed form, with four standard errors.
*
* DGP: y = 5 + (1 + z) x + z + e,  z ~ U(0,1) observed, e ~ N(0,1); the
* effect of x for household i is 1 + z_i.  Object: E[1 + z | rank y = tau].
*
* Step 1 (interaction model): y = a + bx x + bz z + bxz x z + e  ->  b_i = bx + bxz z_i
* Step 2: theta(tau) = sum w_i b_i / sum w_i, w_i = kernel on the rank of y
*       = bx + bxz * zbar_w(tau)          (closed form: effect x composition)
*
* Stacked GMM (Newey-McFadden), exactly identified:
*   E[(1, x, z, xz)' (y - a - bx x - bz z - bxz xz)] = 0        (4 moments)
*   E[ w (bx + bxz z - theta) ] = 0                            (1 moment)
* gmm's robust variance carries the step-1 estimation error into theta.  What
* it does not carry is the estimation of the ranks inside w; the influence
* function of gepwreg carries that term as well, and the pairs bootstrap
* carries everything.
*
* The four standard errors are therefore ordered a priori:
*   delta (composition and ranks fixed)  <  gmm (step 1 only)
*                                        <  influence function (ranks)
*                                       ~=  pairs bootstrap.
* Every one of them is computed at the SAME bandwidth hs, which is why
* band(`hs') is passed to gepwreg: at two different bandwidths the four
* numbers would not be comparable and their ordering would mean nothing.
*
* Expected (n = 20000, at band(hs) = 0.0358).  theta from nlcom, gmm, gepwreg
* and the bootstrap agree to four decimals: 1.4753 / 1.4734 / 1.6136 at
* tau = 0.1 / 0.5 / 0.9, against a truth of 1.461 / 1.461 / 1.597.  The four
* standard errors come out in the order the theory predicts:
*
*   delta               0.0072  0.0072  0.0075
*   gmm robust          0.0086  0.0087  0.0086
*   influence function  0.0088  0.0088  0.0091
*   pairs bootstrap     0.0091  0.0090  0.0094     (+3.0  +2.8  +3.1 percent)
*
* WHY 2000 REPLICATIONS AND NOT 200.  The standard error of a bootstrap
* standard error is about 1/sqrt(2B), which is five percent at B = 200 -- the
* size of the gap being measured.  At B = 200 this file returned 0.0094,
* 0.0097 and 0.0101, and the excess over the influence function looked like
* six, ten and eleven percent; at B = 2000 it is three percent at all three
* quantiles.  Two thirds of the apparent disagreement was the instrument.
* Comparing a bootstrap standard error with an analytical one at the percent
* level needs B in the thousands, whatever the point estimate costs.

version 16
clear all
set more off
which gepwreg

set seed 20260916
set obs 20000
gen double x  = rnormal()
gen double z  = runiform()
gen double e  = rnormal()
gen double y  = 5 + (1+z)*x + z + e
gen double xz = x*z
gen double imp = 1 + z
gen byte one = 1

sort y
gen double pc = _n/_N
quietly summarize pc, detail
local hs = 0.9*min(r(sd), (r(p75)-r(p25))/1.34)*_N^(-0.2)
di as text _n "bandwidth used everywhere below: " as result %6.4f `hs'

capture program drop twostep_z
program define twostep_z, rclass
    syntax , tau(real) h(real)
    tempvar pc w bi
    sort y
    gen double `pc' = _n/_N
    gen double `w'  = exp(-0.25*((`pc' - `tau')/`h')^2)
    quietly regress y x z xz
    gen double `bi' = _b[x] + _b[xz]*z
    quietly summarize `bi' [aw=`w']
    return scalar theta = r(mean)
end

foreach t in 0.1 0.5 0.9 {
    di as text _n "{hline 72}" _n "tau = `t'" _n "{hline 72}"
    quietly summarize imp if abs(pc - `t') <= 0.05
    di as text "true E[1+z | rank y ~ tau]      = " as result %8.4f r(mean)

    * closed form + delta method (ranks and composition treated as fixed)
    quietly gen double w = exp(-0.25*((pc - `t')/`hs')^2)
    quietly summarize z [aw=w]
    local zbar = r(mean)
    quietly regress y x z xz
    quietly nlcom _b[x] + `zbar'*_b[xz]
    matrix b = r(b)
    matrix V = r(V)
    di as text "closed form  bx + bxz*zbar_w    = " as result %8.4f b[1,1] as text "   s.e. (delta, ranks fixed)   = " as result %7.4f sqrt(V[1,1])

    * stacked GMM, ranks fixed
    capture noisily gmm (y - {a} - {bx}*x - {bz}*z - {bxz}*xz) (w*({bx} + {bxz}*z - {theta})), ///
        instruments(1: x z xz) instruments(2: one, noconstant) winitial(identity) onestep
    if _rc == 0 {
        di as text "stacked GMM  theta              = " as result %8.4f _b[/theta] as text "   s.e. (gmm robust)           = " as result %7.4f _se[/theta]
    }
    else di as error "gmm failed with r(`=_rc'); continuing"

    * the two-step of gepwreg, het(z), with the influence-function standard
    * error that carries the rank term
    quietly gepwreg y x, per(`t') het(z) band(`hs')
    di as text "gepwreg het(z)  theta           = " as result %8.4f _b[x] as text "   s.e. (influence function)   = " as result %7.4f _se[x]

    * pairs bootstrap of the whole two-step (ranks, step 1, step 2);
    * 2000 replications, see the header
    quietly bootstrap theta = r(theta), reps(2000) seed(7) nowarn: twostep_z, tau(`t') h(`hs')
    di as text "two-step     theta              = " as result %8.4f _b[theta] as text "   s.e. (bootstrap, 2000 reps) = " as result %7.4f _se[theta]

    * the outcome-ranked slope of version 1.3, kept as rankdep and labelled
    * descriptive: it is NOT an estimate of the same object
    quietly gepwreg y x, per(`t') rankdep band(`hs')
    di as text "rankdep (descriptive slope)     = " as result %8.4f _b[x]
    drop w
}
