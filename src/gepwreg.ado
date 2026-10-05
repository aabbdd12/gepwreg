*! gepwreg.ado  v1.5.0  05oct2026  Araar A.
*! 1.4 (September 2026).  One release since 1.3; the numbers in between were
*! development states that were never distributed.
*! Corrected 04oct2026, results unchanged: every Mata function and structure
*! carries the prefix _gepwreg_, so that a package defining the same names
*! (msat's private PWR engine did: _wls, _se_IF...) cannot replace them.
*! Corrected 04oct2026, results unchanged for every weight it accepts:
*! fweights are refused (they were treated as sampling weights, so their
*! standard errors were not those of frequency weights); the bootstrap
*! under a survey design is labelled as what it is, PSUs within strata;
*! an fpc() declared by svyset is said to be ignored, and a stratum with a
*! single PSU to be left out of the Taylor variance.
*! 1.5.0.  Review C8: a BINARY regressor that is also a
*! heterogeneity variable gets the 1-vs-0 difference as its household effect.
*! The derivative 2 d x of the term d x^2 is not identified when x^2 = x (the
*! coefficients of x, z and x z are then a generalized-inverse choice); the
*! difference, whose own term is d, is.  In the estimate, the influence
*! function, the decomposition and initial().
*! ref(): the effect for a REFERENCE household -- a real type by default
*! (continuous regressors at their mean, a 0/1 variable at its most frequent
*! value, a factor at its most frequent level) or values given, ref(size=5
*! urban=1).  Under het(z) it is M(xref) gamma, with its exact influence
*! function, and does not depend on tau.  Under het(qr) it is beta(u*), the
*! slopes of the conditional quantile regression at the rank u* that puts the
*! reference household at q_tau (the grid path inverted); its standard error
*! is the delta method (with a warning: it rests on the density of the
*! conditional quantiles) or the bootstrap, the reference re-computed in
*! every draw.
*!
*! THE ESTIMATOR.  The effect at the tau-quantile of the outcome is estimated
*! by a two-step estimator: unit-level effects from a heterogeneity model
*! fitted on the whole sample, averaged with the percentile weights on the
*! rank of the outcome.  The kernel-weighted regression on outcome-ranked
*! units of 1.3 converges to a bandwidth-dependent fraction of the effect as
*! soon as the outcome has an error term (Goldberger 1981; see the help file
*! and Araar 2026, version 3); it is kept as rankdep and labelled
*! descriptive.  New options het(), initial(), xref(), qgrid().  rankvar() is
*! unchanged, and now says why it refuses several variables or a
*! factor-variable term.
*!
*! BANDWIDTH.  The MSE-optimal plug-in of Araar (2026) is the default for
*! every estimator, the second step of the two-step included; the second step
*! is a Nadaraya-Watson regression on a uniform rank, so only the curvature
*! of the profile enters the bias, and one bandwidth per call minimises the
*! average relative MSE over the coefficients.  silverman and band() override
*! it.
*!
*! MEASUREMENT ERROR.  merr corrects the step-1 fit of het(varlist) for
*! classical measurement error in the heterogeneity variables.  One variance
*! per variable, identified from the third moments of the step-1 residual and
*! applied only where t = sigma2/se exceeds tcrit(), default 2.  The influence
*! function accounts for the estimated variance, so the analytical and
*! survey-design standard errors stay valid where the correction is active,
*! and boot() redoes the correction and its threshold in every resample.  The
*! command refuses what it cannot identify: a heterogeneity variable taking
*! two values (z^2 is then affine in z, so the third-moment conditions reduce
*! to moments step 1 has already set to zero), a singular step-1 design, and
*! a corrected moment matrix that is not positive definite.  The table
*! reports, per variable, the implied reliability, the skewness of z~ and the
*! R2 of z~ on the regressors, so a variable that is unidentifiable can be
*! told from one this sample merely fails to resolve; e(merr_diag) holds them.
*!
*! STANDARD ERRORS.  Taylor linearisation by default under svyset, over PSUs
*! within strata; the influence-function variance otherwise, always stored in
*! e(V_IF); a bootstrap that re-runs the whole procedure.  The scores of
*! rankvar and rankdep are centred (the constant of the indirect term was
*! omitted up to 1.3, a difference of second order).  Standard errors no
*! longer depend on the order of the observations: the rank term of the
*! influence function is summed over the set of units above i rather than
*! over the rows below it, which the estimates never needed and the variances
*! did.
*!
*! WEIGHTS.  One rule, stated in the header of the output.  A weight declared
*! by svyset is captured first and used for everything -- coefficients,
*! bandwidth, correction, rank, variance -- and a weight on the command line
*! is then ignored with a notice; the command weight is used only when svyset
*! declares none.  e(wtype), e(wexp) and e(wsrc) record the source.
*!
*! ALSO.  Factor variables: every non-base level gets its own indicator.
*! gepwreg_setable, a post-estimation comparison of the standard errors, in
*! its own file and available as the setable option.  A rewritten dialog box
*! covering every option.  Progress dots on all four bootstrap paths, with
*! nodots to suppress them.  Execution time in the header and in e(etime).
*! Percentile Weights Regression
*! Araar (2016, 2023, 2026) ; Firpo, Fortin & Lemieux (2009) ; Goldberger (1981) ;
*! Deville (1999) ; Newey & McFadden (1994)
/* ─────────────────────────────────────────────────────────────────────────────
   Syntax:
     gepwreg depvar indepvars [fw/aw/pw] [if] [in],
           [ PERcentile(#) HET(varlist | qr) INITial(varlist | _all) XREF(numlist) QGRID(numlist)
             RANKvar(varname) RANKdep CBAND(#) BAND(#) SILverman OPTbw
             BOOT(#) SEED(#) noDOTs VCE(svy|if) SETable Level(#)
             noConstant ]
   Modes:
     (default)      two-step estimator of the effect at the tau-quantile of y,
                    heterogeneity model = conditional quantile regression (het(qr))
     het(z1 z2 ...) two-step, heterogeneity model = interactions of x with z
     initial(vars)  two-step, effect at the tau-quantile of y0 = y - b_i (x - xref),
                    x restricted to the listed regressors (_all: every regressor)
     rankvar(z)     profile of the effect along the rank of a covariate z
     rankdep        the one-step weighted regression on outcome-ranked units
                    (version 1.3): descriptive, not the effect
   ───────────────────────────────────────────────────────────────────────── */
#delimit ;

/* ── gepwreg: the command (1.5.0).  Without graph() it hands its arguments,
   unchanged, to _gepwreg_est, the estimator, so every result is that of 1.4
   (tests/lock/test_lock_1_4.do) except where review C8 corrects it (a binary
   regressor also in het(): specification 9 of that lock, exempted there).  With graph(coefficients [, percentiles])
   it runs the estimation at per() as usual, then again, quietly, at every
   percentile of the grid (default 0.05(0.025)0.95), collects the named
   coefficients with their confidence intervals in e(profile), gives the
   estimation at per() back in e(), and draws one panel per coefficient,
   combined into one figure (gepwreg_profile).                           */
capture program drop gepwreg ;
program define gepwreg, eclass properties(svyb) ;
version 16.0 ;
syntax varlist(min=1 numeric fv) [fw aw pw] [if] [in] [,
    PERcentile(real 0.5) GRAPH(string asis) Level(cilevel) * ] ;
if `"`graph'"' == "" {;
    _gepwreg_est `0' ;
    exit ;
} ;
local wt ;
if "`weight'" != "" local wt "[`weight'`exp']" ;

/* the coefficients, the grid and the keywords of graph(): ref draws the curve
   of the reference household (the ref() of the command, or ref(default)),
   uniform adds the simultaneous band */
_parse comma glist grest : graph ;
local grest = trim(subinstr(`"`grest'"', ",", "", 1)) ;
local g_ref 0 ;
local g_uni 0 ;
local gnum ;
foreach tok of local grest {;
    if      lower("`tok'") == "ref"     local g_ref 1 ;
    else if lower("`tok'") == "uniform" local g_uni 1 ;
    else local gnum `gnum' `tok' ;
} ;
if `"`gnum'"' == "" local gnum "0.05(0.025)0.95" ;
capture numlist `"`gnum'"', range(>0 <1) sort ;
if _rc {;
    di as err "graph(): the percentiles must be a numlist strictly between 0 and 1," ;
    di as err "        for example graph(size male, 0.1(0.05)0.9 ref uniform)" ;
    exit 198 ;
} ;
local grid `r(numlist)' ;
local ng : word count `grid' ;
local refopt ;
if `g_ref' & !regexm(`" `options'"', "[ ]ref[(]") local refopt "ref(default)" ;

/* the estimation at per(), displayed as usual, with the profile across the
   grid computed inside it (two-step estimators: one bootstrap loop for every
   percentile, or the joint influence function) */
_gepwreg_est `varlist' `wt' `if' `in', percentile(`percentile') level(`level') `options' `refopt' profile(`grid') ;
capture confirm matrix e(prof_b) ;
local engine = (_rc == 0) ;
if !`engine' & (`g_ref' | `g_uni') {;
    di as text "Note: graph(): ref and uniform need the two-step estimator without initial() or merr --" ;
    di as text "      the figure shows the pointwise intervals only." ;
    local g_ref 0 ;
    local g_uni 0 ;
} ;
local cn : colnames e(b) ;
local coefs ;
foreach t of local glist {;
    if substr("`t'", 1, 2) == "i." {;
        local v = substr("`t'", 3, .) ;
        local nf 0 ;
        foreach c of local cn {;
            if regexm("`c'", "^[0-9]+\.`v'$") {;
                local coefs `coefs' `c' ;
                local ++nf ;
            } ;
        } ;
        if !`nf' {;
            di as err "graph(): `t' has no level among the coefficients: `cn'" ;
            exit 111 ;
        } ;
    } ;
    else {;
        local here : list t in cn ;
        if !`here' {;
            di as err "graph(): `t' is not a coefficient of the model; the coefficients are:" ;
            di as err "        `cn'" ;
            exit 111 ;
        } ;
        local coefs `coefs' `t' ;
    } ;
} ;
local coefs : list uniq coefs ;
local k : word count `coefs' ;
if `k' == 0 {;
    di as err "graph(): name the coefficients to draw, for example graph(size male)" ;
    exit 198 ;
} ;

/* the profile.  Engine: e(prof_*) of the estimation, columns tau, then for
   each coefficient b se lb ub (pointwise), ulb uub (simultaneous band),
   rb rse rlb rub rulb ruub (the reference household)                    */
tempname hold P zq ;
scalar `zq' = invnormal((100 + `level')/200) ;
local nb = 12 ;
matrix `P' = J(`ng', 1 + `nb'*`k', .) ;
local names "tau" ;
forvalues j = 1/`k' {;
    local names `names' b`j' se`j' lb`j' ub`j' ulb`j' uub`j' rb`j' rse`j' rlb`j' rub`j' rulb`j' ruub`j' ;
} ;
matrix colnames `P' = `names' ;
local noV 0 ;
local hasR 0 ;
local hasRS 0 ;
local hasRF 0 ;
if `engine' {;
    tempname PB PS PC PF RB RS RC RF ;
    matrix `PB' = e(prof_b) ;
    capture matrix `PS' = e(prof_se) ;
    if _rc local noV 1 ;
    else if missing(`PS'[1, 1]) local noV 1 ;
    capture matrix `PC' = e(prof_cu) ;
    capture matrix `PF' = e(prof_flat) ;
    if `g_ref' {;
        capture matrix `RB' = e(prof_rb) ;
        local hasR = (_rc == 0) ;
        capture matrix `RS' = e(prof_rse) ;
        local hasRS = (_rc == 0) ;
        if `hasRS' local hasRS = !missing(`RS'[1, 1]) ;
        capture matrix `RC' = e(prof_rcu) ;
        capture matrix `RF' = e(prof_rflat) ;
        local hasRF = (_rc == 0) ;
    } ;
    local j 0 ;
    foreach c of local coefs {;
        local ++j ;
        local pc = colnumb(`PB', "`c'") ;
        local flat`j' = . ;
        local rflat`j' = . ;
        if !`noV' local flat`j' = `PF'[2, `pc'] ;
        if `hasR' & `hasRF' local rflat`j' = `RF'[2, `pc'] ;
        forvalues i = 1/`ng' {;
            local o = 1 + `nb'*(`j' - 1) ;
            matrix `P'[`i', 1] = `: word `i' of `grid'' ;
            matrix `P'[`i', `o' + 1] = `PB'[`i', `pc'] ;
            if !`noV' {;
                local sv = `PS'[`i', `pc'] ;
                matrix `P'[`i', `o' + 2] = `sv' ;
                matrix `P'[`i', `o' + 3] = `PB'[`i', `pc'] - `zq'*`sv' ;
                matrix `P'[`i', `o' + 4] = `PB'[`i', `pc'] + `zq'*`sv' ;
                matrix `P'[`i', `o' + 5] = `PB'[`i', `pc'] - `PC'[1, `pc']*`sv' ;
                matrix `P'[`i', `o' + 6] = `PB'[`i', `pc'] + `PC'[1, `pc']*`sv' ;
            } ;
            if `hasR' {;
                matrix `P'[`i', `o' + 7] = `RB'[`i', `pc'] ;
                if `hasRS' {;
                    local rv = `RS'[`i', `pc'] ;
                    matrix `P'[`i', `o' + 8]  = `rv' ;
                    matrix `P'[`i', `o' + 9]  = `RB'[`i', `pc'] - `zq'*`rv' ;
                    matrix `P'[`i', `o' + 10] = `RB'[`i', `pc'] + `zq'*`rv' ;
                    matrix `P'[`i', `o' + 11] = `RB'[`i', `pc'] - `RC'[1, `pc']*`rv' ;
                    matrix `P'[`i', `o' + 12] = `RB'[`i', `pc'] + `RC'[1, `pc']*`rv' ;
                } ;
            } ;
        } ;
    } ;
} ;
else {;
_estimates hold `hold', copy restore ;
di as text _n "Profile across `ng' percentiles " _c ;
local i 0 ;
foreach t of local grid {;
    local ++i ;
    capture quietly _gepwreg_est `varlist' `wt' `if' `in', percentile(`t') level(`level') `options' ;
    if _rc {;
        local rc = _rc ;
        di as err _n "graph(): the estimation at percentile `t' failed, r(`rc')" ;
        exit `rc' ;
    } ;
    matrix `P'[`i', 1] = `t' ;
    capture confirm matrix e(V) ;
    local hasV = (_rc == 0) ;
    if !`hasV' local noV 1 ;
    local j 0 ;
    foreach c of local coefs {;
        local ++j ;
        local o = 1 + `nb'*(`j' - 1) ;
        matrix `P'[`i', `o' + 1] = _b[`c'] ;
        if `hasV' {;
            matrix `P'[`i', `o' + 2] = _se[`c'] ;
            matrix `P'[`i', `o' + 3] = _b[`c'] - `zq'*_se[`c'] ;
            matrix `P'[`i', `o' + 4] = _b[`c'] + `zq'*_se[`c'] ;
        } ;
    } ;
    di as text "." _c ;
} ;
di "" ;
_estimates unhold `hold' ;
forvalues j = 1/`k' {;
    local flat`j' = . ;
    local rflat`j' = . ;
} ;
} ;
ereturn matrix profile = `P', copy ;
ereturn local profile_coefs "`coefs'" ;
ereturn local cmdline "gepwreg `0'" ;

/* one panel per coefficient, combined; the titles are read on the data */
local dep : word 1 of `varlist' ;
local j 0 ;
foreach c of local coefs {;
    local ++j ;
    local ttl`j' "`c'" ;
    if regexm("`c'", "^([0-9]+)\.(.+)$") {;
        local lev = regexs(1) ;
        local v   = regexs(2) ;
        local vlab : variable label `v' ;
        if "`vlab'" == "" local vlab "`v'" ;
        local vl : value label `v' ;
        local llab "`lev'" ;
        if "`vl'" != "" local llab : label `vl' `lev' ;
        local ttl`j' "`v' = `llab'" ;
    } ;
    else {;
        capture local vlab : variable label `c' ;
        if _rc == 0 & "`vlab'" != "" local ttl`j' "`vlab'" ;
    } ;
    if length(`"`ttl`j''"') > 45 {;
        local ttl`j' = substr(`"`ttl`j''"', 1, 42) + "..." ;
    } ;
} ;
local gl ;
preserve ;
quietly {;
    clear ;
    svmat double `P', names(col) ;
    generate double zero = 0 ;
} ;
local j 0 ;
foreach c of local coefs {;
    local ++j ;
    local band ;
    local show b`j' ;
    if !`noV' {;
        if `g_uni' {;
            local band `band' (rarea ulb`j' uub`j' tau, color(navy%12) lwidth(none)) ;
            local show `show' ulb`j' uub`j' ;
        } ;
        local band `band' (rarea lb`j' ub`j' tau, color(navy%28) lwidth(none)) ;
        local show `show' lb`j' ub`j' ;
    } ;
    local rline ;
    if `hasR' {;
        /* the reference household: its curve dashed, its interval as two
           dotted lines (no second shaded band over the first)               */
        local rline (line rb`j' tau, lcolor(maroon) lpattern(dash)) ;
        local show `show' rb`j' ;
        quietly count if !missing(rse`j') ;
        if r(N) > 0 {;
            local rl = cond(`g_uni', "rulb`j'", "rlb`j'") ;
            local ru = cond(`g_uni', "ruub`j'", "rub`j'") ;
            local rline `rline' (line `rl' tau, lcolor(maroon%70) lpattern(shortdash) lwidth(thin))
                                (line `ru' tau, lcolor(maroon%70) lpattern(shortdash) lwidth(thin)) ;
            local show `show' `rl' `ru' ;
        } ;
    } ;
    /* the zero line, a plotted series rather than yline() (an out-of-range
       yline() is drawn at a wrong height once the panels are combined),
       and only where zero is inside the panel, so as not to flatten a
       profile that stays away from it */
    local ylo = . ;
    local yhi = . ;
    foreach sv of local show {;
        quietly summarize `sv' ;
        if r(N) > 0 {;
            local ylo = min(`ylo', r(min)) ;
            local yhi = max(`yhi', r(max)) ;
        } ;
    } ;
    local zero ;
    if `ylo' <= 0 & `yhi' >= 0 local zero (line zero tau, lpattern(dash) lcolor(gs10)) ;
    /* the p-values of the flatness tests, under the title */
    local sub ;
    if `flat`j'' < . local sub = "constant across tau: p = " + strtrim(string(`flat`j'', "%5.3f")) ;
    if `rflat`j'' < . local sub = "`sub'" + "    reference: p = " + strtrim(string(`rflat`j'', "%5.3f")) ;
    twoway `band' `zero' (line b`j' tau, lcolor(navy) lpattern(solid)) (scatter b`j' tau, msymbol(square) msize(vsmall) mcolor(navy)) `rline',
        xline(`percentile', lpattern(shortdash) lcolor(gs8))
        title(`"`ttl`j''"', size(medsmall)) subtitle(`"`sub'"', size(vsmall))
        xtitle("Percentile {&tau}", size(small)) ytitle("")
        ylabel(, angle(0) labsize(small)) xlabel(0(0.2)1, labsize(small))
        legend(off) name(_gepwreg_p`j', replace) nodraw ;
    local gl `gl' _gepwreg_p`j' ;
} ;
restore ;
local cols = cond(`k' <= 2, `k', cond(`k' <= 4, 2, 3)) ;
local nt "Effect at the percentiles of `dep': `ng' percentiles from `: word 1 of `grid'' to `: word `ng' of `grid''." ;
/* the notes, one short line each */
local notes `"`"`nt'"'"' ;
if `noV' local notes `"`notes' `"No confidence interval: this estimator has no analytical standard errors -- add boot(#)."'"' ;
else {;
    local notes `"`notes' `"Band: pointwise `level'% intervals, `e(SE_type)', the same rule at every percentile."'"' ;
    if `g_uni' local notes `"`notes' `"Light band: the simultaneous `level'% band over the percentiles (sup-t)."'"' ;
} ;
if `engine' & !`noV' local notes `"`notes' `"p: the test that the effect is constant across the percentiles of the figure (sup-t)."'"' ;
if `hasR' {;
    local which = cond(`g_uni', "its simultaneous band", "its pointwise interval") ;
    if `hasRS' local notes `"`notes' `"Dashed (maroon): the reference household, dotted: `which'."'"' ;
    else       local notes `"`notes' `"Dashed (maroon): the reference household."'"' ;
    /* the reference household, cut at its separators */
    local rest `"`e(ref_text)'"' ;
    local line ;
    while `"`rest'"' != "" {;
        local p = strpos(`"`rest'"', char(59)) ;
        if `p' == 0 {;
            local piece `"`rest'"' ;
            local rest ;
        } ;
        else {;
            local piece = substr(`"`rest'"', 1, `p') ;
            local rest  = strtrim(substr(`"`rest'"', `p' + 1, .)) ;
        } ;
        if strlen(`"`line' `piece'"') > 95 & `"`line'"' != "" {;
            local notes `"`notes' `"  `line'"'"' ;
            local line `"`piece'"' ;
        } ;
        else local line = strtrim(`"`line' `piece'"') ;
    } ;
    if `"`line'"' != "" local notes `"`notes' `"  `line'"'"' ;
    /* het(qr): where the reference household would need a conditional rank
       outside the grid of quantile regressions, its curve is held at the
       end of the grid -- say at which percentiles (the rows of the figure's
       grid; the last row of e(prof_ru) is per())                         */
    tempname RU ;
    capture matrix `RU' = e(prof_ru) ;
    if !_rc {;
        foreach s in 1 -1 {;
            local n0 0 ;
            local i0 ;
            local i1 ;
            forvalues i = 1/`ng' {;
                if `RU'[`i', 2] == `s' {;
                    local ++n0 ;
                    if "`i0'" == "" local i0 `i' ;
                    local i1 `i' ;
                    local uend = strtrim(string(`RU'[`i', 1], "%9.0g")) ;
                } ;
            } ;
            if `n0' {;
                local end = cond(`s' == 1, "the last quantile of the grid", "the first quantile of the grid") ;
                local t0 : word `i0' of `grid' ;
                local t1 : word `i1' of `grid' ;
                if `n0' == 1                    local at "at tau = `t0'" ;
                else if `i1' - `i0' + 1 == `n0' local at "for tau from `t0' to `t1'" ;
                else                            local at "at `n0' percentiles from `t0' to `t1'" ;
                local notes `"`notes' `"  Its curve is held at u = `uend', `end', `at': its rank would leave the grid."'"' ;
            } ;
        } ;
    } ;
} ;
local notes `"`notes' `"Dotted vertical line: per()."'"' ;
graph combine `gl', cols(`cols') name(gepwreg_profile, replace)
    note(`notes', size(vsmall)) ;
graph drop `gl' ;
end ;


/* ── _gepwreg_est: the estimator (gepwreg up to 1.4) ──────────────────── */
capture program drop _gepwreg_est ;
program define _gepwreg_est, eclass sortpreserve ;
version 16.0 ;

syntax varlist(min=1 numeric fv) [fw aw pw] [if] [in] [,
    PERcentile(real 0.5)
    HET(string)
    INITial(string)
    XREF(numlist)
    QGRID(numlist >0 <1 sort)
    RANKdep
    RANKvar(string)
    CBAND(real 0.9)
    BAND(real 0)
    OPTbw
    MERR
    TCRIT(real 2)
    SILverman
    BOOT(integer 0)
    SEED(integer 12345)
    REF(string)
    PROFile(numlist >0 <1 sort)
    noDOTs
    Level(cilevel)
    noConstant
    VCE(string)
    SETable
    /* Internal options passed by svy: prefix */
    SVY
    noHEader
    noTABle
] ;

/* (Corrected 04oct2026) fweights are refused: a unit of the survey stands
   for its sampling weight, it is not a replicated record -- the influence-
   function variance and the bootstrap both count the units.              */
if "`weight'" == "fweight" {;
    di as err "fweights are not allowed: sampling weights go in {bf:[pweight=]};" ;
    di as err "a frequency weight built from a sampling weight, such as int(pw*10000), is a pweight;" ;
    di as err "replicated records (true frequencies) can be expanded first: {bf:expand} {it:fwvar}" ;
    exit 101 ;
} ;

/* ── 0a. Estimation mode ────────────────────────────────────────── */
local mode "twostep" ;
local het_type "qr" ;
local hetvars "" ;
if "`rankdep'" != "" local mode "rankdep" ;
if "`rankvar'" != "" {;
    /* (1.4) rankvar() was declared varname, so anything else gave a bare
       "invalid syntax".  The profile is drawn along the RANK of this
       variable, which needs one ordered scale: several variables have no
       joint rank, and a factor-variable term such as i.zone is a set of
       indicators, not an order.  Say so.                              */
    capture unab rankvar : `rankvar' ;
    local rv_n : word count `rankvar' ;
    capture confirm numeric variable `rankvar' ;
    if _rc | `rv_n' != 1 {;
        di as error "rankvar() takes one existing numeric variable." ;
        di as error "The profile is drawn along the RANK of that variable, so it needs" ;
        di as error "a single ordered scale.  Several variables have no joint rank, and" ;
        di as error "a factor-variable term (i.var) is a set of indicators, not an order." ;
        di as error "For a categorical variable, use its numeric code: rankvar(var)." ;
        exit 198 ;
    } ;
    if "`mode'" == "rankdep" {;
        di as error "rankdep and rankvar() cannot be combined" ;
        exit 198 ;
    } ;
    local mode "rankvar" ;
} ;
if "`het'" != "" {;
    if "`mode'" != "twostep" {;
        di as error "het() cannot be combined with rankdep or rankvar()" ;
        exit 198 ;
    } ;
    if lower(trim("`het'")) == "qr" local het_type "qr" ;
    else {;
        local het_type "z" ;
        fvunab hetvars : `het' ;
    } ;
} ;
local do_initial = ("`initial'" != "") ;
if `do_initial' & "`mode'" != "twostep" {;
    di as error "initial requires the two-step estimator (no rankdep, no rankvar())" ;
    exit 198 ;
} ;
if "`xref'" != "" & !`do_initial' {;
    di as text "Note: xref() is used only with initial; ignored." ;
} ;
/* (1.5.0) ref(): the effect for a reference household */
local do_ref = (`"`ref'"' != "") ;
if `do_ref' & "`mode'" != "twostep" {;
    di as error "ref() requires the two-step estimator, het(qr) or het(varlist)" ;
    exit 198 ;
} ;
if `do_ref' & `do_initial' {;
    di as error "ref() cannot be combined with initial()" ;
    exit 198 ;
} ;
/* (1.5.0) profile(): the estimates across the percentiles of graph(),
   computed inside this one estimation (an internal option of gepwreg).
   Two-step estimators without initial() or merr, otherwise graph() loops. */
local do_prof = ("`profile'" != "") ;
if `do_prof' & ("`mode'" != "twostep" | `do_initial' | "`merr'" != "") local do_prof = 0 ;
local nuse : word count `profile' ;

/* ── 0. Mark sample ──────────────────────────────────────────────────────── */
marksample touse ;
if "`hetvars'" != "" {;
    fvrevar `hetvars' if `touse' ;
    markout `touse' `r(varlist)' ;
} ;
qui count if `touse' ;
if r(N) == 0 error 2000 ;

/* ── 1. Parse varlist ────────────────────────────────────────────────────── */
local depvar    : word 1 of `varlist' ;
local indepvars : list varlist - depvar ;

/* Step 1 : fvexpand, then create the temporary variables from the FULL
   expanded list, base levels included.  (1.4) Passing only the non-base
   levels to fvrevar -- 2.gse 3.gse ... 7.gse -- makes Stata take the lowest
   listed level as the base of that set, so its indicator came out as a
   column of zeros: level 2 was silently pooled with level 1 and reported as
   "omitted".  With 1b.gse in the list every level keeps its own indicator;
   the base and omitted terms are then dropped from BOTH lists in parallel,
   so that the display names and the Mata columns always line up.          */
if "`indepvars'" != "" {;
    fvexpand `indepvars' if `touse' ;
    local indepvars_exp_raw `r(varlist)' ;
    fvrevar `indepvars_exp_raw' if `touse' ;
    local indepvars_mata_raw `r(varlist)' ;
    if (`: word count `indepvars_mata_raw'' != `: word count `indepvars_exp_raw'') {;
        di as error "gepwreg: fvrevar returned a list of another length than fvexpand" ;
        exit 498 ;
    } ;
    local indepvars_exp ;
    local indepvars_mata ;
    local j = 0 ;
    foreach v of local indepvars_exp_raw {;
        local ++j ;
        local t : word `j' of `indepvars_mata_raw' ;
        /* a base (Nb.) or omitted (No.) level anywhere in the term, an
           interaction included; Nbn. (no base) is kept                      */
        if !regexm("`v'","[0-9]+b\.") & !regexm("`v'","[0-9]+o\.") {;
            local indepvars_exp  `indepvars_exp' `v' ;
            local indepvars_mata `indepvars_mata' `t' ;
        } ;
    } ;
} ;
else {;
    local indepvars_exp "" ;
    local indepvars_mata "" ;
} ;

/* ── 1b. Ranking variable ───────────────────────────────────────────────── */
local rank_var "`depvar'" ;
local rank_label "`depvar' (outcome)" ;
if "`mode'" == "rankvar" {;
    local rank_var "`rankvar'" ;
    local rank_label "`rankvar'" ;
    markout `touse' `rankvar' ;
    di as text "Ranking variable : `rankvar' (profile of the effect along its rank)" ;
} ;

/* ── 2. Survey weights ───────────────────────────────────────────────────── */
/* Priority rule, fixed 18sep2026 and stated to the user :                    */
/*   1. svyset declares a weight  ->  it is captured HERE, before anything   */
/*      is computed, and used for EVERYTHING: the step-1 coefficients, the   */
/*      bandwidth, the measurement-error correction, the rank, the variance. */
/*      A weight expression on the command is then ignored, with a notice.   */
/*      Two declarations that can disagree would otherwise estimate with one */
/*      weight and build the design variance for the plan of another.        */
/*   2. svyset declares no weight ->  the command weight is used.            */
/*   3. Neither                   ->  unweighted (fw = 1).                   */
tempvar fw_var ;
qui gen double `fw_var' = 1 if `touse' ;
local wgt_name "(unweighted)" ;
local wgt_type "" ;
local wgt_exp  "" ;
local wexpr    "" ;
if "`exp'" != "" {;
    local wexpr = subinstr("`exp'","=","",1) ;
    local wexpr = trim("`wexpr'") ;
} ;
qui svyset ;
local svywvar = trim("`r(wvar)'") ;
local wsrc "none" ;
if "`svywvar'" != "" {;
    qui replace `fw_var' = `svywvar' if `touse' ;
    local wsrc     "svy" ;
    local wgt_type "pweight" ;
    local wgt_exp  "`svywvar'" ;
    local wgt_name "PW[`svywvar'] (svyset)" ;
    di as text "Note: svyset weights [`svywvar'] are used for the estimation and the s.e." ;
    if "`wexpr'" != "" & "`wexpr'" != "`svywvar'" {;
        di as text "      The weight given on the command [`wexpr'] is ignored." ;
    } ;
} ;
if "`wsrc'" == "none" & "`wexpr'" != "" {;
    qui replace `fw_var' = `wexpr' if `touse' ;
    local wsrc     "cmd" ;
    local wgt_type = lower("`weight'") ;
    local wgt_exp  "`wexpr'" ;
    local wgt_ab = upper(substr("`weight'",1,2)) ;
    local wgt_name "`wgt_ab'[`wexpr']" ;
} ;

/* ── 3. Validate options ─────────────────────────────────────────────────── */
if `percentile' <= 0 | `percentile' >= 1 {;
    di as error "percentile() must be strictly between 0 and 1" ;
    exit 198 ;
} ;
if `boot' < 0 {;
    di as error "boot() must be a non-negative integer" ;
    exit 198 ;
} ;
/* (1.4) Execution time in the header, so that the cost of anything --
   the bootstrap, the dots, a large qgrid -- can be read rather than guessed.
   Stata's timers are a shared numbered resource and there is no private
   high-resolution clock, so gepwreg uses timer 99 and says so; a user who
   keeps something in timer 99 across a gepwreg call will find it cleared. */
timer clear 99 ;
timer on 99 ;
local do_dots = ("`dots'" == "") ;
scalar _gepwreg_dots = `do_dots' ;
local do_merr = ("`merr'" != "") ;
if (`do_merr' & "`het_type'" != "z") {;
    di as error "merr requires het(varlist): the correction is defined for" ;
    di as error "the step-1 design (X, Z, X x Z), not for het(qr)." ;
    exit 198 ;
} ;
local n_bwspec = (`band' != 0) + ("`silverman'" != "") + ("`optbw'" != "") ;
if `n_bwspec' > 1 {;
    di as error "specify at most one of band(), silverman, optbw" ;
    exit 198 ;
} ;
if `band' != 0 {;
    local do_optbw = 0 ;
    local bw_label "Manual" ;
} ;
else if "`silverman'" != "" {;
    local do_optbw = 0 ;
    local bw_label "Silverman" ;
} ;
else {;
    local do_optbw = 1 ;
    local bw_label "MSE-optimal" ;
} ;
/* (1.4) The kernel of the two-step smooths unit-level effects along the rank
   of y, which is a Nadaraya-Watson regression with a uniform design density.
   Its MSE-optimal bandwidth is derived in _gepwreg_hopt2() and is the default
   here as it is for the one-step estimators; silverman or band() override.  */
local addcons = ("`constant'" == "") ;

/* ── 3b. Parse vce() option for survey design ───────────────────────────── */
local do_svy = 0 ;
local svy_psu    "" ;
local svy_strata "" ;
if "`svy'" != "" {;
    local vce "linearized" ;
} ;
if "`vce'" == "" {;
    qui svyset ;
    if "`r(su1)'" != "" | "`r(strata1)'" != "" {;
        local vce "svy" ;
        local svy_auto = 1 ;
    } ;
} ;
else if lower("`vce'")=="if" | lower("`vce'")=="analytic" {;
    local vce "" ;
} ;
if "`vce'" != "" {;
    if lower("`vce'")=="svy" | lower("`vce'")=="survey" | lower("`vce'")=="linearized" {;
        qui svyset ;
        if "`r(wvar)'" == "" & "`r(su1)'" == "" {;
            di as error "vce(svy) requires svyset to be called first" ;
            error 119 ;
        } ;
        local do_svy    = 1 ;
        local svy_psu    "`r(su1)'" ;
        local svy_strata "`r(strata1)'" ;
        local svy_fpc    "`r(fpc1)'" ;
        if "`exp'" == "" & "`r(wvar)'" != "" {;
            qui replace `fw_var' = `r(wvar)' if `touse' ;
            local wgt_name "PW[`r(wvar)'] (svyset)" ;
        } ;
        if ("`svy_auto'" == "1") di as text "Survey design detected (svyset); specify vce(if) for the IF-corrected s.e." ;
        else                      di as text "Survey design detected:" ;
        if "`svy_psu'" == "" di as text "  PSU    : (not declared — variance may be conservative)" ;
        else                 di as text "  PSU    : `svy_psu'" ;
        di as text "  Strata : `svy_strata'" ;
        if "`svy_fpc'" != "" {;
            di as text "  FPC    : `svy_fpc' is not used -- no finite-population correction" ;
            di as text "           (the Taylor variance is conservative)" ;
        } ;
        if "`svy_strata'" != "" & "`svy_psu'" != "" {;
            qui {;
                tempvar _first _npsu ;
                bysort `touse' `svy_strata' `svy_psu' : gen byte `_first' = (_n == 1) & `touse' ;
                bysort `touse' `svy_strata' (`svy_psu') : gen `_npsu' = sum(`_first') ;
                bysort `touse' `svy_strata' : replace `_npsu' = `_npsu'[_N] ;
            } ;
            qui sum `_npsu' if `touse' ;
            local min_psu = r(min) ;
            if `min_psu' < 2 {;
                di as text "Warning: a stratum has a single PSU in the estimation sample;" ;
                di as text "  its contribution to the Taylor variance is not defined, and it" ;
                di as text "  is left out of it (it contributes zero)." ;
            } ;
            else if `min_psu' < 5 {;
                di as text "Warning: minimum PSU/stratum = `min_psu'." ;
                di as text "  Taylor SE may be unreliable. Recommend n_h >= 10." ;
            } ;
        } ;
    } ;
    else {;
        di as error "vce(`vce') not supported. Use vce(svy) or vce(if)." ;
        exit 198 ;
    } ;
} ;

/* ═════════════════════════════════════════════════════════════════════════
   T. TWO-STEP ESTIMATOR: het(z) or het(qr), and initial
   ═════════════════════════════════════════════════════════════════════════ */
if "`mode'" == "twostep" {;
    foreach m in _gepwreg_b _gepwreg_batq _gepwreg_V _gepwreg_Vi _gepwreg_Vs _gepwreg_Vb
                 _gepwreg_g _gepwreg_zbar _gepwreg_zpop _gepwreg_decomp _gepwreg_G _gepwreg_ugrid
                 _gepwreg_same _gepwreg_iflag
                 _gepwreg_rk _gepwreg_rfac _gepwreg_rx _gepwreg_rz _gepwreg_xfix _gepwreg_zfix
                 _gepwreg_bref _gepwreg_refx _gepwreg_refu _gepwreg_Vref_IF _gepwreg_Vref_svy
                 _gepwreg_Vref_boot _gepwreg_Erow
                 _gepwreg_pgrid _gepwreg_pb _gepwreg_pse _gepwreg_pcu _gepwreg_pflat
                 _gepwreg_prb _gepwreg_prse _gepwreg_prcu _gepwreg_prflat _gepwreg_Prow
                 _gepwreg_pru {;
        capture matrix drop `m' ;
    } ;
    /* the grid of the profile, per() appended: its bootstrap gives the
       standard errors of the table at per() as well                       */
    if `do_prof' mata: st_matrix("_gepwreg_pgrid", (strtoreal(tokens("`profile'")), `percentile')) ;
    local kx : word count `indepvars_mata' ;
    if `kx' == 0 {;
        di as error "the two-step estimator needs at least one regressor" ;
        exit 198 ;
    } ;
    if !`addcons' di as text "Note: noconstant does not apply to the two-step estimator (step 1 keeps its constant)." ;

    /* reference values of x for the initial-quantile object */
    if "`xref'" == "" {;
        local xref_list ;
        forvalues j = 1/`kx' {; local xref_list `xref_list' 0 ; } ;
    } ;
    else {;
        local kfull_x : word count `indepvars_exp_raw' ;
        if `: word count `xref'' == `kfull_x' & `kfull_x' > `kx' {;
            /* one value per term of the coefficient table, base levels included:
               drop the values that sit on base or omitted levels             */
            local xref_list ;
            local j = 0 ;
            foreach v of local indepvars_exp_raw {;
                local ++j ;
                if !regexm("`v'","[0-9]+b\.") & !regexm("`v'","[0-9]+o\.") {;
                    local xref_list `xref_list' `: word `j' of `xref'' ;
                } ;
            } ;
        } ;
        else if `: word count `xref'' == `kx' local xref_list `xref' ;
        else {;
            di as error "xref() must give one value per regressor: `kx' (non-base terms) or `kfull_x' (all terms of the table)" ;
            exit 198 ;
        } ;
    } ;
    mata: st_matrix("_gepwreg_xref", strtoreal(tokens("`xref_list'"))) ;
    /* which regressors define y0: initial(_all) or a list of regressors      */
    local iflag_list ;
    if `do_initial' {;
        if lower(trim("`initial'")) == "_all" | lower(trim("`initial'")) == "all" {;
            foreach xv of local indepvars_exp {; local iflag_list `iflag_list' 1 ; } ;
        } ;
        else {;
            fvexpand `initial' if `touse' ;
            local init_exp `r(varlist)' ;
            local nfound = 0 ;
            foreach xv of local indepvars_exp {;
                if `: list xv in init_exp' {; local iflag_list `iflag_list' 1 ; local ++nfound ; } ;
                else                        local iflag_list `iflag_list' 0 ;
            } ;
            if `nfound' == 0 {;
                di as error "initial(): none of the listed variables is a regressor of the model" ;
                exit 198 ;
            } ;
        } ;
    } ;
    else {;
        foreach xv of local indepvars_exp {; local iflag_list `iflag_list' 0 ; } ;
    } ;
    mata: st_matrix("_gepwreg_iflag", strtoreal(tokens("`iflag_list'"))) ;

    if `do_initial' & `boot' == 0 {;
        local boot = 200 ;
        di as text "Note: initial() has no analytical standard errors, so the bootstrap" ;
        di as text "      of the whole procedure is used (PSUs within strata under a design);" ;
        di as text "      boot(200) assumed -- specify boot(#) to change it." ;
    } ;

    scalar _gepwreg_h_tmp = 0 ;
    local h0 = . ;
    local neff0 = . ;
    local n_grid = 0 ;

    /* ── T1. Heterogeneity model on observables: closed form + IF ─────── */
    if "`het_type'" == "z" {;
        fvexpand `hetvars' if `touse' ;
        local het_exp_raw `r(varlist)' ;
        fvrevar `het_exp_raw' if `touse' ;
        local het_mata_raw `r(varlist)' ;
        local het_exp ;
        local het_mata ;
        local j = 0 ;
        foreach v of local het_exp_raw {;
            local ++j ;
            local t : word `j' of `het_mata_raw' ;
            if !regexm("`v'","[0-9]+b\.") & !regexm("`v'","[0-9]+o\.") {;
                local het_exp  `het_exp' `v' ;
                local het_mata `het_mata' `t' ;
            } ;
        } ;
        local mz : word count `het_mata' ;
        /* (1.4) same[j,l] = 1 when the heterogeneity variable z_l is the
           regressor x_j itself: the effect of x_j is then the derivative of
           the step-1 model, which also runs through the coefficient of z_l
           and through the interactions of the other regressors with z_l   */
        local same_list ;
        foreach xv of local indepvars_exp {;
            foreach zv of local het_exp {;
                if "`xv'" == "`zv'" local same_list `same_list' 1 ;
                else                local same_list `same_list' 0 ;
            } ;
        } ;
        /* (1.4) merr refuses a heterogeneity variable with two values.
           If z takes only a and b then z^2 = (a+b) z - ab, so the moments
           that identify sigma^2 -- E[e z~^2] and E[e x z~^2] -- are linear
           combinations of E[e], E[e z], E[e x] and E[e x z], each of which
           step 1 has already set to zero.  sigma^2 is not identified: the
           search returns a numerical residual, the reported t is a ratio of
           two such residuals (t = 5.11 was observed on a sigma^2 printed as
           0.000000), and the subtraction of Omega can take M - Omega out of
           positive definiteness, which reaches the user as "estimates post:
           matrix has missing values".  Error in a categorical variable is
           misclassification, which is not classical, so this is a scope
           limit and not a tolerance to widen.                            */
        if `do_merr' {;
            /* (1.5.0) a function of heterogeneity variables in het() -- an
               interaction or a power, c.z#c.z -- is not a variable observed
               with its own classical error: refuse with the reason, rather
               than let the correction fail on its moment matrix              */
            foreach t of local het_exp {;
                if strpos("`t'", "#") {;
                    di as error "merr: het() contains `t', a function of heterogeneity variables." ;
                    di as error "      The correction treats each het() variable as observed with its own" ;
                    di as error "      classical error, which a product or a power of them is not." ;
                    di as error "      Fit this enriched model without merr: it is the check of merr" ;
                    di as error "      (help gepwreg, option merr)." ;
                    exit 198 ;
                } ;
            } ;
            local bad_bin ;
            local jj = 0 ;
            foreach tv of local het_mata {;
                local ++jj ;
                qui su `tv' if `touse', meanonly ;
                local _mn = r(min) ;
                local _mx = r(max) ;
                qui count if `touse' & `tv' > `_mn' & `tv' < `_mx' ;
                if r(N) == 0 {;
                    local nm : word `jj' of `het_exp' ;
                    local bad_bin `bad_bin' `nm' ;
                } ;
            } ;
            if "`bad_bin'" != "" {;
                di as error "merr: the heterogeneity variable(s) `bad_bin' take at most two" ;
                di as error "values.  sigma^2 is not identified for such a variable: z^2 is" ;
                di as error "then affine in z, so the third-moment conditions reduce to" ;
                di as error "moments that step 1 has already set to zero." ;
                di as error "Measurement error in a categorical variable is misclassification," ;
                di as error "whose error is correlated with the true value and is not covered" ;
                di as error "by this correction.  Drop merr, or keep only continuous" ;
                di as error "variables in het()." ;
                exit 198 ;
            } ;
        } ;
        mata: st_matrix("_gepwreg_same", rowshape(strtoreal(tokens("`same_list'")), `kx')) ;
        if `do_ref' {;
            _gepwreg_refrule, spec(`"`ref'"') touse(`touse')
                xexp(`indepvars_exp') xmata(`indepvars_mata') xraw(`indepvars_exp_raw')
                zexp(`het_exp') zmata(`het_mata') zraw(`het_exp_raw') ;
            local ref_nu = r(nu) ;
            local ref_unames `r(unames)' ;
            local ref_kinds  `r(kinds)' ;
            forvalues u = 1/`ref_nu' {; local ref_lev`u' `r(lev`u')' ; } ;
        } ;
        mata: _gepwreg_ts_z("`depvar'", "`indepvars_mata'", "`het_mata'", "`fw_var'", "`touse'",
                            `percentile', `cband', `band', `boot', `do_svy',
                            "`svy_psu'", "`svy_strata'", `do_initial', "_gepwreg_xref", `seed',
                            "_gepwreg_same", "_gepwreg_iflag", `do_optbw',
                            `do_merr', `tcrit', `do_ref', `do_prof', `nuse', `level') ;
        local het_label "interactions of the regressors with: `het_exp'" ;
        /* names of the step-1 coefficients */
        local g_names `indepvars_exp' `het_exp' ;
        foreach xv of local indepvars_exp {;
            foreach zv of local het_exp {;
                local g_names `g_names' `xv'#`zv' ;
            } ;
        } ;
        local g_names `g_names' _cons ;
        capture matrix colnames _gepwreg_g = `g_names' ;
        matrix colnames _gepwreg_zbar = `het_exp' ;
        matrix colnames _gepwreg_zpop = `het_exp' ;
        matrix rownames _gepwreg_decomp = `het_exp' ;
        local dnames "mean_at_tau mean_pop" ;
        foreach xv of local indepvars_exp {; local dnames `dnames' contrib:`xv' ; } ;
        matrix colnames _gepwreg_decomp = `dnames' ;
    } ;

    /* ── T2. Heterogeneity model = conditional quantile regression ─────── */
    else {;
        if "`qgrid'" == "" {;
            local qgrid ;
            forvalues g = 5(5)95 {; local qgrid `qgrid' `=`g'/100' ; } ;
        } ;
        local n_grid : word count `qgrid' ;
        if `n_grid' < 3 {;
            di as error "qgrid() needs at least three quantiles" ;
            exit 198 ;
        } ;
        mata: st_matrix("_gepwreg_ugrid", strtoreal(tokens("`qgrid'"))) ;
        _gepwreg_qrG `depvar' `indepvars_mata' if `touse' [pw=`fw_var'], grid(`qgrid') ;
        mata: _gepwreg_ts_qr("`depvar'", "`indepvars_mata'", "`fw_var'", "`touse'",
                             `percentile', `cband', `band', `do_initial',
                             "_gepwreg_xref", "_gepwreg_G", "_gepwreg_ugrid", "_gepwreg_iflag",
                             `do_optbw') ;
        local het_label "conditional quantile regression, `n_grid' quantiles" ;
        /* (1.5.0) ref(): the estimator E -- beta(u*) with xref' beta(u*) = q_tau
           on the interpolated grid path -- and its delta-method influence
           function: the quantile regressions at the two grid quantiles around
           u* and the quantile q_tau, each with qreg's Powell-kernel density
           (vce(robust, kernel)), and the sampling of the reference means.     */
        if `do_ref' {;
            _gepwreg_refrule, spec(`"`ref'"') touse(`touse')
                xexp(`indepvars_exp') xmata(`indepvars_mata') xraw(`indepvars_exp_raw') ;
            local ref_nu = r(nu) ;
            local ref_unames `r(unames)' ;
            local ref_kinds  `r(kinds)' ;
            forvalues u = 1/`ref_nu' {; local ref_lev`u' `r(lev`u')' ; } ;
            tempvar refrq refkq refr1 refk1 refr2 refk2 ;
            /* q_tau: the tau-quantile of the weighted distribution of y, the
               same rule in the draws and along the profile, and qreg on a
               constant gives the bandwidth of its density                    */
            mata: st_numscalar("_gepwreg_qref", _gepwreg_qtau(st_data(., "`depvar'", "`touse'"), st_data(., "`fw_var'", "`touse'"), `percentile')) ;
            local qref = _gepwreg_qref ;
            qui qreg `depvar' if `touse' [pw=`fw_var'], quantile(`percentile') vce(robust, kernel) ;
            local kbq  = e(kbwidth) ;
            qui gen double `refrq' = `depvar' - `qref' if `touse' ;
            qui gen double `refkq' = 0.75*(1 - 0.2*(`refrq'/`kbq')^2)/sqrt(5) * (abs(`refrq'/`kbq') < sqrt(5)) if `touse' ;
            mata: _gepwreg_qr_ref1("`indepvars_mata'", "`fw_var'", "`touse'", `qref', "_gepwreg_G", "_gepwreg_ugrid") ;
            local ustar = _gepwreg_ustar ;
            local rseg  = _gepwreg_seg ;
            local rclip = _gepwreg_clip ;
            local ru1 : word `rseg' of `qgrid' ;
            local ru2 : word `=`rseg' + 1' of `qgrid' ;
            local refif = 1 ;
            foreach s in 1 2 {;
                capture qui qreg `depvar' `indepvars_mata' if `touse' [pw=`fw_var'], quantile(`ru`s'') vce(robust, kernel) ;
                if _rc {;
                    local refif = 0 ;
                    continue, break ;
                } ;
                local kb`s' = e(kbwidth) ;
                qui predict double `refr`s'' if `touse', residuals ;
                qui gen double `refk`s'' = 0.75*(1 - 0.2*(`refr`s''/`kb`s'')^2)/sqrt(5) * (abs(`refr`s''/`kb`s'') < sqrt(5)) if `touse' ;
            } ;
            if `refif' {;
                mata: _gepwreg_qr_ref2("`indepvars_mata'", "`fw_var'", "`touse'", "_gepwreg_G", "_gepwreg_ugrid", `percentile',
                                       "`refr1'", "`refk1'", `kb1', "`refr2'", "`refk2'", `kb2', "`refrq'", "`refkq'", `kbq',
                                       `do_svy', "`svy_psu'", "`svy_strata'") ;
                local seu = _gepwreg_seu ;
            } ;
        } ;
        /* (1.5.0) the profile: theta and E at every percentile of the grid */
        if `do_prof' {;
            mata: _gepwreg_qr_prof_pt("`depvar'", "`indepvars_mata'", "`fw_var'", "`touse'", "_gepwreg_pgrid",
                                      "_gepwreg_G", "_gepwreg_ugrid", `cband', `band', `do_optbw', `do_ref') ;
        } ;
        /* (1.5.0) the profile's bootstrap: ONE loop of draws for every
           percentile -- in each draw the grid of quantile regressions, the
           ranks and the reference household once, the bandwidth and q_tau at
           every percentile -- and the table at per() from the same draws    */
        if `boot' > 0 & `do_prof' {;
            local ncol = colsof(_gepwreg_pgrid) * `kx' * (1 + `do_ref') ;
            tempname TP ;
            matrix `TP' = J(`boot', `ncol', .) ;
            local rngsave "`c(rngstate)'" ;
            set seed `seed' ;
            local nfail = 0 ;
            if `do_dots' _dots 0, title(Bootstrap replications, every percentile) reps(`boot') ;
            forvalues bb = 1/`boot' {;
                preserve ;
                qui keep if `touse' ;
                if `do_svy' & "`svy_psu'" != "" {;
                    if "`svy_strata'" != "" qui bsample, strata(`svy_strata') cluster(`svy_psu') ;
                    else                    qui bsample, cluster(`svy_psu') ;
                } ;
                else qui bsample ;
                capture {;
                    _gepwreg_qrG `depvar' `indepvars_mata' if `touse' [pw=`fw_var'], grid(`qgrid') ;
                    mata: st_matrix("_gepwreg_Prow", _gepwreg_qr_prof("`depvar'", "`indepvars_mata'", "`fw_var'", "`touse'",
                                    "_gepwreg_pgrid", "_gepwreg_G", "_gepwreg_ugrid", `cband', `band', `do_optbw', `do_ref')) ;
                } ;
                local brc = _rc ;
                if `brc' == 0 matrix `TP'[`bb', 1] = _gepwreg_Prow ;
                else {;
                    if `nfail' == 0 local rcfirst = `brc' ;
                    local ++nfail ;
                } ;
                if `do_dots' _dots `bb' `brc' ;
                restore ;
            } ;
            set rngstate `rngsave' ;
            if `nfail' > 0 di as text "Note: `nfail' bootstrap replication(s) failed and were dropped (first error r(`rcfirst'))." ;
            /* the standard errors, the simultaneous bands and the flatness
               tests of the profile, and the variances at per()            */
            mata: _gepwreg_qr_prof_boot("`TP'", `kx', `do_ref', `level', `nuse') ;
            /* the point estimates of the full sample back in place */
            _gepwreg_qrG `depvar' `indepvars_mata' if `touse' [pw=`fw_var'], grid(`qgrid') ;
            mata: _gepwreg_ts_qr("`depvar'", "`indepvars_mata'", "`fw_var'", "`touse'",
                                 `percentile', `cband', `band', `do_initial',
                                 "_gepwreg_xref", "_gepwreg_G", "_gepwreg_ugrid", "_gepwreg_iflag",
                                 `do_optbw') ;
            matrix _gepwreg_V = _gepwreg_Vb ;
        } ;
        /* pairs (or stratified cluster) bootstrap of the whole procedure */
        else if `boot' > 0 {;
            tempname TT ;
            matrix `TT' = J(`boot', `kx', .) ;
            if `do_ref' {;
                tempname TTr ;
                matrix `TTr' = J(`boot', `kx', .) ;
            } ;
            tempname bpoint ;
            matrix `bpoint' = _gepwreg_b ;
            local rngsave "`c(rngstate)'" ;
            set seed `seed' ;
            local nfail = 0 ;
            /* (1.4) this loop runs in Stata, not in Mata, so it had none of
               the dots added to the other three.  It is also the only one
               that can FAIL a replication, and _dots prints an x rather than
               a dot when it is handed a non-zero return code.            */
            if `do_dots' _dots 0, title(Bootstrap replications) reps(`boot') ;
            forvalues bb = 1/`boot' {;
                preserve ;
                qui keep if `touse' ;
                if `do_svy' & "`svy_psu'" != "" {;
                    if "`svy_strata'" != "" qui bsample, strata(`svy_strata') cluster(`svy_psu') ;
                    else                    qui bsample, cluster(`svy_psu') ;
                } ;
                else qui bsample ;
                capture {;
                    _gepwreg_qrG `depvar' `indepvars_mata' if `touse' [pw=`fw_var'], grid(`qgrid') ;
                    mata: _gepwreg_ts_qr("`depvar'", "`indepvars_mata'", "`fw_var'", "`touse'",
                                         `percentile', `cband', `band', `do_initial',
                                         "_gepwreg_xref", "_gepwreg_G", "_gepwreg_ugrid", "_gepwreg_iflag",
                                         `do_optbw') ;
                    /* (1.5.0) ref(): the whole estimator of E in the draw --
                       q_tau, the reference household of the draw, u*, beta(u*) */
                    if `do_ref' {;
                        mata: st_matrix("_gepwreg_Erow", _gepwreg_qr_refb("`depvar'", "`indepvars_mata'", "`fw_var'", "`touse'", `percentile', "_gepwreg_G", "_gepwreg_ugrid")) ;
                    } ;
                } ;
                local brc = _rc ;
                if `brc' == 0 {;
                    forvalues j = 1/`kx' {; matrix `TT'[`bb', `j'] = _gepwreg_b[1, `j'] ; } ;
                    if `do_ref' {;
                        forvalues j = 1/`kx' {; matrix `TTr'[`bb', `j'] = _gepwreg_Erow[1, `j'] ; } ;
                    } ;
                } ;
                else {;
                    if `nfail' == 0 local rcfirst = `brc' ;
                    local ++nfail ;
                } ;
                if `do_dots' _dots `bb' `brc' ;
                restore ;
            } ;
            set rngstate `rngsave' ;
            if `nfail' > 0 di as text "Note: `nfail' bootstrap replication(s) failed and were dropped (first error r(`rcfirst'))." ;
            mata: st_matrix("_gepwreg_Vb", _gepwreg_bootvar(st_matrix("`TT'"))) ;
            if `do_ref' mata: st_matrix("_gepwreg_Vref_boot", _gepwreg_bootvar(st_matrix("`TTr'"))) ;
            /* restore the point estimates of the full sample */
            _gepwreg_qrG `depvar' `indepvars_mata' if `touse' [pw=`fw_var'], grid(`qgrid') ;
            mata: _gepwreg_ts_qr("`depvar'", "`indepvars_mata'", "`fw_var'", "`touse'",
                                 `percentile', `cband', `band', `do_initial',
                                 "_gepwreg_xref", "_gepwreg_G", "_gepwreg_ugrid", "_gepwreg_iflag",
                                 `do_optbw') ;
            matrix _gepwreg_V = _gepwreg_Vb ;
        } ;
    } ;

    /* ── T3. Retrieve scalars ──────────────────────────────────────────── */
    local h     = _gepwreg_h ;
    local n_eff = _gepwreg_neff ;
    local n_obs = _gepwreg_n ;
    local tau   = `percentile' ;
    if `do_initial' {;
        local h0    = _gepwreg_h0 ;
        local neff0 = _gepwreg_neff0 ;
    } ;
    local hasV = 1 ;
    capture confirm matrix _gepwreg_V ;
    if _rc local hasV = 0 ;
    /* (1.5.0) ref(): the reference household, as text */
    if `do_ref' {;
        local reftext ;
        local sep ;
        forvalues u = 1/`ref_nu' {;
            local v  : word `u' of `ref_unames' ;
            local kd : word `u' of `ref_kinds' ;
            local val = _gepwreg_refu[1, `u'] ;
            if `kd' == 4 local txt "`v' at the shares of its levels" ;
            else if "`ref_lev`u''" != "" {;
                local li = `val' + 1 ;
                local lev : word `li' of `ref_lev`u'' ;
                local lab "`lev'" ;
                local vl : value label `v' ;
                if "`vl'" != "" & "`lev'" != "." local lab : label `vl' `lev' ;
                local txt "`v' = `lab'" ;
            } ;
            else local txt = "`v' = " + strtrim(string(`val', "%10.4g")) ;
            if      `kd' == 0 local how "given" ;
            else if `kd' == 1 local how "mean" ;
            else if `kd' == 2 local how "most frequent value" ;
            else if `kd' == 3 local how "most frequent level" ;
            else              local how "mean" ;
            local reftext `"`reftext'`sep'`txt' (`how')"' ;
            local sep = char(59) + " " ;
        } ;
    } ;

    /* ── T4. Names, base levels, post ──────────────────────────────────── */
    local boot_label "Pairs bootstrap" ;
    if `do_svy' & "`svy_psu'" != "" {;
        if "`svy_strata'" != "" local boot_label "Cluster bootstrap (PSUs within strata)" ;
        else                       local boot_label "Cluster bootstrap (PSUs)" ;
    } ;
    local vnames `indepvars_exp' ;
    local vfull  `indepvars_exp_raw' ;
    local kfull : word count `vfull' ;
    local kest  : word count `vnames' ;
    if (`kfull' > `kest') {;
        tempname E ;
        matrix `E' = J(`kest', `kfull', 0) ;
        local i = 0 ;
        local j = 0 ;
        foreach v of local vfull {;
            local ++j ;
            if !regexm("`v'","[0-9]+b\.") & !regexm("`v'","[0-9]+o\.") {;
                local ++i ;
                matrix `E'[`i', `j'] = 1 ;
            } ;
        } ;
        matrix _gepwreg_b = _gepwreg_b * `E' ;
        capture matrix _gepwreg_batq = _gepwreg_batq * `E' ;
        capture matrix _gepwreg_iflag = _gepwreg_iflag * `E' ;
        if `hasV' matrix _gepwreg_V  = `E'' * _gepwreg_V  * `E' ;
        capture matrix _gepwreg_Vi = `E'' * _gepwreg_Vi * `E' ;
        capture matrix _gepwreg_Vs = `E'' * _gepwreg_Vs * `E' ;
        capture matrix _gepwreg_Vb = `E'' * _gepwreg_Vb * `E' ;
        capture matrix _gepwreg_bref = _gepwreg_bref * `E' ;
        capture matrix _gepwreg_refx = _gepwreg_refx * `E' ;
        foreach nm in IF svy boot {;
            capture matrix _gepwreg_Vref_`nm' = `E'' * _gepwreg_Vref_`nm' * `E' ;
        } ;
        local vnames `vfull' ;
    } ;
    matrix colnames _gepwreg_b = `vnames' ;
    timer off 99 ;
    qui timer list 99 ;
    local etime = r(t99) ;
    if `hasV' {;
        matrix rownames _gepwreg_V = `vnames' ;
        matrix colnames _gepwreg_V = `vnames' ;
        ereturn post _gepwreg_b _gepwreg_V, esample(`touse') obs(`n_obs')
                depname(`depvar') buildfvinfo findomitted ;
    } ;
    else {;
        ereturn post _gepwreg_b, esample(`touse') obs(`n_obs')
                depname(`depvar') buildfvinfo findomitted ;
    } ;
    capture confirm matrix _gepwreg_batq ;
    if !_rc {;
        matrix colnames _gepwreg_batq = `vnames' ;
        ereturn matrix b_atq = _gepwreg_batq ;
    } ;
    capture confirm matrix _gepwreg_Vi ;
    if !_rc {;
        matrix rownames _gepwreg_Vi = `vnames' ;
        matrix colnames _gepwreg_Vi = `vnames' ;
        ereturn matrix V_IF = _gepwreg_Vi ;
    } ;
    capture confirm matrix _gepwreg_Vs ;
    if !_rc {;
        matrix rownames _gepwreg_Vs = `vnames' ;
        matrix colnames _gepwreg_Vs = `vnames' ;
        ereturn matrix V_svy = _gepwreg_Vs ;
    } ;
    capture confirm matrix _gepwreg_Vb ;
    if !_rc {;
        matrix rownames _gepwreg_Vb = `vnames' ;
        matrix colnames _gepwreg_Vb = `vnames' ;
        ereturn matrix V_boot = _gepwreg_Vb ;
    } ;
    if "`het_type'" == "z" {;
        ereturn matrix b_step1  = _gepwreg_g ;
        ereturn matrix zbar_tau = _gepwreg_zbar ;
        ereturn matrix zbar_pop = _gepwreg_zpop ;
        ereturn matrix decomp   = _gepwreg_decomp ;
    } ;
    else {;
        ereturn matrix b_qgrid  = _gepwreg_G ;
        ereturn local  qgrid "`qgrid'" ;
    } ;
    ereturn matrix xref = _gepwreg_xref ;
    /* (1.5.0) the profile across the percentiles (graph()) */
    if `do_prof' {;
        foreach nm in b se cu flat rb rse rcu rflat {;
            capture confirm matrix _gepwreg_p`nm' ;
            if !_rc {;
                matrix colnames _gepwreg_p`nm' = `indepvars_exp' ;
                ereturn matrix prof_`nm' = _gepwreg_p`nm' ;
            } ;
        } ;
        ereturn matrix prof_grid = _gepwreg_pgrid ;
        /* het(qr) with ref(): the reference household's conditional rank u*
           at every percentile, and where it is held at an end of the grid  */
        capture confirm matrix _gepwreg_pru ;
        if !_rc {;
            matrix colnames _gepwreg_pru = ustar clip ;
            ereturn matrix prof_ru = _gepwreg_pru ;
        } ;
    } ;
    /* (1.5.0) ref(): the effect for the reference household */
    if `do_ref' {;
        if "`het_type'" == "z" {;
            local vref = cond(`do_svy', "_gepwreg_Vref_svy", "_gepwreg_Vref_IF") ;
            local ref_se_type = cond(`do_svy', "Taylor linearisation (survey design)", "Analytical (Influence Function)") ;
        } ;
        else if `boot' > 0 {;
            local vref "_gepwreg_Vref_boot" ;
            local ref_se_type "`boot_label'" ;
            if "`ref_se_type'" == "" local ref_se_type "Bootstrap" ;
        } ;
        else {;
            local vref = cond(`do_svy', "_gepwreg_Vref_svy", "_gepwreg_Vref_IF") ;
            local ref_se_type = cond(`do_svy', "Delta method, Taylor linearisation (approximate: see the warning)", "Delta method (approximate: see the warning)") ;
        } ;
        matrix colnames _gepwreg_bref = `vnames' ;
        matrix colnames _gepwreg_refx = `vnames' ;
        capture confirm matrix `vref' ;
        if !_rc {;
            tempname VR ;
            matrix `VR' = `vref' ;
            matrix rownames `VR' = `vnames' ;
            matrix colnames `VR' = `vnames' ;
            ereturn matrix V_ref = `VR' ;
        } ;
        else local ref_se_type "none" ;
        foreach nm in IF svy boot {;
            capture confirm matrix _gepwreg_Vref_`nm' ;
            if !_rc {;
                matrix rownames _gepwreg_Vref_`nm' = `vnames' ;
                matrix colnames _gepwreg_Vref_`nm' = `vnames' ;
                ereturn matrix V_ref_`nm' = _gepwreg_Vref_`nm' ;
            } ;
        } ;
        ereturn matrix b_ref = _gepwreg_bref ;
        ereturn matrix ref_x = _gepwreg_refx ;
        ereturn local  ref         `"`ref'"' ;
        ereturn local  ref_text    `"`reftext'"' ;
        ereturn local  ref_SE_type "`ref_se_type'" ;
        if "`het_type'" == "qr" {;
            ereturn scalar u_ref = `ustar' ;
            if "`seu'" != "" ereturn scalar u_ref_se = `seu' ;
        } ;
    } ;
    if `do_initial' {;
        matrix colnames _gepwreg_iflag = `vnames' ;
        ereturn matrix y0_regressors = _gepwreg_iflag ;
        ereturn local  initial "`initial'" ;
    } ;
    ereturn scalar tau    = `tau' ;
    ereturn scalar h      = `h' ;
    ereturn scalar N_eff  = `n_eff' ;
    if `do_initial' {;
        ereturn scalar h0     = `h0' ;
        ereturn scalar N_eff0 = `neff0' ;
    } ;
    ereturn scalar boot   = `boot' ;
    if (`do_merr') {;
        ereturn matrix merr_s2   = _gepwreg_s2 ;
        ereturn matrix merr_t    = _gepwreg_tm ;
        ereturn matrix merr_keep = _gepwreg_keep ;
        capture matrix colnames _gepwreg_mdiag = reliability skewness R2_on_X ;
        capture ereturn matrix merr_diag = _gepwreg_mdiag ;
        ereturn scalar merr      = 1 ;
        ereturn scalar tcrit     = `tcrit' ;
    }; 
    ereturn scalar do_svy = `do_svy' ;
    ereturn scalar etime  = `etime' ;
    ereturn local  bw_method "`bw_label'" ;
    if `do_initial' ereturn local method "initial" ;
    else            ereturn local method "twostep" ;
    if "`het_type'" == "z" ereturn local het "`het_exp'" ;
    else                   ereturn local het "qr" ;
    ereturn local  cmd      "gepwreg" ;
    ereturn local  rankvar  "`depvar'" ;
    ereturn local  ranklabel "`rank_label'" ;
    ereturn local  wgt      "`wgt_name'" ;
    ereturn local  wsrc     "`wsrc'" ;
    if "`wgt_type'" != "" {;
        ereturn local wtype "`wgt_type'" ;
        ereturn local wexp  "= `wgt_exp'" ;
    } ;
    ereturn local  cmdline  "gepwreg `0'" ;
    ereturn local  title    "Percentile Weights Regression (two-step)" ;
    if `do_initial' {;
        ereturn local SE_type "`boot_label'" ;
        ereturn local vce     "bootstrap" ;
    } ;
    else if "`het_type'" == "z" {;
        if `do_svy' {;
            ereturn local SE_type "Taylor linearisation (survey design)" ;
            ereturn local vce     "linearized" ;
        } ;
        else {;
            ereturn local SE_type "Analytical (Influence Function)" ;
            ereturn local vce     "IF" ;
        } ;
    } ;
    else {;
        if `boot' > 0 {;
            ereturn local SE_type "`boot_label'" ;
            ereturn local vce     "bootstrap" ;
        } ;
        else {;
            ereturn local SE_type "none (specify boot())" ;
            ereturn local vce     "none" ;
        } ;
    } ;

    /* ── T5. Display ───────────────────────────────────────────────────── */
    di "" ;
    di as text "Percentile Weights Regression: two-step estimator" ;
    di as text "{hline 72}" ;
    if `do_initial' {;
        di as text %28s "Object" " = " as result "effect of x for units at the tau-quantile of y0" ;
        di as text %28s "" "   " as result "y0 = y - b_i (x - xref), the outcome before `initial'" ;
    } ;
    else {;
        di as text %28s "Object" " = " as result "effect of x for units at the tau-quantile of `depvar'" ;
    } ;
    di as text %28s "Target percentile (tau)" " = " as result %6.4f `tau' ;
    di as text %28s "Heterogeneity model" " = " as result "`het_label'" ;
    di as text %28s "Bandwidth (h)"           " = " as result %9.6f `h' as text " (`bw_label', rank of `depvar')" ;
    if `do_optbw' & `band' == 0 di as text _col(33) "(chosen at this tau; silverman or band() to fix it)" ;
    if `do_initial' di as text %28s "Bandwidth on rank of y0" " = " as result %9.6f `h0' ;
    di as text %28s "Weights"                 " = " as result "`wgt_name'" ;
    di as text %28s "Observations"            " = " as result %7.0f `n_obs' ;
    if `do_initial' di as text %28s "Kernel N_eff (rank of y0)" " = " as result %7.1f `neff0' ;
    else            di as text %28s "Kernel N_eff"             " = " as result %7.1f `n_eff' ;
    di as text %28s "Execution time"          " = " as result %9.3f `etime' as text " seconds" ;
    if `boot' > 0 di as text %28s "Bootstrap replications" " = " as result %7.0f `boot' ;
    if `do_merr' {;
        di as text %28s "Measurement-error corr." " = " as text
            "sigma2 per heterogeneity variable, t > " as result %3.1f `tcrit' ;
        tempname mS mT mK mD ;
        matrix `mS' = e(merr_s2) ;
        matrix `mT' = e(merr_t) ;
        matrix `mK' = e(merr_keep) ;
        capture matrix `mD' = e(merr_diag) ;
        local hasD = (_rc == 0) ;
        /* (1.4) reliability, skewness and R2(z|X) are shown beside sigma2
           and t.  "t too small" alone covered two situations that call for
           different conduct: a variable that is identifiable but not resolved
           by this sample, where more data would help, and one that is
           symmetric and uncorrelated with the regressors, where no sample size
           ever will.  The three numbers separate them without a threshold. */
        if `hasD' {;
            di as text _col(21) %10s "sigma2" %7s "t" %8s "reliab."
                                %7s "skew" %8s "R2(z|X)" %8s "applied" ;
        } ;
        local nlow  = 0 ;
        local nkeep = 0 ;
        local jm = 0 ;
        foreach zv of local het_exp {;
            local ++jm ;
            local s2v = `mS'[1, `jm'] ;
            local tv  = `mT'[1, `jm'] ;
            local kv  = `mK'[1, `jm'] ;
            local vrd = cond(`kv', "corrected", "not corrected (t too small)") ;
            if `kv' local nkeep = `nkeep' + 1 ;
            if `kv' & `tv' < 1.5 * `tcrit' local nlow = 1 ;
            if `hasD' {;
                di as text %20s abbrev("`zv'", 20) as result %10.6f `s2v'
                   %7.2f `tv' %8.2f `mD'[`jm',1] %7.2f `mD'[`jm',2]
                   %8.2f `mD'[`jm',3] as text %8s cond(`kv', "yes", "no") ;
            } ;
            else {;
                di as text %28s "`zv'" " = " as result %9.6f `s2v'
                   as text "   t = " as result %6.2f `tv' as text "   `vrd'" ;
            } ;
        }; 
        /* Two notes, both conditional.  A note that prints every time stops
           being read; it is the condition that keeps its force.            */
        if `nlow' {;
            di as text _col(3)
               "A t sits just above tcrit(): there the estimator is a pre-test" ;
            di as text _col(3)
               "estimator and no standard error is valid uniformly.  Report the" ;
            di as text _col(3)
               "corrected and the uncorrected profile side by side." ;
        } ;
        if `nkeep' > 0 {;
            di as text _col(3)
               "The third moments cannot separate a missing functional form from" ;
            di as text _col(3)
               "measurement error.  Fit the enriched model WITHOUT merr: if it" ;
            di as text _col(3)
               "lands where merr did, merr was doing the work of the missing" ;
            di as text _col(3)
               "term.  See {help gepwreg##options:help gepwreg}." ;
        } ;
    }; 
    di as text "{hline 72}" ;
    di as text "SE: `e(SE_type)'" ;
    di as text "{hline 72}" ;
    di "" ;
    if `hasV' ereturn display, level(`level') ;
    else {;
        di as text %14s "" "{c |}" %12s "Coefficient" ;
        di as text "{hline 14}{c +}{hline 14}" ;
        local j = 0 ;
        foreach nm of local vnames {;
            local ++j ;
            if regexm("`nm'","[0-9]+b\.") | regexm("`nm'","[0-9]+o\.") continue ;
            di as text %14s abbrev("`nm'",14) "{c |}" as result %12.5f _b[`nm'] ;
        } ;
        di as text "{hline 14}{c BT}{hline 14}" ;
        di as text "Standard errors: specify boot(#) (bootstrap of the two steps)." ;
    } ;
    if "`het_type'" == "z" {;
        di "" ;
        if `do_initial' di as text "Composition of the group at the tau-quantile of y0 and its contribution to the effect" ;
        else            di as text "Composition of the group at tau and its contribution to the effect" ;
        di as text "(effect at tau = coefficient of x in step 1 + sum of the contributions)" ;
        matlist e(decomp), format(%10.4f) ;
    } ;
    if `do_initial' {;
        di "" ;
        di as text "Effect at the tau-quantile of `depvar' itself stored in e(b_atq)." ;
    } ;
    /* (1.5.0) ref(): the effect for the reference household */
    if `do_ref' {;
        di "" ;
        di as text "Effect for the reference household" ;
        di as text "{p 2 4 2}`reftext'{p_end}" ;
        if "`het_type'" == "qr" {;
            di as text "  conditional rank of this household at q_tau: u* = " as result %6.4f e(u_ref)
               as text "  (grid quantiles " as result %4.2f `ru1' as text " and " as result %4.2f `ru2' as text ")" ;
        } ;
        di as text "  SE: `e(ref_SE_type)'" ;
        capture confirm matrix e(V_ref) ;
        if !_rc _gepwreg_reftab e(b_ref) e(V_ref) `level' ;
        else    _gepwreg_reftab e(b_ref) "" `level' ;
        if "`het_type'" == "z" {;
            di as text "{p 2 2 2}Under het(z) the effect of a given household does not depend on tau: the"
               " profile of e(b) across tau is the composition of the group at tau.{p_end}" ;
        } ;
        else {;
            if `rclip' {;
                di as text "{p 2 2 2}Warning: q_tau lies beyond the conditional quantiles of the reference"
                   " household on the grid, and its effect is taken at the end of the grid.{p_end}" ;
            } ;
            if `boot' == 0 {;
                di as text "{p 2 2 2}Warning: these standard errors are approximate.  They rest on the"
                   " density of the conditional quantiles at u* (Powell kernel), which is imprecise in"
                   " finite samples -- off by tens of percent in the tests of the package.  Use boot(#)"
                   " for inference.{p_end}" ;
                if "`seu'" != "" {;
                    local dgrid = min(abs(`ustar' - `ru1'), abs(`ru2' - `ustar')) ;
                    if `dgrid' < 2*`seu' {;
                        di as text "{p 2 2 2}Warning: u* is within two standard errors ("
                           %6.4f `seu' ") of a grid quantile, where the estimate has a kink."
                           "  The delta method can be off there -- use boot(#).{p_end}" ;
                    } ;
                } ;
            } ;
        } ;
    } ;
    di "" ;
    di as text "The one-step weighted regression on outcome-ranked units (version 1.3)" ;
    di as text "is available with rankdep; it does not estimate the effect (help gepwreg)." ;
    if "`setable'" != "" gepwreg_setable ;
    exit ;
} ;

/* ═════════════════════════════════════════════════════════════════════════
   R. ONE-STEP WEIGHTED REGRESSION: rankvar(z) (profile along z) and rankdep
   ═════════════════════════════════════════════════════════════════════════ */
if "`mode'" == "rankdep" {;
    di as text "Note: rankdep -- weighted regression on outcome-ranked units. Its slope is" ;
    di as text "      a descriptive quantity, not the effect of x (see help gepwreg)." ;
} ;

/* ── 4. Call Mata ────────────────────────────────────────────────────────── */
scalar _gepwreg_h_tmp = 0 ;
mata: _gepwreg_main("`depvar'", "`indepvars_mata'", "`fw_var'", "`touse'",
                  `percentile', `cband', `band', `boot', `addcons',
                  `do_svy', "`svy_psu'", "`svy_strata'",
                  "`rank_var'", `do_optbw') ;

/* ── 5. Retrieve scalars ─────────────────────────────────────────────────── */
local h     = _gepwreg_h ;
local n_eff = _gepwreg_neff ;
local n_obs = _gepwreg_n ;
local tau   = `percentile' ;

/* ── 6. Name matrices ────────────────────────────────────────────────────── */
if `addcons' {;
    local vnames `indepvars_exp' _cons ;
    local vfull  `indepvars_exp_raw' _cons ;
} ;
else {;
    local vnames `indepvars_exp' ;
    local vfull  `indepvars_exp_raw' ;
} ;
local kfull : word count `vfull' ;
local kest  : word count `vnames' ;
if (`kfull' > `kest') {;
    tempname E ;
    matrix `E' = J(`kest', `kfull', 0) ;
    local i = 0 ;
    local j = 0 ;
    foreach v of local vfull {;
        local ++j ;
        if !regexm("`v'","[0-9]+b\.") & !regexm("`v'","[0-9]+o\.") {;
            local ++i ;
            matrix `E'[`i', `j'] = 1 ;
        } ;
    } ;
    matrix _gepwreg_b  = _gepwreg_b  * `E' ;
    matrix _gepwreg_V  = `E'' * _gepwreg_V  * `E' ;
    matrix _gepwreg_Vn = `E'' * _gepwreg_Vn * `E' ;
    matrix _gepwreg_Vi = `E'' * _gepwreg_Vi * `E' ;
    if `do_svy'   matrix _gepwreg_Vs = `E'' * _gepwreg_Vs * `E' ;
    if `boot' > 0 matrix _gepwreg_Vb = `E'' * _gepwreg_Vb * `E' ;
    local vnames `vfull' ;
} ;

matrix colnames _gepwreg_b  = `vnames' ;
matrix rownames _gepwreg_V  = `vnames' ;
matrix colnames _gepwreg_V  = `vnames' ;
matrix rownames _gepwreg_Vn = `vnames' ;
matrix colnames _gepwreg_Vn = `vnames' ;
matrix rownames _gepwreg_Vi = `vnames' ;
matrix colnames _gepwreg_Vi = `vnames' ;

/* ── 7. Post results ─────────────────────────────────────────────────────── */
timer off 99 ;
qui timer list 99 ;
local etime = r(t99) ;
ereturn post _gepwreg_b _gepwreg_V, esample(`touse') obs(`n_obs')
        depname(`depvar') buildfvinfo findomitted ;
ereturn matrix V_naive = _gepwreg_Vn ;
ereturn matrix V_IF    = _gepwreg_Vi ;
if `do_svy' {;
    matrix rownames _gepwreg_Vs = `vnames' ;
    matrix colnames _gepwreg_Vs = `vnames' ;
    ereturn matrix V_svy = _gepwreg_Vs ;
} ;
if `boot' > 0 {;
    matrix rownames _gepwreg_Vb = `vnames' ;
    matrix colnames _gepwreg_Vb = `vnames' ;
    ereturn matrix V_boot = _gepwreg_Vb ;
} ;
ereturn scalar tau   = `tau' ;
ereturn scalar h     = `h' ;
ereturn local  bw_method "`bw_label'" ;
ereturn scalar N_eff = `n_eff' ;
ereturn scalar boot   = `boot' ;
ereturn scalar do_svy = `do_svy' ;
ereturn scalar etime  = `etime' ;
ereturn local  method   "`mode'" ;
ereturn local  cmd      "gepwreg" ;
ereturn local  rankvar  "`rank_var'" ;
ereturn local  ranklabel "`rank_label'" ;
ereturn local  wgt      "`wgt_name'" ;
ereturn local  wsrc     "`wsrc'" ;
if "`wgt_type'" != "" {;
    ereturn local wtype "`wgt_type'" ;
    ereturn local wexp  "= `wgt_exp'" ;
} ;
ereturn local  cmdline  "gepwreg `0'" ;
ereturn local  title    "Percentile Weights Regression" ;
if `do_svy' {;
    ereturn local  SE_type  "Taylor linearisation (survey design)" ;
    ereturn local  vce      "linearized" ;
} ;
else {;
    ereturn local  SE_type  "Analytical (Influence Function)" ;
    ereturn local  vce      "IF-corrected" ;
} ;

/* ── 8. Display ──────────────────────────────────────────────────────────── */
di "" ;
if "`mode'" == "rankvar" di as text "Percentile Weights Regression: profile along the rank of `rank_var'" ;
else                     di as text "Percentile Weights Regression: one-step weighted regression (rankdep)" ;
di as text "{hline 72}" ;
if "`mode'" == "rankvar" di as text %28s "Object" " = " as result "effect of x for units at the tau-quantile of `rank_var'" ;
else                     di as text %28s "Object" " = " as result "descriptive slope among units at the tau-quantile of `depvar'" ;
di as text %28s "Target percentile (tau)" " = " as result %6.4f `tau' ;
di as text %28s "Bandwidth (h)"           " = " as result %9.6f `h' as text " (`bw_label')" ;
di as text %28s "Weights"                 " = " as result "`wgt_name'" ;
di as text %28s "Ranking variable"        " = " as result "`rank_label'" ;
di as text %28s "Observations"            " = " as result %7.0f `n_obs' ;
di as text %28s "Kernel N_eff"            " = " as result %7.1f `n_eff' ;
if `boot' > 0 {;
    di as text %28s "Bootstrap replications" " = " as result %7.0f `boot' ;
} ;
di as text %28s "Execution time" " = " as result %9.3f `etime' as text " seconds" ;
di as text "{hline 72}" ;
if `do_svy' di as text "SE: Taylor linearisation under survey design (PSU + strata)" ;
else        di as text "SE: IF-corrected analytical (Deville 1999)" ;
di as text "{hline 72}" ;
di "" ;
ereturn display, level(`level') ;
di "" ;
di as text "Naive WLS SE stored in e(V_naive)." ;
if `boot' > 0 di as text "Bootstrap SE stored in e(V_boot)." ;
if "`setable'" != "" gepwreg_setable ;
else di as text "Use {cmd:gepwreg_setable} (or the setable option) to compare all SE types." ;

end ;


/* -----------------------------------------------------------------------------
   _gepwreg_qrG : step 1 of het(qr) -- conditional quantile regressions on a
   grid of quantiles; the coefficients (regressors, then the constant) go to
   the Stata matrix _gepwreg_G, one row per grid quantile.
   ------------------------------------------------------------------------- */
capture program drop _gepwreg_qrG ;
program define _gepwreg_qrG ;
    version 16.0 ;
    syntax varlist(min=1 numeric) [if] [pw], GRID(numlist) ;
    local depvar : word 1 of `varlist' ;
    local xvars  : list varlist - depvar ;
    local k : word count `xvars' ;
    local G : word count `grid' ;
    matrix _gepwreg_G = J(`G', `k' + 1, .) ;
    /* (1.5.0) qreg's own path to its coefficients -- the weighted least-squares
       start (_qregwls), qreg's stable sort on its residual, the solver _qreg,
       a pweight passed as an iweight -- without the variance that qreg then
       computes and gepwreg never uses (two more quantile regressions at
       u -/+ h for a pweight).  The coefficients are qreg's, to the bit, at a
       quarter of the time.                                                 */
    local wq ;
    if "`weight'" != "" local wq "[iweight`exp']" ;
    /* qreg gives the data their order back after each call, and its stable
       sort breaks the ties of the residual by that order: every quantile
       regression of the grid starts from the order on entry, here too, or a
       degenerate fit (ties, frequent in a bootstrap draw) could land on
       another of its solutions                                             */
    tempvar qr_r qr_i qr_o ;
    qui gen `c(obs_t)' `qr_o' = _n ;
    local g = 0 ;
    foreach u of local grid {;
        local ++g ;
        capture drop `qr_r' `qr_i' ;
        sort `qr_o' ;
        qui gen `c(obs_t)' `qr_i' = _n ;
        capture qui _qregwls `depvar' `xvars' `wq' `if', r(`qr_r') iterate(1) quant(`u') ;
        if !_rc {;
            sort `qr_r' `qr_i' ;
            capture qui _qreg `depvar' `xvars' `if' `wq', quant(`u') ;
            if !_rc & r(convcode) != 1 local rc0 = 430 ;
            else local rc0 = _rc ;
        } ;
        else local rc0 = _rc ;
        if `rc0' {;
            /* a second attempt with more iterations; then give up with the code */
            local rc1 = `rc0' ;
            capture qui _qreg `depvar' `xvars' `if' `wq', quant(`u') iterate(5000) ;
            if _rc | r(convcode) != 1 {;
                di as error "_gepwreg_qrG: the quantile regression failed at quantile `u' (r(`rc1'))" ;
                exit `rc1' ;
            } ;
        } ;
        local j = 0 ;
        foreach xv of local xvars {;
            local ++j ;
            matrix _gepwreg_G[`g', `j'] = _b[`xv'] ;
        } ;
        matrix _gepwreg_G[`g', `k' + 1] = _b[_cons] ;
    } ;
    sort `qr_o' ;
end ;


/* -----------------------------------------------------------------------------
   (1.5.0) _gepwreg_refrule : the reference household of ref(), as a RULE that
   Mata applies to any data, the sample or a bootstrap draw.  One row per
   underlying variable of the regressor terms (and of the het() terms):
     _gepwreg_rk    0 a value given in ref(), 1 the weighted mean (continuous,
                    or every variable with ref(mean)), 2 the most frequent
                    value (a 0/1 variable), 3 the most frequent level (a factor
                    variable, base level included), 4 a factor at the shares of
                    its levels, each indicator at its mean (ref(mean))
     _gepwreg_rfac  1 for a factor variable
     _gepwreg_rx    incidence on the regressor terms, _gepwreg_rz on the het()
                    terms (het(z) only)
     _gepwreg_xfix, _gepwreg_zfix  the values given in ref().
   r(nu), r(unames), r(kinds), r(lev#) the levels of factor # (base first).
   ------------------------------------------------------------------------- */
capture program drop _gepwreg_refrule ;
program define _gepwreg_refrule, rclass ;
    version 16.0 ;
    syntax , SPEC(string asis) TOUSE(varname) XEXP(string) XMATA(string) XRAW(string)
             [ ZEXP(string) ZMATA(string) ZRAW(string) ] ;
    local spec `spec' ;
    local k : word count `xexp' ;
    local m : word count `zexp' ;
    local unames ;
    foreach sd in x z {;
        local nn = cond("`sd'" == "x", `k', `m') ;
        forvalues c = 1/`nn' {;
            local t : word `c' of ``sd'exp' ;
            if strpos("`t'", "#") {;
                di as error "ref(): `t' is an interaction term -- ref() takes models with main effects only" ;
                exit 198 ;
            } ;
            local lev ;
            if regexm("`t'", "^([0-9]+)[a-z]*\.(.+)$") {;
                local lev = regexs(1) ;
                local v   = regexs(2) ;
            } ;
            else local v = regexr("`t'", "^c\.", "") ;
            local u : list posof "`v'" in unames ;
            if `u' == 0 {;
                local unames `unames' `v' ;
                local u : word count `unames' ;
            } ;
            local `sd'pos`u' ``sd'pos`u'' `c' ;
            local `sd'col`u' ``sd'col`u'' `: word `c' of ``sd'mata'' ;
            if "`lev'" != "" {;
                local isfac`u' 1 ;
                local `sd'lev`u' ``sd'lev`u'' `lev' ;
            } ;
            else local ispla`u' 1 ;
        } ;
    } ;
    local nu : word count `unames' ;
    /* the values given: default, or var=value pairs */
    local spec : subinstr local spec " =" "=", all ;
    local spec : subinstr local spec "= " "=", all ;
    local allmean 0 ;
    foreach tok of local spec {;
        if lower("`tok'") == "default" continue ;
        if lower("`tok'") == "mean" {;
            local allmean 1 ;
            continue ;
        } ;
        if !regexm("`tok'", "^([^=]+)=(.+)$") {;
            di as error "ref(): `tok' -- write ref(default), ref(mean), or var=value pairs, as in ref(size=5 urban=1)" ;
            exit 198 ;
        } ;
        local gv   = regexs(1) ;
        local gval = regexs(2) ;
        capture confirm number `gval' ;
        if _rc {;
            di as error "ref(): the value of `gv' must be a number, not `gval'" ;
            exit 198 ;
        } ;
        local u : list posof "`gv'" in unames ;
        if `u' == 0 {;
            di as error "ref(): `gv' is not a variable of the model (the variables are: `unames')" ;
            exit 198 ;
        } ;
        local given`u' `gval' ;
    } ;
    matrix _gepwreg_rk   = J(`nu', 1, .) ;
    matrix _gepwreg_rfac = J(`nu', 1, 0) ;
    matrix _gepwreg_rx   = J(`nu', `k', 0) ;
    matrix _gepwreg_xfix = J(1, `k', 0) ;
    if `m' > 0 {;
        matrix _gepwreg_rz   = J(`nu', `m', 0) ;
        matrix _gepwreg_zfix = J(1, `m', 0) ;
    } ;
    local kinds ;
    forvalues u = 1/`nu' {;
        local v : word `u' of `unames' ;
        if "`isfac`u''" == "1" & "`ispla`u''" == "1" {;
            di as error "ref(): `v' enters both as a factor variable and as a variable, so ref() cannot give it one value" ;
            exit 198 ;
        } ;
        foreach c of local xpos`u' {; matrix _gepwreg_rx[`u', `c'] = 1 ; } ;
        foreach c of local zpos`u' {; matrix _gepwreg_rz[`u', `c'] = 1 ; } ;
        if "`isfac`u''" == "1" {;
            matrix _gepwreg_rfac[`u', 1] = 1 ;
            local base ;
            foreach t in `xraw' `zraw' {;
                if regexm("`t'", "^([0-9]+)b\.`v'$") local base = regexs(1) ;
            } ;
            if "`base'" == "" local base "." ;
            local levs `xlev`u'' ;
            if "`levs'" == "" local levs `zlev`u'' ;
            return local lev`u' "`base' `levs'" ;
        } ;
        if "`given`u''" != "" {;
            matrix _gepwreg_rk[`u', 1] = 0 ;
            local kinds `kinds' 0 ;
            if "`isfac`u''" == "1" {;
                if "`given`u''" != "`base'" {;
                    local p : list posof "`given`u''" in levs ;
                    if `p' == 0 {;
                        di as error "ref(): `given`u'' is not a level of `v' (levels: `base' `levs')" ;
                        exit 198 ;
                    } ;
                    local cx : word `p' of `xpos`u'' ;
                    if "`cx'" != "" matrix _gepwreg_xfix[1, `cx'] = 1 ;
                    local cz : word `p' of `zpos`u'' ;
                    if "`cz'" != "" matrix _gepwreg_zfix[1, `cz'] = 1 ;
                } ;
            } ;
            else {;
                foreach c of local xpos`u' {; matrix _gepwreg_xfix[1, `c'] = `given`u'' ; } ;
                foreach c of local zpos`u' {; matrix _gepwreg_zfix[1, `c'] = `given`u'' ; } ;
            } ;
        } ;
        else if "`isfac`u''" == "1" {;
            local kd = cond(`allmean', 4, 3) ;
            matrix _gepwreg_rk[`u', 1] = `kd' ;
            local kinds `kinds' `kd' ;
        } ;
        else {;
            local col : word 1 of `xcol`u'' ;
            if "`col'" == "" local col : word 1 of `zcol`u'' ;
            qui count if `touse' & !inlist(`col', 0, 1) ;
            local kd = cond(r(N) == 0 & !`allmean', 2, 1) ;
            matrix _gepwreg_rk[`u', 1] = `kd' ;
            local kinds `kinds' `kd' ;
        } ;
    } ;
    return scalar nu = `nu' ;
    return local unames "`unames'" ;
    return local kinds  "`kinds'" ;
end ;


/* (1.5.0) _gepwreg_reftab : the table of the effect for the reference household */
capture program drop _gepwreg_reftab ;
program define _gepwreg_reftab ;
    version 16.0 ;
    args bname vname level ;
    tempname b V ;
    matrix `b' = `bname' ;
    local hasV = (`"`vname'"' != "") ;
    if `hasV' matrix `V' = `vname' ;
    local zc = invnormal((100 + `level')/200) ;
    local names : colnames `b' ;
    di as text "{hline 13}{c TT}{hline 65}" ;
    di as text %12s "" " {c |}" %12s "Coefficient" %11s "Std. err." %8s "z" %8s "P>|z|"
       %26s "[`level'% conf. interval]" ;
    di as text "{hline 13}{c +}{hline 65}" ;
    local j 0 ;
    foreach nm of local names {;
        local ++j ;
        if regexm("`nm'", "[0-9]+b\.") | regexm("`nm'", "[0-9]+o\.") continue ;
        local bj = `b'[1, `j'] ;
        local sj = . ;
        if `hasV' local sj = sqrt(`V'[`j', `j']) ;
        if `sj' < . & `sj' > 0 {;
            local zj = `bj'/`sj' ;
            di as text %12s abbrev("`nm'", 12) " {c |}" as result "  " %10.0g `bj' "  " %9.0g `sj'
               "  " %6.2f `zj' "  " %6.3f 2*normal(-abs(`zj')) "    " %10.0g `bj' - `zc'*`sj'
               "  " %10.0g `bj' + `zc'*`sj' ;
        } ;
        else di as text %12s abbrev("`nm'", 12) " {c |}" as result "  " %10.0g `bj' ;
    } ;
    di as text "{hline 13}{c BT}{hline 65}" ;
end ;


/* ─────────────────────────────────────────────────────────────────────────────
   MATA ENGINE
   #delimit cr required: Mata uses { } and ; internally — must not be
   intercepted by Stata's #delimit ; mode
   ───────────────────────────────────────────────────────────────────────── */
#delimit cr

/* Drop previous definitions so the ado can be sourced more than once */
capture mata: mata drop _gepwreg_wp()
capture mata: mata drop _gepwreg_hopt()
capture mata: mata drop _gepwreg_wls()
capture mata: mata drop _gepwreg_revCumSum()
capture mata: mata drop _gepwreg_revCumSumT()
capture mata: mata drop _gepwreg_se_IF()
capture mata: mata drop _gepwreg_se_naive()
capture mata: mata drop _gepwreg_se_boot()
capture mata: mata drop _gepwreg_se_svy()
capture mata: mata drop _gepwreg_main()
capture mata: mata drop _gepwreg_taylor_psi()
capture mata: mata drop _gepwreg_drawidx()
capture mata: mata drop _gepwreg_interp()
capture mata: mata drop _gepwreg_kw()
capture mata: mata drop _gepwreg_hopt2()
capture mata: mata drop _gepwreg_bootvar()
capture mata: mata drop _gepwreg_Wmat()
capture mata: mata drop _gepwreg_betai()
capture mata: mata drop _gepwreg_gmw_ix()
capture mata: mata drop _gepwreg_gmw_mean()
capture mata: mata drop _gepwreg_gmw_diag()
capture mata: mata drop _gepwreg_gmw_dots_on()
capture mata: mata drop _gepwreg_gmw_Omega()
capture mata: mata drop _gepwreg_gmw_q()
capture mata: mata drop _gepwreg_gmw_m1()
capture mata: mata drop _gepwreg_gmw_m2()
capture mata: mata drop _gepwreg_gmw_pre()
capture mata: mata drop _gepwreg_gmw_gfast()
capture mata: mata drop _gepwreg_gmw_critf()
capture mata: mata drop _gepwreg_gmw_hi()
capture mata: mata drop _gepwreg_gmw_vcov()
capture mata: mata drop _gepwreg_gmw_u()
capture mata: mata drop _gepwreg_gmw_inf()
capture mata: mata drop _gepwreg_gmw_search()
capture mata: mata drop _gepwreg_gmw_fit_t()
capture mata: mata drop _gepwreg_gmw_moms()
capture mata: mata drop _gepwreg_ts_z()
capture mata: mata drop _gepwreg_ts_qr()
capture mata: mata drop _gepwreg_refrow()
capture mata: mata drop _gepwreg_mrow()
capture mata: mata drop _gepwreg_ref_z()
capture mata: mata drop _gepwreg_E_pt()
capture mata: mata drop _gepwreg_qrpsi()
capture mata: mata drop _gepwreg_qr_ref1()
capture mata: mata drop _gepwreg_qr_ref2()
capture mata: mata drop _gepwreg_qr_refb()
capture mata: mata drop _gepwreg_qtau()
capture mata: mata drop _gepwreg_ifmean()
capture mata: mata drop _gepwreg_supt()
capture mata: mata drop _gepwreg_flat()
capture mata: mata drop _gepwreg_prof_stats()
capture mata: mata drop _gepwreg_prof_bstats()
capture mata: mata drop _gepwreg_qr_prof()
capture mata: mata drop _gepwreg_qr_prof_pt()
capture mata: mata drop _gepwreg_qr_prof_boot()

mata:

/* ─────────────────────────────────────────────────────────────────────────────
   _gepwreg_wp()
   Compute gepwe kernel weights and percentile ranks.
   Returns n x 2 matrix : col 1 = w_i,  col 2 = pc_i
   Bandwidth h is written to Stata scalar _gepwreg_h_tmp

   Kernel (exact gepwe.ado formula):
     w_i = exp(-0.25 * ((pc_i - tau)/h)^2) / (h * sqrt(2*pi) * n)

   Silverman bandwidth on the percentile scale:
     h = cband * min(sd(pc), IQR(pc)/1.34) * n^(-1/5)
   Override with band_in > 0
   ───────────────────────────────────────────────────────────────────────── */
real matrix _gepwreg_wp(real colvector y,
                      real colvector fw,
                      real scalar    tau,
                      real scalar    cband,
                      real scalar    band_in,
                      real colvector zrank)
{
    /* zrank : alternative ranking variable z (rows=n) or empty (rows=0).
       If empty, rank on y (standard PWR).
       If provided, rank on z (generalised rankvar option).           */
    real scalar    n, h, tmp, q25, q75, sd_pc, i
    real colvector ord, fw_s, pc_s, pc_sorted, u, w_s, w, pc, z_use, z_s

    n         = rows(y)

    /* Use z if provided, else fall back to y */
    z_use     = (rows(zrank) == n) ? zrank : y

    ord       = order(z_use, 1)       /* sort on z (or y if default) */
    fw_s      = fw[ord]
    pc_s      = runningsum(fw_s) :/ sum(fw_s)

    /* Ties: assign each group of equal z the SAME rank = F_n(z) (paper eq 3),
       so p-hat depends only on the value, not on the sort order of ties.
       Removes run-to-run drift in beta and the MSE-optimal bandwidth.       */
    z_s       = z_use[ord]
    i         = n - 1
    while (i >= 1) {
        if (z_s[i] == z_s[i+1]) pc_s[i] = pc_s[i+1]
        i--
    }

    pc_sorted = sort(pc_s, 1)
    q25       = pc_sorted[ceil(0.25 * n)]
    q75       = pc_sorted[ceil(0.75 * n)]
    tmp       = (q75 - q25) / 1.34
    sd_pc     = sqrt(variance(pc_s))
    if (tmp > sd_pc) tmp = sd_pc
    h         = (band_in == 0) ? cband * tmp * n^(-0.2) : band_in

    u         = (pc_s :- tau) :/ h
    w_s       = exp(-0.25 :* u:^2) :/ (h * sqrt(2 * pi()) * n)

    w         = J(n, 1, 0)
    pc        = J(n, 1, 0)
    w[ord]    = w_s
    pc[ord]   = pc_s

    st_numscalar("_gepwreg_h_tmp", h)
    return((w, pc))
}


/* ─────────────────────────────────────────────────────────────────────────────
   _gepwreg_hopt()
   MSE-optimal bandwidth via two-step plug-in (see paper, eq. h*).

     h*(tau) = [ sigma2(tau) / (8 * n * sqrt(2*pi) * [mu''(tau)]^2) ]^(1/5)

   Algorithm (paper Algorithm 1):
     1. Pilot h_pilot = 2 * h_Silverman
     2. Estimate local-constant fit of y on percentile rank at
        tau-d, tau, tau+d using h_pilot  (d = 0.02)
     3. mu''(tau) ~ [mu(tau+d) - 2*mu(tau) + mu(tau-d)] / d^2
     4. sigma2(tau) = weighted residual variance at tau (from step 2)
     5. h* clipped to [0.3*h_Silverman, 5*h_Silverman] for stability
        (curvature near zero -> h* would otherwise diverge)

   NOTE: h* is computed using y (or zrank) curvature only -- a single
   scalar bandwidth is then applied to ALL regressors, consistent with
   how Silverman's rule is also depvar-based in this module.
   ───────────────────────────────────────────────────────────────────────── */
real scalar _gepwreg_hopt(real colvector y,
                           real colvector fw,
                           real scalar    tau,
                           real scalar    cband,
                           real scalar    n)
{
    real scalar    h_sil, h_pilot, d, mu_lo, mu_mid, mu_hi, curv
    real scalar    sigma2, h_star, lo_clip, hi_clip
    real matrix    wp_sil, wp_lo, wp_mid, wp_hi
    real colvector w_sil, w_lo, w_mid, w_hi

    /* Step 0: pilot Silverman bandwidth (band_in=0 triggers Silverman) */
    wp_sil = _gepwreg_wp(y, fw, tau, cband, 0, J(0,1,.))
    h_sil  = st_numscalar("_gepwreg_h_tmp")
    h_pilot = 2 * h_sil

    d = 0.02

    /* Step 1: local-constant fits (weighted mean of y) at tau-d, tau, tau+d */
    wp_lo  = _gepwreg_wp(y, fw, tau - d, cband, h_pilot, J(0,1,.))
    w_lo   = wp_lo[., 1]
    mu_lo  = sum(w_lo :* y) / sum(w_lo)

    wp_mid = _gepwreg_wp(y, fw, tau,     cband, h_pilot, J(0,1,.))
    w_mid  = wp_mid[., 1]
    mu_mid = sum(w_mid :* y) / sum(w_mid)

    wp_hi  = _gepwreg_wp(y, fw, tau + d, cband, h_pilot, J(0,1,.))
    w_hi   = wp_hi[., 1]
    mu_hi  = sum(w_hi :* y) / sum(w_hi)

    /* Step 2: curvature (second derivative) of local mean profile */
    curv = (mu_hi - 2*mu_mid + mu_lo) / d^2

    /* Step 3: local residual variance at tau (using pilot weights) */
    sigma2 = sum(w_mid :* (y :- mu_mid):^2) / sum(w_mid)

    /* Step 4: MSE-optimal h* (paper eq. h*); guard against curv ~ 0 */
    if (abs(curv) < 1e-8) curv = 1e-8
    h_star = (sigma2 / (8 * n * sqrt(2*pi()) * curv^2))^(0.2)

    /* Step 5: clip for stability (paper Algorithm 1, step 5) */
    lo_clip = 0.3 * h_sil
    hi_clip = 5.0 * h_sil
    if (h_star < lo_clip) h_star = lo_clip
    if (h_star > hi_clip) h_star = hi_clip

    return(h_star)
}


/* ─────────────────────────────────────────────────────────────────────────────
   _gepwreg_wls() : beta = (X'WX)^{-1} X'Wy
   ───────────────────────────────────────────────────────────────────────── */
real colvector _gepwreg_wls(real matrix X, real colvector y, real colvector w)
{
    /* invsym() is more stable than lusolve() for near-singular X'WX:
       returns generalised inverse (zeros for collinear cols) instead
       of missing values, which would crash ereturn post.             */
    real matrix XtW, A
    XtW = (X :* w)'
    A   = XtW * X
    return(invsym(A) * (XtW * y))
}


/* ─────────────────────────────────────────────────────────────────────────────
   _gepwreg_revCumSum()
   Reverse cumulative row sum of matrix M :
     R[j,.] = M[j,.] + M[j+1,.] + ... + M[n,.]

   Algorithm :
     Forward cumsum column-by-column (quadrunningsum requires a vector).
     revCumSum[i,.] = totalSum - forwardCumSum[i-1,.]
     where forwardCumSum[0,.] = 0 (prepended row of zeros).
   Complexity : O(n*k), one while-loop over k columns
   ───────────────────────────────────────────────────────────────────────── */
real matrix _gepwreg_revCumSum(real matrix M)
{
    real scalar n, k, j
    real matrix fcs, total_row

    n         = rows(M)
    k         = cols(M)
    fcs       = J(n, k, 0)

    /* forward cumulative sum, one column at a time */
    j = 1
    while (j <= k) {
        fcs[., j] = quadrunningsum(M[., j])
        j++
    }

    /* total row sum (1 x k) */
    total_row = fcs[n, .]

    /* revCumSum[i,.] = total - fcs[i-1,.]  (fcs[0,.] = 0 by convention) */
    return(total_row :- (J(1, k, 0) \ fcs[1..(n-1), .]))
}


/* ---------------------------------------------------------------------------
   _gepwreg_revCumSumT()   (1.4)
   The same reverse cumulative sum, but over the SET { i : p_i >= p_j } rather
   than over the positions below j.  The two coincide only when the ranking
   variable has no ties.

   Why it matters.  _gepwreg_wp already gives every observation of a tie group
   the same rank p (paper eq. 3), so the kernel weights, and with them the
   point estimates, do not depend on the order of the rows.  The influence
   function did: it is a sum of per-observation scores, and a positional
   cumulative sum splits a tie group into unequal partial sums according to
   which member the sort happened to put first.  Every standard error then
   moved when the rows were permuted -- by 8e-04 on a ranking variable with
   seven distinct values over 8,478 observations, and by 7e-07 on a
   continuous one.  An estimate must not depend on a permutation of the data.

   p is exactly constant inside a tie group and strictly increasing across
   groups, so equality of p identifies the groups and the ranking variable
   itself need not be passed in.  p must arrive in the sorted order of M.
   --------------------------------------------------------------------------- */
real matrix _gepwreg_revCumSumT(real matrix M, real colvector p)
{
    real scalar n, i
    real matrix R
    R = _gepwreg_revCumSum(M)
    n = rows(M)
    for (i = 2; i <= n; i++) {
        if (p[i] == p[i-1]) R[i, .] = R[i-1, .]
    }
    return(R)
}


/* ─────────────────────────────────────────────────────────────────────────────
   _gepwreg_se_IF()
   Influence-function variance (Deville 1999 / Newey-McFadden 1994)

   Total influence of observation j :
     psi_j = w_j * x_j * e_j
             + (1/n) * sum_{i: y_i >= y_j}  dw_i * x_i * e_i

   Kernel derivative :
     dw_i = -0.5 * (pc_i - tau) / h^2 * w_i

   Variance :
     V = A^{-1} * (psi'psi / n) * A^{-1} / n,   A = X'WX / n
   ───────────────────────────────────────────────────────────────────────── */
real matrix _gepwreg_se_IF(real matrix  X,
                   real colvector y,
                   real colvector w,
                   real colvector pc,
                   real scalar    h,
                   real scalar    tau,
                   real colvector beta,
                   real colvector z_ord)
{
    /* z_ord: ordering vector based on ranking variable z.
       If rows(z_ord)=0, fall back to ordering by y (default). */
    real scalar    n, k
    real colvector e, dw, ord
    real matrix    A, Ainv, S_dir, dS, revCS, psi, V

    n          = rows(y)
    k          = cols(X)
    e          = y - X * beta
    A          = (X :* w)' * X / n
    Ainv       = invsym(A)   /* more stable than luinv for near-singular A */
    dw         = (-0.5 / h^2) :* (pc :- tau) :* w
    S_dir      = X :* (w :* e)
    dS         = X :* (dw :* e)
    /* Sort on z if provided, else sort on y */
    ord        = (rows(z_ord) == n) ? z_ord : order(y, 1)
    revCS      = J(n, k, 0)
    revCS[ord,] = _gepwreg_revCumSumT(dS[ord,], pc[ord])
    psi        = S_dir + revCS :/ n
    /* (1.4) centred scores: the indirect term carries the constant
       -(1/n) sum_i dw_i x_i e_i F(y_i), omitted up to 1.3             */
    psi        = psi :- mean(psi)
    V          = Ainv * (psi' * psi / n) * Ainv / n
    /* Warn if any diagonal of A is near zero (sparse category at tau) */
    if (min(diagonal(A)) < 1e-12) {
        printf("{txt}Warning: near-singular X'WX at tau=%g -- SE unreliable for sparse categories.\n", tau)
    }
    return(V)
}


/* ─────────────────────────────────────────────────────────────────────────────
   _gepwreg_se_naive()
   Standard WLS variance assuming fixed weights.
   INCONSISTENT for PWR -- stored in e(V_naive) for reference only.
   ───────────────────────────────────────────────────────────────────────── */
real matrix _gepwreg_se_naive(real matrix  X,
                      real colvector y,
                      real colvector w,
                      real colvector beta)
{
    real scalar    n, k, s2
    real colvector e
    real matrix    Ainv

    n    = rows(y)
    k    = cols(X)
    e    = y - X * beta
    Ainv = invsym((X :* w)' * X / n)   /* stable for sparse categories */
    s2   = sum(w :* e:^2) / (n - k)
    return(s2 * Ainv / n)
}


/* ─────────────────────────────────────────────────────────────────────────────
   _gepwreg_se_boot()
   Bootstrap standard errors.

   If survey design is declared (do_svy=1):
     Cluster (PSU) bootstrap stratifié — equivalent to:
       bootstrap, strata(strata) cluster(psu) reps(B): gepwreg ...
     Algorithm:
       For each stratum h with n_h PSUs:
         Draw n_h PSUs WITH REPLACEMENT from stratum h
         Stack ALL observations from drawn PSUs (PSU drawn k times → k copies)
       Run gepwreg on the stacked sample with RECOMPUTED kernel weights
     This exactly replicates what Stata's bootstrap prefix does.

   If no survey design (do_svy=0):
     Pairs bootstrap — observations drawn with replacement individually.
     Kernel weights RECONSTRUCTED at every draw.
   ─────────────────────────────────────────────────────────────────────────── */
real matrix _gepwreg_se_boot(real matrix   X,
                     real colvector y,
                     real colvector fw,
                     real colvector w0,
                     real scalar    tau,
                     real scalar    cband,
                     real scalar    band_in,
                     real colvector beta0,
                     real scalar    B,
                     real scalar    do_svy,
                     string scalar  touse,
                     string scalar  svy_psu,
                     string scalar  svy_strata)
{
    real scalar    n, k, b, h_val, p_val, i_h, i_p, n_h, n_b
    real colvector idx, y_b, fw_b, w_b, drawn
    real colvector psu_vec, strata_vec, strata_ids, psu_in_h
    real colvector mask_h, mask_p, obs_idx, idx_b
    real matrix    BB, wp_b, X_b

    n      = rows(y)
    k      = rows(beta0)
    BB     = J(B, k, .)
    rseed(12345)

    if (do_svy & svy_psu != "" & svy_strata != "") {
        /* ── Cluster (PSU) bootstrap within strata ──────────────────────
           Mirrors: bootstrap, strata(strata) cluster(psu): gepwreg ...
           For each replication:
             1. For each stratum: draw n_h PSUs with replacement
             2. Stack observations from selected PSUs
             3. Recompute kernel weights on stacked sample
             4. Run WLS on stacked sample                                */
        st_view(psu_vec,    ., svy_psu,    touse)
        st_view(strata_vec, ., svy_strata, touse)
        strata_ids = uniqrows(strata_vec)

        /* (1.4) the standard progress dots: a bootstrap of 200 replications
           with nothing on the screen looks like a hung session.  _dots is
           silenced automatically under -quietly-, so no guard is needed.  */
        if (_gepwreg_gmw_dots_on()) {
            stata("_dots 0, title(Bootstrap replications) reps(" +
                  strofreal(B) + ")")
            displayflush()
        }
        b = 1
        while (b <= B) {
            /* Build index of selected observations for this replication */
            idx_b = J(0, 1, .)

            i_h = 1
            while (i_h <= rows(strata_ids)) {
                h_val    = strata_ids[i_h]
                mask_h   = (strata_vec :== h_val)
                psu_in_h = uniqrows(select(psu_vec, mask_h))
                n_h      = rows(psu_in_h)

                /* Draw n_h PSUs with replacement */
                drawn = ceil(n_h :* runiform(n_h, 1))

                /* Append observations from each drawn PSU */
                i_p = 1
                while (i_p <= n_h) {
                    p_val   = psu_in_h[drawn[i_p]]
                    mask_p  = mask_h :& (psu_vec :== p_val)
                    obs_idx = select((1::n), mask_p)
                    idx_b   = idx_b \ obs_idx
                    i_p++
                }
                i_h++
            }

            /* Bootstrap sample */
            n_b  = rows(idx_b)
            y_b  = y[idx_b]
            fw_b = fw[idx_b]
            X_b  = X[idx_b, .]

            /* Recompute kernel weights on bootstrap sample */
            wp_b    = _gepwreg_wp(y_b, fw_b, tau, cband, band_in, J(0,1,.))
            w_b     = wp_b[., 1]
            BB[b, ] = _gepwreg_wls(X_b, y_b, w_b)'
            if (_gepwreg_gmw_dots_on()) {
                stata("_dots " + strofreal(b) + " 0")
                displayflush()
            }
            b++
        }
    }
    else {
        /* ── Pairs bootstrap (no survey design declared) ─────────────── */
        if (_gepwreg_gmw_dots_on()) {
            stata("_dots 0, title(Bootstrap replications) reps(" +
                  strofreal(B) + ")")
            displayflush()
        }
        b = 1
        while (b <= B) {
            idx     = ceil(n :* runiform(n, 1))
            y_b     = y[idx]
            fw_b    = fw[idx]
            wp_b    = _gepwreg_wp(y_b, fw_b, tau, cband, band_in, J(0,1,.))
            w_b     = wp_b[., 1]
            BB[b, ] = _gepwreg_wls(X[idx, ], y_b, w_b)'
            if (_gepwreg_gmw_dots_on()) {
                stata("_dots " + strofreal(b) + " 0")
                displayflush()
            }
            b++
        }
    }
    return(variance(BB))
}


/* ─────────────────────────────────────────────────────────────────────────────
   _gepwreg_se_svy()
   Taylor linearisation variance under stratified cluster design.
     V_svy = A^{-1} [ sum_h n_h/(n_h-1) sum_i (z_hi-z_bar_h)(z_hi-z_bar_h)' ] A^{-1} / n^2
   where z_hi = sum_{j in PSU(h,i)} psi_j  (k-vector of IF scores).
   ───────────────────────────────────────────────────────────────────────── */
real matrix _gepwreg_se_svy(real matrix  X,
                    real colvector y,
                    real colvector w,
                    real colvector pc,
                    real scalar    h,
                    real scalar    tau,
                    real colvector beta,
                    string scalar  touse,
                    string scalar  svy_psu,
                    string scalar  svy_strata,
                    real scalar    n,
                    real colvector z_ord)
{
    real scalar    k, n_h
    real matrix    psi, A, Ainv, meat, dev, Z_h, S_dir, dS, revCS
    real colvector psu_vec, strata_vec, strata_ids, psu_in_h
    real colvector z_bar, mask_h, mask_p, e, dw, ord
    real scalar    h_val, p_val, i_h, i_p

    k    = cols(X)
    A    = (X :* w)' * X / n
    Ainv = invsym(A)   /* stable for near-singular A (sparse PSU or category) */

    /* ── Raw IF scores psi (n x k) ───────────────────────────────────── */
    e     = y - X * beta
    dw    = (-0.5 / h^2) :* (pc :- tau) :* w
    S_dir = X :* (w :* e)
    dS    = X :* (dw :* e)
    ord   = (rows(z_ord) == n) ? z_ord : order(y, 1)
    revCS = J(n, k, 0)
    revCS[ord,] = _gepwreg_revCumSumT(dS[ord,], pc[ord])
    psi   = S_dir + revCS :/ n
    psi   = psi :- mean(psi)          /* (1.4) centred scores */

    /* ── Load PSU and strata identifiers ─────────────────────────────── */
    if (svy_psu != "") {
        st_view(psu_vec, ., svy_psu, touse)
    }
    else {
        /* No PSU declared: treat each observation as its own PSU */
        psu_vec = (1::n)
        printf("{txt}Note: no PSU declared in svyset; treating each obs as its own PSU.\n")
        printf("{txt}      Variance will be conservative. Recommend: svyset psu [pw=fw], strata(strata)\n")
    }

    if (svy_strata != "") {
        st_view(strata_vec, ., svy_strata, touse)
    }
    else {
        /* No strata: single stratum */
        strata_vec = J(n, 1, 1)
    }

    /* ── Taylor meat: sum over strata and PSUs ───────────────────────── */
    meat     = J(k, k, 0)
    strata_ids = uniqrows(strata_vec)

    i_h = 1
    while (i_h <= rows(strata_ids)) {
        h_val   = strata_ids[i_h]
        mask_h  = (strata_vec :== h_val)
        psu_in_h = uniqrows(select(psu_vec, mask_h))
        n_h     = rows(psu_in_h)

        if (n_h < 2) {
            /* Certainty stratum (single PSU): skip */
            i_h++
            continue
        }

        /* Collect z_hi = sum_{j in PSU i, stratum h} psi_j */
        Z_h = J(n_h, k, 0)
        i_p = 1
        while (i_p <= n_h) {
            p_val    = psu_in_h[i_p]
            mask_p   = mask_h :& (psu_vec :== p_val)
            Z_h[i_p,] = colsum(select(psi, mask_p))
            i_p++
        }

        /* Within-stratum deviation */
        z_bar = colsum(Z_h) :/ n_h       /* 1 x k */
        dev   = Z_h :- z_bar                 /* broadcast 1xk to n_h x k */
        meat  = meat + (n_h / (n_h - 1)) * (dev' * dev)
        i_h++
    }

    return(Ainv * meat * Ainv / n^2)
}

void _gepwreg_main(string scalar depvar,
                 string scalar indepvars,
                 string scalar fw_var,
                 string scalar touse,
                 real scalar   tau,
                 real scalar   cband,
                 real scalar   band_in,
                 real scalar   B,
                 real scalar   addcons,
                 real scalar   do_svy,
                 string scalar svy_psu,
                 string scalar svy_strata,
                 string scalar rank_var,
                 real scalar   do_optbw)
{
    real colvector  y, fw, w, pc, beta, zrank, z_ord_v
    real matrix     X, wp, V_IF, V_naive, V_boot, V_svy
    real scalar     n, k, h, n_eff, band_use, z_for_h
    real colvector  zrank_for_h

    st_view(y,  ., depvar,  touse)
    st_view(fw, ., fw_var,  touse)
    n = rows(y)

    /* Load ranking variable (empty if same as depvar) */
    zrank = J(0, 1, .)
    if (rank_var != depvar & rank_var != "") {
        st_view(zrank, ., rank_var, touse)
    }

    if (indepvars != "") {
        st_view(X, ., tokens(indepvars), touse)
        if (addcons) X = (X, J(n, 1, 1))
    }
    else X = J(n, 1, 1)
    k = cols(X)

    /* MSE-optimal bandwidth: compute h* and override band_in           */
    /* Curvature/variance are based on the ranking variable (y or zrank) */
    band_use = band_in
    if (do_optbw) {
        zrank_for_h = (rows(zrank) == n) ? zrank : y
        band_use = _gepwreg_hopt(zrank_for_h, fw, tau, cband, n)
    }

    wp    = _gepwreg_wp(y, fw, tau, cband, band_use, zrank)
    w     = wp[., 1]
    pc    = wp[., 2]
    h     = st_numscalar("_gepwreg_h_tmp")
    n_eff = sum(w)^2 / sum(w:^2)
    beta  = _gepwreg_wls(X, y, w)

    /* z-based ordering for IF indirect term */
    z_ord_v = (rows(zrank) == n) ? order(zrank, 1) : order(y, 1)

    /* Standard IF variance */
    V_IF    = _gepwreg_se_IF(X, y, w, pc, h, tau, beta, z_ord_v)
    V_naive = _gepwreg_se_naive(X, y, w, beta)

    st_matrix("_gepwreg_Vi", V_IF)          /* (1.4) e(V_IF) always stored */
    /* Survey Taylor variance — overrides V_IF as main SE */
    if (do_svy) {
        V_svy = _gepwreg_se_svy(X, y, w, pc, h, tau, beta,
                        touse, svy_psu, svy_strata, n, z_ord_v)
        st_matrix("_gepwreg_Vs", V_svy)
        V_IF = V_svy
    }

    /* Bootstrap variance */
    if (B > 0) {
        V_boot = _gepwreg_se_boot(X, y, fw, w, tau, cband, band_in, beta, B,
                          do_svy, touse, svy_psu, svy_strata)
        st_matrix("_gepwreg_Vb", V_boot)
    }

    st_matrix("_gepwreg_b",  beta')
    st_matrix("_gepwreg_V",  V_IF)
    st_matrix("_gepwreg_Vn", V_naive)
    st_numscalar("_gepwreg_h",    h)
    st_numscalar("_gepwreg_neff", n_eff)
    st_numscalar("_gepwreg_n",    n)
}


/* ═════════════════════════════════════════════════════════════════════════════
   TWO-STEP ESTIMATOR

   theta_j(tau) = sum_i fw_i w_i b_ij / sum_i fw_i w_i,
   where b_ij is the effect of regressor j for unit i from a heterogeneity
   model estimated on the whole sample, and w_i the percentile weights on
   the rank of y (or of y0 for the initial-quantile object).  The weights are
   the estimator of omega(x) = f(q_tau | x)/f(q_tau) in Firpo, Fortin and
   Lemieux's (2009) representation of the unconditional quantile partial
   effect; b_ij are the CQPEs or their observable-driven counterpart.
   ═══════════════════════════════════════════════════════════════════════════ */

/* Taylor (stratified cluster) variance of a statistic whose linearised
   scores are the rows of psi (n x k): meat / n^2, no bread.                */
real matrix _gepwreg_taylor_psi(real matrix psi, string scalar touse,
                                string scalar svy_psu, string scalar svy_strata,
                                real scalar n)
{
    real scalar    k, n_h, h_val, p_val, i_h, i_p
    real matrix    meat, dev, Z_h
    real colvector psu_vec, strata_vec, strata_ids, psu_in_h, mask_h, mask_p
    real rowvector z_bar

    k = cols(psi)
    if (svy_psu != "") st_view(psu_vec, ., svy_psu, touse)
    else               psu_vec = (1::n)
    if (svy_strata != "") st_view(strata_vec, ., svy_strata, touse)
    else                  strata_vec = J(n, 1, 1)
    meat       = J(k, k, 0)
    strata_ids = uniqrows(strata_vec)
    i_h = 1
    while (i_h <= rows(strata_ids)) {
        h_val    = strata_ids[i_h]
        mask_h   = (strata_vec :== h_val)
        psu_in_h = uniqrows(select(psu_vec, mask_h))
        n_h      = rows(psu_in_h)
        if (n_h < 2) {
            i_h++
            continue
        }
        Z_h = J(n_h, k, 0)
        i_p = 1
        while (i_p <= n_h) {
            p_val     = psu_in_h[i_p]
            mask_p    = mask_h :& (psu_vec :== p_val)
            Z_h[i_p,] = colsum(select(psi, mask_p))
            i_p++
        }
        z_bar = colsum(Z_h) :/ n_h
        dev   = Z_h :- z_bar
        meat  = meat + (n_h / (n_h - 1)) * (dev' * dev)
        i_h++
    }
    return(meat / n^2)
}


/* Bootstrap draw: pairs, or PSUs within strata when clus = 1.              */
real colvector _gepwreg_drawidx(real scalar n, real scalar clus,
                                real colvector psu_vec, real colvector strata_vec)
{
    real colvector idx, strata_ids, psu_in_h, mask_h, drawn
    real scalar    i_h, i_p, n_h, h_val, p_val

    if (!clus) return(ceil(n :* runiform(n, 1)))
    idx        = J(0, 1, .)
    strata_ids = uniqrows(strata_vec)
    i_h = 1
    while (i_h <= rows(strata_ids)) {
        h_val    = strata_ids[i_h]
        mask_h   = (strata_vec :== h_val)
        psu_in_h = uniqrows(select(psu_vec, mask_h))
        n_h      = rows(psu_in_h)
        drawn    = ceil(n_h :* runiform(n_h, 1))
        i_p = 1
        while (i_p <= n_h) {
            p_val = psu_in_h[drawn[i_p]]
            idx   = idx \ select((1::n), mask_h :& (psu_vec :== p_val))
            i_p++
        }
        i_h++
    }
    return(idx)
}


/* Linear interpolation of the rows of G (grid x k) at the conditional
   ranks u (n x 1) on the grid ug (1 x grid); u is clipped to the grid.    */
real matrix _gepwreg_interp(real matrix G, real colvector u, real rowvector ug)
{
    real scalar    n, k, ng, g
    real colvector uc, t, sel
    real matrix    Bi

    n  = rows(u)
    k  = cols(G)
    ng = cols(ug)
    uc = rowmin((rowmax((u, J(n, 1, ug[1]))), J(n, 1, ug[ng])))
    Bi = J(n, k, .)
    g  = 1
    while (g < ng) {
        if (g < ng - 1) sel = selectindex((uc :>= ug[g]) :& (uc :<  ug[g+1]))
        else            sel = selectindex((uc :>= ug[g]) :& (uc :<= ug[g+1]))
        if (rows(sel) > 0) {
            t         = (uc[sel] :- ug[g]) :/ (ug[g+1] - ug[g])
            Bi[sel, .] = J(rows(sel), 1, G[g, .]) + t * (G[g+1, .] - G[g, .])
        }
        g++
    }
    return(Bi)
}


/* Variance of bootstrap replicates, rows with missing values dropped.      */
/* ─────────────────────────────────────────────────────────────────────────────
   _gepwreg_kw()  percentile weights at (tau, h) from ranks already computed
   ───────────────────────────────────────────────────────────────────────── */
real colvector _gepwreg_kw(real colvector pc, real scalar tau,
                           real scalar h, real scalar n)
{
    return(exp(-0.25 :* ((pc :- tau) :/ h):^2) :/ (h * sqrt(2*pi()) * n))
}

/* ─────────────────────────────────────────────────────────────────────────────
   _gepwreg_hopt2()
   MSE-optimal bandwidth for the SECOND step of the two-step estimator.

   Step 2 is a Nadaraya-Watson regression of the household effects b_ij on the
   rank of the outcome, evaluated at tau.  The rank is uniform on [0,1], so the
   design density is flat, the usual g' f'/f term of the local constant bias
   vanishes, and only the curvature of the profile survives:

       MSE_j(h) = h^4 theta_j''(tau)^2 + s_j^2(tau) / (2 sqrt(2 pi) n h),
       s_j^2(tau) = Var(b_ij | p_i = tau).

   The coefficients do not share units, so the average RELATIVE mean squared
   error is minimised -- each term divided by s_j^2 -- which gives one
   bandwidth for the call and reduces to the single-coefficient formula when
   k = 1:

       h* = [ k / (8 sqrt(2 pi) n sum_j (theta_j'' / s_j)^2 ) ]^(1/5).

   theta_j'' is obtained by second differences at a pilot bandwidth of
   2 h_Silverman, s_j^2 by the weighted variance of b_ij in the pilot group.
   Because b_ij does not depend on h in either heterogeneity model, the pilot
   costs three weighted means and no refitting.

   h* is clipped to [0.5, 3] x h_Silverman.  The plug-in is pulled towards
   zero when the profile is locally straight, since
   E[thetahat''^2] = theta''^2 + Var(thetahat'') > 0 inflates the denominator,
   while in population h* diverges there; the clipping covers both.
   ───────────────────────────────────────────────────────────────────────── */
real scalar _gepwreg_hopt2(real matrix Bi, real colvector fw,
                           real colvector pc, real scalar tau,
                           real scalar hsil, real scalar n)
{
    real scalar    k, kk, hp, d, agg, hstar, j, S
    real rowvector th0, thp, thm, d2, s2
    real colvector ww

    k  = cols(Bi)
    hp = 2 * hsil
    d  = min((0.05, tau/2, (1-tau)/2))
    if (d <= 0 | hsil <= 0) return(hsil)

    ww  = fw :* _gepwreg_kw(pc, tau, hp, n)
    S   = sum(ww)
    th0 = (ww' * Bi) / S
    s2  = J(1, k, 0)
    for (j = 1; j <= k; j++) s2[j] = sum(ww :* (Bi[., j] :- th0[j]):^2) / S

    ww  = fw :* _gepwreg_kw(pc, tau + d, hp, n)
    thp = (ww' * Bi) / sum(ww)
    ww  = fw :* _gepwreg_kw(pc, tau - d, hp, n)
    thm = (ww' * Bi) / sum(ww)
    d2  = (thp - 2*th0 + thm) :/ d^2

    agg = 0
    kk  = 0
    for (j = 1; j <= k; j++) {
        if (s2[j] > 1e-12) {
            agg = agg + (d2[j]^2) / s2[j]
            kk  = kk + 1
        }
    }
    if (kk == 0 | agg <= 0) return(3 * hsil)

    hstar = (kk / (8 * sqrt(2*pi()) * n * agg))^0.2
    if (hstar < 0.5*hsil) hstar = 0.5*hsil
    if (hstar > 3  *hsil) hstar = 3  *hsil
    return(hstar)
}

real matrix _gepwreg_bootvar(real matrix TT)
{
    real colvector ok
    ok = (rowmissing(TT) :== 0)
    if (sum(ok) < 3) return(J(cols(TT), cols(TT), .))
    return(variance(select(TT, ok)))
}


/* ─────────────────────────────────────────────────────────────────────────────
   (1.5.0) THE REFERENCE HOUSEHOLD, option ref()
   _gepwreg_refrow(): the reference household on any data -- the sample or a
   bootstrap draw -- from the rule of _gepwreg_refrule (see the ado part).
   Returns xr, zr (the values of the regressor and het() terms) and uval (per
   variable: the value, or for a factor the index of the level, 0 = base).
   ───────────────────────────────────────────────────────────────────────── */
void _gepwreg_refrow(real matrix X, real matrix Z, real colvector fw,
                     real colvector rk, real colvector isfac,
                     real matrix RX, real matrix RZ,
                     real rowvector xfix, real rowvector zfix,
                     real rowvector xr, real rowvector zr, real colvector uval)
{
    real scalar    nu, u, val, S, nl, c, best, i, fx
    real rowvector cx, cz, F
    real colvector v
    nu   = rows(rk)
    xr   = xfix
    zr   = zfix
    uval = J(nu, 1, .)
    S    = sum(fw)
    for (u = 1; u <= nu; u++) {
        cx = J(1, 0, .)
        cz = J(1, 0, .)
        if (cols(RX)) cx = selectindex(RX[u, .])
        if (cols(RZ)) cz = selectindex(RZ[u, .])
        nl = cols(cx)
        if (nl == 0) nl = cols(cz)
        if (rk[u] == 0) {
            if (isfac[u]) {
                uval[u] = 0
                for (c = 1; c <= nl; c++) {
                    if (cols(cx)) fx = xfix[cx[c]]
                    else          fx = zfix[cz[c]]
                    if (fx == 1) uval[u] = c
                }
            }
            else {
                if (cols(cx)) uval[u] = xfix[cx[1]]
                else          uval[u] = zfix[cz[1]]
            }
        }
        else if (rk[u] == 4) {
            /* ref(mean): every indicator of the factor at its mean */
            for (c = 1; c <= nl; c++) {
                if (cols(cx)) val = sum(fw :* X[., cx[c]]) / S
                else          val = sum(fw :* Z[., cz[c]]) / S
                if (cols(cx)) xr[cx[c]] = val
                if (cols(cz)) zr[cz[c]] = val
            }
            uval[u] = .
        }
        else if (rk[u] == 3) {
            F = J(1, nl + 1, 0)
            for (c = 1; c <= nl; c++) {
                if (cols(cx)) F[c + 1] = sum(fw :* X[., cx[c]])
                else          F[c + 1] = sum(fw :* Z[., cz[c]])
            }
            F[1] = S - sum(F)
            best = 1
            for (i = 2; i <= nl + 1; i++) {
                if (F[i] > F[best]) best = i
            }
            if (cols(cx)) xr[cx] = J(1, cols(cx), 0)
            if (cols(cz)) zr[cz] = J(1, cols(cz), 0)
            if (best > 1) {
                if (cols(cx)) xr[cx[best - 1]] = 1
                if (cols(cz)) zr[cz[best - 1]] = 1
            }
            uval[u] = best - 1
        }
        else {
            if (cols(cx)) v = X[., cx[1]]
            else          v = Z[., cz[1]]
            val = sum(fw :* v) / S
            if (rk[u] == 2) val = (val > 0.5)
            if (cols(cx)) xr[cx] = J(1, cols(cx), val)
            if (cols(cz)) zr[cz] = J(1, cols(cz), val)
            uval[u] = val
        }
    }
}

/* het(z): the derivative of the step-1 design row with respect to x_j at
   one household (x, z) -- or, for a binary x_j that is also z_l, the 1-vs-0
   difference (review C8) -- so that the effect is M gam'.                */
real rowvector _gepwreg_mrow(real rowvector x, real rowvector z, real scalar j,
                             real matrix same, real scalar isbin)
{
    real scalar    k, m, l, jj, c
    real rowvector M
    k = cols(x)
    m = cols(z)
    M = J(1, k + m + k*m + 1, 0)
    M[j] = 1
    for (l = 1; l <= m; l++) {
        if (same[j, l]) M[k + l] = 1
        for (jj = 1; jj <= k; jj++) {
            c = k + m + (jj-1)*m + l
            if (jj == j) {
                if (same[j, l]) M[c] = (isbin ? 1 : 2 * x[j])
                else            M[c] = z[l]
            }
            else if (same[j, l]) M[c] = x[jj]
        }
    }
    return(M)
}

/* het(z): the effect of the reference household, theta_ref_j = M_j(xr, zr)
   gam', and (dopsi) its influence function: psi_g M_j' plus, for each
   reference value that is a sample mean, d theta_ref / d mean times the
   influence of that mean -- theta_ref is affine in the reference values, so
   the derivative is the difference over one unit, exactly.  Does not depend
   on tau.                                                                 */
real rowvector _gepwreg_ref_z(real matrix X, real matrix Z, real colvector fw,
                              real rowvector gam, real matrix psi_g,
                              real matrix same, real rowvector isbin,
                              real colvector rk, real colvector isfac,
                              real matrix RX, real matrix RZ,
                              real rowvector xfix, real rowvector zfix,
                              real matrix psir, real rowvector xr,
                              real colvector uval, real scalar dopsi)
{
    real scalar    k, n, j, u, d, vbar, S, nl, c
    real rowvector zr, xe, ze, thr, cx, cz
    real colvector v
    k  = cols(X)
    n  = rows(X)
    zr = J(1, 0, .)
    _gepwreg_refrow(X, Z, fw, rk, isfac, RX, RZ, xfix, zfix, xr, zr, uval)
    thr = J(1, k, .)
    for (j = 1; j <= k; j++) thr[j] = _gepwreg_mrow(xr, zr, j, same, isbin[j]) * gam'
    if (!dopsi) return(thr)
    S    = sum(fw)
    psir = J(n, k, 0)
    for (j = 1; j <= k; j++) psir[., j] = psi_g * _gepwreg_mrow(xr, zr, j, same, isbin[j])'
    for (u = 1; u <= rows(rk); u++) {
        if (rk[u] != 1 & rk[u] != 4) continue
        cx = selectindex(RX[u, .])
        cz = J(1, 0, .)
        if (cols(RZ)) cz = selectindex(RZ[u, .])
        /* one mean (rk 1: every position of the variable moves together),
           or one mean per indicator of a factor (rk 4)                     */
        if (rk[u] == 1) nl = 1
        else            nl = max((cols(cx), cols(cz)))
        for (c = 1; c <= nl; c++) {
            xe = xr
            ze = zr
            if (rk[u] == 1) {
                if (cols(cx)) xe[cx] = xe[cx] :+ 1
                if (cols(cz)) ze[cz] = ze[cz] :+ 1
            }
            else {
                if (cols(cx)) xe[cx[c]] = xe[cx[c]] + 1
                if (cols(cz)) ze[cz[c]] = ze[cz[c]] + 1
            }
            if (cols(cx)) v = X[., cx[c]]
            else          v = Z[., cz[c]]
            vbar = sum(fw :* v) / S
            for (j = 1; j <= k; j++) {
                d = _gepwreg_mrow(xe, ze, j, same, isbin[j]) * gam' - thr[j]
                psir[., j] = psir[., j] + (d * n / S) :* fw :* (v :- vbar)
            }
        }
    }
    return(thr)
}

/* het(qr): the estimator E -- beta(u*) with xt' beta(u*) = q on the grid
   path interpolated linearly in u (the constant included in beta and xt);
   seg, tt: the segment [u_seg, u_seg+1] and the position in it; clip = 1
   when q lies beyond the reference's conditional quantiles of the grid.  */
real rowvector _gepwreg_E_pt(real matrix G, real rowvector ug, real scalar q,
                             real rowvector xt, real scalar ustar, real scalar seg,
                             real scalar tt, real scalar clip)
{
    real scalar    ng
    real rowvector Qr
    ng   = cols(ug)
    Qr   = xt * G'
    clip = 0
    if (q <= Qr[1]) {
        seg  = 1
        tt   = 0
        clip = 1
    }
    else if (q >= Qr[ng]) {
        seg  = ng - 1
        tt   = 1
        clip = 1
    }
    else {
        seg = max(selectindex(Qr :<= q))
        tt  = (q - Qr[seg]) / (Qr[seg+1] - Qr[seg])
    }
    ustar = ug[seg] + tt * (ug[seg+1] - ug[seg])
    return((1 - tt) * G[seg, .] + tt * G[seg+1, .])
}

/* the influence function (sum form) of a weighted quantile regression at u,
   with qreg's Powell-kernel density: H xt_i f_i s_i,
   H = ( sum_i f_i K(e_i/kb) xt_i xt_i' / kb )^-1, s_i = u - 1{e_i < 0}    */
real matrix _gepwreg_qrpsi(real matrix Xt, real colvector fw, real colvector r,
                           real colvector kv, real scalar kb, real scalar u)
{
    real matrix    H
    real colvector s
    H = invsym(quadcross(Xt, fw :* kv, Xt) / kb)
    s = (u :- (r :< 0)) :* (abs(r) :>= 1e-10)
    return((Xt :* (fw :* s)) * H)
}

/* het(qr), ref(): the point -- the reference household of the sample, u*, E */
void _gepwreg_qr_ref1(string scalar xvars, string scalar fw_var, string scalar touse,
                      real scalar q, string scalar G_name, string scalar ug_name)
{
    real matrix    X, G, RX
    real colvector fw, rk, isfac, uval
    real rowvector ug, xfix, xr, zr, bs
    real scalar    k, ustar, seg, tt, clip
    st_view(X,  ., tokens(xvars), touse)
    st_view(fw, ., fw_var, touse)
    k     = cols(X)
    G     = st_matrix(G_name)
    ug    = st_matrix(ug_name)
    rk    = st_matrix("_gepwreg_rk")
    isfac = st_matrix("_gepwreg_rfac")
    RX    = st_matrix("_gepwreg_rx")
    xfix  = st_matrix("_gepwreg_xfix")
    xr    = J(1, 0, .)
    zr    = J(1, 0, .)
    uval  = J(0, 1, .)
    _gepwreg_refrow(X, J(rows(X), 0, .), fw, rk, isfac, RX, J(rows(rk), 0, .),
                    xfix, J(1, 0, .), xr, zr, uval)
    ustar = .
    seg   = .
    tt    = .
    clip  = .
    bs = _gepwreg_E_pt(G, ug, q, (xr, 1), ustar, seg, tt, clip)
    st_matrix("_gepwreg_bref", bs[1..k])
    st_matrix("_gepwreg_refx", xr)
    st_matrix("_gepwreg_refu", uval')
    st_numscalar("_gepwreg_ustar", ustar)
    st_numscalar("_gepwreg_seg",   seg)
    st_numscalar("_gepwreg_tt",    tt)
    st_numscalar("_gepwreg_clip",  clip)
}

/* het(qr), ref(): the delta-method influence function of E (n form),
     psi_j = psi(beta*_j) + (dbeta_j / D) [ psi(q) - xt' psi(beta*)
             - sum over the reference means of beta*_x psi(mean) ],
   beta* = (1-t) beta_seg + t beta_seg+1, D = xt'(beta_seg+1 - beta_seg);
   the variance, its Taylor version under a design, and the SE of u*.     */
void _gepwreg_qr_ref2(string scalar xvars, string scalar fw_var, string scalar touse,
                      string scalar G_name, string scalar ug_name, real scalar tau,
                      string scalar r1, string scalar k1, real scalar kb1,
                      string scalar r2, string scalar k2, real scalar kb2,
                      string scalar rq, string scalar kq, real scalar kbq,
                      real scalar do_svy, string scalar svy_psu, string scalar svy_strata)
{
    real matrix    X, Xt, G, RX, P1, P2, Pstar, psi, V
    real colvector fw, rk, core, v, rr, kk
    real rowvector ug, xr, xt, bs, cx
    real scalar    n, k, seg, tt, clip, D, j, u, S, vbar, du
    st_view(X,  ., tokens(xvars), touse)
    st_view(fw, ., fw_var, touse)
    n    = rows(X)
    k    = cols(X)
    Xt   = (X, J(n, 1, 1))
    G    = st_matrix(G_name)
    ug   = st_matrix(ug_name)
    rk   = st_matrix("_gepwreg_rk")
    RX   = st_matrix("_gepwreg_rx")
    xr   = st_matrix("_gepwreg_refx")
    xt   = (xr, 1)
    seg  = st_numscalar("_gepwreg_seg")
    tt   = st_numscalar("_gepwreg_tt")
    clip = st_numscalar("_gepwreg_clip")
    bs   = (1 - tt) * G[seg, .] + tt * G[seg + 1, .]
    st_view(rr, ., r1, touse)
    st_view(kk, ., k1, touse)
    P1 = n * _gepwreg_qrpsi(Xt, fw, rr, kk, kb1, ug[seg])
    st_view(rr, ., r2, touse)
    st_view(kk, ., k2, touse)
    P2 = n * _gepwreg_qrpsi(Xt, fw, rr, kk, kb2, ug[seg + 1])
    st_view(rr, ., rq, touse)
    st_view(kk, ., kq, touse)
    core  = n * _gepwreg_qrpsi(J(n, 1, 1), fw, rr, kk, kbq, tau)
    Pstar = (1 - tt) * P1 + tt * P2
    core  = core - Pstar * xt'
    S = sum(fw)
    for (u = 1; u <= rows(rk); u++) {
        if (rk[u] == 1) {
            cx   = selectindex(RX[u, .])
            v    = X[., cx[1]]
            vbar = sum(fw :* v) / S
            core = core - ((n / S) * sum(bs[cx])) :* fw :* (v :- vbar)
        }
        else if (rk[u] == 4) {
            /* ref(mean): each indicator of the factor at its own mean */
            cx = selectindex(RX[u, .])
            for (j = 1; j <= cols(cx); j++) {
                v    = X[., cx[j]]
                vbar = sum(fw :* v) / S
                core = core - ((n / S) * bs[cx[j]]) :* fw :* (v :- vbar)
            }
        }
    }
    D = xt * (G[seg + 1, .] - G[seg, .])'
    if (clip) core = J(n, 1, 0)
    psi = J(n, k, .)
    for (j = 1; j <= k; j++) psi[., j] = Pstar[., j] + ((G[seg + 1, j] - G[seg, j]) / D) * core
    psi = psi :- mean(psi)
    V   = quadcross(psi, psi) / n^2
    st_matrix("_gepwreg_Vref_IF", V)
    if (do_svy) st_matrix("_gepwreg_Vref_svy", _gepwreg_taylor_psi(psi, touse, svy_psu, svy_strata, n))
    /* the SE of u*: du* = (u_seg+1 - u_seg) dt, dt = core / D */
    du = (ug[seg + 1] - ug[seg]) / D
    v  = du :* core
    v  = v :- mean(v)
    st_numscalar("_gepwreg_seu", sqrt(quadcross(v, v)) / n)
}

/* het(qr), ref(): E in a bootstrap draw -- the reference household of the
   draw, u*, beta(u*) (values only)                                      */
real rowvector _gepwreg_qr_refb(string scalar depvar, string scalar xvars, string scalar fw_var,
                                string scalar touse, real scalar tau, string scalar G_name,
                                string scalar ug_name)
{
    real matrix    X, G, RX
    real colvector fw, rk, isfac, uval, y
    real rowvector ug, xfix, xr, zr, bs
    real scalar    k, ustar, seg, tt, clip, q
    st_view(y,  ., depvar, touse)
    st_view(X,  ., tokens(xvars), touse)
    st_view(fw, ., fw_var, touse)
    q     = _gepwreg_qtau(y, fw, tau)
    k     = cols(X)
    G     = st_matrix(G_name)
    ug    = st_matrix(ug_name)
    rk    = st_matrix("_gepwreg_rk")
    isfac = st_matrix("_gepwreg_rfac")
    RX    = st_matrix("_gepwreg_rx")
    xfix  = st_matrix("_gepwreg_xfix")
    xr    = J(1, 0, .)
    zr    = J(1, 0, .)
    uval  = J(0, 1, .)
    _gepwreg_refrow(X, J(rows(X), 0, .), fw, rk, isfac, RX, J(rows(rk), 0, .),
                    xfix, J(1, 0, .), xr, zr, uval)
    ustar = .
    seg   = .
    tt    = .
    clip  = .
    bs = _gepwreg_E_pt(G, ug, q, (xr, 1), ustar, seg, tt, clip)
    return(bs[1..k])
}


/* Step-1 design of het(z): (X, Z, X x Z, 1); the interaction of x_j with
   z_l sits in column k + m + (j-1)*m + l.                                  */
real matrix _gepwreg_Wmat(real matrix X, real matrix Z)
{
    real scalar n, k, m, j, l
    real matrix W
    n = rows(X)
    k = cols(X)
    m = cols(Z)
    W = (X, Z, J(n, k*m, .), J(n, 1, 1))
    j = 1
    while (j <= k) {
        l = 1
        while (l <= m) {
            W[., k + m + (j-1)*m + l] = X[., j] :* Z[., l]
            l++
        }
        j++
    }
    return(W)
}


/* Unit-level effects of het(z): the derivative of the step-1 model with
   respect to x_j,
     b_ij = g_j + sum_l d_jl z_il
          + sum_{l: z_l is x_j} [ c_l + sum_j' d_j'l x_ij' ]
   (n x k).  The second line is the part of the derivative that runs through
   z_l when z_l is the regressor itself: its own coefficient c_l and the
   interactions of every regressor with it (for j' = j this doubles the
   quadratic term, d/dx of d x^2 = 2 d x).
   (1.5.0, review C8) For a BINARY x_j that is also z_l, x_j^2 = x_j: the
   doubled term 2 d_jl x_ij is not identified (the coefficients of x_j, z_l
   and x_j z_l are a generalized-inverse choice), and the effect is the
   1-vs-0 difference, whose own term is d_jl: 2 d_jl x_ij + d_jl (1 - 2 x_ij)
   = d_jl.  isbin[j] = 1 for a regressor taking only the values 0 and 1.   */
real matrix _gepwreg_betai(real matrix X, real matrix Z, real rowvector gam,
                           real scalar k, real matrix same, real rowvector isbin)
{
    real scalar    n, m, j, l, jj, i0
    real matrix    Bi
    real colvector dcol, add
    n  = rows(Z)
    m  = cols(Z)
    Bi = J(n, k, .)
    j = 1
    while (j <= k) {
        i0   = k + m + (j-1)*m
        dcol = gam[(i0+1)..(i0+m)]'
        Bi[., j] = J(n, 1, gam[j]) + Z * dcol
        l = 1
        while (l <= m) {
            if (same[j, l] == 1) {
                add = J(n, 1, gam[k + l])
                jj = 1
                while (jj <= k) {
                    add = add + X[., jj] * gam[k + m + (jj-1)*m + l]
                    jj++
                }
                Bi[., j] = Bi[., j] + add
                if (isbin[j]) Bi[., j] = Bi[., j] + gam[i0 + l] :* (1 :- 2 :* X[., j])
            }
            l++
        }
        j++
    }
    return(Bi)
}



/* ─────────────────────────────────────────────────────────────────────────────
   (1.5.0) THE PROFILE ACROSS PERCENTILES (graph()): helpers
   ───────────────────────────────────────────────────────────────────────── */
/* q_tau: the tau-quantile of the weighted distribution of y, the smallest y
   whose weighted ECDF reaches tau (ties share the rank, as in _gepwreg_wp) */
real scalar _gepwreg_qtau(real colvector y, real colvector fw, real scalar tau)
{
    real colvector o, ys, cs
    o  = order(y, 1)
    ys = y[o]
    cs = quadrunningsum(fw[o]) :/ sum(fw)
    return(min(select(ys, cs :>= tau - 1e-12)))
}

/* the influence function (n form) of the kernel mean of v at tau: the
   direct term and the term of the estimated ranks (as for psi_mu in
   _gepwreg_ts_z); ord = order(y, 1)                                       */
real colvector _gepwreg_ifmean(real colvector v, real colvector fw, real colvector w,
                               real colvector pc, real scalar tau, real scalar h,
                               real colvector ord)
{
    real scalar    n, S, mu
    real colvector dw, a, rc
    n   = rows(v)
    S   = sum(fw :* w)
    mu  = sum(fw :* w :* v) / S
    dw  = (-0.5 / h^2) :* (pc :- tau) :* w
    a   = dw :* fw :* (v :- mu)
    rc  = J(n, 1, 0)
    rc[ord] = _gepwreg_revCumSumT(a[ord], pc[ord])
    return((n :* fw :* w :* (v :- mu)) :/ S + (fw :* (rc :- sum(a :* pc))) :/ S)
}

/* sup-t critical value of a curve with covariance V across the percentiles:
   the level quantile of max_t |Z_t| / sd_t, Z ~ N(0, V), M draws (V may be
   singular: symmetric square root)                                       */
real scalar _gepwreg_supt(real matrix V, real scalar level, real scalar M)
{
    real matrix    U, Z
    real rowvector lam, mx
    real colvector sd
    real scalar    T
    T  = rows(V)
    sd = sqrt(diagonal(V))
    if (min(sd) <= 0 | hasmissing(sd)) return(.)
    symeigensystem(V, U, lam)
    lam = lam :* (lam :> 0)
    Z   = (U :* sqrt(lam)) * rnormal(T, M, 0, 1)
    mx  = colmax(abs(Z :/ sd))
    mx  = sort(mx', 1)'
    return(mx[ceil(M * level / 100)])
}

/* the flatness test of a curve th (T x 1) with covariance V: H0 th constant
   across the percentiles.  D = th - mean(th), stat = max_t |D_t| / sd(D_t),
   its null distribution N(0, A V A'), A = I - J/T, M draws.  (stat, p)    */
real rowvector _gepwreg_flat(real colvector th, real matrix V, real scalar M)
{
    real scalar    T, stat
    real matrix    A, VD, U, Z
    real colvector D, sd
    real rowvector lam, mx
    T  = rows(th)
    A  = I(T) :- J(T, T, 1/T)
    D  = A * th
    VD = A * V * A'
    sd = sqrt(abs(diagonal(VD)))
    if (min(sd) <= 0 | hasmissing(sd)) return((., .))
    stat = max(abs(D :/ sd))
    symeigensystem(VD, U, lam)
    lam = lam :* (lam :> 0)
    Z   = (U :* sqrt(lam)) * rnormal(T, M, 0, 1)
    mx  = colmax(abs(Z :/ sd))
    return((stat, (1 + sum(mx :>= stat)) / (M + 1)))
}

/* analytic: the profile Pb (T x k), the covariance VP of its Tk scores
   ((g-1)k + j), the first nuse percentiles for the band and the test.
   Stores pre+b, pre+se, pre+cu, pre+flat.  The simulation runs on its own
   seed and gives the random-number state back.                           */
void _gepwreg_prof_stats(real matrix Pb, real matrix VP, real scalar level,
                         real scalar seed, real scalar nuse, string scalar pre)
{
    real scalar    T, k, j
    real matrix    se, Vj, fl
    real rowvector cu, idx, u
    string scalar  rs
    T  = rows(Pb)
    k  = cols(Pb)
    se = J(T, k, .)
    cu = J(1, k, .)
    fl = J(2, k, .)
    u  = 1..nuse
    rs = rngstate()
    rseed(seed)
    for (j = 1; j <= k; j++) {
        idx = (0..(T-1)) :* k :+ j
        Vj  = VP[idx, idx]
        se[., j] = sqrt(abs(diagonal(Vj)))
        cu[j]    = _gepwreg_supt(Vj[u, u], level, 10000)
        fl[., j] = _gepwreg_flat(Pb[u, j], Vj[u, u], 10000)'
    }
    rngstate(rs)
    st_matrix(pre + "b",    Pb)
    st_matrix(pre + "se",   se)
    st_matrix(pre + "cu",   cu)
    st_matrix(pre + "flat", fl)
}

/* bootstrap: the profile Pb (T x k) and its draws TT (B x Tk, (g-1)k + j).
   Stores pre+se (the SD of the draws), pre+cu (the level quantile of
   max_t |th_b - th| / se over the first nuse percentiles) and pre+flat
   (stat, p of the flatness test: the draws of the deviations from their
   own mean, centred at the estimate's).                                  */
void _gepwreg_prof_bstats(real matrix Pb, real matrix TT, real scalar level,
                          real scalar nuse, string scalar pre)
{
    real scalar    T, k, j, B, stat
    real matrix    se, TH, Db, fl
    real rowvector cu, idx, u, sdD
    real colvector mx, sb, thc
    T  = rows(Pb)
    k  = cols(Pb)
    TT = select(TT, rowmissing(TT) :== 0)
    B  = rows(TT)
    se = J(T, k, .)
    cu = J(1, k, .)
    fl = J(2, k, .)
    u  = 1..nuse
    for (j = 1; j <= k; j++) {
        idx = (0..(T-1)) :* k :+ j
        TH  = TT[., idx]
        if (B < 3) continue
        se[., j] = sqrt(diagonal(variance(TH)))
        if (min(se[u, j]) <= 0) continue
        mx = rowmax(abs(TH[., u] :- Pb[u, j]') :/ se[u, j]')
        mx = sort(mx, 1)
        cu[j] = mx[ceil(B * level / 100)]
        thc = Pb[u, j] :- mean(Pb[u, j])
        Db  = (TH[., u] :- rowsum(TH[., u]) :/ nuse) :- thc'
        sdD = sqrt(diagonal(variance(Db)))'
        if (min(sdD) <= 0) continue
        stat = max(abs(thc') :/ sdD)
        sb   = rowmax(abs(Db) :/ sdD)
        fl[., j] = (stat \ (1 + sum(sb :>= stat)) / (B + 1))
    }
    st_matrix(pre + "se",   se)
    st_matrix(pre + "cu",   cu)
    st_matrix(pre + "flat", fl)
}

/* het(qr): the profile on the data in memory (the sample or a draw): the
   households' effects b_i = beta(u_i) once, the ranks once, then at every
   percentile of the grid the bandwidth (the plug-in recomputed) and theta,
   and (do_ref) the reference household once and E at every percentile.
   Returns (vec theta', vec E'): (g-1)k + j.                               */
real rowvector _gepwreg_qr_prof(string scalar depvar, string scalar xvars, string scalar fw_var,
                                string scalar touse, string scalar grid_name, string scalar G_name,
                                string scalar ug_name, real scalar cband, real scalar band_in,
                                real scalar do_optbw, real scalar do_ref, | real matrix U)
{
    real colvector y, fw0, fw, c, ext, u, pc, w, ww, rk, isfac, uval
    real matrix    X, G, Q, Bi, Pb, Pr, RX, wp
    real rowvector ug, tg, umid, xfix, xr, zr, bs
    real scalar    n, k, ng, T, g, hsil, hg, q, ustar, seg, tt, clip
    st_view(y,   ., depvar,        touse)
    st_view(X,   ., tokens(xvars), touse)
    st_view(fw0, ., fw_var,        touse)
    n  = rows(y)
    k  = cols(X)
    fw = fw0 :/ mean(fw0)
    G  = st_matrix(G_name)
    ug = st_matrix(ug_name)
    ng = cols(ug)
    tg = st_matrix(grid_name)
    T  = cols(tg)
    Q    = (X, J(n, 1, 1)) * G'
    c    = rowsum(Q :< y)
    umid = (ug[1..(ng-1)] + ug[2..ng]) :/ 2
    ext  = (ug[1], umid, ug[ng])'
    u    = ext[c :+ 1]
    Bi   = _gepwreg_interp(G[., 1..k], u, ug)
    wp   = _gepwreg_wp(y, fw, tg[1], cband, band_in, J(0,1,.))
    pc   = wp[., 2]
    hsil = st_numscalar("_gepwreg_h_tmp")
    Pb   = J(T, k, .)
    for (g = 1; g <= T; g++) {
        if (do_optbw & band_in <= 0) hg = _gepwreg_hopt2(Bi, fw, pc, tg[g], hsil, n)
        else                         hg = hsil
        w  = _gepwreg_kw(pc, tg[g], hg, n)
        ww = fw :* w
        Pb[g, .] = (ww' * Bi) / sum(ww)
    }
    if (!do_ref) return(vec(Pb')')
    rk    = st_matrix("_gepwreg_rk")
    isfac = st_matrix("_gepwreg_rfac")
    RX    = st_matrix("_gepwreg_rx")
    xfix  = st_matrix("_gepwreg_xfix")
    xr    = J(1, 0, .)
    zr    = J(1, 0, .)
    uval  = J(0, 1, .)
    _gepwreg_refrow(X, J(n, 0, .), fw, rk, isfac, RX, J(rows(rk), 0, .),
                    xfix, J(1, 0, .), xr, zr, uval)
    Pr = J(T, k, .)
    /* U (the point only, not the draws): u* and the clip, +1 when q_tau is
       above the reference's conditional quantiles of the grid (the curve
       held at the last quantile), -1 below (held at the first), 0 inside   */
    if (args() == 12) U = J(T, 2, .)
    for (g = 1; g <= T; g++) {
        q     = _gepwreg_qtau(y, fw, tg[g])
        ustar = .
        seg   = .
        tt    = .
        clip  = .
        bs    = _gepwreg_E_pt(G, ug, q, (xr, 1), ustar, seg, tt, clip)
        Pr[g, .] = bs[1..k]
        if (args() == 12) U[g, .] = (ustar, clip * (tt == 1 ? 1 : -1))
    }
    return((vec(Pb')', vec(Pr')'))
}

/* het(qr): the profile at the point, stored */
void _gepwreg_qr_prof_pt(string scalar depvar, string scalar xvars, string scalar fw_var,
                         string scalar touse, string scalar grid_name, string scalar G_name,
                         string scalar ug_name, real scalar cband, real scalar band_in,
                         real scalar do_optbw, real scalar do_ref)
{
    real rowvector v
    real matrix    U
    real scalar    T, k
    v = _gepwreg_qr_prof(depvar, xvars, fw_var, touse, grid_name, G_name, ug_name,
                         cband, band_in, do_optbw, do_ref, U)
    T = cols(st_matrix(grid_name))
    k = cols(tokens(xvars))
    st_matrix("_gepwreg_pb", rowshape(v[1..(T*k)], T))
    if (do_ref) {
        st_matrix("_gepwreg_prb", rowshape(v[(T*k + 1)..(2*T*k)], T))
        st_matrix("_gepwreg_pru", U)
    }
}

/* het(qr): after the draws -- the standard errors, the simultaneous bands
   and the flatness tests of the profile (and of the reference curve), and
   the bootstrap variances at per(), the last percentile of the grid      */
void _gepwreg_qr_prof_boot(string scalar TP_name, real scalar k, real scalar do_ref,
                           real scalar level, real scalar nuse)
{
    real matrix    TP, Pb, Pr
    real scalar    T
    real rowvector last
    TP = st_matrix(TP_name)
    TP = select(TP, rowmissing(TP) :== 0)
    Pb = st_matrix("_gepwreg_pb")
    T  = rows(Pb)
    _gepwreg_prof_bstats(Pb, TP[., 1..(T*k)], level, nuse, "_gepwreg_p")
    last = (T - 1) * k :+ (1..k)
    st_matrix("_gepwreg_Vb", _gepwreg_bootvar(TP[., last]))
    if (do_ref) {
        Pr = st_matrix("_gepwreg_prb")
        _gepwreg_prof_bstats(Pr, TP[., (T*k + 1)..(2*T*k)], level, nuse, "_gepwreg_pr")
        st_matrix("_gepwreg_Vref_boot", _gepwreg_bootvar(TP[., T*k :+ last]))
    }
}


/* ---------------------------------------------------------------------------
   Measurement-error correction for het(z) -- see notes section 15.
   W = (X, Z, X x Z, 1); one sigma^2 per heterogeneity variable, errors
   independent across them.  fw is normalised to mean 1 upstream.
   ------------------------------------------------------------------------- */
real scalar _gepwreg_gmw_ix(real scalar k, real scalar m, real scalar j, real scalar l)
{
    return(k + m + (j - 1) * m + l)
}

/* (1.4) nodots.  The progress dots cost about 26 microseconds a
   replication, which is nothing, but a user who wants to measure that for
   himself should be able to switch them off -- and so should anyone piping
   the output somewhere that does not want them.                          */
real scalar _gepwreg_gmw_dots_on()
{
    real scalar v
    v = st_numscalar("_gepwreg_dots")
    if (v == .) return(1)
    return(v != 0)
}

/* weighted column means of A */
real rowvector _gepwreg_gmw_mean(real matrix A, real colvector fw)
{
    return((fw' * A) / sum(fw))
}

/* (1.4) The two a priori conditions of the correction, one row per
   heterogeneity variable, so that the user can read WHY a variable was left
   alone instead of being told only that its t was small.

   col 1  implied reliability, 1 - sigma2 / Var(z~).  More interpretable than
          sigma2 itself, which has no scale: it says how much work the
          correction is doing, hence how much variance it costs.
   col 2  skewness of z~.  Classical error DILUTES skewness and never creates
          it -- with u symmetric, skew(z~) = skew(z) * reliability^(3/2) -- so
          a visible skewness is a real one and a lower bound on the latent.
          With everything Gaussian the model is not identified from moments at
          all, whatever the sample size.
   col 3  R2 of z~ on the regressors.  The identifying moments have -2 sigma2
          E[q z] on their right-hand side, q = c + sum_j d_j x_j the
          derivative of the step-1 index in z: for a variable uncorrelated
          with the regressors and centred, E[q z] = 0, and it carries no
          information about its own error, however large that error is.

   No threshold is applied to either.  The diagnosis is shown; the judgement
   stays with the reader.                                                   */
real matrix _gepwreg_gmw_diag(real matrix X, real matrix Zt, real colvector fw,
                      real colvector s2)
{
    real scalar n, k, m, l, sw, mu, v, m3, rss
    real matrix D, Xc, XtX
    real colvector cf, r, zc
    n = rows(Zt) ; k = cols(X) ; m = cols(Zt) ; sw = sum(fw)
    Xc = (X, J(n, 1, 1))
    XtX = cross(Xc, fw, Xc)
    D = J(m, 3, .)
    for (l = 1; l <= m; l++) {
        mu = (fw' * Zt[., l]) / sw
        zc = Zt[., l] :- mu
        v  = (fw' * (zc:^2)) / sw
        m3 = (fw' * (zc:^3)) / sw
        cf = lusolve(XtX, cross(Xc, fw, Zt[., l]))
        r  = Zt[., l] - Xc * cf
        rss = (fw' * (r:^2)) / sw
        if (v > 0) {
            if (s2[l] > 0) D[l, 1] = 1 - s2[l] / v
            D[l, 2] = m3 / v^1.5
            D[l, 3] = 1 - rss / v
        }
    }
    return(D)
}

real matrix _gepwreg_gmw_Omega(real matrix X, real colvector fw, real colvector s2,
                       real scalar k, real scalar m)
{
    real scalar p, l, j, jj, zl, a, b
    real matrix O, Exx
    real rowvector Ex
    p   = k + m + k * m + 1
    O   = J(p, p, 0)
    Ex  = _gepwreg_gmw_mean(X, fw)
    Exx = cross(X, fw, X) / sum(fw)
    for (l = 1; l <= m; l++) {
        zl = k + l
        O[zl, zl] = O[zl, zl] + s2[l]
        for (j = 1; j <= k; j++) {
            a = _gepwreg_gmw_ix(k, m, j, l)
            O[zl, a] = O[zl, a] + Ex[j] * s2[l]
            O[a, zl] = O[zl, a]
            for (jj = 1; jj <= k; jj++) {
                b = _gepwreg_gmw_ix(k, m, jj, l)
                O[a, b] = O[a, b] + Exx[j, jj] * s2[l]
            }
        }
    }
    return(O)
}

real matrix _gepwreg_gmw_q(real matrix X, real colvector bb, real scalar k,
                   real scalar m)
{
    real scalar n, l, j
    real matrix Q
    n = rows(X)
    Q = J(n, m, 0)
    for (l = 1; l <= m; l++) {
        Q[., l] = J(n, 1, bb[k + l])
        for (j = 1; j <= k; j++) {
            Q[., l] = Q[., l] + bb[_gepwreg_gmw_ix(k, m, j, l)] :* X[., j]
        }
    }
    return(Q)
}

real matrix _gepwreg_gmw_m1(real matrix W, real matrix X, real colvector y,
                    real colvector bb, real colvector s2, real scalar k,
                    real scalar m)
{
    real scalar l, j, a
    real matrix out, Q
    real colvector e
    e   = y - W * bb
    Q   = _gepwreg_gmw_q(X, bb, k, m)
    out = W :* e
    for (l = 1; l <= m; l++) {
        out[., k + l] = out[., k + l] + s2[l] :* Q[., l]
        for (j = 1; j <= k; j++) {
            a = _gepwreg_gmw_ix(k, m, j, l)
            out[., a] = out[., a] + s2[l] :* X[., j] :* Q[., l]
        }
    }
    return(out)
}

real matrix _gepwreg_gmw_m2(real matrix W, real matrix X, real matrix Zt,
                    real colvector y, real colvector bb, real colvector s2,
                    real scalar k, real scalar m)
{
    real scalar n, l, j, c
    real matrix out, Q
    real colvector e, z2
    n   = rows(y)
    e   = y - W * bb
    Q   = _gepwreg_gmw_q(X, bb, k, m)
    out = J(n, m * (k + 1), 0)
    c   = 0
    for (l = 1; l <= m; l++) {
        z2 = Zt[., l]:^2
        c++
        out[., c] = z2 :* e + 2 * s2[l] :* Zt[., l] :* Q[., l]
        for (j = 1; j <= k; j++) {
            c++
            out[., c] = X[., j] :* z2 :* e + 2 * s2[l] :* X[., j] :* Zt[., l] :* Q[., l]
        }
    }
    return(out)
}

struct _gepwreg_gmw_moms {
    real matrix A
    real matrix C
    real matrix E1
    real matrix E2
}

struct _gepwreg_gmw_moms scalar _gepwreg_gmw_pre(real matrix X, real matrix Zt,
                                 real colvector y, real colvector fw,
                                 real matrix W, real scalar k, real scalar m)
{
    struct _gepwreg_gmw_moms scalar P
    real scalar n, l, jj, j2, r, sw
    real colvector xj, z2, wv
    n  = rows(y)
    sw = sum(fw)
    P.A  = J(m, k + 1, 0)
    P.C  = J(m * (k + 1), cols(W), 0)
    P.E1 = J(m, k + 1, 0)
    P.E2 = J(m * (k + 1), k, 0)
    for (l = 1; l <= m; l++) {
        z2 = Zt[., l]:^2
        for (jj = 0; jj <= k; jj++) {
            xj = (jj == 0 ? J(n, 1, 1) : X[., jj])
            wv = xj :* z2
            r  = (l - 1) * (k + 1) + jj + 1
            P.A[l, jj + 1]  = (fw' * (wv :* y)) / sw
            P.C[r, .]       = (fw' * (W :* wv)) / sw
            P.E1[l, jj + 1] = (fw' * (xj :* Zt[., l])) / sw
            for (j2 = 1; j2 <= k; j2++) {
                P.E2[r, j2] = (fw' * (xj :* X[., j2] :* Zt[., l])) / sw
            }
        }
    }
    return(P)
}

real rowvector _gepwreg_gmw_gfast(struct _gepwreg_gmw_moms scalar P, real colvector bb,
                          real colvector s2, real scalar k, real scalar m)
{
    real scalar l, jj, j2, r, acc
    real rowvector g
    g = J(1, m * (k + 1), 0)
    for (l = 1; l <= m; l++) {
        for (jj = 0; jj <= k; jj++) {
            r   = (l - 1) * (k + 1) + jj + 1
            acc = bb[k + l] * P.E1[l, jj + 1]
            for (j2 = 1; j2 <= k; j2++) {
                acc = acc + bb[_gepwreg_gmw_ix(k, m, j2, l)] * P.E2[r, j2]
            }
            g[r] = P.A[l, jj + 1] - (P.C[r, .] * bb)[1,1] + 2 * s2[l] * acc
        }
    }
    return(g)
}

real scalar _gepwreg_gmw_critf(struct _gepwreg_gmw_moms scalar P, real matrix X,
                       real colvector fw, real matrix M, real colvector my,
                       real matrix A, real colvector s2, real scalar k,
                       real scalar m)
{
    real rowvector g
    real colvector bb
    real matrix Om
    Om = _gepwreg_gmw_Omega(X, fw, s2, k, m)
    if (min(symeigenvalues(M - Om)) <= 1e-12) return(1e12)
    bb = lusolve(M - Om, my)
    g  = _gepwreg_gmw_gfast(P, bb, s2, k, m)
    return((g * A * g')[1,1])
}

real colvector _gepwreg_gmw_hi(real matrix X, real matrix Zt, real colvector fw,
                       real scalar k, real scalar m)
{
    real scalar n, l, sw
    real matrix Xc, XtX
    real colvector hi, cf, r
    n   = rows(X)
    sw  = sum(fw)
    Xc  = (X, J(n, 1, 1))
    XtX = cross(Xc, fw, Xc)
    hi  = J(m, 1, .)
    for (l = 1; l <= m; l++) {
        cf    = lusolve(XtX, cross(Xc, fw, Zt[., l]))
        r     = Zt[., l] - Xc * cf
        hi[l] = 0.97 * (fw' * (r:^2)) / sw
    }
    return(hi)
}

/* weighted covariance of the moment columns */
real matrix _gepwreg_gmw_vcov(real matrix mm, real colvector fw)
{
    real matrix mc
    mc = mm :- _gepwreg_gmw_mean(mm, fw)
    return(cross(mc, fw, mc) / sum(fw))
}

real matrix _gepwreg_gmw_u(real matrix X, real scalar k, real scalar m, real scalar l)
{
    real scalar n, p, j
    real matrix u
    n = rows(X)
    p = k + m + k * m + 1
    u = J(n, p, 0)
    u[., k + l] = J(n, 1, 1)
    for (j = 1; j <= k; j++) u[., _gepwreg_gmw_ix(k, m, j, l)] = X[., j]
    return(u)
}

void _gepwreg_gmw_inf(real matrix X, real matrix Zt, real colvector y,
              real colvector fw, real colvector s2, real colvector bb,
              real matrix psib, real matrix psis)
{
    real scalar n, k, m, p, q, l, j, r, sw
    real matrix W, M, Q, m1, m2, G1b, G1s, G2b, G2s, G1bi, mtil, dbds, D, A, U
    real colvector xj, z2, sc
    n = rows(y) ; k = cols(X) ; m = cols(Zt)
    sw = sum(fw)
    p = k + m + k * m + 1
    q = m * (k + 1)
    W  = _gepwreg_Wmat(X, Zt)
    M  = cross(W, fw, W) / sw
    Q  = _gepwreg_gmw_q(X, bb, k, m)
    m1 = _gepwreg_gmw_m1(W, X, y, bb, s2, k, m)
    m2 = _gepwreg_gmw_m2(W, X, Zt, y, bb, s2, k, m)
    G1b = -M + _gepwreg_gmw_Omega(X, fw, s2, k, m)
    G1s = J(p, m, 0)
    G2b = J(q, p, 0)
    G2s = J(q, m, 0)
    r = 0
    for (l = 1; l <= m; l++) {
        U = _gepwreg_gmw_u(X, k, m, l)
        G1s[., l] = _gepwreg_gmw_mean(U :* Q[., l], fw)'
        z2 = Zt[., l]:^2
        for (j = 0; j <= k; j++) {
            xj = (j == 0 ? J(n, 1, 1) : X[., j])
            r++
            G2b[r, .] = -_gepwreg_gmw_mean(W :* (xj :* z2), fw) + 2 * s2[l] :* _gepwreg_gmw_mean(U :* (xj :* Zt[., l]), fw)
            G2s[r, l] = 2 * (fw' * (xj :* Zt[., l] :* Q[., l])) / sw
        }
    }
    G1bi = luinv(G1b)
    mtil = m2 - m1 * (G1bi' * G2b')
    dbds = -(G1bi * G1s)
    D    = G2b * dbds + G2s
    A    = invsym(_gepwreg_gmw_vcov(m2, fw))
    psis = -(mtil * (A * D) * luinv(D' * A * D))
    psib = -(m1 * G1bi') + psis * dbds'
    sc   = (n / sw) :* fw
    psib = psib :* sc
    psis = psis :* sc
}

real colvector _gepwreg_gmw_search(struct _gepwreg_gmw_moms scalar P, real matrix X,
                           real colvector fw, real matrix M, real colvector my,
                           real matrix A, real colvector s2, real colvector free,
                           real colvector hi, real scalar k, real scalar m,
                           real scalar nsweep)
{
    real scalar sweep, l, gr, a, b, c, d, fc, fd, it
    real colvector t
    gr = (sqrt(5) - 1) / 2
    for (sweep = 1; sweep <= nsweep; sweep++) {
        for (l = 1; l <= m; l++) {
            if (free[l] == 0) continue
            a = 0
            b = hi[l]
            c = b - gr * (b - a)
            d = a + gr * (b - a)
            t = s2 ; t[l] = c
            fc = _gepwreg_gmw_critf(P, X, fw, M, my, A, t, k, m)
            t[l] = d
            fd = _gepwreg_gmw_critf(P, X, fw, M, my, A, t, k, m)
            for (it = 1; it <= 70; it++) {
                if (fc < fd) {
                    b = d ; d = c ; fd = fc
                    c = b - gr * (b - a)
                    t = s2 ; t[l] = c
                    fc = _gepwreg_gmw_critf(P, X, fw, M, my, A, t, k, m)
                }
                else {
                    a = c ; c = d ; fc = fd
                    d = a + gr * (b - a)
                    t = s2 ; t[l] = d
                    fd = _gepwreg_gmw_critf(P, X, fw, M, my, A, t, k, m)
                }
            }
            s2[l] = (a + b) / 2
        }
    }
    return(s2)
}

real colvector _gepwreg_gmw_fit_t(real matrix X, real matrix Zt, real colvector y,
                          real colvector fw, real scalar tcrit,
                          real colvector bb, real colvector tstat,
                          real colvector keep, real matrix psib,
                          real matrix psis)
{
    real scalar n, k, m, rep, l, sw
    real matrix W, M, A, Mc
    real colvector my, s2, hi, free, sc, se, ev
    struct _gepwreg_gmw_moms scalar P
    n = rows(y) ; k = cols(X) ; m = cols(Zt) ; sw = sum(fw)
    W  = _gepwreg_Wmat(X, Zt)
    M  = cross(W, fw, W) / sw
    my = cross(W, fw, y) / sw
    P  = _gepwreg_gmw_pre(X, Zt, y, fw, W, k, m)
    hi = _gepwreg_gmw_hi(X, Zt, fw, k, m)
    /* (1.4) A singular step-1 design is not an over-correction, and the
       two must not share a message.  lusolve() on a singular matrix returns
       missing values, which used to travel to ereturn post and surface as
       "estimates post: matrix has missing values", r(504).  The common
       cause is a factor expansion of a variable that is also a
       heterogeneity variable: with i.gse among the regressors and gse in
       het(), the interaction 2.gse x gse equals 2 x 2.gse exactly.  The
       existing same[] guard compares names, so "gse" and "2.gse" slip past
       it.                                                              */
    ev = symeigenvalues(makesymmetric(M))
    if (min(ev) <= 1e-12 * max(ev)) {
        errprintf("merr: the step-1 design is singular -- its columns are collinear.\n")
        errprintf("      A frequent cause is a factor expansion of a variable that is\n")
        errprintf("      also a heterogeneity variable: with i.z among the regressors\n")
        errprintf("      and z in het(), the interaction of an indicator of z with z\n")
        errprintf("      is proportional to the indicator itself.  Use i.z in het() or\n")
        errprintf("      z among the regressors, not both.\n")
        exit(198)
    }
    s2 = J(m, 1, 0)
    free = J(m, 1, 1)
    sc = abs(_gepwreg_gmw_gfast(P, lusolve(M, my), s2, k, m))' :+ 1e-10
    A  = diag(1 :/ sc:^2)
    for (rep = 1; rep <= 2; rep++) {
        s2 = _gepwreg_gmw_search(P, X, fw, M, my, A, s2, free, hi, k, m, 8)
        bb = lusolve(M - _gepwreg_gmw_Omega(X, fw, s2, k, m), my)
        A  = invsym(_gepwreg_gmw_vcov(_gepwreg_gmw_m2(W, X, Zt, y, bb, s2, k, m), fw))
    }
    Mc = M - _gepwreg_gmw_Omega(X, fw, s2, k, m)
    /* (1.4) The correction subtracts Omega from the moment matrix.  Too
       large a sigma^2 takes M - Omega out of positive definiteness, the
       solve returns missing values, and they used to travel all the way to
       ereturn post, where Stata says only "matrix has missing values".
       The notes of 18sep2026 asked for this guard; here it is.         */
    if (min(symeigenvalues(makesymmetric(Mc))) <= 0 | hasmissing(Mc)) {
        errprintf("merr: the corrected moment matrix is not positive definite.\n")
        errprintf("      The estimated error variances are too large for this design;\n")
        errprintf("      the correction is refused rather than returned as missing.\n")
        exit(198)
    }
    bb = lusolve(Mc, my)
    if (hasmissing(bb)) {
        errprintf("merr: the corrected step-1 solve returned missing values.\n")
        exit(198)
    }
    _gepwreg_gmw_inf(X, Zt, y, fw, s2, bb, psib, psis)
    se = J(m, 1, .)
    tstat = J(m, 1, 0)
    for (l = 1; l <= m; l++) {
        se[l] = sqrt(psis[., l]' * psis[., l]) / n
        /* (1.4) se is the scale of the information the moments carry
           about sigma^2_l.  When it collapses relative to hi[l], the
           widest admissible variance, they carry none and s2/se is a
           ratio of two numerical residuals.  Declare the variable not
           identified instead of reporting a t that means nothing.    */
        if (se[l] > 1e-8 * hi[l]) tstat[l] = s2[l] / se[l]
        else {
            tstat[l] = 0
            s2[l]    = 0
        }
    }
    keep = (tstat :> tcrit)
    if (sum(keep) < m) {
        for (l = 1; l <= m; l++) {
            if (keep[l] == 0) s2[l] = 0
        }
        free = keep
        if (sum(keep) > 0) {
            s2 = _gepwreg_gmw_search(P, X, fw, M, my, A, s2, free, hi, k, m, 8)
        }
        bb = lusolve(M - _gepwreg_gmw_Omega(X, fw, s2, k, m), my)
        _gepwreg_gmw_inf(X, Zt, y, fw, s2, bb, psib, psis)
    }
    return(s2)
}

/* theta_j(tau), its i.i.d. standard error, and the per-observation scores.
   The rank is the WEIGHTED one: a household's place in the population.      */
/* Design-based variance from arbitrary per-observation scores.
   meat = sum_h (n_h/(n_h-1)) sum_p (Z_ph - Zbar_h)'(Z_ph - Zbar_h), / n^2.
   psi already carries the bread, so no sandwich is applied here.            */

/* ─────────────────────────────────────────────────────────────────────────────
   _gepwreg_ts_z()  two-step estimator, heterogeneity through observables z

   Step 1: y = W gamma + e, W = (X, Z, X x Z, 1), weighted by fw.
           b_ij = gamma_j + sum_l gamma_jl z_il.
   Step 2: theta_j = sum fw w b_ij / sum fw w  (w: kernel on the rank of y)
                   = gamma_j + sum_l gamma_jl zbar_l(tau)      (closed form)
   Variance (analytical): influence function of theta_j,
       psi_j(theta) = psi(gamma_j) + sum_l [ psi(gamma_jl) zbar_l
                                            + gamma_jl psi(zbar_l) ],
     with psi(gamma) = (W'FW/n)^-1 W_i fw_i e_i the OLS influence, and
     psi(zbar_l) = direct term  n fw_i w_i (z_il - zbar_l) / S
                 + rank term    fw_i [ sum_{i': y_i' >= y_i} a_i' - sum_i' a_i' pc_i' ] / S,
     a_i = dw_i fw_i (z_il - zbar_l),  dw_i = -0.5 (pc_i - tau) w_i / h^2,
     S = sum fw w.  V = psi'psi / n^2 (scores centred).  Under a survey
     design the Taylor formula is applied to the same scores.
   Bootstrap: pairs, or PSUs within strata under a design, the whole
     procedure re-run at each draw.
   Initial-quantile object: y0_i = y_i - sum_j b_ij (x_ij - xref_j), then
     the same kernel mean on the rank of y0 (bootstrap variance).
   ───────────────────────────────────────────────────────────────────────── */
void _gepwreg_ts_z(string scalar depvar, string scalar xvars, string scalar zvars,
                   string scalar fw_var, string scalar touse,
                   real scalar tau, real scalar cband, real scalar band_in,
                   real scalar B, real scalar do_svy,
                   string scalar svy_psu, string scalar svy_strata,
                   real scalar do_initial, string scalar xref_name, real scalar seed,
                   string scalar same_name, string scalar iflag_name,
                   real scalar do_optbw, real scalar do_merr, real scalar tcrit,
                   real scalar do_ref, real scalar do_prof, real scalar nuse,
                   real scalar level)
{
    real rowvector tg, xg, zg
    real matrix    Pb, PSIP, VP
    real colvector wg, wwg, ordp
    real scalar    T, g, hg, Sg, hsil
    real colvector rk, isfac, uval, uvb
    real matrix    RX, RZ, psir, Vr, Vrs, TTr, dumP
    real rowvector xfix, zfix, xrr, xrb, thr
    real matrix    same, psi_xb, psib_m, psis_m, pbb_m, pss_m
    real colvector s2_m, bb_m, t_m, keep_m
    real colvector s2b_m, bbb_m, tb_m, kb_m
    real rowvector xbar, iflag, mu_d, xbar_d, isbin
    real scalar    jj, i_dl
    real colvector y, fw0, fw, w, pc, e, ww, dw, a, rc, ord, y0, w0, ww0, idx
    real colvector yb, fwb, eb, wb, wwb, y0b, w0b, ww0b, psu_vec, strata_vec
    real matrix    X, Z, W, A, Ainv, Bi, wp, psi_g, psi_mu, psi_th, V, Vs, Vb
    real matrix    TT, decomp, wp0, Xb, Zb, Wb, Bib, wpb, wp0b
    real rowvector gam, theta, mu, zpop, xref, theta0, gamb, thb
    real scalar    n, k, m, p, h, S, neff, j, l, idx_jl, h0, neff0, b, clus, Sb
    real scalar    hb

    st_view(y,   ., depvar,        touse)
    st_view(X,   ., tokens(xvars), touse)
    st_view(Z,   ., tokens(zvars), touse)
    st_view(fw0, ., fw_var,        touse)
    n  = rows(y)
    k  = cols(X)
    m  = cols(Z)
    fw = fw0 :/ mean(fw0)
    same = st_matrix(same_name)
    /* (1.5.0, review C8) the regressors taking only the values 0 and 1 */
    isbin = J(1, k, 0)
    for (jj = 1; jj <= k; jj++) isbin[jj] = all((X[., jj] :== 0) :| (X[., jj] :== 1))

    /* ── step 1 ──────────────────────────────────────────────────────── */
    W    = _gepwreg_Wmat(X, Z)
    p    = cols(W)
    A    = quadcross(W, fw, W) / n
    Ainv = invsym(A)
    gam  = (Ainv * (quadcross(W, fw, y) / n))'
    if (do_merr) {
        /* A heterogeneity variable that is ALSO a regressor puts the
           measurement error in the X block too, which Omega does not cover:
           that is the x = z degeneracy where the estimand is not defined.
           Refuse rather than return numbers that look right.               */
        if (sum(same) > 0) {
            errprintf("merr: a heterogeneity variable is also a regressor;")
            errprintf(" the correction does not cover that case.\n")
            exit(198)
        }
        bb_m = J(0, 1, .) ; t_m = J(0, 1, .) ; keep_m = J(0, 1, .)
        psib_m = J(0, 0, .) ; psis_m = J(0, 0, .)
        s2_m = _gepwreg_gmw_fit_t(X, Z, y, fw, tcrit, bb_m, t_m, keep_m, psib_m, psis_m)
        gam  = bb_m'
        st_matrix("_gepwreg_s2",   s2_m')
        st_matrix("_gepwreg_tm",   t_m')
        st_matrix("_gepwreg_keep", keep_m')
        st_matrix("_gepwreg_mdiag", _gepwreg_gmw_diag(X, Z, fw, s2_m))
    }
    e    = y - W * gam'
    Bi   = _gepwreg_betai(X, Z, gam, k, same, isbin)

    /* ── step 2 ──────────────────────────────────────────────────────── */
    wp    = _gepwreg_wp(y, fw, tau, cband, band_in, J(0,1,.))
    pc    = wp[., 2]
    h     = st_numscalar("_gepwreg_h_tmp")
    hsil  = h
    if (do_optbw & band_in <= 0) {
        h = _gepwreg_hopt2(Bi, fw, pc, tau, h, n)
        st_numscalar("_gepwreg_h_tmp", h)
        w = _gepwreg_kw(pc, tau, h, n)
    }
    else w = wp[., 1]
    ww    = fw :* w
    S     = sum(ww)
    theta = (ww' * Bi) / S
    mu    = (ww' * Z)  / S
    xbar  = (ww' * X)  / S
    zpop  = (fw' * Z)  / n
    neff  = S^2 / sum(ww:^2)

    /* ── influence function ──────────────────────────────────────────── */
    if (do_merr) psi_g = psib_m
    else         psi_g = (W :* (fw :* e)) * Ainv
    dw     = (-0.5 / h^2) :* (pc :- tau) :* w
    ord    = order(y, 1)
    psi_mu = J(n, m, 0)
    l = 1
    while (l <= m) {
        a       = dw :* fw :* (Z[., l] :- mu[l])
        rc      = J(n, 1, 0)
        rc[ord] = _gepwreg_revCumSumT(a[ord], pc[ord])
        psi_mu[., l] = (n :* fw :* w :* (Z[., l] :- mu[l])) :/ S +
                       (fw :* (rc :- sum(a :* pc))) :/ S
        l++
    }
    /* influence of the weighted means of the regressors in the group, needed
       when a regressor is also a heterogeneity variable                     */
    psi_xb = J(n, k, 0)
    jj = 1
    while (jj <= k) {
        a       = dw :* fw :* (X[., jj] :- xbar[jj])
        rc      = J(n, 1, 0)
        rc[ord] = _gepwreg_revCumSumT(a[ord], pc[ord])
        psi_xb[., jj] = (n :* fw :* w :* (X[., jj] :- xbar[jj])) :/ S +
                        (fw :* (rc :- sum(a :* pc))) :/ S
        jj++
    }
    psi_th = J(n, k, 0)
    j = 1
    while (j <= k) {
        psi_th[., j] = psi_g[., j]
        l = 1
        while (l <= m) {
            idx_jl = k + m + (j-1)*m + l
            psi_th[., j] = psi_th[., j] + psi_g[., idx_jl] :* mu[l] + gam[idx_jl] :* psi_mu[., l]
            if (same[j, l]) {
                psi_th[., j] = psi_th[., j] + psi_g[., k + l]
                jj = 1
                while (jj <= k) {
                    i_dl = k + m + (jj-1)*m + l
                    psi_th[., j] = psi_th[., j] + psi_g[., i_dl] :* xbar[jj] + gam[i_dl] :* psi_xb[., jj]
                    jj++
                }
                /* (1.5.0, review C8) binary x_j: + d_jl (1 - 2 xbar_j(tau)) */
                if (isbin[j]) psi_th[., j] = psi_th[., j] + psi_g[., idx_jl] :* (1 - 2*xbar[j]) - (2*gam[idx_jl]) :* psi_xb[., j]
            }
            l++
        }
        j++
    }
    psi_th = psi_th :- mean(psi_th)
    V = quadcross(psi_th, psi_th) / n^2
    if (do_svy) Vs = _gepwreg_taylor_psi(psi_th, touse, svy_psu, svy_strata, n)

    /* ── (1.5.0) the effect of the reference household, ref() ────────── */
    if (do_ref) {
        rk    = st_matrix("_gepwreg_rk")
        isfac = st_matrix("_gepwreg_rfac")
        RX    = st_matrix("_gepwreg_rx")
        RZ    = st_matrix("_gepwreg_rz")
        xfix  = st_matrix("_gepwreg_xfix")
        zfix  = st_matrix("_gepwreg_zfix")
        psir  = .
        xrr   = .
        uval  = .
        thr   = _gepwreg_ref_z(X, Z, fw, gam, psi_g, same, isbin, rk, isfac, RX, RZ,
                               xfix, zfix, psir, xrr, uval, 1)
        psir  = psir :- mean(psir)
        Vr    = quadcross(psir, psir) / n^2
        if (do_svy) Vrs = _gepwreg_taylor_psi(psir, touse, svy_psu, svy_strata, n)
    }

    /* ── (1.5.0) the profile across the percentiles of graph(): theta at
          every percentile and its influence function, psi_g Mbar(tau)' +
          the IF of the kernel mean of b (Mbar(tau) = M at the weighted means,
          M being affine in x and z) -- the joint covariance gives the
          pointwise standard errors, the simultaneous band and the flatness
          test.  The reference household's effect does not depend on tau.  */
    if (do_prof) {
        tg   = st_matrix("_gepwreg_pgrid")
        T    = cols(tg)
        ordp = order(y, 1)
        Pb   = J(T, k, .)
        PSIP = J(n, T*k, .)
        for (g = 1; g <= T; g++) {
            if (do_optbw & band_in <= 0) hg = _gepwreg_hopt2(Bi, fw, pc, tg[g], hsil, n)
            else                         hg = hsil
            wg  = _gepwreg_kw(pc, tg[g], hg, n)
            wwg = fw :* wg
            Sg  = sum(wwg)
            xg  = (wwg' * X) / Sg
            zg  = (wwg' * Z) / Sg
            for (j = 1; j <= k; j++) {
                Pb[g, j] = (wwg' * Bi[., j]) / Sg
                PSIP[., (g-1)*k + j] = psi_g * _gepwreg_mrow(xg, zg, j, same, isbin[j])' +
                                       _gepwreg_ifmean(Bi[., j], fw, wg, pc, tg[g], hg, ordp)
            }
        }
        PSIP = PSIP :- mean(PSIP)
        if (do_svy) VP = _gepwreg_taylor_psi(PSIP, touse, svy_psu, svy_strata, n)
        else        VP = quadcross(PSIP, PSIP) / n^2
        _gepwreg_prof_stats(Pb, VP, level, seed, nuse, "_gepwreg_p")
        if (do_ref) {
            if (do_svy) st_matrix("_gepwreg_prse", J(T, 1, sqrt(diagonal(Vrs))'))
            else        st_matrix("_gepwreg_prse", J(T, 1, sqrt(diagonal(Vr))'))
            st_matrix("_gepwreg_prb",   J(T, 1, thr))
            st_matrix("_gepwreg_prcu",  J(1, k, invnormal((100 + level)/200)))
        }
    }

    /* ── initial-quantile object ─────────────────────────────────────── */
    if (do_initial) {
        xref   = st_matrix(xref_name)
        iflag  = st_matrix(iflag_name)
        y0     = y - rowsum(((X :- xref) :* Bi) :* iflag)
        wp0    = _gepwreg_wp(y0, fw, tau, cband, band_in, J(0,1,.))
        h0     = st_numscalar("_gepwreg_h_tmp")
        if (do_optbw & band_in <= 0) {
            h0 = _gepwreg_hopt2(Bi, fw, wp0[., 2], tau, h0, n)
            w0 = _gepwreg_kw(wp0[., 2], tau, h0, n)
        }
        else w0 = wp0[., 1]
        ww0    = fw :* w0
        theta0 = (ww0' * Bi) / sum(ww0)
        neff0  = sum(ww0)^2 / sum(ww0:^2)
        /* composition of the group at the tau-quantile of y0, for the
           decomposition of the reported effect theta0 */
        mu_d   = (ww0' * Z) / sum(ww0)
        xbar_d = (ww0' * X) / sum(ww0)
        st_numscalar("_gepwreg_h_tmp", h)
    }
    else {
        mu_d   = mu
        xbar_d = xbar
    }

    /* ── bootstrap ───────────────────────────────────────────────────── */
    if (B > 0) {
        clus = (do_svy & svy_psu != "")
        if (clus) {
            st_view(psu_vec, ., svy_psu, touse)
            if (svy_strata != "") st_view(strata_vec, ., svy_strata, touse)
            else                  strata_vec = J(n, 1, 1)
        }
        else {
            psu_vec    = J(0, 1, .)
            strata_vec = J(0, 1, .)
        }
        rseed(seed)
        if (_gepwreg_gmw_dots_on()) {
            stata("_dots 0, title(Bootstrap replications) reps(" +
                  strofreal(B) + ")")
            displayflush()
        }
        TT = J(B, k, .)
        if (do_ref) TTr = J(B, k, .)
        b  = 1
        while (b <= B) {
            idx  = _gepwreg_drawidx(n, clus, psu_vec, strata_vec)
            yb   = y[idx]
            fwb  = fw[idx]
            Xb   = X[idx, .]
            Zb   = Z[idx, .]
            Wb   = _gepwreg_Wmat(Xb, Zb)
            if (do_merr) {
                /* the resample must redo the WHOLE procedure, the correction
                   and its t rule included, or the bootstrap would describe a
                   different estimator from the point estimate               */
                bbb_m = J(0, 1, .) ; tb_m = J(0, 1, .) ; kb_m = J(0, 1, .)
                pbb_m = J(0, 0, .) ; pss_m = J(0, 0, .)
                s2b_m = _gepwreg_gmw_fit_t(Xb, Zb, yb, fwb, tcrit, bbb_m, tb_m, kb_m,
                                   pbb_m, pss_m)
                gamb  = bbb_m'
            }
            else gamb = (invsym(quadcross(Wb, fwb, Wb)) * quadcross(Wb, fwb, yb))'
            Bib  = _gepwreg_betai(Xb, Zb, gamb, k, same, isbin)
            if (do_initial) {
                y0b   = yb - rowsum(((Xb :- xref) :* Bib) :* iflag)
                wp0b  = _gepwreg_wp(y0b, fwb, tau, cband, band_in, J(0,1,.))
                if (do_optbw & band_in <= 0) {
                    hb   = _gepwreg_hopt2(Bib, fwb, wp0b[., 2], tau,
                                          st_numscalar("_gepwreg_h_tmp"), n)
                    ww0b = fwb :* _gepwreg_kw(wp0b[., 2], tau, hb, n)
                }
                else ww0b = fwb :* wp0b[., 1]
                thb   = (ww0b' * Bib) / sum(ww0b)
            }
            else {
                wpb   = _gepwreg_wp(yb, fwb, tau, cband, band_in, J(0,1,.))
                if (do_optbw & band_in <= 0) {
                    hb  = _gepwreg_hopt2(Bib, fwb, wpb[., 2], tau,
                                         st_numscalar("_gepwreg_h_tmp"), n)
                    wwb = fwb :* _gepwreg_kw(wpb[., 2], tau, hb, n)
                }
                else wwb = fwb :* wpb[., 1]
                thb   = (wwb' * Bib) / sum(wwb)
            }
            TT[b, .] = thb
            /* (1.5.0) ref(): the reference household of the draw */
            if (do_ref) {
                dumP = .
                xrb  = .
                uvb  = .
                TTr[b, .] = _gepwreg_ref_z(Xb, Zb, fwb, gamb, J(0, 0, .), same, isbin, rk, isfac,
                                           RX, RZ, xfix, zfix, dumP, xrb, uvb, 0)
            }
            if (_gepwreg_gmw_dots_on()) {
                stata("_dots " + strofreal(b) + " 0")
                displayflush()
            }
            b++
        }
        Vb = _gepwreg_bootvar(TT)
        st_matrix("_gepwreg_Vb", Vb)
        if (do_ref) st_matrix("_gepwreg_Vref_boot", _gepwreg_bootvar(TTr))
        st_numscalar("_gepwreg_h_tmp", h)
    }

    /* ── decomposition table: rows z_l; mean at tau, population mean,
          contribution gamma_jl * zbar_l to each theta_j ────────────── */
    decomp = J(m, 2 + k, .)
    l = 1
    while (l <= m) {
        decomp[l, 1] = mu_d[l]
        decomp[l, 2] = zpop[l]
        j = 1
        while (j <= k) {
            decomp[l, 2 + j] = gam[k + m + (j-1)*m + l] * mu_d[l]
            if (same[j, l]) {
                decomp[l, 2 + j] = decomp[l, 2 + j] + gam[k + l]
                for (jj = 1; jj <= k; jj++) decomp[l, 2 + j] = decomp[l, 2 + j] + gam[k + m + (jj-1)*m + l] * xbar_d[jj]
                /* (1.5.0, review C8) binary x_j */
                if (isbin[j]) decomp[l, 2 + j] = decomp[l, 2 + j] + gam[k + m + (j-1)*m + l] * (1 - 2*xbar_d[j])
            }
            j++
        }
        l++
    }

    /* ── store ───────────────────────────────────────────────────────── */
    st_matrix("_gepwreg_batq", theta)
    st_matrix("_gepwreg_Vi",   V)
    if (do_svy) st_matrix("_gepwreg_Vs", Vs)
    if (do_initial) {
        st_matrix("_gepwreg_b", theta0)
        st_matrix("_gepwreg_V", Vb)
        st_numscalar("_gepwreg_h0",    h0)
        st_numscalar("_gepwreg_neff0", neff0)
    }
    else {
        st_matrix("_gepwreg_b", theta)
        if (do_svy) st_matrix("_gepwreg_V", Vs)
        else        st_matrix("_gepwreg_V", V)
    }
    if (do_ref) {
        st_matrix("_gepwreg_bref",    thr)
        st_matrix("_gepwreg_refx",    xrr)
        st_matrix("_gepwreg_refu",    uval')
        st_matrix("_gepwreg_Vref_IF", Vr)
        if (do_svy) st_matrix("_gepwreg_Vref_svy", Vrs)
    }
    st_matrix("_gepwreg_g",      gam)
    st_matrix("_gepwreg_zbar",   mu_d)
    st_matrix("_gepwreg_zpop",   zpop)
    st_matrix("_gepwreg_decomp", decomp)
    st_numscalar("_gepwreg_h",    h)
    st_numscalar("_gepwreg_neff", neff)
    st_numscalar("_gepwreg_n",    n)
}


/* ─────────────────────────────────────────────────────────────────────────────
   _gepwreg_ts_qr()  two-step estimator, heterogeneity = conditional quantile
   regression (Firpo, Fortin and Lemieux 2009: the UQPE as the average of
   the conditional quantile partial effects at the units located at q_tau)

   G (grid x (k+1)) holds the qreg coefficients (regressors, constant) on
   the grid ug.  The conditional rank of unit i is the midpoint of the two
   grid quantiles between which y_i falls given x_i (clipped to the grid);
   its effect b_ij is the slope interpolated at that rank.  Step 2 as in
   _gepwreg_ts_z.  Variance by the bootstrap of the whole procedure, run
   from the Stata side.
   ───────────────────────────────────────────────────────────────────────── */
void _gepwreg_ts_qr(string scalar depvar, string scalar xvars, string scalar fw_var,
                    string scalar touse, real scalar tau, real scalar cband,
                    real scalar band_in, real scalar do_initial,
                    string scalar xref_name, string scalar G_name, string scalar ug_name,
                    string scalar iflag_name, real scalar do_optbw)
{
    real colvector y, fw0, fw, w, ww, c, u, ext, y0, w0, ww0
    real matrix    X, G, Q, Bi, wp, wp0
    real rowvector ug, theta, theta0, xref, umid, iflag
    real scalar    n, k, ng, h, S, neff, h0, neff0
    real colvector pc

    st_view(y,   ., depvar,        touse)
    st_view(X,   ., tokens(xvars), touse)
    st_view(fw0, ., fw_var,        touse)
    n  = rows(y)
    k  = cols(X)
    fw = fw0 :/ mean(fw0)
    G  = st_matrix(G_name)
    ug = st_matrix(ug_name)
    ng = cols(ug)

    /* conditional rank of each unit */
    Q    = (X, J(n, 1, 1)) * G'
    c    = rowsum(Q :< y)
    umid = (ug[1..(ng-1)] + ug[2..ng]) :/ 2
    ext  = (ug[1], umid, ug[ng])'
    u    = ext[c :+ 1]
    Bi   = _gepwreg_interp(G[., 1..k], u, ug)

    /* step 2 */
    wp    = _gepwreg_wp(y, fw, tau, cband, band_in, J(0,1,.))
    pc    = wp[., 2]
    h     = st_numscalar("_gepwreg_h_tmp")
    if (do_optbw & band_in <= 0) {
        h = _gepwreg_hopt2(Bi, fw, pc, tau, h, n)
        st_numscalar("_gepwreg_h_tmp", h)
        w = _gepwreg_kw(pc, tau, h, n)
    }
    else w = wp[., 1]
    ww    = fw :* w
    S     = sum(ww)
    theta = (ww' * Bi) / S
    neff  = S^2 / sum(ww:^2)

    if (do_initial) {
        xref   = st_matrix(xref_name)
        iflag  = st_matrix(iflag_name)
        y0     = y - rowsum(((X :- xref) :* Bi) :* iflag)
        wp0    = _gepwreg_wp(y0, fw, tau, cband, band_in, J(0,1,.))
        h0     = st_numscalar("_gepwreg_h_tmp")
        if (do_optbw & band_in <= 0) {
            h0 = _gepwreg_hopt2(Bi, fw, wp0[., 2], tau, h0, n)
            w0 = _gepwreg_kw(wp0[., 2], tau, h0, n)
        }
        else w0 = wp0[., 1]
        ww0    = fw :* w0
        theta0 = (ww0' * Bi) / sum(ww0)
        neff0  = sum(ww0)^2 / sum(ww0:^2)
        st_matrix("_gepwreg_b", theta0)
        st_numscalar("_gepwreg_h0",    h0)
        st_numscalar("_gepwreg_neff0", neff0)
    }
    else st_matrix("_gepwreg_b", theta)
    st_matrix("_gepwreg_batq", theta)
    st_numscalar("_gepwreg_h",    h)
    st_numscalar("_gepwreg_neff", neff)
    st_numscalar("_gepwreg_n",    n)
}

end
