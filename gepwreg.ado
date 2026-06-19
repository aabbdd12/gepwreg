*! gepwreg.ado  v1.2  Araar A.
*! Percentile Weights Regression
*! Araar (2016, 2023) ; Deville (1999) ; Newey & McFadden (1994)
/* ─────────────────────────────────────────────────────────────────────────────
   Syntax:
     gepwreg depvar [indepvars] [fw/aw/pw] [if] [in],
           [ PERcentile(#) CBAND(#) BAND(#) OPTbw BOOT(#) Level(#) noConstant ]
   ───────────────────────────────────────────────────────────────────────── */
#delimit ;

capture program drop gepwreg ;
program define gepwreg, eclass properties(svyb) sortpreserve ;
version 16.0 ;

syntax varlist(min=1 numeric fv) [fw aw pw] [if] [in] [,
    PERcentile(real 0.5)
    CBAND(real 0.9)
    BAND(real 0)
    OPTbw
    BOOT(integer 0)
    Level(cilevel)
    noConstant
    VCE(string)
    RANKvar(varname)
    /* Internal options passed by svy: prefix */
    SVY
    noHEader
    noTABle
] ;

/* ── 0. Mark sample ──────────────────────────────────────────────────────── */
marksample touse ;
qui count if `touse' ;
if r(N) == 0 error 2000 ;

/* ── 1. Parse varlist ────────────────────────────────────────────────────── */
local depvar    : word 1 of `varlist' ;
local indepvars : list varlist - depvar ;

/* Step 1 : fvexpand then filter base (Nb.) and omitted (No.) categories
   0b.rururb → removed    1.rururb → kept
   This gives the display names AND the clean list for fvrevar             */
if "`indepvars'" != "" {;
    fvexpand `indepvars' if `touse' ;
    local indepvars_exp_raw `r(varlist)' ;
    local indepvars_exp ;
    foreach v of local indepvars_exp_raw {;
        if !regexm("`v'","^[0-9]+b\.") & !regexm("`v'","^[0-9]+o\.") {;
            local indepvars_exp `indepvars_exp' `v' ;
        } ;
    } ;
} ;
else local indepvars_exp "" ;

/* Step 2 : fvrevar on the FILTERED list → real temp vars for Mata
   fvrevar on  hhsize 1.rururb  →  hhsize __000002  (no base var)
   Dimensions of indepvars_mata always equal dimensions of indepvars_exp  */
if "`indepvars_exp'" != "" {;
    fvrevar `indepvars_exp' if `touse' ;
    local indepvars_mata `r(varlist)' ;
} ;
else local indepvars_mata "" ;

/* ── 1b. Ranking variable ───────────────────────────────────────────────── */
/* Default: rank on depvar. If rankvar() specified, rank on that variable.  */
local rank_var "`depvar'" ;
local rank_label "y (outcome)" ;
if "`rankvar'" != "" {;
    local rank_var "`rankvar'" ;
    local rank_label "`rankvar'" ;
    markout `touse' `rankvar' ;
    di as text "Ranking variable : `rankvar' (instead of `depvar')" ;
} ;

/* ── 2. Survey weights ───────────────────────────────────────────────────── */
/* Priority rule (highest to lowest) :                                       */
/*   1. Weights declared in the command  [pw=weight]  → used as-is          */
/*   2. Weights declared in svyset       [pw=weight]  → used automatically   */
/*   3. No weights anywhere              → unweighted (fw=1)                 */
/* If both command and svyset declare weights, command weights take priority  */
/* and a note is displayed so the user is aware.                             */
tempvar fw_var ;
qui gen double `fw_var' = 1 if `touse' ;
local wgt_name "(unweighted)" ;
if "`exp'" != "" {;
    /* ── Priority 1: user declared weights in the command ── */
    local wexpr = subinstr("`exp'","=","",1) ;
    /* Trim spaces from wexpr for clean comparison */
    local wexpr = trim("`wexpr'") ;
    qui replace `fw_var' = `wexpr' if `touse' ;
    local wgt_type = upper(substr("`weight'",1,2)) ;
    local wgt_name "`wgt_type'[`wexpr']" ;
    /* Check svyset — compare with command weights */
    qui svyset ;
    if "`r(wvar)'" != "" {;
        local svywvar = trim("`r(wvar)'") ;
        if "`svywvar'" == "`wexpr'" {;
            /* Same variable declared in both places — silently use it once */
            /* No message needed : consistent declaration */
        } ;
        else {;
            /* Different variables — command takes priority, warn user */
            di as text "Note: command weights [pw=`wexpr'] take priority over" ;
            di as text "      svyset weights [`svywvar']. Using [pw=`wexpr']." ;
        } ;
    } ;
} ;
else {;
    /* ── Priority 2: no command weights — check svyset ── */
    qui svyset ;
    if "`r(wvar)'" != "" {;
        qui replace `fw_var' = `r(wvar)' if `touse' ;
        local wgt_name "PW[`r(wvar)'] (svyset)" ;
    } ;
    /* else: Priority 3 — unweighted, wgt_name stays "(unweighted)" */
} ;

/* ── 3. Validate options ─────────────────────────────────────────────────── */
if `percentile' <= 0 | `percentile' >= 1 {;
    di as error "percentile() must be strictly between 0 and 1" ;
    error 198 ;
} ;
if `boot' < 0 {;
    di as error "boot() must be a non-negative integer" ;
    error 198 ;
} ;
local do_optbw = ("`optbw'" != "") ;
if `do_optbw' & `band' != 0 {;
    di as error "optbw and band() are mutually exclusive" ;
    error 198 ;
} ;
local addcons = ("`constant'" == "") ;

/* ── 3b. Parse vce() option for survey design ───────────────────────────── */
/* Detect svy: prefix (Stata passes vce(linearized) or sets svy option)    */
local do_svy = 0 ;
local svy_psu    "" ;
local svy_strata "" ;
if "`svy'" != "" {;
    local vce "linearized" ;
} ;
if "`vce'" != "" {;
    if lower("`vce'")=="svy" | lower("`vce'")=="survey" | lower("`vce'")=="linearized" {;
        /* Check svyset has been called */
        qui svyset ;
        if "`r(wvar)'" == "" & "`r(su1)'" == "" {;
            di as error "vce(svy) requires svyset to be called first" ;
            error 119 ;
        } ;
        local do_svy    = 1 ;
        local svy_psu    "`r(su1)'" ;     /* svyset stores PSU in r(su1)    */
        local svy_strata "`r(strata1)'" ;
        /* If no weights supplied by user, take them from svyset            */
        if "`exp'" == "" & "`r(wvar)'" != "" {;
            qui replace `fw_var' = `r(wvar)' if `touse' ;
            local wgt_name "PW[`r(wvar)'] (svyset)" ;
        } ;
        di as text "Survey design detected:" ;
        if "`svy_psu'" == "" {;
            di as text "  PSU    : (not declared — variance may be conservative)" ;
        } ;
        else {;
            di as text "  PSU    : `svy_psu'" ;
        } ;
        di as text "  Strata : `svy_strata'" ;
        /* Check minimum PSU per stratum */
        if "`svy_strata'" != "" & "`svy_psu'" != "" {;
            qui tab `svy_strata' if `touse' ;
            qui {;
                tempvar _npsu ;
                bysort `svy_strata' (`svy_psu') : gen `_npsu' = (_n==1) ;
                bysort `svy_strata' : replace `_npsu' = sum(`_npsu') ;
                bysort `svy_strata' : replace `_npsu' = `_npsu'[_N] ;
            } ;
            qui sum `_npsu' if `touse' ;
            local min_psu = r(min) ;
            if `min_psu' < 5 {;
                di as text "Warning: minimum PSU/stratum = `min_psu'." ;
                di as text "  Taylor SE may be unreliable. Recommend n_h >= 10." ;
            } ;
        } ;
    } ;
    else {;
        di as error "vce(`vce') not supported. Use vce(svy) for survey design." ;
        error 198 ;
    } ;
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
if `addcons' local vnames `indepvars_exp' _cons ;
else         local vnames `indepvars_exp' ;

matrix colnames _gepwreg_b  = `vnames' ;
matrix rownames _gepwreg_V  = `vnames' ;
matrix colnames _gepwreg_V  = `vnames' ;
matrix rownames _gepwreg_Vn = `vnames' ;
matrix colnames _gepwreg_Vn = `vnames' ;

/* ── 7. Post results ─────────────────────────────────────────────────────── */
ereturn post _gepwreg_b _gepwreg_V, esample(`touse') obs(`n_obs') ;
ereturn matrix V_naive = _gepwreg_Vn ;
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
ereturn local  bw_method = cond(`do_optbw', "MSE-optimal", "Silverman") ;
ereturn scalar N_eff = `n_eff' ;
ereturn scalar boot   = `boot' ;
ereturn scalar do_svy = `do_svy' ;
ereturn local  cmd      "gepwreg" ;
ereturn local  rankvar  "`rank_var'" ;
ereturn local  ranklabel "`rank_label'" ;
ereturn local  wgt      "`wgt_name'" ;
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
di as text "Percentile Weights Regression" ;
di as text "{hline 60}" ;
di as text %28s "Target percentile (tau)" " = " as result %6.4f `tau' ;
if `do_optbw' {;
    di as text %28s "Bandwidth (h)"           " = " as result %9.6f `h' as text " (MSE-optimal)" ;
} ;
else {;
    di as text %28s "Bandwidth (h)"           " = " as result %9.6f `h' ;
} ;
di as text %28s "Weights"                 " = " as result "`wgt_name'" ;
di as text %28s "Ranking variable"        " = " as result "`rank_label'" ;
di as text %28s "Observations"            " = " as result %7.0f `n_obs' ;
di as text %28s "Kernel N_eff"            " = " as result %7.1f `n_eff' ;
if `boot' > 0 {;
    di as text %28s "Bootstrap replications" " = " as result %7.0f `boot' ;
} ;
di as text "{hline 60}" ;
if `do_svy' {;
    di as text "SE: Taylor linearisation under survey design (PSU + strata)" ;
} ;
else {;
    di as text "SE: IF-corrected analytical (Deville 1999)" ;
} ;
di as text "{hline 60}" ;
di "" ;
ereturn display, level(`level') ;
di "" ;
di as text "Naive WLS SE stored in e(V_naive)." ;
if `boot' > 0 di as text "Bootstrap SE stored in e(V_boot)." ;
di as text "Use {cmd:gepwreg_setable} to compare all SE types." ;

end ;


/* ─────────────────────────────────────────────────────────────────────────────
   MATA ENGINE
   #delimit cr required: Mata uses { } and ; internally — must not be
   intercepted by Stata's #delimit ; mode
   ───────────────────────────────────────────────────────────────────────── */
#delimit cr

/* Drop previous definitions so the ado can be sourced more than once */
capture mata: mata drop _gepwe_wp()
capture mata: mata drop _gepwreg_hopt()
capture mata: mata drop _wls()
capture mata: mata drop _revCumSum()
capture mata: mata drop _se_IF()
capture mata: mata drop _se_naive()
capture mata: mata drop _se_boot()
capture mata: mata drop _se_svy()
capture mata: mata drop _gepwreg_main()

mata:

/* ─────────────────────────────────────────────────────────────────────────────
   _gepwe_wp()
   Compute gepwe kernel weights and percentile ranks.
   Returns n x 2 matrix : col 1 = w_i,  col 2 = pc_i
   Bandwidth h is written to Stata scalar _gepwreg_h_tmp

   Kernel (exact gepwe.ado formula):
     w_i = exp(-0.25 * ((pc_i - tau)/h)^2) / (h * sqrt(2*pi) * n)

   Silverman bandwidth on the percentile scale:
     h = cband * min(sd(pc), IQR(pc)/1.34) * n^(-1/5)
   Override with band_in > 0
   ───────────────────────────────────────────────────────────────────────── */
real matrix _gepwe_wp(real colvector y,
                      real colvector fw,
                      real scalar    tau,
                      real scalar    cband,
                      real scalar    band_in,
                      real colvector zrank)
{
    /* zrank : alternative ranking variable z (rows=n) or empty (rows=0).
       If empty, rank on y (standard PWR).
       If provided, rank on z (generalised rankvar option).           */
    real scalar    n, h, tmp, q25, q75, sd_pc
    real colvector ord, fw_s, pc_s, pc_sorted, u, w_s, w, pc, z_use

    n         = rows(y)

    /* Use z if provided, else fall back to y */
    z_use     = (rows(zrank) == n) ? zrank : y

    ord       = order(z_use, 1)       /* sort on z (or y if default) */
    fw_s      = fw[ord]
    pc_s      = runningsum(fw_s) :/ sum(fw_s)

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
    wp_sil = _gepwe_wp(y, fw, tau, cband, 0, J(0,1,.))
    h_sil  = st_numscalar("_gepwreg_h_tmp")
    h_pilot = 2 * h_sil

    d = 0.02

    /* Step 1: local-constant fits (weighted mean of y) at tau-d, tau, tau+d */
    wp_lo  = _gepwe_wp(y, fw, tau - d, cband, h_pilot, J(0,1,.))
    w_lo   = wp_lo[., 1]
    mu_lo  = sum(w_lo :* y) / sum(w_lo)

    wp_mid = _gepwe_wp(y, fw, tau,     cband, h_pilot, J(0,1,.))
    w_mid  = wp_mid[., 1]
    mu_mid = sum(w_mid :* y) / sum(w_mid)

    wp_hi  = _gepwe_wp(y, fw, tau + d, cband, h_pilot, J(0,1,.))
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
   _wls() : beta = (X'WX)^{-1} X'Wy
   ───────────────────────────────────────────────────────────────────────── */
real colvector _wls(real matrix X, real colvector y, real colvector w)
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
   _revCumSum()
   Reverse cumulative row sum of matrix M :
     R[j,.] = M[j,.] + M[j+1,.] + ... + M[n,.]

   Algorithm :
     Forward cumsum column-by-column (quadrunningsum requires a vector).
     revCumSum[i,.] = totalSum - forwardCumSum[i-1,.]
     where forwardCumSum[0,.] = 0 (prepended row of zeros).
   Complexity : O(n*k), one while-loop over k columns
   ───────────────────────────────────────────────────────────────────────── */
real matrix _revCumSum(real matrix M)
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


/* ─────────────────────────────────────────────────────────────────────────────
   _se_IF()
   Influence-function variance (Deville 1999 / Newey-McFadden 1994)

   Total influence of observation j :
     psi_j = w_j * x_j * e_j
             + (1/n) * sum_{i: y_i >= y_j}  dw_i * x_i * e_i

   Kernel derivative :
     dw_i = -0.5 * (pc_i - tau) / h^2 * w_i

   Variance :
     V = A^{-1} * (psi'psi / n) * A^{-1} / n,   A = X'WX / n
   ───────────────────────────────────────────────────────────────────────── */
real matrix _se_IF(real matrix  X,
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
    revCS[ord,] = _revCumSum(dS[ord,])
    psi        = S_dir + revCS :/ n
    V          = Ainv * (psi' * psi / n) * Ainv / n
    /* Warn if any diagonal of A is near zero (sparse category at tau) */
    if (min(diagonal(A)) < 1e-12) {
        printf("{txt}Warning: near-singular X'WX at tau=%g -- SE unreliable for sparse categories.\n", tau)
    }
    return(V)
}


/* ─────────────────────────────────────────────────────────────────────────────
   _se_naive()
   Standard WLS variance assuming fixed weights.
   INCONSISTENT for PWR -- stored in e(V_naive) for reference only.
   ───────────────────────────────────────────────────────────────────────── */
real matrix _se_naive(real matrix  X,
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
   _se_boot()
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
real matrix _se_boot(real matrix   X,
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
            wp_b    = _gepwe_wp(y_b, fw_b, tau, cband, band_in, J(0,1,.))
            w_b     = wp_b[., 1]
            BB[b, ] = _wls(X_b, y_b, w_b)'
            b++
        }
    }
    else {
        /* ── Pairs bootstrap (no survey design declared) ─────────────── */
        b = 1
        while (b <= B) {
            idx     = ceil(n :* runiform(n, 1))
            y_b     = y[idx]
            fw_b    = fw[idx]
            wp_b    = _gepwe_wp(y_b, fw_b, tau, cband, band_in, J(0,1,.))
            w_b     = wp_b[., 1]
            BB[b, ] = _wls(X[idx, ], y_b, w_b)'
            b++
        }
    }
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

    wp    = _gepwe_wp(y, fw, tau, cband, band_use, zrank)
    w     = wp[., 1]
    pc    = wp[., 2]
    h     = st_numscalar("_gepwreg_h_tmp")
    n_eff = sum(w)^2 / sum(w:^2)
    beta  = _wls(X, y, w)

    /* z-based ordering for IF indirect term */
    z_ord_v = (rows(zrank) == n) ? order(zrank, 1) : order(y, 1)

    /* Standard IF variance */
    V_IF    = _se_IF(X, y, w, pc, h, tau, beta, z_ord_v)
    V_naive = _se_naive(X, y, w, beta)

    /* Survey Taylor variance — overrides V_IF as main SE */
    if (do_svy) {
        V_svy = _se_svy(X, y, w, pc, h, tau, beta,
                        touse, svy_psu, svy_strata, n, z_ord_v)
        st_matrix("_gepwreg_Vs", V_svy)
        V_IF = V_svy
    }

    /* Bootstrap variance */
    if (B > 0) {
        V_boot = _se_boot(X, y, fw, w, tau, cband, band_in, beta, B,
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
note: variable y may be used before set.ile Weights Regression
*! Araar (2016, 2023) ; Deville (1999) ; Newey & McFadden (1994)
/* ─────────────────────────────────────────────────────────────────────────────
   Syntax:
     gepwreg depvar [indepvars] [fw/aw/pw] [if] [in],
           [ PERcentile(#) CBAND(#) BAND(#) OPTbw BOOT(#) Level(#) noConstant ]
   ───────────────────────────────────────────────────────────────────────── */
#delimit ;

capture program drop gepwreg ;
#delimit ;
/* ─────────────────────────────────────────────────────────────────────────────
   gepwreg_setable : side-by-side SE comparison table
   ───────────────────────────────────────────────────────────────────────── */
capture program drop gepwreg_setable ;
program define gepwreg_setable ;
    if "`e(cmd)'" != "gepwreg" {;
        di as error "gepwreg_setable: last estimation must be gepwreg" ;
        error 301 ;
    } ;
    matrix b    = e(b) ;
    matrix V_IF = e(V) ;
    matrix V_n  = e(V_naive) ;
    local hasBoot = (e(boot) > 0) ;
    if `hasBoot' matrix V_bt = e(V_boot) ;
    di "" ;
    di as text "SE comparison  (tau=" as result %5.3f e(tau)
       as text "  h=" as result %8.6f e(h) as text ")" ;
    di as text "{hline 72}" ;
    if `hasBoot' {;
        di as text %18s "Variable"
           %11s "Coeff."
           %11s "SE_naive"
           %11s "SE_IF"
           %11s "SE_boot"
           %9s  "IF/Boot" ;
    } ;
    else {;
        di as text %18s "Variable"
           %11s "Coeff."
           %11s "SE_naive"
           %11s "SE_IF"
           %10s "IF/Naive" ;
    } ;
    di as text "{hline 72}" ;
    local names : colnames e(b) ;
    local j = 1 ;
    foreach nm of local names {;
        local bj = b[1,`j'] ;
        local sn = sqrt(V_n[`j',`j']) ;
        local si = sqrt(V_IF[`j',`j']) ;
        if `hasBoot' {;
            local sb = sqrt(V_bt[`j',`j']) ;
            local r  = `si' / `sb' ;
            di as text %18s abbrev("`nm'",18)
               as result %11.5f `bj' %11.5f `sn'
                         %11.5f `si' %11.5f `sb' %9.3f `r' ;
        } ;
        else {;
            local r  = `si' / `sn' ;
            di as text %18s abbrev("`nm'",18)
               as result %11.5f `bj' %11.5f `sn'
                         %11.5f `si' %10.3f `r' ;
        } ;
        local ++j ;
    } ;
    di as text "{hline 72}" ;
    di as text "SE_naive : WLS weights-fixed (inconsistent)" ;
    di as text "SE_IF    : influence-function corrected (Deville 1999)" ;
    if `hasBoot' {;
        di as text "SE_boot  : pairs bootstrap B=" as result e(boot)
           as text " (weights reconstructed each draw)" ;
    } ;
    if e(do_svy) {;
        di as text "SE_svy   : Taylor linearisation (PSU + strata)" ;
        di as text "           stored in e(V_svy)" ;
    } ;
    di "" ;
end ;
