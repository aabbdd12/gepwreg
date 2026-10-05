* p13a_E_draws.do -- called by p13_profile_delta.do (Paper, version 5,
* Section 7.5): the estimator E of the effect for a given household under
* het(qr), its delta-method standard error with the Powell kernel density
* (as qreg vce(robust, kernel)) on the sample, and B = 200 bootstrap draws in
* which the whole estimator is re-run (the 19 quantile regressions, q_tau,
* the reference type, u*, E), with no standard error in the draws.
* Argument: what [chunk]
*   dgp         the perfect case, beta_x1(u) = 0.4 + 2u, beta_x2(u) = 0.5 + 3u,
*               n = 2000, tau = 0.25: point, the 200 draws and the report;
*   bkf c       Burkina (bkf98I.dta), lexp on size male urban i.gse,
*               unweighted, tau = 0.25: chunk c of 3 of the 200 draws (seed
*               777 + c); chunk 1 also computes the point;
*   bkfreport   combines the Burkina chunks and reports.
* Reference household: the default real type (continuous regressors at their
* mean, binary ones at their most frequent value, a factor at its most
* frequent level), recomputed in every draw; in the DGP also the type (5, 1).
* It does not call gepwreg: the estimator is programmed here, with official
* quantile regressions.  The draws are saved in results/ for p13b.
version 16
clear all
args what chunk
set more off
set varabbrev off
local OUT "results"
capture mkdir "results"
local tau 0.25
local B   200
local nch 3
global TAU `tau'
numlist "0.05(0.05)0.95"
global UG `r(numlist)'

mata:
real matrix _qr_psi(real matrix Xt, real colvector r, real colvector kv,
                    real scalar kb, real scalar u)
{
    real matrix    H
    real colvector s
    H = invsym(quadcross(Xt, kv, Xt) / kb)
    s = (u :- (r :< 0)) :* (abs(r) :>= 1e-10)
    return((Xt :* s) * H)
}

real rowvector _E_seg(real matrix G, real rowvector ug, real scalar q,
                      real rowvector xt, real scalar ustar, real scalar seg,
                      real scalar tt)
{
    real scalar    ng
    real rowvector Qr
    ng = cols(ug)
    Qr = xt * G'
    if (q <= Qr[1]) {
        seg = 1
        tt  = 0
    }
    else if (q >= Qr[ng]) {
        seg = ng - 1
        tt  = 1
    }
    else {
        seg = max(selectindex(Qr :<= q))
        tt  = (q - Qr[seg]) / (Qr[seg+1] - Qr[seg])
    }
    ustar = ug[seg] + tt * (ug[seg+1] - ug[seg])
    return((1 - tt) * G[seg, .] + tt * G[seg+1, .])
}

// block of one reference: E (k), SE (k), slope t (k), u*, segment
real rowvector _E_pt(real matrix G, real rowvector ug, real scalar q,
                     real rowvector xr)
{
    real scalar    k, ustar, seg, tt
    real rowvector bs
    k     = cols(xr)
    ustar = .
    seg   = .
    tt    = .
    bs = _E_seg(G, ug, q, (xr, 1), ustar, seg, tt)
    return((bs[1..k], J(1, 2*k, .), ustar, seg))
}

real rowvector _E_ref(real matrix G, real rowvector ug, real scalar q,
                      real colvector PSIQ, real matrix X, real rowvector xr,
                      real rowvector flag, real colvector KB)
{
    real scalar    n, k, ustar, seg, tt, D, j, clip
    real rowvector xt, bs, xbar, out, Qr
    real matrix    Xt, P1, P2, Pstar
    real colvector core, psi, dd
    n  = rows(X)
    k  = cols(X)
    Xt = (X, J(n, 1, 1))
    xt = (xr, 1)
    ustar = .
    seg   = .
    tt    = .
    bs = _E_seg(G, ug, q, xt, ustar, seg, tt)
    Qr = xt * G'
    clip = (q <= Qr[1] | q >= Qr[cols(Qr)])
    D  = xt * (G[seg+1, .] - G[seg, .])'
    P1 = _qr_psi(Xt, st_data(., "_r" + strofreal(seg)),     st_data(., "_k" + strofreal(seg)),     KB[seg],     ug[seg])
    P2 = _qr_psi(Xt, st_data(., "_r" + strofreal(seg + 1)), st_data(., "_k" + strofreal(seg + 1)), KB[seg + 1], ug[seg + 1])
    Pstar = (1 - tt) * P1 + tt * P2
    core  = PSIQ - Pstar * xt'
    xbar  = mean(X)
    for (j = 1; j <= k; j++) {
        if (flag[j]) core = core - ((X[., j] :- xbar[j]) / n) * bs[j]
    }
    if (clip) core = J(n, 1, 0)
    out = J(1, 3*k + 2, .)
    for (j = 1; j <= k; j++) {
        psi = Pstar[., j] + ((G[seg+1, j] - G[seg, j]) / D) * core
        psi = psi :- mean(psi)
        out[j]     = bs[j]
        out[k + j] = sqrt(quadcross(psi, psi))
        dd  = P2[., j] - P1[., j]
        dd  = dd :- mean(dd)
        out[2*k + j] = (G[seg+1, j] - G[seg, j]) / sqrt(quadcross(dd, dd))
    }
    out[3*k + 1] = ustar
    out[3*k + 2] = seg
    return(out)
}

// _Eres = block of the default reference [, block of the fixed reference]
void _E_run(real scalar qh, real scalar qkb, real scalar se)
{
    real matrix    X, G
    real rowvector ug, xr, fl, xf, out
    real colvector KB, PSIQ
    real scalar    n, tau
    G   = st_matrix("_Gr")
    ug  = strtoreal(tokens(st_global("UG")))
    tau = strtoreal(st_global("TAU"))
    xr  = st_matrix("_XR")
    fl  = st_matrix("_FL")
    xf  = strtoreal(tokens(st_global("XFIX")))
    if (se) {
        X    = st_data(., tokens(st_global("XV")))
        n    = rows(X)
        KB   = st_matrix("_KB")
        PSIQ = _qr_psi(J(n, 1, 1), st_data(., "_rq"), st_data(., "_kq"), qkb, tau)
        out  = _E_ref(G, ug, qh, PSIQ, X, xr, fl, KB)
        if (cols(xf)) out = (out, _E_ref(G, ug, qh, PSIQ, X, xf, J(1, cols(xf), 0), KB))
    }
    else {
        out = _E_pt(G, ug, qh, xr)
        if (cols(xf)) out = (out, _E_pt(G, ug, qh, xf))
    }
    st_matrix("_Eres", out)
}

// the default reference type of the DGP: x1 at its mean, x2 at its mode
void _ref_dgp()
{
    real colvector x1, x2
    x1 = st_data(., "x1")
    x2 = st_data(., "x2")
    st_matrix("_XR", (mean(x1), mean(x2) > 0.5))
    st_matrix("_FL", (1, 0))
}

// Burkina: size at its mean, male and urban at their mode, gse at its most
// frequent level (dummies g2..g7)
void _ref_bkf()
{
    real colvector s, m, u, g
    real rowvector cnt, xr, top
    real scalar    l, gm
    s = st_data(., "size")
    m = st_data(., "male")
    u = st_data(., "urban")
    g = st_data(., "gse")
    cnt = J(1, 7, 0)
    for (l = 1; l <= 7; l++) cnt[l] = sum(g :== l)
    top = selectindex(cnt :== max(cnt))
    gm  = top[1]
    xr = (mean(s), mean(m) > 0.5, mean(u) > 0.5)
    for (l = 2; l <= 7; l++) xr = (xr, gm == l)
    st_matrix("_XR", xr)
    st_matrix("_FL", (1, J(1, 8, 0)))
    st_numscalar("_gmode", gm)
}

real colvector _sq(real scalar q, real colvector x1, real colvector x2,
                   real scalar s1, real scalar s2)
{
    real colvector lo, hi, mid, up
    real scalar    it
    lo = J(rows(x1), 1, 1e-12)
    hi = J(rows(x1), 1, 1 - 1e-12)
    for (it = 1; it <= 70; it++) {
        mid = (lo + hi) :/ 2
        up  = (invnormal(mid) + (0.4 :+ s1 :* mid) :* x1 + (0.5 :+ s2 :* mid) :* x2) :< q
        lo  = up :* mid + (1 :- up) :* lo
        hi  = up :* hi + (1 :- up) :* mid
    }
    return((lo + hi) :/ 2)
}

// truth of the DGP, deterministic: plim E for the default type (5, mode)
// then the fixed type (5, 1); also q_tau, u* and beta(u*)
real rowvector _truth(real scalar tau, real scalar s1, real scalar s2)
{
    real scalar    M, it, lo, hi, mid, F, q, mode, ud, uf, us, sg, tt
    real colvector x1, p, z0, z1
    real rowvector ug, pd, pf
    real matrix    Gt
    M  = 4000
    x1 = 1 :+ 8 :* ((1::M) :- 0.5) :/ M
    p  = invlogit(-1 :+ 0.1 :* x1)
    z0 = J(M, 1, 0)
    z1 = J(M, 1, 1)
    mode = (mean(p) > 0.5)
    lo = -20
    hi = 60
    for (it = 1; it <= 70; it++) {
        mid = (lo + hi) / 2
        F   = mean((1 :- p) :* _sq(mid, x1, z0, s1, s2) + p :* _sq(mid, x1, z1, s1, s2))
        if (F < tau) lo = mid
        else         hi = mid
    }
    q  = (lo + hi) / 2
    ud = _sq(q, 5, mode, s1, s2)
    uf = _sq(q, 5, 1, s1, s2)
    ug = (1..19) :* 0.05
    Gt = ((0.4 :+ s1 :* ug)', (0.5 :+ s2 :* ug)', invnormal(ug)')
    us = .
    sg = .
    tt = .
    pd = _E_seg(Gt, ug, q, (5, mode, 1), us, sg, tt)
    pf = _E_seg(Gt, ug, q, (5, 1, 1), us, sg, tt)
    return((pd[1], pd[2], pf[1], pf[2], q, ud, uf, 0.4 + s1*ud, 0.5 + s2*ud, 0.4 + s1*uf, 0.5 + s2*uf))
}

// report: E0 the point row, TT the draws, names of the coefficients, nref
// blocks, plim (empty if unknown), labels of the references
void _report(real rowvector E0, real matrix TT, string rowvector names,
             real scalar nref, real rowvector plim, string rowvector rlab)
{
    real scalar    k, r, j, B, base, ie, sd, rob, cr
    real colvector e, srt
    k = cols(names)
    B = rows(TT)
    for (r = 1; r <= nref; r++) {
        base = (r - 1) * (3*k + 2)
        printf("\n  reference: %s -- u* %6.4f (segment %g), draws whose segment changed: %4.1f%%\n",
               rlab[r], E0[base + 3*k + 1], E0[base + 3*k + 2], 100 * mean(TT[., base + 3*k + 2] :!= E0[base + 3*k + 2]))
        if (cols(plim)) printf("  coefficient     point     plim      t   SE_anal   SD_boot  robust SD  anal/SD  anal/rob  slope t\n")
        else            printf("  coefficient     point   SE_anal   SD_boot  robust SD  anal/SD  anal/rob  slope t\n")
        for (j = 1; j <= k; j++) {
            ie  = base + j
            e   = TT[., ie]
            sd  = sqrt(variance(e))
            srt = sort(e, 1)
            rob = (srt[ceil(0.75*B)] - srt[ceil(0.25*B)]) / 1.349
            if (cols(plim)) {
                printf("  %-12s %8.4f %8.4f %6.2f  %8.4f  %8.4f  %8.4f   %6.3f   %6.3f   %6.2f\n",
                       names[j], E0[ie], plim[(r-1)*k + j], (E0[ie] - plim[(r-1)*k + j]) / E0[ie + k],
                       E0[ie + k], sd, rob, E0[ie + k] / sd, E0[ie + k] / rob, E0[ie + 2*k])
            }
            else {
                printf("  %-12s %8.4f  %8.4f  %8.4f  %8.4f   %6.3f   %6.3f   %6.2f\n",
                       names[j], E0[ie], E0[ie + k], sd, rob, E0[ie + k] / sd, E0[ie + k] / rob, E0[ie + 2*k])
            }
        }
    }
}
end

* every step of E on the data in memory -> matrix _Eres; se: 1 = with the
* analytic SE, 0 = the value only (the draws)
capture program drop _E_full
program define _E_full
    args se
    quietly {
        local k : word count $XV
        local ng : word count $UG
        matrix _Gr = J(`ng', `k' + 1, .)
        matrix _KB = J(`ng', 1, .)
        local g 0
        foreach u of global UG {
            local ++g
            if `se' {
                capture qreg $Y $XV, quantile(`u') vce(robust, kernel)
                if _rc qreg $Y $XV, quantile(`u') vce(robust, kernel) iterate(5000)
            }
            else {
                capture _qreg $Y $XV, quant(`u')
                if _rc _qreg $Y $XV, quant(`u') iterate(5000)
            }
            local j 0
            local xb 0
            foreach v of global XV {
                local ++j
                matrix _Gr[`g', `j'] = _b[`v']
            }
            matrix _Gr[`g', `k' + 1] = _b[_cons]
            if `se' {
                local kb = e(kbwidth)
                matrix _KB[`g', 1] = `kb'
                capture drop _r`g' _k`g'
                predict double _r`g', residuals
                * qreg's Epanechnikov kernel (_qreg_kernel, private to qreg.ado)
                generate double _k`g' = 0.75*(1 - 0.2*(_r`g'/`kb')^2)/sqrt(5) * (abs(_r`g'/`kb') < sqrt(5))
            }
        }
        local qkb 0
        if `se' {
            capture qreg $Y, quantile($TAU) vce(robust, kernel)
            if _rc qreg $Y, quantile($TAU) vce(robust, kernel) iterate(5000)
            local qh  = _b[_cons]
            local qkb = e(kbwidth)
            capture drop _rq _kq
            generate double _rq = $Y - `qh'
            generate double _kq = 0.75*(1 - 0.2*(_rq/`qkb')^2)/sqrt(5) * (abs(_rq/`qkb') < sqrt(5))
        }
        else {
            capture _qreg $Y, quant($TAU)
            if _rc _qreg $Y, quant($TAU) iterate(5000)
            local qh = _b[_cons]
        }
        * the default reference type, from this sample
        mata: _ref_$REFTYPE()
        mata: _E_run(`qh', `qkb', `se')
        if `se' drop _r* _k*
    }
end

capture program drop _dgp
program define _dgp
    args nobs s1 s2
    quietly {
        set obs `nobs'
        generate double x1 = 1 + 8*runiform()
        generate double x2 = runiform() < invlogit(-1 + 0.1*x1)
        generate double u  = runiform()
        generate double y  = invnormal(u) + (0.4 + `s1'*u)*x1 + (0.5 + `s2'*u)*x2
        drop u
    }
end

capture program drop _bkf_data
program define _bkf_data
    use bkf98I, clear
    quietly {
        generate double lexp = ln(exppc)
        generate byte male  = (sex == 1)
        generate byte urban = (zone == 2)
        forvalues l = 2/7 {
            generate byte g`l' = (gse == `l')
        }
        keep lexp size male urban gse g2-g7
        keep if !missing(lexp, size, male, urban, gse)
    }
end

* ---- the perfect case ------------------------------------------------------------------
if "`what'" == "dgp" {
    global Y y
    global XV x1 x2
    global XFIX 5 1
    global REFTYPE dgp
    mata: TR = _truth(`tau', 2, 3)
    mata: printf("\ntruth: q_tau %8.4f; default type (5, 0): u* %6.4f, beta(u*) (%7.4f, %7.4f), plim E (%7.4f, %7.4f)\n", TR[5], TR[6], TR[8], TR[9], TR[1], TR[2])
    mata: printf("       fixed type (5, 1): u* %6.4f, beta(u*) (%7.4f, %7.4f), plim E (%7.4f, %7.4f)\n", TR[7], TR[10], TR[11], TR[3], TR[4])
    set seed 1
    _dgp 2000 2 3
    _E_full 1
    matrix E0 = _Eres
    keep y x1 x2
    matrix TT = J(`B', colsof(E0), .)
    set seed 777
    timer clear 1
    timer on 1
    forvalues b = 1/`B' {
        preserve
        bsample
        _E_full 0
        restore
        matrix TT[`b', 1] = _Eres
    }
    timer off 1
    quietly timer list 1
    di as text "bootstrap time: " as result %6.1f r(t1) as text " s"
    di as text _n "==== perfect case: beta_x1(u) = 0.4 + 2u, beta_x2(u) = 0.5 + 3u, n = 2000, tau = `tau', B = `B'"
    mata: _report(st_matrix("E0"), st_matrix("TT"), ("x1", "x2"), 2, TR[1..4], ("default type (x1 mean, x2 mode)", "fixed type (5, 1)"))
    preserve
    clear
    quietly svmat double TT, names(c)
    save "`OUT'/p13_dgp_draws.dta", replace
    restore
    exit
}

* ---- Burkina: one chunk of the draws ---------------------------------------------------
if "`what'" == "bkf" {
    global Y lexp
    global XV size male urban g2 g3 g4 g5 g6 g7
    global XFIX
    global REFTYPE bkf
    _bkf_data
    if `chunk' == 1 {
        timer clear 2
        timer on 2
        _E_full 1
        timer off 2
        quietly timer list 2
        di as text "point with its SE: " as result %6.1f r(t2) as text " s; most frequent gse level " as result _gmode
        matrix E0 = _Eres
        preserve
        clear
        quietly svmat double E0, names(c)
        save "`OUT'/p13_bkf_point.dta", replace
        restore
    }
    local b1 = floor((`chunk' - 1) * `B' / `nch') + 1
    local b2 = floor(`chunk' * `B' / `nch')
    local nb = `b2' - `b1' + 1
    matrix TT = J(`nb', 29, .)
    set seed `= 777 + `chunk''
    timer clear 1
    timer on 1
    forvalues b = 1/`nb' {
        preserve
        bsample
        _E_full 0
        restore
        matrix TT[`b', 1] = _Eres
        if mod(`b', 20) == 0 di as text "  chunk `chunk': `b' of `nb'"
    }
    timer off 1
    quietly timer list 1
    di as text "chunk `chunk' (`nb' draws): " as result %6.1f r(t1) as text " s"
    clear
    quietly svmat double TT, names(c)
    save "`OUT'/p13_bkf_draws_c`chunk'.dta", replace
    exit
}

* ---- Burkina: combine and report ---------------------------------------------------------
if "`what'" == "bkfreport" {
    use "`OUT'/p13_bkf_point.dta", clear
    mkmat c*, matrix(E0)
    use "`OUT'/p13_bkf_draws_c1.dta", clear
    forvalues c = 2/`nch' {
        append using "`OUT'/p13_bkf_draws_c`c'.dta"
    }
    mkmat c*, matrix(TT)
    di as text _n "==== Burkina (bkf98I), lexp on size male urban i.gse, unweighted, n = 8478, tau = `tau', B = " rowsof(TT)
    mata: _report(st_matrix("E0"), st_matrix("TT"), ("size", "male", "urban", "2.gse", "3.gse", "4.gse", "5.gse", "6.gse", "7.gse"), 1, J(1, 0, .), ("default type (size mean; male, urban, gse at their mode)"))
    exit
}
di as error "usage: do p13a_E_draws.do dgp | bkf c | bkfreport"
exit 198
