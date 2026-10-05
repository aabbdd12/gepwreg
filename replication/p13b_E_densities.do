* p13b_E_densities.do -- called by p13_profile_delta.do (Paper, version 5,
* Section 7.5): the delta-method standard error of E with two densities for
* the quantile-regression influence function, judged by the SAME draws of
* p13a_E_draws.do (the draws do not depend on how the standard error is
* computed):
*   Powell kernel     : one bandwidth for every household (qreg vce(robust, kernel));
*   Hendricks-Koenker : each household's own density 2h/(x'b(u+h) - x'b(u-h))
*                       (qreg's default vce(robust)); Siddiqui quotient for q_tau.
* Argument: dgp | bkf.  Checks: the Stata form u(1-u) H X'X H of both
* reproduces e(V) of qreg; the point E equals p13a's.  Report: the two
* standard errors against the standard deviation of the draws -- the numbers
* of Section 7.5.  It does not call gepwreg.
version 16
clear all
args what
set more off
set varabbrev off
local OUT "results"
capture mkdir "results"
local tau 0.25
global TAU `tau'
numlist "0.05(0.05)0.95"
global UG `r(numlist)'

mata:
real matrix _psiK(real matrix Xt, real colvector r, real colvector kv,
                  real scalar kb, real scalar u)
{
    real matrix    H
    real colvector s
    H = invsym(quadcross(Xt, kv, Xt) / kb)
    s = (u :- (r :< 0)) :* (abs(r) :>= 1e-10)
    return((Xt :* s) * H)
}

real matrix _psiF(real matrix Xt, real colvector r, real colvector f, real scalar u)
{
    real matrix    H
    real colvector s
    H = invsym(quadcross(Xt, f, Xt))
    s = (u :- (r :< 0)) :* (abs(r) :>= 1e-10)
    return((Xt :* s) * H)
}

// check (1): both Stata forms against e(V) of the two qreg runs; row g of _CHK
void _chk(real scalar g, real scalar u, string scalar vk, string scalar vf,
          real scalar kb, string scalar rname, string scalar kname,
          string scalar fname, string scalar xv)
{
    real matrix    X, Xt, HK, HF, C
    real colvector kv, f
    real scalar    n
    X  = (strlen(xv) ? st_data(., tokens(xv)) : J(st_nobs(), 0, .))
    n  = rows(X)
    Xt = (X, J(n, 1, 1))
    kv = st_data(., kname)
    f  = st_data(., fname)
    HK = invsym(quadcross(Xt, kv, Xt) / kb)
    HF = invsym(quadcross(Xt, f, Xt))
    C  = st_matrix("_CHK")
    C[g, 1] = max(abs(sqrt(diagonal(u*(1-u) * HK * quadcross(Xt, Xt) * HK)) :/ sqrt(diagonal(st_matrix(vk))) :- 1))
    C[g, 2] = max(abs(sqrt(diagonal(u*(1-u) * HF * quadcross(Xt, Xt) * HF)) :/ sqrt(diagonal(st_matrix(vf))) :- 1))
    st_matrix("_CHK", C)
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

// block of one reference: E (k), SE kernel (k), SE Hendricks-Koenker (k), u*, segment
real rowvector _E_ref2(real matrix G, real rowvector ug, real scalar q,
                       real colvector PQK, real colvector PQF, real matrix X,
                       real rowvector xr, real rowvector flag, real colvector KB)
{
    real scalar    n, k, ustar, seg, tt, D, j, clip, m
    real rowvector xt, bs, xbar, out, Qr
    real matrix    Xt, P1, P2, Pstar
    real colvector core, psi, r1, r2
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
    r1 = st_data(., "_r" + strofreal(seg))
    r2 = st_data(., "_r" + strofreal(seg + 1))
    xbar = mean(X)
    out  = J(1, 3*k + 2, .)
    out[1..k] = bs[1..k]
    for (m = 1; m <= 2; m++) {
        if (m == 1) {
            P1 = _psiK(Xt, r1, st_data(., "_k" + strofreal(seg)),     KB[seg],     ug[seg])
            P2 = _psiK(Xt, r2, st_data(., "_k" + strofreal(seg + 1)), KB[seg + 1], ug[seg + 1])
            core = PQK
        }
        else {
            P1 = _psiF(Xt, r1, st_data(., "_f" + strofreal(seg)),     ug[seg])
            P2 = _psiF(Xt, r2, st_data(., "_f" + strofreal(seg + 1)), ug[seg + 1])
            core = PQF
        }
        Pstar = (1 - tt) * P1 + tt * P2
        core  = core - Pstar * xt'
        for (j = 1; j <= k; j++) {
            if (flag[j]) core = core - ((X[., j] :- xbar[j]) / n) * bs[j]
        }
        if (clip) core = J(n, 1, 0)
        for (j = 1; j <= k; j++) {
            psi = Pstar[., j] + ((G[seg+1, j] - G[seg, j]) / D) * core
            psi = psi :- mean(psi)
            out[m*k + j] = sqrt(quadcross(psi, psi))
        }
    }
    out[3*k + 1] = ustar
    out[3*k + 2] = seg
    return(out)
}

void _E_run2(real scalar qh, real scalar qkb)
{
    real matrix    X, G
    real rowvector ug, xr, fl, xf, out
    real colvector KB, PQK, PQF
    real scalar    n, tau
    G   = st_matrix("_Gr")
    ug  = strtoreal(tokens(st_global("UG")))
    tau = strtoreal(st_global("TAU"))
    xr  = st_matrix("_XR")
    fl  = st_matrix("_FL")
    xf  = strtoreal(tokens(st_global("XFIX")))
    X   = st_data(., tokens(st_global("XV")))
    n   = rows(X)
    KB  = st_matrix("_KB")
    PQK = _psiK(J(n, 1, 1), st_data(., "_rq"), st_data(., "_kq"), qkb, tau)
    PQF = _psiF(J(n, 1, 1), st_data(., "_rq"), st_data(., "_fq"), tau)
    out = _E_ref2(G, ug, qh, PQK, PQF, X, xr, fl, KB)
    if (cols(xf)) out = (out, _E_ref2(G, ug, qh, PQK, PQF, X, xf, J(1, cols(xf), 0), KB))
    st_matrix("_Eres", out)
}

void _ref_dgp()
{
    real colvector x1, x2
    x1 = st_data(., "x1")
    x2 = st_data(., "x2")
    st_matrix("_XR", (mean(x1), mean(x2) > 0.5))
    st_matrix("_FL", (1, 0))
}

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
}

void _report4(real rowvector E0, real matrix TT, string rowvector names,
              real scalar nref, string rowvector rlab)
{
    real scalar    k, r, j, B, base, ie, sd, rob
    real colvector e, srt
    k = cols(names)
    B = rows(TT)
    for (r = 1; r <= nref; r++) {
        base = (r - 1) * (3*k + 2)
        printf("\n  reference: %s -- u* %6.4f (segment %g)\n", rlab[r], E0[base + 3*k + 1], E0[base + 3*k + 2])
        printf("  coefficient     point   SE_kern     SE_HK   SD_boot  robust SD  kern/SD    HK/SD   HK/rob\n")
        for (j = 1; j <= k; j++) {
            ie  = base + j
            e   = TT[., ie]
            sd  = sqrt(variance(e))
            srt = sort(e, 1)
            rob = (srt[ceil(0.75*B)] - srt[ceil(0.25*B)]) / 1.349
            printf("  %-12s %8.4f  %8.4f  %8.4f  %8.4f   %8.4f   %6.3f   %6.3f   %6.3f\n",
                   names[j], E0[ie], E0[ie + k], E0[ie + 2*k], sd, rob,
                   E0[ie + k] / sd, E0[ie + 2*k] / sd, E0[ie + 2*k] / rob)
        }
    }
}
end

* the point of E with both SEs -> _Eres; _CHK holds check (1)
capture program drop _E_point2
program define _E_point2
    quietly {
        local k  : word count $XV
        local ng : word count $UG
        matrix _Gr  = J(`ng', `k' + 1, .)
        matrix _KB  = J(`ng', 1, .)
        matrix _CHK = J(`ng' + 1, 2, .)
        tempname VK VF
        local g 0
        foreach u of global UG {
            local ++g
            qreg $Y $XV, quantile(`u') vce(robust, kernel)
            matrix `VK' = e(V)
            local kb = e(kbwidth)
            matrix _KB[`g', 1] = `kb'
            capture drop _r`g' _k`g' _f`g'
            predict double _r`g', residuals
            generate double _k`g' = 0.75*(1 - 0.2*(_r`g'/`kb')^2)/sqrt(5) * (abs(_r`g'/`kb') < sqrt(5))
            qreg $Y $XV, quantile(`u') vce(robust)
            matrix `VF' = e(V)
            local j 0
            foreach v of global XV {
                local ++j
                matrix _Gr[`g', `j'] = _b[`v']
            }
            matrix _Gr[`g', `k' + 1] = _b[_cons]
            local h = e(bwidth)
            * qreg's fitted method (VCE_fitted in qreg.ado)
            _qreg $Y $XV, quant(`=`u' - `h'')
            predict double _xl, xb
            _qreg $Y $XV, quant(`=`u' + `h'')
            predict double _xu, xb
            generate double _f`g' = cond(_xu - _xl > sqrt(c(epsdouble)), 2*`h'/(_xu - _xl), 0)
            drop _xl _xu
            mata: _chk(`g', `u', "`VK'", "`VF'", `kb', "_r`g'", "_k`g'", "_f`g'", "$XV")
        }
        local g = `ng' + 1
        qreg $Y, quantile($TAU) vce(robust, kernel)
        matrix `VK' = e(V)
        local qh  = _b[_cons]
        local qkb = e(kbwidth)
        capture drop _rq _kq _fq
        generate double _rq = $Y - `qh'
        generate double _kq = 0.75*(1 - 0.2*(_rq/`qkb')^2)/sqrt(5) * (abs(_rq/`qkb') < sqrt(5))
        qreg $Y, quantile($TAU) vce(robust)
        matrix `VF' = e(V)
        local hq = e(bwidth)
        _qreg $Y, quant(`=$TAU - `hq'')
        local ql = _b[_cons]
        _qreg $Y, quant(`=$TAU + `hq'')
        local qu = _b[_cons]
        generate double _fq = cond(`qu' - `ql' > sqrt(c(epsdouble)), 2*`hq'/(`qu' - `ql'), 0)
        mata: _chk(`g', $TAU, "`VK'", "`VF'", `qkb', "_rq", "_kq", "_fq", "")
        mata: _ref_$REFTYPE()
        mata: _E_run2(`qh', `qkb')
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

if "`what'" == "dgp" {
    global Y y
    global XV x1 x2
    global XFIX 5 1
    global REFTYPE dgp
    set seed 1
    _dgp 2000 2 3
    _E_point2
    matrix E0 = _Eres
    mata: C = st_matrix("_CHK"); printf("\n(1) Stata forms against e(V), largest relative SE gap over 19 grid points and q: kernel %9.2e, fitted %9.2e\n", max(C[., 1]), max(C[., 2]))
    di as text "(2) p13a's point (printed by p13a above): default (0.8016, 1.5653), fixed (0.6198, 1.1077)"
    mata: E0 = st_matrix("E0"); printf("    p13b's point: default (%6.4f, %6.4f), fixed (%6.4f, %6.4f)\n", E0[1], E0[2], E0[9], E0[10])
    use "`OUT'/p13_dgp_draws.dta", clear
    mkmat c*, matrix(TT)
    di as text _n "(3) perfect case, n = 2000, tau = `tau', the 200 draws of p13a"
    mata: _report4(st_matrix("E0"), st_matrix("TT"), ("x1", "x2"), 2, ("default type (x1 mean, x2 mode)", "fixed type (5, 1)"))
    exit
}
if "`what'" == "bkf" {
    global Y lexp
    global XV size male urban g2 g3 g4 g5 g6 g7
    global XFIX
    global REFTYPE bkf
    _bkf_data
    _E_point2
    matrix E0 = _Eres
    mata: C = st_matrix("_CHK"); printf("\n(1) Stata forms against e(V), largest relative SE gap over 19 grid points and q: kernel %9.2e, fitted %9.2e\n", max(C[., 1]), max(C[., 2]))
    preserve
    use "`OUT'/p13_bkf_point.dta", clear
    mkmat c*, matrix(E3P)
    restore
    mata: a = st_matrix("E0"); b = st_matrix("E3P"); printf("(2) point E against p13a's, largest gap over the 9 coefficients: %9.2e\n", max(abs(a[1..9] - b[1..9])))
    use "`OUT'/p13_bkf_draws_c1.dta", clear
    append using "`OUT'/p13_bkf_draws_c2.dta"
    append using "`OUT'/p13_bkf_draws_c3.dta"
    mkmat c*, matrix(TT)
    di as text _n "(3) Burkina, unweighted, n = 8478, tau = `tau', the 200 draws of p13a"
    mata: _report4(st_matrix("E0"), st_matrix("TT"), ("size", "male", "urban", "2.gse", "3.gse", "4.gse", "5.gse", "6.gse", "7.gse"), 1, ("default type (size mean; male, urban, gse at their mode)"))
    exit
}
di as error "usage: do p13b_E_densities.do dgp | bkf"
exit 198
