* p12_profile_dgp.do -- Paper (version 5), Section 7.4, Table 9: the linear
* correction (L) against the exact estimator (E) of the effect for a given
* household under het(qr), where the truth is known.  tau = 0.25.
*   L: the households' effects b_i = beta(u_i) (gepwreg's het(qr): conditional
*      rank by the count of grid predictions below y, midpoint, interpolation),
*      regressed on (x, 1) with the percentile weights at tau, predicted at
*      the reference:   L_j = (xref, 1) (Xt'WXt)^-1 Xt'W b_j.
*   E: beta(u*) with xt_ref'beta(u*) = q_tau on the interpolated grid.
* References: the type (5, 1) and the mean household.
*   (1) bias: L and E on one large sample (n = 200,000) against the truth
*       (2,000,000 draws);
*   (2) variance: pairs bootstrap at n = 2000, B = 400, the whole procedure
*       re-run, the bandwidth held at the point estimate's.
* DGP: y = Phi^-1(u) + (0.4 + 0.4u) x1 + (0.5 + 0.6u) x2, u ~ U(0,1)
* independent of (x1, x2), x1 ~ U(1,9), Pr[x2 = 1 | x1] = invlogit(-1 + 0.1 x1).
* Needs p12_core_mata.do (beside this file).
* Output: results/p12_profile_dgp.csv (rows x1, x2: the truth, L and E for
* the two references, and the bootstrap standard errors of L and E).
* Run: do p12_profile_dgp   (about ten minutes: 400 x 19 quantile regressions)

* ---------------------------------------------------------------------
* Replication, Araar (2026), version 5.  Run this file from the folder
* that contains it; the results are written to the results/ folder
* beside it, which is created if it does not exist.
*
* It needs gepwreg 1.5.0:
*
*   net install gepwreg, from("https://raw.githubusercontent.com/aabbdd12/gepwreg/v15") replace
*
* The which gepwreg below must report v1.5.0.
* ---------------------------------------------------------------------
version 16
clear all
set more off
set varabbrev off
which gepwreg
capture mkdir "results"
run "p12_core_mata.do"

local tau  0.25
local B    400
local xr1  5
local xr2  1

capture program drop _dgp
program define _dgp
    args nobs
    quietly {
        set obs `nobs'
        generate double x1 = 1 + 8*runiform()
        generate double x2 = runiform() < invlogit(-1 + 0.1*x1)
        generate double u  = runiform()
        generate double y  = invnormal(u) + (0.4 + 0.4*u)*x1 + (0.5 + 0.6*u)*x2
        drop u
    }
end

mata:
real scalar _ustar(real scalar x1, real scalar x2, real scalar q)
{
    real scalar lo, hi, mid, it
    lo = 1e-9
    hi = 1 - 1e-9
    for (it = 1; it <= 200; it++) {
        mid = (lo + hi) / 2
        if (invnormal(mid) + (0.4 + 0.4*mid)*x1 + (0.5 + 0.6*mid)*x2 < q) lo = mid
        else hi = mid
    }
    return((lo + hi) / 2)
}

// E: beta(u*) at the reference xt (with its constant); also u*
real rowvector _E_pt(real matrix G, real rowvector ug, real scalar q,
                     real rowvector xt, real scalar ustar)
{
    real scalar    ng, seg, tt
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

// theta(tau), L (ref fixed, ref mean), E (ref fixed, ref mean), u* (fixed, mean)
real rowvector _L_E(real colvector y, real matrix X, real colvector fw,
                    real matrix G, real rowvector ug, real scalar tau,
                    real scalar h, real scalar q, real rowvector xf)
{
    real scalar    n, k, ng, uf, um
    real colvector c, ext, u, pc, w, ww
    real matrix    Bi, Xt, Pi
    real rowvector umid, th, xm, Lf, Lm, Ef, Em
    n    = rows(y)
    k    = cols(X)
    ng   = cols(ug)
    c    = rowsum(((X, J(n, 1, 1)) * G') :< y)
    umid = (ug[1..(ng-1)] + ug[2..ng]) :/ 2
    ext  = (ug[1], umid, ug[ng])'
    u    = ext[c :+ 1]
    Bi   = _rfc_interp(G[., 1..k], u, ug)
    pc   = _rfc_wp(y, fw, tau, 0.9, h, J(0, 1, .))[., 2]
    w    = _rfc_kw(pc, tau, h, n)
    ww   = fw :* w
    th   = (ww' * Bi) / sum(ww)
    Xt   = (X, J(n, 1, 1))
    Pi   = invsym(quadcross(Xt, ww, Xt)) * quadcross(Xt, ww, Bi)
    xm   = (fw' * X) / sum(fw)
    Lf   = (xf, 1) * Pi
    Lm   = (xm, 1) * Pi
    uf   = .
    um   = .
    Ef   = _E_pt(G, ug, q, (xf, 1), uf)
    Em   = _E_pt(G, ug, q, (xm, 1), um)
    return((th, Lf, Lm, Ef[1..k], Em[1..k], uf, um))
}
end

* ---- the truth ---------------------------------------------------------------------
set seed 20261005
_dgp 2000000
_pctile y, percentiles(`=100*`tau'')
local qpop = r(r1)
summarize x2, meanonly
local m2 = r(mean)
mata: uf = _ustar(`xr1', `xr2', `qpop'); um = _ustar(5, `m2', `qpop')
mata: TR = (0.4 + 0.4*uf, 0.5 + 0.6*uf, 0.4 + 0.4*um, 0.5 + 0.6*um)
mata: printf("\ntruth: u*(fixed) %6.4f  u*(mean) %6.4f   fixed (x1, x2) = (%6.4f, %6.4f)   mean (x1, x2) = (%6.4f, %6.4f)\n", uf, um, TR[1], TR[2], TR[3], TR[4])

* ---- (1) bias: one large sample ---------------------------------------------------
clear
set seed 2
_dgp 200000
generate double one = 1
quietly gepwreg y x1 x2, per(`tau')
local hL = e(h)
tempname GL bL
matrix `GL' = e(b_qgrid)
matrix `bL' = e(b)
local ugL "`e(qgrid)'"
quietly _qreg y, quant(`tau')
local qL = _b[_cons]
mata: yL = st_data(., "y"); XL = st_data(., ("x1", "x2")); fL = st_data(., "one")
mata: ugrid = strtoreal(tokens("`ugL'"))
mata: RL = _L_E(yL, XL, fL, st_matrix("`GL'"), ugrid, `tau', `hL', `qL', (`xr1', `xr2'))
mata: printf("\n(1) n = 200,000, h = %6.4f: theta vs gepwreg %9.2e\n", `hL', max(abs(RL[1..2] - st_matrix("`bL'"))))
mata: printf("    u*: fixed %6.4f (truth %6.4f)   mean %6.4f (truth %6.4f)\n", RL[11], uf, RL[12], um)
mata: printf("                 theta(tau)   L fixed     L mean      E fixed     E mean\n")
mata: printf("    x1        %10.4f  %10.4f  %10.4f  %10.4f  %10.4f\n", RL[1], RL[3], RL[5], RL[7], RL[9])
mata: printf("    x2        %10.4f  %10.4f  %10.4f  %10.4f  %10.4f\n", RL[2], RL[4], RL[6], RL[8], RL[10])
mata: printf("    truth x1               %10.4f  %10.4f  %10.4f  %10.4f\n", TR[1], TR[3], TR[1], TR[3])
mata: printf("    truth x2               %10.4f  %10.4f  %10.4f  %10.4f\n", TR[2], TR[4], TR[2], TR[4])
mata: BIAS = (RL[3], RL[4], RL[5], RL[6], RL[7], RL[8], RL[9], RL[10]) :- (TR[1], TR[2], TR[3], TR[4], TR[1], TR[2], TR[3], TR[4])

* ---- (2) variance: n = 2000, pairs bootstrap ------------------------------------
clear
set seed 1
_dgp 2000
generate double one = 1
quietly gepwreg y x1 x2, per(`tau')
local h = e(h)
tempname G0 b0
matrix `G0' = e(b_qgrid)
matrix `b0' = e(b)
local ug "`e(qgrid)'"
quietly _qreg y, quant(`tau')
local q0 = _b[_cons]
mata: y0 = st_data(., "y"); X0 = st_data(., ("x1", "x2")); f0 = st_data(., "one"); ug = strtoreal(tokens("`ug'"))
mata: R0 = _L_E(y0, X0, f0, st_matrix("`G0'"), ug, `tau', `h', `q0', (`xr1', `xr2'))
mata: printf("\n(2) n = 2000, h = %6.4f: theta vs gepwreg %9.2e\n", `h', max(abs(R0[1..2] - st_matrix("`b0'"))))

keep y x1 x2 one
tempname TT Gb
matrix `TT' = J(`B', 12, .)
set seed 777
forvalues b = 1/`B' {
    preserve
    bsample
    matrix `Gb' = J(19, 3, .)
    local g 0
    foreach u of local ug {
        local ++g
        quietly _qreg y x1 x2, quant(`u')
        matrix `Gb'[`g', 1] = _b[x1]
        matrix `Gb'[`g', 2] = _b[x2]
        matrix `Gb'[`g', 3] = _b[_cons]
    }
    quietly _qreg y, quant(`tau')
    local qb = _b[_cons]
    mata: st_matrix("_rowb", _L_E(st_data(., "y"), st_data(., ("x1", "x2")), st_data(., "one"), st_matrix("`Gb'"), ug, `tau', `h', `qb', (`xr1', `xr2')))
    restore
    matrix `TT'[`b', 1] = _rowb
    if mod(`b', 100) == 0 di as text "  `b' replications"
}
mata:
    TT  = st_matrix("`TT'")
    sb  = sqrt(diagonal(variance(TT)))'
    printf("\n    bootstrap SE (B = `B')      theta(tau)   L fixed     L mean      E fixed     E mean\n")
    printf("    x1                      %10.4f  %10.4f  %10.4f  %10.4f  %10.4f\n", sb[1], sb[3], sb[5], sb[7], sb[9])
    printf("    x2                      %10.4f  %10.4f  %10.4f  %10.4f  %10.4f\n", sb[2], sb[4], sb[6], sb[8], sb[10])
    printf("    SE ratio L / E          x1 fixed %6.3f   x2 fixed %6.3f   x1 mean %6.3f   x2 mean %6.3f\n", sb[3]/sb[7], sb[4]/sb[8], sb[5]/sb[9], sb[6]/sb[10])
    printf("    stability of the L SEs: B = 100 / 200 / 400\n")
    for (j = 3; j <= 6; j++) printf("      col %g: %7.4f %7.4f %7.4f\n", j, sqrt(variance(TT[1..100, j])), sqrt(variance(TT[1..200, j])), sb[j])
    sbe = (sb[3], sb[4], sb[5], sb[6], sb[7], sb[8], sb[9], sb[10])
    printf("\n    large-sample bias / bootstrap SE at n = 2000:\n")
    printf("                 L fixed   L mean    E fixed   E mean\n")
    printf("    x1          %7.3f   %7.3f   %7.3f   %7.3f\n", BIAS[1]/sbe[1], BIAS[3]/sbe[3], BIAS[5]/sbe[5], BIAS[7]/sbe[7])
    printf("    x2          %7.3f   %7.3f   %7.3f   %7.3f\n", BIAS[2]/sbe[2], BIAS[4]/sbe[4], BIAS[6]/sbe[6], BIAS[8]/sbe[8])
    printf("\n    the n = 2000 sample against the truth (t with the bootstrap SE):\n")
    printf("    L fixed  x1 %7.4f (t %5.2f)  x2 %7.4f (t %5.2f)   L mean  x1 %7.4f (t %5.2f)  x2 %7.4f (t %5.2f)\n", R0[3], (R0[3]-TR[1])/sb[3], R0[4], (R0[4]-TR[2])/sb[4], R0[5], (R0[5]-TR[3])/sb[5], R0[6], (R0[6]-TR[4])/sb[6])
    printf("    E fixed  x1 %7.4f (t %5.2f)  x2 %7.4f (t %5.2f)   E mean  x1 %7.4f (t %5.2f)  x2 %7.4f (t %5.2f)\n", R0[7], (R0[7]-TR[1])/sb[7], R0[8], (R0[8]-TR[2])/sb[8], R0[9], (R0[9]-TR[3])/sb[9], R0[10], (R0[10]-TR[4])/sb[10])
end

* ---- the results ---------------------------------------------------------------
mata: st_matrix("_p12", ((TR[1], RL[3], RL[7], TR[3], RL[5], RL[9], sb[3], sb[7], sb[5], sb[9]) \ (TR[2], RL[4], RL[8], TR[4], RL[6], RL[10], sb[4], sb[8], sb[6], sb[10])))
matrix colnames _p12 = truth_fixed L_fixed E_fixed truth_mean L_mean E_mean se_L_fixed se_E_fixed se_L_mean se_E_mean
preserve
clear
quietly svmat double _p12, names(col)
generate str2 coef = cond(_n == 1, "x1", "x2")
order coef
export delimited using "results/p12_profile_dgp.csv", replace
restore
