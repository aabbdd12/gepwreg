* p11_profile_burkina.do -- Paper (version 5), Section 7.4, Table 8: the
* effect of the mean household under het(qr) along the distribution, on the
* Burkina Faso 1998 survey (bkf98I.dta): the effect for the group at tau,
* theta(tau); the linear correction to the population means,
* L = theta(tau) + pi(tau)'(Xbar_pop - Xbar(tau)), pi(tau) from the
* percentile-weighted regression of the household effects on X; and the exact
* effect E = beta(u*), u* the root of xbar'beta(u) = q_tau on the grid of
* quantile regressions interpolated linearly in u (the matching function of
* Alejo et al. 2025 evaluated for the mean household).
* Regressors: size, male head, urban, socio-economic group (i.gse);
* sampling weights; het(qr), grid 0.05(0.05)0.95; tau = 0.05(0.025)0.95.
* The household effects are rebuilt from the grid and their kernel mean is
* checked against gepwreg's theta(tau) (gap printed, of the order 1e-15).
* Output: results/p11_profile_burkina.csv (37 rows: tau, q_tau, u* and the
* three effects for size, male, urban and crop farmers); the table printed
* at tau = 0.1, 0.3, ..., 0.9 is Table 8.
* Run: do p11_profile_burkina   (about two minutes)

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
generate double f = weight
forvalues g = 2/7 {
    generate byte g`g' = (gse == `g')
}
local XV "size male urban g2 g3 g4 g5 g6 g7"

* ---- step 1 of het(qr): the grid of quantile regressions ---------------------
numlist "0.05(0.05)0.95"
local ug `r(numlist)'
local ng : word count `ug'
tempname G
matrix `G' = J(`ng', 10, .)
local g 0
foreach u of local ug {
    local ++g
    quietly qreg lexp `XV' [pw = f], quantile(`u')
    local j 0
    foreach v of local XV {
        local ++j
        matrix `G'[`g', `j'] = _b[`v']
    }
    matrix `G'[`g', 10] = _b[_cons]
}

* ---- the household effects as gepwreg builds them (for the linear correction)
mata:
    X  = st_data(., tokens("`XV'"))
    y  = st_data(., "lexp")
    G  = st_matrix("`G'")
    ug = strtoreal(tokens("`ug'"))
    n  = rows(y); k = cols(X); ng = cols(ug)
    Q  = (X, J(n, 1, 1)) * G'
    c  = rowsum(Q :< y)
    umid = (ug[1..(ng-1)] + ug[2..ng]) :/ 2
    ext  = (ug[1], umid, ug[ng])'
    u  = ext[c :+ 1]
    uc = rowmin((rowmax((u, J(n, 1, ug[1]))), J(n, 1, ug[ng])))
    Bi = J(n, k, .)
    for (g = 1; g < ng; g++) {
        if (g < ng - 1) sel = selectindex((uc :>= ug[g]) :& (uc :<  ug[g+1]))
        else            sel = selectindex((uc :>= ug[g]) :& (uc :<= ug[g+1]))
        if (rows(sel) > 0) {
            t = (uc[sel] :- ug[g]) :/ (ug[g+1] - ug[g])
            Bi[sel, .] = J(rows(sel), 1, G[g, 1..k]) + t * (G[g+1, 1..k] - G[g, 1..k])
        }
    }
    st_store(., st_addvar("double", ("q_size", "q_male", "q_urban", "q_crop")), Bi[., (1, 2, 3, 7)])
end

* ---- ranks as gepwreg; population means of X -----------------------------------
sort lexp
quietly generate double pc = sum(f)
quietly replace pc = pc / pc[_N]
quietly by lexp: replace pc = pc[_N]
tempname XB
matrix `XB' = J(1, 10, 1)
local j 0
foreach v of local XV {
    local ++j
    quietly summarize `v' [aw = f], meanonly
    local pop_`v' = r(mean)
    matrix `XB'[1, `j'] = r(mean)
}

* ---- the conditional quantile function of the mean household on the grid ------
mata:
    xb  = st_matrix("`XB'")
    Qb  = xb * G'
    printf("\nconditional quantiles of the mean household on the grid (must increase):\n")
    for (g = 1; g <= ng; g++) printf("  u = %4.2f  xbar'beta(u) = %8.4f\n", ug[g], Qb[g])
    printf("non-increasing steps: %f\n", sum(Qb[2..ng] :<= Qb[1..(ng-1)]))
end

* the exact effect: (us, ul, beta at us [size male urban crop], beta at ul [...]);
* us = gepwreg's count rule (midpoint), ul = the linear interpolation in u (E)
mata:
real rowvector _xact(real scalar qt, real rowvector Qb, real colvector ext,
                     real rowvector ug, real matrix G, real scalar k)
{
    real scalar ng, cnt, us, ul, g, m, uu, w
    real rowvector b, out
    ng  = cols(ug)
    // (a) gepwreg's rule: count of the grid predictions below, midpoint
    cnt = sum(Qb :< qt)
    us  = ext[cnt + 1]
    // (b) linear interpolation of xbar'beta(u_g) in u
    if (qt <= Qb[1])       ul = ug[1]
    else if (qt >= Qb[ng]) ul = ug[ng]
    else {
        g  = max(selectindex(Qb :<= qt))
        ul = ug[g] + (qt - Qb[g]) / (Qb[g+1] - Qb[g]) * (ug[g+1] - ug[g])
    }
    out = (us, ul)
    for (m = 1; m <= 2; m++) {
        uu  = (m == 1 ? us : ul)
        g   = max((1, min((ng - 1, max(selectindex(ug :<= uu + 1e-12))))))
        w   = (uu - ug[g]) / (ug[g+1] - ug[g])
        b   = G[g, 1..k] + w * (G[g+1, 1..k] - G[g, 1..k])
        out = out, b[(1, 2, 3, 7)]
    }
    return(out)
}
end

* ---- along the grid of tau -----------------------------------------------------
numlist "0.05(0.025)0.95"
local grid `r(numlist)'
local nt : word count `grid'
tempname R XR
matrix `R' = J(`nt', 21, .)
matrix colnames `R' = tau qtau u_step u_lin                          ///
    th_size th_male th_urban th_crop  lc_size lc_male lc_urban lc_crop ///
    es_size es_male es_urban es_crop  ex_size ex_male ex_urban ex_crop chk
local i 0
foreach t of local grid {
    local ++i
    matrix `R'[`i', 1] = `t'
    * q_tau as the weighted quantile of the ranks
    quietly summarize lexp if pc >= `t' - 1e-12, meanonly
    local qt = r(min)
    matrix `R'[`i', 2] = `qt'
    * gepwreg (het(qr)) and its bandwidth
    quietly gepwreg lexp size male urban i.gse, per(`t')
    local h = e(h)
    matrix `R'[`i', 5] = _b[size]
    matrix `R'[`i', 6] = _b[male]
    matrix `R'[`i', 7] = _b[urban]
    matrix `R'[`i', 8] = _b[5.gse]
    capture drop kw
    quietly generate double kw = f * exp(-0.25*((pc - `t')/`h')^2)
    local c 8
    local mx 0
    foreach e in size male urban crop {
        local ++c
        quietly summarize q_`e' [aw = kw], meanonly
        local mx = max(`mx', abs(r(mean) - `R'[`i', `c' - 4]))
        quietly regress q_`e' `XV' [aw = kw]
        local pred = _b[_cons]
        foreach v of local XV {
            local pred = `pred' + _b[`v']*`pop_`v''
        }
        matrix `R'[`i', `c'] = `pred'
    }
    matrix `R'[`i', 21] = `mx'
    * the exact effect of the mean household at q_tau
    mata: st_matrix("`XR'", _xact(`qt', Qb, ext, ug, G, k))
    matrix `R'[`i', 3] = `XR'[1, 1]
    matrix `R'[`i', 4] = `XR'[1, 2]
    forvalues c = 1/8 {
        matrix `R'[`i', 12 + `c'] = `XR'[1, 2 + `c']
    }
}

* ---- checks and the table --------------------------------------------------------
mata: R = st_matrix("`R'")
mata: printf("\nrebuilt het(qr) effects vs gepwreg, max |gap|: %9.2e\n", max(R[., 21]))
di as text _n "Table 8: theta = the group at tau; L = the linear correction; E = the exact effect"
mata: printf("\n  tau   q_tau  u*(lin) u*(step) | size: theta   lin   exact | male: theta   lin  exact | urban: theta  lin  exact | crop: theta   lin  exact\n")
mata: for (i = 3; i <= rows(R); i = i + 8) printf("%5.3f %7.3f %6.3f %6.3f | %8.4f %7.4f %7.4f | %6.3f %6.3f %6.3f | %6.3f %6.3f %6.3f | %6.3f %6.3f %6.3f\n", R[i,1], R[i,2], R[i,4], R[i,3], R[i,5], R[i,9], R[i,17], R[i,6], R[i,10], R[i,18], R[i,7], R[i,11], R[i,19], R[i,8], R[i,12], R[i,20])

* ---- save --------------------------------------------------------------------------
preserve
clear
quietly svmat double `R', names(col)
quietly summarize tau if u_lin >= 0.95 - 1e-9, meanonly
di as text "the mean household is held at u = 0.95 from tau = " as result r(min)
keep tau qtau u_lin th_* lc_* ex_*
export delimited using "results/p11_profile_burkina.csv", replace
restore
