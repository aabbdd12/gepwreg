*! gepwreg_setable.ado  1.5.0  17sep2026  Araar A.
*! Post-estimation command of gepwreg: the consistent standard errors
*! available after the last estimation side by side -- analytical (influence
*! function), Taylor (svyset design) and bootstrap -- with the ratio of the two
*! most informative ones, as a check of the asymptotic approximation.  Kept in
*! its own file so that Stata finds it whether or not gepwreg.ado is in
*! memory.
#delimit ;

capture program drop gepwreg_setable ;
program define gepwreg_setable ;
    version 16 ;
    syntax [, NAIVE] ;
    if "`e(cmd)'" != "gepwreg" {;
        di as error "gepwreg_setable: last estimation must be gepwreg" ;
        error 301 ;
    } ;
    tempname b Vn Vi Vs Vb ;
    matrix `b' = e(b) ;
    local hasN = 0 ;
    local hasI = 0 ;
    local hasS = 0 ;
    local hasB = 0 ;
    /* an undefined e() matrix evaluates to a 1 x 1 missing, not an error:
       test the names in e(matrices) instead                                */
    local emats : e(matrices) ;
    /* e(V_naive), the WLS variance with fixed weights, is inconsistent and kept
       for the record; displayed only on request (naive), for the replication
       of the comparisons of the paper                                        */
    if "`naive'" != "" & `: list posof "V_naive" in emats' > 0 {; matrix `Vn' = e(V_naive) ; local hasN = 1 ; } ;
    if `: list posof "V_IF"    in emats' > 0 {; matrix `Vi' = e(V_IF)    ; local hasI = 1 ; } ;
    if `: list posof "V_svy"   in emats' > 0 {; matrix `Vs' = e(V_svy)   ; local hasS = 1 ; } ;
    if `: list posof "V_boot"  in emats' > 0 {; matrix `Vb' = e(V_boot)  ; local hasB = 1 ; } ;
    if `hasI' + `hasS' + `hasB' + `hasN' == 0 {;
        di as text "gepwreg_setable: no standard errors stored by the last estimation" ;
        di as text "(het(qr) without boot(): specify boot(#))." ;
        exit ;
    } ;
    /* the ratio: bootstrap against the best analytical one when both exist */
    if `hasB' & `hasS'      local rlab "svy/Boot" ;
    else if `hasB' & `hasI' local rlab "IF/Boot" ;
    else if `hasS' & `hasI' local rlab "svy/IF" ;
    else if `hasI' & `hasN' local rlab "IF/naive" ;
    else                    local rlab "" ;
    local method "`e(method)'" ;
    if "`method'" == "" local method "rankdep" ;
    di "" ;
    di as text "SE comparison  (" as result "`method'" as text ", tau=" as result %5.3f e(tau)
       as text "  h=" as result %8.6f e(h) as text ")" ;
    di as text "{hline 76}" ;
    local hdr : display %14s "Variable" %10s "Coeff." ;
    if `hasN' {; local h2 : display %10s "SE_naive" ; local hdr "`hdr'`h2'" ; } ;
    if `hasI' {; local h2 : display %10s "SE_IF"    ; local hdr "`hdr'`h2'" ; } ;
    if `hasS' {; local h2 : display %10s "SE_svy"   ; local hdr "`hdr'`h2'" ; } ;
    if `hasB' {; local h2 : display %10s "SE_boot"  ; local hdr "`hdr'`h2'" ; } ;
    if "`rlab'" != "" {; local h2 : display %10s "`rlab'" ; local hdr "`hdr'`h2'" ; } ;
    di as text "`hdr'" ;
    di as text "{hline 76}" ;
    local names : colnames e(b) ;
    local j = 0 ;
    foreach nm of local names {;
        local ++j ;
        if regexm("`nm'","[0-9]+b\.") | regexm("`nm'","[0-9]+o\.") continue ;
        local bj = `b'[1,`j'] ;
        local line : display %14s abbrev("`nm'",14) ;
        di as text "`line'" as result %10.5f `bj' _continue ;
        local sn = . ; local si = . ; local ss = . ; local sb = . ;
        if `hasN' {; local sn = sqrt(`Vn'[`j',`j']) ; di as result %10.5f `sn' _continue ; } ;
        if `hasI' {; local si = sqrt(`Vi'[`j',`j']) ; di as result %10.5f `si' _continue ; } ;
        if `hasS' {; local ss = sqrt(`Vs'[`j',`j']) ; di as result %10.5f `ss' _continue ; } ;
        if `hasB' {; local sb = sqrt(`Vb'[`j',`j']) ; di as result %10.5f `sb' _continue ; } ;
        if "`rlab'" == "svy/Boot"      local r = cond(`sb' > 0, `ss' / `sb', .) ;
        else if "`rlab'" == "IF/Boot"  local r = cond(`sb' > 0, `si' / `sb', .) ;
        else if "`rlab'" == "svy/IF"   local r = cond(`si' > 0, `ss' / `si', .) ;
        else if "`rlab'" == "IF/naive" local r = cond(`sn' > 0, `si' / `sn', .) ;
        if "`rlab'" != "" di as result %10.3f `r' ;
        else              di "" ;
    } ;
    di as text "{hline 76}" ;
    di as text "The table is a check of the asymptotic approximation, not a choice: the" ;
    di as text "main table already reports the appropriate standard error. A ratio far" ;
    di as text "from 1 (small N_eff, few PSUs per stratum, a rare category) says to read" ;
    di as text "the bootstrap." ;
    if `hasN' di as text "SE_naive : WLS with the kernel weights taken as fixed (inconsistent; shown on request)" ;
    if `hasI' {;
        if "`method'" == "twostep" | "`method'" == "initial"
            di as text "SE_IF    : influence function of the two-step (step 1, composition, ranks)" ;
        else
            di as text "SE_IF    : linearisation, kernel weights estimated (Deville 1999)" ;
    } ;
    if `hasS' di as text "SE_svy   : Taylor linearisation under the svyset design (PSU + strata)" ;
    if `hasB' {;
        di as text "SE_boot  : pairs bootstrap B=" as result e(boot)
           as text " (the whole procedure re-run at each draw)" ;
    } ;
    di "" ;
end ;
