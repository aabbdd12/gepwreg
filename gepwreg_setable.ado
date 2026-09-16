*! gepwreg_setable.ado  1.4.0  16sep2026  Araar A.
*! Post-estimation command of gepwreg: the naive, IF-corrected, Taylor (when
*! svyset) and bootstrap standard errors of the last gepwreg estimation side
*! by side.  Kept in its own file so that Stata finds it whether or not
*! gepwreg.ado is in memory.
#delimit ;

/* -----------------------------------------------------------------------------
   gepwreg_setable : side-by-side SE comparison table
   Columns: SE_naive (e(V_naive)), SE_IF (e(V_IF)), SE_svy (e(V_svy), when
   the design was used), SE_boot (e(V_boot), when boot() > 0).  The last
   column is the ratio of the two most informative ones: IF/Boot when a
   bootstrap was run, svy/IF under a survey design, IF/naive otherwise.
   ------------------------------------------------------------------------- */
capture program drop gepwreg_setable ;
program define gepwreg_setable ;
    version 16 ;
    if "`e(cmd)'" != "gepwreg" {;
        di as error "gepwreg_setable: last estimation must be gepwreg" ;
        error 301 ;
    } ;
    tempname b Vn Vi Vs Vb ;
    matrix `b'  = e(b) ;
    matrix `Vn' = e(V_naive) ;
    capture matrix `Vi' = e(V_IF) ;
    if _rc matrix `Vi' = e(V) ;
    local hasSvy  = (e(do_svy) == 1) ;
    local hasBoot = (e(boot) > 0) ;
    if `hasSvy'  matrix `Vs' = e(V_svy) ;
    if `hasBoot' matrix `Vb' = e(V_boot) ;
    local rlab = cond(`hasBoot', "IF/Boot", cond(`hasSvy', "svy/IF", "IF/naive")) ;
    di "" ;
    di as text "SE comparison  (tau=" as result %5.3f e(tau)
       as text "  h=" as result %8.6f e(h) as text ")" ;
    di as text "{hline 74}" ;
    local hdr : display %14s "Variable" %10s "Coeff." %10s "SE_naive" %10s "SE_IF" ;
    if `hasSvy' {;
        local h2 : display %10s "SE_svy" ;
        local hdr "`hdr'`h2'" ;
    } ;
    if `hasBoot' {;
        local h3 : display %10s "SE_boot" ;
        local hdr "`hdr'`h3'" ;
    } ;
    local h4 : display %9s "`rlab'" ;
    di as text "`hdr'`h4'" ;
    di as text "{hline 74}" ;
    local names : colnames e(b) ;
    local j = 0 ;
    foreach nm of local names {;
        local ++j ;
        /* base and omitted levels sit in e(b) as zeros: skip them */
        if regexm("`nm'","[0-9]+b\.") | regexm("`nm'","[0-9]+o\.") continue ;
        local bj = `b'[1,`j'] ;
        local sn = sqrt(`Vn'[`j',`j']) ;
        local si = sqrt(`Vi'[`j',`j']) ;
        local line : display %14s abbrev("`nm'",14) ;
        local line "`line'" ;
        di as text "`line'" as result %10.5f `bj' %10.5f `sn' %10.5f `si' _continue ;
        if `hasSvy' {;
            local ss = sqrt(`Vs'[`j',`j']) ;
            di as result %10.5f `ss' _continue ;
        } ;
        if `hasBoot' {;
            local sb = sqrt(`Vb'[`j',`j']) ;
            di as result %10.5f `sb' _continue ;
        } ;
        if `hasBoot'      local r = cond(`sb' > 0, `si' / `sb', .) ;
        else if `hasSvy'  local r = cond(`si' > 0, `ss' / `si', .) ;
        else              local r = cond(`sn' > 0, `si' / `sn', .) ;
        di as result %9.3f `r' ;
    } ;
    di as text "{hline 74}" ;
    di as text "SE_naive : WLS with the kernel weights taken as fixed (inconsistent)" ;
    di as text "SE_IF    : linearisation, kernel weights estimated (Deville 1999)" ;
    if `hasSvy' {;
        di as text "SE_svy   : Taylor linearisation under the svyset design (PSU + strata)" ;
    } ;
    if `hasBoot' {;
        di as text "SE_boot  : pairs bootstrap B=" as result e(boot)
           as text " (weights re-estimated at each draw)" ;
    } ;
    di "" ;
end ;
