*! gepwreg_setable.ado  1.3.1  16sep2026  Araar A.
*! Post-estimation command of gepwreg: the naive, IF-corrected and bootstrap
*! standard errors of the last gepwreg estimation side by side.  Kept in its
*! own file so that Stata finds it whether or not gepwreg.ado is in memory.
#delimit ;

/* -----------------------------------------------------------------------------
   gepwreg_setable : side-by-side SE comparison table
   ------------------------------------------------------------------------- */
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
