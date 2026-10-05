*! gepwreg_examples 1.5.0  05oct2026  A. Araar (Universite Laval / PEP)
*! The examples of help gepwreg, run from their links.
*!   gepwreg_examples #          run example # in the command window
*!   gepwreg_examples #, db      open the dialog box of gepwreg filled in for example #
*!   gepwreg_examples #, db(k)   the same, with the k-th estimation of an example that
*!                               has several (examples 1 and 2: two; 7: four); db = db(1)
*!   gepwreg_examples #, do      open example # as a do-file in the Do-file Editor
*! The data in memory are never lost: the run keeps them (preserve) and gives
*! them back at the end, even after an error or a Break; the do-file does the
*! same; the dialog box, which needs the example data in memory, refuses to
*! replace data that have unsaved changes.
*! Examples 1-6 use bkf98I.dta, the Burkina Faso 1998 survey extract, an
*! ancillary file of the package, read from the current folder (where
*! "ssc install gepwreg, all" or "net get gepwreg" copies it, possibly as
*! bkf98i.dta), else from the SSC archive, else from GitHub (_gepwreg_exload);
*! nothing is written.  Example 7 simulates its data (set seed 20260919).
*! 1.5.0 (05oct2026): example 8, the profile in a figure
*! against the mean household (het(z), ref(mean), graph() with the reference
*! curve and the simultaneous band); example 9, the profile for a given
*! household under het(qr) with the bootstrap of the whole profile; db(k),
*! one link of the dialog box per estimation of an example.
program define gepwreg_examples
    version 16
    syntax anything(name=ex id="example number") [, DO NOEDIT *]
    capture confirm integer number `ex'
    if _rc | !inrange(`ex', 1, 9) {
        di as err "gepwreg_examples: the examples are numbered 1 to 9 (see help gepwreg)"
        exit 198
    }
    * db or db(k): the dialog box, filled in with the k-th estimation (db = db(1))
    local db
    local dbk 1
    if `"`options'"' != "" {
        if !ustrregexm(`"`options'"', "^\s*db\s*(\(\s*([0-9]+)\s*\))?\s*$") {
            di as err `"gepwreg_examples: option `options' not allowed (db, db(#) or do)"'
            exit 198
        }
        local db db
        if ustrregexs(2) != "" local dbk = real(ustrregexs(2))
    }
    if "`db'" != "" & "`do'" != "" {
        di as err "gepwreg_examples: db and do cannot be combined"
        exit 198
    }
    * the estimations an example offers to the dialog box
    local ndb 1
    if inlist(`ex', 1, 2) local ndb 2
    if `ex' == 7          local ndb 4
    if "`db'" != "" & !inrange(`dbk', 1, `ndb') {
        if `ndb' == 1 di as err "gepwreg_examples: example `ex' has one estimation for the dialog box: db (or db(1))"
        else          di as err "gepwreg_examples: example `ex' has `ndb' estimations for the dialog box: db(1) to db(`ndb')"
        exit 198
    }
    local X "lexp size male i.gse"
    local RIF "rifhdreg lexp size male i.gse [pw=weight], rif(q(25))"
    local data "bkf98I"
    * the setup of the survey data (examples 1-6)
    local s1 "generate lexp  = ln(exppc)"
    local s2 "generate male  = (sex == 1)"
    local s3 "generate urban = (zone == 2)"
    local ns 3
    local n 0
    local rif 0
    local loop 0
    if `ex' == 1 {
        local title "The effect of the covariates for the households at the first quartile of expenditure"
        local c1 "gepwreg `X', per(0.25) het(urban size)"
        local c2 "gepwreg `X', per(0.25) boot(50)"
        local n 2
        local rif 1
    }
    else if `ex' == 2 {
        local title "Survey design: Taylor standard errors, the bootstrap of the design, the three side by side"
        local c1 "svyset psu [pw=weight], strata(strata)"
        local c2 "gepwreg `X', per(0.25) het(urban size)"
        local c3 "gepwreg `X', per(0.25) het(urban size) boot(200)"
        local c4 "gepwreg_setable"
        local n 4
    }
    else if `ex' == 3 {
        local title "The profile of the effects across percentiles"
        local loop 1
        local n 0
        foreach t in 0.10 0.25 0.50 0.75 0.90 {
            local q = round(100*`t')
            local ++n
            local c`n' "quietly gepwreg `X', per(`t') het(urban size)"
            local ++n
            local c`n' "estimates store q`q'"
        }
        local ++n
        local c`n' "estimates table q10 q25 q50 q75 q90, se"
    }
    else if `ex' == 4 {
        local title "The effect at the initial quantile: the households poor before the urban premium"
        local c1 "gepwreg lexp urban size male i.gse, per(0.25) het(size) initial(urban)"
        local n 1
    }
    else if `ex' == 5 {
        local title "The profile of the effects along household size"
        local c1 "gepwreg lexp male urban i.gse, per(0.25) rankvar(size)"
        local n 1
    }
    else if `ex' == 6 {
        local title "The one-step regression of version 1.3, for comparison (descriptive only)"
        local c1 "gepwreg `X', per(0.25) rankdep"
        local n 1
    }
    else if `ex' == 8 {
        local title "The profile across percentiles against the mean household: composition (het(z))"
        local c1 "gepwreg lexp size male urban i.gse, per(0.25) het(urban size) ref(mean) graph(size male 5.gse urban, ref uniform)"
        local n 1
    }
    else if `ex' == 9 {
        local title "The profile for a given household: is the effect constant along the distribution? (het(qr))"
        local c1 "gepwreg lexp size male urban i.gse, per(0.25) ref(default) boot(50) graph(size male 5.gse urban, ref uniform)"
        local n 1
    }
    else if `ex' == 7 {
        local title "Measurement error in a heterogeneity variable: the correction against a known truth"
        local data ""
        local s1 "clear"
        local s2 "set seed 20260919"
        local s3 "set obs 8000"
        local s4 "generate double g   = rnormal()"
        local s5 "generate double z1  = exp(0.45*g)"
        local s6 "generate double z2  = runiform()"
        local s7 "generate double x   = 1 + 0.45*g + sqrt(1-0.45^2)*rnormal()"
        local s8 "generate double b   = 0.5 + 1.5*z1 + 1.0*z2"
        local s9 "generate double y   = 5 + b*x + 3*z1 + 2*z2 + rnormal()"
        local s10 "generate double zt1 = z1 + 0.50*rnormal()"
        local s11 "generate double zt2 = z2 + 0.25*rnormal()"
        local ns 11
        local c1 "gepwreg y x, het(z1 z2)   per(0.75)"
        local c2 "gepwreg y x, het(zt1 zt2) per(0.75)"
        local c3 "gepwreg y x, het(zt1 zt2) per(0.75) merr"
        local c4 "gepwreg y x, het(zt1 c.zt1#c.zt1 zt2) per(0.75)"
        local n 4
    }

    * ---- as a do-file, in Stata's temporary folder ----
    if "`do'" != "" {
        local T = c(tmpdir)
        local T : subinstr local T "\" "/", all
        if substr("`T'", -1, 1) != "/" local T "`T'/"
        local fn "`T'gepwreg_example_`ex'.do"
        tempname fh
        file open `fh' using "`fn'", write text replace
        file write `fh' "* gepwreg, example `ex': `title'" _n
        file write `fh' "* Written by gepwreg_examples in Stata's temporary folder; save it elsewhere to keep it." _n
        file write `fh' "* preserve keeps the data in memory and gives them back when this do-file ends;" _n
        file write `fh' "* delete that line to keep working on the example data." _n
        file write `fh' "preserve" _n
        if "`data'" != "" {
            file write `fh' "* the example data: the current folder (where ssc install gepwreg, all copies" _n
            file write `fh' "* them), else the SSC archive, else GitHub" _n
            file write `fh' "capture use bkf98I, clear" _n
            file write `fh' "if _rc capture use bkf98i, clear" _n
            file write `fh' `"if _rc capture use "http://fmwww.bc.edu/repec/bocode/b/bkf98i.dta", clear"' _n
            file write `fh' `"if _rc capture use "http://fmwww.bc.edu/repec/bocode/b/bkf98I.dta", clear"' _n
            file write `fh' `"if _rc use "https://raw.githubusercontent.com/aabbdd12/gepwreg/main/examples/bkf98I.dta", clear"' _n
        }
        forvalues i = 1/`ns' {
            file write `fh' `"`s`i''"' _n
        }
        if `loop' {
            file write `fh' "foreach tau in 0.10 0.25 0.50 0.75 0.90 {" _n
            file write `fh' "    gepwreg `X', per(\`tau') het(urban size)" _n
            file write `fh' "    estimates store q\`=round(100*\`tau')'" _n
            file write `fh' "}" _n
            file write `fh' "estimates table q10 q25 q50 q75 q90, se" _n
        }
        else {
            forvalues i = 1/`n' {
                file write `fh' `"`c`i''"' _n
            }
        }
        if `rif' {
            file write `fh' "* RIF regression for comparison: rifhdreg is part of the rif package (ssc install rif)" _n
            file write `fh' `"capture noisily `RIF'"' _n
        }
        file close `fh'
        if "`noedit'" == "" doedit "`fn'"
        di as txt "(example `ex' written to " as res `"`fn'"' as txt ")"
        exit
    }

    * ---- in the dialog box: needs the example data in memory ----
    if "`db'" != "" {
        * the example data loaded for a dialog box carry a mark (_dta[gepwreg_example],
        * or the mark of the examples of another package): they are replaced without
        * asking; the user's own data with changes are not
        local isex 0
        local cl : char _dta[]
        foreach c of local cl {
            if substr("`c'", -8, .) == "_example" & `"`: char _dta[`c']'"' == "1" local isex 1
        }
        if c(changed) & !`isex' {
            di as err "gepwreg_examples, db: the data in memory have changes not saved;"
            di as err "save them (or clear) first: the dialog box needs the example data in memory"
            exit 4
        }
        if "`data'" != "" _gepwreg_exload
        forvalues i = 1/`ns' {
            quietly `s`i''
        }
        char _dta[gepwreg_example] "1"
        if "`data'" != "" di as txt "(example data `data' loaded, with lexp, male and urban, for the dialog box)"
        else              di as txt "(simulated data of example 7 in memory for the dialog box)"
        if `ex' == 3 di as txt "(the dialog box fits one percentile; the profile is the loop of the do-file: gepwreg_examples 3, do)"
        if `ex' == 1 & `dbk' == 2 di as txt "(het(qr) with boot(50): a little over a minute)"
        if `ex' == 9 di as txt "(example 9 runs in under two minutes: 50 bootstrap draws of the whole procedure, one loop for every percentile)"
        db gepwreg
        * Stata keeps the state of the dialog between two openings: every control
        * an example may set is first put back to its default
        .gepwreg_dlg.main.vn_dv.setvalue ""
        .gepwreg_dlg.main.vl_iv.setvalue ""
        .gepwreg_dlg.main.ed_per.setvalue "0.5"
        .gepwreg_dlg.main.ck_nc.setoff
        .gepwreg_dlg.main.ed_in.setvalue ""
        .gepwreg_dlg.main.ed_xr.setvalue ""
        .gepwreg_dlg.het.vl_hv.setvalue ""
        .gepwreg_dlg.het.ed_qg.setvalue ""
        .gepwreg_dlg.het.vn_rv.setvalue ""
        .gepwreg_dlg.het.ck_me.setoff
        .gepwreg_dlg.het.ed_tc.setvalue "2"
        .gepwreg_dlg.bw.rb_opt.seton
        .gepwreg_dlg.bw.ed_bv.setvalue ""
        .gepwreg_dlg.se.ed_bt.setvalue "0"
        .gepwreg_dlg.se.ed_sd.setvalue "12345"
        .gepwreg_dlg.wgt.rb_wnone.seton
        .gepwreg_dlg.wgt.ed_wexp.setvalue ""
        .gepwreg_dlg.rep.ck_st.setoff
        .gepwreg_dlg.rep.ed_gv.setvalue ""
        .gepwreg_dlg.rep.ed_gg.setvalue ""
        .gepwreg_dlg.rep.ck_gu.setoff
        .gepwreg_dlg.rep.ck_gr.setoff
        .gepwreg_dlg.rep.rb_rn.seton
        .gepwreg_dlg.rep.ed_rv.setvalue ""
        * the example
        if `ex' == 7 {
            .gepwreg_dlg.main.vn_dv.setvalue "y"
            .gepwreg_dlg.main.vl_iv.setvalue "x"
            .gepwreg_dlg.main.ed_per.setvalue "0.75"
        }
        else {
            .gepwreg_dlg.main.vn_dv.setvalue "lexp"
            .gepwreg_dlg.main.ed_per.setvalue "0.25"
            if `ex' == 4      .gepwreg_dlg.main.vl_iv.setvalue "urban size male i.gse"
            else if inlist(`ex', 8, 9) .gepwreg_dlg.main.vl_iv.setvalue "size male urban i.gse"
            else if `ex' == 5 .gepwreg_dlg.main.vl_iv.setvalue "male urban i.gse"
            else              .gepwreg_dlg.main.vl_iv.setvalue "size male i.gse"
        }
        local qr = (`ex' == 9 | (`ex' == 1 & `dbk' == 2))
        if `ex' == 5      _gepwreg_exdbmode rv
        else if `ex' == 6 _gepwreg_exdbmode rd
        else if `qr'      _gepwreg_exdbmode qr
        else              _gepwreg_exdbmode z
        if inlist(`ex', 1, 2, 3, 8) & !`qr' .gepwreg_dlg.het.vl_hv.setvalue "urban size"
        if `ex' == 1 & `dbk' == 2 .gepwreg_dlg.se.ed_bt.setvalue "50"
        if inlist(`ex', 8, 9) {
            .gepwreg_dlg.rep.ed_gv.setvalue "size male 5.gse urban"
            .gepwreg_dlg.rep.ck_gu.seton
            .gepwreg_dlg.rep.ck_gr.seton
        }
        if `ex' == 8 .gepwreg_dlg.rep.rb_rm.seton
        if `ex' == 9 {
            .gepwreg_dlg.rep.rb_rd.seton
            .gepwreg_dlg.se.ed_bt.setvalue "50"
        }
        if `ex' == 4 {
            .gepwreg_dlg.het.vl_hv.setvalue "size"
            .gepwreg_dlg.main.ed_in.setvalue "urban"
        }
        if `ex' == 5 .gepwreg_dlg.het.vn_rv.setvalue "size"
        if `ex' == 2 & `dbk' == 2 {
            .gepwreg_dlg.se.ed_bt.setvalue "200"
            .gepwreg_dlg.rep.ck_st.seton
        }
        * example 7: the oracle, the naive fit, the correction, the enriched model
        if `ex' == 7 {
            if `dbk' == 1      .gepwreg_dlg.het.vl_hv.setvalue "z1 z2"
            else if `dbk' == 4 .gepwreg_dlg.het.vl_hv.setvalue "zt1 c.zt1#c.zt1 zt2"
            else               .gepwreg_dlg.het.vl_hv.setvalue "zt1 zt2"
            if `dbk' == 3 {
                .gepwreg_dlg.het.ck_me.seton
                .gepwreg_dlg.het.tx_tc.enable
                .gepwreg_dlg.het.ed_tc.enable
            }
        }
        exit
    }

    * ---- in the command window: the data in memory are kept ----
    if `ex' == 9 di as txt "(under two minutes: 50 bootstrap draws of the whole procedure, one loop for every percentile)"
    preserve
    di as txt _n "{hline 78}" _n "gepwreg, example `ex': " as res "`title'" _n as txt "{hline 78}"
    if "`data'" != "" {
        _gepwreg_exload
        di as txt "(example data `data', read from `r(source)'; the data in memory come back at the end)"
    }
    else di as txt "(simulated data; the data in memory come back at the end)"
    forvalues i = 1/`ns' {
        di as txt `". `s`i''"'
        quietly `s`i''
    }
    forvalues i = 1/`n' {
        di as txt _n `". `c`i''"'
        `c`i''
    }
    if `rif' {
        di as txt _n `". `RIF'"'
        capture which rifhdreg
        if _rc {
            di as txt "(rifhdreg is part of the rif package, which is not installed: " _c
            di as txt `"{stata "ssc install rif":ssc install rif}"' as txt ")"
        }
        else `RIF'
    }
end

* ============================================================================
* the example data: bkf98I.dta, an ancillary file of the package, from the
* current folder (where "ssc install gepwreg, all" or "net get gepwreg" copies
* it; net get may write it as bkf98i.dta), else beside the command or in the
* examples/ folder of a copy of the repository (src/../examples), else from
* the SSC archive, else from GitHub; nothing is written to disk
program define _gepwreg_exload, rclass
    foreach f in bkf98I bkf98i {
        capture confirm file "`f'.dta"
        if !_rc {
            quietly use "`f'.dta", clear
            return local source "the current folder"
            exit
        }
    }
    capture findfile gepwreg_examples.ado
    if !_rc {
        local here = subinstr("`r(fn)'", "\", "/", .)
        local here = substr("`here'", 1, strrpos("`here'", "/") - 1)
        foreach d in "`here'" "`here'/../examples" {
            foreach f in bkf98I bkf98i {
                capture confirm file "`d'/`f'.dta"
                if !_rc {
                    quietly use "`d'/`f'.dta", clear
                    return local source "`d'"
                    exit
                }
            }
        }
    }
    foreach f in bkf98i bkf98I {
        capture quietly use "http://fmwww.bc.edu/repec/bocode/b/`f'.dta", clear
        if !_rc {
            return local source "the SSC archive"
            exit
        }
    }
    capture quietly use "https://raw.githubusercontent.com/aabbdd12/gepwreg/main/examples/bkf98I.dta", clear
    if !_rc {
        return local source "GitHub"
        exit
    }
    di as err "gepwreg_examples: bkf98I.dta is not in the current folder (`c(pwd)'),"
    di as err "  and it could not be read from the SSC archive or from GitHub."
    di as txt "  Copy the example data into the current folder with"
    di as txt `"  {stata "ssc install gepwreg, all replace"} (from SSC), or"'
    di as txt `"  {stata "net get gepwreg, from(https://raw.githubusercontent.com/aabbdd12/gepwreg/main)"} (from GitHub)."'
    exit 601
end

* ============================================================================
* the dialog box: the estimator chosen, and the controls that go with it
* enabled or disabled explicitly (as the dialog's own script het_ck does), in
* case the radio button does not run its script when set from here
program define _gepwreg_exdbmode
    args m
    if "`m'" == "z"  .gepwreg_dlg.het.rb_z.seton
    if "`m'" == "qr" .gepwreg_dlg.het.rb_qr.seton
    if "`m'" == "rv" .gepwreg_dlg.het.rb_rv.seton
    if "`m'" == "rd" .gepwreg_dlg.het.rb_rd.seton
    local hv  = cond("`m'" == "z",  "enable", "disable")
    local qg  = cond("`m'" == "qr", "enable", "disable")
    local rv  = cond("`m'" == "rv", "enable", "disable")
    local in  = cond(inlist("`m'", "z", "qr"), "enable", "disable")
    .gepwreg_dlg.het.tx_hv.`hv'
    .gepwreg_dlg.het.vl_hv.`hv'
    .gepwreg_dlg.het.ck_me.`hv'
    .gepwreg_dlg.het.tx_qg.`qg'
    .gepwreg_dlg.het.ed_qg.`qg'
    .gepwreg_dlg.het.tx_rv.`rv'
    .gepwreg_dlg.het.vn_rv.`rv'
    .gepwreg_dlg.het.tx_tc.disable
    .gepwreg_dlg.het.ed_tc.disable
    .gepwreg_dlg.main.tx_in.`in'
    .gepwreg_dlg.main.ed_in.`in'
    .gepwreg_dlg.main.tx_xr.`in'
    .gepwreg_dlg.main.ed_xr.`in'
end
