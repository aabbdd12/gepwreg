* p12_core_mata.do -- used by p12_profile_dgp.do: copies of the Mata functions
* of gepwreg 1.5.0 that it calls (the percentile weights, the kernel, the
* interpolation of the grid of quantile regressions), renamed _gepwreg_ ->
* _rfc_, because the functions of an ado-file are private to it.  Generated
* from src/gepwreg.ado; not to be run on its own.
capture mata: mata drop _rfc_wp()
capture mata: mata drop _rfc_revCumSum()
capture mata: mata drop _rfc_revCumSumT()
capture mata: mata drop _rfc_kw()
capture mata: mata drop _rfc_Wmat()
capture mata: mata drop _rfc_interp()
mata:
real matrix _rfc_wp(real colvector y,
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

    st_numscalar("_rfc_h_tmp", h)
    return((w, pc))
}

real matrix _rfc_revCumSum(real matrix M)
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

real matrix _rfc_revCumSumT(real matrix M, real colvector p)
{
    real scalar n, i
    real matrix R
    R = _rfc_revCumSum(M)
    n = rows(M)
    for (i = 2; i <= n; i++) {
        if (p[i] == p[i-1]) R[i, .] = R[i-1, .]
    }
    return(R)
}

real colvector _rfc_kw(real colvector pc, real scalar tau,
                           real scalar h, real scalar n)
{
    return(exp(-0.25 :* ((pc :- tau) :/ h):^2) :/ (h * sqrt(2*pi()) * n))
}

real matrix _rfc_Wmat(real matrix X, real matrix Z)
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

real matrix _rfc_interp(real matrix G, real colvector u, real rowvector ug)
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

end
