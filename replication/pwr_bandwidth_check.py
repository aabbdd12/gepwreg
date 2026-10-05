"""Does the MSE-optimal bandwidth of gepwreg 1.4 beat the Silverman rule?

This reproduces the two-step \\texttt{het(z)} estimator and the plug-in
bandwidth of version 1.4 in numpy, on the DGPs of
p4_twostep_vs_rif.do at three sample sizes, against the EXACT
theta(tau) of pwr_exact_truth.py and against the ORACLE bandwidth -- the one
that minimises the realised root mean squared error, which no rule can beat.

HOW TO READ IT, AND HOW NOT TO.  The natural summary, the mean RMSE over the
cells, is the wrong one: the cells differ in scale by a factor of fifty, and
averaging their levels lets one DGP decide the answer.  (On a first pass it
said the rule was worth one percent -- it is not.)  The scale-free summary is
each cell's RMSE divided by the oracle's, which is the same argument that
produced the aggregate rule itself: the coefficients do not share units, so
what is minimised is the average RELATIVE mean squared error.

U1 is excluded here.  Its heterogeneity variable is irrelevant by
construction, the bias of the estimator is 0.2 to 0.3 at every bandwidth, and
no bandwidth rule can be judged on it.

RESULT (45 cells: 3 DGPs x 5 quantiles x n = 2000, 6000, 20000; 300 seeds).
RMSE relative to the oracle bandwidth:

    Silverman          mean 1.326   median 1.096   worst 3.161   22/45 > 10%
    plug-in 1.4      mean 1.131   median 1.063   worst 1.524   17/45 > 10%
    direct evaluation  mean 1.090   median 1.064   worst 1.458   11/45 > 10%
    shifted stencil    mean 1.431   median 1.098   worst 3.797   22/45 > 10%
    one-sided stencil  mean 1.153   median 1.058   worst 1.929   18/45 > 10%

The plug-in cuts the average distance to the oracle from 33 to 13 percent and
the worst case from 216 to 52 percent.  Where Silverman fails it fails badly:
D1 at tau = 0.1 and 0.9 is a steep, nearly straight profile at the boundary,
where the local constant smoother at the Silverman bandwidth is almost all
bias (0.0087 of an RMSE of 0.0092 at n = 6000) and the plug-in halves the
bandwidth and removes it.

THE ONE FAILURE, AND TWO ATTEMPTS THAT DID NOT REPAIR IT.  The plug-in is
worse than Silverman in one recurring cell, H1 at tau = 0.9, at all three
sample sizes (1.47, 1.40, 1.46 against the oracle, where Silverman is 1.04,
1.16, 1.38).  The cause is identified: the pilot second difference evaluates
theta at tau + 0.05 = 0.95, where the pilot kernel -- standard deviation
2 h_Sil sqrt(2) = 0.13 at n = 6000 -- is truncated by the support of the
rank, so the estimated curvature is too small and h* comes out too wide.
Two repairs were tried and BOTH ARE WORSE than doing nothing: recentring the
stencil inward destroys the D1 boundary cells the rule handles best (1.02 ->
1.99), and a one-sided second difference fixes H1 at tau = 0.9 but breaks D3
at tau = 0.75 (1.23 -> 1.93).  They are kept in the file as a negative
result.  The variant that does help is `direct', which computes the bias by
evaluating the actual kernel on the actual sample instead of through the
asymptotic expansion, so the truncation needs no theory; it is better on
average and on the count of bad cells, and has an erratic failure of its own
where the pilot profile is locally flat (H1, tau = 0.75, n = 20000): it sees
no bias, so it takes the widest bandwidth on offer and the true bias then
bites.  It is a candidate for a later version, not a reason to hold this one.

Run from python/:  python pwr_bandwidth_check.py     (about five minutes)
"""
import numpy as np

R = 300
NS = (2000, 6000, 20000)
DGPS = ("D3", "H1", "D1")
TAUS = (0.1, 0.25, 0.5, 0.75, 0.9)
MULTS = np.array([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.25,
                  1.4, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0])
PROFILE_GRID = np.linspace(0.005, 0.995, 100)
SQ2 = np.sqrt(2.0)

# exact theta(tau), from pwr_exact_truth.py
TRUE = {
    "D3": {0.1: 1.231125, 0.25: 1.319275, 0.5: 1.492782,
           0.75: 1.672133, 0.9: 1.780590},
    "H1": {0.1: 1.464797, 0.25: 1.447873, 0.5: 1.462897,
           0.75: 1.513874, 0.9: 1.588059},
    "D1": {0.1: 1.100110, 0.25: 1.250125, 0.5: 1.500150,
           0.75: 1.750175, 0.9: 1.900190},
    "U1": {0.1: 1.182857, 0.25: 1.289697, 0.5: 1.532903,
           0.75: 1.683244, 0.9: 1.799876},
}


def draw(dgp, rng, n):
    if dgp == "D3":
        z = rng.random(n); x = rng.standard_normal(n)
        y = 5 * z + (1 + z) * x
    elif dgp == "H1":
        z = rng.random(n); x = rng.standard_normal(n)
        y = 5 + (1 + z) * x + z + rng.standard_normal(n)
    elif dgp == "D1":
        z = rng.random(n); x = rng.standard_normal(n)
        y = 100 * z + (1 + z) * x
    else:                                     # U1: z is irrelevant
        u = rng.random(n); x = 2 * rng.random(n)
        z = rng.standard_normal(n)
        y = 5 + (1 + u) * x + 3 * u
    return y, x, z


def bhat(y, x, z):
    """Step 1: b_i = gamma + d z_i from the fit of y on [1, x, z, x z].
    It does not depend on the bandwidth, which is what makes the plug-in
    cost three weighted means instead of a second fit."""
    W = np.column_stack([np.ones(len(y)), x, z, x * z])
    g = np.linalg.lstsq(W, y, rcond=None)[0]
    return g[1] + g[3] * z


def ranks(y):
    n = len(y)
    o = np.argsort(y, kind="stable")
    p = np.empty(n)
    p[o] = (np.arange(n) + 1.0) / n
    return p


def silverman(p):
    iqr = np.subtract(*np.quantile(p, [0.75, 0.25]))
    return 0.9 * min(p.std(ddof=1), iqr / 1.34) * len(p) ** -0.2


def nw(b, p, tau, h):
    """Step 2: the local constant average of b on the rank of y at tau."""
    w = np.exp(-0.25 * ((p - tau) / h) ** 2)
    return float(w @ b / w.sum())


def _h(s2, d2, n, hsil):
    if s2 <= 1e-12 or d2 == 0:
        return 3 * hsil
    h = (s2 / (8 * np.sqrt(2 * np.pi) * n * d2 ** 2)) ** 0.2
    return float(np.clip(h, 0.5 * hsil, 3 * hsil))


def v_plugin(b, p, tau, hsil, s2, th0, n, prof):
    """gepwreg 1.4: pilot at 2 h_Sil, theta'' by a centred second
    difference, h* = [s^2 / (8 sqrt(2 pi) n theta''^2)]^(1/5)."""
    hp = 2 * hsil
    d = min(0.05, tau / 2, (1 - tau) / 2)
    d2 = (nw(b, p, tau + d, hp) - 2 * th0 + nw(b, p, tau - d, hp)) / d ** 2
    return _h(s2, d2, n, hsil)


def v_shifted(b, p, tau, hsil, s2, th0, n, prof, c=2.0):
    """NEGATIVE RESULT: recentre the stencil so that tau' +- d stays c pilot
    standard deviations inside [0,1].  Worse than doing nothing."""
    hp = 2 * hsil
    d = min(0.05, tau / 2, (1 - tau) / 2)
    m = c * SQ2 * hp + d
    t = float(np.clip(tau, min(m, 0.5), max(1 - m, 0.5)))
    a = nw(b, p, t, hp) if t != tau else th0
    d2 = (nw(b, p, t + d, hp) - 2 * a + nw(b, p, t - d, hp)) / d ** 2
    return _h(s2, d2, n, hsil)


def v_onesided(b, p, tau, hsil, s2, th0, n, prof, c=2.0):
    """NEGATIVE RESULT: one-sided second difference when the outward point is
    within c pilot standard deviations of the edge.  Also worse."""
    hp = 2 * hsil
    d = min(0.05, tau / 2, (1 - tau) / 2)
    m = c * SQ2 * hp + d
    if tau > 1 - m and tau >= 0.5:
        d2 = (th0 - 2 * nw(b, p, tau - d, hp) + nw(b, p, tau - 2 * d, hp)) / d ** 2
    elif tau < m and tau < 0.5:
        d2 = (th0 - 2 * nw(b, p, tau + d, hp) + nw(b, p, tau + 2 * d, hp)) / d ** 2
    else:
        d2 = (nw(b, p, tau + d, hp) - 2 * th0 + nw(b, p, tau - d, hp)) / d ** 2
    return _h(s2, d2, n, hsil)


def v_direct(b, p, tau, hsil, s2, th0, n, prof):
    """Minimise the mean squared error by DIRECT evaluation: the bias of the
    local constant average is the kernel average of a pilot profile minus the
    pilot at tau, computed with the actual kernel on the actual sample, so
    the truncation at the boundary is handled without an expansion."""
    tprof = np.interp(p, PROFILE_GRID, prof)
    t_tau = float(np.interp(tau, PROFILE_GRID, prof))
    best, bh = np.inf, hsil
    for m in MULTS:
        w = np.exp(-0.25 * ((p - tau) / (m * hsil)) ** 2)
        S = w.sum()
        mse = ((w @ tprof) / S - t_tau) ** 2 + s2 * (w @ w) / S ** 2
        if mse < best:
            best, bh = mse, m * hsil
    return bh


RULES = [("plug-in 1.4", v_plugin), ("direct evaluation", v_direct),
         ("shifted stencil", v_shifted), ("one-sided stencil", v_onesided)]


def main():
    summ = {k: [] for k, _ in RULES}
    summ["Silverman"] = []
    head = f"{'n':>6}{'DGP':>4}{'tau':>6}|{'oracle':>8}{'@':>6}{'Silv':>7}|"
    print(head + "".join(f"{k[:9]:>22}" for k, _ in RULES))
    for n in NS:
        for dgp in DGPS:
            G = {(t, m): [] for t in TAUS for m in MULTS}
            V = {(t, k): [] for t in TAUS for k, _ in RULES}
            RA = {(t, k): [] for t in TAUS for k, _ in RULES}
            for s in range(R):
                rng = np.random.default_rng(1000 + s)
                y, x, z = draw(dgp, rng, n)
                b = bhat(y, x, z)
                p = ranks(y)
                hs = silverman(p)
                dd = (p[None, :] - PROFILE_GRID[:, None]) / (2 * hs)
                K = np.exp(-0.25 * dd * dd)
                prof = (K @ b) / K.sum(1)
                for t in TAUS:
                    w = np.exp(-0.25 * ((p - t) / (2 * hs)) ** 2)
                    S = w.sum()
                    th0 = w @ b / S
                    s2 = (w @ (b - th0) ** 2) / S
                    for m in MULTS:
                        G[(t, m)].append(nw(b, p, t, m * hs))
                    for k, f in RULES:
                        h = f(b, p, t, hs, s2, th0, n, prof)
                        V[(t, k)].append(nw(b, p, t, h))
                        RA[(t, k)].append(h / hs)
            for t in TAUS:
                tt = TRUE[dgp][t]

                def rmse(v):
                    return float(np.sqrt(((np.array(v) - tt) ** 2).mean()))

                rg = {m: rmse(G[(t, m)]) for m in MULTS}
                mb = min(rg, key=rg.get)
                rb = rg[mb]
                summ["Silverman"].append(rg[1.0] / rb)
                line = (f"{n:>6}{dgp:>4}{t:6.2f}|{rb:8.4f}{mb:6.2f}"
                        f"{rg[1.0] / rb:7.3f}|")
                for k, _ in RULES:
                    rr = rmse(V[(t, k)])
                    summ[k].append(rr / rb)
                    line += f"{rr:10.4f}{rr / rb:6.3f}{np.mean(RA[(t, k)]):6.2f}"
                print(line)
            print()
    ncell = len(summ["Silverman"])
    print(f"RMSE relative to the oracle bandwidth, {ncell} cells")
    for k in ["Silverman"] + [x for x, _ in RULES]:
        a = np.array(summ[k])
        print(f"  {k:<20} mean {a.mean():5.3f}  median {np.median(a):5.3f}"
              f"  worst {a.max():5.3f}  > 10% off {(a > 1.10).sum():2d}/{ncell}")


if __name__ == "__main__":
    main()
