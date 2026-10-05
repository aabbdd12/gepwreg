"""What choosing the bandwidth costs in variance, and what recovers it.

The influence-function variance of the two-step estimator is derived with the
bandwidth held fixed.  Since version 1.4 the bandwidth is chosen on the
sample, so the estimate inherits the variability of that choice, and the
influence function does not see it.  The Monte Carlo of p5 shows where it
matters: on DGP H1 at n = 4000 the ratio of the mean influence-function
standard error to the Monte Carlo standard deviation is 1.038 at tau = 0.1,
0.962 at tau = 0.5 and 0.835 at tau = 0.9, where the coverage of the nominal
95 percent interval falls to 0.824.  Under the Silverman rule, which is
effectively non-random because the rank is uniform, the same ratio was 1.00
at all three quantiles.

The mechanism is a first-order expansion: theta(h_hat) - theta(h_bar) is
approximately theta'(h_bar) (h_hat - h_bar), so the extra variance is
theta'(h)^2 Var(h_hat).  It is large exactly where the root mean squared
error curve is steep, which is where the plug-in is least reliable anyway:
near a boundary of the rank, with a curved profile.

This file measures the two pieces separately.  "sd at h*" re-chooses the
bandwidth on every sample, as the command does.  "sd at hbar" uses the SAME
average width, fixed in advance, so the only difference between the two
columns is the act of choosing.  It then checks whether the bootstrap, which
re-chooses the bandwidth inside every resample as gepwreg does, recovers the
right dispersion.

RESULT (H1, n = 4000, 500 seeds; bootstrap 150 seeds x 200 resamples).

    tau    sd at h*   sd at hbar   excess     boot/MC
    0.10    0.01685     0.01707     -1.3%      1.046
    0.50    0.01841     0.01771     +4.0%      0.994
    0.90    0.02023     0.01647    +22.8%      1.030

0.01647 / 0.02023 = 0.814 against the 0.835 observed in p5: the whole gap is
the bandwidth choice, and none of it is an error in the influence function.
The bootstrap covers it at all three quantiles.

CONSEQUENCE.  The plug-in for the point estimate, which is what
pwr_bandwidth_check.py justifies; boot() for the standard error wherever the
estimate is sensitive to the bandwidth -- near a boundary with a curved
profile, or whenever e(h) sits near the top of its clipping range; or
silverman, which restores the influence function by making the bandwidth
effectively deterministic, at the cost of a worse point estimate.

Run from python/:  python pwr_bandwidth_variance.py     (about two minutes)
"""
import numpy as np

from pwr_bandwidth_check import (TRUE, bhat, draw, nw, ranks, silverman,
                                 v_plugin)

N, R, B, RB = 4000, 500, 200, 150
TAUS = (0.1, 0.5, 0.9)


def fit(y, x, z):
    b = bhat(y, x, z)
    p = ranks(y)
    return b, p, silverman(p)


def hstar(b, p, tau, hs, n):
    w = np.exp(-0.25 * ((p - tau) / (2 * hs)) ** 2)
    S = w.sum()
    th0 = w @ b / S
    s2 = (w @ (b - th0) ** 2) / S
    return v_plugin(b, p, tau, hs, s2, th0, n, None)


def main():
    print(f"H1, n = {N}, {R} seeds.  Where does the extra variance come from?\n")
    print(f"{'tau':>5}{'true':>9}{'mean h*/hS':>12}{'sd(h*/hS)':>11}|"
          f"{'sd at h*':>10}{'sd at hbar':>12}{'sd at hSil':>12}|{'excess':>9}")
    mc = {}
    for tau in TAUS:
        TH, HR = [], []
        for s in range(R):
            rng = np.random.default_rng(2000 + s)
            y, x, z = draw("H1", rng, N)
            b, p, hs = fit(y, x, z)
            h = hstar(b, p, tau, hs, N)
            TH.append(nw(b, p, tau, h))
            HR.append(h / hs)
        hbar = float(np.mean(HR))
        TF, TS = [], []
        for s in range(R):
            rng = np.random.default_rng(2000 + s)
            y, x, z = draw("H1", rng, N)
            b, p, hs = fit(y, x, z)
            TF.append(nw(b, p, tau, hbar * hs))
            TS.append(nw(b, p, tau, hs))
        sa = float(np.std(TH, ddof=1))
        sf = float(np.std(TF, ddof=1))
        ss = float(np.std(TS, ddof=1))
        mc[tau] = sa
        print(f"{tau:5.2f}{TRUE['H1'][tau]:9.4f}{hbar:12.2f}"
              f"{np.std(HR, ddof=1):11.2f}|{sa:10.5f}{sf:12.5f}{ss:12.5f}|"
              f"{100 * (sa / sf - 1):+8.1f}%")
    print("\n  sd at h*   : the plug-in as used, bandwidth re-chosen on every"
          " sample")
    print("  sd at hbar : the SAME average width, fixed in advance -- no choice")
    print("  excess     : what choosing the bandwidth costs in variance\n")

    print(f"Does the bootstrap see it?  {RB} seeds x B = {B} pairs resamples\n")
    print(f"{'tau':>5}{'MC sd':>10}{'mean boot sd':>14}{'ratio':>8}")
    for tau in TAUS:
        BS = []
        for s in range(RB):
            rng = np.random.default_rng(2000 + s)
            y, x, z = draw("H1", rng, N)
            rb = np.random.default_rng(90000 + s)
            tb = np.empty(B)
            for k in range(B):
                i = rb.integers(0, N, N)
                bb, pb, hsb = fit(y[i], x[i], z[i])
                tb[k] = nw(bb, pb, tau, hstar(bb, pb, tau, hsb, N))
            BS.append(tb.std(ddof=1))
        print(f"{tau:5.2f}{mc[tau]:10.5f}{np.mean(BS):14.5f}"
              f"{np.mean(BS) / mc[tau]:8.3f}")


if __name__ == "__main__":
    main()
