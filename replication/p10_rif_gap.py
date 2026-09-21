"""Where RIF regression and theta(tau) coincide, and where they part.

Section 4.6 of the paper.  Writes results/p10_rif_gap.csv.

THE DESIGN.  Take the simplest model that exists,

    y = b x + u,     u independent of x,     b constant,

with no heterogeneity, no endogeneity and nothing misspecified.  Shifting x
by delta shifts every quantile of y by b*delta, so the unconditional quantile
partial effect and theta(tau) = E[b | y = q_tau] both equal b at every tau,
exactly and with no approximation of any kind.  Whatever a RIF regression
reports other than b is therefore its own projection error, isolated from
every other question.

THE TWO QUANTITIES.  The population slope of a RIF regression is a linear
projection of the conditional cdf at the quantile,

    gamma_RIF(tau) = - Cov( F(q_tau | x), x ) / ( f(q_tau) Var(x) ),

while, since dF(q|x)/dx = -f(q|x) b, the effect itself is an average
derivative of the same function,

    theta(tau) = E[ -dF(q_tau|X)/dx ] / f(q_tau) = b.

So the gap between them is exactly the gap between a linear projection of
F(q_tau|.) and its average derivative, and nothing else.

TWO SUFFICIENT CONDITIONS FOR EQUALITY.  The projection equals the average
derivative when F(q_tau|.) is affine on the support of x -- which needs the
cdf of u to be affine over the whole range q_tau - b*supp(x), so u uniform
AND no part of the support running off its edge, a knife-edge -- or when x is
jointly Gaussian, by Stein's lemma, Cov(g(X),X) = Var(X) E[g'(X)], whatever
the law of u.  The second is the one that matters in practice, and it fails
for every regressor a household survey carries.

WHAT THE TABLE SHOWS.  The law of the REGRESSOR decides, not the law of the
error: with a Gaussian x the projection is exact even when u is strongly
skewed, and with a skewed x it is far from exact even when u is Gaussian.

Everything below is quadrature, not simulation: there is no Monte Carlo noise
and no density is estimated anywhere.  The convergence table printed at the
end is the evidence that the numbers have settled.

    python p10_rif_gap.py

Needs numpy and scipy.  Takes about a minute.
"""
import csv
import os
import numpy as np
from scipy import stats, optimize

TAUS = (0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)
NODES = 400_000          # midpoint rule on the probability scale
RATIO = 0.5              # b * sd(x) / sd(u), held fixed across the table


# --------------------------------------------------------------------------
# laws, each standardised to mean 0 and variance 1
# --------------------------------------------------------------------------
def law(dist, n=NODES):
    """Quadrature nodes and weights for a standardised version of dist."""
    p = (np.arange(n) + 0.5) / n
    x = dist.ppf(p)
    w = np.full(n, 1.0 / n)
    m = w @ x
    s = np.sqrt(w @ (x - m) ** 2)
    return (x - m) / s, w


def lognormal_error(s):
    """cdf and pdf of a lognormal(0, s) standardised to variance one."""
    mean = np.exp(s ** 2 / 2)
    sd = np.sqrt((np.exp(s ** 2) - 1) * np.exp(s ** 2))
    skew = (np.exp(s ** 2) + 2) * np.sqrt(np.exp(s ** 2) - 1)
    F = lambda t: stats.lognorm.cdf(t * sd + mean, s)
    f = lambda t: stats.lognorm.pdf(t * sd + mean, s) * sd
    return F, f, skew


ERRORS = (("normal", stats.norm.cdf, stats.norm.pdf, 0.0),
          ("lognormal", *lognormal_error(0.9)))

X_LAWS = (("normal", stats.norm),
          ("uniform", stats.uniform),
          ("Student t(5)", stats.t(5)),
          ("chi-square(3)", stats.chi2(3)),
          ("exponential", stats.expon),
          ("lognormal", stats.lognorm(1.0)))


# --------------------------------------------------------------------------
def limits(x, w, b, F_u, f_u, tau):
    """Population limits of the two estimators at tau, by the same quadrature.

    Returns (RIF-regression slope, two-step limit).  The two-step limit is
    computed, not asserted: it is the percentile-weighted mean of the unit
    effects, sum_i omega_tau(x_i) b / sum_i omega_tau(x_i), with the weights
    omega_tau(x) = f(q_tau|x)/f(q_tau) that Proposition 2 shows the percentile
    weights converge to.  In this design b does not vary, so the weighted mean
    returns b whatever the weights are, and the figure the programme prints is
    a check on the quadrature as much as on the estimator.
    """
    mx = w @ x
    vx = w @ (x - mx) ** 2
    cdf_y = lambda q: w @ F_u(q - b * x)
    q = optimize.brentq(lambda q: cdf_y(q) - tau, -200, 200,
                        xtol=1e-14, rtol=8.9e-16)
    f_y = w @ f_u(q - b * x)

    Fc = F_u(q - b * x)
    cov = w @ ((Fc - w @ Fc) * (x - mx))
    gamma_rif = -cov / (f_y * vx)

    omega = f_u(q - b * x) / f_y                 # f(q_tau|x)/f(q_tau)
    theta_pwr = (w @ (omega * b)) / (w @ omega)

    return gamma_rif, theta_pwr


def main():
    os.makedirs("results", exist_ok=True)
    rows = []
    print("departure of each estimator's population limit from the true\n"
          f"effect, in per cent;  b sd(x)/sd(u) = {RATIO}\n")
    for u_name, F_u, f_u, skew in ERRORS:
        print(f"u {u_name} (skewness {skew:.2f})")
        print(f"{'law of x':<16}" + "".join(f"{t:>9.2f}" for t in TAUS))
        for x_name, dist in X_LAWS:
            x, w = law(dist)
            line, line2 = [], []
            for tau in TAUS:
                g, t2 = limits(x, w, RATIO, F_u, f_u, tau)
                gap = 100.0 * (g / RATIO - 1.0)
                gap2 = 100.0 * (t2 / RATIO - 1.0)
                line.append(f"{gap:>9.2f}")
                line2.append(f"{gap2:>9.2f}")
                rows.append(dict(u_law=u_name, u_skewness=f"{skew:.4f}",
                                 x_law=x_name, tau=tau, b=RATIO,
                                 gamma_rif=f"{g:.10f}",
                                 theta_twostep=f"{t2:.10f}",
                                 truth=f"{RATIO:.10f}",
                                 rif_gap_pct=f"{gap:.4f}",
                                 twostep_gap_pct=f"{gap2:.4f}"))
            print(f"  RIF      {x_name:<15}" + "".join(line))
            print(f"  two-step {x_name:<15}" + "".join(line2))
        print()

    with open("results/p10_rif_gap.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)
    print("written: results/p10_rif_gap.csv")

    # ---- the two things the table is claimed to show, checked separately --
    F_ln, f_ln, _ = ERRORS[1][1], ERRORS[1][2], None
    print("\nconvergence: x normal, u lognormal, tau = 0.25 "
          "(the exactness claimed by Stein's lemma)")
    for n in (25_000, 100_000, 400_000, 1_600_000):
        x, w = law(stats.norm, n)
        g, _ = limits(x, w, RATIO, F_ln, f_ln, 0.25)
        print(f"   {n:>9,} nodes   gap {100*(g/RATIO-1):+10.5f} %")

    print("\nthe gap against the signal-to-noise ratio, x lognormal, "
          "u lognormal")
    x, w = law(stats.lognorm(1.0))
    print(f"{'b sd(x)/sd(u)':>14}" + "".join(f"{t:>9.2f}" for t in TAUS))
    for b in (0.02, 0.05, 0.1, 0.25, 0.5, 1.0):
        line = "".join(f"{100*(limits(x, w, b, F_ln, f_ln, t)[0]/b-1):>9.2f}"
                       for t in TAUS)
        print(f"{b:>14.2f}" + line)
    print("\nThe gap does shrink as the signal shrinks, since the projection "
          "and the\naverage derivative agree in the limit, but slowly when x "
          "is strongly\nskewed: at b sd(x)/sd(u) = 0.02, a ratio at which the "
          "regressor explains\nalmost nothing, it is still -36 per cent at "
          "tau = 0.05.  A weak model is\nno protection.")


if __name__ == "__main__":
    main()
