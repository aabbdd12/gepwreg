"""The exact value of theta_j(tau) = E[b_ij | y = q_tau] in the Monte Carlo
DGPs of the paper, by numerical integration.

WHY THIS FILE EXISTS.  The Monte Carlo programs p4, p5 and p6 used, as their
"truth", the mean effect of the units whose rank in y falls within +-0.05 of
tau.  That reference is not the estimand: it is the estimand smoothed over a
window of half-width 0.05.  It was harmless while the bandwidth was fixed at
the Silverman rule and only biases of order 0.1 were at stake; it stopped
being harmless the moment the bandwidth became the object of study.  The
kernel of gepwe has standard deviation h*sqrt(2) on the rank scale, which at
the Silverman rule with n = 4000 to 20000 is 0.050 to 0.059 -- the half-width
of the window itself.  A reference smoothed over the same width as the
estimator's own kernel rewards that kernel: it flatters the Silverman rule by
construction, and any bandwidth comparison run against it is rigged.

A Monte Carlo reference on a very large sample is no substitute.  On 400,000
draws with a +-0.005 window the error still reaches 0.007, against an RMSE of
0.008 in the experiment it is meant to arbitrate.

Every DGP of the paper has a latent variable u that carries the effect, and a
conditional law of y given u that is known in closed form.  Then

    theta(tau) = E_u[ b(u) f(q_tau | u) ] / E_u[ f(q_tau | u) ],
    q_tau solving E_u[ F(q | u) ] = tau,

a one-dimensional integral, computed here by the midpoint rule on 2,000,000
points (error O(N^-2), below 1e-12 for these integrands) and a bisection on
q.  U1 also has a closed form, used as a check.

  D3  y = 5 z + (1+z) x,          z ~ U(0,1), x ~ N(0,1)
      y | z ~ N(5z, (1+z)^2),                             b = 1 + z
  H1  y = 5 + (1+z) x + z + e,    z ~ U(0,1), x, e ~ N(0,1)
      y | z ~ N(5+z, (1+z)^2 + 1),                        b = 1 + z
  D1  y = 100 z + (1+z) x
      y | z ~ N(100z, (1+z)^2),                           b = 1 + z
  U1  y = 5 + (1+U) x + 3U,       U ~ U(0,1), x ~ U(0,2)
      y | U ~ U(5+3U, 7+5U),                              b = 1 + U
  P6  the H1 design under a survey plan: 10 strata x 40 PSUs x 20 households,
      y = 5 + s_h + (1+z) x + z + a_c + e,  x = c_c + N(0,1),
      a_c, c_c ~ N(0, 0.6^2), s_h = 0.1 (h - 5.5), weights w_h = 1 + 0.1 h.
      a and c integrate out in closed form: given z and the stratum,
          y | z, h ~ N(5 + s_h + z, 1.36 [(1+z)^2 + 1]),
      because Var[(1+z) c + a + e | z] = 0.36 (1+z)^2 + 0.36 + 1.  The
      population is the WEIGHTED one, so stratum h enters with weight w_h.

Output: the values, a Stata block ready to paste, and a check against a
large simulation with two window widths extrapolated to zero.

Run from python/:  python pwr_exact_truth.py
"""
import numpy as np
from scipy.stats import norm

N = 2_000_000
TAUS4 = (0.1, 0.25, 0.5, 0.75, 0.9)
TAUS5 = (0.1, 0.5, 0.9)
TAUS6 = (0.25, 0.75)


def midpoints(n=N):
    return (np.arange(n) + 0.5) / n


# --- the normal-mixture DGPs ------------------------------------------------
NORMAL = {                       # name -> (mu(u), sd(u)), b(u) = 1 + u
    "D3": (lambda u: 5.0 * u, lambda u: 1.0 + u),
    "H1": (lambda u: 5.0 + u, lambda u: np.sqrt((1.0 + u) ** 2 + 1.0)),
    "D1": (lambda u: 100.0 * u, lambda u: 1.0 + u),
}
BRACKET = {"D3": (-10.0, 20.0), "H1": (-10.0, 25.0), "D1": (-15.0, 120.0)}


def normal_ref(name, tau, u=None):
    mu_f, sd_f = NORMAL[name]
    if u is None:
        u = midpoints()
    mu, sd = mu_f(u), sd_f(u)
    lo, hi = BRACKET[name]
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if norm.cdf((mid - mu) / sd).mean() < tau:
            lo = mid
        else:
            hi = mid
    q = 0.5 * (lo + hi)
    d = norm.pdf((q - mu) / sd) / sd
    return q, float(((1.0 + u) * d).mean() / d.mean())


# --- U1: uniform conditional law, closed form and quadrature ----------------
def u1_ref(tau):
    """y | U ~ U(5+3U, 7+5U).  The set of U compatible with q is
    [max(0,(q-7)/5), min(1,(q-5)/3)], the conditional density is 1/(2+2u) and
    b(u)/(2+2u) = 1/2 exactly, so
        theta = |A| / ln((1+u_hi)/(1+u_lo)).
    """
    def cdf(q):
        u = midpoints()
        return np.clip((q - 5.0 - 3.0 * u) / (2.0 + 2.0 * u), 0.0, 1.0).mean()

    lo, hi = 5.0, 12.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if cdf(mid) < tau:
            lo = mid
        else:
            hi = mid
    q = 0.5 * (lo + hi)
    ulo, uhi = max(0.0, (q - 7.0) / 5.0), min(1.0, (q - 5.0) / 3.0)
    closed = (uhi - ulo) / np.log((1.0 + uhi) / (1.0 + ulo))
    u = midpoints()
    d = np.where((u >= ulo) & (u <= uhi), 1.0 / (2.0 + 2.0 * u), 0.0)
    quad = float(((1.0 + u) * d).mean() / d.mean())
    return q, closed, quad


# --- P6: the survey design --------------------------------------------------
def p6_ref(tau, n=200_000):
    h = np.arange(1, 11)
    w = 1.0 + 0.1 * h
    w = w / w.sum()
    s = 0.1 * (h - 5.5)
    z = midpoints(n)
    sd = np.sqrt(1.36 * ((1.0 + z) ** 2 + 1.0))
    lo, hi = -10.0, 25.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        F = sum(wk * norm.cdf((mid - 5.0 - sk - z) / sd).mean()
                for wk, sk in zip(w, s))
        if F < tau:
            lo = mid
        else:
            hi = mid
    q = 0.5 * (lo + hi)
    num = den = 0.0
    for wk, sk in zip(w, s):
        d = norm.pdf((q - 5.0 - sk - z) / sd) / sd
        num += wk * ((1.0 + z) * d).mean()
        den += wk * d.mean()
    return q, float(num / den)


# --- the window reference, for the record -----------------------------------
def window_check(name, tau, n=4_000_000, seed=11):
    """E[b | |rank(y) - tau| <= w] on a large sample, for two widths, plus the
    linear extrapolation to w -> 0.  Shows how far the +-0.05 window sits from
    the estimand."""
    rng = np.random.default_rng(seed)
    u = rng.random(n)
    if name == "U1":
        y = 5.0 + (1.0 + u) * (2.0 * rng.random(n)) + 3.0 * u
    else:
        mu_f, sd_f = NORMAL[name]
        y = mu_f(u) + sd_f(u) * rng.standard_normal(n)
    b = 1.0 + u
    o = np.argsort(y, kind="stable")
    pc = np.empty(n)
    pc[o] = (np.arange(n) + 1.0) / n
    out = []
    for wd in (0.05, 0.02, 0.01):
        m = np.abs(pc - tau) <= wd
        out.append(b[m].mean())
    # E[b|window] is even in w to first order, so extrapolate in w^2
    a = (out[2] * 0.02 ** 2 - out[1] * 0.01 ** 2) / (0.02 ** 2 - 0.01 ** 2)
    return out, a


def main():
    print(__doc__.split("Run from")[0].strip()[:0] or "", end="")
    print("EXACT theta(tau) = E[b | y = q_tau], by quadrature\n")
    print(f"{'dgp':>5}{'tau':>7}{'q_tau':>12}{'theta':>12}")
    stata = {}
    for name in ("D3", "H1", "D1"):
        u = midpoints()
        for tau in TAUS4:
            q, th = normal_ref(name, tau, u)
            print(f"{name:>5}{tau:7.2f}{q:12.5f}{th:12.6f}")
            stata[(name, tau)] = th
    for tau in TAUS4:
        q, closed, quad = u1_ref(tau)
        print(f"{'U1':>5}{tau:7.2f}{q:12.5f}{closed:12.6f}"
              f"   (quadrature {quad:.6f}, gap {abs(closed-quad):.2e})")
        stata[("U1", tau)] = closed

    print("\nP6, the survey design (weighted population)")
    print(f"{'tau':>7}{'q_tau':>12}{'theta':>12}")
    p6 = {}
    for tau in TAUS6:
        q, th = p6_ref(tau)
        print(f"{tau:7.2f}{q:12.5f}{th:12.6f}")
        p6[tau] = th

    print("\nHow far the +-0.05 window sits from the estimand"
          " (n = 4,000,000)")
    print(f"{'dgp':>5}{'tau':>7}{'w=.05':>11}{'w=.02':>11}{'w=.01':>11}"
          f"{'w->0':>11}{'exact':>11}{'bias of .05':>13}")
    for name in ("D3", "H1", "D1", "U1"):
        for tau in TAUS4:
            o, a = window_check(name, tau)
            ex = stata[(name, tau)]
            print(f"{name:>5}{tau:7.2f}{o[0]:11.5f}{o[1]:11.5f}{o[2]:11.5f}"
                  f"{a:11.5f}{ex:11.5f}{o[0]-ex:13.5f}")

    print("\n--- Stata block for p4 -------------------------------------")
    for name in ("D3", "H1", "D1", "U1"):
        line = f'    if "`dgp\'" == "{name}" {{'
        print(line)
        for tau in TAUS4:
            k = str(tau).replace(".", "p")
            print(f"        local TRUE_{k} = {stata[(name, tau)]:.6f}")
        print("    }")
    print("\n--- Stata block for p5 (H1) --------------------------------")
    for tau in TAUS5:
        k = str(tau).replace(".", "p")
        print(f"local true_z_{k} = {stata[('H1', tau)]:.6f}")
    print("\n--- Stata block for p6 -------------------------------------")
    for tau in TAUS6:
        k = str(tau).replace(".", "p")
        print(f"local true_{k} = {p6[tau]:.6f}")


if __name__ == "__main__":
    main()
