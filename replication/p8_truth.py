"""theta(0.75) for the two measurement-error designs of Section 5.

WHY THIS FILE EXISTS.  The design of p8_merr.do has 8,000 observations, and
the effect at a quantile cannot be read off it: a window of half-width 0.02
on the rank holds 320 units, so the Monte Carlo standard error of a window
mean is around 0.06 -- larger than every difference the section discusses.
Estimating the truth from the same small sample that is being corrected would
compare two noisy numbers and call the noise a bias.

The truth is a population quantity, theta(tau) = E[b | y = q_tau], and it
needs only a large enough draw.  The table below is its own evidence: the
estimate must stop moving both as the sample grows and as the window
narrows, and it does: 2.976 for the first design, 2.714 for the second.

    python p8_truth.py

Design A is that of p8_merr.do (Section 5.3): two heterogeneity variables,
both observed with error, the effect linear in them.  Design B is that of
p9_merr_curvature.do (Section 5.5): one heterogeneity variable observed
EXACTLY, with an effect quadratic in it.

Needs numpy.  Takes about two minutes.
"""
import numpy as np

TAU = 0.75


def theta(tau, n, halfwidth, seed):
    rng = np.random.default_rng(seed)
    g = rng.standard_normal(n)
    z1 = np.exp(0.45 * g)
    z2 = rng.random(n)
    x = 1 + 0.45 * g + np.sqrt(1 - 0.45 ** 2) * rng.standard_normal(n)
    b = 0.5 + 1.5 * z1 + 1.0 * z2
    y = 5 + b * x + 3 * z1 + 2 * z2 + rng.standard_normal(n)
    b_sorted = b[np.argsort(y)]
    p = np.arange(1, n + 1) / n
    keep = np.abs(p - tau) <= halfwidth
    return b_sorted[keep].mean(), int(keep.sum())


def theta_b(tau, n, halfwidth, seed):
    """Design B: z measured exactly, effect quadratic in z."""
    rng = np.random.default_rng(seed)
    g = rng.standard_normal(n)
    z = np.exp(0.45 * g)
    x = 1 + 0.45 * g + np.sqrt(1 - 0.45 ** 2) * rng.standard_normal(n)
    b = 0.5 + 1.5 * z + 0.8 * (z ** 2 - 1.5)
    y = 5 + b * x + 3 * z + rng.standard_normal(n)
    b_sorted = b[np.argsort(y)]
    p = np.arange(1, n + 1) / n
    keep = np.abs(p - tau) <= halfwidth
    return b_sorted[keep].mean(), int(keep.sum())


for name, fn, used in (("A  (p8_merr.do, Section 5.3)", theta, "2.976"),
                       ("B  (p9_merr_curvature.do, Section 5.5)", theta_b, "2.714")):
    print(f"Design {name}:  theta({TAU}) = E[b | y = q_tau]\n")
    print(f"{'n':>12} {'half-width':>11} {'units':>9} {'estimate':>10}")
    for n in (2_000_000, 10_000_000):
        for hw in (0.02, 0.01, 0.005, 0.002):
            v, k = fn(TAU, n, hw, seed=7 if fn is theta else 777)
            print(f"{n:>12,} {hw:>11} {k:>9,} {v:>10.4f}")
        print()
    print(f"    value used in the paper: {used}\n")

print("A window that is too narrow trades the smoothing bias for Monte Carlo")
print("noise and starts to wander, which is why several widths are reported")
print("rather than one: the number is established by the double stability, in")
print("the sample size and in the width, not by a single figure.")
