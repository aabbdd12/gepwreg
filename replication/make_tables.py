"""make_tables.py -- builds the LaTeX table bodies of the paper (version 3)
from the CSV files the do-files of this folder write into results/.

Run it from the folder that contains it:  python make_tables.py
It writes tables/t_*.tex beside itself, each one a tabular environment that
the paper includes with \\input.  Nothing in the paper's numbers is typed by
hand.

If a result the paper needs is absent, the script stops and names the
do-file that writes it, rather than producing a half-built set of tables.
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
OUT = os.path.join(HERE, "tables")
os.makedirs(OUT, exist_ok=True)

# which file writes which result, so that a missing one names its own cause
SOURCE = {
    "p1_onestep.csv": "p1_onestep_attenuation.do",
    "p1_sweep.csv": "p1_onestep_attenuation.do",
    "p1_goldberger.csv": "p1_onestep_attenuation.do",
    "p2_closed_doors.csv": "p2_closed_doors.do",
    "p2_band.csv": "p2_closed_doors.do",
    "p3_note2016.csv": "p3_replications.do",
    "p3_mc2023.csv": "p3_replications.do",
    "p4_twostep.csv": "p4_twostep_vs_rif.do",
    "p5_coverage.csv": "p5_coverage_mc.do",
    "p6_svy.csv": "p6_svy_mc.do",
    "p7_tableA.csv": "p7_burkina.do",
    "p7_tableB.csv": "p7_burkina.do",
    "p7_tableC.csv": "p7_burkina.do",
    "p7_tableD.csv": "p7_burkina.do or p7d_tableD.do",
}
_missing = [c for c in SOURCE if not os.path.exists(os.path.join(RES, c))]
if _missing:
    print("missing in results/, run the file that writes it:")
    for c in sorted(_missing):
        print(f"  {c:22s} <- {SOURCE[c]}")
    raise SystemExit(1)


def read(name):
    with open(os.path.join(RES, name), newline="") as f:
        return list(csv.DictReader(f))


def fl(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def f3(v):
    x = fl(v)
    return "." if x != x else f"{x:.3f}"


def f4(v):
    x = fl(v)
    return "." if x != x else f"{x:.4f}"


def f2(v):
    x = fl(v)
    return "." if x != x else f"{x:.2f}"


def f1(v):
    x = fl(v)
    return "." if x != x else f"{x:.1f}"


def fi(v):
    x = fl(v)
    return "." if x != x else f"{x:,.0f}"


def write(name, body):
    with open(os.path.join(OUT, name), "w", newline="\n") as f:
        f.write(body)
    print("wrote", name)


# ---------------------------------------------------------------------------
# Table: the one-step against the truth (p1_onestep.csv)
# ---------------------------------------------------------------------------
rows = read("p1_onestep.csv")
dgp_label = {"A": "A: $y=5+0.5x+e$",
             "B": "B: $y=5+0.5x+e^{0.5x}e$",
             "C": "C: $y=5+0.5x+e^{-0.5x}e$",
             "D": "D: $y=5+(1+U)x+U+e$"}
lines = [r"\begin{tabular}{lcrrrrrr}", r"\toprule",
         r"DGP & $\tau$ & Truth & One-step & s.e.\ (IF) & $h$ & $n_{\mathrm{eff}}$ & OLS \\",
         r"\midrule"]
for dgp in ["A", "B", "C", "D"]:
    first = True
    for r in [x for x in rows if x["dgp"] == dgp]:
        lab = dgp_label[dgp] if first else ""
        first = False
        lines.append(f"{lab} & {r['tau']} & {f3(r['true'])} & {f4(r['onestep'])} & {f4(r['se_IF'])} & "
                     f"{f3(r['h'])} & {fi(r['N_eff'])} & {f3(r['ols'])} \\\\")
    lines.append(r"\addlinespace")
lines += [r"\bottomrule", r"\end{tabular}"]
write("t_p1_onestep.tex", "\n".join(lines) + "\n")

# ---------------------------------------------------------------------------
# Table: the bandwidth sweep (p1_sweep.csv), A and D side by side
# ---------------------------------------------------------------------------
rows = read("p1_sweep.csv")
hs = []
for r in rows:
    if r["h"] not in hs:
        hs.append(r["h"])
A = {r["h"]: r for r in rows if r["dgp"] == "A"}
D = {r["h"]: r for r in rows if r["dgp"] == "D"}
lines = [r"\begin{tabular}{rrrrrr}", r"\toprule",
         r"& & \multicolumn{2}{c}{DGP A (truth 0.500)} & \multicolumn{2}{c}{DGP D (truth "
         + f3(D[hs[0]]["true"]) + r")} \\",
         r"\cmidrule(lr){3-4}\cmidrule(lr){5-6}",
         r"$h$ & $n_{\mathrm{eff}}$ & One-step & OLS & One-step & OLS \\", r"\midrule"]
for h in hs:
    lines.append(f"{h} & {fi(A[h]['N_eff'])} & {f4(A[h]['onestep'])} & {f3(A[h]['ols'])} & "
                 f"{f4(D[h]['onestep'])} & {f3(D[h]['ols'])} \\\\")
lines += [r"\bottomrule", r"\end{tabular}"]
write("t_p1_sweep.tex", "\n".join(lines) + "\n")

# ---------------------------------------------------------------------------
# Table: the Goldberger identity (p1_goldberger.csv)
# ---------------------------------------------------------------------------
rows = read("p1_goldberger.csv")
lines = [r"\begin{tabular}{rrrrrr}", r"\toprule",
         r"$h$ & $n_{\mathrm{eff}}$ & $\Var_w(y)/\Var(y)$ & $\Var_w(x)/\Var(x)$ & Ratio & $b_w/\beta$ \\",
         r"\midrule"]
for r in rows:
    lines.append(f"{r['h']} & {fi(r['N_eff'])} & {f4(r['Vw_y_over_V_y'])} & {f4(r['Vw_x_over_V_x'])} & "
                 f"{f4(r['ratio'])} & {f4(r['b_over_beta'])} \\\\")
lines += [r"\bottomrule", r"\end{tabular}"]
write("t_p1_goldberger.tex", "\n".join(lines) + "\n")

# ---------------------------------------------------------------------------
# Table: closed doors (p2_closed_doors.csv)
# ---------------------------------------------------------------------------
rows = read("p2_closed_doors.csv")
lines = [r"\begin{tabular}{lcrrrrr}", r"\toprule",
         r"DGP & $\tau$ & Truth & One-step & Retained-variance ratio & Variance-corrected & Double-weighted \\",
         r"\midrule"]
for dgp in ["A", "B", "C", "D"]:
    first = True
    for r in [x for x in rows if x["dgp"] == dgp]:
        lab = dgp if first else ""
        first = False
        lines.append(f"{lab} & {r['tau']} & {f3(r['true'])} & {f4(r['onestep'])} & {f4(r['ratio'])} & "
                     f"{f3(r['corrected'])} & {f4(r['double_weighted'])} \\\\")
    lines.append(r"\addlinespace")
lines += [r"\bottomrule", r"\end{tabular}"]
write("t_p2_closed.tex", "\n".join(lines) + "\n")

# ---------------------------------------------------------------------------
# Table: the hard band (p2_band.csv)
# ---------------------------------------------------------------------------
rows = read("p2_band.csv")
lab = {"OLS full sample": "OLS, full sample",
       "OLS in the band of y": "OLS on the band $0.35 \\le F(y) \\le 0.45$",
       "truncreg on the band of y": "Truncated-normal MLE on the same band of $y$ (\\texttt{truncreg})",
       "OLS in the band 0.35-0.45 of x": "OLS on the band $0.35 \\le F(x) \\le 0.45$",
       "OLS in the band 0.25-0.75 of x": "OLS on the band $0.25 \\le F(x) \\le 0.75$"}
lines = [r"\begin{tabular}{lrr}", r"\toprule", r"Estimator & Slope & s.e. \\", r"\midrule"]
for r in rows:
    lines.append(f"{lab[r['estimator']]} & {f4(r['slope'])} & {f4(r['se'])} \\\\")
lines += [r"\bottomrule", r"\end{tabular}"]
write("t_p2_band.tex", "\n".join(lines) + "\n")

# ---------------------------------------------------------------------------
# Table: the 2016 note (p3_note2016.csv)
# ---------------------------------------------------------------------------
rows = read("p3_note2016.csv")
A1 = [r for r in rows if r["example"] == "A1"]
A3 = [r for r in rows if r["example"] == "A3"]
lines = [r"\begin{tabular}{crrr@{\hspace{2em}}crrr}", r"\toprule",
         r"\multicolumn{4}{c}{A1: $y=1000+p\,x_1$, $x_1$ the rank} & \multicolumn{4}{c}{A3: $y=2000+p\,x_1+60u$, $x_1\perp p$} \\",
         r"\cmidrule(lr){1-4}\cmidrule(lr){5-8}",
         r"$\tau$ & Total derivative $2\tau$ & $h_{\mathrm{Sil}}/3$ & $h_{\mathrm{Sil}}$ & $\tau$ & Coefficient $\tau$ & $h_{\mathrm{Sil}}/3$ & $h_{\mathrm{Sil}}$ \\",
         r"\midrule"]
n = max(len(A1), len(A3))
for i in range(n):
    a = A1[i] if i < len(A1) else None
    b = A3[i] if i < len(A3) else None
    la = f"{a['tau']} & {f2(a['true'])} & {f3(a['onestep_h3'])} & {f3(a['onestep_hSil'])}" if a else " & & & "
    lb = f"{b['tau']} & {f2(b['true'])} & {f3(b['onestep_h3'])} & {f3(b['onestep_hSil'])}" if b else " & & & "
    lines.append(f"{la} & {lb} \\\\")
lines += [r"\bottomrule", r"\end{tabular}"]
write("t_p3_2016.tex", "\n".join(lines) + "\n")

# ---------------------------------------------------------------------------
# Table: the 2023 Monte Carlo (p3_mc2023.csv)
# ---------------------------------------------------------------------------
rows = read("p3_mc2023.csv")
lines = [r"\begin{tabular}{lcrrrrrr}", r"\toprule",
         r"$\sigma_\varepsilon$ ($R^2$; $\mathrm{sd}(x_1|y)/\mathrm{sd}(x_1)$) & $\tau$ & Truth & One-step on $y$ & Ranked on $x_1$ & RIF-OLS & OLS \\",
         r"\midrule"]
for sig in ["1", "10", "25"]:
    sub = [r for r in rows if r["sd_eps"] == sig]
    first = True
    for r in sub:
        lab = f"{sig} ({f2(r['R2'])}; {f2(r['sdx1_given_y'])})" if first else ""
        first = False
        lines.append(f"{lab} & {r['tau']} & {f1(r['true'])} & {f1(r['onestep_mean'])} [{f2(r['onestep_sd'])}] & "
                     f"{f1(r['rankvar_mean'])} [{f2(r['rankvar_sd'])}] & {f1(r['rif_mean'])} [{f2(r['rif_sd'])}] & "
                     f"{f1(r['ols_mean'])} \\\\")
    lines.append(r"\addlinespace")
lines += [r"\bottomrule", r"\end{tabular}"]
write("t_p3_2023.tex", "\n".join(lines) + "\n")

# ---------------------------------------------------------------------------
# Table: two-step against RIF and the truth (p4_twostep.csv)
# ---------------------------------------------------------------------------
rows = read("p4_twostep.csv")
dgp4 = {"D3": "D3: $y=5r+(1+r)x$",
        "H1": "H1: $y=5+(1+z)x+z+e$",
        "D1": "D1: $y=100r+(1+r)x$",
        "U1": "U1: $y=5+(1+U)x+3U$, $z$ irrelevant"}
lines = [r"\begin{tabular}{lcrrrrrr}", r"\toprule",
         r"DGP & $\tau$ & Truth & PWR-1 & PWR-2(z) & s.e.\ IF & PWR-2(qr) & RIF-OLS \\",
         r"\midrule"]
for dgp in ["D3", "H1", "D1", "U1"]:
    first = True
    for r in [x for x in rows if x["dgp"] == dgp]:
        lab = dgp4[dgp] if first else ""
        first = False
        lines.append(f"{lab} & {r['tau']} & {f3(r['true'])} & {f3(r['onestep_m'])} [{f3(r['onestep_sd'])}] & "
                     f"{f3(r['hetz_m'])} [{f3(r['hetz_sd'])}] & {f4(r['hetz_seIF'])} & "
                     f"{f3(r['hetqr_m'])} [{f3(r['hetqr_sd'])}] & {f3(r['rif_m'])} [{f3(r['rif_sd'])}] \\\\")
    lines.append(r"\addlinespace")
lines += [r"\bottomrule", r"\end{tabular}"]
write("t_p4_twostep.tex", "\n".join(lines) + "\n")

# ---------------------------------------------------------------------------
# Table: the MSE-optimal bandwidth against the Silverman rule, same DGPs and
# the same samples (p4_twostep.csv).  The truth being exact, the root mean
# squared error is sqrt(bias^2 + variance) and can be read off the table.
# ---------------------------------------------------------------------------


def _rmse(mean, sd, true):
    return ((float(mean) - float(true)) ** 2 + float(sd) ** 2) ** 0.5


lines = [r"\begin{tabular}{lcrrrrrr}", r"\toprule",
         r"& & & \multicolumn{2}{c}{PWR-2(z), $h^\ast$} "
         r"& \multicolumn{2}{c}{PWR-2(z), Silverman} & \\",
         r"\cmidrule(lr){4-5}\cmidrule(lr){6-7}",
         r"DGP & $\tau$ & $h^\ast/h_{\mathrm{Sil}}$ & Mean & RMSE "
         r"& Mean & RMSE & Gain \\",
         r"\midrule"]
for dgp in ["D3", "H1", "D1", "U1"]:
    first = True
    for r in [x for x in rows if x["dgp"] == dgp]:
        lab = dgp4[dgp] if first else ""
        first = False
        ro = _rmse(r["hetz_m"], r["hetz_sd"], r["true"])
        rs = _rmse(r["hetzsil_m"], r["hetzsil_sd"], r["true"])
        gain = 100.0 * (rs - ro) / rs
        lines.append(f"{lab} & {r['tau']} & {f2(r['h_ratio'])} & "
                     f"{f3(r['hetz_m'])} & {ro:.4f} & "
                     f"{f3(r['hetzsil_m'])} & {rs:.4f} & {gain:+.1f}\\% \\\\")
    lines.append(r"\addlinespace")
lines += [r"\bottomrule", r"\end{tabular}"]
write("t_p4_bandwidth.tex", "\n".join(lines) + "\n")

# ---------------------------------------------------------------------------
# Table: IF variance and coverage (p5_coverage.csv)
# ---------------------------------------------------------------------------
rows = read("p5_coverage.csv")
lines = [r"\begin{tabular}{lcrrrrrrr}", r"\toprule",
         r"Estimator & $\tau$ & Object & Mean & MC s.d. & Mean s.e.\ IF & Ratio & Cov.\ (MC mean) & Cov.\ (object) \\",
         r"\midrule"]
for est, pwr in [("het(z)", "PWR-2(z)"), ("rankvar(z)", "PWR-cov(z)")]:
    first = True
    for r in [x for x in rows if x["estimator"] == est]:
        lab = pwr if first else ""
        first = False
        lines.append(f"{lab} & {r['tau']} & {f4(r['true_object'])} & {f4(r['mean_est'])} & {f4(r['mc_sd'])} & "
                     f"{f4(r['mean_seIF'])} & {f3(r['ratio_seIF_mcsd'])} & {f3(r['coverage95_mcmean'])} & "
                     f"{f3(r['coverage95_true'])} \\\\")
    lines.append(r"\addlinespace")
lines += [r"\bottomrule", r"\end{tabular}"]
write("t_p5_coverage.tex", "\n".join(lines) + "\n")

# ---------------------------------------------------------------------------
# Table: survey design (p6_svy.csv)
# ---------------------------------------------------------------------------
rows = read("p6_svy.csv")
lines = [r"\begin{tabular}{crrrrrrrrrr}", r"\toprule",
         r"& & & & \multicolumn{3}{c}{Taylor (PSU, strata)} & \multicolumn{3}{c}{IF ignoring the design} & Cluster boot. \\",
         r"\cmidrule(lr){5-7}\cmidrule(lr){8-10}",
         r"$\tau$ & Truth & Mean & MC s.d. & Mean s.e. & Ratio & Cov. & Mean s.e. & Ratio & Cov. & one sample \\",
         r"\midrule"]
for r in rows:
    lines.append(f"{r['tau']} & {f4(r['true_exact'])} & {f4(r['mean_est'])} & {f4(r['mc_sd'])} & "
                 f"{f4(r['mean_se_taylor'])} & {f3(r['ratio_taylor'])} & {f3(r['cov95_taylor'])} & "
                 f"{f4(r['mean_se_IF_nodesign'])} & {f3(r['ratio_IF'])} & {f3(r['cov95_IF'])} & "
                 f"{f4(r['boot_cluster_se_one_sample'])} \\\\")
lines += [r"\bottomrule", r"\end{tabular}"]
write("t_p6_svy.tex", "\n".join(lines) + "\n")

# ---------------------------------------------------------------------------
# Tables: Burkina Faso, Table A (p7_tableA.csv)
# ---------------------------------------------------------------------------
rows = read("p7_tableA.csv")
est_lab = {"het(urban size) Taylor": r"PWR-2(z), \texttt{het(urban size)}, Taylor",
           "het(qr) boot50": r"PWR-2(qr), cluster bootstrap",
           "RIF-OLS": r"RIF-OLS (\texttt{rifhdreg})",
           "one-step rankdep (<=1.4)": r"PWR-1, rank of $y$ (\texttt{rankdep})"}
ests = list(est_lab)
taus = ["0.1", "0.25", "0.5", "0.75", "0.9"]


def cell(t, e, v, what):
    r = [x for x in rows if x["tau"] == t and x["estimator"] == e and x["variable"] == v]
    if not r:
        return "."
    r = r[0]
    return f4(r["coef"]) if what == "coef" else f"({f4(r['se'])})"


# household size and male head, one table with a panel each
lines = [r"\begin{tabular}{lrrrrr}", r"\toprule",
         r"Estimator & $\tau=0.1$ & $\tau=0.25$ & $\tau=0.5$ & $\tau=0.75$ & $\tau=0.9$ \\"]
for v, label in [("size", r"\textit{Household size}"),
                 ("male", r"\textit{Male head}")]:
    lines.append(r"\midrule")
    lines.append(r"\multicolumn{6}{l}{" + label + r"} \\[2pt]")
    for e in ests:
        lines.append(est_lab[e] + " & " + " & ".join(cell(t, e, v, "coef") for t in taus) + r" \\")
        lines.append(" & " + " & ".join(cell(t, e, v, "se") for t in taus) + r" \\[2pt]")
ne = [x for x in rows if x["estimator"] == ests[0] and x["variable"] == "size"]
lines.append(r"\midrule")
lines.append(r"$n_{\mathrm{eff}}$ (two-step) & " + " & ".join(fi(x["N_eff"]) for x in ne) + r" \\")
lines += [r"\bottomrule", r"\end{tabular}"]
write("t_p7_sizemale.tex", "\n".join(lines) + "\n")

# socio-economic groups, three estimators, tau = 0.1, 0.5, 0.9
gse_lab = {"2.gse": "Private wage-earner", "3.gse": "Artisan, trader", "4.gse": "Other earner",
           "5.gse": "Crop farmer", "6.gse": "Subsistence farmer", "7.gse": "Inactive"}
lines = [r"\begin{tabular}{lrrrrrrrrr}", r"\toprule",
         r"& \multicolumn{3}{c}{PWR-2(z)} & \multicolumn{3}{c}{PWR-2(qr)} & \multicolumn{3}{c}{RIF-OLS} \\",
         r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}\cmidrule(lr){8-10}",
         r"Group (ref.\ public wage-earner) & 0.1 & 0.5 & 0.9 & 0.1 & 0.5 & 0.9 & 0.1 & 0.5 & 0.9 \\",
         r"\midrule"]
for v in ["2.gse", "3.gse", "4.gse", "5.gse", "6.gse", "7.gse"]:
    cells = []
    for e in ests[:3]:
        for t in ["0.1", "0.5", "0.9"]:
            r = [x for x in rows if x["tau"] == t and x["estimator"] == e and x["variable"] == v][0]
            cells.append(f3(r["coef"]))
    lines.append(gse_lab[v] + " & " + " & ".join(cells) + r" \\")
lines += [r"\bottomrule", r"\end{tabular}"]
write("t_p7_gse.tex", "\n".join(lines) + "\n")

# ---------------------------------------------------------------------------
# Table B: composition (p7_tableB.csv)
# ---------------------------------------------------------------------------
rows = read("p7_tableB.csv")
lines = [r"\begin{tabular}{crrrrrr}", r"\toprule",
         r"& \multicolumn{2}{c}{Urban share} & \multicolumn{2}{c}{Household size} & \multicolumn{2}{c}{Contribution of size to} \\",
         r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7}",
         r"$\tau$ & at $\tau$ & population & at $\tau$ & population & effect of size & effect of male \\",
         r"\midrule"]
for t in taus:
    u = [x for x in rows if x["tau"] == t and x["zvar"] == "urban"][0]
    s = [x for x in rows if x["tau"] == t and x["zvar"] == "size"][0]
    lines.append(f"{t} & {f3(u['mean_at_tau'])} & {f3(u['mean_pop'])} & {f2(s['mean_at_tau'])} & {f2(s['mean_pop'])} & "
                 f"{f4(s['contrib_size'])} & {f4(s['contrib_male'])} \\\\")
lines += [r"\bottomrule", r"\end{tabular}"]
write("t_p7_comp.tex", "\n".join(lines) + "\n")

# ---------------------------------------------------------------------------
# Table C: profile along household size (p7_tableC.csv)
# ---------------------------------------------------------------------------
rows = read("p7_tableC.csv")
lines = [r"\begin{tabular}{lrrrrr}", r"\toprule",
         r"Variable & $\tau_{\mathrm{size}}=0.1$ & $0.25$ & $0.5$ & $0.75$ & $0.9$ \\", r"\midrule"]
for v, lab in [("male", "Male head"), ("urban", "Urban"), ("5.gse", "Crop farmer"),
               ("6.gse", "Subsistence farmer"), ("7.gse", "Inactive")]:
    cells = []
    for t in taus:
        r = [x for x in rows if x["tau_of_size"] == t and x["variable"] == v][0]
        cells.append(f"{f3(r['coef'])} ({f3(r['se'])})")
    lines.append(lab + " & " + " & ".join(cells) + r" \\")
ne = [[x for x in rows if x["tau_of_size"] == t][0]["N_eff"] for t in taus]
lines.append(r"\addlinespace")
lines.append(r"$n_{\mathrm{eff}}$ & " + " & ".join(fi(x) for x in ne) + r" \\")
lines += [r"\bottomrule", r"\end{tabular}"]
write("t_p7_profile.tex", "\n".join(lines) + "\n")

# ---------------------------------------------------------------------------
# Table D: the initial quantile (p7_tableD.csv) -- new layout (with het column)
# or old layout (without), whichever the file has
# ---------------------------------------------------------------------------
rows = read("p7_tableD.csv")
if rows and "het" in rows[0]:
    lines = [r"\begin{tabular}{lcrrrrr}", r"\toprule",
             r"Heterogeneity model & $\tau$ & At the initial quantile & s.e.\ (boot.) & At the quantile of $y$ & Size at $\tau$ of $y_0$ & Size at $\tau$ of $y$ \\",
             r"\midrule"]
    for het, pwr in [("het(size)", "PWR-2(z), \\texttt{het(size)}"), ("het(qr)", "PWR-2(qr)")]:
        first = True
        for r in [x for x in rows if x["het"] == het]:
            lab = pwr if first else ""
            first = False
            lines.append(f"{lab} & {r['tau']} & {f4(r['urban_at_initial_quantile'])} & {f4(r['se_boot'])} & "
                         f"{f4(r['urban_at_quantile_of_lexp'])} & {f2(r['size_at_initial_quantile'])} & "
                         f"{f2(r['size_at_quantile_of_lexp'])} \\\\")
        lines.append(r"\addlinespace")
else:
    lines = [r"\begin{tabular}{crrr}", r"\toprule",
             r"$\tau$ & At the initial quantile & s.e.\ (boot.) & At the quantile of $y$ \\", r"\midrule"]
    for r in rows:
        lines.append(f"{r['tau']} & {f4(r['urban_at_initial_quantile'])} & {f4(r['se_boot'])} & "
                     f"{f4(r['urban_at_quantile_of_lexp'])} \\\\")
lines += [r"\bottomrule", r"\end{tabular}"]
write("t_p7_initial.tex", "\n".join(lines) + "\n")

# ---------------------------------------------------------------------------
# Table: the dispersion of the group at tau, Burkina (p8_dispersion.csv)
# ---------------------------------------------------------------------------
# It belongs to a section that is not in version 3 of the paper, and the
# script that wrote the CSV is not part of this package; skip it if absent.
_p8 = os.path.join(RES, "p8_dispersion.csv")
if os.path.exists(_p8):
    rows = read("p8_dispersion.csv")
    lines = [r"\begin{tabular}{crrrrrr}", r"\toprule",
             r"& \multicolumn{3}{c}{Dispersion the group would have} & \multicolumn{2}{c}{Observed in the group} & \\",
             r"\cmidrule(lr){2-4}\cmidrule(lr){5-6}",
             r"$\tau$ & within types & between types & total & $\mathrm{sd}_w(y)$ & retained & $n_{\mathrm{eff}}$ \\",
             r"\midrule"]
    for r in rows:
        lines.append(f"{f2(r['tau'])} & {f4(r['sd_within'])} & {f4(r['sd_between'])} & {f4(r['sd_natural'])} & "
                     f"{f4(r['sd_observed'])} & {f4(r['retained'])} & {fi(r['n_eff'])} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    write("t_p8_dispersion.tex", "\n".join(lines) + "\n")
else:
    print("skipped t_p8_dispersion.tex: p8_dispersion.csv is not produced\n        by this package (the table is not in version 3 of the paper)")
