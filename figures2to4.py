"""
make_results_figures.py
Figures 2-4 for the JPKPD adaptation, built directly from the verified modules
(structural_identifiability.py, analysis.py, diagnostics_S4.py) -- no numbers here
are re-derived independently; this script only plots what those modules compute.
"""
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import chi2

from model import CNBIO, DrugParams
from analysis import generate_synthetic_data, fit_mechanistic_model, practical_identifiability
import structural_identifiability as si
import diagnostics_S4 as d4

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(OUT, exist_ok=True)

mpl.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 8, "axes.linewidth": 0.6,
    "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "legend.fontsize": 7, "pdf.fonttype": 42,
    "mathtext.fontset": "custom", "mathtext.rm": "Liberation Sans", "mathtext.it": "Liberation Sans:italic",
})
W = 174 / 25.4

# ============================================================== Fig 2: structural degeneracy
drug_true = DrugParams(CLint_uL_per_min=4.5, CLd_uL_per_min=8.0, Kp=5.0, kb_per_min=0.0008)
t_eval = np.linspace(1, 2000, 300)
drug_alt, Cm_true, Cm_alt = si.demonstrate_degeneracy(CNBIO, drug_true, kb_alt=0.003, t_eval=t_eval)
rel_diff = np.abs((Cm_true - Cm_alt) / Cm_true)

fig, ax = plt.subplots(1, 2, figsize=(W * 0.72, 2.6))
a = ax[0]
a.semilogy(t_eval / 60, Cm_true, color="black", lw=1.8,
          label=f"$CL_d$={drug_true.CLd_uL_per_min}, $K_p$={drug_true.Kp}, $CL_{{int}}$={drug_true.CLint_uL_per_min}, $k_b$={drug_true.kb_per_min}")
a.semilogy(t_eval / 60, Cm_alt, color="#D55E00", lw=1.0, ls=(0, (4, 2)),
          label=f"$CL_d$={drug_alt.CLd_uL_per_min:.2f}, $K_p$={drug_alt.Kp:.2f}, $CL_{{int}}$={drug_alt.CLint_uL_per_min:.2f}, $k_b$={drug_alt.kb_per_min}")
a.set_xlabel("Time (h)"); a.set_ylabel(r"$C_m(t)$ ($\mu$M)")
a.legend(frameon=False, loc="upper right", fontsize=6.2)
a.text(-0.18, 1.04, "a", transform=a.transAxes, fontsize=10, fontweight="bold")

b = ax[1]
b.semilogy(t_eval / 60, rel_diff, color="black", lw=1.2)
b.set_xlabel("Time (h)"); b.set_ylabel("Relative difference in " + r"$C_m(t)$")
b.set_ylim(1e-9, 1e-5)
b.text(-0.18, 1.04, "b", transform=b.transAxes, fontsize=10, fontweight="bold")
fig.tight_layout(w_pad=1.3)
fig.savefig(f"{OUT}/Fig2.pdf"); fig.savefig(f"{OUT}/Fig2.png", dpi=300)
plt.close(fig)
print(f"Fig2: max relative difference = {rel_diff.max():.2e}")

# ============================================================== Fig 3: FIM instability
scenarios = {
    "low-$CL_{int}$, dense (n=12)": (0.3, 20.6 * 60, 12),
    "high-$CL_{int}$, dense (n=12)": (15.0, 6.2 * 60, 12),
}
rows = []
for label, (CLint_total, duration_min, n) in scenarios.items():
    t_sample = d4.sampling_schedule(duration_min, n)
    drug = d4.make_drug(CLint_total)
    Cm_obs, _ = generate_synthetic_data(CNBIO, drug, d4.DOSE_UM, t_sample, d4.NOISE_CV, seed=d4.SEED)
    theta0 = [CLint_total * 1.3, d4.GROUND_TRUTH["CLd"] * 0.8, d4.GROUND_TRUTH["Kp"] * 1.2,
              d4.GROUND_TRUTH["kb"] * 1.5, CNBIO.e_uL_per_min * 1.2 + 1e-6]
    res = fit_mechanistic_model(t_sample, Cm_obs, d4.DOSE_UM, CNBIO, theta0)
    J = res.jac; JTJ = J.T @ J
    cond = np.linalg.cond(JTJ)
    n_obs = len(res.fun); dof = max(n_obs - len(res.x), 1)
    sigma2 = (res.fun @ res.fun) / dof
    se_plain = np.sqrt(np.abs(np.diag(sigma2 * np.linalg.inv(JTJ))))
    se_pinv = np.sqrt(np.abs(np.diag(sigma2 * np.linalg.pinv(JTJ, rcond=1e-10))))
    rows.append((label, cond, se_plain[0] / abs(res.x[0]), se_pinv[0] / abs(res.x[0])))
    print(f"Fig3 [{label}]: cond={cond:.3e}, relSE plain={se_plain[0]/abs(res.x[0]):.3e}, "
          f"relSE pinv={se_pinv[0]/abs(res.x[0]):.3e}")

fig, ax = plt.subplots(1, 2, figsize=(W * 0.72, 2.6))
labels = [r[0] for r in rows]
x = np.arange(len(labels))
ax[0].bar(x, [np.log10(r[1]) for r in rows], color="#7f7f7f", width=0.5)
ax[0].set_xticks(x); ax[0].set_xticklabels(labels, fontsize=6.5, rotation=15, ha="right")
ax[0].set_ylabel(r"log$_{10}$ condition number of $J^TJ$")
ax[0].text(-0.22, 1.04, "a", transform=ax[0].transAxes, fontsize=10, fontweight="bold")

w = 0.35
ax[1].bar(x - w / 2, [np.log10(r[2]) for r in rows], width=w, color="#D55E00", label="plain inverse")
ax[1].bar(x + w / 2, [np.log10(r[3]) for r in rows], width=w, color="#0072B2", label="pseudoinverse")
ax[1].set_xticks(x); ax[1].set_xticklabels(labels, fontsize=6.5, rotation=15, ha="right")
ax[1].set_ylabel(r"log$_{10}$ relative SE of $CL_{int}$")
ax[1].legend(frameon=False, loc="upper left", fontsize=6.5)
ax[1].text(-0.22, 1.04, "b", transform=ax[1].transAxes, fontsize=10, fontweight="bold")
fig.tight_layout(w_pad=1.6)
fig.savefig(f"{OUT}/Fig3.pdf"); fig.savefig(f"{OUT}/Fig3.png", dpi=300)
plt.close(fig)

# ============================================================== Fig 4: single-start vs multi-start profile likelihood
duration_min, n = 20.6 * 60, 5
t_sample = d4.sampling_schedule(duration_min, n)
drug_low = d4.make_drug(d4.CLINT_LOW)
Cm_obs, _ = generate_synthetic_data(CNBIO, drug_low, d4.DOSE_UM, t_sample, d4.NOISE_CV, seed=d4.SEED)
ratios_wide = np.array([0.2, 0.48, 0.96, 1.92, 4.8, 9.6, 19.2])
deltas_single = d4.profile_clint_single_start(t_sample, Cm_obs, drug_low, fixed={}, ratios=ratios_wide)

ratios_local = np.array([0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0])
ms_scenarios = [
    ("sparse (n=5), $k_b$, $e$ unknown", 20.6 * 60, 5, {}),
    ("sparse (n=5), $k_b$, $e$ known", 20.6 * 60, 5, {"kb": d4.GROUND_TRUTH["kb"], "e": CNBIO.e_uL_per_min}),
    ("dense (n=12), $k_b$, $e$ known", 20.6 * 60, 12, {"kb": d4.GROUND_TRUTH["kb"], "e": CNBIO.e_uL_per_min}),
]
ms_results = {}
for label, dur, nn, fixed_base in ms_scenarios:
    ts = d4.sampling_schedule(dur, nn)
    Cm_o, _ = generate_synthetic_data(CNBIO, drug_low, d4.DOSE_UM, ts, d4.NOISE_CV, seed=d4.SEED)
    deltas, sigma2 = d4.multistart_profile(ts, Cm_o, d4.CLINT_LOW, fixed_base, ratios_local)
    ms_results[label] = (deltas, sigma2)
    threshold = chi2.ppf(0.95, df=1) * sigma2
    n_inside = np.sum(deltas <= threshold)
    print(f"Fig4 [{label}]: {n_inside}/{len(ratios_local)} inside 95% CI")

fig, ax = plt.subplots(1, 2, figsize=(W, 2.7))
a = ax[0]
order = np.argsort(ratios_wide)
a.plot(ratios_wide[order], deltas_single[order], "o-", color="black", ms=4, lw=1.0)
a.set_xscale("log")
a.set_xlabel(r"$CL_{int}$ / fitted $CL_{int}$")
a.set_ylabel(r"$\Delta$SSE (profile $-$ optimum)")
a.set_title("Single-start (naive)", fontsize=8)
a.text(-0.18, 1.06, "a", transform=a.transAxes, fontsize=10, fontweight="bold")

b = ax[1]
colors = ["#D55E00", "#0072B2", "#009E73"]
for (label, (deltas, sigma2)), col in zip(ms_results.items(), colors):
    threshold = chi2.ppf(0.95, df=1) * sigma2
    b.plot(ratios_local, deltas, "o-", color=col, ms=4, lw=1.0, label=label)
    b.axhline(threshold, color=col, ls=":", lw=0.7)
b.set_xscale("log")
b.set_xlabel(r"$CL_{int}$ / true $CL_{int}$")
b.set_ylabel(r"$\Delta$SSE (profile $-$ optimum)")
b.set_title("Multi-start (6 restarts)", fontsize=8)
b.legend(frameon=False, loc="upper center", fontsize=6.0)
b.text(-0.18, 1.06, "b", transform=b.transAxes, fontsize=10, fontweight="bold")
fig.tight_layout(w_pad=1.6)
fig.savefig(f"{OUT}/Fig4.pdf"); fig.savefig(f"{OUT}/Fig4.png", dpi=300)
plt.close(fig)
print("saved Fig2, Fig3, Fig4")
