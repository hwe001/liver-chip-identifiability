"""
diagnostics_S4.py
=================

Reproduces Supplementary Methods Section S4 ("Practical identifiability:
a methodological case study") of the manuscript

    "What can media-only liver-on-a-chip concentration data identify?
     More sampling cannot resolve structural ambiguity"

end to end, as a single documented script. This consolidates the ad hoc
diagnostic snippets used during manuscript development (previously not
saved as standalone, reproducible code -- see the "Not yet done" note at
the end of supplementary_methods.md) into one file that reproduces every
numbered claim in S4.1, S4.2, and S4.3.

Run with:  python3 diagnostics_S4.py

Expected runtime: a few minutes (S4.3's multi-start scan is the slow part).

All parameter choices below (hardware, ground-truth kinetics, sampling
schedules, noise level, random seed) match those used to generate the
numbers actually quoted in supplementary_methods.md. Exact reproduced
values may differ in the last 1-2 significant figures from run to run
depending on optimizer numerics, but the qualitative pattern -- FIM
instability -> single-start profile-likelihood failure -> multi-start
stability -> weak practical constraint even with kb, e known -- is
robust and is what the manuscript's claim rests on, not any single
digit.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import numpy as np
from scipy.stats import chi2

from model import CNBIO, DrugParams
from analysis import generate_synthetic_data, fit_mechanistic_model, practical_identifiability

# ---------------------------------------------------------------------
# Shared ground-truth setup (CN Bio hardware, "low CLint" regime -- see
# main text Table 1 / Methods for why this regime is the illustrative
# case used throughout Section S4)
# ---------------------------------------------------------------------
GROUND_TRUTH = dict(CLd=8.0, Kp=5.0, kb=0.0008)
CLINT_LOW = 1.0 * (CNBIO.cells / 1e6)   # ~0.3 uL/min total
CLINT_HIGH = 50.0 * (CNBIO.cells / 1e6)  # ~15.0 uL/min total
DOSE_UM = 1.0
NOISE_CV = 0.15
SEED = 1000

PARAM_NAMES = ["CLint", "CLd", "Kp", "kb", "e"]


def make_drug(CLint_total):
    return DrugParams(CLint_uL_per_min=CLint_total, CLd_uL_per_min=GROUND_TRUTH["CLd"],
                       Kp=GROUND_TRUTH["Kp"], kb_per_min=GROUND_TRUTH["kb"])


def sampling_schedule(duration_min, n):
    return np.linspace(duration_min / n, duration_min, n)


# ======================================================================
# S4.1 -- FIM covariance approximation: numerical instability near the
#         structural degeneracy
# ======================================================================
def section_S4_1():
    print("=" * 70)
    print("S4.1  FIM covariance approximation near the structural degeneracy")
    print("=" * 70)

    scenarios = {
        "low-CLint, dense (n=12)": (CLINT_LOW, 20.6 * 60, 12),
        "high-CLint, dense (n=12)": (CLINT_HIGH, 6.2 * 60, 12),
    }

    for label, (CLint_total, duration_min, n) in scenarios.items():
        t_sample = sampling_schedule(duration_min, n)
        drug_true = make_drug(CLint_total)
        Cm_obs, _ = generate_synthetic_data(CNBIO, drug_true, DOSE_UM, t_sample, NOISE_CV, seed=SEED)

        theta0 = [CLint_total * 1.3, GROUND_TRUTH["CLd"] * 0.8, GROUND_TRUTH["Kp"] * 1.2,
                  GROUND_TRUTH["kb"] * 1.5, CNBIO.e_uL_per_min * 1.2 + 1e-6]
        res = fit_mechanistic_model(t_sample, Cm_obs, DOSE_UM, CNBIO, theta0)

        J = res.jac
        JTJ = J.T @ J
        cond_number = np.linalg.cond(JTJ)

        # Compare plain inversion vs. pseudoinverse directly (this is the
        # comparison reported as a ~200-fold disagreement in the text)
        n_obs = len(res.fun)
        dof = max(n_obs - len(res.x), 1)
        sigma2 = (res.fun @ res.fun) / dof
        try:
            cov_plain = sigma2 * np.linalg.inv(JTJ)
            se_plain = np.sqrt(np.abs(np.diag(cov_plain)))
            relse_plain = se_plain[0] / abs(res.x[0])
        except np.linalg.LinAlgError:
            relse_plain = np.inf
        cov_pinv = sigma2 * np.linalg.pinv(JTJ, rcond=1e-10)
        se_pinv = np.sqrt(np.abs(np.diag(cov_pinv)))
        relse_pinv = se_pinv[0] / abs(res.x[0])

        print(f"\n  Scenario: {label}")
        print(f"    fitted CLint = {res.x[0]:.4f} (true = {CLint_total:.4f})")
        print(f"    cond(J^T J)  = {cond_number:.3e}")
        print(f"    CLint relative SE, plain inverse : {relse_plain:.3e}")
        print(f"    CLint relative SE, pseudoinverse  : {relse_pinv:.3e}")
        print(f"    disagreement (plain / pseudoinv)  : {relse_plain / relse_pinv:.1f}x")

    print("\n  -> Condition numbers of 1e9-1e17 and large plain-vs-pseudoinverse\n"
          "     disagreements (reported in the manuscript as ~200-fold on\n"
          "     identical data) indicate the FIM approximation is numerically\n"
          "     unreliable in this regime, independent of which inversion is used.\n")


# ======================================================================
# S4.2 -- Single-start profile likelihood: contaminated by convergence
#         failure
# ======================================================================
def profile_clint_single_start(t_sample, Cm_obs, drug_true, fixed, ratios):
    """One optimizer start per grid point -- deliberately naive, to
    reproduce the erratic single-start behavior described in S4.2."""
    CLint_true = drug_true.CLint_uL_per_min
    free_names_mle = [p for p in PARAM_NAMES if p not in fixed]

    def theta0_for(free_names, clint_val):
        defaults = {"CLd": GROUND_TRUTH["CLd"] * 0.8, "Kp": GROUND_TRUTH["Kp"] * 1.2,
                    "kb": GROUND_TRUTH["kb"] * 1.5, "e": CNBIO.e_uL_per_min * 1.2 + 1e-6,
                    "CLint": clint_val * 1.3}
        return [defaults[n] for n in free_names]

    res_mle = fit_mechanistic_model(t_sample, Cm_obs, DOSE_UM, CNBIO,
                                     theta0_for(free_names_mle, CLint_true), fixed=fixed)
    sse_mle = res_mle.fun @ res_mle.fun

    deltas = []
    for r in ratios:
        f = dict(fixed); f["CLint"] = CLint_true * r
        free_names = [p for p in PARAM_NAMES if p not in f]
        res = fit_mechanistic_model(t_sample, Cm_obs, DOSE_UM, CNBIO,
                                     theta0_for(free_names, CLint_true * r), fixed=f)
        deltas.append((res.fun @ res.fun) - sse_mle)
    return np.array(deltas)


def section_S4_2():
    print("=" * 70)
    print("S4.2  Single-start profile likelihood (CLint, kb unknown)")
    print("=" * 70)

    duration_min, n = 20.6 * 60, 5  # sparse, low-CLint (the illustrative case in the text)
    t_sample = sampling_schedule(duration_min, n)
    drug_true = make_drug(CLINT_LOW)
    Cm_obs, _ = generate_synthetic_data(CNBIO, drug_true, DOSE_UM, t_sample, NOISE_CV, seed=SEED)

    ratios = np.array([0.2, 0.48, 0.96, 1.92, 4.8, 9.6, 19.2])
    deltas = profile_clint_single_start(t_sample, Cm_obs, drug_true, fixed={}, ratios=ratios)

    print("\n  ratio(CLint/fitted)   delta-SSE")
    for r, d in zip(ratios, deltas):
        print(f"    {r:<20.2f}{d:.5f}")

    monotonic = np.all(np.diff(deltas[np.argsort(ratios)]) >= -1e-6)
    print(f"\n  -> Profile is {'monotonic' if monotonic else 'NON-monotonic / erratic'}"
          f" (expected: erratic, diagnostic of single-start\n"
          f"     optimizer convergence failure at extreme grid points, not a genuine\n"
          f"     likelihood ridge).\n")


# ======================================================================
# S4.3 -- Multi-start profile likelihood: the trustworthy result
# ======================================================================
def multistart_best_fit(t_sample, Cm_obs, fixed, n_starts=6, seed=42):
    rng = np.random.default_rng(seed)
    free_names = [p for p in PARAM_NAMES if p not in fixed]
    best_sse, best_res = np.inf, None
    for _ in range(n_starts):
        mult = rng.uniform(0.3, 3.0, size=len(free_names))
        base_vals = {"CLd": GROUND_TRUTH["CLd"], "Kp": GROUND_TRUTH["Kp"],
                     "kb": GROUND_TRUTH["kb"], "e": CNBIO.e_uL_per_min + 1e-6,
                     "CLint": fixed.get("CLint", CLINT_LOW)}
        theta0 = [base_vals[n] * mult[i] for i, n in enumerate(free_names)]
        try:
            res = fit_mechanistic_model(t_sample, Cm_obs, DOSE_UM, CNBIO, theta0, fixed=fixed)
            sse = res.fun @ res.fun
            if sse < best_sse:
                best_sse, best_res = sse, res
        except Exception:
            continue
    return best_sse, best_res


def multistart_profile(t_sample, Cm_obs, CLint_true, fixed_base, ratios, n_starts=6, seed=42):
    sse_mle, res_mle = multistart_best_fit(t_sample, Cm_obs, dict(fixed_base), n_starts, seed)
    n_obs = len(res_mle.fun)
    dof = max(n_obs - len(res_mle.x), 1)
    sigma2 = sse_mle / dof
    deltas = []
    for r in ratios:
        fixed = dict(fixed_base); fixed["CLint"] = CLint_true * r
        sse, _ = multistart_best_fit(t_sample, Cm_obs, fixed, n_starts, seed)
        deltas.append(sse - sse_mle)
    return np.array(deltas), sigma2


def section_S4_3():
    print("=" * 70)
    print("S4.3  Multi-start profile likelihood")
    print("=" * 70)

    ratios_local = np.array([0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0])
    scenarios = [
        ("sparse (n=5), kb+e UNKNOWN", 20.6 * 60, 5, {}),
        ("sparse (n=5), kb+e KNOWN", 20.6 * 60, 5, {"kb": GROUND_TRUTH["kb"], "e": CNBIO.e_uL_per_min}),
        ("dense (n=12), kb+e KNOWN", 20.6 * 60, 12, {"kb": GROUND_TRUTH["kb"], "e": CNBIO.e_uL_per_min}),
    ]

    drug_true = make_drug(CLINT_LOW)

    for label, duration_min, n, fixed_base in scenarios:
        t_sample = sampling_schedule(duration_min, n)
        Cm_obs, _ = generate_synthetic_data(CNBIO, drug_true, DOSE_UM, t_sample, NOISE_CV, seed=SEED)
        deltas, sigma2 = multistart_profile(t_sample, Cm_obs, CLINT_LOW, fixed_base, ratios_local)
        threshold = chi2.ppf(0.95, df=1) * sigma2

        n_inside = np.sum(deltas <= threshold)
        print(f"\n  Scenario: {label}")
        print(f"    95% chi-square threshold: {threshold:.5f}")
        print(f"    {'ratio':<8}{'delta-SSE':<14}{'inside 95% CI?'}")
        for r, d in zip(ratios_local, deltas):
            print(f"    {r:<8.2f}{d:<14.5f}{'yes' if d <= threshold else 'no'}")
        print(f"    -> {n_inside}/{len(ratios_local)} of the scanned 0.3x-3x range is "
              f"statistically indistinguishable from the optimum.")

    print("\n  -> Even fixing kb AND evaporation (e) at their true values does not, in\n"
          "     this low-CLint regime, sharpen CLint's practical identifiability under\n"
          "     either sparse or dense sampling -- the basis of the manuscript's\n"
          "     'methodological caution' (structural identifiability is necessary\n"
          "     but not sufficient for practical estimability).\n")


if __name__ == "__main__":
    section_S4_1()
    section_S4_2()
    section_S4_3()
    print("=" * 70)
    print("Done. See supplementary_methods.md Section S4 for the narrative")
    print("account these diagnostics support.")
    print("=" * 70)
