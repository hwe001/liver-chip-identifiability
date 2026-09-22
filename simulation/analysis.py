"""
Synthetic-data generation, naive one-compartment fitting (mimicking the
conventional literature approach), mechanistic-model fitting, and
practical identifiability analysis (Fisher Information Matrix) for the
liver-chip model in model.py.

Also includes IVIVE scaling to human hepatic clearance (well-stirred
model), using the same scaling equations as DigiLoCS (Aravindakshan et
al. 2025), so bias in on-chip CLint can be translated into a clinically
interpretable error in predicted human clearance.
"""

import numpy as np
from scipy.optimize import least_squares
from model import Platform, DrugParams, simulate_media_concentration

RNG = np.random.default_rng(20260922)

# ---- IVIVE / well-stirred model constants (DigiLoCS, Eq. 14-15) ----
HC = 120e6      # human hepatocellularity, cells / g liver
LW = 25.7       # g liver / kg body weight
QH = 20.7       # mL/min/kg, average hepatic blood flow
FU_INC = 1.0    # fraction unbound in incubation medium (simplification: unbound drug)
FU_B = 1.0      # fraction unbound in blood (simplification)


def clint_uL_per_min_per_1e6cells_to_human_CLH(clint_per_1e6, fu_inc=FU_INC, fu_b=FU_B):
    """Scale an in-vitro intrinsic clearance (per million cells) to a
    predicted human hepatic clearance (mL/min/kg) via the well-stirred model."""
    CLint_H = (clint_per_1e6 * HC * LW) / fu_inc  # uL/min/kg -> convert below
    CLint_H_mL = CLint_H / 1000.0                  # uL -> mL
    CLH_pred = (QH * fu_b * CLint_H_mL) / (QH + fu_b * CLint_H_mL)
    return CLH_pred


def generate_synthetic_data(platform: Platform, drug: DrugParams, dose_uM,
                             t_sample, noise_cv, seed=None):
    """Simulate the true chip trajectory and add proportional (CV%) noise
    to mimic LC-MS/MS assay variability."""
    rng = RNG if seed is None else np.random.default_rng(seed)
    Cm_true = simulate_media_concentration(t_sample, dose_uM, platform, drug)
    Cm_obs = Cm_true * (1 + rng.normal(0, noise_cv, size=Cm_true.shape))
    Cm_obs = np.clip(Cm_obs, 1e-4 * dose_uM, None)  # LC-MS/MS lower limit of quantification, relative to dose
    return Cm_obs, Cm_true


def fit_naive_one_compartment(t_sample, Cm_obs, Vm0_uL):
    """Conventional literature approach: single well-mixed compartment,
    constant volume, log-linear regression -> CLc (uL/min)."""
    logC = np.log(Cm_obs)
    A = np.vstack([t_sample, np.ones_like(t_sample)]).T
    slope, intercept = np.linalg.lstsq(A, logC, rcond=None)[0]
    CLc = -slope * Vm0_uL  # uL/min
    return CLc


def _residuals(theta, t_sample, Cm_obs, dose_uM, platform, fixed=None):
    """fixed: optional dict of {param_name: value} to hold constant
    (e.g. {'kb': known_value} when NSB has been pre-characterized via
    a cell-free control, as in Docci et al. 2022 and Rajan et al. 2023).
    theta then contains only the remaining free parameters, in the
    fixed PARAM_NAMES order with the fixed ones skipped."""
    from run_scenarios import PARAM_NAMES
    fixed = fixed or {}
    free_names = [p for p in PARAM_NAMES if p not in fixed]
    if min(theta) <= 0:
        return np.full_like(Cm_obs, 1e3)
    values = dict(fixed)
    values.update(zip(free_names, theta))
    p = DrugParams(CLint_uL_per_min=values["CLint"], CLd_uL_per_min=values["CLd"],
                   Kp=values["Kp"], kb_per_min=values["kb"])
    plat = Platform(platform.name, platform.Vm0_uL, platform.Vc_uL, platform.cells,
                     values["e"], platform.e_source)
    try:
        Cm_pred = simulate_media_concentration(t_sample, dose_uM, plat, p)
    except Exception:
        return np.full_like(Cm_obs, 1e3)
    return (np.log(Cm_obs) - np.log(np.clip(Cm_pred, 1e-9, None)))


def fit_mechanistic_model(t_sample, Cm_obs, dose_uM, platform, theta0, fixed=None):
    """Fit the free parameters (all 5, or fewer if `fixed` pre-specifies
    some, e.g. kb from a cell-free NSB control) by nonlinear least
    squares on log-concentration residuals. Returns fitted params +
    Jacobian at the optimum (for the FIM / practical identifiability
    analysis). theta0 must match the free-parameter order/length."""
    from run_scenarios import PARAM_NAMES
    fixed = fixed or {}
    free_names = [p for p in PARAM_NAMES if p not in fixed]
    bounds_lb = {"CLint": 1e-4, "CLd": 1e-4, "Kp": 1e-2, "kb": 1e-6, "e": 0.0}
    bounds_ub = {"CLint": 1e5, "CLd": 1e5, "Kp": 1e3, "kb": 1.0, "e": platform.Vm0_uL / t_sample[-1] * 0.9}
    lb = [bounds_lb[p] for p in free_names]
    ub = [bounds_ub[p] for p in free_names]
    res = least_squares(
        _residuals, theta0, args=(t_sample, Cm_obs, dose_uM, platform, fixed),
        bounds=(lb, ub), method="trf", xtol=1e-12, ftol=1e-12,
    )
    res.free_names = free_names
    return res


def practical_identifiability(res, noise_cv, rcond=1e-10):
    """Fisher Information Matrix (Gauss-Newton approx) -> covariance ->
    relative standard errors and parameter correlation matrix.

    IMPORTANT: near a structural (or near-structural) non-identifiability,
    J^T J is numerically singular (condition numbers of 1e12-1e17 were
    observed in this model). Plain matrix inversion (np.linalg.inv) on
    such a matrix is numerically unstable and can disagree with a
    pseudoinverse by two orders of magnitude on the SAME data -- i.e. the
    instability is a property of the inversion algorithm, not of the
    underlying parameter uncertainty. We therefore use a truncated SVD
    pseudoinverse (np.linalg.pinv with an explicit rcond) throughout, and
    report the condition number alongside every relSE so a near-singular
    result is visible rather than silently amplified."""
    J = res.jac  # d(residual)/d(theta), residuals already in log space
    n = len(res.fun)
    dof = max(n - len(res.x), 1)
    sigma2 = (res.fun @ res.fun) / dof
    JTJ = J.T @ J
    cond_number = np.linalg.cond(JTJ) if J.size else np.inf
    try:
        cov = sigma2 * np.linalg.pinv(JTJ, rcond=rcond)
        se = np.sqrt(np.abs(np.diag(cov)))
        rel_se = se / np.abs(res.x)
        with np.errstate(invalid="ignore", divide="ignore"):
            corr = cov / np.outer(se, se)
    except np.linalg.LinAlgError:
        rel_se = np.full(len(res.x), np.inf)
        corr = np.full((len(res.x), len(res.x)), np.nan)
    return rel_se, corr, cond_number
