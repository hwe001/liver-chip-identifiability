"""
Profile-likelihood practical identifiability for CLint.

Motivation: a Fisher-Information-Matrix (FIM) / Gauss-Newton covariance
approximation was found to be numerically unreliable in this model near
the structural degeneracy identified in structural_identifiability.py.
J^T J condition numbers of 1e9-1e17 were observed; plain matrix inversion
and a truncated-SVD pseudoinverse disagreed on the resulting relative
standard errors by up to ~200-fold on IDENTICAL data (see
diagnostic output referenced in the manuscript). A direct profile-
likelihood check (below) showed the pseudoinverse-based numbers were
misleadingly optimistic: for a representative low-clearance, sparse-
sampling, kb-unknown scenario, the profile is essentially FLAT (delta-SSE
~1e-3) across a 100-fold range of CLint, i.e. genuinely unidentifiable in
practice -- despite a pseudoinverse relative SE of ~1 suggesting good
identification.

Profile likelihood does not require inverting anything and is therefore
immune to this failure mode: for each candidate CLint, all remaining free
parameters are re-optimized, and the resulting best-attainable SSE is
compared to the global MLE. A 95% confidence interval on CLint is the set
of values where SSE(CLint) - SSE(MLE) <= chi2.ppf(0.95, df=1) * sigma^2.
"""

import numpy as np
from scipy.stats import chi2
from model import Platform, DrugParams
from analysis import generate_synthetic_data, fit_mechanistic_model


def profile_clint(platform, drug_true, dose_uM, t_sample, Cm_obs, fixed_base=None,
                   ratios=None, theta0_builder=None):
    """Return (ratios, delta_sse array, mle_ratio, sigma2) for a profile
    likelihood scan of CLint at the given fixed ratios relative to the
    true CLint. fixed_base: dict of any other parameters to hold fixed
    (e.g. {'kb': known_value}) throughout the whole profile, including
    at the MLE itself, for a fair like-for-like comparison."""
    fixed_base = dict(fixed_base or {})
    ratios = ratios if ratios is not None else np.geomspace(0.05, 30, 25)
    CLint_true = drug_true.CLint_uL_per_min

    # MLE with CLint free (but fixed_base still applied, e.g. kb known)
    from run_scenarios import PARAM_NAMES
    free_names_mle = [p for p in PARAM_NAMES if p not in fixed_base]
    theta0_mle = theta0_builder(free_names_mle, CLint_true)
    res_mle = fit_mechanistic_model(t_sample, Cm_obs, dose_uM, platform, theta0_mle, fixed=fixed_base)
    sse_mle = res_mle.fun @ res_mle.fun
    n = len(res_mle.fun)
    dof = max(n - len(res_mle.x), 1)
    sigma2 = sse_mle / dof
    mle_ratio = res_mle.x[free_names_mle.index("CLint")] / CLint_true if "CLint" in free_names_mle else 1.0

    delta_sse = []
    for r in ratios:
        clint_fixed = CLint_true * r
        fixed = dict(fixed_base)
        fixed["CLint"] = clint_fixed
        free_names = [p for p in PARAM_NAMES if p not in fixed]
        theta0 = theta0_builder(free_names, clint_fixed)
        res = fit_mechanistic_model(t_sample, Cm_obs, dose_uM, platform, theta0, fixed=fixed)
        sse = res.fun @ res.fun
        delta_sse.append(sse - sse_mle)
    return ratios, np.array(delta_sse), mle_ratio, sigma2


def ci_width_from_profile(ratios, delta_sse, sigma2, conf=0.95):
    """95% profile-likelihood CI width (in log10(ratio) units) -- the
    range of CLint/CLint_true ratios where delta_sse stays below the
    chi-square threshold. Returns None if the CI extends beyond the
    scanned range (i.e. is at least this wide -- a lower bound)."""
    threshold = chi2.ppf(conf, df=1) * sigma2
    inside = ratios[delta_sse <= threshold]
    if len(inside) == 0:
        return 0.0, True  # nothing inside -- immediately rejects; CI essentially a point
    lo, hi = inside.min(), inside.max()
    hit_boundary = (lo == ratios.min()) or (hi == ratios.max())
    return np.log10(hi / lo), hit_boundary
