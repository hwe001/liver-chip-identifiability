"""
phase_map.py
Maps practical identifiability of CLint across the dimensionless parameter
space (mu = Kp*CLint/CLd, beta = kb*Vm/CLd), with kb and e fixed at their
true values (i.e. the structural degeneracy already resolved via a correct
NSB control) -- directly answering the reviewer's request to show *when*
media-only sampling is practically sufficient, not only when it fails.

Hardware held fixed at CN Bio values (Vm=1600 uL, Vc=3 uL); Kp=5.0 and
CLd=8.0 uL/min held fixed as the reference scale, so mu and beta are varied
by choosing CLint = mu*CLd/Kp and kb = beta*CLd/Vm. Dense sampling (n=12)
only, to keep runtime bounded; this is noted explicitly as a scope
limitation in the manuscript.
"""
import numpy as np
from scipy.stats import chi2

from model import CNBIO, DrugParams
from analysis import generate_synthetic_data
from run_scenarios import adaptive_duration
import diagnostics_S4 as d4

CLD, KP, VM = 8.0, 5.0, CNBIO.Vm0_uL
MU_GRID = np.array([0.03, 0.1, 0.3, 1.0, 3.0])
BETA_GRID = np.array([0.01, 0.03, 0.1, 0.3, 1.0])
N_RESTARTS = 8
RATIOS = np.array([0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0])
N_SAMPLES = 12
SEED = 1000


def ci_width(deltas, ratios, sigma2):
    threshold = chi2.ppf(0.95, df=1) * sigma2
    rejected = deltas > threshold
    # An isolated rejection surrounded by accepted points on both sides is
    # inconsistent with a well-behaved (unimodal, monotonically increasing
    # away from the optimum) profile, and is the known multi-start
    # convergence-failure artifact documented in Fig. 4a/S4.2. Reclassify
    # such isolated points as accepted, and flag them.
    spurious = np.zeros_like(rejected)
    for k in range(1, len(rejected) - 1):
        if rejected[k] and not rejected[k - 1] and not rejected[k + 1]:
            spurious[k] = True
    cleaned_rejected = rejected & ~spurious
    inside = ratios[~cleaned_rejected]
    if len(inside) == 0:
        return 0.0, False, int(spurious.sum())
    lo, hi = inside.min(), inside.max()
    hit_edge = (lo == ratios.min()) or (hi == ratios.max())
    return float(np.log10(hi / lo)), hit_edge, int(spurious.sum())


def run_row(i):
    """Compute one row (fixed mu, all beta) and save it to data/row_{i}.npz.
    Run as: python phase_map.py <row_index>  -- keeps each invocation short
    enough to finish well within a single tool call's time budget."""
    mu = MU_GRID[i]
    widths, hit_edges, spurious = [], [], []
    for beta in BETA_GRID:
        CLint = mu * CLD / KP
        kb = beta * CLD / VM
        drug_true = DrugParams(CLint_uL_per_min=CLint, CLd_uL_per_min=CLD, Kp=KP, kb_per_min=kb)
        duration = adaptive_duration(CNBIO, drug_true, 1.0, max_duration_min=96 * 60, target_frac=0.15)
        t_sample = np.linspace(duration / N_SAMPLES, duration, N_SAMPLES)
        Cm_obs, _ = generate_synthetic_data(CNBIO, drug_true, 1.0, t_sample, 0.15, seed=SEED)
        fixed_base = {"kb": kb, "e": CNBIO.e_uL_per_min}
        deltas, sigma2 = d4.multistart_profile(t_sample, Cm_obs, CLint, fixed_base, RATIOS,
                                                n_starts=N_RESTARTS, seed=42)
        width, hit_edge, n_spurious = ci_width(deltas, RATIOS, sigma2)
        widths.append(width); hit_edges.append(hit_edge); spurious.append(n_spurious)
        print(f"mu={mu:<6.3g} beta={beta:<6.3g}  CI width (log10)={width:.3f}  "
              f"{'(reaches scan edge)' if hit_edge else ''}"
              f"{f'  [{n_spurious} spurious rejection(s) reclassified]' if n_spurious else ''}")
    import os
    os.makedirs("data", exist_ok=True)
    np.savez(f"data/row_{i}.npz", mu=mu, beta=BETA_GRID, width=np.array(widths),
              hit_edge=np.array(hit_edges), spurious=np.array(spurious), ratios=RATIOS)
    print(f"saved data/row_{i}.npz")


def combine_rows():
    """Combine all data/row_*.npz into the final data/phase_map_data.npz."""
    results = np.full((len(MU_GRID), len(BETA_GRID)), np.nan)
    hit_edge_map = np.zeros_like(results, dtype=bool)
    spurious_count = np.zeros_like(results, dtype=int)
    for i in range(len(MU_GRID)):
        d = np.load(f"data/row_{i}.npz")
        results[i, :] = d["width"]
        hit_edge_map[i, :] = d["hit_edge"]
        spurious_count[i, :] = d["spurious"]
    np.savez("data/phase_map_data.npz", mu=MU_GRID, beta=BETA_GRID, width=results,
              hit_edge=hit_edge_map, spurious=spurious_count, ratios=RATIOS)
    print("saved data/phase_map_data.npz (combined)")
    return results, hit_edge_map, spurious_count


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "combine":
        combine_rows()
    elif len(sys.argv) > 1:
        run_row(int(sys.argv[1]))
    else:
        for i in range(len(MU_GRID)):
            run_row(i)
        combine_rows()
