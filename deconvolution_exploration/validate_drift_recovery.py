import sys
sys.path.insert(0, '..')
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from tv_model import simulate_tv_media_concentration, recover_instantaneous_CLint
from model import CNBIO

Vm, Vc, Kp, kb = CNBIO.Vm0_uL, CNBIO.Vc_uL, 5.0, 0.0008
CLd = 8.0
NOISE_CV = 0.15
DURATION_MIN = 60 * 60  # 60 h -- long incubation, where drift is most plausible
N_SAMPLES = 30
RNG = np.random.default_rng(2026)

SCENARIOS = {
    "flat (control)": lambda t: 0.5 * np.ones_like(t),
    "linear decline (50% loss over 60h)": lambda t: 0.5 * (1 - 0.5 * t / DURATION_MIN),
    "step change at 24h": lambda t: np.where(t < 24*60, 0.6, 0.3),
}


def run_scenario(name, CLint_func, n_reps=20):
    t = np.linspace(1, DURATION_MIN, N_SAMPLES)
    Cm_true = simulate_tv_media_concentration(t, 1.0, CNBIO, CLint_func, CLd, Kp, kb)
    true_vals = CLint_func(t)

    recovered_reps = []
    slopes = []
    for rep in range(n_reps):
        noise = RNG.normal(0, NOISE_CV, size=Cm_true.shape)
        Cm_obs = np.clip(Cm_true * (1 + noise), 1e-6, None)
        CLint_rec, _ = recover_instantaneous_CLint(t, Cm_obs, Vm, Vc, Kp, kb,
                                                      smooth_window=11, smooth_poly=2)
        recovered_reps.append(CLint_rec)
        # simple linear-trend test on the recovered profile (excluding edge artifacts)
        interior = slice(3, -3)
        slope, intercept, r, p, se = stats.linregress(t[interior], CLint_rec[interior])
        slopes.append(slope)

    recovered_reps = np.array(recovered_reps)
    slopes = np.array(slopes)
    # detection test: is zero outside the 95% CI of the recovered slope distribution?
    slope_mean, slope_se = slopes.mean(), slopes.std(ddof=1) / np.sqrt(n_reps)
    ci_lo, ci_hi = slope_mean - 1.96*slope_se, slope_mean + 1.96*slope_se
    detected = not (ci_lo <= 0 <= ci_hi)

    return dict(t=t, true_vals=true_vals, recovered_median=np.median(recovered_reps, axis=0),
                recovered_lo=np.percentile(recovered_reps, 5, axis=0),
                recovered_hi=np.percentile(recovered_reps, 95, axis=0),
                slope_mean=slope_mean, slope_ci=(ci_lo, ci_hi), detected=detected)


if __name__ == "__main__":
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
    print(f"{'Scenario':<38}{'slope (per min)':<20}{'95% CI':<28}{'drift detected?'}")
    for ax, (name, func) in zip(axes, SCENARIOS.items()):
        res = run_scenario(name, func)
        ax.plot(res["t"]/60, res["true_vals"], "k--", lw=1.5, label="true CLint(t)")
        ax.plot(res["t"]/60, res["recovered_median"], color="#b00000", lw=1.8, label="recovered (median)")
        ax.fill_between(res["t"]/60, res["recovered_lo"], res["recovered_hi"],
                          color="#b00000", alpha=0.2, label="5-95th pct (20 reps)")
        ax.set_title(name, fontsize=9)
        ax.set_xlabel("time (h)")
        ax.legend(fontsize=7)
        ax.set_ylim(0, 1.0)
        print(f"{name:<38}{res['slope_mean']:<20.6f}"
              f"[{res['slope_ci'][0]:.6f}, {res['slope_ci'][1]:.6f}]".ljust(28) +
              f"{res['detected']}")
    axes[0].set_ylabel("CLint (uL/min)")
    fig.suptitle("Recovering time-varying CLint(t) from noisy media concentration data\n"
                 "(quasi-steady-state method, 15% noise, 30 samples over 60h, 20 Monte Carlo reps)")
    fig.tight_layout()
    fig.savefig("/mnt/user-data/outputs/deconv_fig1_drift_recovery.png", dpi=160)
    print("\nsaved deconv_fig1_drift_recovery.png")
