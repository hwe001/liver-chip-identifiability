import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from model import Platform, DrugParams, PLATFORMS, simulate_media_concentration
from analysis import (
    generate_synthetic_data, fit_naive_one_compartment, fit_mechanistic_model,
    practical_identifiability, clint_uL_per_min_per_1e6cells_to_human_CLH,
)

PARAM_NAMES = ["CLint", "CLd", "Kp", "kb", "e"]

# ---- Scenario grid ----
CLINT_LEVELS = {"low": 1.0, "high": 50.0}          # uL/min/1e6 cells (Rajan-anchored range)
SAMPLING = {
    "sparse": None,   # filled per-platform below (typical reported schedules)
    "dense": None,
}
NOISE_CV = 0.15
DOSE_UM = 1.0
GROUND_TRUTH = dict(CLd=8.0, Kp=5.0, kb=0.0008)  # shared mechanistic constants (illustrative)

N_REPLICATES = 15  # Monte Carlo replicates per scenario for bias/CI estimates


def sampling_schedule(platform_name, density, duration_min):
    if density == "sparse":
        n = 5
    else:
        n = 12
    # avoid t=0 (log undefined issues at exact zero for naive fit robustness)
    return np.linspace(duration_min / n, duration_min, n)


def adaptive_duration(platform, drug, dose_uM, max_duration_min, target_frac=0.15):
    """Pick an incubation window the way an experimentalist would: long
    enough to see substantial depletion, capped at a platform-realistic
    maximum (mirrors Rajan et al. 2023 using shorter windows for
    high-clearance and longer windows for low-clearance compounds)."""
    t_grid = np.geomspace(max_duration_min / 500, max_duration_min, 300)
    Cm = simulate_media_concentration(t_grid, dose_uM, platform, drug)
    frac = Cm / Cm[0]
    below = np.where(frac < target_frac)[0]
    if len(below) == 0:
        return max_duration_min
    return float(t_grid[below[0]])


def max_duration_for_platform(name):
    # platform-realistic ceilings (Rajan et al. 2023 ran up to 196 h;
    # CnBio/spheroid studies typically run up to ~4 days)
    return {"CnBio": 96 * 60, "Javelin": 196 * 60, "Spheroid": 96 * 60}[name]


def run_all():
    rows = []
    detail_rows = []

    for plat_name, platform in PLATFORMS.items():
        max_duration = max_duration_for_platform(plat_name)
        for clint_label, clint_per_1e6 in CLINT_LEVELS.items():
            CLint_total = clint_per_1e6 * (platform.cells / 1e6)
            drug_true = DrugParams(
                CLint_uL_per_min=CLint_total,
                CLd_uL_per_min=GROUND_TRUTH["CLd"],
                Kp=GROUND_TRUTH["Kp"],
                kb_per_min=GROUND_TRUTH["kb"],
            )
            duration = adaptive_duration(platform, drug_true, DOSE_UM, max_duration)
            for density in ["sparse", "dense"]:
                t_sample = sampling_schedule(plat_name, density, duration)

                naive_bias, mech_bias, mech_bias_kbknown = [], [], []
                rel_se_accum, rel_se_kbknown_accum = [], []
                cond_accum, cond_kbknown_accum = [], []

                for rep in range(N_REPLICATES):
                    Cm_obs, Cm_true = generate_synthetic_data(
                        platform, drug_true, DOSE_UM, t_sample, NOISE_CV, seed=1000 + rep
                    )

                    # naive one-compartment fit
                    CLc_naive = fit_naive_one_compartment(t_sample, Cm_obs, platform.Vm0_uL)
                    naive_bias.append(CLc_naive / CLint_total)

                    # mechanistic fit, ALL 5 parameters jointly free
                    theta0_full = [
                        CLint_total * 1.3, GROUND_TRUTH["CLd"] * 0.8, GROUND_TRUTH["Kp"] * 1.2,
                        GROUND_TRUTH["kb"] * 1.5, platform.e_uL_per_min * 1.2 + 1e-6,
                    ]
                    res = fit_mechanistic_model(t_sample, Cm_obs, DOSE_UM, platform, theta0_full)
                    mech_bias.append(res.x[0] / CLint_total)
                    rel_se, corr, cond_n = practical_identifiability(res, NOISE_CV)
                    rel_se_accum.append(rel_se)
                    cond_accum.append(cond_n)

                    # mechanistic fit with kb PRE-CHARACTERIZED (fixed at truth,
                    # mirroring a cell-free non-specific-binding control) --
                    # structural identifiability predicts this should resolve
                    # CLint, CLd and Kp uniquely (see structural_identifiability.py)
                    fixed = {"kb": GROUND_TRUTH["kb"]}
                    theta0_kbknown = [CLint_total * 1.3, GROUND_TRUTH["CLd"] * 0.8, GROUND_TRUTH["Kp"] * 1.2,
                                      platform.e_uL_per_min * 1.2 + 1e-6]
                    res_kb = fit_mechanistic_model(t_sample, Cm_obs, DOSE_UM, platform, theta0_kbknown, fixed=fixed)
                    mech_bias_kbknown.append(res_kb.x[0] / CLint_total)
                    rel_se_kb, _, cond_kb = practical_identifiability(res_kb, NOISE_CV)
                    rel_se_kbknown_accum.append(rel_se_kb)
                    cond_kbknown_accum.append(cond_kb)

                rel_se_mean = np.nanmean(np.array(rel_se_accum), axis=0)
                rel_se_kbknown_mean = np.nanmean(np.array(rel_se_kbknown_accum), axis=0)
                cond_median = np.nanmedian(cond_accum)
                cond_kbknown_median = np.nanmedian(cond_kbknown_accum)

                naive_bias = np.array(naive_bias)
                mech_bias = np.array(mech_bias)
                mech_bias_kbknown = np.array(mech_bias_kbknown)

                # translate CLint bias into human CLH prediction error
                true_CLH = clint_uL_per_min_per_1e6cells_to_human_CLH(clint_per_1e6)
                naive_CLH_fold_error = np.median(naive_bias)  # proportional bias carries through linearly-ish
                mech_CLH_fold_error = np.median(mech_bias)

                rows.append(dict(
                    platform=plat_name, CLint_level=clint_label, sampling=density,
                    n_samples=len(t_sample), duration_h=round(duration / 60, 1),
                    evaporation_source=platform.e_source[:40] + "...",
                    naive_bias_median=np.median(naive_bias),
                    naive_bias_cv=np.std(naive_bias) / np.median(naive_bias),
                    mech_bias_median=np.median(mech_bias),
                    mech_bias_cv=np.std(mech_bias) / np.median(mech_bias),
                    mech_kbknown_bias_median=np.median(mech_bias_kbknown),
                    mech_kbknown_bias_cv=np.std(mech_bias_kbknown) / np.median(mech_bias_kbknown),
                    CLint_relSE=rel_se_mean[0],
                    CLd_relSE=rel_se_mean[1],
                    Kp_relSE=rel_se_mean[2],
                    kb_relSE=rel_se_mean[3],
                    e_relSE=rel_se_mean[4],
                    CLint_relSE_kbknown=rel_se_kbknown_mean[0],
                    CLd_relSE_kbknown=rel_se_kbknown_mean[1],
                    Kp_relSE_kbknown=rel_se_kbknown_mean[2],
                    e_relSE_kbknown=rel_se_kbknown_mean[3],
                    cond_JTJ_median=cond_median,
                    cond_JTJ_kbknown_median=cond_kbknown_median,
                ))

                for nb, mb, mbk in zip(naive_bias, mech_bias, mech_bias_kbknown):
                    detail_rows.append(dict(
                        platform=plat_name, CLint_level=clint_label, sampling=density,
                        naive_bias=nb, mech_bias=mb, mech_bias_kbknown=mbk,
                    ))

    return pd.DataFrame(rows), pd.DataFrame(detail_rows)


def make_figures(summary_df, detail_df, outdir):
    plats = list(PLATFORMS.keys())
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
    for ax, plat in zip(axes, plats):
        sub = detail_df[detail_df.platform == plat]
        labels, data = [], []
        for clint_label in ["low", "high"]:
            for density in ["sparse", "dense"]:
                for approach, col in [("naive", "naive_bias"), ("mechanistic", "mech_bias")]:
                    s = sub[(sub.CLint_level == clint_label) & (sub.sampling == density)][col]
                    labels.append(f"{clint_label}\n{density}\n{approach}")
                    data.append(s.values)
        bp = ax.boxplot(data, labels=labels, showmeans=True)
        ax.axhline(1.0, color="green", linestyle="--", linewidth=1, label="unbiased")
        ax.axhline(3.0, color="red", linestyle=":", linewidth=1)
        ax.axhline(1/3, color="red", linestyle=":", linewidth=1, label="3-fold band")
        ax.set_title(plat)
        ax.set_yscale("log")
        ax.tick_params(axis="x", labelsize=7, rotation=0)
    axes[0].set_ylabel("Recovered CLint / True CLint")
    axes[0].legend(fontsize=8, loc="lower right")
    fig.suptitle("On-chip CLint recovery bias: naive one-compartment vs mechanistic fit\n"
                 "(synthetic data, 15% proportional noise, 15 Monte Carlo replicates per scenario)")
    fig.tight_layout()
    fig.savefig(f"{outdir}/fig1_clint_bias_by_platform.png", dpi=160)
    plt.close(fig)

    # Identifiability heatmap: log10(relative SE) of each parameter
    param_cols = ["CLint_relSE", "CLd_relSE", "Kp_relSE", "kb_relSE", "e_relSE"]
    fig2, axes2 = plt.subplots(1, 2, figsize=(12, 4.8))
    from matplotlib.colors import Normalize
    norm = Normalize(vmin=-1, vmax=4)  # log10(relSE): -1 (SE=10% of estimate) to 4 (SE=10000x estimate)
    for ax, density in zip(axes2, ["sparse", "dense"]):
        mat = []
        row_labels = []
        for plat in plats:
            for clint_label in ["low", "high"]:
                r = summary_df[(summary_df.platform == plat) &
                                (summary_df.CLint_level == clint_label) &
                                (summary_df.sampling == density)]
                if len(r):
                    vals = r[param_cols].values[0].astype(float)
                    vals = np.clip(vals, 1e-3, 1e4)
                    mat.append(np.log10(vals))
                    row_labels.append(f"{plat}-{clint_label}")
        mat = np.array(mat)
        im = ax.imshow(mat, aspect="auto", cmap="RdYlGn_r", norm=norm)
        ax.set_xticks(range(len(param_cols)))
        ax.set_xticklabels([c.replace("_relSE", "") for c in param_cols])
        ax.set_yticks(range(len(row_labels)))
        ax.set_yticklabels(row_labels, fontsize=8)
        ax.set_title(f"{density} sampling")
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                ax.text(j, i, f"{mat[i, j]:.1f}", ha="center", va="center", fontsize=7)
        cbar = fig2.colorbar(im, ax=ax, label="log10(relative SE)")
        cbar.ax.axhline((1 - norm.vmin) / (norm.vmax - norm.vmin), color="black", lw=1)
    fig2.suptitle("Practical identifiability: log10(relative standard error) of each fitted parameter\n"
                  "(green = well identified, SE << estimate; red = poorly identified, SE >> estimate; "
                  "black line on colorbar marks SE = estimate)")
    fig2.tight_layout()
    fig2.savefig(f"{outdir}/fig2_identifiability_heatmap.png", dpi=160)
    plt.close(fig2)


def make_kbknown_figure(summary_df, outdir):
    """Compare relative SE of CLint, CLd, Kp when kb is jointly fit
    (unknown) vs pre-characterized (known, e.g. from a cell-free NSB
    control) -- the empirical counterpart of the structural
    identifiability result in structural_identifiability.py."""
    plats = list(PLATFORMS.keys())
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, density in zip(axes, ["sparse", "dense"]):
        rows_unknown, rows_known, labels = [], [], []
        for plat in plats:
            for clint_label in ["low", "high"]:
                r = summary_df[(summary_df.platform == plat) &
                                (summary_df.CLint_level == clint_label) &
                                (summary_df.sampling == density)]
                if len(r):
                    rows_unknown.append(np.log10(np.clip(r["CLint_relSE"].values[0], 1e-3, 1e4)))
                    rows_known.append(np.log10(np.clip(r["CLint_relSE_kbknown"].values[0], 1e-3, 1e4)))
                    labels.append(f"{plat}-{clint_label}")
        x = np.arange(len(labels))
        w = 0.35
        ax.bar(x - w / 2, rows_unknown, width=w, label="kb jointly fit (unknown)", color="firebrick")
        ax.bar(x + w / 2, rows_known, width=w, label="kb pre-characterized (known)", color="seagreen")
        ax.axhline(0, color="black", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
        ax.set_ylabel("log10(relative SE of CLint)")
        ax.set_title(f"{density} sampling")
        ax.legend(fontsize=8)
    fig.suptitle("Effect of pre-characterizing non-specific binding (kb) on CLint identifiability\n"
                 "(structural result: media-only sampling cannot separate CLint, CLd, Kp, kb "
                 "unless kb is fixed independently)")
    fig.tight_layout()
    fig.savefig(f"{outdir}/fig3_kbknown_vs_unknown.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    import os
    outdir = "/mnt/user-data/outputs"
    os.makedirs(outdir, exist_ok=True)
    summary_df, detail_df = run_all()
    summary_df.to_csv(f"{outdir}/liver_chip_simulation_summary.csv", index=False)
    make_figures(summary_df, detail_df, outdir)
    make_kbknown_figure(summary_df, outdir)
    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 25)
    print(summary_df.round(3).to_string(index=False))
