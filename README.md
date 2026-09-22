# What can media-only liver-on-a-chip concentration data identify?

Code and manuscript supporting the CPT:PSP Perspective:

> **What can media-only liver-on-a-chip concentration data identify? More sampling cannot resolve structural ambiguity**

## Repository structure

```
manuscript/
    manuscript.md              - Main text (Introduction, Structural identifiability,
                                  Exact bias limits, Practical identifiability caution,
                                  Reporting standard, Limitations, Conclusion, Table 1,
                                  Figure legend, references)
    supplementary_methods.md   - Full derivations, diagnostic details (FIM condition
                                  numbers, profile-likelihood methodology), and scope notes
    plain_language_summary.md  - Non-technical summary of the paper's argument

simulation/
    model.py                   - Two-compartment ODE model, platform hardware parameters
    analysis.py                - Synthetic data generation, naive/mechanistic fitting,
                                  Fisher-information practical-identifiability diagnostics
    structural_identifiability.py
                                - Symbolic + numerical verification of the exact
                                  structural non-identifiability result (Section:
                                  "Structural identifiability of media-only sampling")
    profile_likelihood.py      - Profile-likelihood scanning utilities used for the
                                  practical-identifiability investigation
    run_scenarios.py           - Driver script for the exploratory cross-platform
                                  synthetic simulations (superseded in the final
                                  manuscript -- see note below)

figures/
    build_figure1.py           - Generates Figure 1 (three-panel schematic) at the
                                  journal's required dimensions (178mm, 600dpi)
    figure1_schematic.pdf/png  - Rendered Figure 1

deconvolution_exploration/
    tv_model.py                - Time-varying-CLint extension and quasi-steady-state
                                  recovery method (see note below)
    validate_drift_recovery.py - Synthetic validation of drift detection
    deconv_fig1_drift_recovery.png
```

## Important note on scope

The **final manuscript's** central claims rest entirely on:
- `simulation/structural_identifiability.py` (the exact structural non-identifiability
  proof and continuous equivalence family)
- The exact asymptotic bias derivations and nondimensionalization (symbolic/analytic;
  see `manuscript/supplementary_methods.md` Sections S1-S3)
- `simulation/profile_likelihood.py` and the multi-start profile-likelihood diagnostic
  described in `supplementary_methods.md` Section S4

`simulation/run_scenarios.py` and `simulation/analysis.py`'s cross-platform synthetic
comparisons were part of the manuscript's development process but were **removed from
the submitted paper** after multi-start profile likelihood revealed the earlier
FIM-based practical-identifiability claims were unreliable near the structural
degeneracy (see `supplementary_methods.md` Section S4 for the full account). This code
is retained here for transparency and reproducibility of that investigation, not
because its quantitative outputs are cited in the final manuscript.

The `deconvolution_exploration/` folder is a **separate, exploratory follow-on idea**
(testing whether CLint is constant over an incubation, using a quasi-steady-state
recovery method) developed after the main manuscript was finalized. It is not part of
the submitted paper and has not been through the same verification process — see the
docstring in `tv_model.py` for its known limitations (validity restricted to the
fast-distribution regime; no correction for the systematic bias documented there).

## Requirements

See `requirements.txt`. Developed and tested with Python 3.11+, NumPy 2.4, SciPy 1.17,
SymPy, Matplotlib.

## Reproducing the core result

```bash
cd simulation
python3 structural_identifiability.py
```

This reproduces the symbolic derivation, the closed-form continuous equivalence
family, and the numerical verification (two structurally-degenerate parameter sets
producing media trajectories agreeing to solver precision).

## Status

Manuscript version corresponds to `manuscript.md` as of the last update in this
repository. Two reference entries (Milani et al. 2022 and Pamies et al. 2024) have
truncated author lists pending final verification against primary sources — see the
"Still to verify before submission" note at the end of `manuscript.md`.

## License

See `LICENSE`. (Code: MIT license suggested as a permissive default -- replace with
your preferred license before making the repository public if you'd rather use
something else.)
