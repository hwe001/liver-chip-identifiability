# What can media-only liver-on-a-chip concentration data identify?

*This repository previously hosted a shorter Perspective-format version of
this work submitted to CPT: Pharmacometrics & Systems Pharmacology, which
was not accepted. It has been replaced with the current, fuller version
prepared for the Journal of Pharmacokinetics and Pharmacodynamics; the
core structural-identifiability result and code are unchanged, but the
manuscript format, several figures, and the practical-identifiability
case study are new. The previous version remains available in this
repository's git history.*

Code, manuscript, and supplementary material for:

> **What can media-only liver-on-a-chip concentration data identify?
> Structural identifiability, a continuous equivalence family, and a
> cautionary case study in practical identifiability**
> (submitted to the *Journal of Pharmacokinetics and Pharmacodynamics*)

## Reproduce every figure in one command

```bash
pip install -r requirements.txt
python make_all_figures.py
```

This regenerates `figures/Fig1.pdf` through `Fig4.pdf` (and `.png` previews)
and prints, to stdout, the exact numerical values quoted in the manuscript's
Results section (condition numbers, relative standard errors, the
3.92e-07 maximum relative difference in the equivalence-family
demonstration, etc.) — so the script's own output is a direct
reproducibility check against the paper's text.

## Repository structure

```
main.tex, abstract.tex, introduction.tex, theoretical.tex,
methods.tex, results.tex, discussion.tex, conclusions.tex   - manuscript source (compile with latexmk -pdf main.tex)
main.pdf                                                     - compiled manuscript

supplementary.tex, supplementary.pdf                         - full proofs (Propositions/Theorem with proofs)
                                                                and the complete practical-identifiability
                                                                diagnostic case study (Sections S1-S5)

model.py                  - ODE right-hand side, platform hardware parameters,
                             media-concentration simulation
analysis.py                - synthetic data generation, naive/mechanistic fitting,
                             Fisher-information practical-identifiability diagnostics
structural_identifiability.py
                            - symbolic + numerical verification of the exact structural
                             non-identifiability result and the continuous equivalence
                             family (Theoretical section, Eq. 5; Supplementary S2)
profile_likelihood.py      - profile-likelihood scanning utilities
diagnostics_S4.py          - consolidated script reproducing every numbered claim in
                             Supplementary Section S4 (FIM instability, single-start
                             vs. multi-start profile likelihood)
run_scenarios.py           - shared parameter-naming constants and the adaptive-duration
                             helper used across the other scripts

figure1_schematic.py       - generates Figure 1 (three-panel schematic summary)
figures2to4.py             - generates Figures 2-4 (equivalence-family demonstration,
                             FIM instability, single- vs. multi-start profile likelihood)
make_all_figures.py        - runs both of the above in sequence; the single
                             entry point for full reproduction

phase_map.py                - EXPLORATORY, NOT a finished analysis (see "Status and known
                             gaps" below). Attempts to map practical identifiability across
                             the dimensionless (mu, beta) parameter space at fixed r =
                             Kp*Vc/Vm. Run separately: python phase_map.py <row_index>,
                             then python phase_map.py combine. Output saved to
                             data/phase_map_data.npz. Not referenced by any figure in the
                             manuscript; retained for transparency of the investigation
                             described in the Discussion's Limitations.

figures/                   - output directory for all generated figures
data/                      - output directory for saved numerical results (e.g. the
                             phase-map grid)
```

## What each figure shows

- **Figure 1** — Schematic summary: media concentration data determine three
  combinations of four mechanistic parameters (panel A); the resulting
  continuous equivalence family collapses to a point once non-specific
  binding (kb) is independently measured (panel B); structural
  identifiability, once achieved, is necessary but not sufficient for
  practical estimability (panel C).
- **Figure 2** — Two parameter quadruples differing up to 4-fold, related
  by the closed-form equivalence family (Eq. 5), produce media
  concentration trajectories agreeing to within solver precision
  (relative difference < 4e-7).
- **Figure 3** — Fisher-information-based practical-identifiability
  diagnostics are themselves numerically unreliable near the structural
  degeneracy: condition numbers of 10^12-10^17, and a plain-inverse vs.
  pseudoinverse disagreement of two to three orders of magnitude in the
  estimated relative standard error of CLint, on identical data.
- **Figure 4** — A single-start profile likelihood for CLint is
  non-monotonic and erratic (panel a), diagnostic of optimizer
  non-convergence; a multi-start profile likelihood (panel b) is stable,
  and shows that CLint remains weakly constrained across three sampling
  scenarios even after the structural degeneracy is resolved.

## Status and known gaps

- **The mu/beta phase map (`phase_map.py`) is an unfinished, exploratory
  investigation, not a validated result, and is not referenced by any
  figure in the manuscript.** A 2D sweep over (mu, beta) at fixed hardware
  (CN Bio, r = Kp*Vc/Vm ≈ 0.009) found the 95% profile-likelihood
  confidence interval for CLint spanning the full tested 100-fold range
  at every grid point — including regimes where Fisher-information-based
  diagnostics had suggested better identifiability (consistent with this
  paper's finding that FIM diagnostics are unreliable near this
  degeneracy). Follow-up diagnosis traced this to the fixed, small value
  of r, not to mu or beta: a zero-noise sensitivity check showed CLint
  sensitivity increases sharply with r (a 100-fold increase in Vc
  increased the fit-degradation sensitivity ratio by roughly eleven
  orders of magnitude). Attempting to confirm this under realistic noise
  by varying r directly did not work cleanly, because changing Vc also
  changes the system's intrinsic timescales, so the fixed sampling
  schedule used elsewhere in this repository becomes mismatched at
  different r and the comparison is confounded. This is reported
  honestly as an open question in the manuscript's Discussion
  (Limitations) and Conclusions, not as a settled finding. A correct
  treatment requires re-deriving an r-appropriate sampling schedule at
  each grid point — a reasonable next step for follow-up work, not
  completed here.
- Full author lists for two references (Milani et al. 2022; Pamies et al.
  2024) have not been independently verified against the primary source
  beyond the lead author; verify before submission.

## Requirements

See `requirements.txt`. Developed and tested with Python 3.11+, NumPy,
SciPy, SymPy, Matplotlib. LaTeX compilation requires a standard TeX
distribution (tested with `latexmk` + `pdflatex`).

## License

Add your preferred license before making this repository public.
