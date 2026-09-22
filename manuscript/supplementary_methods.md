# Supplementary Methods

## S1. Model equations

Two-compartment model: media compartment (concentration C_m, volume V_m) and a lumped cell-associated compartment (concentration C_c, volume V_c), linked by distributional clearance CL_d and cell:media partition coefficient K_p, first-order non-specific binding (NSB) loss from media at rate k_b, and metabolism (intrinsic clearance CLint) confined to the cell compartment:

dC_m/dt = −(CL_d/V_m)C_m + (CL_d/(K_pV_m))C_c − k_bC_m
dC_c/dt = (CL_d/V_c)C_m − (CL_d/(K_pV_c) + CLint/V_c)C_c

Matrix form dx/dt = Ax with x=[C_m, C_c]ᵀ and

A = [[a11, a12], [a21, a22]],  a11=−(CL_d/V_m+k_b),  a12=CL_d/(K_pV_m),  a21=CL_d/V_c,  a22=−(CL_d/(K_pV_c)+CLint/V_c)

Media volume with evaporation: V_m(t) = max(V_m0 − e·t, floor), where e is a first-order-in-time (linear) evaporation rate. Evaporation is treated as a known, externally measured input to the model (via residual-volume weighing), not fitted jointly with the kinetic parameters, throughout the main-text structural identifiability derivation (Section S2). The full simulation code (model.py) allows e to be either fixed or fitted; the distinction matters for the practical-identifiability investigation in Section S5.

Initial condition: dose in media only, C_m(0)=C_0, C_c(0)=0.

## S2. Laplace-transform derivation of structural identifiability (media-only observation)

Taking the Laplace transform of the linear system above with y(t)=C_m(t) observed and x(0)=[C_0,0]ᵀ:

Y(s)/C_0 = (s−a22) / (s²−(a11+a22)s+(a11a22−a12a21))

The denominator and numerator coefficients are the only quantities determined by data of arbitrary density and duration:
- Numerator constant term: −a22
- Denominator s¹ coefficient: −(a11+a22), hence a11 (given a22)
- Denominator s⁰ coefficient: a11a22−a12a21, hence a12a21 ≡ P (given a11, a22)

So exactly three real numbers — a11, a22, P=a12a21=CL_d²/(K_pV_mV_c) — are identifiable from C_m(t), regardless of sampling. With V_m, V_c known (fixed hardware), this gives 3 equations in the 4 unknowns (CL_d, K_p, CLint, k_b): a11=−(CL_d/V_m+k_b), a22=−(CL_d/(K_pV_c)+CLint/V_c), P=CL_d²/(K_pV_mV_c).

### S2.1 Continuous equivalence family

Solving these three equations for (CL_d, K_p, CLint) in terms of a free parameter k_b:

CL_d(k_b) = V_m(−a11−k_b)
K_p(k_b) = CL_d(k_b)² / (P·V_m·V_c)
CLint(k_b) = −a22·V_c − CL_d(k_b)/K_p(k_b)

Any k_b in the admissible range (yielding positive CL_d, K_p, CLint) generates a full parameter quadruple consistent with the same (a11, a22, P), hence the same C_m(t) trajectory exactly. This was verified two ways:
1. **Symbolically** (SymPy): substituting the closed-form solution back into the definitions of a11, a22, P recovers the assumed values identically for arbitrary k_b, confirmed via `sympy.solve` on the three-equation system with numeric coefficients (see `structural_identifiability.py`).
2. **Numerically, at the ODE level**: two parameter sets differing by design — (CLint=4.5, CL_d=8.0, K_p=5.0, k_b=0.0008) and (CLint=3.243, CL_d=4.48, K_p=1.568, k_b=0.003), both consistent with the same (a11,a22,P) by construction — were each simulated through the full nonlinear ODE solver (constant-volume case, e=0) and their C_m(t) trajectories compared. Maximum absolute difference: 3.76×10⁻⁸; maximum relative difference: 3.92×10⁻⁷ (attributable to `solve_ivp` numerical tolerance, not a real difference — the trajectories are analytically identical).

### S2.2 Resolution via independent k_b

If k_b is measured independently (e.g., a cell-free NSB control), the three equations above become 3 equations in the 3 remaining unknowns (CL_d, K_p, CLint), with a unique solution given by the same closed-form expressions evaluated at the known k_b. This was confirmed by direct substitution and by refitting synthetic data with k_b held fixed (Section S5).

## S3. Exact asymptotic bias limits

Define the conventional (naive) apparent clearance as CL_naive = V_m·k_dep, where k_dep is the monoexponential media-depletion rate constant obtained by log-linear regression of C_m(t) assuming a single well-mixed compartment at constant V_m. In the true two-compartment system, the "terminal" decay rate an analyst would fit corresponds to the eigenvalue of A with smaller magnitude (λ_slow), so CL_naive = −λ_slow·V_m.

### S3.1 Fast-distribution limit (CL_d/CLint → ∞)

Substituting CL_d → 1/ε and taking ε→0 (equivalently, examining λ_slow as CL_d grows without bound with CLint, K_p, k_b, V_m, V_c fixed), confirmed numerically via direct eigenvalue evaluation across CL_d ∈ {8, 80, 800, ..., 8×10¹²} (Table S1): λ_slow converges to

λ_slow → −(K_pCLint+k_bV_m) / (V_m+K_pV_c)

so CL_naive → V_m(K_pCLint+k_bV_m)/(V_m+K_pV_c). Confirmed to 5 significant figures by CL_d=8×10⁸ (23.5585 vs. asymptotic 23.5591; floating-point noise appears by CL_d=8×10¹²).

**Table S1.** CL_naive vs. CL_d (CLint=4.5, K_p=5.0, k_b=0.0008, V_m=1600, V_c=3.0; asymptotic value 23.5591)

| CL_d | CL_naive |
|---|---|
| 8 | 7.177 |
| 80 | 18.734 |
| 800 | 22.961 |
| 8,000 | 23.498 |
| 8×10⁵ | 23.5585 |
| 8×10⁸ | 23.5591 |

### S3.2 Fast-metabolism limit (CLint/CL_d → ∞)

As CLint → ∞ with other parameters fixed, the cell compartment approaches a quasi-steady-state sink (C_c → 0 relative to C_m), and λ_slow → −(CL_d/V_m + k_b), giving CL_naive → CL_d + k_bV_m, independent of CLint. Confirmed numerically across CLint ∈ {4.5, 45, ..., 4.5×10¹⁰} (CL_d=8.0 fixed): CL_naive converges to 9.28 (=8.0+1600×0.0008) by CLint=4500, matching to 4 significant figures.

### S3.3 Nondimensionalization

With τ=(CL_d/V_m)t, x=C_m/C_m(0), y=C_c/(K_pC_m(0)):

dx/dτ = −(1+β)x + y
dy/dτ = (1/r)[x − (1+μ)y]

where r=K_pV_c/V_m, β=k_bV_m/CL_d, μ=K_pCLint/CL_d. Derived and verified symbolically (SymPy, chain-rule substitution of C_m=C_m0·x(τ), C_c=K_pC_m0·y(τ), t=(V_m/CL_d)τ into the original ODEs; see inline derivation in the analysis notebook). Rewriting the two asymptotic bias-ratio limits in these dimensionless terms:

CL_naive/CLint → [K_p/(1+r)](1+β/μ)  as CL_d/CLint→∞ (μ→0 along this physical path)
CL_naive/CLint → (K_p/μ)(1+β)  as CLint/CL_d→∞ (μ→∞ along this physical path)

Note μ→0 or μ→∞ alone does not uniquely specify which physical limit is meant (μ=K_pCLint/CL_d can approach these limits via CL_d→∞, CLint→0, or K_p→0/∞ separately); the two boxed asymptotic formulas above are specifically derived along the CL_d→∞ and CLint→∞ paths respectively, with the other parameters held fixed, and should not be read as general statements about arbitrary routes to μ→0 or μ→∞.

## S4. Practical identifiability: a methodological case study

This section documents, in full, the sequence of diagnostic attempts summarized qualitatively in the main text.

### S4.1 FIM covariance approximation: numerical instability near the structural degeneracy

For a representative scenario (CN Bio hardware, V_m=1600 µL, V_c=3 µL, low-CLint regime, sparse [n=5] and dense [n=12] sampling, 15% proportional noise), nonlinear least-squares fits (`scipy.optimize.least_squares`, trust-region-reflective) were performed with all 5 parameters (CLint, CL_d, K_p, k_b, e) free. The Gauss-Newton Fisher Information Matrix J^TJ was computed at each fitted optimum.

Observed condition numbers of J^TJ: 1.377×10¹⁷ (low-CLint, dense sampling) and 1.447×10¹² (high-CLint, dense sampling) — both indicating numerical singularity at or beyond double-precision limits. Comparing plain matrix inversion (`np.linalg.inv`) against a truncated-SVD pseudoinverse (`np.linalg.pinv`, rcond=1e-10) on the identical fitted Jacobian: relative standard errors for CLint of 1.404×10³ (plain inverse) vs. 6.885 (pseudoinverse) in one scenario — a ~200-fold disagreement from the choice of inversion algorithm alone, on identical data and an identical point estimate.

### S4.2 Single-start profile likelihood: contaminated by convergence failure

A profile likelihood for CLint (fixing CLint at a grid of values relative to its fitted value, re-optimizing all other free parameters, recording the resulting SSE) was computed with a single optimizer start per grid point. The resulting profile was non-monotonic and erratic: near-zero delta-SSE (as good as the unconstrained optimum) recurred at extreme CLint ratios (e.g., 0.05× and 23× the fitted value) interspersed with much larger delta-SSE at intermediate ratios — a pattern inconsistent with a well-behaved likelihood profile and diagnostic of the optimizer failing to locate the true conditional optimum at those grid points (poor local-minimum convergence from a fixed, non-adaptive starting guess).

### S4.3 Multi-start profile likelihood

Repeating the profile at each grid point with 6 random restarts (starting multipliers drawn uniformly from [0.3, 3.0]× the nominal parameter values) and retaining the best (lowest-SSE) fit at each point gave materially more stable profiles. Two findings:

1. **kb+e unknown, sparse sampling (n=5):** delta-SSE remained near zero (all points "inside" a 95% profile-likelihood confidence region, chi-square threshold with 1 d.f.) across the full scanned range of CLint ratios (0.3×–3.0×), i.e. CLint is not practically constrained by this design.
2. **kb+e known (fixed at true values), sparse AND dense (n=5, n=12) sampling:** delta-SSE also remained near zero across the same 0.3×–3.0× range in both cases — i.e., fixing k_b and e (removing the structural degeneracy) did not, in this regime, materially sharpen the practical identifiability of CLint, even under dense sampling.

This is the basis for the main text's "methodological caution": resolving the structural degeneracy (S2.2) does not by itself guarantee practical estimability under realistic finite, noisy sampling designs, in the parameter regime examined (low CLint relative to CL_d, i.e. small μ).

### S4.4 Scope and caveats of the practical-identifiability investigation

- All numerical values above are for one illustrative parameter regime (CN Bio hardware, GROUND_TRUTH CL_d=8.0, K_p=5.0, k_b=0.0008, low-CLint level) and should not be read as general claims about all regimes or platforms.
- The multi-start profile likelihood used 6 restarts per grid point with uniform-multiplier initialization; a more exhaustive search (more restarts, alternative initialization strategies, or a dedicated global optimizer) might further change the profile shape, particularly at the extremes of the scanned range.
- We did not perform a DAISY-based or other formal structural identifiability check on the complete 5-parameter, time-varying-volume model (only the reduced constant-volume 4-parameter model in S2); it is possible additional structural degeneracies exist between k_b and e in the full model, which would be a natural extension of this work.

## S5. Code and reproducibility

All derivations, symbolic verification, and numerical checks were performed in Python 3 (NumPy 2.4, SciPy 1.17, SymPy). Code files (to be deposited in a public repository prior to submission, per journal policy):

- `model.py` — ODE right-hand side, platform hardware parameters, media-concentration simulation
- `analysis.py` — synthetic data generation, naive and mechanistic fitting, FIM-based practical identifiability (plain inverse and pseudoinverse)
- `structural_identifiability.py` — symbolic and numerical verification of the equivalence family (Section S2)
- `profile_likelihood.py` — profile-likelihood scanning and confidence-interval utilities (Section S4)
- Figure-generation script for Figure 1 (`build_figure1.py`)

**Not yet done:** consolidating the ad hoc diagnostic scripts used for Sections S4.1–S4.3 (currently run as standalone snippets during manuscript development) into a single, clean, documented reproducibility script or notebook. This should be completed before the repository is made public, per the journal's requirement that "example model code and datasets" be provided "so that the key modeling/simulation steps can be replicated."

## S6. Related but not exhaustively re-verified for this submission

Full DAISY-based structural identifiability of the 5-parameter, time-varying-volume model; a global (rather than local, multi-start) practical identifiability analysis; and extension of the equivalence-family argument to three-compartment models (e.g., DigiLoCS) are all natural follow-on work, not undertaken here.
