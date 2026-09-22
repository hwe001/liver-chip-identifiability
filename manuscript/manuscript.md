# Working title
**What can media-only liver-on-a-chip concentration data identify? More sampling cannot resolve structural ambiguity**

*Perspective — CPT: Pharmacometrics & Systems Pharmacology*
*Format (confirmed): 1,600-word main text + 75-word Introduction (Introduction is additional to, not part of, the 1,600 words) | ≤10 references | combined Figure/Table limit of 2*
*(No separate Abstract for this article type. Article type confirmed as Perspective; APC-free status confirmed.)*

---

### Introduction

A liver-chip well looks simple from the outside: drug goes in, fluid is sampled over time, concentration falls. But that fall reflects several processes at once — distribution into cells, metabolism, and non-specific loss — tangled together in a single curve. Previous studies have examined identifiability for particular chip models and designs. Here, with a minimal two-compartment model, we ask what the curve alone can tell us, independent of how densely it is sampled.

---

### Structural identifiability of media-only sampling

This is not a new concern. Conventional human liver microsome and hepatocyte assays have historically underpredicted in vivo metabolic intrinsic clearance, with Hallifax et al.¹ reporting average biases of approximately five- and four-fold for microsomes and hepatocytes, respectively, in their retrospective dataset. Liver-chip platforms such as the CN Bio PhysioMimix²,³ and the Javelin liver tissue chip⁴ were developed partly to improve on this, and several groups have already applied formal identifiability analysis to specific configurations. Docci et al.² performed a priori structural identifiability, via DAISY, for three candidate models of their system. Milani et al.⁵ assessed structural and practical identifiability alongside global sensitivity analysis for a mycophenolate mofetil Gut-Liver OoC. Milani et al.⁶ extended this into a general experimental-design tutorial demonstrated on midazolam, noting explicitly that structural identifiability depends on the model and the compartments sampled, and does not by itself guarantee precise estimation.

We ask the same question at its simplest: for a minimal two-compartment model — media C_m (volume V_m) exchanging with cells C_c (volume V_c) via distributional clearance CL_d and partition coefficient K_p, first-order non-specific binding (NSB, rate k_b), and metabolism CLint confined to the cells — what can watching only the media compartment ever tell us? This reduced model captures the processes common to more detailed liver-chip models, such as the three-compartment DigiLoCS digital twin⁷, without their added complexity. Observing y(t)=C_m(t) with dose in media only, the Laplace transform gives

Y(s)/C_0 = (s−a22) / (s²−(a11+a22)s+(a11a22−a12a21))

with a11=−(CL_d/V_m+k_b), a22=−(CL_d/(K_pV_c)+CLint/V_c), and a12a21 = CL_d²/(K_pV_mV_c) ≡ P. However long or densely you sample, only a11, a22, and P can ever be recovered from this curve — three identifiable combinations for four mechanistic parameters. The shortfall is not abstract. Choosing any admissible k_b defines a corresponding CL_d, K_p, and CLint through

CL_d=V_m(−a11−k_b),  K_p=CL_d²/(PV_mV_c),  CLint=−a22V_c−CL_d/K_p,

so an entire curve of parameter combinations — not just two or three isolated alternatives — reproduces exactly the same media trajectory.

This degeneracy can be resolved. Independently measuring k_b — for example via a cell-free NSB control, already run by Docci et al.² and Rajan et al.⁴ — collapses that curve to a single point, making CL_d, K_p, and CLint identifiable through the same three relations. It is not the only way: sampling the cells directly, or independently fixing K_p or CL_d, would work equally well (cf. Milani et al.⁶). This fix has a limit worth stating: it resolves the first-order NSB representation assumed here, but real adsorption can be reversible and depends on the compound (Rajan et al.⁴ found >85% recovery for 17 of 21 compounds, but appreciably lower recovery for some lipophilic acidic drugs), so a cell-free control characterizes k_b for that representation, not necessarily every loss mechanism at play. Evaporation is more tractable: because it shows up as a measurable V_m(t), it can in principle be supplied to the model as a known input rather than estimated jointly, provided the volume trajectory and any sampling withdrawals are properly tracked.

### Exact bias limits and a dimensionless reduction

Define the conventional apparent clearance as CL_naive=V_mk_dep, where k_dep is the monoexponential media-depletion rate constant an analyst would ordinarily fit. Two extremes make clear how this conventional analysis can mislead. When distribution is much faster than metabolism, drug in the media and the cells reaches equilibrium almost immediately, and what the media curve actually reflects is metabolism and NSB acting together across the whole accessible volume — not metabolism alone. When metabolism is instead much faster than distribution, the cells behave like an absorbing sink: no matter how much faster CLint becomes, media depletion cannot outrun the rate at which drug arrives at the cells in the first place, so the apparent clearance simply saturates. Formally, in the fast-distribution limit (CL_d/CLint→∞), CL_naive → V_m(K_pCLint+k_bV_m)/(V_m+K_pV_c); in the fast-metabolism limit (CLint/CL_d→∞), CL_naive saturates exactly at CL_d+k_bV_m.

These two regimes, and everything between them, collapse neatly if we stop tracking six dimensional quantities and instead track three ratios. Writing τ=(CL_d/V_m)t, x=C_m/C_m(0), y=C_c/(K_pC_m(0)) turns the system into dx/dτ=−(1+β)x+y, dy/dτ=(1/r)[x−(1+μ)y], where r=K_pV_c/V_m, β=k_bV_m/CL_d, and μ=K_pCLint/CL_d. The concentration dynamics x(τ) depend on nothing but these three numbers: two devices with entirely different volumes, cell counts, and absolute clearances will trace out the identical curve if r, β, and μ agree (Figure 1). The bias in CL_naive/CLint is not quite as clean, since it compares apparent clearance directly against CLint, pulling K_p back in on its own: the fast-distribution limit becomes CL_naive/CLint → [K_p/(1+r)](1+β/μ), simplifying to K_p/(1+r) once NSB is small relative to metabolism (β≪μ), while the fast-metabolism limit is CL_naive/CLint → (K_p/μ)(1+β). One consequence is asymmetric: as CLint/CL_d grows, the apparent-to-true ratio falls toward zero no matter what — an underestimate that only worsens. In the fast-distribution regime, by contrast, whether the conventional analysis over- or underestimates CLint depends on K_p, r, and how much NSB contributes.

### A methodological caution on practical identifiability

Near the structural degeneracy, even the tools used to assess parameter uncertainty became unreliable. Fisher-information covariance estimates had condition numbers of 10⁹–10¹⁷, while single-start profile likelihoods were vulnerable to optimizer failure. Multi-start profiling was more stable, but still indicated weak constraint on CLint under the sampling designs examined, even with k_b and evaporation fixed at their true values. We offer this as a caution rather than a number: near a structural degeneracy, standard local diagnostics can mislead, so practical estimability needs to be checked directly, not assumed once structural identifiability is established — consistent with the staged approach Milani et al.⁶ advocate. Full diagnostic detail is in the Supplementary Material.

### Toward a minimum reporting standard

These results have a direct consequence for how liver-chip PK studies should report their methods: a reader cannot judge whether a reported CLint is identifiable and adequately constrained unless the study also reports the information needed to assess it. Existing MPS quality-management recommendations⁸ are focused primarily on biological reproducibility rather than on this kind of identifiability assessment. Table 1 sets out what we think a liver-chip PK study should report if its CLint estimate is meant to support IVIVE.

### Limitations

K_p, CL_d, and the dimensionless-group values used in our numerical checks are illustrative rather than tied to any specific compound. The closed-form structural result holds exactly for the constant-volume reduction; we did not run a DAISY-based structural analysis of the full five-parameter, time-varying-volume model. Nor do we claim to reproduce the bias of any particular published dataset, or that this reduced model stands in for every liver-chip analysis in use today.

### Conclusion

A media concentration curve, on its own, identifies three combinations of the four mechanistic parameters in this model, leaving a continuous family of equally valid parameterizations. Measuring NSB independently is one way to collapse that family and make CL_d, K_p, and CLint identifiable. The same simple model also gives exact limits on how badly conventional analysis can be biased, and shows that normalized concentration dynamics reduce to just three ratios, regardless of a platform's absolute dimensions. But identifiability in principle is only half the story: CLint can still be poorly constrained by real, finite, noisy data even once the ambiguity above is resolved. We think liver-chip studies aiming to support IVIVE should report what is needed to check both — including compound-specific NSB data and the raw concentration-time curves themselves.

---

## Methods (brief)

All derivations and numerical/symbolic verification were performed in Python (NumPy, SciPy, SymPy); code and full diagnostic detail (Supplementary Material) provided as supplementary material.

## Table 1. Minimum reporting checklist for liver-chip PK studies

| # | Report | Item |
|---|--------|------|
| 1 | ✓ | Compartment volumes and sampling/replacement volumes |
| 2 | ✓ | Cell number or measured functional cell mass |
| 3 | ✓ | Flow/recirculation rate and geometry |
| 4 | ✓ | Evaporation trajectory |
| 5 | ✓ | Compound-specific cell-free recovery/NSB data |
| 6 | ✓ | Raw concentration-time data with below-LOQ values identified |
| 7 | ✓ | Actual sampling times and analytical LOQ |
| 8 | ✓ | Model equations, parameter definitions, units, externally fixed parameters |
| 9 | recommend | Structural identifiability assessment for the specified observation scheme |
| 10 | recommend | Practical-identifiability/uncertainty analysis appropriate near any known structural degeneracy, not point estimates alone |

## Figure 1

Built (see figure1_schematic.pdf / figure1_schematic.png, 178 mm width, 600 dpi). Three-panel schematic: (A) information bottleneck, C_m(t) → {a11,a22,P} vs. {CL_d,K_p,CLint,k_b}; (B) continuous equivalence family parametrized by k_b, collapsing to a point once k_b is independently known; (C) structural identifiability ⇒ (necessary, not sufficient) ⇒ practical estimability. Full legend below, after References, per journal convention.

## References (numbered, in order of first appearance — superscript in-text citation style per journal convention)

1. Hallifax D, Foster JA, Houston JB. Prediction of human metabolic clearance from in vitro systems: retrospective analysis and prospective view. Pharm Res. 2010;27(10):2150-2161.
2. Docci L, Milani N, Ramp T, Romeo AA, Godoy P, Franyuti DO, et al. Exploration and application of a liver-on-a-chip device in combination with modelling and simulation for quantitative drug metabolism studies. Lab Chip. 2022;22(6):1187-1205.
3. Tsamandouras N, Kostrzewski T, Stokes CL, Griffith LG, Hughes DJ, Cirit M. Quantitative assessment of population variability in hepatic drug metabolism using a perfused three-dimensional human liver microphysiological system. J Pharmacol Exp Ther. 2017;360(1):95-105.
4. Rajan SAP, Sherfey J, Ohri S, Nichols L, Smith JT, Parekh P, et al. A novel milli-fluidic liver tissue chip with continuous recirculation for predictive pharmacokinetics applications. AAPS J. 2023;25(6):102. Correction: AAPS J. 2024;26(1):2.
5. Milani N, Parrott N, Ortiz Franyuti D, et al. Application of a gut-liver-on-a-chip device and mechanistic modelling to the quantitative in vitro pharmacokinetic study of mycophenolate mofetil. Lab Chip. 2022;22(15):2853-2868.
6. Milani N, Parrott N, Galetin A, Fowler S, Gertz M. In silico modeling and simulation of organ-on-a-chip systems to support data analysis and a priori experimental design. CPT Pharmacometrics Syst Pharmacol. 2024;13:524-543.
7. Aravindakshan MR, Mandal C, Pothen A, Schaller S, Maass C. DigiLoCS: A leap forward in predictive organ-on-chip simulations. PLoS One. 2025;20(1):e0314083.
8. Pamies D, Ekert J, Zurich MG, et al. Recommendations on fit-for-purpose criteria to establish quality management for microphysiological systems and for monitoring their reproducibility. Stem Cell Reports. 2024;19(5):604-617.

**NOT YET VERIFIED — full author lists for refs 5 and 8 are incomplete (shown as "et al." beyond the lead author) because I have not independently confirmed every co-author name against the primary source. Complete these from the actual papers before submission; do not submit with a truncated author list where the journal requires all authors.**

---

## Figure Legend

**Figure 1. Media-only concentration data underdetermine the mechanistic model; resolving this does not by itself guarantee practical estimability.**
**(A)** The observed media concentration-time trajectory Cm(t) determines only three combinations of the model's rate constants (a11, a22, and P = a12a21), while the mechanistic model has four unknown parameters (CLd, Kp, CLint, kb) — an information bottleneck of exactly one degree of freedom.
**(B)** Schematic illustration (not to scale) of the resulting one-parameter continuous equivalence family: for any admissible choice of kb, corresponding values of CLd, Kp, and CLint exist that reproduce an identical Cm(t). Independent measurement of kb (e.g., a cell-free non-specific-binding control) selects a single point on this curve, corresponding to the true parameter values.
**(C)** Structural identifiability (a unique noise-free solution, achieved once kb is independently fixed) is necessary but not sufficient for practical estimability (reliable parameter recovery from finite, noisy data), which must be assessed separately and explicitly.

---

## Change log vs. previous draft (for tracking; delete before submission)

**v10 → this version:**
1. Title changed to "What can media-only liver-on-a-chip concentration data identify? More sampling cannot resolve structural ambiguity" — sharper and directly quotes the paper's most memorable finding, at the cost of no longer signaling the bias-asymptotics and reporting-checklist contributions in the title itself (both remain fully present in the body).
2. Figure 1 regenerated with panel text stripped to the essentials (box labels only: variable sets, "3 vs. 4", curve + star + "k_b known", and the three-box necessary/sufficient chain) — all explanatory prose that previously lived inside the panels now lives only in the Figure Legend, where it belongs.

**v10 changes (hybrid pass: ~70% v9 narrative + 30% v8 restraint, per reviewer's explicit sentence-by-sentence guidance):**
1. Cut the "Picture what actually happens inside the well..." paragraph entirely — it restated what the Introduction already established, and removing it let the section open faster with "This is not a new concern."
2. Fixed "three equations, four unknowns" → "three identifiable combinations for four mechanistic parameters" (the original undersold what was actually proved).
3. Softened three accumulating conversational phrases: "There is a natural way out" → "This degeneracy can be resolved"; "And the fix has a limit worth stating plainly" → "This fix has a limit worth stating"; "Evaporation is more forgiving" → "Evaporation is more tractable" (reverting to the more precise v7/v8 wording).
4. "Two extremes make clear how this everyday analysis can mislead" → "...this conventional analysis..." (removed "everyday").
5. Practical-identifiability section: restored v8's direct, observation-first opening ("Near the structural degeneracy, even the tools used to assess parameter uncertainty became unreliable...") in place of v9's "we tried X, X failed" narrative framing, per reviewer's specific critique that this section's story should be about the phenomenon, not about what we tried.
6. Reporting-standard section: "trustworthy" → "identifiable and adequately constrained" (trustworthy oversells beyond what identifiability analysis alone can support); softened the categorical claim about MPS recommendations ("speak to...but not to this" → "focused primarily on...rather than on this kind of identifiability assessment").
7. Conclusion: removed the rhetorically absolute "never the parameters themselves," replaced with "the four mechanistic parameters in this model, leaving a continuous family of equally valid parameterizations" — same meaning, but properly scoped to the stated model rather than reading as a universal claim.
8. Word count: 1,210/1,600 main text (down from v9's 1,329, mostly from cutting the redundant paragraph), Introduction unchanged at 75/75.

**Kept from v9 (per reviewer's explicit "keep" list):** the vivid Introduction; "We ask the same question at its simplest... what can watching only the media compartment ever tell us?"; "The shortfall is not abstract"; "However long or densely you sample..."; the intuitive physical explanation of both asymptotic limits before the formulas; the motivating sentence for nondimensionalization; the "two devices... trace out the identical curve" interpretation; the accessible Conclusion structure.

## Still to verify before submission

- Complete author lists for references 5 (Milani et al. 2022) and 8 (Pamies et al. 2024).
- Consolidate diagnostic scripts into one clean reproducibility notebook.
- This version has now been reviewed for science (multiple rounds), format/citations, and prose/voice (two rounds, converging on this hybrid) — a final cold read is still worthwhile before submission, but no further major revision axis is expected.
