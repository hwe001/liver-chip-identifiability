"""
Structural identifiability of the 2-compartment liver-chip model
(constant media volume, i.e. evaporation e treated separately -- see
note at the end of this file) when only the MEDIA compartment is
sampled (the standard practice: cell-associated/intracellular
concentrations are rarely sampled serially, since that requires
lysing the tissue).

State-space form (Cm, Cc), constant Vm:
    dCm/dt = a11*Cm + a12*Cc
    dCc/dt = a21*Cm + a22*Cc
where
    a11 = -(CLd/Vm + kb)
    a12 =  CLd/(Kp*Vm)
    a21 =  CLd/Vc
    a22 = -(CLd/(Kp*Vc) + CLint/Vc)

Observing y(t) = Cm(t) from an impulse-type initial condition
(Cm(0)=C0, Cc(0)=0), the Laplace transform is

    Y(s)/C0 = (s - a22) / (s^2 - (a11+a22) s + (a11*a22 - a12*a21))

Only three quantities are recoverable from noise-free data of arbitrary
duration and density: a22, (a11+a22) [hence a11], and the product
a12*a21. These are the *identifiable combinations*:

    A11 = a11                              = -(CLd/Vm + kb)
    A22 = a22                              = -(CLd/(Kp*Vc) + CLint/Vc)
    P   = a12*a21 = CLd^2 / (Kp*Vm*Vc)

That is 3 equations in 4 unknowns (CLint, CLd, Kp, kb) once Vm, Vc are
fixed (known) hardware constants -- the model is STRUCTURALLY
UNIDENTIFIABLE from media-only sampling: a one-parameter family of
(CLint, CLd, Kp, kb) combinations produces exactly the same observed
trajectory (verified numerically in the accompanying script to
solver precision, ~1e-7).

If kb is independently known (e.g. from a cell-free non-specific-
binding control -- standard practice in Docci et al. 2022 and Rajan
et al. 2023), all three remaining parameters become uniquely
determined:
    CLd   = Vm * (-A11 - kb)
    Kp    = CLd^2 / (P * Vm * Vc)
    CLint = -A22 * Vc - CLd / Kp

This file provides both the symbolic derivation (reproducible) and a
solver for the two identifiability regimes.

Note on evaporation (e): e enters the model as a known, directly
measurable time-varying volume Vm(t) = Vm0 - e*t (from residual-volume
weighing), not as a kinetic parameter competing for information in the
same way as kb. Treating e as measured-and-fixed (as the field already
does) sidesteps this identifiability problem for e specifically; the
structural result above concerns the four *kinetic* parameters at
fixed, known Vm(t).
"""

import sympy as sp
import numpy as np
from model import Platform, DrugParams, simulate_media_concentration

CLint, CLd, Kp, kb, Vm, Vc, s = sp.symbols("CLint CLd Kp kb Vm Vc s", positive=True)

a11 = -(CLd / Vm + kb)
a12 = CLd / (Kp * Vm)
a21 = CLd / Vc
a22 = -(CLd / (Kp * Vc) + CLint / Vc)

A11_expr = sp.simplify(a11)
A22_expr = sp.simplify(a22)
P_expr = sp.simplify(a12 * a21)


def identifiable_combinations(params: DrugParams, Vm0_uL: float, Vc_uL: float):
    """Return the 3 structurally identifiable combinations (A11, A22, P)
    for a given (true) parameter set and fixed hardware volumes."""
    subs = {CLint: params.CLint_uL_per_min, CLd: params.CLd_uL_per_min,
            Kp: params.Kp, kb: params.kb_per_min, Vm: Vm0_uL, Vc: Vc_uL}
    return (float(A11_expr.subs(subs)), float(A22_expr.subs(subs)), float(P_expr.subs(subs)))


def is_structurally_identifiable_without_kb(A11_val=None, A22_val=None, P_val=None, Vm0_uL=None, Vc_uL=None):
    """3 equations, 4 unknowns (CLint, CLd, Kp, kb) at fixed Vm, Vc:
    always underdetermined by exactly 1 degree of freedom. Returns the
    1-parameter family (free variable: kb) as explicit expressions,
    numeric if A11_val/A22_val/P_val/Vm0_uL/Vc_uL are given, else fully symbolic."""
    CLd_s, Kp_s, CLint_s, kb_s = sp.symbols("CLd_s Kp_s CLint_s kb_s", positive=True)
    A11_sym, A22_sym, P_sym = sp.symbols("A11 A22 P")
    A11_use = A11_val if A11_val is not None else A11_sym
    A22_use = A22_val if A22_val is not None else A22_sym
    P_use = P_val if P_val is not None else P_sym
    Vm_use = Vm0_uL if Vm0_uL is not None else Vm
    Vc_use = Vc_uL if Vc_uL is not None else Vc

    eqs = [
        sp.Eq(-(CLd_s / Vm_use + kb_s), A11_use),
        sp.Eq(-(CLd_s / (Kp_s * Vc_use) + CLint_s / Vc_use), A22_use),
        sp.Eq(CLd_s ** 2 / (Kp_s * Vm_use * Vc_use), P_use),
    ]
    family = sp.solve(eqs, [CLd_s, Kp_s, CLint_s], dict=True)
    return family, kb_s  # family[i][var] is an expression in kb_s (and A11/A22/P if symbolic)


def resolve_with_known_kb(A11_val, A22_val, P_val, Vm0_uL, Vc_uL, kb_known):
    """Closed-form unique solution for CLd, Kp, CLint once kb is fixed
    (e.g. measured via a cell-free non-specific-binding control)."""
    CLd_hat = Vm0_uL * (-A11_val - kb_known)
    Kp_hat = CLd_hat ** 2 / (P_val * Vm0_uL * Vc_uL)
    CLint_hat = -A22_val * Vc_uL - CLd_hat / Kp_hat
    return CLd_hat, Kp_hat, CLint_hat


def demonstrate_degeneracy(platform: Platform, drug_true: DrugParams, kb_alt: float, t_eval):
    """Build the alternate (CLd, Kp, CLint) that is structurally
    indistinguishable from drug_true when kb is unknown, then simulate
    both and return their trajectories for direct comparison."""
    A11_val, A22_val, P_val = identifiable_combinations(drug_true, platform.Vm0_uL, platform.Vc_uL)
    CLd_alt, Kp_alt, CLint_alt = resolve_with_known_kb(A11_val, A22_val, P_val,
                                                        platform.Vm0_uL, platform.Vc_uL, kb_alt)
    drug_alt = DrugParams(CLint_uL_per_min=CLint_alt, CLd_uL_per_min=CLd_alt, Kp=Kp_alt, kb_per_min=kb_alt)

    flat_platform = Platform(platform.name, platform.Vm0_uL, platform.Vc_uL, platform.cells,
                              0.0, "constant-volume (e=0) for structural-identifiability check")
    Cm_true = simulate_media_concentration(t_eval, 1.0, flat_platform, drug_true)
    Cm_alt = simulate_media_concentration(t_eval, 1.0, flat_platform, drug_alt)
    return drug_alt, Cm_true, Cm_alt


if __name__ == "__main__":
    print("Symbolic identifiable combinations (media-only observation):")
    print("  A11 =", A11_expr)
    print("  A22 =", A22_expr)
    print("  P   =", P_expr)
    print()
    print("With kb UNKNOWN (numeric example, CnBio hardware, true params below):")
    from model import CNBIO
    drug_true0 = DrugParams(CLint_uL_per_min=4.5, CLd_uL_per_min=8.0, Kp=5.0, kb_per_min=0.0008)
    A11v, A22v, Pv = identifiable_combinations(drug_true0, CNBIO.Vm0_uL, CNBIO.Vc_uL)
    family, kb_s = is_structurally_identifiable_without_kb(A11v, A22v, Pv, CNBIO.Vm0_uL, CNBIO.Vc_uL)
    print(f"  (CLd, Kp, CLint) as a function of free parameter kb={kb_s}:")
    for sol in family:
        for k, v in sol.items():
            print(f"    {k} = {sp.simplify(v)}")
    print()

    from model import CNBIO
    drug_true = DrugParams(CLint_uL_per_min=4.5, CLd_uL_per_min=8.0, Kp=5.0, kb_per_min=0.0008)
    t_eval = np.linspace(1, 2000, 30)
    drug_alt, Cm_true, Cm_alt = demonstrate_degeneracy(CNBIO, drug_true, kb_alt=0.003, t_eval=t_eval)
    print(f"True params:      CLint={drug_true.CLint_uL_per_min}, CLd={drug_true.CLd_uL_per_min}, "
          f"Kp={drug_true.Kp}, kb={drug_true.kb_per_min}")
    print(f"Degenerate params: CLint={drug_alt.CLint_uL_per_min:.4f}, CLd={drug_alt.CLd_uL_per_min:.4f}, "
          f"Kp={drug_alt.Kp:.4f}, kb={drug_alt.kb_per_min}")
    print(f"Max relative difference in Cm(t) between the two: {np.max(np.abs((Cm_true - Cm_alt) / Cm_true)):.2e}")
