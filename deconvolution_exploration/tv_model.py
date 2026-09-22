"""
Exploring whether CLint is actually constant during a liver-chip incubation.

Core idea: Paper 1 derived (fast-distribution limit, CL_d >> CLint):
    CL_naive(t) = Vm*(Kp*CLint + kb*Vm) / (Vm + Kp*Vc)

If CLint is CONSTANT, log(Cm(t)) is a straight line and CL_naive extracted
from its slope is one fixed number. If CLint(t) VARIES over the incubation,
the same relation holds locally at each instant (quasi-steady-state
approximation: Cc(t) ~= Kp*Cm(t) tracks Cm(t) faster than CLint(t) changes),
which lets us invert the *local* slope of ln(Cm(t)) to recover an estimate
of CLint(t) directly -- no assumption of constancy required.

Rearranging the mass-balance under quasi-steady-state:
    (Vm + Kp*Vc) * dCm/dt = -(kb*Vm + CLint(t)*Kp) * Cm(t)
    =>  CLint(t) = [ -dln(Cm)/dt * (Vm + Kp*Vc) - kb*Vm ] / Kp

This requires Vm, Vc, Kp, kb already known (e.g. from the NSB-control +
structural-identifiability approach in Paper 1) -- this piece assumes that
problem is already solved and asks a different question: is CLint itself
constant?

Validity: quasi-steady-state requires CL_d >> CLint(t) at all times (the
same "fast-distribution" regime as Paper 1's low-mu asymptote). This is a
real, stated limitation -- for CLint(t) comparable to or exceeding CL_d,
recovering a time-varying rate needs a full LTV (linear time-varying)
state-estimation approach (e.g. Kalman filtering treating CLint(t) as a
random-walk state), which is future work, not implemented here.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import savgol_filter
from model import Platform, CNBIO


def simulate_tv_media_concentration(t_eval, dose_uM, platform: Platform,
                                     CLint_func, CLd, Kp, kb):
    """Simulate the two-compartment model with a TIME-VARYING CLint(t),
    given by CLint_func(t) (total intrinsic clearance, uL/min, not
    per-cell). Otherwise identical to model.py's constant-CLint version."""
    Vm0, Vc = platform.Vm0_uL, platform.Vc_uL

    def rhs(t, y):
        Am, Ac = y
        Cm = Am / Vm0  # constant volume for this exploration (no evaporation)
        Cc = Ac / Vc
        CLint_t = CLint_func(t)
        dAm = -CLd * Cm + (CLd / Kp) * Cc - kb * Am
        dAc = CLd * Cm - (CLd / Kp) * Cc - CLint_t * Cc
        return [dAm, dAc]

    y0 = [dose_uM * Vm0, 0.0]
    sol = solve_ivp(rhs, (0, t_eval[-1]), y0, t_eval=t_eval,
                     method="LSODA", rtol=1e-9, atol=1e-12)
    Cm = sol.y[0] / Vm0
    return Cm


def recover_instantaneous_CLint(t, Cm_obs, Vm, Vc, Kp, kb, smooth_window=9, smooth_poly=2):
    """Quasi-steady-state recovery of CLint(t) from noisy Cm(t) data.

    Smooths ln(Cm) with a Savitzky-Golay filter (this IS the
    regularization step -- differentiating raw noisy data directly is
    the classic ill-posed-deconvolution failure mode) before taking its
    derivative, then inverts the quasi-steady-state relation."""
    logCm = np.log(np.clip(Cm_obs, 1e-9, None))
    # Savitzky-Golay requires odd window <= len(t); clip if necessary
    w = min(smooth_window, len(t) - (1 - len(t) % 2))
    if w % 2 == 0:
        w -= 1
    w = max(w, smooth_poly + 2 + (1 - (smooth_poly + 2) % 2))  # ensure valid odd window
    logCm_smooth = savgol_filter(logCm, window_length=w, polyorder=smooth_poly)
    dlogCm_dt = savgol_filter(logCm, window_length=w, polyorder=smooth_poly, deriv=1, delta=(t[1] - t[0]))

    CLint_t = (-dlogCm_dt * (Vm + Kp * Vc) - kb * Vm) / Kp
    return CLint_t, logCm_smooth
