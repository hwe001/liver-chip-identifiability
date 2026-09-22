"""
Mechanistic 2-compartment model of a recirculating liver-on-a-chip.

State variables (amounts, not concentrations):
    Am(t) : drug amount in the media compartment
    Ac(t) : drug amount in the cell-associated compartment

Media volume shrinks linearly due to evaporation:
    Vm(t) = max(Vm0 - e*t, Vm_floor)

Concentrations (what an experimentalist actually measures by sampling
the well at time t) are:
    Cm(t) = Am(t) / Vm(t)
    Cc(t) = Ac(t) / Vc          (Vc assumed constant, small cell-compartment volume)

ODEs:
    dAm/dt = -CLd*Cm + (CLd/Kp)*Cc - kb*Am
    dAc/dt =  CLd*Cm - (CLd/Kp)*Cc - CLint*Cc

Parameters:
    CLint : intrinsic metabolic clearance in the cell compartment (uL/min)
    CLd   : distributional (permeability x surface area) clearance between
            media and cell compartment (uL/min)
    Kp    : cell:media partition coefficient at equilibrium (dimensionless)
    kb    : first-order non-specific binding / adsorption loss rate from
            media (1/min)
    e     : evaporation rate (uL/min), Vm(t) = Vm0 - e*t
"""

import numpy as np
from scipy.integrate import solve_ivp
from dataclasses import dataclass


@dataclass
class Platform:
    """Fixed hardware parameters for a chip platform (not fitted)."""
    name: str
    Vm0_uL: float      # initial media volume, microliters
    Vc_uL: float       # nominal cell-compartment volume, microliters
    cells: float       # number of cells (or cell-equivalent) on chip
    e_uL_per_min: float  # evaporation rate, microliters/min (literature-anchored or assumed)
    e_source: str      # provenance note for the evaporation rate


@dataclass
class DrugParams:
    """Fitted / ground-truth drug- and system-specific parameters."""
    CLint_uL_per_min: float   # total (not per-cell) intrinsic clearance on-chip
    CLd_uL_per_min: float
    Kp: float
    kb_per_min: float


def rhs(t, y, Vm0, e, Vc, p: DrugParams):
    Am, Ac = y
    Vm = max(Vm0 - e * t, 0.05 * Vm0)  # floor to avoid division blow-up
    Cm = Am / Vm
    Cc = Ac / Vc
    dAm = -p.CLd_uL_per_min * Cm + (p.CLd_uL_per_min / p.Kp) * Cc - p.kb_per_min * Am
    dAc = p.CLd_uL_per_min * Cm - (p.CLd_uL_per_min / p.Kp) * Cc - p.CLint_uL_per_min * Cc
    return [dAm, dAc]


def simulate_media_concentration(t_eval, dose_uM, platform: Platform, drug: DrugParams):
    """Return the media concentration Cm(t) (what gets sampled) at t_eval."""
    y0 = [dose_uM * platform.Vm0_uL, 0.0]  # Am0 = dose * Vm0 (amount in concentration*volume units)
    sol = solve_ivp(
        rhs, (0, t_eval[-1]), y0,
        t_eval=t_eval,
        args=(platform.Vm0_uL, platform.e_uL_per_min, platform.Vc_uL, drug),
        method="LSODA", rtol=1e-8, atol=1e-10,
    )
    Vm_t = np.maximum(platform.Vm0_uL - platform.e_uL_per_min * t_eval, 0.05 * platform.Vm0_uL)
    Cm = sol.y[0] / Vm_t
    return Cm


# ---- Platform anchors (see references in accompanying notes) ----

CNBIO = Platform(
    name="CnBio/PhysioMimix",
    Vm0_uL=1600.0,       # 1.6 mL, DigiLoCS Table 1 (reading Docci et al. 2022)
    Vc_uL=3.0,           # nominal small cell-compartment volume (assumption, not literature)
    cells=3e5,           # DigiLoCS Table 1
    e_uL_per_min=0.07,   # Docci et al. 2022 evaporation sensitivity grid: kev in {0, 0.05, 0.07, 0.1} uL/min
    e_source="Docci et al. 2022, Fig. evaporation sensitivity analysis (kev grid, mid value chosen)",
)

JAVELIN = Platform(
    name="Javelin LTC",
    Vm0_uL=1700.0,       # 1.7 mL, Rajan et al. 2023 protocol
    Vc_uL=3.0,
    cells=2.25e5,        # 200,000-250,000 PHH, midpoint, Rajan et al. 2023
    e_uL_per_min=1700.0 * 0.004 / (24 * 60),  # <0.4%/day -> uL/min, Rajan et al. 2023
    e_source="Rajan et al. 2023: <0.4%/day evaporative loss, converted to uL/min",
)

SPHEROID = Platform(
    name="3D spheroid (Hurel-type)",
    Vm0_uL=75.0,         # 0.05-0.1 mL range midpoint, DigiLoCS Table 1 (Kanebratt 2021, Bonn 2016)
    Vc_uL=0.5,
    cells=3e4,           # DigiLoCS Table 1 (order-of-magnitude midpoint)
    e_uL_per_min=0.07 * (75.0 / 1600.0) ** (2 / 3),  # ASSUMPTION: scaled from CnBio by surface/volume^(2/3)
    e_source="ASSUMED (not from literature): scaled from CnBio evaporation by (Vm0 ratio)^(2/3) "
             "to reflect higher surface-to-volume ratio in small static wells. Flag as a sensitivity "
             "case, not a validated number.",
)

PLATFORMS = {"CnBio": CNBIO, "Javelin": JAVELIN, "Spheroid": SPHEROID}
