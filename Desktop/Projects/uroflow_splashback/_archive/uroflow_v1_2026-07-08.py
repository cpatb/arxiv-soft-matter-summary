"""
uroflow.py -- Splashback physics of a urine stream: Rayleigh-Weber jet breakup
plus projectile stream-tracking, framed as a standing-distance / aim-angle
optimisation (a sibling of the bidet directed-jet problem).

Core chain of physics
---------------------
1. Exit conditions. Peak uroflow  Q ~ 15-25 mL/s through an effective meatus
   diameter d ~ 3-6 mm gives an exit speed  v = Q / (pi d^2 / 4) ~ 1-4 m/s.

2. Plateau-Rayleigh breakup. The jet is unstable to varicose perturbations; the
   fastest-growing mode sets a breakup LENGTH L_b beyond which the coherent stream
   has fragmented into drops. We compute the max growth rate from Weber's viscous
   long-wave dispersion relation (reduces to Rayleigh's inviscid result) and use
   L_b = v/omega_max * ln(r/delta0), with ln(r/delta0) ~ 12 (Weber 1931).

3. Splash. The Splash Lab (Hurd & Truscott, APS-DFD 2013) showed a COHERENT
   stream splashes far less than a train of drops. So the controlling question is
   whether the stream reaches the target surface BEFORE it breaks up: compare the
   trajectory path length s to L_b. For post-breakup impact, the drop Weber number
   and impact angle set the splash severity (grazing impact and pre-breakup impact
   both suppress splashback).

4. Trajectory. The stream centreline is itself a projectile launched from height h
   at angle theta and speed v; we integrate it (gravity; drag negligible for a
   coherent cm-scale stream over <0.5 m) to the target surface to get path length
   and impact angle.

Deliverable: a splashback map over (standing distance, aim angle) and concrete,
quantitative recommendations. SI units throughout.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize_scalar

# --- Fluid properties: urine is close to water, slightly higher solute load. ---
RHO = 1020.0        # kg/m^3
MU = 1.0e-3         # Pa.s   (~1.0-1.1 mPa.s, near water)
SIGMA = 0.066       # N/m    urine-air surface tension (~66 mN/m)
G = 9.81
LN_R_DELTA0 = 12.0  # Weber's empirical noise-amplification factor for L_b


# --------------------------------------------------------------------------- #
# 1-2. Jet exit + Rayleigh-Weber breakup.
# --------------------------------------------------------------------------- #
def exit_velocity(Q, d):
    """Exit speed (m/s) from volumetric flow Q (m^3/s) and diameter d (m)."""
    return Q / (np.pi * d**2 / 4.0)


def ohnesorge(d):
    """Oh = mu / sqrt(rho sigma d)  (diameter-based)."""
    return MU / np.sqrt(RHO * SIGMA * d)


def max_growth_rate(d):
    """Fastest-growing varicose-mode growth rate omega_max (1/s) from Weber's
    viscous long-wave dispersion relation on a jet of radius r = d/2:

        omega(x) = -1.5 nu k^2 + sqrt( (1.5 nu k^2)^2 + (sigma/(2 rho r^3)) x^2 (1-x^2) )

    with x = k r in (0,1), k = x/r, nu = mu/rho. Returns (omega_max, x_opt).
    Inviscid limit -> omega_max = 0.354 sqrt(sigma/(rho r^3)) at x = 1/sqrt(2),
    recovering Rayleigh's classical result.
    """
    r = d / 2.0
    nu = MU / RHO
    base = SIGMA / (2.0 * RHO * r**3)

    def neg_omega(x):
        k = x / r
        visc = 1.5 * nu * k**2
        disc = visc**2 + base * x**2 * (1.0 - x**2)
        if disc < 0:
            return 0.0
        return -(-visc + np.sqrt(disc))

    res = minimize_scalar(neg_omega, bounds=(1e-4, 0.999), method="bounded")
    return -res.fun, res.x


def breakup_length(Q, d):
    """Intact (coherent) length L_b (m) of the jet before it fragments."""
    v = exit_velocity(Q, d)
    omega, _ = max_growth_rate(d)
    return v / omega * LN_R_DELTA0


def weber(v, d):
    return RHO * v**2 * d / SIGMA


def reynolds(v, d):
    return RHO * v * d / MU


# --------------------------------------------------------------------------- #
# 4. Stream trajectory (projectile) to a target surface.
# --------------------------------------------------------------------------- #
def trajectory(v, theta, h, horizontal_target, n=400):
    """Integrate the stream centreline from launch (height h, angle theta below
    horizontal is negative; aiming downward uses theta<0) until it reaches a
    vertical target plane at x = horizontal_target OR hits the floor z=0.

    Drag on a coherent ~cm stream over sub-metre distances is negligible (the
    stream momentum flux dwarfs air drag), so pure ballistic.

    Returns dict with path length s, impact point, impact angle to the surface.
    """
    # Analytic parabola sampled densely; robust and drag-free.
    vx = v * np.cos(theta)
    vz = v * np.sin(theta)          # theta<0 => downward launch
    if vx <= 1e-6:
        vx = 1e-6
    t_target = horizontal_target / vx
    ts = np.linspace(0, t_target, n)
    x = vx * ts
    z = h + vz * ts - 0.5 * G * ts**2
    # Truncate if it hits the floor before the target.
    below = np.where(z <= 0)[0]
    if below.size:
        i = below[0]
        ts, x, z = ts[:i + 1], x[:i + 1], z[:i + 1]
    # Path length via arc length.
    s = np.sum(np.hypot(np.diff(x), np.diff(z)))
    # Impact velocity direction at the end.
    t_end = ts[-1]
    vz_end = vz - G * t_end
    speed_end = np.hypot(vx, vz_end)
    # Impact angle to a VERTICAL wall surface (angle between velocity and wall):
    # wall normal is horizontal (x); grazing (small) angle = better.
    impact_angle_wall = np.arctan2(abs(vz_end), abs(vx))  # from horizontal
    return dict(x=x, z=z, s=s, x_impact=x[-1], z_impact=z[-1],
                speed_impact=speed_end, angle_from_horizontal=impact_angle_wall,
                hit_floor=bool(below.size))


# --------------------------------------------------------------------------- #
# 3. Splash severity metric.
# --------------------------------------------------------------------------- #
def splash_index(Q, d, theta, h, horizontal_target):
    """Dimensionless splashback severity in [0, ~], combining:
      - coherence penalty: 1 if impact occurs AFTER breakup (drops), <1 if before,
      - Weber factor: sqrt(We_impact) (splash grows with impact momentum),
      - angle factor: sin(impact angle) (normal impact splashes most; grazing least).

    This is a physically-motivated composite, not a measured calibration; it ranks
    configurations consistently for optimisation. Lower = drier.
    """
    v = exit_velocity(Q, d)
    traj = trajectory(v, theta, h, horizontal_target)
    L_b = breakup_length(Q, d)
    # Coherence: smooth transition around s = L_b.
    frac = traj["s"] / L_b
    coherence_penalty = 1.0 / (1.0 + np.exp(-(frac - 1.0) / 0.15))  # ~0 pre-breakup, ~1 post
    coherence_penalty = 0.15 + 0.85 * coherence_penalty            # coherent still splashes a little
    We_imp = weber(traj["speed_impact"], d)
    angle_factor = np.sin(traj["angle_from_horizontal"])           # normal (vertical wall) worst
    return dict(index=coherence_penalty * np.sqrt(We_imp) * angle_factor,
                L_b=L_b, path=traj["s"], frac_of_Lb=frac,
                We_impact=We_imp, angle_deg=np.degrees(traj["angle_from_horizontal"]),
                coherent=frac < 1.0, hit_floor=traj["hit_floor"], traj=traj)
