"""
bpa.py -- Bloodstain Pattern Analysis: ballistic region-of-origin reconstruction.

Physical premise
----------------
A bloodstain on a surface is the fingerprint of a single drop-impact event.
Three geometric measurements per stain are used:

  1. Ellipse aspect ratio  ->  impact angle   alpha = arcsin(W / L)      [Balthazard 1939]
  2. Ellipse major-axis azimuth (+ tail)      ->  in-surface direction of travel
  3. Ellipse size (minor axis W) + drop size  ->  impact SPEED, via a drop-impact
                                                   spread-factor correlation.

The first two are standard forensic practice. The third is the bridge to
drop-impact physics: on a smooth surface the maximum spread factor
xi = D_stain / D_0 is a monotonic function of the (normal) Weber number
We_n = rho * D_0 * v_n^2 / sigma  [Pasandideh-Fard et al., Phys. Fluids 8, 650 (1996)].
Measuring xi therefore *measures the normal impact speed*, which -- together with
the impact angle -- gives the full impact velocity VECTOR at each stain.

Given a velocity vector at a known landing point we can integrate the equations of
motion BACKWARD in time (gravity + aerodynamic drag) to recover the flight path.
The common point through which all back-integrated paths pass is the region of
origin. The classical "straight-line stringing" method is the special case that
ignores gravity; it is known to OVERESTIMATE the height of origin. This module
quantifies that bias and shows the ballistic method removes it.

Everything here is validated by round-trip: emit droplets from a KNOWN origin,
forward-integrate to synthetic stains, then reconstruct and check we recover the
origin. Units are SI throughout (metres, seconds, kilograms) unless noted.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize

# --------------------------------------------------------------------------- #
# Physical constants (whole blood at ~body/room temperature; air at 20 C).
# Blood is weakly shear-thinning; at the high shear of impact mu ~ 3-5 mPa.s,
# so a Newtonian value is used for the (drag-only, non-impact) flight phase.
# --------------------------------------------------------------------------- #
RHO_BLOOD = 1060.0      # kg/m^3   whole-blood density
SIGMA_BLOOD = 0.062     # N/m      blood-air surface tension
MU_BLOOD = 4.0e-3       # Pa.s     high-shear apparent viscosity (flight-phase drag negligible internally)
RHO_AIR = 1.20          # kg/m^3
MU_AIR = 1.81e-5        # Pa.s
G = 9.81                # m/s^2

# Advancing contact angle of blood on a smooth non-porous surface (deg).
# Enters the spread model weakly; blood on glass/painted drywall ~ 100-140 deg.
THETA_C_DEG = 120.0


# --------------------------------------------------------------------------- #
# Drop-impact spread factor  xi(We_n)  -- the tie to drop-impact physics.
# --------------------------------------------------------------------------- #
def spread_factor(we_n: np.ndarray, re_n: np.ndarray,
                  theta_c_deg: float = THETA_C_DEG) -> np.ndarray:
    """Maximum spread factor xi = D_stain / D_0 (Pasandideh-Fard et al. 1996).

        xi = sqrt( (We + 12) / ( 3(1 - cos theta_c) + 4 We / sqrt(Re) ) )

    Monotonically increasing in We_n, which is what makes stain size invertible
    for impact speed. We_n and Re_n use the surface-normal velocity component,
    the classical driver of splat spreading.
    """
    we_n = np.asarray(we_n, dtype=float)
    re_n = np.asarray(re_n, dtype=float)
    theta = np.radians(theta_c_deg)
    denom = 3.0 * (1.0 - np.cos(theta)) + 4.0 * we_n / np.sqrt(re_n)
    return np.sqrt((we_n + 12.0) / denom)


def n_spines(we_n: float) -> float:
    """Expected number of rim spines/fingers ~ 1.14 * We^{1/2} (Rayleigh rim
    instability; Marmanis & Thoroddsen 1996). A secondary size-independent
    estimator of We -- useful when the parent drop diameter D_0 is unknown."""
    return 1.14 * np.sqrt(max(we_n, 0.0))


def _we_re_from_vn(v_n: float, d0: float):
    we = RHO_BLOOD * d0 * v_n**2 / SIGMA_BLOOD
    re = RHO_BLOOD * d0 * v_n / MU_BLOOD
    return we, re


def vn_from_spread(xi: float, d0: float,
                   theta_c_deg: float = THETA_C_DEG) -> float:
    """Invert the spread model: given measured spread factor xi and parent drop
    diameter d0, solve for the normal impact speed v_n (m/s).

    Monotonicity of xi(We_n) guarantees a unique root; solved with brentq.
    """
    if xi <= 1.0:
        # No measurable spreading beyond the drop footprint -> ~zero normal speed.
        return 0.0

    def resid(v_n):
        we, re = _we_re_from_vn(v_n, d0)
        return spread_factor(we, re, theta_c_deg) - xi

    lo, hi = 1e-3, 50.0
    if resid(lo) > 0:          # xi below the model's minimum -> clamp
        return lo
    if resid(hi) < 0:          # xi beyond model at 50 m/s -> clamp (unphysical for BPA)
        return hi
    return brentq(resid, lo, hi, xtol=1e-6)


# --------------------------------------------------------------------------- #
# Ballistic flight: gravity + aerodynamic drag on a spherical drop.
# --------------------------------------------------------------------------- #
def drag_coefficient(re: np.ndarray) -> np.ndarray:
    """Sphere drag coefficient. Schiller-Naumann for Re < 1000, constant 0.44
    beyond -- standard for mm-scale drops at Re ~ 10^2-10^3."""
    re = np.asarray(re, dtype=float)
    cd = np.where(re > 1e-9, 24.0 / np.maximum(re, 1e-9) * (1.0 + 0.15 * re**0.687), 1e9)
    return np.where(re > 1000.0, 0.44, cd)


def _accel(vel: np.ndarray, d0: float, drag: bool) -> np.ndarray:
    """Acceleration of a drop of diameter d0 moving at velocity `vel` (3-vector).
    a = g - (3 rho_air Cd |v|)/(4 rho_blood d0) * v  (per-unit-mass drag)."""
    g_vec = np.array([0.0, 0.0, -G])
    if not drag:
        return g_vec
    speed = np.linalg.norm(vel)
    if speed < 1e-12:
        return g_vec
    re = RHO_AIR * speed * d0 / MU_AIR
    cd = float(drag_coefficient(re))
    a_drag = -(3.0 * RHO_AIR * cd * speed) / (4.0 * RHO_BLOOD * d0) * vel
    return g_vec + a_drag


def integrate_flight(x0, v0, d0, direction=+1, drag=True,
                     t_max=5.0, stop_plane_z=None, max_step=2e-3,
                     rtol=1e-8, atol=1e-10):
    """Integrate a drop's trajectory. direction=+1 forward in time, -1 backward.

    Returns (t, X, V) arrays. If stop_plane_z is given (forward only), integration
    stops when the drop crosses that horizontal plane from above.
    """
    x0 = np.asarray(x0, float)
    v0 = np.asarray(v0, float)

    def rhs(t, y):
        vel = y[3:]
        acc = _accel(vel, d0, drag)
        return np.concatenate([vel, acc]) * direction

    events = None
    if stop_plane_z is not None and direction == +1:
        def hit(t, y):
            return y[2] - stop_plane_z
        hit.terminal = True
        hit.direction = -1          # crossing downward
        events = hit

    sol = solve_ivp(rhs, (0.0, t_max), np.concatenate([x0, v0]),
                    method="RK45", max_step=max_step, rtol=rtol, atol=atol,
                    events=events, dense_output=False)
    X = sol.y[:3].T
    V = sol.y[3:].T
    return sol.t, X, V


# --------------------------------------------------------------------------- #
# Stain data structure.
# --------------------------------------------------------------------------- #
@dataclass
class Stain:
    x: float            # landing position on the target surface (m)
    y: float
    z: float            # surface height at the stain (0 for floor; wall stains carry z)
    L: float            # ellipse major axis (m)
    W: float            # ellipse minor axis (m)
    azimuth: float      # in-surface direction of travel, radians (points AWAY from origin)
    d0: float           # parent drop diameter (m); measured/estimated in practice

    @property
    def alpha(self) -> float:
        """Impact angle from Balthazard: alpha = arcsin(W / L)."""
        return np.arcsin(np.clip(self.W / self.L, 0.0, 1.0))


# --------------------------------------------------------------------------- #
# FORWARD MODEL: known origin -> synthetic stains (for validation).
# --------------------------------------------------------------------------- #
def emit_and_land(origin, n_drops, speed_range, d0_range, cone_half_angle,
                  floor_z=0.0, seed=0, drag=True):
    """Emit n_drops from `origin` into a downward-ish cone, forward-integrate to
    the floor, and synthesise the resulting elliptical stains.

    Returns list[Stain]. Droplets that miss the floor within t_max are dropped.
    """
    rng = np.random.default_rng(seed)
    origin = np.asarray(origin, float)
    stains = []
    for _ in range(n_drops):
        speed = rng.uniform(*speed_range)
        d0 = rng.uniform(*d0_range)
        # Launch direction: azimuth uniform, polar within a cone about horizontal,
        # biased slightly downward (typical of impact/cast-off spatter).
        az = rng.uniform(0.0, 2.0 * np.pi)
        elev = rng.uniform(-cone_half_angle, cone_half_angle)  # above/below horizontal
        dir_vec = np.array([np.cos(elev) * np.cos(az),
                            np.cos(elev) * np.sin(az),
                            np.sin(elev)])
        v0 = speed * dir_vec
        t, X, V = integrate_flight(origin, v0, d0, direction=+1, drag=drag,
                                   stop_plane_z=floor_z)
        if X[-1, 2] > floor_z + 1e-6:      # never reached floor
            continue
        p = X[-1]
        v_imp = V[-1]
        speed_imp = np.linalg.norm(v_imp)
        # Impact angle: between velocity and the (horizontal) surface.
        v_n = abs(v_imp[2])                       # normal component (floor)
        alpha = np.arcsin(np.clip(v_n / speed_imp, 0.0, 1.0))
        if np.degrees(alpha) < 3.0:               # grazing: ellipse method invalid
            continue
        # Synthesize stain size from drop-impact spread of the NORMAL velocity.
        we_n, re_n = _we_re_from_vn(v_n, d0)
        xi = float(spread_factor(we_n, re_n))
        W = d0 * xi                               # minor axis ~ normal-driven spread
        L = W / np.sin(alpha)                     # major axis from impact angle
        # In-surface direction of travel = horizontal projection of impact velocity.
        azimuth = np.arctan2(v_imp[1], v_imp[0])
        stains.append(Stain(x=p[0], y=p[1], z=floor_z, L=L, W=W,
                            azimuth=azimuth, d0=d0))
    return stains


# --------------------------------------------------------------------------- #
# RECONSTRUCTION 1: classical straight-line "stringing".
# --------------------------------------------------------------------------- #
def convergence_2d(stains):
    """Least-squares intersection of the in-surface direction-of-travel lines.

    Each stain contributes a line through (x,y) with in-surface direction
    (cos az, sin az). The convergence point minimises the sum of squared
    perpendicular distances to all back-directed lines (closed-form SVD-free
    normal-equations solution). Returns (x0, y0).
    """
    A = np.zeros((2, 2))
    b = np.zeros(2)
    for s in stains:
        d = np.array([np.cos(s.azimuth), np.sin(s.azimuth)])
        n = np.array([-d[1], d[0]])               # unit normal to the line
        p = np.array([s.x, s.y])
        A += np.outer(n, n)
        b += np.outer(n, n) @ p
    x0y0 = np.linalg.solve(A, b)
    return float(x0y0[0]), float(x0y0[1])


def straight_line_height(stains, x0, y0):
    """Classic tangent method: z = d * tan(alpha), averaged over stains, where d
    is the horizontal distance from stain to the 2-D convergence point. Ignores
    gravity -> biased high."""
    zs = []
    for s in stains:
        d = np.hypot(s.x - x0, s.y - y0)
        zs.append(s.z + d * np.tan(s.alpha))
    zs = np.asarray(zs)
    return float(np.mean(zs)), float(np.std(zs)), zs


# --------------------------------------------------------------------------- #
# RECONSTRUCTION 2: ballistic back-projection (gravity + drag).
# --------------------------------------------------------------------------- #
def impact_velocity_vector(stain: Stain):
    """Recover the full 3-D impact velocity vector at a stain from its geometry
    and the drop-impact spread model. Points in the direction of travel (into
    the surface). Returns (v_vector, speed, alpha)."""
    alpha = stain.alpha
    xi = stain.W / stain.d0
    v_n = vn_from_spread(xi, stain.d0)            # normal speed from spread factor
    if v_n <= 0.0:
        return None, 0.0, alpha
    speed = v_n / np.sin(alpha)                   # total speed from impact angle
    v_horiz = speed * np.cos(alpha)
    vx = v_horiz * np.cos(stain.azimuth)
    vy = v_horiz * np.sin(stain.azimuth)
    vz = -v_n                                     # downward onto floor
    return np.array([vx, vy, vz]), speed, alpha


def _point_to_polyline_dist2(p, poly):
    """Min squared distance from point p to a polyline (segment-aware)."""
    a = poly[:-1]
    b = poly[1:]
    ab = b - a
    ap = p - a
    denom = np.einsum("ij,ij->i", ab, ab)
    t = np.clip(np.einsum("ij,ij->i", ap, ab) / np.maximum(denom, 1e-18), 0.0, 1.0)
    proj = a + t[:, None] * ab
    d2 = np.einsum("ij,ij->i", p - proj, p - proj)
    return float(d2.min())


def ballistic_origin(stains, drag=True, guess=None, t_back=0.9):
    """Back-integrate each stain's trajectory and find the 3-D point of closest
    approach to all of them (the region of origin).

    t_back (s) need only exceed the true flight time; 0.9 s covers heights up to
    ~2 m. Reconstruction is at the cm scale, so relaxed ODE tolerances suffice.

    Returns (origin_xyz, per_stain_min_dist, trajectories).
    """
    trajs = []
    for s in stains:
        v_imp, speed, alpha = impact_velocity_vector(s)
        if v_imp is None:
            continue
        # Backward integration: reverse the impact velocity as the "initial"
        # backward velocity; integrate_flight(direction=-1) undoes gravity+drag.
        x_imp = np.array([s.x, s.y, s.z])
        t, X, V = integrate_flight(x_imp, v_imp, s.d0, direction=-1, drag=drag,
                                   t_max=t_back, max_step=6e-3,
                                   rtol=1e-6, atol=1e-9)
        trajs.append(X)

    if guess is None:
        x0, y0 = convergence_2d(stains)
        z0, _, _ = straight_line_height(stains, x0, y0)
        guess = np.array([x0, y0, z0])

    def cost(p):
        return sum(_point_to_polyline_dist2(p, tr) for tr in trajs)

    res = minimize(cost, np.asarray(guess, float), method="Nelder-Mead",
                   options={"xatol": 1e-4, "fatol": 1e-8, "maxiter": 20000})
    origin = res.x
    dists = np.sqrt([_point_to_polyline_dist2(origin, tr) for tr in trajs])
    return origin, dists, trajs
