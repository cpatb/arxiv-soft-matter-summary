"""
uroflow.py -- Splashback physics of a directed liquid jet (uroflow), v2.

WHAT CHANGED FROM v1 (2026-07-08 -> 2026-08-07).  v1 is kept in _archive/.
Five defects were found and fixed; each is marked [FIX n] at the point of use.

  [FIX 1] TARGET GEOMETRY WAS SELF-INCONSISTENT.  v1 integrated the trajectory to a
          VERTICAL plane at x = horizontal_target, then scored splash with
          sin(angle-from-horizontal).  For a vertical wall the normal is horizontal, so
          sin(angle-from-horizontal) is the TANGENTIAL fraction: v1 scored a purely
          normal impact (the worst case) as the driest.  The target is now an explicit
          object carrying a unit normal, and every splash number is built on the
          resolved NORMAL velocity component.

  [FIX 2] THE JET WAS TREATED AS PARALLEL.  A falling jet accelerates and thins by
          continuity: r(s) = r0*sqrt(v0/v(s)).  Over a 0.75 m fall the diameter drops
          by ~37%, which changes both the local Rayleigh growth rate (omega ~ r^-3/2)
          and the impact Weber number (We ~ d).  v1 used the exit diameter everywhere
          and therefore over-stated We at impact by ~1.6x.

  [FIX 3] BREAKUP LENGTH WAS UNCALIBRATED AND ~3-4x TOO LONG.  v1 printed L_b = 63 cm
          next to its own note that real streams break at 15-20 cm, then used the 63 cm
          value to draw the "coherent" region that the headline recommendation depends
          on.  Breakup is now the solution of an accumulated-disturbance integral along
          the thinning jet, and the disturbance amplification Lambda = ln(r0/delta0) is
          CALIBRATED against the observed breakup distance instead of assumed.

  [FIX 4] POST-BREAKUP DROP SIZE WAS WRONG.  v1 used the jet diameter as the drop
          diameter in the impact Weber number.  Rayleigh breakup produces drops of
          D = (9/x_opt)^(1/3) * r, i.e. ~1.89 d_jet, so v1 under-stated the drop We by
          ~1.9x in exactly the regime it called "worst".

  [FIX 5] THE SPLASH INDEX WAS INVENTED.  sqrt(We)*sin(theta) with an ad-hoc logistic
          is monotone in the right variables but has no threshold, so it could rank but
          never predict.  Replaced by published criteria: Mundo et al. (1995) for a dry
          solid and Cossali et al. (1997) for a wetted surface, both evaluated on the
          normal component.  The model now answers "does it splash at all", not only
          "which is worse".

PHYSICS CHAIN
-------------
1. Exit.  Peak uroflow Q ~ 15-25 mL/s through an effective meatus diameter d ~ 3-6 mm
   gives v0 = Q/(pi d^2/4) ~ 1-4 m/s, We ~ 100-600, Oh ~ 0.002 (essentially inviscid).

2. Falling jet.  Energy along the centreline gives v(s); continuity gives r(s).  Both
   feed the local stability calculation and the impact conditions.

3. Plateau-Rayleigh-Weber breakup.  The local fastest-growing varicose mode has growth
   rate omega_max(r) from Weber's viscous long-wave dispersion relation.  The jet breaks
   where the accumulated amplification integral G(s) = int_0^s omega(s')/v(s') ds'
   reaches Lambda = ln(r0/delta0).

4. Trajectory.  Ballistic centreline for the coherent column (drag on a coherent
   cm-scale stream is negligible, and the model now reports the number rather than
   asserting it).  After breakup the fragments are integrated with Schiller-Naumann
   drag, so the drag correction is computed instead of assumed.

5. Impact.  The target carries a unit normal.  Velocity is resolved into normal and
   tangential parts and the splash criterion is evaluated on the normal part.

ASSUMPTIONS, STATED
-------------------
* Axisymmetric jet.  The meatus is a slit and real streams split; the effective circular
  diameter lumps that in.  This is the largest unmodelled geometric effect.
* Newtonian, water-like fluid.  Urine is 1.00-1.10 mPa.s and 66 mN/m at body temperature.
* No air entrainment, no evaporation, isothermal.
* Stretching stabilisation of the falling jet is NOT modelled.  A gravitationally
  stretched jet is more stable than a parallel one (Frankel & Weihs 1985), so the
  computed breakup length is a LOWER bound on the ideal value.  Since the calibrated
  Lambda absorbs the difference, this matters for interpreting Lambda, not for the
  predicted breakup point.

SI units throughout.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

# --------------------------------------------------------------------------- #
# Material properties, with provenance.
# --------------------------------------------------------------------------- #
RHO = 1020.0          # kg/m^3   urine, 1.005-1.030 specific gravity
MU = 1.0e-3           # Pa.s     1.0-1.1 mPa.s, near water at 37 C
SIGMA = 0.066         # N/m      66 mN/m, slightly below water (urea + salts)
G = 9.81              # m/s^2

RHO_AIR = 1.20        # kg/m^3   20 C, 1 atm
MU_AIR = 1.81e-5      # Pa.s

# Weber's noise-amplification factor Lambda = ln(r0/delta0).  The classical textbook
# value for a quiescent laboratory jet is ~12.  A meatus jet is strongly perturbed, so
# the default here is CALIBRATED (see calibrate_lambda) rather than assumed.  12.0 is
# retained only as the "ideal laboratory jet" reference case.
LAMBDA_IDEAL = 12.0

# Observed breakup distance for a real human stream, used as the calibration anchor.
# Hurd & Truscott (Splash Lab), "Urinal Dynamics", APS-DFD 2013: coherent length
# ~15-20 cm.  Stored as a range so the calibration can report a band, not a point.
L_BREAKUP_OBSERVED = (0.15, 0.20)   # m

# Watson region-I / region-II transition, reused from the validated bidet wall-jet
# model (Projects/bidet/model/bidet_model.py, r0_over_a): r0/a = 0.330 Re_a^(1/3).
C_DELTA = np.sqrt(280.0 / 39.0)


def _check_positive(**kwargs):
    """Fail loudly and early on non-physical inputs."""
    for name, value in kwargs.items():
        arr = np.asarray(value, dtype=float)
        if not np.all(np.isfinite(arr)):
            raise ValueError(f"{name} must be finite, got {value!r}")
        if np.any(arr <= 0.0):
            raise ValueError(f"{name} must be strictly positive, got {value!r}")


# --------------------------------------------------------------------------- #
# 1. Exit conditions and dimensionless groups.
# --------------------------------------------------------------------------- #
def exit_velocity(Q, d):
    """Exit speed (m/s) from volumetric flow Q (m^3/s) and diameter d (m)."""
    _check_positive(Q=Q, d=d)
    return Q / (np.pi * d**2 / 4.0)


def weber(v, d, rho=RHO, sigma=SIGMA):
    """We = rho v^2 d / sigma."""
    return rho * np.asarray(v, dtype=float) ** 2 * d / sigma


def reynolds(v, d, rho=RHO, mu=MU):
    """Re = rho v d / mu."""
    return rho * np.asarray(v, dtype=float) * d / mu


def ohnesorge(d, rho=RHO, mu=MU, sigma=SIGMA):
    """Oh = mu / sqrt(rho sigma d) = sqrt(We)/Re."""
    _check_positive(d=d)
    return mu / np.sqrt(rho * sigma * d)


# --------------------------------------------------------------------------- #
# 2. Falling-jet kinematics: energy + continuity.  [FIX 2]
# --------------------------------------------------------------------------- #
def jet_speed_after_drop(v0, fall):
    """Centreline speed (m/s) after falling a height `fall` (m, positive downward).

    Energy along a free streamline, neglecting air drag on the column and viscous
    losses (Oh ~ 2e-3, so viscous dissipation over a sub-metre fall is negligible).
    `fall` may be negative if the stream rises; the speed then decreases.
    """
    v0 = np.asarray(v0, dtype=float)
    v_sq = v0**2 + 2.0 * G * np.asarray(fall, dtype=float)
    if np.any(v_sq <= 0.0):
        raise ValueError("jet decelerates to rest before the target (rise too large)")
    return np.sqrt(v_sq)


def jet_radius(r0, v0, v):
    """Local jet radius (m) from continuity, Q = pi r^2 v = const."""
    _check_positive(r0=r0, v0=v0, v=v)
    return r0 * np.sqrt(np.asarray(v0, dtype=float) / np.asarray(v, dtype=float))


# --------------------------------------------------------------------------- #
# 3. Local stability and the accumulated-disturbance breakup criterion.
# --------------------------------------------------------------------------- #
def max_growth_rate(d, rho=RHO, mu=MU, sigma=SIGMA):
    """Fastest-growing varicose growth rate omega_max (1/s) and its x = k r.

    Weber's viscous long-wave dispersion relation on a jet of radius r = d/2:

        omega(x) = -1.5 nu k^2 + sqrt[ (1.5 nu k^2)^2 + (sigma/(2 rho r^3)) x^2 (1-x^2) ]

    CLOSED FORM.  Writing A = 1.5 nu / r^2, B = sigma/(2 rho r^3) and y = x^2,

        omega(y) = -A y + sqrt[ (A^2 - B) y^2 + B y ],

    and d(omega)/dy = 0 gives, after squaring once and cancelling a factor of B,

        4 (A^2 - B) y^2 + 4 B y - B = 0   ->   y_opt = sqrt(B) / (2 (A + sqrt(B))).

    The inviscid limit A -> 0 returns y = 1/2, i.e. x = 1/sqrt(2), which is Rayleigh's
    classical optimum; omega_max is then sqrt(sigma/(8 rho r^3)).  The root was verified
    symbolically and is cross-checked against the numerical optimiser in the tests.

    Vectorised over d, which matters: the falling-jet integral evaluates this at every
    station along the jet, and the scalar optimiser made that step ~1000x slower.
    """
    d = np.asarray(d, dtype=float)
    _check_positive(d=d)
    r = d / 2.0
    nu = mu / rho
    A = 1.5 * nu / r**2
    B = sigma / (2.0 * rho * r**3)
    sqrtB = np.sqrt(B)
    y = sqrtB / (2.0 * (A + sqrtB))
    omega = -A * y + np.sqrt(np.maximum((A**2 - B) * y**2 + B * y, 0.0))
    x = np.sqrt(y)
    if omega.ndim == 0:
        return float(omega), float(x)
    return omega, x


def max_growth_rate_numeric(d, rho=RHO, mu=MU, sigma=SIGMA):
    """Brute-force optimisation of the same dispersion relation, for cross-checking the
    closed form above.  Not used in the hot path."""
    _check_positive(d=d)
    r = float(d) / 2.0
    nu = mu / rho
    base = sigma / (2.0 * rho * r**3)

    def neg_omega(x):
        k = x / r
        visc = 1.5 * nu * k**2
        disc = visc**2 + base * x**2 * (1.0 - x**2)
        return -(-visc + np.sqrt(max(disc, 0.0)))

    res = minimize_scalar(neg_omega, bounds=(1e-6, 1.0 - 1e-9), method="bounded",
                          options={"xatol": 1e-12})
    return -res.fun, res.x


def rayleigh_drop_diameter(d_jet, x_opt=None):
    """Diameter (m) of the drops produced by varicose breakup of a jet of diameter d_jet.

    One wavelength lambda = 2 pi r / x_opt of jet collapses into one sphere, so
    (pi/6) D^3 = pi r^2 lambda  ->  D = (12 pi / x_opt)^(1/3) * r * ... expanded below.
    For the inviscid optimum x_opt = 1/sqrt(2) this gives D = 1.89 d_jet.  [FIX 4]
    """
    _check_positive(d_jet=d_jet)
    if x_opt is None:
        _, x_opt = max_growth_rate(d_jet)
    r = d_jet / 2.0
    wavelength = 2.0 * np.pi * r / x_opt
    volume = np.pi * r**2 * wavelength
    return (6.0 * volume / np.pi) ** (1.0 / 3.0)


def disturbance_integral(v0, r0, path_s, fall_of_s, n=400, rho=RHO, mu=MU, sigma=SIGMA):
    """Accumulated logarithmic disturbance growth G(s) along the jet.  [FIX 3]

    G(s) = int_0^s omega_max(r(s')) / v(s') ds'

    where v and r follow the falling-jet relations of section 2.  The jet breaks where
    G reaches Lambda = ln(r0/delta0).  For a parallel, non-accelerating jet this reduces
    exactly to the textbook L_b = v/omega * Lambda, which the tests check.

    Parameters
    ----------
    path_s : float          arclength to integrate over (m)
    fall_of_s : callable    fall(s) -> vertical drop (m, positive down) at arclength s

    Returns
    -------
    s_grid, G_grid : ndarray
    """
    _check_positive(v0=v0, r0=r0, path_s=path_s)
    s = np.linspace(0.0, float(path_s), int(n))
    fall = np.asarray([float(fall_of_s(si)) for si in s])
    v = jet_speed_after_drop(v0, fall)
    r = jet_radius(r0, v0, v)
    omega, _ = max_growth_rate(2.0 * r, rho, mu, sigma)
    integrand = omega / v
    G = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(s))])
    return s, G


def breakup_arclength(v0, r0, fall_of_s, lam, s_max=2.0, **kw):
    """Arclength (m) at which the jet fragments, or np.inf if it survives to s_max."""
    s, G = disturbance_integral(v0, r0, s_max, fall_of_s, **kw)
    if G[-1] < lam:
        return np.inf
    return float(np.interp(lam, G, s))


def breakup_length_parallel(Q, d, lam=LAMBDA_IDEAL):
    """Textbook parallel-jet breakup length L_b = v/omega_max * Lambda (m).

    Retained for comparison and for the tests; the falling-jet version above supersedes
    it for anything with a vertical drop.
    """
    v = exit_velocity(Q, d)
    omega, _ = max_growth_rate(d)
    return v / omega * lam


def calibrate_lambda(Q, d, launch_height, observed_range=L_BREAKUP_OBSERVED):
    """Fit Lambda = ln(r0/delta0) so the model breaks where real streams are seen to.

    Returns (lam_lo, lam_hi, delta_over_r_lo, delta_over_r_hi) for the observed band.
    A large Lambda means a quiet jet; a small Lambda means a strongly perturbed exit.
    Reporting the implied delta0/r is the honest way to say how disturbed the model
    thinks the meatus is.
    """
    _check_positive(Q=Q, d=d)
    v0 = exit_velocity(Q, d)
    r0 = d / 2.0

    # Straight-down stream, so arclength equals fall.  This is the geometry the
    # Splash Lab anchor was measured in.
    def fall_of_s(s):
        return s

    out = []
    for L_obs in observed_range:
        s, Gc = disturbance_integral(v0, r0, max(2.0 * L_obs, 0.05), fall_of_s)
        out.append(float(np.interp(L_obs, s, Gc)))
    lam_lo, lam_hi = min(out), max(out)
    return lam_lo, lam_hi, float(np.exp(-lam_hi)), float(np.exp(-lam_lo))


# --------------------------------------------------------------------------- #
# 4. Target geometry.  [FIX 1]
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Target:
    """An impact surface in the (x, z) plane, with an explicit outward unit normal.

    kind : 'pool'  horizontal free surface at z = z_plane, normal (0, 1)
           'wall'  vertical surface at x = x_plane, normal (-1, 0)
    film : liquid depth on the surface (m).  0.0 means dry.  For 'pool' this is the
           standing water depth, which is the dominant design lever in the real problem
           and was absent from v1 entirely.
    """
    kind: str
    plane: float
    film: float = 0.0

    def __post_init__(self):
        if self.kind not in ("pool", "wall"):
            raise ValueError(f"kind must be 'pool' or 'wall', got {self.kind!r}")
        if self.film < 0.0:
            raise ValueError("film depth cannot be negative")

    @property
    def normal(self):
        """Outward unit normal as (nx, nz)."""
        return np.array([0.0, 1.0]) if self.kind == "pool" else np.array([-1.0, 0.0])


def resolve_impact(vx, vz, target):
    """Split an impact velocity into normal and tangential magnitudes (m/s).

    Splash thresholds for oblique impact are conventionally evaluated on the NORMAL
    component; the tangential part slides the lamella rather than ejecting it.
    """
    v = np.array([float(vx), float(vz)])
    n = target.normal
    v_n = abs(float(np.dot(v, n)))
    v_t = float(np.linalg.norm(v - np.dot(v, n) * n))
    return v_n, v_t


# --------------------------------------------------------------------------- #
# 5. Trajectory.
# --------------------------------------------------------------------------- #
class StreamMissesTarget(ValueError):
    """The stream does not reach the requested target: it falls short, hits the floor
    first, or is aimed away from it.  A distinct type so callers can skip these
    configurations without swallowing genuine input errors."""


def ballistic_to_target(v0, theta, h, target, n=2000, floor=0.0):
    """Ballistic centreline from (0, h) at speed v0 and angle theta (rad, negative =
    aimed downward) until it meets `target`.

    `floor` is the ground plane.  A wall-bound stream that drops below the floor before
    reaching the wall has landed short, which is a different configuration, not a
    grazing impact.  v2.0 of this function omitted that check and reported nonsensical
    steep-aim results; the guard below is the fix.

    Returns a dict with the arclength travelled, the impact point and velocity, and the
    vertical fall at impact.  Raises StreamMissesTarget if the stream never gets there.
    """
    _check_positive(v0=v0)
    vx0 = v0 * np.cos(theta)
    vz0 = v0 * np.sin(theta)
    if vx0 <= 1e-9 and target.kind == "wall":
        raise StreamMissesTarget("a stream with no horizontal velocity never reaches a wall")

    if target.kind == "wall":
        t_end = target.plane / vx0
        # Time at which the parabola crosses the floor, if it does.
        disc_floor = vz0 * vz0 + 2.0 * G * (h - floor)
        if disc_floor > 0.0:
            t_floor = (vz0 + np.sqrt(disc_floor)) / G
            if t_floor < t_end:
                raise StreamMissesTarget(
                    f"stream lands {target.plane - vx0 * t_floor:.3f} m short of the wall")
    else:
        # Solve h + vz0 t - g t^2/2 = z_plane for the later root.
        a, b, c = -0.5 * G, vz0, h - target.plane
        disc = b * b - 4 * a * c
        if disc < 0:
            raise StreamMissesTarget("stream never reaches the pool plane")
        t_end = max((-b + np.sqrt(disc)) / (2 * a), (-b - np.sqrt(disc)) / (2 * a))
    if t_end <= 0:
        raise StreamMissesTarget("target is behind the launch point")

    t = np.linspace(0.0, t_end, int(n))
    x = vx0 * t
    z = h + vz0 * t - 0.5 * G * t**2
    s = float(np.sum(np.hypot(np.diff(x), np.diff(z))))
    vz_end = vz0 - G * t_end

    # Arclength as a function of time, so a breakup arclength can be converted back.
    s_of_t = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(x), np.diff(z)))])

    return dict(t=t, x=x, z=z, s=s, s_of_t=s_of_t,
                vx_impact=vx0, vz_impact=float(vz_end),
                speed_impact=float(np.hypot(vx0, vz_end)),
                fall_at_impact=float(h - z[-1]),
                fall_of_s=lambda ss: float(np.interp(ss, s_of_t, h - z)))


def drop_flight_with_drag(vx, vz, x0, z0, D, target, max_t=2.0):
    """Integrate a single fragment of diameter D from (x0, z0) to the target with
    Schiller-Naumann sphere drag, so the drag correction is COMPUTED, not asserted.

    Returns the impact velocity components and the fractional speed change relative to
    the drag-free result.
    """
    _check_positive(D=D)

    def cd(re):
        re = max(float(re), 1e-6)
        return 24.0 / re * (1.0 + 0.15 * re**0.687) if re < 1000.0 else 0.44

    def rhs(_t, y):
        _x, _z, u, w = y
        speed = np.hypot(u, w)
        if speed < 1e-9:
            return [u, w, 0.0, -G]
        re = RHO_AIR * speed * D / MU_AIR
        # a = -(3 rho_air Cd /(4 rho_l D)) |v| v
        k = 3.0 * RHO_AIR * cd(re) / (4.0 * RHO * D)
        return [u, w, -k * speed * u, -G - k * speed * w]

    def hit(_t, y):
        return y[0] - target.plane if target.kind == "wall" else y[1] - target.plane
    hit.terminal = True
    hit.direction = 0

    sol = solve_ivp(rhs, (0.0, max_t), [x0, z0, vx, vz], events=hit,
                    rtol=1e-8, atol=1e-10, dense_output=True)
    if not sol.t_events[0].size:
        raise StreamMissesTarget("fragment never reached the target within max_t")
    y_end = sol.y_events[0][0]

    # Drag-free comparison over the same flight time.
    t_end = float(sol.t_events[0][0])
    vx_free, vz_free = vx, vz - G * t_end
    speed_free = np.hypot(vx_free, vz_free)
    speed_drag = np.hypot(y_end[2], y_end[3])
    return dict(vx=float(y_end[2]), vz=float(y_end[3]),
                speed=float(speed_drag), t_flight=t_end,
                drag_speed_deficit=float(1.0 - speed_drag / speed_free))


# --------------------------------------------------------------------------- #
# 6. Splash criteria: published thresholds, evaluated on the normal component. [FIX 5]
# --------------------------------------------------------------------------- #
K_DRY_CRITICAL = 57.7        # Mundo, Sommerfeld & Tropea (1995), Int. J. Multiphase Flow


def splash_dry(v_n, D):
    """Mundo et al. (1995) dry-surface splash parameter and its threshold.

    K = We^(1/2) Re^(1/4); deposition below K ~ 57.7, splash above.  Roughness and gas
    pressure move this by order-unity factors, so it is an order-of-magnitude criterion
    and is reported as an excess ratio rather than a yes/no.
    """
    _check_positive(D=D)
    if v_n <= 0.0:
        return dict(K=0.0, K_c=K_DRY_CRITICAL, excess=-1.0, splashes=False)
    We = weber(v_n, D)
    Re = reynolds(v_n, D)
    K = float(np.sqrt(We) * Re**0.25)
    return dict(K=K, K_c=K_DRY_CRITICAL, excess=K / K_DRY_CRITICAL - 1.0,
                splashes=bool(K > K_DRY_CRITICAL), We=float(We), Re=float(Re))


def splash_film(v_n, D, film):
    """Cossali, Coghe & Marengo (1997) wetted-surface splash criterion.

    K = We * Oh^(-0.4) = We^0.8 Re^0.4, with threshold
        K_c = 2100 + 5880 * delta^1.44,  delta = film/D,
    correlated over roughly 0.1 <= delta <= 1.  Outside that range the threshold is
    extrapolated and flagged, because a silently extrapolated correlation is the kind of
    thing that turns into a wrong headline.
    """
    _check_positive(D=D)
    if film <= 0.0:
        raise ValueError("splash_film requires a positive film depth; use splash_dry")
    if v_n <= 0.0:
        return dict(K=0.0, K_c=np.inf, excess=-1.0, splashes=False, delta=film / D,
                    extrapolated=False)
    delta = film / D
    We = weber(v_n, D)
    Oh = ohnesorge(D)
    K = float(We * Oh**-0.4)
    K_c = 2100.0 + 5880.0 * delta**1.44
    return dict(K=K, K_c=float(K_c), excess=K / K_c - 1.0, splashes=bool(K > K_c),
                delta=float(delta), We=float(We), Oh=float(Oh),
                extrapolated=bool(delta < 0.1 or delta > 1.0))


def watson_transition_radius(a, v):
    """Radius (m) where the viscous boundary layer fills the radial film of an impinging
    coherent jet: r0/a = 0.330 Re_a^(1/3).

    Reused from the validated bidet wall-jet model (Projects/bidet/model/bidet_model.py).
    Inside r0 the film thins as mass conservation dictates; outside it thickens.
    """
    _check_positive(a=a, v=v)
    nu = MU / RHO
    Re_a = v * a / nu
    return a * (Re_a / (4.0 * C_DELTA**2)) ** (1.0 / 3.0)


def coherent_film_thickness(Q, r, v):
    """Inviscid radial-film thickness h(r) = Q/(2 pi r v) for a coherent impinging jet.

    Valid inside the Watson transition radius, where the film has not yet thickened
    viscously.  Used to report what a coherent column actually deposits, instead of
    scoring it with a drop criterion that does not apply to a continuous column.
    """
    _check_positive(Q=Q, r=r, v=v)
    return Q / (2.0 * np.pi * r * v)


# --------------------------------------------------------------------------- #
# 7. The composite evaluation.
# --------------------------------------------------------------------------- #
def evaluate(Q, d, theta, h, target, lam, use_drag=True):
    """Full chain for one configuration.

    Returns a dict with the regime ('coherent' or 'fragmented'), the impact conditions
    resolved on the target normal, and the applicable published splash criterion.

    `excess` is the headline number: K/K_c - 1.  Negative means no splash is predicted;
    positive is a splash whose magnitude is above threshold by that fraction.  Unlike
    v1's index it is comparable across geometries and has a meaningful zero.
    """
    _check_positive(Q=Q, d=d, lam=lam)
    v0 = exit_velocity(Q, d)
    r0 = d / 2.0
    traj = ballistic_to_target(v0, theta, h, target)

    s_break = breakup_arclength(v0, r0, traj["fall_of_s"], lam,
                                s_max=max(traj["s"] * 1.5, 0.05))
    coherent = s_break >= traj["s"]

    v_local = jet_speed_after_drop(v0, traj["fall_at_impact"])
    d_local = 2.0 * jet_radius(r0, v0, v_local)

    if coherent:
        v_n, v_t = resolve_impact(traj["vx_impact"], traj["vz_impact"], target)
        a_jet = d_local / 2.0
        r_watson = watson_transition_radius(a_jet, max(v_n, 1e-6))
        h_film = coherent_film_thickness(Q, r_watson, max(v_n, 1e-6))
        crit = dict(regime="coherent", K=np.nan, K_c=np.nan, excess=np.nan,
                    splashes=False, watson_radius=float(r_watson),
                    film_at_watson=float(h_film))
        D_drop = np.nan
    else:
        # Fragment where the jet broke, then fly the drop the rest of the way.
        d_break = 2.0 * jet_radius(
            r0, v0, jet_speed_after_drop(v0, traj["fall_of_s"](s_break)))
        D_drop = rayleigh_drop_diameter(d_break)
        i_break = int(np.searchsorted(traj["s_of_t"], s_break))
        i_break = min(max(i_break, 0), len(traj["t"]) - 1)
        vz_at_break = v0 * np.sin(theta) - G * traj["t"][i_break]
        if use_drag:
            fl = drop_flight_with_drag(v0 * np.cos(theta), vz_at_break,
                                       traj["x"][i_break], traj["z"][i_break],
                                       D_drop, target)
            vx_i, vz_i, drag_deficit = fl["vx"], fl["vz"], fl["drag_speed_deficit"]
        else:
            vx_i, vz_i = traj["vx_impact"], traj["vz_impact"]
            drag_deficit = 0.0
        v_n, v_t = resolve_impact(vx_i, vz_i, target)
        crit = (splash_film(v_n, D_drop, target.film) if target.film > 0.0
                else splash_dry(v_n, D_drop))
        crit["regime"] = "fragmented"
        crit["drag_speed_deficit"] = drag_deficit

    return dict(regime=crit["regime"], excess=crit.get("excess", np.nan),
                splashes=crit.get("splashes", False),
                v_normal=float(v_n), v_tangential=float(v_t),
                path=traj["s"], breakup_arclength=s_break,
                coherent=bool(coherent), d_at_impact=float(d_local),
                drop_diameter=float(D_drop) if np.isfinite(D_drop) else np.nan,
                incidence_deg=float(np.degrees(np.arctan2(v_t, max(v_n, 1e-12)))),
                criterion=crit, trajectory=traj)
