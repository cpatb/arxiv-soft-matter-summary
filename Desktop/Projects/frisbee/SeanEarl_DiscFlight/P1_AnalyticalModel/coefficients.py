"""
Aerodynamic-coefficient sets for the Ultrastar 175 g flying disc.

PROVENANCE
----------
The functional forms and coefficient values here are taken directly from the
primary source:

    Hummel, S. A. (2003). "Frisbee Flight Simulation and Throw Biomechanics."
    MS Thesis, University of California, Davis.
    - Eqs. (3.3)-(3.8), p. 35        : force/moment definitions
    - Table 3-1, p. 35               : mass, inertia, geometry, rho, g
    - Table 3-2, p. 39               : fitted coefficient values
    - Appendix B, `discfltEOM.m`     : the reference implementation that was
                                       used to FIT the Table 3-2 values

The coefficient values in Table 3-2 are only meaningful inside the exact
functional forms of `discfltEOM.m`, because they were estimated by fitting that
code to measured flights. We therefore reproduce those forms literally rather
than substituting a "tidier" non-dimensionalisation. Hummel's relevant lines:

    CL  = CLo + CLa*alpha;
    alphaeq = -CLo/CLa;                    % angle of attack at zero lift
    CD  = CDo + CDa*(alpha-alphaeq)*(alpha-alphaeq);
    CM  = CMo + CMa*alpha;
    Mvp  = Adp*d*(sqrt(d/g)*CRr*omegaspin + CRp*omegavp)*uvp;   % roll
    Mlat = Adp*d*(CM + CMq*omegalat)*ulat;                      % pitch
    Mspin = [0 0 +CNr*(omegaspin)];                             % spin-down

Note three consequences that are easy to get wrong:
  1. alpha0 is NOT a free parameter. It is the zero-lift angle -CL0/CLa, and it
     appears ONLY in the drag parabola. It does not shift CL or CM.
  2. The rate-dependent terms (CMq, CRp, CRr) multiply RAW angular rates in
     rad/s. They are NOT non-dimensionalised by d/(2V); the CRr term carries an
     explicit sqrt(d/g) time scale instead. Non-dimensionalising them by d/(2V)
     under-damps pitch and roll by roughly two orders of magnitude.
  3. The spin-down moment `Mspin` is returned directly in N*m and is NOT scaled
     by the dynamic-pressure factor Adp*d, unlike the roll and pitch moments.

CROSS-CHECKS against independent statements in the same thesis (see
`validation_hummel.py` for the executable versions):
  - Sec. 2.4.1: "when alpha = 0 the Frisbee must be thrown at 16 m/s to produce
    horizontal flight."  P&C CL(0) = CL0 = 0.20 gives v = 15.6 m/s.
  - Sec. 2.4.2: pitching moment vanishes at alpha = 9 deg.  P&C
    CM = -0.02 + 0.13*alpha = 0 at alpha = 8.81 deg.
  - Sec. 2.4.2: at that trim, CL = 0.67.  P&C gives CL(8.81 deg) = 0.665.
These three independent checks only hold with CL and CM taken as functions of
alpha itself (not alpha - alpha0), which is why the forms above are used.

SIGN / FRAME CONVENTION
-----------------------
Hummel works in an inertial frame with z pointing DOWN and body axis c3 through
the disc's underside. This module and `disc_6dof.py` use z UP with body axis b3
through the disc's top. The angle of attack has the same sense in both (alpha>0
means the relative wind strikes the underside). The pitch axis, however, is
ulat = c3 x uvp in Hummel's convention and b3 x b1_aero in ours, and these two
axes are anti-parallel. `disc_6dof.py` therefore applies the pitching and
rolling moments about the z-up axes with the sign flip handled explicitly and
documented there; the coefficient values below are Hummel's, unmodified.

CONVENTIONS USED BELOW
----------------------
    alpha : angle of attack [rad]; alpha > 0 = flow strikes bottom of disc.
    p     : roll rate about the projected-velocity axis   [rad/s]
    q     : pitch rate about the in-plane lateral axis    [rad/s]
    r     : spin rate about the disc normal               [rad/s]

    CL(alpha)      = CL0 + CLa*alpha
    CD(alpha)      = CD0 + CDa*(alpha - alpha0)**2,   alpha0 = -CL0/CLa
    CM(alpha, q)   = CM0 + CMa*alpha + CMq*q
    CR(p, r)       = sqrt(d/g)*CRr*r + CRp*p
    spin_down(r)   = CNr*r                            [N*m, unscaled]
"""

from dataclasses import dataclass
from math import pi, sqrt

import numpy as np

# Standard gravity used by Hummel (Table 3-1). Kept as a named constant because
# it enters the sqrt(d/g) time scale in the roll-moment term, not just weight.
G_HUMMEL = 9.794            # m/s^2

DEG = pi / 180.0


@dataclass(frozen=True)
class DiscCoefficients:
    """A linearised coefficient set in Hummel's (2003) functional form."""

    name: str
    # lift
    CL0: float
    CLa: float
    # drag
    CD0: float
    CDa: float
    # pitching moment
    CM0: float
    CMa: float
    CMq: float
    # rolling moment
    CRr: float
    CRp: float
    # spin-axis (spin-down) moment
    CNr: float
    # Range of angle of attack over which this linear fit is defensible.
    # Outside it the linear CL/CM extrapolation is not supported by the source
    # data (the real disc stalls near 25-30 deg). Reported by validation.
    alpha_valid_deg: tuple = (-10.0, 25.0)
    # Provenance note for the paper's methods section.
    source: str = ""

    @property
    def alpha0(self) -> float:
        """Zero-lift angle of attack [rad], Hummel's `alphaeq = -CLo/CLa`.

        This is a derived quantity, not a fitted one. It sets the vertex of the
        drag parabola and appears nowhere else.
        """
        return -self.CL0 / self.CLa

    def CL(self, alpha: float) -> float:
        """Lift coefficient. Hummel: CL = CLo + CLa*alpha."""
        return self.CL0 + self.CLa * alpha

    def CD(self, alpha: float) -> float:
        """Drag coefficient. Hummel: CD = CDo + CDa*(alpha-alphaeq)^2."""
        da = alpha - self.alpha0
        return self.CD0 + self.CDa * da * da

    def CM(self, alpha: float, q: float) -> float:
        """Pitching-moment coefficient, including rate damping.

        Hummel: CM = CMo + CMa*alpha, then Mlat = Adp*d*(CM + CMq*omegalat).
        `q` is the RAW pitch rate in rad/s about the in-plane lateral axis.
        """
        return self.CM0 + self.CMa * alpha + self.CMq * q

    def CR(self, p: float, r: float, d: float, g: float = G_HUMMEL) -> float:
        """Rolling-moment coefficient about the projected-velocity axis.

        Hummel: sqrt(d/g)*CRr*omegaspin + CRp*omegavp, with RAW rates in rad/s.
        The sqrt(d/g) factor (~0.166 s for an Ultrastar) is the time scale that
        makes the spin-coupling term dimensionally consistent.
        """
        return sqrt(d / g) * self.CRr * r + self.CRp * p

    def spin_down_moment(self, r: float) -> float:
        """Spin-down moment about the disc normal, in N*m.

        Hummel: Mspin = [0 0 CNr*omegaspin]. Unlike the roll and pitch moments
        this is NOT multiplied by the dynamic-pressure factor Adp*d, so it is
        returned as a moment rather than as a coefficient.
        """
        return self.CNr * r

    def alpha_out_of_range(self, alpha_deg: float) -> bool:
        lo, hi = self.alpha_valid_deg
        return alpha_deg < lo or alpha_deg > hi


# ---------------------------------------------------------------------------
# Published coefficient sets (Hummel 2003, Table 3-2)
# ---------------------------------------------------------------------------

# Row "short flights": all ten coefficients estimated from three atypical
# flights (v = 3 m/s, spin = 8.91 rad/s, distance = 2 m; thesis Sec. 3.3.2).
#
# Values are the full-precision ones from the thesis Appendix A comment block
# (p. 82), which Table 3-2 reports rounded to two significant figures.
#
# CRr WARNING: Table 3-2 lists CRr = 1.4e-2, but Appendix A carries
#   CRr = 0.00171;  % for nondimensionalization = sqrt(d/g), magnitude of CRr
#                   % changes with nondimensionalization
# The two numbers describe the same physics under different scalings. Because
# `CR()` above implements the sqrt(d/g) form that `discfltEOM.m` actually
# integrates, the matching value is 0.00171. Using Table 3-2's 1.4e-2 inside the
# sqrt(d/g) form overstates the spin-induced roll moment by roughly 8x.
HUMMEL_SHORT = DiscCoefficients(
    name="Hummel 2003 short flights",
    CL0=0.3331, CLa=1.9124,
    CD0=0.1769, CDa=0.685,
    CM0=-0.0821, CMa=0.4338, CMq=-5.0e-3,
    CRr=0.00171, CRp=-5.5e-3,
    CNr=-7.1e-6,
    source=("Hummel (2003) Appendix A (full precision); Table 3-2 row 'short "
            "flights' (rounded). CRr is the sqrt(d/g)-scaled value 0.00171."),
)

# Row "long flights": ONLY the three damping coefficients CMq, CRp, CNr were
# re-optimised (v = 14 m/s, spin = 50 rad/s, distance = 18 m). Hummel Sec.
# 3.3.4: "the seven [other] coefficients were specified according to the values
# reported by Potts and Crowther (2002)." The published set is therefore the
# P&C force/pitch coefficients PLUS these three damping values -- not the
# short-flight force coefficients. Getting this wrong produces a hybrid that
# corresponds to no published measurement.
HUMMEL_LONG = DiscCoefficients(
    name="Hummel 2003 long flights (P&C forces + refitted damping)",
    CL0=0.20, CLa=2.96,
    CD0=0.08, CDa=2.72,
    CM0=-0.02, CMa=0.13, CMq=-1.44e-2,
    CRr=0.00171, CRp=-1.25e-2,
    CNr=-3.41e-5,
    source=("Hummel (2003) Appendix A, 'long flight f2302' for CMq/CRp/CNr "
            "(Table 3-2 row 'long flights', rounded); force and pitch "
            "coefficients from Potts & Crowther (2002) per Sec. 3.3.4. "
            "CRr retained at Hummel's sqrt(d/g)-scaled 0.00171: Table 3-2 "
            "lists -3e-3 for P&C, but under an unstated scaling, so the "
            "internally consistent value is used and flagged here."),
)

# Row "P & C": Potts & Crowther (2002), AIAA 2002-3150, as tabulated by Hummel.
# The angular-rate coefficients were NOT measured in the wind tunnel (Hummel
# Sec. 3.4: "wind tunnel measurement techniques cannot measure a moment about an
# axis that is spinning"), so CMq/CRp/CNr are borrowed from the short-flight
# estimate and flagged. Table 3-2 lists CRr = -3e-3 for P&C.
POTTS_CROWTHER = DiscCoefficients(
    name="Potts & Crowther 2002 (wind tunnel)",
    CL0=0.20, CLa=2.96,
    CD0=0.08, CDa=2.72,
    CM0=-0.02, CMa=0.13, CMq=-5.0e-3,     # CMq borrowed: not measured by P&C
    CRr=-3.0e-3, CRp=-5.5e-3,             # CRp borrowed: not measured by P&C
    CNr=-7.1e-6,                          # borrowed: not measured by P&C
    source=("Potts & Crowther (2002) AIAA 2002-3150 as tabulated in Hummel "
            "(2003) Table 3-2. CMq, CRp, CNr are NOT from P&C (unmeasurable in "
            "a wind tunnel) and are carried over from Hummel's short flights."),
)

# ---------------------------------------------------------------------------
# Tabulated Potts & Crowther (2002) wind-tunnel curves
# ---------------------------------------------------------------------------
#
# Transcribed verbatim from Hummel (2003) Appendix A (pp. 82-83), where they are
# used by the `CoefUsed==2` branch of `discfltEOM.m`. These are the measured
# curves, so they replace the linear fits and remain valid through stall; they
# also carry the advance-ratio dependence of the spin-induced rolling moment,
# which the linear sets cannot represent.
#
# Columns: angle of attack [rad], coefficient [-].

_PC_CL = np.array([
    [-0.17450, -0.2250],
    [-0.05236,  0.0000],
    [ 0.00000,  0.1500],
    [ 0.08727,  0.4500],
    [ 0.17453,  0.7250],
    [ 0.26180,  0.9750],
    [ 0.34907,  1.2000],
    [ 0.43633,  1.4500],
    [ 0.52360,  1.6750],
])

_PC_CD = np.array([
    [-0.17450, 0.1500],
    [-0.05236, 0.0800],
    [ 0.00000, 0.1000],
    [ 0.08727, 0.1500],
    [ 0.17450, 0.2600],
    [ 0.26180, 0.3900],
    [ 0.34910, 0.5700],
    [ 0.43630, 0.7500],
    [ 0.52360, 0.9200],
])

_PC_CM = np.array([
    [-0.174532925, -0.0380],
    [-0.087266463, -0.0220],
    [-0.052359878, -0.0140],
    [ 0.000000000, -0.0060],
    [ 0.052359878, -0.0060],
    [ 0.104719755, -0.0020],
    [ 0.157079633,  0.0000],
    [ 0.209439510,  0.0100],
    [ 0.261799388,  0.0220],
    [ 0.349065850,  0.0440],
    [ 0.401425728,  0.0600],
    [ 0.453785606,  0.0840],
    [ 0.523598776,  0.1100],
])

# Spin-induced rolling-moment coefficient CRr as a function of angle of attack
# and advance ratio AdvR = d*omega_spin/(2*V). Rows are advance ratios, columns
# angles of attack.
_PC_CRR_ALPHA = np.array([
    -0.0873, -0.0698, -0.0524, -0.0349, -0.0175, 0.0000, 0.0175, 0.0349,
     0.0524,  0.0698,  0.0873,  0.1047,  0.1222, 0.1396, 0.1571, 0.1745,
     0.1920,  0.2094,  0.2269,  0.2443,  0.2618, 0.5236,
])

_PC_CRR_ADVR = np.array([2.0, 1.04, 0.69, 0.35, 0.17, 0.0])

_PC_CRR = np.array([
    [-0.0172, -0.0192, -0.0180, -0.0192, -0.0180, -0.0172, -0.0172, -0.0168,
     -0.0188, -0.0164, -0.0136, -0.0100, -0.0104, -0.0108, -0.0084, -0.0080,
     -0.0080, -0.0060, -0.0048, -0.0064, -0.0080, -0.0030],
    [-0.0112, -0.0132, -0.0120, -0.0132, -0.0120, -0.0112, -0.0112, -0.0108,
     -0.0128, -0.0104, -0.0096, -0.0068, -0.0072, -0.0076, -0.0052, -0.0048,
     -0.0048, -0.0028, -0.0032, -0.0048, -0.0064, -0.0030],
    [-0.0056, -0.0064, -0.0064, -0.0068, -0.0064, -0.0064, -0.0052, -0.0064,
     -0.0028, -0.0028, -0.0040, -0.0020, -0.0040, -0.0020, -0.0016,  0.0000,
      0.0000,  0.0000,  0.0000, -0.0020, -0.0048, -0.0030],
    [-0.0012, -0.0016, -0.0004, -0.0028, -0.0016, -0.0016, -0.0004,  0.0004,
      0.0004,  0.0008,  0.0004,  0.0008,  0.0012,  0.0008,  0.0020,  0.0028,
      0.0032,  0.0024,  0.0028,  0.0004, -0.0012, -0.0030],
    [-0.0012, -0.0012, -0.0016, -0.0016, -0.0012, -0.0004,  0.0004,  0.0008,
      0.0008,  0.0016,  0.0004,  0.0020,  0.0004,  0.0016,  0.0020,  0.0020,
      0.0020,  0.0012,  0.0012,  0.0000, -0.0012, -0.0030],
    [-0.0012, -0.0012, -0.0004, -0.0008, -0.0008, -0.0008,  0.0004,  0.0004,
      0.0004,  0.0008,  0.0004,  0.0008, -0.0004,  0.0000,  0.0000,  0.0004,
      0.0000,  0.0000,  0.0004, -0.0020, -0.0012, -0.0030],
])


class TabulatedCoefficients:
    """Measured Potts & Crowther (2002) curves, drop-in for DiscCoefficients.

    Interpolation is cubic (matching MATLAB's 'spline') in angle of attack and
    linear in advance ratio, with clamping to the tabulated range so that a
    transient excursion cannot produce a wild extrapolated coefficient. The
    clamping is recorded so validation can report how often it fired.

    The rate-damping coefficients CMq and CRp are NOT part of the P&C wind-tunnel
    data (a tunnel cannot measure a moment about a spinning axis, Hummel Sec.
    3.4) and are supplied from a linear set, by default the long-flight fit.
    """

    def __init__(self, damping: DiscCoefficients = None):
        from scipy.interpolate import CubicSpline, RegularGridInterpolator

        damping = damping if damping is not None else HUMMEL_LONG
        self.name = "Potts & Crowther 2002 (tabulated)"
        self.source = ("Potts & Crowther (2002) AIAA 2002-3150, curves as "
                       "tabulated in Hummel (2003) Appendix A. Damping terms "
                       f"CMq/CRp/CNr from: {damping.name}.")
        self.alpha_valid_deg = (-10.0, 30.0)

        self.CMq = damping.CMq
        self.CRp = damping.CRp
        self.CNr = damping.CNr

        self._alpha_lo = float(_PC_CL[0, 0])
        self._alpha_hi = float(_PC_CL[-1, 0])
        self._cl = CubicSpline(_PC_CL[:, 0], _PC_CL[:, 1])
        self._cd = CubicSpline(_PC_CD[:, 0], _PC_CD[:, 1])
        self._cm = CubicSpline(_PC_CM[:, 0], _PC_CM[:, 1])
        # RegularGridInterpolator needs ascending axes; AdvR is descending.
        self._crr = RegularGridInterpolator(
            (_PC_CRR_ADVR[::-1], _PC_CRR_ALPHA),
            _PC_CRR[::-1, :],
            method="linear", bounds_error=False, fill_value=None,
        )
        self.n_clamped = 0

    def _clamp(self, alpha: float) -> float:
        if alpha < self._alpha_lo:
            self.n_clamped += 1
            return self._alpha_lo
        if alpha > self._alpha_hi:
            self.n_clamped += 1
            return self._alpha_hi
        return alpha

    @property
    def alpha0(self) -> float:
        """Zero-lift angle from the measured CL curve (tabulated at -3 deg)."""
        return -0.05236

    def CL(self, alpha: float) -> float:
        return float(self._cl(self._clamp(alpha)))

    def CD(self, alpha: float) -> float:
        return float(self._cd(self._clamp(alpha)))

    def CM(self, alpha: float, q: float) -> float:
        return float(self._cm(self._clamp(alpha))) + self.CMq * q

    def CR(self, p: float, r: float, d: float, g: float = G_HUMMEL,
           V: float = None) -> float:
        """Rolling moment coefficient.

        With `V` supplied the tabulated CRr(alpha, AdvR) surface would be used;
        the trajectory integrator instead calls `CR_tabulated` directly because
        it already holds alpha and V. This signature exists so the class stays
        interface-compatible with `DiscCoefficients`.
        """
        return sqrt(d / g) * 0.00171 * r + self.CRp * p

    def CR_tabulated(self, alpha: float, p: float, r: float, d: float,
                     V: float) -> float:
        """CRr interpolated on (alpha, advance ratio), plus roll damping."""
        advr = abs(d * r / (2.0 * V)) if V > 1e-6 else 0.0
        advr = min(advr, float(_PC_CRR_ADVR[0]))
        crr = float(self._crr((advr, self._clamp(alpha))))
        # Table is for positive spin; the moment reverses with spin direction.
        return np.sign(r) * crr + self.CRp * p

    def spin_down_moment(self, r: float) -> float:
        return self.CNr * r

    def alpha_out_of_range(self, alpha_deg: float) -> bool:
        lo, hi = self.alpha_valid_deg
        return alpha_deg < lo or alpha_deg > hi


# The tabulated Potts & Crowther curves are the DEFAULT because they are the
# only set that reproduces Hummel's own published flights: they recover the
# Sec. 2.4.2 trim glide to three significant figures (alpha_trim 9.00 vs 9 deg,
# gamma -19.3 vs -19.3 deg, CL 0.672 vs 0.67) and the Sec. 2.4.3 25 m flight to
# within 8% in range, 4% in apex height and 1% in flight time. The linear sets
# extrapolate badly: a realistic throw reaches alpha > 30 deg late in the
# flight, well outside the +/-15 deg band where the linear fits were made.
# See validation_hummel.py for the executable comparison.
DEFAULT_COEFFS = TabulatedCoefficients()

COEFFICIENT_SETS = {
    "hummel_short": HUMMEL_SHORT,
    "hummel_long": HUMMEL_LONG,
    "potts_crowther": POTTS_CROWTHER,
}


def get_coefficients(name: str):
    """Look up a coefficient set by name, including the tabulated option."""
    if name in ("potts_tabulated", "tabulated"):
        return TabulatedCoefficients()
    if name not in COEFFICIENT_SETS:
        raise KeyError(
            f"unknown coefficient set {name!r}; available: "
            f"{sorted(list(COEFFICIENT_SETS) + ['potts_tabulated'])}"
        )
    return COEFFICIENT_SETS[name]


# ---------------------------------------------------------------------------
# Geometry / mass
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class DiscGeometry:
    """Discraft Ultrastar 175 g, using Hummel (2003) Table 3-1 MEASURED values.

    The moments of inertia matter: a thin uniform disc of this mass and diameter
    would have Is = 0.5*m*R^2 = 1.59e-3 and Ip = 0.25*m*R^2 = 7.97e-4 kg*m^2,
    which are 32% and 35% BELOW the measured 2.35e-3 and 1.22e-3. A real disc
    concentrates mass in the rim. Since the spin-axis inertia Is sets the
    gyroscopic precession rate (Omega = M / (Is*omega_z)), the uniform-disc
    idealisation overestimates turn and fade by roughly 50%. Use the measured
    values.
    """

    mass: float = 0.175                # kg      (Table 3-1)
    diameter: float = 0.269            # m       (Table 3-1)
    area: float = 0.057                # m^2     (Table 3-1, measured planform)
    Is: float = 0.00235                # kg*m^2  (Table 3-1, I_zz, spin axis)
    Ip: float = 0.00122                # kg*m^2  (Table 3-1, I_xx = I_yy)
    rho_air: float = 1.23              # kg/m^3  (Table 3-1)
    g: float = G_HUMMEL                # m/s^2   (Table 3-1)

    @property
    def radius(self) -> float:
        return self.diameter / 2.0

    @property
    def area_geometric(self) -> float:
        """pi*R^2 for reference; differs from the measured `area` by ~0.3%."""
        return pi * self.radius ** 2


ULTRASTAR = DiscGeometry()
