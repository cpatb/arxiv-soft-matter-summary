"""
Six-degree-of-freedom rigid-body integrator for a flying disc.

State vector (13 components):
    r       : position in inertial frame      (3)  [m]
    v       : velocity in inertial frame      (3)  [m/s]
    q       : quaternion B -> N, (qw, qx, qy, qz) (4) [-]
              q rotates a body-frame vector into the inertial frame.
    omega_B : body-frame angular velocity     (3)  [rad/s]
              omega_B = (omega_x, omega_y, omega_z); omega_z is the spin rate.

Body frame B is material (rotates with the disc); b3_B is the disc normal
(spin axis), pointing through the convex "top" of the disc.

Wind field U(r, t) is supplied as a callable. Defaults to still air.

Reference: Hummel (2003), Crowther & Potts (2007), Hubbard & Hummel (2000).
See Project1_AnalyticalTreatment.tex for the derivation of the EOMs.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Callable, Optional
import numpy as np
from scipy.integrate import solve_ivp

from coefficients import (
    DEFAULT_COEFFS, DiscCoefficients, DiscGeometry, HUMMEL_LONG, HUMMEL_SHORT,
    ULTRASTAR,
)


# ---------------------------------------------------------------------------
# Quaternion utilities
# ---------------------------------------------------------------------------

def q_normalize(q: np.ndarray) -> np.ndarray:
    return q / np.linalg.norm(q)


def q_to_R(q: np.ndarray) -> np.ndarray:
    """Convert quaternion (qw, qx, qy, qz) to a rotation matrix R s.t.
    v_N = R @ v_B."""
    qw, qx, qy, qz = q
    R = np.array([
        [1 - 2*(qy*qy + qz*qz), 2*(qx*qy - qz*qw),     2*(qx*qz + qy*qw)],
        [2*(qx*qy + qz*qw),     1 - 2*(qx*qx + qz*qz), 2*(qy*qz - qx*qw)],
        [2*(qx*qz - qy*qw),     2*(qy*qz + qx*qw),     1 - 2*(qx*qx + qy*qy)],
    ])
    return R


def q_kinematics(q: np.ndarray, omega_B: np.ndarray) -> np.ndarray:
    """dq/dt = 0.5 * q * (0, omega_B)."""
    qw, qx, qy, qz = q
    wx, wy, wz = omega_B
    dq = 0.5 * np.array([
        -qx*wx - qy*wy - qz*wz,
         qw*wx + qy*wz - qz*wy,
         qw*wy - qx*wz + qz*wx,
         qw*wz + qx*wy - qy*wx,
    ])
    return dq


def axis_angle_to_q(axis: np.ndarray, angle: float) -> np.ndarray:
    a = axis / np.linalg.norm(axis)
    s = np.sin(angle / 2.0)
    return np.array([np.cos(angle / 2.0), s * a[0], s * a[1], s * a[2]])


# ---------------------------------------------------------------------------
# Aerodynamic force and moment
# ---------------------------------------------------------------------------

def aero_force_moment(
    v_N: np.ndarray,
    R_BN: np.ndarray,
    omega_B: np.ndarray,
    U_N: np.ndarray,
    coeffs: DiscCoefficients,
    geom: DiscGeometry,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """
    Compute aerodynamic force in inertial frame and moment in body frame.

    Returns
    -------
    F_N      : aerodynamic force in inertial frame, shape (3,)
    M_B      : aerodynamic moment in body frame,    shape (3,)
    diag     : dict with intermediate quantities (alpha, V, V_par, CL, CD, ...)
    """
    # Relative wind in inertial frame, then transformed to body
    V_rel_N = v_N - U_N
    V = np.linalg.norm(V_rel_N)
    if V < 1e-4:
        return np.zeros(3), np.zeros(3), {"V": V, "alpha": 0.0}

    # Body-frame components of relative wind
    V_rel_B = R_BN.T @ V_rel_N
    Vx_B, Vy_B, Vz_B = V_rel_B
    V_par = float(np.hypot(Vx_B, Vy_B))

    # Angle of attack: alpha > 0 when wind strikes underside (Vz_B < 0)
    alpha = np.arctan2(-Vz_B, V_par) if V_par > 1e-6 else (
        np.pi / 2.0 * np.sign(-Vz_B)
    )

    # Aerodynamic-frame unit vectors expressed in body frame
    if V_par > 1e-6:
        b1_aero = np.array([Vx_B / V_par, Vy_B / V_par, 0.0])
        b2_aero = np.array([-Vy_B / V_par, Vx_B / V_par, 0.0])
    else:
        # Edge case: wind purely along spin axis; pick arbitrary in-plane axes
        b1_aero = np.array([1.0, 0.0, 0.0])
        b2_aero = np.array([0.0, 1.0, 0.0])

    # Angular rates resolved in aero frame: p, q, r
    omega_aero = np.array([
        omega_B @ b1_aero,   # p, about V_par direction (roll)
        omega_B @ b2_aero,   # q, about in-plane perpendicular (pitch)
        omega_B[2],          # r, about disc normal (spin)
    ])
    p_a, q_a, r_a = omega_aero

    # ---- Frame reconciliation with Hummel (2003) --------------------------
    # The coefficient sets were fitted in Hummel's inertial frame, which has z
    # pointing DOWN and disc normal c3 through the disc's UNDERSIDE. We use z UP
    # with b3 through the disc's TOP, so
    #     c3   = -b3
    #     ulat = c3 x uvp = -b2_aero      (pitch axis is anti-parallel)
    #     uvp  = +b1_aero                 (roll axis is parallel)
    # Angular rates therefore transform as below. Both the rate and the axis
    # flip for pitch, so the DAMPING term keeps its sign while the STATIC
    # pitching moment reverses; getting only one of the two right silently
    # inverts the disc's turn-and-fade behaviour.
    omega_spin_H = -r_a          # spin rate about c3
    omega_lat_H = -q_a           # pitch rate about ulat
    omega_vp_H = +p_a            # roll rate about uvp

    # Aerodynamic coefficients, evaluated in Hummel's convention
    CL = coeffs.CL(alpha)
    CD = coeffs.CD(alpha)
    CM_H = coeffs.CM(alpha, omega_lat_H)
    if hasattr(coeffs, "CR_tabulated"):
        CR_H = coeffs.CR_tabulated(alpha, omega_vp_H, omega_spin_H,
                                   geom.diameter, V)
    else:
        CR_H = coeffs.CR(omega_vp_H, omega_spin_H, geom.diameter, geom.g)
    # Spin-down moment comes back as a moment in N*m, not a coefficient, and is
    # NOT scaled by the dynamic-pressure factor (Hummel `discfltEOM.m`).
    M_spin_H = coeffs.spin_down_moment(omega_spin_H)

    qbar = 0.5 * geom.rho_air * V * V
    qbarA = qbar * geom.area
    qbarAd = qbarA * geom.diameter

    # Aero-frame unit vectors of lift and drag, expressed in body frame
    v_hat_B = V_rel_B / V
    # b3_B = (0,0,1) in body frame; project out v_hat component
    b3_B = np.array([0.0, 0.0, 1.0])
    L_dir_B = b3_B - (b3_B @ v_hat_B) * v_hat_B
    L_dir_norm = np.linalg.norm(L_dir_B)
    if L_dir_norm > 1e-6:
        L_hat_B = L_dir_B / L_dir_norm
    else:
        L_hat_B = np.zeros(3)

    F_B = qbarA * (CL * L_hat_B - CD * v_hat_B)
    F_N = R_BN @ F_B

    # Moments, converted back into our z-up body basis:
    #   roll  about uvp  = +b1_aero  -> coefficient unchanged
    #   pitch about ulat = -b2_aero  -> coefficient negated
    #   spin  about c3   = -b3       -> moment negated
    M_B = (qbarAd * (CR_H * b1_aero + (-CM_H) * b2_aero)
           + (-M_spin_H) * b3_B)

    diag = dict(
        V=V, V_par=V_par, Vz_B=Vz_B,
        alpha=alpha, alpha_deg=np.degrees(alpha),
        CL=CL, CD=CD, CM=CM_H, CR=CR_H, M_spin=M_spin_H,
        p_a=p_a, q_a=q_a, r_a=r_a,
        qbarA=qbarA,
    )
    return F_N, M_B, diag


# ---------------------------------------------------------------------------
# Equations of motion
# ---------------------------------------------------------------------------

def rhs(
    t: float,
    y: np.ndarray,
    coeffs: DiscCoefficients,
    geom: DiscGeometry,
    wind: Callable[[np.ndarray, float], np.ndarray],
    g: float | None = None,
) -> np.ndarray:
    """Right-hand side of the 6-DOF system. y is the 13-element state."""
    # Default to the gravity used when the coefficients were fitted (9.794),
    # so the model is self-consistent with its own calibration.
    if g is None:
        g = geom.g
    r = y[0:3]
    v = y[3:6]
    q = y[6:10]
    omega_B = y[10:13]

    # Renormalize quaternion (light, no integrator-side projection needed)
    q = q_normalize(q)
    R_BN = q_to_R(q)

    U_N = wind(r, t)
    F_N, M_B, _ = aero_force_moment(v, R_BN, omega_B, U_N, coeffs, geom)

    # Translation
    dr = v
    dv = F_N / geom.mass + np.array([0.0, 0.0, -g])

    # Orientation
    dq = q_kinematics(q, omega_B)

    # Rotation (Euler's equations for axisymmetric body):
    # Ip dwx + (Is - Ip) wy wz = Mx
    # Ip dwy - (Is - Ip) wx wz = My
    # Is dwz = Mz
    wx, wy, wz = omega_B
    Is, Ip = geom.Is, geom.Ip
    Mx, My, Mz = M_B
    domega = np.array([
        (Mx - (Is - Ip) * wy * wz) / Ip,
        (My + (Is - Ip) * wx * wz) / Ip,
        Mz / Is,
    ])

    dy = np.empty(13)
    dy[0:3]   = dr
    dy[3:6]   = dv
    dy[6:10]  = dq
    dy[10:13] = domega
    return dy


# ---------------------------------------------------------------------------
# Convenience wrappers
# ---------------------------------------------------------------------------

def still_air(r: np.ndarray, t: float) -> np.ndarray:
    return np.zeros(3)


def uniform_wind(U_vec: np.ndarray):
    U_vec = np.asarray(U_vec, dtype=float)
    def _f(r, t):
        return U_vec
    return _f


@dataclass
class ThrowResult:
    t: np.ndarray
    r: np.ndarray            # (N, 3)
    v: np.ndarray            # (N, 3)
    q: np.ndarray            # (N, 4)
    omega_B: np.ndarray      # (N, 3)
    diag: list[dict]         # per-timestep aero diagnostics
    range_m: float
    flight_time: float
    landing_idx: int


def simulate(
    r0: np.ndarray,
    v0: np.ndarray,
    q0: np.ndarray,
    omega_B0: np.ndarray,
    coeffs: DiscCoefficients = DEFAULT_COEFFS,
    geom: DiscGeometry = ULTRASTAR,
    wind: Callable[[np.ndarray, float], np.ndarray] = still_air,
    t_max: float = 8.0,
    rtol: float = 1e-6,
    atol: float = 1e-8,
    z_ground: float = 0.0,
) -> ThrowResult:
    """Integrate one throw, stopping when the disc hits the ground."""
    def hit_ground(t, y):
        return y[2] - z_ground
    hit_ground.terminal = True
    hit_ground.direction = -1

    y0 = np.concatenate([r0, v0, q_normalize(q0), omega_B0])
    sol = solve_ivp(
        lambda t, y: rhs(t, y, coeffs, geom, wind),
        t_span=(0.0, t_max),
        y0=y0,
        method="LSODA",
        events=hit_ground,
        rtol=rtol, atol=atol,
        dense_output=False,
        max_step=0.05,
    )

    t = sol.t
    r = sol.y[0:3].T
    v = sol.y[3:6].T
    q = sol.y[6:10].T
    omega = sol.y[10:13].T

    # Reconstruct per-timestep diagnostics
    diag = []
    for k in range(len(t)):
        R_BN = q_to_R(q_normalize(q[k]))
        _, _, d = aero_force_moment(v[k], R_BN, omega[k],
                                     wind(r[k], t[k]), coeffs, geom)
        diag.append(d)

    landing_idx = len(t) - 1
    landing_xy = r[landing_idx, 0:2]
    range_m = float(np.linalg.norm(landing_xy - r0[0:2]))

    return ThrowResult(
        t=t, r=r, v=v, q=q, omega_B=omega, diag=diag,
        range_m=range_m, flight_time=float(t[-1]), landing_idx=landing_idx,
    )


# ---------------------------------------------------------------------------
# Initial-condition helper for canonical throws
# ---------------------------------------------------------------------------

def build_initial_conditions(
    launch_pos: np.ndarray,
    speed: float,
    heading_deg: float,
    elevation_deg: float,
    release_alpha_deg: float,
    roll_deg: float = 0.0,
    spin_rad_s: float = 100.0,
    handedness: str = "right_backhand",
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build (r0, v0, q0, omega_B0) for a throw with intuitive parameters.

    Parameters
    ----------
    launch_pos : (3,) array, m
    speed      : release speed, m/s
    heading_deg: azimuth of release velocity, measured from +x axis
                 counterclockwise (i.e. standard math convention).
    elevation_deg: pitch of release velocity above horizontal
    release_alpha_deg: angle between disc plane and ground at release
                      (NOT angle of attack — this is the body orientation;
                      α arises from the velocity relative to the disc).
                      For a "flat" release this is 0; for hyzer-style
                      (right-handed backhand standard), negative.
    roll_deg   : bank angle of the disc at release (hyzer/anhyzer);
                 positive = right-hand-side-down (anhyzer for RH-BH)
    spin_rad_s : magnitude of spin rate (positive number)
    handedness : "right_backhand" or "right_forehand" (or "left_*")
                 sets the sign of spin about the disc normal.
                 Right backhand: disc spins clockwise viewed from above
                                 (angular velocity vector points DOWN
                                 through the disc), so omega_z (body) < 0
                                 if body z points up through disc top.
    """
    r0 = np.asarray(launch_pos, dtype=float).copy()

    az = np.radians(heading_deg)
    el = np.radians(elevation_deg)
    v_dir = np.array([
        np.cos(el) * np.cos(az),
        np.cos(el) * np.sin(az),
        np.sin(el),
    ])
    v0 = speed * v_dir

    # Disc orientation: start with body frame aligned with inertial,
    # then rotate to apply pitch (release_alpha_deg about y-body), then
    # bank (roll_deg about x-body, after the pitch), then yaw to align
    # heading. We construct a rotation matrix R_BN and convert to quat.

    # Build orientation matrix from sequence: yaw -> pitch -> roll
    # so disc x_B initially points in the throw heading direction.
    # Convention: release_alpha_deg = angle of attack at release.
    # Geometry: disc-plane pitch above horizontal = elevation + alpha
    # (when both alpha and elevation are positive: disc is more nose-up
    # than the velocity by exactly alpha).
    # In our z-up frame, my Ry takes body z to (sin theta, 0, cos theta) for
    # positive theta. A NOSE-UP disc has body z tilted BACKWARD (opposite
    # to heading), so we pass a NEGATIVE pitch to Ry to obtain a nose-up disc.
    yaw = az
    pitch = -np.radians(elevation_deg + release_alpha_deg)
    roll = np.radians(roll_deg)

    cy, sy = np.cos(yaw), np.sin(yaw)
    cp, sp = np.cos(pitch), np.sin(pitch)
    cr, sr = np.cos(roll), np.sin(roll)

    # ZYX intrinsic rotations (yaw, then pitch, then roll)
    Rz = np.array([[cy, -sy, 0], [sy, cy, 0], [0, 0, 1]])
    Ry = np.array([[cp, 0, sp], [0, 1, 0], [-sp, 0, cp]])
    Rx = np.array([[1, 0, 0], [0, cr, -sr], [0, sr, cr]])
    R_BN = Rz @ Ry @ Rx

    # Convert R_BN to quaternion (qw, qx, qy, qz)
    tr = R_BN[0, 0] + R_BN[1, 1] + R_BN[2, 2]
    if tr > 0:
        S = 2.0 * np.sqrt(1.0 + tr)
        qw = 0.25 * S
        qx = (R_BN[2, 1] - R_BN[1, 2]) / S
        qy = (R_BN[0, 2] - R_BN[2, 0]) / S
        qz = (R_BN[1, 0] - R_BN[0, 1]) / S
    elif R_BN[0, 0] > R_BN[1, 1] and R_BN[0, 0] > R_BN[2, 2]:
        S = 2.0 * np.sqrt(1.0 + R_BN[0, 0] - R_BN[1, 1] - R_BN[2, 2])
        qw = (R_BN[2, 1] - R_BN[1, 2]) / S
        qx = 0.25 * S
        qy = (R_BN[0, 1] + R_BN[1, 0]) / S
        qz = (R_BN[0, 2] + R_BN[2, 0]) / S
    elif R_BN[1, 1] > R_BN[2, 2]:
        S = 2.0 * np.sqrt(1.0 + R_BN[1, 1] - R_BN[0, 0] - R_BN[2, 2])
        qw = (R_BN[0, 2] - R_BN[2, 0]) / S
        qx = (R_BN[0, 1] + R_BN[1, 0]) / S
        qy = 0.25 * S
        qz = (R_BN[1, 2] + R_BN[2, 1]) / S
    else:
        S = 2.0 * np.sqrt(1.0 + R_BN[2, 2] - R_BN[0, 0] - R_BN[1, 1])
        qw = (R_BN[1, 0] - R_BN[0, 1]) / S
        qx = (R_BN[0, 2] + R_BN[2, 0]) / S
        qy = (R_BN[1, 2] + R_BN[2, 1]) / S
        qz = 0.25 * S
    q0 = np.array([qw, qx, qy, qz])

    # Spin direction
    sign = {
        "right_backhand": -1.0,   # spin angular velocity points DOWN
        "left_forehand":  -1.0,
        "right_forehand": +1.0,   # spin angular velocity points UP
        "left_backhand":  +1.0,
    }.get(handedness, -1.0)
    omega_B0 = np.array([0.0, 0.0, sign * spin_rad_s])

    return r0, v0, q0, omega_B0
