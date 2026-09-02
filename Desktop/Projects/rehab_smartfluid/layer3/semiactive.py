"""Layer 3: semi-active control of a strain-clocked smart-fluid element.

The plant is layer 2 wrapped in a device geometry. The control input is the
onset stress ``sigma*``, which is what a magnetic field actually moves in an MR
fluid: more field means frictional contacts form at lower stress. Modelling the
field this way rather than as a directly commanded force is the whole point of
what follows.

Three properties of this plant, none of which a Bouc-Wen damper model has
--------------------------------------------------------------------------
**1. The actuator has no instantaneous authority.** At any instant the force is
``F = eta(f) * (v/h) * S``. The control appears nowhere in that expression. It
enters only through ``df/dt``, so the relative degree from field to force is one
and every command is served through a lag. Clipped-optimal and skyhook laws are
written for an actuator whose damping coefficient switches immediately; here it
cannot.

**2. The lag is a strain, so the bandwidth depends on speed.** Linearising about
a steady operating point,

    d(delta f)/dt = -(gamma_dot / gamma_0) delta f + (gamma_dot / gamma_0) k u

which is first order with time constant ``tau_c = gamma_0 / gamma_dot`` and
bandwidth ``omega_c = gamma_dot / gamma_0``. The actuator gets faster the harder
you shear it.

**3. At rest the actuator is dead.** ``gamma_dot -> 0`` sends the bandwidth to
zero: with no relative motion the microstructure cannot be driven anywhere at
all, whatever the field. A controller cannot pre-position this device before an
impact. It has to earn its state during the event.

The unification worth noticing
------------------------------
The number of actuator time constants available during an event of duration
``tau`` is ``tau / tau_c = gamma_dot tau / gamma_0``, which is the same
activation number ``A`` that governs the passive shortfall in layer 2.

**A is simultaneously the strain budget of the material and the control-authority
budget of the actuator.** When ``A < 1`` the passive material has not finished
thickening and the active controller has not completed one response, for the same
reason and by the same arithmetic. Layers 2 and 3 fail together.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "layer2"))

from activation import STFMaterial, jamming_fraction  # noqa: E402

__all__ = [
    "SemiActiveElement",
    "ControlResult",
    "bandwidth",
    "time_constant",
    "simulate",
    "quasi_steady_command",
]


@dataclass(frozen=True)
class SemiActiveElement:
    """A shear-cell element: gap ``h``, shear area ``area``, controllable onset stress.

    ``sigma_star_lo`` is the strong-field limit (contacts form easily, most
    thickening) and ``sigma_star_hi`` the zero-field limit. The material's own
    ``sigma_star`` is ignored; the control supplies it.
    """

    material: STFMaterial
    h: float
    area: float
    sigma_star_lo: float
    sigma_star_hi: float

    def __post_init__(self) -> None:
        if self.h <= 0 or self.area <= 0:
            raise ValueError("gap and area must be positive")
        if not 0.0 < self.sigma_star_lo < self.sigma_star_hi:
            raise ValueError("need 0 < sigma_star_lo < sigma_star_hi")

    def clamp(self, sigma_star: float | np.ndarray) -> np.ndarray:
        """The actuator saturation. Commands outside the field range are unachievable."""
        return np.clip(sigma_star, self.sigma_star_lo, self.sigma_star_hi)

    def shear_rate(self, velocity: float | np.ndarray) -> np.ndarray:
        return np.abs(np.asarray(velocity, dtype=float)) / self.h

    def force(self, f: float | np.ndarray, velocity: float | np.ndarray) -> np.ndarray:
        """Force at frictional fraction ``f`` and relative velocity. Sign opposes motion."""
        eta = self.material.viscosity(f)
        v = np.asarray(velocity, dtype=float)
        return -np.sign(v) * eta * self.shear_rate(v) * self.area

    def f_target(self, f: float, velocity: float, sigma_star: float) -> float:
        """Steady frictional fraction the current stress and field would drive towards."""
        eta = float(self.material.viscosity(f))
        if not np.isfinite(eta):
            return 1.0
        sigma = eta * float(self.shear_rate(velocity))
        with np.errstate(divide="ignore", over="ignore"):
            return float(np.exp(-((sigma_star / max(sigma, 1e-300)) ** self.material.beta)))


def time_constant(velocity: float, element: SemiActiveElement) -> float:
    """tau_c = gamma_0 / gamma_dot. Infinite at rest: the actuator cannot act."""
    gd = float(element.shear_rate(velocity))
    return float("inf") if gd == 0.0 else element.material.gamma_0 / gd


def bandwidth(velocity: float, element: SemiActiveElement) -> float:
    """omega_c = gamma_dot / gamma_0, rad/s. Zero at rest."""
    return float(element.shear_rate(velocity)) / element.material.gamma_0


def quasi_steady_command(
    f_desired: float, f: float, velocity: float, element: SemiActiveElement
) -> float:
    """The field a quasi-steady controller asks for: whatever makes f_ss = f_desired now.

    Inverting f_ss = exp(-(sigma*/sigma)^beta) for sigma*. This is the command a
    designer writes when they believe the material tracks its flow curve, and it
    is what the study in scripts/run_control_study.py holds up against reality.
    """
    eta = float(element.material.viscosity(f))
    sigma = eta * float(element.shear_rate(velocity))
    fd = float(np.clip(f_desired, 1e-12, 1.0 - 1e-12))
    if sigma <= 0.0:
        return element.sigma_star_hi          # no motion: ask for no thickening
    return float(element.clamp(sigma * (-np.log(fd)) ** (1.0 / element.material.beta)))


@dataclass(frozen=True)
class ControlResult:
    t: np.ndarray
    f: np.ndarray
    force: np.ndarray
    force_desired: np.ndarray
    sigma_star: np.ndarray
    velocity: np.ndarray
    saturated: np.ndarray

    @property
    def tracking_rms(self) -> float:
        """RMS force error as a fraction of the peak demanded force."""
        peak = np.max(np.abs(self.force_desired))
        if peak <= 0:
            return float("nan")
        return float(np.sqrt(np.mean((self.force - self.force_desired) ** 2)) / peak)

    @property
    def saturated_fraction(self) -> float:
        """Fraction of the event during which the command was outside the field range."""
        return float(np.mean(self.saturated))


def simulate(
    element: SemiActiveElement,
    velocity_of_t: Callable[[float], float],
    force_desired_of_t: Callable[[float], float],
    controller: Callable[[float, float, float, SemiActiveElement], float],
    t_end: float,
    *,
    f0: float = 0.0,
    n_out: int = 1001,
) -> ControlResult:
    """Closed-loop simulation of the element under a semi-active control law.

    ``controller(f_desired, f, velocity, element) -> sigma_star``. The desired
    frictional fraction is recovered from the demanded force by inverting the
    viscosity, which is exact and instantaneous; the point of the study is that
    reaching it is neither.
    """
    mat = element.material

    def f_for_force(F_d: float, velocity: float) -> float:
        """Frictional fraction that would deliver F_d right now, or 1 if unreachable."""
        gd = float(element.shear_rate(velocity))
        if gd <= 0.0:
            return 0.0
        eta_needed = abs(F_d) / (gd * element.area)
        if eta_needed <= mat.eta_s:
            return 0.0
        # eta = eta_s (1 - phi/phi_J)^-n  ->  phi_J, then invert phi_J(f).
        phi_j = mat.phi / (1.0 - (eta_needed / mat.eta_s) ** (-1.0 / mat.n))
        f = (mat.phi_0 - phi_j) / (mat.phi_0 - mat.phi_m)
        return float(np.clip(f, 0.0, 1.0))

    def rhs(t: float, y: np.ndarray) -> list[float]:
        f = float(np.clip(y[0], 0.0, 1.0))
        v = float(velocity_of_t(t))
        gd = float(element.shear_rate(v))
        if gd == 0.0:
            return [0.0]                        # dead actuator, property 3
        u = controller(f_for_force(force_desired_of_t(t), v), f, v, element)
        return [gd * (element.f_target(f, v, u) - f) / mat.gamma_0]

    t_eval = np.linspace(0.0, t_end, n_out)
    sol = solve_ivp(rhs, (0.0, t_end), [f0], t_eval=t_eval,
                    rtol=1e-8, atol=1e-10, method="LSODA", max_step=t_end / 400.0)
    if not sol.success:
        raise RuntimeError(f"closed-loop integration failed: {sol.message}")

    t = sol.t
    f = np.clip(sol.y[0], 0.0, 1.0)
    v = np.array([float(velocity_of_t(ti)) for ti in t])
    Fd = np.array([float(force_desired_of_t(ti)) for ti in t])
    F = element.force(f, v)

    us, sat = [], []
    for fi, vi, Fdi in zip(f, v, Fd):
        raw_target = f_for_force(Fdi, vi)
        u_unclamped = quasi_steady_command(raw_target, float(fi), float(vi), element)
        u = controller(raw_target, float(fi), float(vi), element)
        us.append(u)
        sat.append(not np.isclose(u, u_unclamped, rtol=1e-9)
                   or np.isclose(u, element.sigma_star_lo, rtol=1e-9))
    assert t.shape == f.shape == F.shape, "closed-loop arrays must align"

    return ControlResult(t=t, f=f, force=np.abs(F), force_desired=np.abs(Fd),
                         sigma_star=np.array(us), velocity=v,
                         saturated=np.array(sat, dtype=bool))


def steady_f_at_field(element: SemiActiveElement, velocity: float,
                      sigma_star: float, *, iters: int = 200) -> float:
    """Fixed point of f = f_ss(eta(f) gamma_dot; sigma*), by damped iteration.

    Used as the analytic operating point the linearisation in the module
    docstring is taken about.
    """
    f = 0.0
    for _ in range(iters):
        nxt = element.f_target(f, velocity, sigma_star)
        if not np.isfinite(nxt) or element.material.phi >= jamming_fraction(
            nxt, element.material.phi_0, element.material.phi_m
        ):
            return float("nan")
        f, prev = 0.5 * (f + nxt), f
        if abs(f - prev) < 1e-14:
            break
    return f
