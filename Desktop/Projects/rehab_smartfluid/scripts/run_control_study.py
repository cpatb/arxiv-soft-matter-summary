"""What a quasi-steady controller actually delivers on a strain-clocked element.

Two panels.

  (a) Tracking error of a quasi-steady semi-active controller against the
      activation number A. The controller asks, at every instant, for the field
      that would produce the demanded force IF the material tracked its flow
      curve. It does not, and the error is measured rather than argued.

  (b) The actuator bandwidth map. omega_c = gamma_dot / gamma_0, converted to Hz,
      over the published gait shear-rate band and the measured strain-clock
      bracket, against the bandwidth a gait event demands. This panel is closed
      form; it is drawn rather than simulated.

Run from the project root:

    python scripts/run_control_study.py
"""

from __future__ import annotations

import hashlib
import sys
from datetime import date
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "layer2"))
sys.path.insert(0, str(ROOT / "layer3"))

from activation import STFMaterial  # noqa: E402
from semiactive import (  # noqa: E402
    SemiActiveElement,
    quasi_steady_command,
    simulate,
    steady_f_at_field,
)

FIGDIR = ROOT / "figures"

# Published critical band for the STF prosthetic foot (Smart Mater. Struct. 2020).
GAIT_BAND = (10.0, 80.0)
# Strain clock measured in dst_fit from Lin et al. (2015) reversal transients.
GAMMA_0_LO, GAMMA_0_HI = 0.5, 1.6
# Event bandwidths: weight acceptance is 100 ms (loading response, 2-12% of a
# 1.0 s cycle). The 30 ms stumble is a stand-in and is NOT sourced.
TAU_WEIGHT, TAU_STUMBLE = 0.100, 0.030

SIGMA_STAR, PHI_M, PHI_0, N_EXP = 50.0, 0.54, 0.64, 2.0
ETA_S = SIGMA_STAR / ((1.0 - 0.40 / PHI_0) ** -N_EXP * 100.0)
GAP, AREA = 1e-3, 1e-3

# Velocity grid for the tabulated steady demand. Numerical, not physical: main()
# checks that doubling it does not move the reported tracking error.
DEMAND_GRID = 241


def build(phi: float, gamma_0: float) -> SemiActiveElement:
    mat = STFMaterial(phi=phi, phi_m=PHI_M, sigma_star=SIGMA_STAR, eta_s=ETA_S,
                      gamma_0=gamma_0, phi_0=PHI_0, n=N_EXP)
    return SemiActiveElement(material=mat, h=GAP, area=AREA,
                             sigma_star_lo=5.0, sigma_star_hi=500.0)


def tracking_error(el: SemiActiveElement, gamma_dot_peak: float, tau: float,
                   *, demand_grid: int = DEMAND_GRID) -> float:
    """RMS force error of a quasi-steady controller over one half-sine event."""
    v_peak = gamma_dot_peak * el.h

    def velocity(t: float) -> float:
        return v_peak * np.sin(np.pi * min(max(t / tau, 0.0), 1.0))

    # The demand: what full field would deliver in STEADY state at this instant.
    # This is exactly the specification a quasi-steady designer writes down.
    #
    # Tabulated on a velocity grid and interpolated. The fixed point is solved by
    # iteration and the ODE solver asks for the demand thousands of times; the
    # demand depends on t only through v, so tabulating in v costs nothing but
    # interpolation error. main() checks that doubling the grid does not move the
    # reported number.
    v_grid = np.linspace(0.0, v_peak, demand_grid)
    F_grid = np.empty_like(v_grid)
    for i, vi in enumerate(v_grid):
        if vi <= 0.0:
            F_grid[i] = 0.0
            continue
        f_ss = steady_f_at_field(el, float(vi), el.sigma_star_lo)
        F_grid[i] = 0.0 if not np.isfinite(f_ss) else float(abs(el.force(f_ss, float(vi))))

    def demand(t: float) -> float:
        return float(np.interp(velocity(t), v_grid, F_grid))

    res = simulate(el, velocity, demand, quasi_steady_command, t_end=tau, n_out=801)
    return res.tracking_rms


def main() -> int:
    FIGDIR.mkdir(exist_ok=True)

    gamma_dot_peak = 100.0
    activations = np.geomspace(0.1, 60.0, 22)
    phis = (0.40, 0.46, 0.50)
    curves: dict[float, np.ndarray] = {}
    for phi in phis:
        el = build(phi, 1.0)
        curves[phi] = np.array([
            tracking_error(el, gamma_dot_peak, A * el.material.gamma_0 / gamma_dot_peak)
            for A in activations
        ])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.0), constrained_layout=True)

    for phi in phis:
        ax1.semilogx(activations, curves[phi], lw=1.4, label=f"{phi:.2f}")
    # Linear y from 0: the error is a fraction of peak demand.
    ax1.set_ylim(0.0, 1.0)
    ax1.set_xlabel(r"$A = \dot\gamma\,\tau/\gamma_0$")
    ax1.set_ylabel("force tracking error (RMS / peak)")
    ax1.legend(title=r"$\phi$", frameon=False, fontsize=7, title_fontsize=7)

    gd = np.geomspace(1.0, 1e3, 200)
    ax2.fill_between(gd, gd / GAMMA_0_HI / (2 * np.pi), gd / GAMMA_0_LO / (2 * np.pi),
                     color="0.85", zorder=0)
    ax2.loglog(gd, gd / GAMMA_0_LO / (2 * np.pi), lw=1.2, color="k")
    ax2.loglog(gd, gd / GAMMA_0_HI / (2 * np.pi), lw=1.2, color="k")
    ax2.axhline(1.0 / TAU_WEIGHT, color="C0", lw=1.0, ls="--")
    ax2.axhline(1.0 / TAU_STUMBLE, color="C3", lw=1.0, ls=":")
    ax2.axvspan(*GAIT_BAND, color="C2", alpha=0.15, zorder=0)
    # Both axes log: each spans three decades.
    ax2.set_xlabel(r"$\dot\gamma$  (s$^{-1}$)")
    ax2.set_ylabel(r"$\omega_c/2\pi$  (Hz)")

    stem = FIGDIR / "control_study_2026-09-01"
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".png"), dpi=200)
    plt.close(fig)

    with stem.with_suffix(".csv").open("w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"# generated {date.today().isoformat()} by scripts/run_control_study.py\n")
        fh.write(f"# gamma_dot_peak={gamma_dot_peak} 1/s, gap={GAP} m, area={AREA} m^2\n")
        fh.write("panel,x,series,y\n")
        for phi in phis:
            for A, e in zip(activations, curves[phi]):
                fh.write(f"a,{A:.6g},{phi:.3f},{e:.6g}\n")

    print("actuator bandwidth over the published gait band:")
    for gd_ in GAIT_BAND:
        for g0, tag in ((GAMMA_0_HI, "slow clock"), (GAMMA_0_LO, "fast clock")):
            print(f"  gamma_dot={gd_:5.1f} 1/s, gamma_0={g0:.1f} ({tag:10s}) -> "
                  f"{gd_/g0:7.2f} rad/s = {gd_/g0/(2*np.pi):6.2f} Hz")
    print(f"\nevent bandwidth demanded: weight acceptance {1/TAU_WEIGHT:.0f} Hz, "
          f"stumble {1/TAU_STUMBLE:.0f} Hz (stumble UNSOURCED)")

    # Grid adequacy: the interpolation must not be setting the answer.
    el_chk, tau_chk = build(0.46, 1.0), 1.0 / gamma_dot_peak
    coarse = tracking_error(el_chk, gamma_dot_peak, tau_chk, demand_grid=DEMAND_GRID)
    fine = tracking_error(el_chk, gamma_dot_peak, tau_chk, demand_grid=2 * DEMAND_GRID - 1)
    print(f"\ndemand-grid check at A=1, phi=0.46: {DEMAND_GRID} -> {coarse:.6f}, "
          f"{2*DEMAND_GRID-1} -> {fine:.6f}, shift {abs(fine-coarse):.2e}")

    print("\nquasi-steady controller tracking error (RMS / peak demand):")
    print("     A   " + "".join(f"{p:>10.2f}" for p in phis))
    for i, A in enumerate(activations):
        if i % 4 == 0 or i == len(activations) - 1:
            print(f"  {A:6.2f} " + "".join(f"{curves[p][i]:10.4f}" for p in phis))

    for path in (stem.with_suffix(".pdf"), stem.with_suffix(".png"), stem.with_suffix(".csv")):
        print(f"{hashlib.md5(path.read_bytes()).hexdigest()[:12]}  {path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
