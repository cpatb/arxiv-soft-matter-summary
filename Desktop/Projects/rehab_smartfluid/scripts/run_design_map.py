"""Phase-0 design map: what a strain-clocked STF element can and cannot deliver.

Produces two panels and the CSV behind them.

  (a) Quasi-steady shortfall. Peak stress reached during an event of duration
      tau at rate gamma_dot, divided by the steady stress the device literature
      would have designed to, against the activation number A = gamma_dot tau /
      gamma_0. The quasi-steady assumption is the right edge of this plot.

  (b) The reachable-rate ceiling. The largest steady shear rate a packing can
      sustain before it shear jams, against packing fraction, with the measured
      gait band marked. This is the constraint that a model without a phi term
      cannot express.

Run from the project root:

    python scripts/run_design_map.py

Every number in the accompanying SCALING note comes from this script. Nothing
here is fitted to anything: the material is placed at the thickening onset by
construction (see MATERIAL below) and the result is reported as a dimensionless
map, so no invented material constants enter.
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

from activation import (  # noqa: E402
    STFMaterial,
    quasi_steady_ratio,
    steady_curve,
)

FIGDIR = ROOT / "figures"

# Reference rate inside the 10-80 1/s critical band measured for the published
# STF prosthetic foot (Smart Mater. Struct. 2020, verified 2026-08-31).
GAIT_BAND = (10.0, 80.0)
GAMMA_DOT_REF = 100.0

# Carrier viscosity solved, not chosen, so the lubricated stress at GAMMA_DOT_REF
# equals sigma_star: eta_s = sigma_star / ((1 - phi/phi_0)^-n gamma_dot). The
# fixture therefore sits at the thickening onset instead of in a dead regime.
SIGMA_STAR = 50.0
PHI_REF, PHI_M, PHI_0, N_EXP = 0.40, 0.54, 0.64, 2.0
ETA_S = SIGMA_STAR / ((1.0 - PHI_REF / PHI_0) ** -N_EXP * GAMMA_DOT_REF)

# Strain clock bracket measured from Lin et al. (2015) shear-reversal transients,
# dst_fit/STRAIN_CLOCK_2026-08-28.md. NOT the 0.2 borrowed from corn starch.
GAMMA_0_BRACKET = (0.5, 1.6)

MATERIAL = STFMaterial(phi=PHI_REF, phi_m=PHI_M, sigma_star=SIGMA_STAR,
                       eta_s=ETA_S, gamma_0=1.0, phi_0=PHI_0, n=N_EXP)


def shortfall_curve(phi: float, activation: np.ndarray) -> np.ndarray:
    """Peak stress delivered / quasi-steady stress, against activation number."""
    mat = STFMaterial(phi=phi, phi_m=PHI_M, sigma_star=SIGMA_STAR,
                      eta_s=ETA_S, gamma_0=1.0, phi_0=PHI_0, n=N_EXP)
    out = np.empty_like(activation)
    for i, A in enumerate(activation):
        tau = A * mat.gamma_0 / GAMMA_DOT_REF
        out[i], _, _ = quasi_steady_ratio(GAMMA_DOT_REF, tau, mat)
    return out


def reachable_rate(phi: float) -> float:
    """Largest steady shear rate this packing can sustain before shear jamming."""
    mat = STFMaterial(phi=phi, phi_m=PHI_M, sigma_star=SIGMA_STAR,
                      eta_s=ETA_S, gamma_0=1.0, phi_0=PHI_0, n=N_EXP)
    sigma = np.geomspace(1e-2, 1e5, 200001)
    rate, eta = steady_curve(sigma, mat)
    ok = np.isfinite(eta) & (rate > 0)
    return float(rate[ok].max()) if ok.any() else float("nan")


def main() -> int:
    FIGDIR.mkdir(exist_ok=True)

    activation = np.geomspace(0.05, 100.0, 60)
    phis_short = (0.40, 0.46, 0.50, 0.52)
    shortfalls = {phi: shortfall_curve(phi, activation) for phi in phis_short}

    phis_rate = np.round(np.arange(0.40, 0.565, 0.005), 4)
    rates = np.array([reachable_rate(p) for p in phis_rate])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.0), constrained_layout=True)

    for phi in phis_short:
        ax1.semilogx(activation, shortfalls[phi], lw=1.4, label=f"{phi:.2f}")
    ax1.axhline(1.0, color="0.6", lw=0.8, ls=":")
    # Linear y from 0: the shortfall is a fraction and suppressing zero would
    # exaggerate it. x is log because A spans more than three decades.
    ax1.set_ylim(0.0, 1.05)
    ax1.set_xlabel(r"$A = \dot\gamma\,\tau/\gamma_0$")
    ax1.set_ylabel(r"$\sigma_{\rm peak}/\sigma_{\rm steady}$")
    ax1.legend(title=r"$\phi$", frameon=False, fontsize=7, title_fontsize=7)

    ax2.semilogy(phis_rate, rates, lw=1.4, color="k")
    ax2.axhspan(*GAIT_BAND, color="0.85", zorder=0)
    ax2.axvline(PHI_M, color="0.6", lw=0.8, ls=":")
    # y is log: the ceiling spans three decades over this phi range.
    ax2.set_xlabel(r"$\phi$")
    ax2.set_ylabel(r"$\dot\gamma_{\rm max}$  (s$^{-1}$)")

    stem = FIGDIR / "design_map_2026-08-31"
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".png"), dpi=200)
    plt.close(fig)

    csv = stem.with_suffix(".csv")
    with csv.open("w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"# generated {date.today().isoformat()} by scripts/run_design_map.py\n")
        fh.write(f"# eta_s={ETA_S:.7g} Pa.s (solved for onset at {GAMMA_DOT_REF} 1/s), "
                 f"sigma_star={SIGMA_STAR} Pa, phi_m={PHI_M}, phi_0={PHI_0}, n={N_EXP}\n")
        fh.write(f"# gamma_0 = 1.0 in the map; measured bracket {GAMMA_0_BRACKET}\n")
        fh.write("panel,x,series,y\n")
        for phi in phis_short:
            for A, r in zip(activation, shortfalls[phi]):
                fh.write(f"a,{A:.6g},{phi:.3f},{r:.6g}\n")
        for p, r in zip(phis_rate, rates):
            fh.write(f"b,{p:.4g},rate_max,{r:.6g}\n")

    print(f"eta_s solved            = {ETA_S:.6g} Pa.s  (PEG-like, not water-like)")
    print(f"gait band               = {GAIT_BAND[0]:g} to {GAIT_BAND[1]:g} 1/s")
    print("\nshortfall at phi = 0.40:")
    for A in (0.3, 1.0, 3.0, 10.0, 30.0):
        r = shortfall_curve(0.40, np.array([A]))[0]
        print(f"  A = {A:5.1f}   sigma_peak/sigma_steady = {r:.4f}")
    print("\nreachable rate ceiling:")
    for p in (0.40, 0.46, 0.50, 0.52, 0.54, 0.56):
        rm = reachable_rate(p)
        flag = "ABOVE band" if rm > GAIT_BAND[1] else "BELOW band"
        print(f"  phi = {p:.2f}   gamma_dot_max = {rm:10.1f} 1/s   {flag}")

    for path in (stem.with_suffix(".pdf"), stem.with_suffix(".png"), csv):
        h = hashlib.md5(path.read_bytes()).hexdigest()[:12]
        print(f"{h}  {path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
