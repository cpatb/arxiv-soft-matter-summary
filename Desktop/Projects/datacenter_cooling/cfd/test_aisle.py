"""Checks on the vented-partition cell, and on the through-flow added 2026-08-31.

`aisle.py`'s docstring has claimed since it was written that the de Vahl Davis reduction
is "checked in test_aisle.py". This file did not exist. It does now.

It does NOT check the de Vahl Davis reduction: that benchmark lives in
`test_buoyant2d.py`, against `Cavity`, and is passing (Nu 1.115 vs 1.118 at Ra = 1e3,
2.203 vs 2.243 at Ra = 1e4). `aisle.py`'s docstring is corrected to point there rather
than to make a claim about this file that this file does not honour.

Each test carries an anchor that does not pass through the code under test: the stored
calibration JSON, an exact analytic flux, a conservation law, or a discriminating negative
case that must give a *different* answer.

Run:  python test_aisle.py        (or: python -m pytest test_aisle.py -q)
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from aisle import AisleCell
from run_provisioning import solve_to_phi

HERE = Path(__file__).resolve().parent
INLET = (0.10, 0.40)
OUTLET = (0.60, 0.90)


# --------------------------------------------------------------------------- #
# The anchor that does not pass through any of tonight's code: the committed
# calibration, measured before the through-flow existed.
# --------------------------------------------------------------------------- #
def test_zero_flux_reproduces_the_stored_calibration():
    stored = json.load(open(HERE / "cfd_recirc.json"))
    for rec in stored[:3]:
        a = AisleCell(N=40, H_room=1.0, eta=rec["eta"], Ra=rec["Ra"])
        solve_to_phi(a)
        # cfd_recirc.json rounds to four decimals, so that is the tolerance the
        # file itself can support -- not a tolerance chosen to make this pass.
        assert abs(a.recirc_fraction() - rec["phi"]) < 5e-5, (
            f"eta={rec['eta']} Ra={rec['Ra']}: {a.recirc_fraction():.6f} "
            f"vs stored {rec['phi']}"
        )
    print("[PASS] test_zero_flux_reproduces_the_stored_calibration")


def test_zero_flux_takes_the_old_code_path_exactly():
    """flux = 0 must skip the Dirichlet lift and the opening BCs entirely."""
    a = AisleCell(N=30, eta=0.4, Ra=4e3)
    assert a.flux == 0.0
    assert np.all(a.psi_bc == 0.0)
    b = AisleCell(N=30, eta=0.4, Ra=4e3, flux=0.0)
    for _ in range(200):
        a.step(1e-4)
        b.step(1e-4)
    assert np.array_equal(a.T, b.T) and np.array_equal(a.psi, b.psi)
    print("[PASS] test_zero_flux_takes_the_old_code_path_exactly")


# --------------------------------------------------------------------------- #
# Conservation: the flux that arrives is the flux that was asked for.
# --------------------------------------------------------------------------- #
def test_net_flux_is_conserved_through_every_horizontal_line():
    """The net upward volumetric flux through any horizontal cut must equal `flux`.

    Integrated from the SOLVED psi across the full width, at five heights spanning
    the domain. Two things make this the right form of the check rather than the
    obvious one:

    * Integrating v over the INLET SPAN ALONE at the first interior row is not the
      inflow. Between the floor and y = dy the entering stream has already spread
      sideways, and that version of the test reads 7.15 against a requested 8.0 --
      a real 11% shortfall in a quantity that is not conserved, not a solver defect.
    * Reading psi at the two walls would be a tautology: those are the Dirichlet
      values that were set. Differencing the interior solution across the width, and
      requiring the SAME answer at five heights, tests that the lift actually carried
      the through-flow into the interior and that the solve conserved it.
    """
    Q = 8.0
    a = AisleCell(N=60, eta=0.4, Ra=4e3, flux=Q, inlet=INLET, outlet=OUTLET)
    solve_to_phi(a)
    x = np.linspace(0.0, 1.0, a.nx)
    rows = [1, a.ny // 4, a.ny // 2, 3 * a.ny // 4, a.ny - 2]
    fluxes = []
    for j in rows:
        v = -np.gradient(a.psi[j, :], x)      # v = -dpsi/dx
        fluxes.append(float(np.trapezoid(v, x)))
    f = np.array(fluxes)
    assert np.all(np.abs(f - Q) / Q < 0.05), (
        f"net flux by row {np.round(f, 3)} against requested {Q}")
    assert (f.max() - f.min()) / Q < 0.02, (
        f"the net flux is not the same at every height: {np.round(f, 3)}")
    print(f"[PASS] test_net_flux_is_conserved_through_every_horizontal_line "
          f"({np.round(f, 3)} against {Q})")


def test_boundary_psi_closes_around_the_domain():
    a = AisleCell(N=30, eta=0.2, Ra=4e3, flux=5.0, inlet=INLET, outlet=OUTLET)
    p = a.psi_bc
    assert np.isclose(p[0, 0], 5.0) and np.isclose(p[-1, 0], 5.0)
    assert np.isclose(p[0, -1], 0.0) and np.isclose(p[-1, -1], 0.0)
    assert np.allclose(p[:, 0], 5.0) and np.allclose(p[:, -1], 0.0)
    print("[PASS] test_boundary_psi_closes_around_the_domain")


def test_a_partition_standing_on_a_nonzero_streamline_is_refused():
    """The solid is baked in at psi = 0, so the geometry must put it there."""
    try:
        AisleCell(N=30, eta=0.4, Ra=4e3, flux=5.0, inlet=(0.10, 0.70), outlet=OUTLET)
    except ValueError as e:
        assert "partition" in str(e)
        print("[PASS] test_a_partition_standing_on_a_nonzero_streamline_is_refused")
        return
    raise AssertionError("an inlet reaching under the partition was accepted")


# --------------------------------------------------------------------------- #
# Physics: the direction of every effect, and one case that must NOT move.
# --------------------------------------------------------------------------- #
def test_recirculation_falls_monotonically_with_provisioning():
    prev = None
    for p in (0.0, 0.5, 1.0, 2.0):
        a = AisleCell(N=40, eta=0.0, Ra=4e3, flux=p * 3.2035,
                      inlet=INLET, outlet=OUTLET)
        solve_to_phi(a)
        phi = a.recirc_fraction()
        if prev is not None:
            assert phi <= prev + 1e-4, f"phi rose from {prev:.4f} to {phi:.4f} at p={p}"
        prev = phi
    print("[PASS] test_recirculation_falls_monotonically_with_provisioning")


def test_advected_signal_is_lewis_independent_and_diffused_one_is_not():
    """The discriminating pair behind the 2026-08-31 floor retraction.

    At eta = 0 the inlet signal is advected, so changing the tracer's diffusivity by
    50x must not move it. At eta = 0.95 the partition shuts the advective path and
    what is left is conduction, so the same change must move it a lot. A version of
    the solver that ignored `tracer_diff` would pass the first and fail the second.
    """
    def tracer(eta, D):
        a = AisleCell(N=40, eta=eta, Ra=1e4, tracer_diff=D)
        solve_to_phi(a)
        return a.recirc_fraction_tracer(), a.recirc_fraction()

    lo0, _ = tracer(0.0, 0.02)
    hi0, _ = tracer(0.0, 1.0)
    assert abs(hi0 - lo0) / hi0 < 0.05, (
        f"eta = 0 moved {abs(hi0 - lo0) / hi0:.1%} with a 50x diffusivity change; "
        "an advected signal must be Lewis-independent"
    )

    lo1, phi_T1 = tracer(0.95, 0.02)
    hi1, _ = tracer(0.95, 1.0)
    assert hi1 / lo1 > 5.0, (
        f"eta = 0.95 moved only {hi1 / lo1:.1f}x with a 50x diffusivity change"
    )
    assert abs(hi1 - phi_T1) / phi_T1 < 0.02, (
        f"at Le = 1 the tracer should reproduce the temperature channel; "
        f"{hi1:.5f} vs {phi_T1:.5f}"
    )
    print(f"[PASS] test_advected_signal_is_lewis_independent_and_diffused_one_is_not "
          f"(eta=0: {lo0:.4f}->{hi0:.4f}; eta=0.95: {lo1:.5f}->{hi1:.5f}, "
          f"phi_T {phi_T1:.5f})")


def test_containment_reduces_recirculation_at_zero_flux():
    a = AisleCell(N=40, eta=0.0, Ra=4e3)
    b = AisleCell(N=40, eta=0.95, Ra=4e3)
    solve_to_phi(a)
    solve_to_phi(b)
    assert b.recirc_fraction() < a.recirc_fraction()
    print("[PASS] test_containment_reduces_recirculation_at_zero_flux")


if __name__ == "__main__":
    test_zero_flux_takes_the_old_code_path_exactly()
    test_boundary_psi_closes_around_the_domain()
    test_a_partition_standing_on_a_nonzero_streamline_is_refused()
    test_zero_flux_reproduces_the_stored_calibration()
    test_net_flux_is_conserved_through_every_horizontal_line()
    test_recirculation_falls_monotonically_with_provisioning()
    test_containment_reduces_recirculation_at_zero_flux()
    test_advected_signal_is_lewis_independent_and_diffused_one_is_not()
    print("\nALL AISLE TESTS PASSED")
