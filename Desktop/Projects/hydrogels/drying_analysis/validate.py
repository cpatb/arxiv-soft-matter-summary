"""
Validation suite: every module checked against synthetic data with KNOWN ground
truth (CLAUDE.md: "compare to a reference", "test the edges"). Run directly:

    python -m drying_analysis.validate      # from Desktop/Projects/hydrogels
    python validate.py                      # from inside the package dir

Deterministic: all randomness uses np.random.default_rng(SEED). Each check prints
PASS/FAIL with the recovered-vs-true numbers. Exit code 0 iff all pass.
"""
from __future__ import annotations

import sys

import numpy as np

try:  # allow running as a script or as a module
    from . import (calibration, crack_metrics, curl_radius,
                   indentation_modulus, kymograph, swelling_kinetics, toughness)
except ImportError:
    import calibration, crack_metrics, curl_radius  # type: ignore
    import indentation_modulus, kymograph, swelling_kinetics, toughness  # type: ignore

SEED = 20260706
_results = []


def check(name, ok, detail=""):
    _results.append(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  -- {detail}" if detail else ""))


# --------------------------------------------------------------------------
def test_calibration():
    # 200 px over 25 mm -> 8 px/mm exactly
    ppm = calibration.px_per_mm_from_fiducial((0, 0), (0, 200), 25.0)
    check("calibration: px/mm from fiducial", np.isclose(ppm, 8.0),
          f"got {ppm:.4f}, expected 8.0")
    # diagonal 3-4-5
    ppm2 = calibration.px_per_mm_from_fiducial((0, 0), (30, 40), 5.0)
    check("calibration: diagonal fiducial", np.isclose(ppm2, 10.0),
          f"got {ppm2:.4f}, expected 10.0")
    for bad in [((0, 0), (0, 0), 5.0), ((0, 0), (0, 10), 0.0)]:
        try:
            calibration.px_per_mm_from_fiducial(*bad)
            check("calibration: rejects degenerate input", False, f"{bad}")
            break
        except ValueError:
            pass
    else:
        check("calibration: rejects degenerate input", True)


def test_swelling_kinetics():
    L = 1.0e-3            # 1 mm slab
    D_true = 3.0e-10      # m^2/s
    t = np.linspace(1.0, 4000.0, 400)
    F = swelling_kinetics.crank_fixed_surface(t, D_true, L)

    early = swelling_kinetics.fit_D_early_sqrt(t, F, L)
    check("swelling: early sqrt(t) D recovers Crank input",
          np.isclose(early.D, D_true, rtol=0.06),
          f"D={early.D:.3e} vs {D_true:.3e}, r2={early.r2:.4f}")

    late = swelling_kinetics.fit_D_longtime(t, F, L)
    check("swelling: long-time exponential D recovers Crank input",
          np.isclose(late.D, D_true, rtol=0.05),
          f"D={late.D:.3e} vs {D_true:.3e}, r2={late.r2:.4f}")

    n = swelling_kinetics.fit_exponent_n(t, F)
    check("swelling: Fickian exponent n ~ 0.5", np.isclose(n, 0.5, atol=0.05),
          f"n={n:.3f}")

    # CONVENTION cross-check (catches a thickness/discretization bug that a pure
    # round-trip cannot): Robin sim with k_m -> inf must reduce to the fixed-
    # surface Crank series for the SAME full thickness L.
    t_lim = np.linspace(1.0, 4000.0, 200)
    F_robin_inf = swelling_kinetics._simulate_robin(t_lim, D_true, 1.0e3, L, n_nodes=81)
    F_crank = swelling_kinetics.crank_fixed_surface(t_lim, D_true, L)
    check("swelling: Robin(k_m->inf) reduces to Crank fixed-surface (same L)",
          np.allclose(F_robin_inf, F_crank, atol=0.02),
          f"max|dF|={np.max(np.abs(F_robin_inf - F_crank)):.4f}")

    # Robin-BC round trip: forward-simulate with known (D, k_m), recover them.
    # Bi = k_m (L/2)/D = 1.0  ->  k_m = 2 D / L.
    D2 = 2.0e-10
    km2 = 2.0 * D2 / L               # gives Bi = 1.0
    t2 = np.linspace(0.0, 8000.0, 120)
    F2 = swelling_kinetics._simulate_robin(t2, D2, km2, L, n_nodes=61)
    fit = swelling_kinetics.fit_D_km_robin(t2, F2, L, n_nodes=61)
    check("swelling: Robin (D,k_m) round-trip -> D",
          np.isclose(fit.D, D2, rtol=0.15), f"D={fit.D:.3e} vs {D2:.3e}")
    check("swelling: Robin round-trip -> Bi ~ 1.0",
          np.isclose(fit.Bi, 1.0, rtol=0.20), f"Bi={fit.Bi:.3f}")
    # Non-zero start time must NOT rescale the model (t=0 referencing).
    t3 = np.linspace(600.0, 8000.0, 100)
    F3 = swelling_kinetics._simulate_robin(t3, D2, km2, L, n_nodes=61)
    check("swelling: Robin sim referenced to t=0 (F>0 at first sampled t>0)",
          F3[0] > 0.02, f"F(t0={t3[0]:.0f}s)={F3[0]:.4f}")


def test_indentation():
    E_true = 5.0e3        # 5 kPa gel
    R = 2.0e-3            # 2 mm sphere radius
    # half-space round trip
    for delta in (1e-4, 3e-4):
        F = indentation_modulus.forward_indentation_force(E_true, R, delta)
        E = indentation_modulus.modulus_from_indentation(F, R, delta)
        check(f"indentation: half-space round-trip (delta={delta*1e3:.1f}mm)",
              np.isclose(E, E_true, rtol=1e-9), f"E={E:.2f} vs {E_true}")
    # finite thickness: forward with correction, invert with correction
    h = 1.0e-3
    delta = 2.0e-4
    Fh = indentation_modulus.forward_indentation_force(E_true, R, delta, h=h)
    Eh = indentation_modulus.modulus_from_indentation(Fh, R, delta, h=h)
    check("indentation: finite-thickness round-trip",
          np.isclose(Eh, E_true, rtol=1e-9), f"E={Eh:.2f} vs {E_true}")
    # ignoring finite thickness over-stiffens (E too high) -- sanity of sign
    E_ignored = indentation_modulus.modulus_from_indentation(Fh, R, delta)
    check("indentation: ignoring finite h inflates E (sign check)",
          E_ignored > E_true, f"E_ignored={E_ignored:.1f} > {E_true}")
    # buoyant weight sanity: steel sphere in water
    Fb = indentation_modulus.buoyant_weight(1e-3, 7850.0, 998.0)
    check("indentation: buoyant weight positive & plausible",
          0 < Fb < 1e-2, f"F={Fb:.3e} N")
    # HAND-COMPUTED reference (independent of the forward model, so it can catch
    # a wrong 16/9 prefactor): E = 9 F / (16 sqrt(R) delta^{3/2}).
    # F=1e-4 N, R=2e-3 m, delta=1e-4 m -> E = 9e-4/(16*sqrt(2e-3)*1e-6) = 1257.79 Pa
    E_ref = indentation_modulus.modulus_from_indentation(1e-4, 2e-3, 1e-4)
    check("indentation: matches hand-computed half-space E",
          np.isclose(E_ref, 1257.79, rtol=1e-3), f"E={E_ref:.2f} vs 1257.79 Pa")


def test_curl_radius():
    rng = np.random.default_rng(SEED)
    R_true = 120.0        # px
    cy, cx = 300.0, 250.0
    ang = np.linspace(np.deg2rad(200), np.deg2rad(340), 40)  # ~140 deg arc
    y = cy + R_true * np.sin(ang) + rng.normal(0, 0.5, ang.size)
    x = cx + R_true * np.cos(ang) + rng.normal(0, 0.5, ang.size)
    fit = curl_radius.fit_curl_radius(np.column_stack([y, x]), px_per_mm=8.0)
    check("curl: radius recovered from noisy arc",
          np.isclose(fit.radius_px, R_true, rtol=0.03),
          f"R={fit.radius_px:.2f}px vs {R_true}, rms={fit.rms_residual_px:.2f}")
    check("curl: mm conversion", np.isclose(fit.radius_mm, R_true / 8.0, rtol=0.03),
          f"{fit.radius_mm:.3f} mm")
    try:
        curl_radius.fit_curl_radius([[0, 0], [1, 1]])
        check("curl: rejects <3 points", False)
    except ValueError:
        check("curl: rejects <3 points", True)


def test_crack_metrics():
    # synthetic 3x3 grid of square fragments separated by 1-px cracks.
    ppm = 10.0            # px/mm
    cell = 50            # px -> 5 mm side -> 25 mm^2 per interior fragment
    gap = 2
    n = 3
    size = n * cell + (n + 1) * gap
    mask = np.zeros((size, size), dtype=bool)
    # draw crack lines (True) on a grid
    for i in range(n + 1):
        r = i * (cell + gap)
        mask[r:r + gap, :] = True
        mask[:, r:r + gap] = True
    m = crack_metrics.crack_metrics(mask, px_per_mm=ppm, min_area_px=50,
                                    exclude_border=True)
    # interior fragments = (n)^2 minus border-touching; border excluded leaves
    # the fully-interior ones. With a 3x3 grid all are border-adjacent except
    # none are fully interior -> allow >=1; check area of any detected cell.
    expected_area = (cell / ppm) ** 2       # 25 mm^2
    check("crack: fragment count > 0", m.n_fragments > 0, f"n={m.n_fragments}")
    if m.n_fragments > 0:
        check("crack: fragment area ~ 25 mm^2",
              np.isclose(m.median_area_mm2, expected_area, rtol=0.08),
              f"{m.median_area_mm2:.2f} vs {expected_area:.2f} mm^2")
    # all-zero mask (no cracks) -> one big fragment, but it touches border ->
    # excluded -> zero fragments; that's the documented edge behavior.
    empty = crack_metrics.crack_metrics(np.zeros((30, 30), bool), ppm,
                                        exclude_border=True)
    check("crack: single border-touching fragment excluded",
          empty.n_fragments == 0, f"n={empty.n_fragments}")


def test_toughness():
    G = 10.0e3            # shear modulus 10 kPa
    lam_c = 1.8
    # sample the neo-Hookean nominal stress sigma = G(lam - 1/lam^3) (pure shear)
    lam = np.linspace(1.0, lam_c, 60)
    sigma = G * (lam - 1.0 / lam**3)
    W_num = toughness.strain_energy_density(lam, sigma)
    W_ana = toughness.neo_hookean_W_pure_shear(G, lam_c)
    check("toughness: W integral matches neo-Hookean analytic",
          np.isclose(W_num, W_ana, rtol=1e-3),
          f"W_num={W_num:.1f} vs W_ana={W_ana:.1f} Pa")
    H0 = 5e-3
    Gam = toughness.toughness_pure_shear(W_num, H0)
    check("toughness: Gamma = W*H0", np.isclose(Gam, W_num * H0),
          f"Gamma={Gam:.2f} J/m^2")
    ell = toughness.network_length(Gam, 3 * G)   # E ~ 3G incompressible
    check("toughness: ell = Gamma/E positive length", ell > 0, f"ell={ell*1e3:.3f} mm")
    # ell_from_lfc cross-check: pure-shear f(lam_c); for brittle lam_c~1.2, f~1 so
    # ell ~ l_fc; hand value f(1.8)=(1.8^2+1/1.8^2-2)/6 = 0.26697...
    f_expected = (1.8**2 + 1.0 / 1.8**2 - 2.0) / 6.0
    ell_x = toughness.ell_from_lfc(1.0, 1.8)
    check("toughness: ell_from_lfc uses pure-shear f(lam_c)",
          np.isclose(ell_x, f_expected, rtol=1e-9), f"f={ell_x:.5f} vs {f_expected:.5f}")
    # Direction check: ell = Gamma/E = l_fc * (W_c/E) with W_c/E = f < 1 for any
    # finite stretch, and f -> 0 as lam_c -> 1. So ell << l_fc for brittle
    # low-stretch gels, and grows toward l_fc for stretchy ones. (NB: hydrogel_MOC
    # section 5 prose says "ell ~ l_fc for brittle low-stretch" -- that is inverted;
    # the formula, which is the standard Gamma/E vs Gamma/W_c ratio, is correct.
    # Flagged in drying_analysis/NOTES_moc_flags.md.)
    f_brittle = toughness.ell_from_lfc(1.0, 1.15)
    f_stretch = toughness.ell_from_lfc(1.0, 2.5)
    check("toughness: ell/l_fc monotonic & < 1 (grows with stretch)",
          f_brittle < f_stretch < 1.0 and f_brittle < 0.1,
          f"f(1.15)={f_brittle:.4f} < f(2.5)={f_stretch:.4f} < 1")
    # strain_energy_density must REJECT a non-monotonic (unloading) trace.
    try:
        toughness.strain_energy_density([1.0, 1.5, 1.2], [0.0, 1.0, 0.5])
        check("toughness: rejects non-monotonic stretch", False)
    except ValueError:
        check("toughness: rejects non-monotonic stretch", True)

    # fractocohesive knee: plateau 1.6 for c<l_fc=2mm, then falling
    l_true = 2.0
    c = np.array([0.3, 0.6, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0])
    y = np.where(c <= l_true, 1.6, 1.6 - 0.15 * (c - l_true))
    y = y + np.random.default_rng(SEED).normal(0, 0.01, c.size)
    knee = toughness.fractocohesive_knee(c, y)
    check("toughness: fractocohesive knee ~ l_fc",
          np.isclose(knee.l_fc, l_true, atol=0.4),
          f"l_fc={knee.l_fc:.2f} vs {l_true}")


def test_kymograph():
    # synthetic frame stack: a crack descends from the surface (row 0) then
    # deflects (depth stops increasing) -> should classify 'crack->peel'.
    T, H, W = 60, 80, 40
    frames = np.full((T, H, W), 200.0)     # bright background
    col = 20
    max_depth = 40
    for t in range(T):
        # depth grows linearly for first half, then holds (deflection)
        if t < T // 2:
            depth = int(max_depth * t / (T // 2))
        else:
            depth = max_depth
        frames[t, :depth + 1, col - 1:col + 2] = 30.0   # dark crack column
    kymo = kymograph.build_kymograph(frames, col, band=1)
    check("kymograph: shape is (H, T)", kymo.shape == (H, T), f"{kymo.shape}")
    traj = kymograph.crack_tip_trajectory(kymo, surface_row=0)
    check("kymograph: onset at surface (z_first ~ 0)",
          abs(traj.z_first) <= 3, f"z_first={traj.z_first}")
    # the tip must actually TRACK the descent, not sit at the surface (this is
    # what the cancelling tip/classifier bugs broke): max detected depth ~ 40.
    check("kymograph: tip tracks the descending leading edge",
          traj.depth.max() >= 0.9 * max_depth,
          f"max depth={traj.depth.max():.0f} vs {max_depth}")
    label = kymograph.classify_ordering(traj)
    check("kymograph: classifies descend-then-flatten as crack->peel",
          label == "crack->peel", f"got '{label}'")
    # crack-ONLY case: descends all the way to the last frame (never settles) ->
    # must NOT be labelled crack->peel (independent guard on the classifier).
    framesC = np.full((T, H, W), 200.0)
    for t in range(T):
        depth = int((max_depth) * t / (T - 1))     # descends until the final frame
        framesC[t, :depth + 1, col - 1:col + 2] = 30.0
    kymoC = kymograph.build_kymograph(framesC, col, band=1)
    trajC = kymograph.crack_tip_trajectory(kymoC, surface_row=0)
    check("kymograph: monotonic descent -> crack-only (not crack->peel)",
          kymograph.classify_ordering(trajC) == "crack-only",
          f"got '{kymograph.classify_ordering(trajC)}'")
    # peel-first case: a subsurface feature APPEARS (after a pristine frame 0),
    # never touching the surface -> peel-first. Frame 0 is the pre-event reference.
    frames2 = np.full((T, H, W), 200.0)
    frames2[5:, 30:33, col - 1:col + 2] = 30.0
    kymo2 = kymograph.build_kymograph(frames2, col, band=1)
    traj2 = kymograph.crack_tip_trajectory(kymo2, surface_row=0, ref_col=0)
    check("kymograph: subsurface onset -> peel-first",
          kymograph.classify_ordering(traj2) == "peel-first",
          f"z_first={traj2.z_first}")
    # surface-at-bottom invariance: same descent with surface_row=H-1 and rows
    # flipped must still read as crack->peel (bug-5 guard).
    frames_flip = frames[:, ::-1, :]
    kymo_flip = kymograph.build_kymograph(frames_flip, col, band=1)
    traj_flip = kymograph.crack_tip_trajectory(kymo_flip, surface_row=H - 1)
    check("kymograph: surface-at-bottom gives same crack->peel verdict",
          kymograph.classify_ordering(traj_flip) == "crack->peel",
          f"got '{kymograph.classify_ordering(traj_flip)}'")


def main():
    print("drying_analysis validation suite  (seed", SEED, ")")
    print("=" * 60)
    test_calibration()
    test_swelling_kinetics()
    test_indentation()
    test_curl_radius()
    test_crack_metrics()
    test_toughness()
    test_kymograph()
    print("=" * 60)
    n_pass = sum(_results)
    n_tot = len(_results)
    print(f"{n_pass}/{n_tot} checks passed")
    return 0 if n_pass == n_tot else 1


if __name__ == "__main__":
    sys.exit(main())
