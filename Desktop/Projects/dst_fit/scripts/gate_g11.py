"""Gate G11: measure the strain clock ``gamma_0``, the constant a flow curve cannot give.

Every other gate in this package fits steady flow curves. At steady state the
relaxation law's left-hand side is zero, so the two clocks that govern every
transient -- Richards et al. Eq. 4's strain clock ``gamma_0`` and time clock
``t_r`` -- drop out of the fit entirely. They are unmeasurable from a flow curve
in the strict sense: not poorly determined, absent.

Lin et al. (2015) deposited the measurement that does contain them. Reversing the
flow destroys the frictional contact network; the instantaneous viscosity then
climbs back to its steady value over accumulated strain, at four rates, for two
materials. Inverting those records through the *same* Wyart-Cates fit that G10
validated turns them into ``gamma_0``.

Why this is worth a gate rather than a script
---------------------------------------------
The Basilisk drop-impact solver runs exactly this law, with ``gamma_0 ~ 0.2``
(close to Han et al. 2018's 0.197 start-up scale; Richards et al. 2019 used
O(0.1) -- the earlier 'Richards' corn starch' attribution was wrong) and ``t_r = 0.24 +/- 0.05 s`` bracketed from
flow-cessation studies on other systems. Its own status card calls both
placeholders that "cannot be pinned without step-down rheometry", and its L9 row
shows the observable is not flat in ``gamma_0``: peak ``jam_frac`` runs
0.912 / 0.997 / 0.671 at ``gamma_0`` = 0.05 / 0.2 / 0.8. A measured value on
silica replaces a corn-starch literature number in the one term that solver
claims as its differentiator.

What was fixed before this gate ran, and what was not
-----------------------------------------------------
Being straight about this, because two of the four conditions below are
confirmations rather than discoveries.

``SOLVER_HANDOFF.md``, written earlier the same day, had already landmarked the
63% recovery point of these records **in viscosity space** and found the recovery
collapsing on strain rather than on time (spread 1.63x against 16.3x), with the
per-rate strain constant *falling* monotonically with rate -- the opposite of what
Richards' additive time term produces. It left one question open: whether that
wrong-sign trend is real or an artefact of ``eta`` being strongly nonlinear in
``f`` near jamming, which landmarking ``eta`` cannot separate. **G11b and G11c
ask whether those two observations survive the model inversion**, which is the
thing that answers the open question, and they extend both to a second material
and a simulation control that were never landmarked at all.

G11a and G11d were not suggested by anything and are genuine pre-registration.
All four tolerances were fixed before this script was first run.

  G11a  What the clock does and does not depend on. Two parts, because the first
        version of this condition **could not fail** and had to be replaced.

        The retracted version swept ``phi_0`` over 0.58 to 0.68, refitting the
        steady curve each time, and asked that ``gamma_eff`` move by less than 1%.
        It returned 7e-8. That is an identity, not a measurement: Wyart-Cates is
        exactly covariant under ``(phi, phi_m, phi_0) -> c (phi, phi_m, phi_0)``
        because ``phi_J`` is homogeneous of degree one, so ``phi/phi_J`` is
        invariant pointwise. ``Model.with_fixed`` moves the bounds with ``phi_0``,
        the fit slides along exactly that ray -- ``phi/phi_0`` held at 0.611341 and
        ``phi_m/phi`` at 1.299132 to six digits across the whole sweep -- and the
        inverted ``f`` is the same array. The gate was re-measuring algebra.

        The real theorem is stronger and is what is gated now. Writing
        ``u = phi_J/phi = 1/(1 - eta_r^(-1/n))``, the inversion is
        ``f = A + B u`` with ``A`` and ``B`` set by the packing triple. So ``f``
        is **affine in u for any triple**, and a relaxation fit with both
        amplitudes free absorbs any affine map of the ordinate. ``gamma_eff``
        therefore does not depend on ``phi``, ``phi_m`` or ``phi_0`` at all,
        separately or jointly, and depends only on ``n``.

        G11a-i  Re-measure at four packing triples imposed outright rather than
                fitted, spanning ``phi/phi_0`` = 0.32 to 0.87 and ``phi_m/phi`` =
                1.09 to 1.75, which is nowhere near the covariance ray. No
                ``gamma_eff`` may move by more than **1e-6**, since the expected
                answer is bit-identical and the tolerance covers only float noise.
        G11a-ii Sweep the exponent ``n`` over 1.8 to 2.5, which is the constant
                that does move it and is not measured for either system. Bar
                **5%**, set from the smallest effect this gate claims rather than
                from what the sweep returns; see ``N_TOL``.

        Together these are what decides whether ``gamma_0`` is quotable at all:
        it should inherit ``sigma*``'s identifiability rather than ``phi_m``'s
        degeneracy. Proved algebraically in ``tests/test_transient.py``; this
        checks it end to end on real fits.

        The retraction is recorded rather than quietly edited because it is the
        same failure mode this project keeps finding elsewhere -- a gate that
        names the wrong quantity -- and it was caught by review, not by the gate.

  G11b  The recovery collapses on strain, not on time: across the measured rate
        span the spread of ``gamma_eff`` is smaller than the spread of
        ``t_eff = gamma_eff / gamma_dot``, for both materials. An ordering claim
        with no threshold, because a threshold here would be arbitrary and the
        ordering is the whole content.

  G11c  Richards Eq. 4's additive form holds: fitting
        ``1/gamma_eff = 1/gamma_0 + 1/(gamma_dot t_r)`` returns ``t_r > 0``, for
        both materials. A negative ``t_r`` is not a large error bar; it means
        ``gamma_eff`` falls with rate, which that functional form cannot produce
        at any ``gamma_0`` and any ``t_r``.

  G11d  The answer is not an artefact of the fit window. Refitting each record
        truncated at ``gamma = 3`` instead of its full length must move
        ``gamma_eff`` by less than a **factor of 1.5**. These records are sampled
        uniformly in strain, so about 90% of the samples sit in the plateau above
        ``gamma = 1`` and carry no information about the clock; a constant that
        moves when that ballast is removed was measuring the ballast. 1.5x is
        set by what the number is for: the solver sweeps ``gamma_0`` over a factor
        of 16, so a 1.5x window is comfortably inside what it resolves, while
        still being far tighter than the effect being claimed.

Run:  python scripts/gate_g11.py
Exit code is 0 only if all four pass, on both materials.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dst_fit import fit_model, thickening_branch  # noqa: E402
from dst_fit.io_literature import read_lin_reversal, read_lin_transients  # noqa: E402
from dst_fit.models import WC_EXPONENT, WYART_CATES, frictional_fraction  # noqa: E402
from dst_fit.transient import (  # noqa: E402
    admissible_window,
    fit_recovery,
    fit_strain_clock,
    frictionless_relative_viscosity,
    half_recovery_strain,
    smooth,
    smoothing_samples,
)

ROOT = Path(__file__).resolve().parents[1]
LIN_ROOT = ROOT / (
    "Hydrodynamic and Contact Contributions to Continuous Shear Thickening "
    "in Colloidal Suspensions"
)

PHI0_DEFAULT = 0.64
PHI0_SWEEP = (0.58, 0.61, 0.64, 0.68)

# G11a-i. Packing triples chosen to sit well OFF the covariance ray: across them
# phi/phi_0 runs 0.32 to 0.87 and phi_m/phi runs 1.09 to 1.75, so a gate that
# passed only because the fit slid along the ray cannot pass here. The clock is
# affine-invariant, so the expected answer is bit-identical and the tolerance only
# has to cover floating-point noise.
PACKING_TRIPLES = (
    (0.40, 0.52, 0.64),
    (0.20, 0.35, 0.62),
    (0.55, 0.60, 0.63),
    (0.30, 0.45, 0.58),
)
PACKING_TOL = 1e-6

# G11a-ii. The exponent is the one fixed constant the clock genuinely depends on,
# and it is not measured for either system: n = 2 is Maron-Pierce, Guy et al. fit
# 2.2(2) to their own upper branch, and values to 2.5 appear in the literature.
#
# 5% is set by what this gate claims, not by what the sweep returns. The smallest
# effect asserted anywhere in G11 is the 1.70x rate trend that survives the
# two-component correction; 5% is more than an order of magnitude below that, so
# an n-sensitivity under the bar cannot account for any conclusion here. It is
# deliberately not set to the measured value.
N_SWEEP = (1.8, 2.0, 2.2, 2.5)
N_TOL = 0.05

TRUNCATION_TOL = 1.5  # G11d, factor
TRUNCATE_AT = 3.0  # G11d, strain

# What the Basilisk drop-impact solver currently runs. NOT Richards' value: Richards
# et al. 2019 used gamma_0 = O(0.1); 0.2 matches Han et al. 2018's 0.197, a START-UP
# scale in a different functional form (corrected 2026-09-27, PRIOR_ART_STRAIN_CLOCK). The point of comparison:
# this gate exists to replace it with a number measured on the solver's own material.
SOLVER_GAMMA_0 = 0.2

# transient panel -> steady panel measured on the same sample
PAIRS = (("Fig_2b", "Fig_2a"), ("Fig_2d", "Fig_2c"))


def steady_fit(system: str, phi_0: float):
    """The G10 fit: Wyart-Cates on relative viscosity against reduced stress."""
    split = read_lin_reversal(LIN_ROOT / "Fig2" / "Fig2", system)
    model = WYART_CATES.with_fixed(eta_s=1.0, phi_0=phi_0)
    br = thickening_branch(split.sweep_eta_r)
    fit = fit_model(model, split.sweep_sigma, split.sweep_eta_r, branch=br, n_starts=24)
    phi, phi_m, sigma_star = (float(v) for v in fit.params)
    return split, phi, phi_m, sigma_star


def fit_recovery_in_viscosity(strain, eta_r, phi, phi_m, phi_0):
    """The same exponential fitted to ``eta_r`` directly, for comparison only.

    This is what landmarking the viscosity measures. It is *not* the ``gamma_0``
    of the relaxation law: ``eta`` is strongly nonlinear in ``f`` near jamming, so
    its e-folding strain differs from ``f``'s by a factor that depends on how
    close the sample sits to its jamming point. Reported so the size of that
    factor is a measured number rather than an argument.
    """
    win = admissible_window(strain, eta_r, phi, phi_m, phi_0, require_admissible=False)
    es = smooth(eta_r, smoothing_samples(strain))
    g = strain[win.i_min :]
    y = es[win.i_min :]
    g0 = float(g[0])

    def residual(p):
        y0, y_ss, g_eff = p
        return y_ss + (y0 - y_ss) * np.exp(-(g - g0) / g_eff) - y

    span = float(g[-1] - g0)
    sol = least_squares(
        residual,
        [float(y[0]), float(np.median(y[-max(1, y.size // 20) :])), span / 5.0],
        bounds=([0.0, 0.0, 1e-4], [1e3, 1e3, 10.0 * span]),
    )
    return float(sol.x[2])


def measure(
    panel: str,
    steady_panel: str,
    phi_0: float,
    gamma_fit_max: float = np.inf,
    n: float = WC_EXPONENT,
):
    """Fit every admitted rate of one transient panel. Returns recoveries and skips."""
    split, phi, phi_m, sigma_star = steady_fit(steady_panel, phi_0)
    trans = read_lin_transients(LIN_ROOT, panel)
    recoveries, skipped = [], []
    for rate, g, eta in trans:
        try:
            recoveries.append(
                fit_recovery(
                    g, eta, rate, phi, phi_m, phi_0, n=n, gamma_fit_max=gamma_fit_max
                )
            )
        except ValueError as exc:
            skipped.append((rate, str(exc)))
    return split, trans, phi, phi_m, sigma_star, recoveries, skipped


def measure_with_packing(panel: str, triple: tuple[float, float, float]):
    """Re-measure every rate at a packing triple supplied outright.

    Deliberately does *not* refit the steady curve. Refitting would slide the
    triple back onto the covariance ray, which is exactly the identity that made
    the first version of G11a unfailable. Here the triple is imposed, so the
    ratios really do change and the invariance really is being tested.
    """
    phi, phi_m, phi_0 = triple
    trans = read_lin_transients(LIN_ROOT, panel)
    out = []
    for rate, g, eta in trans:
        try:
            out.append(fit_recovery(g, eta, rate, phi, phi_m, phi_0))
        except ValueError:
            continue
    return out


def landmarks(trans, phi, phi_m, phi_0) -> dict[float, tuple[float, float]]:
    """Fit-free half-recovery strain per rate, skipping rates the model cannot reach."""
    out: dict[float, tuple[float, float]] = {}
    for rate, g, eta in trans:
        try:
            g50, f_plateau = half_recovery_strain(g, eta, phi, phi_m, phi_0)
        except ValueError:
            continue
        if np.isfinite(g50):
            out[rate] = (g50, f_plateau)
    return out


def spread(values) -> float:
    """Max over min. The scale-free way to compare two collapses."""
    v = np.asarray(list(values), dtype=float)
    return float(v.max() / v.min())


def main() -> int:
    failures: list[str] = []
    print("Gate G11: the strain clock of the frictional network")
    print("Lin et al. 2015, PRL 115, 228304; DataShare 10283/886")
    print("Richards et al. 2019, PRL 123, 038004, Eq. 4 is the law being tested")
    print()
    print("Every +/- printed below is linearised and assumes independent samples.")
    print("These records carry 9 000 to 23 000 samples smoothed over ~100 of them, so")
    print("consecutive residuals are strongly correlated and those intervals are severe")
    print("underestimates. The honest uncertainty on the clock is the spread across the")
    print("estimators in G11a and G11d, not the +/-.\n")
    delivered: dict[str, list[float]] = {}

    for panel, steady_panel in PAIRS:
        split, trans, phi, phi_m, sigma_star, recs, skipped = measure(
            panel, steady_panel, PHI0_DEFAULT
        )
        floor = frictionless_relative_viscosity(phi, PHI0_DEFAULT)
        print(f"=== {panel}: {trans.label}   (steady fit from {steady_panel})")
        print(
            f"  steady fit  phi = {phi:.4f}  phi_m = {phi_m:.4f}  "
            f"sigma* = {sigma_star:.4g}  (phi_0 = {PHI0_DEFAULT} fixed)"
        )
        print(
            f"  the model's frictionless floor is eta_r = {floor:.3f}; "
            f"nothing below it is invertible at any f"
        )

        # ---------------- the headline observation, not gated ----------------
        print("\n  how far the reversal dip falls below the model's own floor:")
        for rate, g, eta in trans:
            try:
                win = admissible_window(g, eta, phi, phi_m, PHI0_DEFAULT)
            except ValueError as exc:
                print(f"    gamma_dot = {rate:<5g}  ENTIRELY below the floor -- {exc}")
                continue
            print(
                f"    gamma_dot = {rate:<5g}  min eta_r = {win.eta_min:6.3f} "
                f"= {win.depth_below_floor:.3f}x the floor   "
                f"admissible from gamma = {win.gamma_start:.4f} "
                f"({100 * win.admissible_fraction:.1f}% of the record)"
            )
        for rate, why in skipped:
            print(f"    gamma_dot = {rate:<5g}  NOT FITTED: {why.split('. ')[0]}")

        if len(recs) < 3:
            failures.append(
                f"{panel}: only {len(recs)} rates admitted, cannot separate two clocks"
            )
            print()
            continue

        # ---------------- the two estimators ----------------
        marks = landmarks(trans, phi, phi_m, PHI0_DEFAULT)
        print("\n  the recovery, measured two ways. gamma_eff is the exponential")
        print("  constant fitted in f-space; gamma_50 is the fit-free strain from the")
        print("  reversal minimum to half the plateau, which needs no window at all.")
        print(f"    {'gdot':>6} {'gamma_eff':>10} {'+/-':>8} {'gamma_50':>9} "
              f"{'f_plateau':>10} {'t_eff [s]':>10}")
        for r in recs:
            g50 = marks.get(r.rate, (np.nan, np.nan))[0]
            print(f"    {r.rate:6g} {r.gamma_eff:10.4f} {r.gamma_eff_err:8.4f} "
                  f"{g50:9.4f} {r.f_plateau:10.3f} {r.t_eff:10.4f}")

        # A record whose inferred plateau in f lies outside [0, 1] is outside the
        # model's range: no frictional fraction reproduces it. Such a record still
        # yields a relaxation constant, because the clock is affine-invariant and
        # is really a property of u = phi_J/phi rather than of f, and it is
        # reported above. But it cannot be used to *test* a law of the model whose
        # range it falls outside, so it is excluded here. The rule is a property
        # of the model, not of the answer, and it excludes exactly one record.
        usable = [r for r in recs if 0.0 <= r.f_plateau <= 1.0]
        dropped = [r for r in recs if r not in usable]
        for r in dropped:
            print(f"\n  gamma_dot = {r.rate:g} is excluded from G11c and from the "
                  f"quoted clock: its plateau in f is {r.f_plateau:.3f}, outside "
                  f"[0, 1], so the model has no state that describes this record.")
        if len(usable) < 3:
            failures.append(
                f"{panel}: only {len(usable)} records lie inside the model's range, "
                f"which is too few to separate two clocks"
            )
            print()
            continue

        rates = np.array([r.rate for r in usable])
        g_eff = np.array([r.gamma_eff for r in usable])
        g_err = np.array([r.gamma_eff_err for r in usable])
        clock = fit_strain_clock(rates, g_eff, g_err)

        # Whether the clock even moves monotonically with rate. G11c fits a
        # straight line through 1/gamma_eff against 1/gamma_dot and reads its sign;
        # if the points are not ordered, that sign is describing scatter. Reported
        # separately so a failure for the wrong reason is visible as one.
        order = np.argsort(rates)
        seq = g_eff[order]
        monotone = bool(np.all(np.diff(seq) < 0) or np.all(np.diff(seq) > 0))

        # ---------------- G11b ----------------
        # Judged on the fit-free landmark. gamma_eff is measured over a window
        # whose *start* moves with rate, because a slower recovery climbs back
        # above the frictionless floor later, so a rate trend in gamma_eff could
        # in principle be a trend in the window. gamma_50 cannot carry that.
        mrates = np.array(sorted(marks))
        g50s = np.array([marks[r][0] for r in mrates])
        t50s = g50s / mrates
        strain_spread, time_spread = spread(g50s), spread(t50s)
        rate_span = float(mrates.max() / mrates.min())
        strain_tighter = strain_spread < time_spread
        print("\n  G11b  which clock governs, on the fit-free landmark")
        print(f"        rate span {rate_span:.0f}x over {mrates.size} rates. A pure "
              f"strain clock gives gamma_50 spread 1x and t_50 spread {rate_span:.0f}x;")
        print("        a pure time clock gives the reverse.")
        print(f"        measured: gamma_50 spread {strain_spread:.2f}x, "
              f"t_50 spread {time_spread:.2f}x  ->  "
              f"{'STRAIN' if strain_tighter else 'TIME'}   "
              f"{'PASS' if strain_tighter else 'FAIL'}")
        print(f"        the fitted constant agrees: gamma_eff spread "
              f"{clock.strain_spread:.2f}x against t_eff {clock.time_spread:.2f}x")
        if not strain_tighter:
            failures.append(
                f"G11b[{panel}]: t_50 collapses tighter ({time_spread:.2f}x) "
                f"than gamma_50 ({strain_spread:.2f}x)"
            )

        # ---------------- G11c ----------------
        print("\n  G11c  Richards Eq. 4 additive form")
        print(clock.report())
        print(f"        gamma_eff against ascending rate: "
              f"{', '.join(f'{v:.4f}' for v in seq)}  -> "
              f"{'monotone' if monotone else 'NOT MONOTONE'}")
        if not monotone:
            print("        With the clocks unordered in rate, the straight line's sign")
            print("        is describing scatter rather than a trend. This is a failure")
            print("        of the test, not a refutation of the law: on this material")
            print("        the additive form is untested, not contradicted.")
        elif not clock.form_holds:
            print("        The clocks fall monotonically with rate, which the additive")
            print("        time term cannot produce at any gamma_0 and any t_r.")
        print(f"        {'PASS' if clock.form_holds else 'FAIL'}")
        if not clock.form_holds:
            reason = ("gamma_eff falls monotonically with rate, which the additive "
                      "time term cannot produce" if monotone else
                      f"gamma_eff is not monotone in rate (r^2 = {clock.r_squared:.2f}), "
                      "so the form is untested here rather than refuted")
            failures.append(f"G11c[{panel}]: t_r = {clock.t_r:.3f} s is negative; {reason}")

        # ---------------- G11a ----------------
        print("\n  G11a  what the clock does and does not depend on")
        print("        (i) the packing triple, varied OFF the covariance ray")
        by_triple = {}
        for triple in PACKING_TRIPLES:
            rr = measure_with_packing(panel, triple)
            by_triple[triple] = {r.rate: r.gamma_eff for r in rr}
            ph, phm, p0 = triple
            print(f"        (phi, phi_m, phi_0) = ({ph:.2f}, {phm:.2f}, {p0:.2f})   "
                  f"phi/phi_0 = {ph / p0:.4f}, phi_m/phi = {phm / ph:.4f}")
        worst_pack, worst_pack_rate = 0.0, None
        for r in recs:
            vals = np.array(
                [by_triple[t].get(r.rate, np.nan) for t in PACKING_TRIPLES], dtype=float
            )
            if not np.all(np.isfinite(vals)):
                failures.append(f"G11a[{panel}]: rate {r.rate:g} dropped out of the sweep")
                continue
            rel = float((vals.max() - vals.min()) / np.mean(vals))
            if rel > worst_pack:
                worst_pack, worst_pack_rate = rel, r.rate
        where = ("identical at every rate" if worst_pack_rate is None
                 else f"at gamma_dot = {worst_pack_rate:g}")
        print("        these ratios span 0.32-0.87 and 1.09-1.75, so this is not the ray.")
        print(f"        worst relative spread in gamma_eff: {worst_pack:.2e} ({where})")
        print(f"        tolerance {PACKING_TOL:.0e}   "
              f"{'PASS' if worst_pack <= PACKING_TOL else 'FAIL'}")
        if worst_pack > PACKING_TOL:
            failures.append(
                f"G11a-i[{panel}]: gamma_eff moves {worst_pack:.2e} across packing triples"
            )

        print("\n        (ii) the divergence exponent n, which is the one constant it "
              "does depend on")
        by_n = {}
        for nn in N_SWEEP:
            _, _, _, _, _, rr, _ = measure(panel, steady_panel, PHI0_DEFAULT, n=nn)
            by_n[nn] = {r.rate: r.gamma_eff for r in rr}
        print(f"        {'gdot':>6}" + "".join(f"{f'n={nn}':>10}" for nn in N_SWEEP)
              + f"{'spread':>9}")
        worst_n, worst_n_rate = 0.0, None
        for r in recs:
            vals = np.array([by_n[nn].get(r.rate, np.nan) for nn in N_SWEEP], dtype=float)
            if not np.all(np.isfinite(vals)):
                continue
            rel = float((vals.max() - vals.min()) / np.mean(vals))
            print(f"        {r.rate:6g}" + "".join(f"{v:10.4f}" for v in vals)
                  + f"{rel:8.2%}")
            if rel > worst_n:
                worst_n, worst_n_rate = rel, r.rate
        print(f"        worst {worst_n:.2%} at gamma_dot = {worst_n_rate}   "
              f"tolerance {N_TOL:.0%}   "
              f"{'PASS' if worst_n <= N_TOL else 'FAIL'}")
        if worst_n > N_TOL:
            failures.append(
                f"G11a-ii[{panel}]: gamma_eff moves {worst_n:.1%} across n = "
                f"{min(N_SWEEP)} to {max(N_SWEEP)}"
            )

        # ---------------- G11d ----------------
        print(f"\n  G11d  stability against truncating the fit at gamma = {TRUNCATE_AT:g}")
        _, _, _, _, _, recs_short, _ = measure(
            panel, steady_panel, PHI0_DEFAULT, gamma_fit_max=TRUNCATE_AT
        )
        short = {r.rate: r.gamma_eff for r in recs_short}
        # Same set G11c uses, and for the same reason: this asks whether *the
        # quoted answer* is an artefact of the window, and a record outside the
        # model's range is not part of the quoted answer. Its ratio is printed
        # below the gate so the exclusion cannot hide a bad number.
        worst_f, worst_r = 1.0, None
        for r in usable:
            if r.rate not in short:
                failures.append(f"G11d[{panel}]: rate {r.rate:g} lost under truncation")
                continue
            ratio = max(short[r.rate], r.gamma_eff) / min(short[r.rate], r.gamma_eff)
            print(f"        gamma_dot = {r.rate:<5g}  full {r.gamma_eff:.4f}  "
                  f"truncated {short[r.rate]:.4f}   {ratio:.3f}x")
            if ratio > worst_f:
                worst_f, worst_r = ratio, r.rate
        for r in dropped:
            if r.rate in short:
                ratio = max(short[r.rate], r.gamma_eff) / min(short[r.rate], r.gamma_eff)
                print(f"        gamma_dot = {r.rate:<5g}  full {r.gamma_eff:.4f}  "
                      f"truncated {short[r.rate]:.4f}   {ratio:.3f}x   "
                      f"(excluded: outside the model's range, not gated)")
        print(f"        worst {worst_f:.3f}x at gamma_dot = {worst_r}   "
              f"tolerance {TRUNCATION_TOL:.1f}x   "
              f"{'PASS' if worst_f <= TRUNCATION_TOL else 'FAIL'}")
        if worst_f > TRUNCATION_TOL:
            failures.append(
                f"G11d[{panel}]: gamma_eff moves {worst_f:.2f}x when the fit window "
                f"is truncated at gamma = {TRUNCATE_AT:g}"
            )

        # ---------------- what the solver gets ----------------
        allvals = (
            list(g_eff)
            + list(g50s)
            + [short[r.rate] for r in usable if r.rate in short]
        )
        delivered[trans.label] = allvals
        print("\n  the number the solver needs, bracketed over every estimator above:")
        print(f"        gamma_0 = {min(allvals):.2f} to {max(allvals):.2f} strain units, "
              f"from {len(allvals)} estimates")
        print(f"        the solver runs gamma_0 = {SOLVER_GAMMA_0}, taken from the "
              f"corn-starch literature (Han et al. 2018's 0.197 start-up scale; Richards et al. "
              f"2019 used O(0.1)). Every estimate here is larger, by "
              f"{min(allvals) / SOLVER_GAMMA_0:.1f}x to {max(allvals) / SOLVER_GAMMA_0:.1f}x.")

        # ---------------- observations ----------------
        print("\n  observations (not gated):")
        print("    the fitted plateau of f against the steady-state f the flow-curve fit")
        print("    predicts at the same rate, f_ss = exp(-sigma*/Sigma). Two different")
        print("    files, one model; they are not forced to agree.")
        for r in recs:
            j = int(np.argmin(np.abs(split.rev_rate - r.rate)))
            if not np.isclose(split.rev_rate[j], r.rate, rtol=1e-6):
                continue
            sigma_red = float(split.rev_total[j] * split.rev_rate[j])
            predicted = float(frictional_fraction(np.array([sigma_red]), sigma_star)[0])
            print(f"      gamma_dot = {r.rate:<5g}  fitted {r.f_plateau:.3f} "
                  f"+/- {r.f_plateau_err:.3f}   predicted {predicted:.3f}   "
                  f"gap {r.f_plateau - predicted:+.3f}")

        print("\n    the same exponential fitted to eta_r instead of f, i.e. what")
        print("    landmarking the viscosity measures. The ratio is the nonlinearity.")
        for r in recs:
            i = trans.rates.index(r.rate)
            g_eta = fit_recovery_in_viscosity(
                trans.strain[i], trans.eta_r[i], phi, phi_m, PHI0_DEFAULT
            )
            print(f"      gamma_dot = {r.rate:<5g}  f-space {r.gamma_eff:.4f}   "
                  f"eta-space {g_eta:.4f}   {g_eta / r.gamma_eff:.2f}x")

        print("\n    a post-hoc collapse, printed so that it is on the record as having")
        print("    been tested. gamma_50 * f_plateau looked constant on the silica; the")
        print("    PMMA is the out-of-sample test of it.")
        prod = np.array([marks[r][0] * marks[r][1] for r in mrates])
        print(f"      products {np.round(prod, 3).tolist()}, spread {spread(prod):.2f}x "
              f"against gamma_50's own {strain_spread:.2f}x -- "
              f"{'tighter' if spread(prod) < strain_spread else 'NOT tighter, so it is'}"
              f"{'' if spread(prod) < strain_spread else ' a coincidence of four points'}")
        print()

    # ----------------------------------------------------------------------
    if delivered:
        print("the strain clock across both materials:")
        for label, vals in delivered.items():
            print(f"  {label:<28s} gamma_0 = {min(vals):.2f} to {max(vals):.2f}")
        allv = [v for vals in delivered.values() for v in vals]
        print(f"  every estimate on both materials exceeds the solver's "
              f"{SOLVER_GAMMA_0}; the smallest is {min(allv) / SOLVER_GAMMA_0:.1f}x it.\n")

    if failures:
        print("FAILURES")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("all G11 conditions pass")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
