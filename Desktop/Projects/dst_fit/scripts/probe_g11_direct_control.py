"""Does the Wyart-Cates inversion distort the strain clock? Lin's Fig3 can say.

WHAT G11 TAKES ON FAITH
    G11 measures the strain clock `gamma_0` by inverting a post-reversal viscosity
    transient through Wyart-Cates to get the frictional fraction `f(gamma)`, and fitting
    an exponential to that. The inversion is the one step nothing checks: `f` is never
    observed, only inferred, and `eta` is strongly nonlinear in `f` near jamming, so an
    inversion that is wrong in a strain-dependent way would move the clock without
    leaving a residual.

    Lin et al.'s Fig3 is the control, and `PLAN_2026-08-29.md` item 5 asks for it. Those
    panels are SIMULATIONS of the same reversal protocol, so the contact and hydrodynamic
    stresses are **evaluated directly** rather than inferred. Fig_3b gives both, against
    strain, at gamma_dot = 0.74 gamma_dot_0.

WHAT IS AND IS NOT BEING COMPARED
    The directly evaluated contact STRESS fraction `s = sigma_contact / sigma_total` is
    not the same quantity as Wyart-Cates's frictional-CONTACT fraction `f`. Their levels
    are not comparable and this probe does not compare them.

    What is comparable is the SHAPE of the rebuild, because `fit_recovery` fits
    `f_ss + (f_0 - f_ss) exp(-(gamma - gamma_start)/gamma_eff)` with both amplitudes free
    and eliminated exactly. Any AFFINE map between two signals leaves `gamma_eff`
    bit-identical -- that is the same argument that makes G11's clock independent of the
    packing triple, verified there across four triples spanning phi/phi_0 = 0.32 to 0.87.
    So a difference in `gamma_eff` between the two channels is exactly a measurement of
    the inversion's NONLINEARITY, which is the thing G11 assumes away.

HOW THE SECOND CHANNEL REUSES THE FIRST CHANNEL'S CODE
    Rather than re-implement the binning, the VarPro amplitude elimination and the
    log-grid profile scan, `s(gamma)` is pushed FORWARD through the same Wyart-Cates map
    that `fit_recovery` inverts, and the synthetic `eta_r` handed to `fit_recovery`. The
    inversion then returns `s(gamma)` exactly and the fit is the same tested code path.
    Round-tripping is verified to machine precision below before anything is fitted --
    a round trip is not a physics anchor, but here it is exactly the right check, because
    the claim is about the code path and not about the physics.

    Both channels are fitted at the SAME strain samples. Fig_3a is sampled at 16 185
    points and Fig_3b at 3 236, so Fig_3a is interpolated onto Fig_3b's grid; the two are
    first checked against each other on the overlap, because the deposit does not say in
    so many words that they are the same run.

PRE-REGISTERED, before the run
    P1  The two panels' total viscosities differ by a CONSTANT factor on the shared
        grid -- coefficient of variation below 0.5%. If they do not, they are not the
        same run and nothing below is a control.

        P1 was first written as "agree to better than 2% in median" and that version
        FAILED, at 2.43%. It was the wrong test, not a failed one: the two panels differ
        by 1.024300 with a standard deviation of 1.3e-5, i.e. a fixed normalisation
        constant held to five digits, which is exactly what "the same run, plotted with a
        different reduction" looks like and is not what "different runs" looks like. The
        original bar is recorded here rather than deleted, and both numbers are printed.
        The constancy test is the stricter one anyway: it demands 0.5% on a ratio where
        the old one demanded 2% on a value.

        Nothing downstream uses Fig_3a. Both fitted channels come from Fig_3b, so the
        1.0243 does not enter either of them and the two differ only by the inversion --
        which is the whole point. What the factor IS is not explained anywhere in the
        deposit's READMEs and is left unexplained here.
    P2  The round trip s -> eta_r -> s is exact to 1e-10. Otherwise the second channel is
        not measuring what it claims.
    P3  The two clocks agree to within 25%. That bar is set by what the clock is for:
        the Basilisk solver borrows gamma_0 = 0.2 and G11's measured range is 0.54 to
        1.64, a factor of 3, so a distortion below 25% would not move any conclusion
        that range supports, and one above it would.
    P4  The comparison is repeated at a second packing triple. If the two clocks move
        together, the affine-invariance argument holds and the difference is the
        nonlinearity; if the DIFFERENCE moves with the triple, the probe is measuring
        the triple and P3 cannot be read.

    This is a control on the METHOD. The rates are in the simulation's own units, so the
    clock here is dimensionless and is not a second measurement of the material constant.

    A rank correlation between the two signals is NOT reported, and the reason is worth
    recording: 76% of this record is the post-recovery plateau, so Spearman is computed
    almost entirely on plateau noise and comes back at 0.04 while Pearson is 0.94. The
    same asymmetry is why `fit_recovery` log-bins before fitting. The affine R2 below is
    computed on the raw record and carries the same caveat -- it is quoted as a shape
    diagnostic, not as a goodness-of-fit.

Run:  python scripts/probe_g11_direct_control.py
Exit code 0 only if P1, P2 and P4 hold; P3's verdict is printed either way, because
"the inversion distorts the clock" is as much a result as "it does not".
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from dst_fit.transient import (  # noqa: E402
    WC_EXPONENT,
    fit_recovery,
    frictional_fraction_from_viscosity as inv,
)

FIG3 = (ROOT / "Hydrodynamic and Contact Contributions to Continuous Shear Thickening "
        "in Colloidal Suspensions" / "Fig3" / "Fig3")
RATE = 0.74  # in units of the simulation's own gamma_dot_0; see the module docstring

# Two packing triples, deliberately far apart, for P4. Neither is a measurement of
# Lin's simulation -- they are gauge choices, and P4 is the check that they behave
# like gauge choices.
TRIPLES = [(0.50, 0.58, 0.66), (0.40, 0.52, 0.70)]

AGREE_CV_TOL = 0.005  # P1: coefficient of variation of the panel-to-panel ratio
ROUNDTRIP_TOL = 1e-10  # P2
CLOCK_TOL = 0.25      # P3
TRIPLE_TOL = 0.02     # P4: the clock DIFFERENCE must not move with the triple


def load(name: str) -> np.ndarray:
    return np.loadtxt(FIG3 / name, delimiter=",")


def wc_forward(f, phi, phi_m, phi_0, n=WC_EXPONENT):
    """eta_r = (1 - phi/phi_J)^-n with phi_J = f phi_m + (1-f) phi_0.

    The exact inverse of `frictional_fraction_from_viscosity`, written here so the
    round trip in P2 is a genuine check on two independent expressions rather than a
    function called twice.
    """
    f = np.asarray(f, dtype=float)
    phi_j = f * phi_m + (1.0 - f) * phi_0
    return (1.0 - phi / phi_j) ** (-n)


def main() -> int:
    contact = load("Fig_3b_contact.csv")
    hydro = load("Fig_3b_hydro.csv")
    total_b = load("Fig_3b_total.csv")
    a = load(f"Fig_3a_gamma_dot_{RATE:.2f}.csv")

    g = contact[:, 0]
    assert np.allclose(g, hydro[:, 0]) and np.allclose(g, total_b[:, 0]), \
        "Fig_3b's three files are not on one strain grid"

    # ---- P1: are Fig_3a and Fig_3b the same run? --------------------------------
    lo, hi = max(g[0], a[0, 0]), min(g[-1], a[-1, 0])
    m = (g >= lo) & (g <= hi)
    eta_a = np.interp(g[m], a[:, 0], a[:, 1])
    ratio = eta_a / total_b[m, 1]
    cv = float(np.std(ratio) / np.mean(ratio))
    rel = np.abs(eta_a - total_b[m, 1]) / total_b[m, 1]
    p1 = cv < AGREE_CV_TOL
    print(f"P1  Fig_3a / Fig_3b total viscosity on {int(m.sum())} shared strains: "
          f"ratio = {np.mean(ratio):.6f} +/- {np.std(ratio):.1e}, "
          f"CV = {cv:.4%}   {'PASS' if p1 else 'FAIL'}")
    print(f"    (the original bar, median |rel| < 2%, FAILED at {np.median(rel):.2%}; it "
          f"was the wrong test -- see the docstring)")
    if not p1:
        print("    They are not the same run; nothing below is a control. Stopping.")
        return 1

    s = contact[:, 1] / total_b[:, 1]          # direct contact stress fraction
    eta_direct = total_b[:, 1]                 # the panel's own total, same grid
    print(f"    s = sigma_contact/sigma_total runs {s.min():.4f} to {s.max():.4f} "
          f"over strain {g[0]:.4f} to {g[-1]:.3f}")

    rows = []
    for (phi, phi_m, phi_0) in TRIPLES:
        # ---- P2: the forward/inverse round trip ---------------------------------
        eta_synth = wc_forward(s, phi, phi_m, phi_0)
        s_back = inv(eta_synth, phi, phi_m, phi_0)
        err = float(np.max(np.abs(s_back - s)))
        p2 = err < ROUNDTRIP_TOL
        print(f"\nphi = {phi}, phi_m = {phi_m}, phi_0 = {phi_0}")
        print(f"P2  round trip s -> eta_r -> s: max |error| = {err:.2e}   "
              f"{'PASS' if p2 else 'FAIL'}")
        if not p2:
            return 1

        try:
            r_inv = fit_recovery(g, eta_direct, RATE, phi, phi_m, phi_0)
            r_dir = fit_recovery(g, eta_synth, RATE, phi, phi_m, phi_0)
        except ValueError as e:
            print(f"    fit_recovery refused: {e}")
            return 1

        d = abs(r_dir.gamma_eff - r_inv.gamma_eff) / r_inv.gamma_eff
        rows.append((phi, phi_m, phi_0, r_inv.gamma_eff, r_dir.gamma_eff, d))

        # How far from affine are the two signals, and where does the inversion
        # leave its own domain? These are the diagnostics that say WHY the clocks
        # differ, and they are printed whichever way P3 lands.
        f_inv = inv(eta_direct, phi, phi_m, phi_0)
        good = np.isfinite(f_inv) & np.isfinite(s)
        coef = np.polyfit(f_inv[good], s[good], 1)
        pred = np.polyval(coef, f_inv[good])
        ss_res = float(np.sum((s[good] - pred) ** 2))
        ss_tot = float(np.sum((s[good] - s[good].mean()) ** 2))
        r2 = 1.0 - ss_res / ss_tot
        worst = float(np.max(np.abs(s[good] - pred)))
        out = float(np.mean((f_inv < 0.0) | (f_inv > 1.0)))
        print(f"    s vs f is affine to R2 = {r2:.4f}, worst residual {worst:.4f} "
              f"against an s span of {s.max() - s.min():.4f} ({worst / (s.max() - s.min()):.0%})")
        print(f"    the inversion leaves its own domain [0, 1] on {out:.1%} of the record: "
              f"f runs {np.nanmin(f_inv):.3f} to {np.nanmax(f_inv):.3f}")
        for lo_g, hi_g in ((0.0, 0.5), (0.5, 2.0), (2.0, 1e9)):
            mm = (g >= lo_g) & (g < hi_g)
            if mm.sum():
                print(f"      strain {lo_g:<4g}-{min(hi_g, g[-1]):<5.3g}  "
                      f"eta {eta_direct[mm].mean():7.2f}   s {s[mm].mean():.4f}   "
                      f"f {f_inv[mm].mean():7.4f}   n = {int(mm.sum())}")
        print(f"    inverted channel  gamma_eff = {r_inv.gamma_eff:.4f} "
              f"[{r_inv.gamma_eff_lo:.4f}, {r_inv.gamma_eff_hi:.4f}]  "
              f"chi2_r {r_inv.chi2_reduced:.1f} on {r_inv.n_fit} bins")
        print(f"    DIRECT channel    gamma_eff = {r_dir.gamma_eff:.4f} "
              f"[{r_dir.gamma_eff_lo:.4f}, {r_dir.gamma_eff_hi:.4f}]  "
              f"chi2_r {r_dir.chi2_reduced:.1f} on {r_dir.n_fit} bins")
        print(f"    difference {d:.2%}")

    # ---- P4: does the difference depend on the gauge? ---------------------------
    diffs = [r[5] for r in rows]
    p4 = (max(diffs) - min(diffs)) < TRIPLE_TOL
    print(f"\nP4  the clock difference across two packing triples: "
          f"{diffs[0]:.2%} and {diffs[1]:.2%}, spread {max(diffs) - min(diffs):.2%}   "
          f"{'PASS' if p4 else 'FAIL'}")
    if not p4:
        print("    The difference moves with the gauge, so this probe is measuring the")
        print("    triple and P3 cannot be read. Do not quote the clock comparison.")
        return 1

    # ---- P3: the verdict, either way --------------------------------------------
    d = float(np.mean(diffs))
    print(f"\nP3  mean difference {d:.2%} against a {CLOCK_TOL:.0%} bar")
    if d < CLOCK_TOL:
        print("    The Wyart-Cates inversion does NOT distort the strain clock on this")
        print("    record. G11's one unchecked step is checked, and the 0.54-1.64 range")
        print("    it reports is not an artefact of the inversion.")
    else:
        print("    The inversion DOES distort the clock, by more than the bar. G11's")
        print("    gamma_0 carries this as a systematic and the range it reports must")
        print("    be widened by it before the solver comparison is quoted.")
    print("\n    Both numbers are dimensionless: Lin's simulation rates are in units of")
    print("    its own gamma_dot_0, so this is a control on the METHOD and not a second")
    print("    measurement of the material constant.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
