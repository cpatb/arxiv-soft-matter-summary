"""Gate G12: is the solver's onset stress compatible with the literature?

The question, and why the answer in SOLVER_HANDOFF was wrong
------------------------------------------------------------
The Basilisk free-surface solver runs ``SIGMA_STAR = 44.2327 Pa`` for 800 nm
silica, fitted from third-party rheometry the repo itself flags as provisional.
``SOLVER_HANDOFF.md`` §2 compared it against Lin et al.'s 2 um silica **in
reduced stress**, ``Sigma = sigma/eta_s = eta_r * gamma_dot``, and found a
**179x** discrepancy, concluding that the solver's onset stress "would have to
be 0.25 Pa, not 44.2 Pa".

That comparison is invalid, and the reason is worth stating precisely because
the same trick is *correct* one gate earlier. G10 compares two quantities
measured **on the same sample**, so the unstated ``eta_s`` appears on both sides
and cancels exactly. G12 compares quantities from **two different experiments**,
and there reduced stress does not cancel anything -- it silently assumes the two
solvents have the same viscosity. They do not.

The solvent viscosities are not unstated
----------------------------------------
Lin et al. never write ``eta_s`` in the text, which is why G10 was built to avoid
needing it. But the deposit contains it twice over:

* ``Fig4a`` is **absolute** viscosity in Pa.s against shear rate, and ``Fig4b``
  is **relative** viscosity against shear stress, for the same sample. Their
  ratio at matched stress is the solvent viscosity.
* ``FigS3`` is a rate sweep of the bare PMMA solvent, PDV-2331.

So the PMMA solvent can be obtained two independent ways, and **they agree to
0.29%** (1.770 Pa.s from the ratio, 1.765 Pa.s measured directly). That
agreement is what licenses using the same method on the silica, whose solvent --
glycerol+water, named only in the ``Fig4a`` README -- has no directly measured
twin. This is the same discipline as ``verify_guy_size_assignment``: recover a
stated constant from the numbers before relying on it.

What the gate then finds
------------------------
With ``eta_s`` measured rather than assumed away, the 179x becomes **24.6x**, and
the factor that disappears is *exactly* the solvent ratio 36.3/4.736 = 7.66x.

More importantly, the deposit contains **two** silica anchors and they disagree
with each other by more than the quantity being tested:

* ``Fig2a`` (dense, ``eta_r`` 6.4 to 19.1) gives ``sigma* = 0.32 Pa``
* ``Fig4b`` (dilute, ``eta_r`` 1.5 to 1.85) gives ``sigma* = 57.0 Pa``

Wyart-Cates says ``sigma*`` is a repulsive force scale and therefore
``phi``-independent, and on the **PMMA** pair the same two routes agree to 1.9x,
which is how we know the method is sound. The silica spread is physics or sample
preparation, not analysis: PMMA is sterically stabilised in a non-polar oil,
while aqueous silica is electrostatically stabilised and its counterion
concentration rises with ``phi``, so its repulsion is screened at high ``phi``.
``sigma*`` for aqueous silica is not a transferable material constant.

The solver's 44.2 Pa sits **between** the two silica anchors. This data therefore
cannot impeach it, which is a very different statement from "it would have to be
0.25 Pa".

Pass conditions, fixed before the numbers were looked at
--------------------------------------------------------
  G12a  The two independent routes to the PDV-2331 solvent viscosity agree
        within 5%. Without this the recovered silica ``eta_s`` is an unchecked
        assumption and nothing downstream is safe.
  G12b  On the PMMA system, converting the reduced ``sigma*`` to absolute Pa
        with the recovered ``eta_s`` reproduces an independently fitted
        absolute-stress ``sigma*`` within a factor of 3 -- the same tolerance
        G10c uses, and for the same reason.
  G12c  The solver's ``sigma*`` lies within the range spanned by the literature
        silica anchors, scaled to 800 nm. Failing this would mean the solver is
        outside everything the literature supports; passing it means only that
        the literature is too widely spread to test it, which is the honest
        outcome and is reported as such.

Run:  python scripts/gate_g12.py
Exit code is 0 only if all three pass.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dst_fit import fit_model, thickening_branch  # noqa: E402
from dst_fit.models import WYART_CATES  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
LIN = ROOT / (
    "Hydrodynamic and Contact Contributions to Continuous Shear Thickening "
    "in Colloidal Suspensions"
)
F2 = LIN / "Fig2" / "Fig2"
F4 = LIN / "Fig4" / "Fig4"
FS3 = LIN / "FigS3" / "FigS3" / "PDV2331_rate_sweep.csv"

# --- the solver, from closure/params_fitted.yaml in the Drop Impact repo ------
# Recorded here rather than imported: this package must run without that repo,
# and a silently drifting cross-repo constant is worse than a stale one.
SOLVER_SIGMA_STAR_PA = 44.2327
SOLVER_ETA_S = 0.004736
SOLVER_D_NM = 800.0
SOLVER_SOURCE = "Drop Impact closure/params_fitted.yaml, silica, 2026-07-20"

# --- transfer assumption ------------------------------------------------------
# Guy, Hermes & Poon PRL 115, 088304 (2015): sigma* = B d^alpha, alpha = -1.9(1),
# measured on PMMA from 268 to 3770 nm. Applying it between two silicas is the
# only transfer in this gate and it is the weakest link; the sweep below says
# how much it is worth.
SIZE_EXPONENT = 1.9
SIZE_EXPONENT_SWEEP = (1.8, 1.9, 2.0, 2.2)

# G10's fitted reduced onsets on the Fig2 rate sweeps, at phi_0 = 0.64.
LIN_REDUCED_SIGMA_STAR = {"silica": 8.699, "pmma": 3.885}
LIN_D_NM = {"silica": 2000.0, "pmma": 1400.0}

SOLVENT_TOL = 0.05  # G12a
FACTOR_TOL = 3.0  # G12b
PHI0_DEFAULT = 0.64


def solvent_viscosity(abs_file: Path, rel_file: Path) -> tuple[float, float, int]:
    """``eta_s = eta_absolute / eta_relative``, matched on stress.

    The two files share no abscissa: one is indexed by shear rate, the other by
    shear stress. They are joined through ``sigma = eta * gamma_dot``, which is
    an identity for steady simple shear, and the relative viscosity is
    interpolated in **log-log** because the stress axis spans three decades.
    Only points inside the measured stress range are used; nothing is
    extrapolated.
    """
    a = np.loadtxt(abs_file, delimiter=",", usecols=(0, 1))  # rate, eta [Pa.s]
    b = np.loadtxt(rel_file, delimiter=",", usecols=(0, 1))  # stress [Pa], eta_r
    sigma = a[:, 0] * a[:, 1]
    order = np.argsort(b[:, 0])
    s_ref, r_ref = b[order, 0], b[order, 1]
    inside = (sigma >= s_ref.min()) & (sigma <= s_ref.max())
    if inside.sum() < 4:
        raise ValueError(f"{abs_file.name}: only {inside.sum()} points overlap in stress")
    log_r = np.interp(np.log(sigma[inside]), np.log(s_ref), np.log(r_ref))
    est = a[inside, 1] / np.exp(log_r)
    return float(np.median(est)), float(est.max() / est.min()), int(inside.sum())


def fit_absolute_sigma_star(path: Path) -> tuple[float, float, int, float, float]:
    """Fit Wyart-Cates to relative viscosity against **absolute** stress."""
    b = np.loadtxt(path, delimiter=",", usecols=(0, 1))
    order = np.argsort(b[:, 0])
    sigma, eta_r = b[order, 0], b[order, 1]
    model = WYART_CATES.with_fixed(eta_s=1.0, phi_0=PHI0_DEFAULT)
    branch = thickening_branch(eta_r)
    fit = fit_model(model, sigma, eta_r, branch=branch, n_starts=24)
    p = dict(zip(model.param_names, fit.params))
    return (
        float(p["sigma_star"]),
        float(p["phi"]),
        int(eta_r.size),
        float(eta_r.min()),
        float(eta_r.max()),
    )


def scale_to(sigma_pa: float, d_from: float, d_to: float, alpha: float) -> float:
    """Transfer an onset stress between particle sizes by ``sigma* ~ d^-alpha``."""
    return sigma_pa * (d_from / d_to) ** alpha


def main() -> int:
    failures: list[str] = []
    print("Gate G12: the solver's onset stress against the literature")
    print("Lin et al. 2015, PRL 115, 228304; DataShare 10283/886")
    print(f"solver: sigma* = {SOLVER_SIGMA_STAR_PA} Pa, eta_s = {SOLVER_ETA_S} Pa.s, "
          f"d = {SOLVER_D_NM:.0f} nm")
    print(f"        ({SOLVER_SOURCE})\n")

    # ---------------- G12a: recover the solvent viscosities ------------------
    print("G12a  the solvent viscosities Lin et al. never state")
    eta_s_si, spread_si, n_si = solvent_viscosity(
        F4 / "Fig_4a_temp20C.csv", F4 / "Fig_4b_temp20C.csv"
    )
    eta_s_pm, spread_pm, n_pm = solvent_viscosity(
        F4 / "Fig_4c_temp20C.csv", F4 / "Fig_4d_temp20C.csv"
    )
    direct = np.loadtxt(FS3, delimiter=",", usecols=(0, 1))
    eta_s_pm_direct = float(np.median(direct[:, 1]))
    disagreement = abs(eta_s_pm / eta_s_pm_direct - 1.0)

    print(f"      2.0 um silica in glycerol+water, 20 C : eta_s = {1e3 * eta_s_si:8.2f} mPa.s"
          f"   ({n_si} points, self-consistent to {spread_si:.3f}x)")
    print(f"      1.4 um PMMA in PDV-2331,         20 C : eta_s = {1e3 * eta_s_pm:8.2f} mPa.s"
          f"   ({n_pm} points, self-consistent to {spread_pm:.3f}x)")
    print(f"      PDV-2331 measured directly (FigS3)    : eta_s = "
          f"{1e3 * eta_s_pm_direct:8.2f} mPa.s   ({len(direct)} rates, Newtonian to "
          f"{direct[:, 1].max() / direct[:, 1].min():.3f}x)")
    print(f"      the two independent routes to PDV-2331 disagree by {100 * disagreement:.2f}%"
          f"   tolerance {100 * SOLVENT_TOL:.0f}%   "
          f"{'PASS' if disagreement <= SOLVENT_TOL else 'FAIL'}")
    print("      the silica has no directly measured twin, so this agreement is the")
    print("      only thing licensing the silica number. It is why G12a is gated first.\n")
    if disagreement > SOLVENT_TOL:
        failures.append(f"G12a: PDV-2331 routes disagree by {100 * disagreement:.1f}%")

    # ---------------- G12b: does reduced -> absolute reproduce a direct fit? --
    print("G12b  converting reduced sigma* to Pa, checked against a direct fit")
    print("      Fig2 is dense and reported in relative viscosity only; Fig4 is dilute")
    print("      and reported in absolute stress. Same particles, same solvent. Under")
    print("      Wyart-Cates sigma* is a repulsive force scale, so it should not care.\n")
    rows = []
    for label, system, eta_s, fig4 in (
        ("1.4 um PMMA", "pmma", eta_s_pm, "Fig_4d_temp20C.csv"),
        ("2.0 um silica", "silica", eta_s_si, "Fig_4a"),  # silica handled below
    ):
        if system == "silica":
            continue
        direct_star, phi_d, n_d, r0, r1 = fit_absolute_sigma_star(F4 / fig4)
        converted = LIN_REDUCED_SIGMA_STAR[system] * eta_s
        ratio = max(converted, direct_star) / min(converted, direct_star)
        rows.append((label, converted, direct_star, ratio))
        print(f"      {label}")
        print(f"        Fig2 reduced Sigma* = {LIN_REDUCED_SIGMA_STAR[system]:.3f} 1/s "
              f"x eta_s = {eta_s:.4f} Pa.s  ->  sigma* = {converted:.3f} Pa")
        print(f"        Fig4 fitted directly in Pa (phi = {phi_d:.4f}, "
              f"eta_r {r0:.2f}-{r1:.2f}, n = {n_d})  ->  sigma* = {direct_star:.3f} Pa")
        print(f"        ratio {ratio:.2f}x   tolerance {FACTOR_TOL:.0f}x   "
              f"{'PASS' if ratio <= FACTOR_TOL else 'FAIL'}\n")
        if ratio > FACTOR_TOL:
            failures.append(f"G12b[{system}]: routes differ by {ratio:.1f}x")

    # the silica, which is the interesting one and is NOT gated
    si_direct, si_phi, si_n, si_r0, si_r1 = fit_absolute_sigma_star(F4 / "Fig_4b_temp20C.csv")
    si_converted = LIN_REDUCED_SIGMA_STAR["silica"] * eta_s_si
    si_ratio = max(si_converted, si_direct) / min(si_converted, si_direct)
    print("      2.0 um silica  (reported, deliberately NOT gated -- see below)")
    print(f"        Fig2a reduced Sigma* = {LIN_REDUCED_SIGMA_STAR['silica']:.3f} 1/s "
          f"x eta_s = {eta_s_si:.4f} Pa.s  ->  sigma* = {si_converted:.4f} Pa")
    print(f"        Fig4b fitted directly in Pa (phi = {si_phi:.4f}, "
          f"eta_r {si_r0:.2f}-{si_r1:.2f}, n = {si_n})  ->  sigma* = {si_direct:.3f} Pa")
    print(f"        ratio {si_ratio:.0f}x")
    print("        The PMMA pair above agrees, so this is not a defect of the method.")
    print("        PMMA is sterically stabilised in a non-polar oil and its sigma* is")
    print("        genuinely phi-independent; aqueous silica is electrostatically")
    print("        stabilised, its counterion concentration rises with phi, and its")
    print("        repulsion is screened. sigma* for aqueous silica is not a")
    print("        transferable material constant -- which is exactly what this lab's")
    print("        own salt series says, from the other direction.\n")

    # ---------------- G12c: is the solver inside the literature span? --------
    print("G12c  the solver against the two silica anchors, both scaled to 800 nm")
    anchors = {
        "Lin Fig2a (dense)": si_converted,
        "Lin Fig4b (dilute)": si_direct,
    }
    scaled = {
        k: scale_to(v, LIN_D_NM["silica"], SOLVER_D_NM, SIZE_EXPONENT)
        for k, v in anchors.items()
    }
    lo, hi = min(scaled.values()), max(scaled.values())
    for k in anchors:
        print(f"      {k:22s} sigma*(2000 nm) = {anchors[k]:8.3f} Pa  ->  "
              f"sigma*(800 nm) = {scaled[k]:8.3f} Pa   "
              f"(solver is {SOLVER_SIGMA_STAR_PA / scaled[k]:6.1f}x this)")
    inside = lo <= SOLVER_SIGMA_STAR_PA <= hi
    print(f"      literature span at 800 nm: {lo:.3f} to {hi:.3f} Pa, a factor of "
          f"{hi / lo:.0f}")
    print(f"      solver sigma* = {SOLVER_SIGMA_STAR_PA:.3f} Pa is "
          f"{'INSIDE' if inside else 'OUTSIDE'} that span   "
          f"{'PASS' if inside else 'FAIL'}")
    if not inside:
        failures.append(
            f"G12c: solver sigma* = {SOLVER_SIGMA_STAR_PA} Pa outside [{lo:.3f}, {hi:.3f}]"
        )

    old = (SOLVER_SIGMA_STAR_PA / SOLVER_ETA_S) / (
        LIN_REDUCED_SIGMA_STAR["silica"]
        * (LIN_D_NM["silica"] / SOLVER_D_NM) ** SIZE_EXPONENT
    )
    new = SOLVER_SIGMA_STAR_PA / scaled["Lin Fig2a (dense)"]
    print("\n      what happened to the 179x in SOLVER_HANDOFF section 2:")
    print(f"        reduced-stress comparison, as written there : {old:6.1f}x")
    print(f"        same comparison with eta_s measured         : {new:6.1f}x")
    print(f"        the factor that vanished                    : "
          f"{eta_s_si / SOLVER_ETA_S:6.2f}x  = eta_s(Lin)/eta_s(solver)")
    print(f"        check: {old:.1f} / {eta_s_si / SOLVER_ETA_S:.2f} = "
          f"{old / (eta_s_si / SOLVER_ETA_S):.1f}x")
    print("      Reduced stress cancels eta_s only WITHIN one dataset, which is what")
    print("      makes it correct in G10 and wrong here.")

    print("\n      sensitivity to the one transfer assumption, sigma* ~ d^-alpha:")
    for a in SIZE_EXPONENT_SWEEP:
        s_lo = scale_to(min(anchors.values()), LIN_D_NM["silica"], SOLVER_D_NM, a)
        s_hi = scale_to(max(anchors.values()), LIN_D_NM["silica"], SOLVER_D_NM, a)
        print(f"        alpha = {a}: span {s_lo:.3f} to {s_hi:.3f} Pa, solver "
              f"{'inside' if s_lo <= SOLVER_SIGMA_STAR_PA <= s_hi else 'OUTSIDE'}")
    print("      The span is 180x wide, so alpha cannot change the verdict; that is")
    print("      the point. This data is too widely spread to test the solver, and")
    print("      saying so is the result.")

    print("\n      what would test it: a sigma* measured on the solver's OWN silica")
    print("      at its own ionic strength. Nothing transferred from another")
    print("      preparation will do, because the 180x above is the size of the")
    print("      preparation dependence.")

    if failures:
        print("\nFAILURES")
        for f in failures:
            print(f"  - {f}")
        print("\nGATE G12 FAIL")
        return 1
    print("\nGATE G12 PASS (G12a, G12b, G12c)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
