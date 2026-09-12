"""analyse_phase8.py — extract phi_c(CV) from the constant-volume shear ladders.

FOUR ESTIMATORS, AND WHY THE PRIMARY IS THE ONE IT IS
----------------------------------------------------
At constant volume every cell reaches a steady state that lies on the granular
`mu(I)` / `phi(I)` rheology, so the ladder is not a set of independent points: it is a
sampling of one curve, parameterised by the inertial number `I = gammadot d50
sqrt(rho/p)`. The critical-state fraction is the quasi-static limit of that curve.

  E2  PRESSURE ROOT, free exponent  [PRIMARY]   p_ss = c (phi - phi_c)^b, fitted in log
      space. The critical state is where the constant-volume grain pressure VANISHES, so
      extrapolating p -> 0 targets it directly. `b` is reported as a diagnostic of the
      contact law, never used to redefine phi_c.
  E3  PRESSURE ROOT, b pinned to 1      the harmonic-contact convention, as a check.
  E4  COORDINATION CROSSING             the phi at which z_nr crosses Z_ISO_FRICTIONAL.
      Independent of the stress entirely, which is its whole value: a defect in the
      stress read-out cannot move it.

  E1  phi(I) EXTRAPOLATION   🔴 **RETRACTED 2026-09-08 AS A phi_c ESTIMATOR**, having
      been the PRIMARY one. It fits phi against the inertial number and extrapolates to
      I = 0, which is the standard construction WHEN the pressure is imposed and the
      shear rate is varied. This protocol does the opposite: `gdot` is IDENTICAL on every
      rung of a ladder (measured: 3.000e-05 across all five at cv = 0, varying by
      1.0000x), so I = gdot d sqrt(rho/p) moves ONLY through p. Therefore

          I -> 0  is the HIGH-PRESSURE, DENSEST end of the ladder,
                  NOT the critical state.

      The extrapolation ran away from phi_c, not toward it. That is why E1 sat 0.074
      ABOVE E2 -- it was reporting the dense end while E2 reported the vanishing-pressure
      end, and the "systematic between estimators" was two different quantities. E1 is
      still computed and printed, relabelled `phi_dense_E1`, because the dense-end value
      is a real number and its CV contrast happens to agree with E2's; it is no longer
      offered as phi_c.

      ⚠ The CV CONTRAST is barely affected (+0.0578 from E1 vs +0.0623 from E2), because
      both ends of the ladder shift together with CV. P1 and P2 therefore survive the
      retraction; what does not survive is any ABSOLUTE phi_c quoted from E1.

The four are reported side by side for every cell group. A CV contrast is quoted only
if it has the same sign under all four; the SPREAD between them is the systematic and
it is printed.

WHAT IS REFUSED
---------------
  * a root outside the ladder it was fitted on (extrapolation dressed as measurement)
  * fitting against `phi_target`; every fit uses the REALISED `phi`, because stage B
    lands on the target only to a tolerance and the first version of it did not

THE QUALITY GATES, AND WHY EVERY CONTRAST IS PRINTED FOUR TIMES
---------------------------------------------------------------
🔴 Until 2026-09-09 this docstring claimed "a group whose cells are not all steady, or
not all quasi-static -- excluded cells are PRINTED with their reason, never dropped
silently". That was FALSE for the primary estimator. `e1_phi_of_I` does filter on
`I <= I_FIT_MAX`, but `_pressure_root` (E2 and E3) and `e4_z_crossing` filtered on
NOTHING except a pressure floor. The flags were counted into `n_steady`/`n_quasistatic`,
printed, and never used to exclude. Found by a numerics reviewer 2026-09-09.

The fix is NOT to switch the gates on, because there is no defensible single answer:

  * the QUASI-STATIC gate removes, at cv = 0 and 0.2, exactly and only the BOTTOM rung
    of each ladder -- the rung NEAREST the root the fit is extrapolating to. Gating
    therefore buys admissibility at the cost of a longer extrapolation. Both costs are
    real and they point opposite ways.
  * the exclusions are CV-ASYMMETRIC (3/3/1 non-quasi-static cells at cv = 0/0.2/0.5),
    which is the direction that inflates the contrast.

So `main` computes every contrast under ALL FOUR gate sets and prints them side by side,
with the per-CV-level `n` beside each, and never reports one alone. Measured on prod2
as of 2026-09-09 (the numbers this module prints, not typed here from memory -- see
out/phase8_verdict.txt):

        ungated                         dphi_c = +0.0623
        quasi-static cells only         dphi_c = +0.0444
        steady cells only               dphi_c = +0.0485
        both gates                      dphi_c = +0.0438

**The headline falls ~30 % under any gate.** Sign and qualitative conclusion survive;
the number does not. **Quote the range +0.044 to +0.062, with the gate named.**

`test_shear.py::test_every_defined_gate_is_applied` fails if a gate is added to `GATES`
and no estimator path honours it.

THE SHEAR RATE IS NOT A CONTROLLED VARIABLE ALONG THE CV LADDER
---------------------------------------------------------------
🔴 Found 2026-09-09, previously unrecorded. `shear.py` sets
`gdot = gamma_max / (shear_steps * O.dt)` with `O.dt = DT_SAFETY * PWaveTimeStep()`,
and the P-wave timestep is set by the SMALLEST particle. Widening the size distribution
shrinks the smallest particle, shrinks `dt`, and therefore RAISES `gdot`. Measured on
prod2: `gdot` = 3.000e-05 at cv = 0 against 1.437e-04 at cv = 0.5, a factor **4.79**.

Every CV contrast on prod2 is therefore confounded with a shear-rate contrast, and
`I = gdot d50 sqrt(rho/p)` carries it directly into the quasi-static gate above. This
module now prints the realised `gdot` per CV level beside every contrast and raises a
🔴 line when it varies by more than `GDOT_MATCH_TOL` across the levels being contrasted.
The control run that removes the confound is tagged `ratem` (analyse it with
`--pattern "phase8_ratem*.jsonl"`), where `shear_steps` was chosen per cell so that
`gdot` lands on the cv = 0.5 value.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "out")

# --- named constants, with origin ---------------------------------------------------
I_FIT_MAX = 1.0e-2      # quasi-static bar; same value as shear.I_MAX, kept in sync by
                        # test_shear.py::test_i_max_agrees
# 🔴 The isostatic contact number is NOT one constant: it is 2d = 6 for FRICTIONLESS
# spheres and d + 1 = 4 for fully frictional ones. Hardcoding 4 was right for the mu_p =
# 0.3 production arm and WRONG for the mu_p = 0 anchor run, where the pack is isostatic
# at 6 and E4 would have hunted for a crossing at 4 that happens far below jamming --
# a constant that is correct for one arm and silently wrong for another, which is the
# trap this project keeps meeting. Chosen per cell from mu_p.
Z_ISO_FRICTIONLESS = 6.0
Z_ISO_FRICTIONAL = 4.0  # the project's A2/A3 anchors measure 6.37-6.92 at mu_p = 0 and
                        # 3.98-4.07 at mu_p = 3, which is what pins both values.
Z_ISO_MU_BAR = 1.0e-9   # below this mu_p the pack is treated as frictionless

# 🔴 AND THE FRICTIONAL CONSTANT ABOVE IS THE WRONG NUMBER FOR THE PRODUCTION ARM.
# d+1 = 4 is the mu_p -> INFINITY limit; at finite friction nothing fixes the count and
# it has to be measured. `measure_z_iso.py` extrapolates z_nr to zero confining pressure
# on the project's OWN isotropic ladders (phase2_*.jsonl, the same records and the same
# construction phase2_phiJ.py uses for phi_J) and gets, at N = 1500:
#
#     mu_p = 0.00   5.859 +- 0.078 (11 cells)   against 6.0 used -- theory-exact, fine
#     mu_p = 0.30   4.535 +- 0.072 (18 cells)   against 4.0 used -- 0.535 CONTACTS OUT
#     mu_p = 3.00   3.961          ( 1 cell )   against 4.0 used -- fine
#
# At the shear ladders' own dz/dphi (28.1 / 23.0 / 21.9 per unit phi at cv 0 / 0.2 / 0.5)
# that is 0.019-0.024 in phi_c, roughly half the CV contrast E4 is used to measure. And
# z_iso DRIFTS with CV (lognormal 4.609 -> 4.377 over CV 0 -> 0.5), so one constant puts
# a CV-dependent bias into E4's own contrast, worth about +0.0045 in the inflating
# direction.
#
# ⚠ Using these values makes E4 depend on the ISOTROPIC measurement, which is also what
# the shear ladders are centred on -- so a corrected E4 is a better estimator of phi_c
# and is NOT an independent check of "phi_c = phi_J". Both readings are reported; the
# default stays the theory constant so that no previously published number moves
# silently. Source: out/z_iso_measured.txt.
#
# ✅ AUDITED 2026-09-12 -- PREREG_ZISO_2026-09-12.md, out/ziso_audit.txt, 8/8 gates.
# The table below is measured under the SHIPPED rung-selection rule, which was found to
# fit two distinct pressures while claiming five rungs (D2). The values are DELIBERATELY
# LEFT as they are: the rung-correct rule moves them by at most 0.052 contacts and moves
# E4m's CV contrast by +0.0003 (0.0497 -> 0.0501), so changing them here would invalidate
# the reproduction of every committed artifact to buy the third decimal place. The
# correction is measured and bounded, which is worth more.
# Two consequences worth carrying:
#   - Delta phi_c is NOT sensitive to this target. A uniform shift cancels out of a
#     contrast exactly; only the CV-dependence survives, and that is 0.0003 in phi. So
#     the E4/E4m split on H-3 (0.0069 vs 0.0126) CANNOT be a target effect of this size.
#   - the CV drift is a LOGNORMAL effect (-0.213); on the zero-skew families it is
#     -0.056 (uniform) and -0.012 (bidisperse). `test_shear.py` recomputes this table
#     with rule="broken" on purpose.
Z_ISO_MEASURED = {
    ("lognormal", 0.0): 4.609, ("lognormal", 0.1): 4.553, ("lognormal", 0.2): 4.513,
    ("lognormal", 0.3): 4.431, ("lognormal", 0.4): 4.446, ("lognormal", 0.5): 4.377,
    ("uniform", 0.0): 4.609, ("uniform", 0.1): 4.612, ("uniform", 0.2): 4.563,
    ("uniform", 0.3): 4.547, ("uniform", 0.4): 4.550, ("uniform", 0.5): 4.416,
    ("bidisperse", 0.0): 4.609, ("bidisperse", 0.1): 4.595, ("bidisperse", 0.2): 4.535,
    ("bidisperse", 0.3): 4.552, ("bidisperse", 0.4): 4.564, ("bidisperse", 0.5): 4.542,
}
Z_ISO_MEASURED_MU = 0.3         # the friction those numbers were measured at


def z_iso_measured(family: str, cv: float, mu_p: float) -> float:
    """The MEASURED isostatic count for this cell, or NaN if none was measured.

    Only mu_p = 0.3 is tabulated, because that is the only friction at which the
    theory constant is wrong and an isotropic ladder exists to replace it.
    """
    if abs(mu_p - Z_ISO_MEASURED_MU) > 1e-9:
        return float("nan")
    return Z_ISO_MEASURED.get((family, round(float(cv), 3)), float("nan"))
P_FLOOR = 1.0e-9        # below this p_hat the cell carries no measurable stress
B_MAX = 8.0             # upper bound on E2's free exponent. A fit landing ON it is
                        # reported as bound-limited and its phi_c is refused.
PHI_MATCH_TOL = 5.0e-3  # how close two rungs from different CV ladders must be to count
                        # as "matched phi". The ladder spacing is 0.020, so this admits
                        # only genuine near-coincidences and never two different rungs.
PHI_GROUP_TOL = 1.0e-3  # how close two SEEDS' realised phi must be to be one rung;
                        # stage B lands within ~1e-4, so this is loose by 10x
GDOT_MATCH_TOL = 1.20   # how far the realised shear rate may vary across the CV levels
                        # of a contrast before the contrast is flagged as confounded
                        # with a rate contrast. The three SEEDS at one CV span at
                        # most 1.10 (dt is set by the smallest particle DRAWN), so 1.20
                        # sits just above seed scatter and any ladder-wide move beyond
                        # it is a systematic. Measured span on prod2: 4.78x.

# --- the quality gates -------------------------------------------------------------
# Every gate is a named predicate on ONE record. `analyse_group` takes a tuple of gate
# names and drops every record failing any of them BEFORE any estimator sees the data,
# and `main` runs every set in GATE_SETS and prints them side by side. Adding an entry
# here without an estimator honouring it is what went wrong before 2026-09-09; the test
# `test_shear.py::test_every_defined_gate_is_applied` now fails if that happens again.
GATES = {
    "quasistatic": lambda r: bool(r.get("quasistatic", False)),
    "steady": lambda r: bool(r.get("steady", False)),
}
GATE_SETS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("ungated", ()),
    ("quasistatic", ("quasistatic",)),
    ("steady", ("steady",)),
    ("both", ("quasistatic", "steady")),
)


def apply_gates(recs: list[dict], gate_names: tuple[str, ...]) -> list[dict]:
    """Records passing EVERY named gate. Unknown names raise rather than pass silently."""
    for g in gate_names:
        if g not in GATES:
            raise KeyError(f"unknown gate {g!r}; defined: {sorted(GATES)}")
    return [r for r in recs if all(GATES[g](r) for g in gate_names)]


def load(pattern: str) -> list[dict]:
    """Records under out/ matching `pattern`, which may be a COMMA-SEPARATED list.

    The list form exists for the rate-matched control, whose arms were run with
    different `shear_steps` per CV level and therefore live under different tags; a
    single glob cannot name them. Every record still carries its own `protocol`,
    `gamma_max`, `shear_steps` and realised `gdot`, all of which `group_key` and the
    contrast block report, so pooling two patterns cannot hide a knob.

    ⚠ Records are DEDUPLICATED on their `key`, which encodes every swept knob including
    the realised-target `phi`. Two patterns can legitimately name the same cell -- the
    down-shifted ladder and the production one share their top two rungs at cv = 0 by
    construction -- and without this the duplicate-state guard fires 18 times on a
    perfectly sound selection. Every dropped duplicate is CHECKED against the record it
    duplicates on `p_ss_hat` and `z_nr_ss`; a disagreement is a REPRODUCIBILITY failure
    and is printed rather than swallowed.
    """
    recs: list[dict] = []
    seen: dict[str, dict] = {}
    paths: list[str] = []
    for pat in (p.strip() for p in pattern.split(",") if p.strip()):
        paths.extend(glob.glob(os.path.join(OUT_DIR, pat)))
    for path in sorted(set(paths)):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                k = r.get("key")
                if isinstance(k, str) and k in seen:
                    prev = seen[k]
                    for fld in ("p_ss_hat", "z_nr_ss"):
                        a, b = prev.get(fld), r.get(fld)
                        if (isinstance(a, float) and isinstance(b, float)
                                and not np.isclose(a, b, rtol=1e-9, atol=0.0)):
                            print(f"🔴 REPRODUCIBILITY: {k} appears twice with "
                                  f"{fld} = {a!r} and {b!r}")
                    continue
                if isinstance(k, str):
                    seen[k] = r
                recs.append(r)
    return recs


def group_key(r: dict) -> tuple:
    """Everything a comparison must hold fixed, in the key.

    🔴 `protocol`, `gamma_max` and `shear_steps` are in here deliberately. The default
    glob `phase8_prod*.jsonl` matches BOTH `phase8_prod*` (protocol cv-shear-v1) and
    `phase8_prod2*` (v2, after the periodic-cell flip fix and with a different ladder,
    strain and step count). Pooling those would be the N-pooling defect of 2026-09-08
    repeated in a new place: a key that does not encode a knob the analysis groups by.
    `main` additionally REFUSES to report contrasts when more than one protocol is
    present, rather than quietly averaging over it.
    """
    return (r.get("family", "lognormal"), round(float(r["cv"]), 3),
            round(float(r["mu_p"]), 3), int(r["seed"]), int(r["n_requested"]),
            r.get("protocol", "?"), round(float(r.get("gamma_max", float("nan"))), 3),
            int(r.get("shear_steps", -1)))


def e1_phi_of_I(phi: np.ndarray, ivals: np.ndarray) -> dict:
    """phi at the DENSE end, lim_{I->0} phi(I), over the quasi-static rungs.

    🔴 NOT phi_c. See the module docstring: on this protocol gdot is fixed, so I moves
    only through p and I -> 0 is the densest rung. Kept and reported because the value
    is real and its CV contrast is informative; renamed so it cannot be read as phi_c.
    """
    ok = np.isfinite(phi) & np.isfinite(ivals) & (ivals <= I_FIT_MAX) & (ivals > 0)
    out = {"phi_dense_E1": np.nan, "slope_E1": np.nan, "n_E1": int(ok.sum()),
           "r2_E1": np.nan, "I_min": np.nan, "I_max": np.nan}
    if ok.sum() < 3:
        return out
    x, y = ivals[ok], phi[ok]
    A = np.vstack([np.ones_like(x), x]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ coef
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    out.update(phi_dense_E1=float(coef[0]), slope_E1=float(coef[1]),
               r2_E1=float(1 - np.sum((y - pred) ** 2) / ss_tot) if ss_tot > 0 else np.nan,
               I_min=float(x.min()), I_max=float(x.max()))
    return out


def _pressure_root(phi: np.ndarray, p: np.ndarray, pin_b: float | None) -> dict:
    """Root of p(phi) = c (phi - phi_c)^b, fitted on log p."""
    from scipy.optimize import curve_fit

    tag = "E3" if pin_b is not None else "E2"
    out = {f"phi_c_{tag}": np.nan, f"b_{tag}": (pin_b if pin_b else np.nan),
           f"n_{tag}": 0}
    ok = np.isfinite(phi) & np.isfinite(p) & (p > P_FLOOR)
    if ok.sum() < (3 if pin_b is not None else 4):
        return out
    x, y = phi[ok], np.log(p[ok])
    lo = float(x.min())

    if pin_b is None:
        def model(xx, phic, logc, b):
            return logc + b * np.log(np.maximum(xx - phic, 1e-12))
        p0 = [lo - 0.02, float(y.mean()), 1.0]
        bounds = ([-np.inf, -np.inf, 0.1], [lo - 1e-6, np.inf, B_MAX])
    else:
        def model(xx, phic, logc):
            return logc + pin_b * np.log(np.maximum(xx - phic, 1e-12))
        p0 = [lo - 0.02, float(y.mean())]
        bounds = ([-np.inf, -np.inf], [lo - 1e-6, np.inf])
    try:
        popt, _ = curve_fit(model, x, y, p0=p0, bounds=bounds, maxfev=40000)
    except Exception as exc:                     # noqa: BLE001 - report, do not hide
        out[f"fit_error_{tag}"] = repr(exc)
        return out
    out[f"phi_c_{tag}"] = float(popt[0])
    out[f"n_{tag}"] = int(ok.sum())
    if pin_b is None:
        out["b_E2"] = float(popt[2])

    # 🔴 A FIT ON A BOUND IS NOT A MEASUREMENT. Both roots are bounded above by
    # `lo - 1e-6` (the root must lie below the lowest rung). On the 2026-09-08 prod2 run
    # E3 returned EXACTLY that bound on 6 of 9 groups -- 0.5745 at cv 0 and 0.5880 at
    # cv 0.2, i.e. the ladder's own bottom rung, with a seed sd of 0.0000 and an
    # apparent significance of 161 sd. That is the bound reporting itself, not phi_c.
    # Flagged and returned as NaN, with the value kept under a `_at_bound` key.
    if abs(popt[0] - (lo - 1e-6)) < 1e-4:
        out[f"phi_c_{tag}_at_bound"] = float(popt[0])
        out[f"phi_c_{tag}"] = np.nan

    # 🔴 AND THE SAME RULE FOR THE EXPONENT. E2's `b` is bounded to [0.1, 8.0]; a fit
    # that lands ON that bound has not found a power law, it has run out of room, and
    # its phi_c is then whatever the bound forced. Measured on the 2026-09-08 P4
    # frictionless anchor: b = 8.00, 7.33, 8.00 on the three seeds, and the phi_c it
    # returned (0.580) missed the anchor by 0.052 while the isostatic crossing hit it to
    # 0.0015. Frictionless jamming is genuinely much sharper than a low-order power law,
    # so this is the fit reporting its own inadequacy -- which is worth keeping, not
    # hiding.
    if pin_b is None and popt[2] >= B_MAX - 1e-6:
        out["b_E2_at_bound"] = float(popt[2])
        out["phi_c_E2_at_bound"] = float(popt[0])
        out["phi_c_E2"] = np.nan
    return out


def e4_z_crossing(phi: np.ndarray, z: np.ndarray, mu_p: float = 1.0,
                  z_iso: float | None = None, tag: str = "E4") -> dict:
    """phi where z_nr crosses the isostatic count, by linear interpolation.

    The target depends on friction: 6 for frictionless spheres, 4 for frictional ones.
    Passing `mu_p` is mandatory in spirit -- the default is the frictional case only so
    that older callers keep their previous behaviour rather than silently changing it.

    `z_iso` overrides the theory constant with a MEASURED count (see Z_ISO_MEASURED and
    measure_z_iso.py). It is passed, not looked up here, so that a caller cannot get a
    measured target by accident.
    """
    if z_iso is None:
        z_iso = (Z_ISO_FRICTIONLESS if mu_p <= Z_ISO_MU_BAR else Z_ISO_FRICTIONAL)
    out = {f"phi_c_{tag}": np.nan, f"z_iso_used_{tag}" if tag != "E4" else "z_iso_used":
           z_iso}
    if not np.isfinite(z_iso):
        return out
    ok = np.isfinite(phi) & np.isfinite(z)
    if ok.sum() < 2:
        return out
    o = np.argsort(phi[ok])
    x, y = phi[ok][o], z[ok][o]
    if y.min() > z_iso or y.max() < z_iso:
        return out
    for i in range(len(x) - 1):
        if (y[i] - z_iso) * (y[i + 1] - z_iso) <= 0:
            if y[i + 1] == y[i]:
                continue
            t = (z_iso - y[i]) / (y[i + 1] - y[i])
            out[f"phi_c_{tag}"] = float(x[i] + t * (x[i + 1] - x[i]))
            break
    return out


def analyse_group(recs: list[dict], gate_names: tuple[str, ...] = ()) -> dict:
    """Every estimator for one (family, cv, mu, seed, ...) ladder.

    `gate_names` is applied to the RECORDS, before any estimator sees them, so that
    every estimator honours the same gate. The unfiltered counts are still reported as
    `n_steady`/`n_quasistatic` so a reader can see what a gate would cost.
    """
    n_before = len(recs)
    recs = apply_gates(recs, gate_names)
    recs = sorted(recs, key=lambda r: r["phi"])
    if not recs:
        return {"n_cells": 0, "n_dropped_by_gate": n_before,
                "gate": "+".join(gate_names) or "ungated",
                "phi_lo": np.nan, "phi_hi": np.nan, "n_steady": 0, "n_quasistatic": 0,
                "p_hat_max": np.nan, "phi_c_E2": np.nan, "phi_c_E3": np.nan,
                "phi_c_E4": np.nan, "phi_c_E4m": np.nan, "phi_dense_E1": np.nan,
                "mu_s": np.nan, "b_E2": np.nan, "gdot_mean": np.nan}
    phi = np.array([r["phi"] for r in recs], float)
    p = np.array([r["p_ss_hat"] for r in recs], float)
    ivals = np.array([r["i_ss"] for r in recs], float)
    z = np.array([r["z_nr_ss"] for r in recs], float)
    mu = np.array([r["mu_macro"] for r in recs], float)

    res = {"n_cells": len(recs), "phi_lo": float(phi.min()), "phi_hi": float(phi.max()),
           "n_steady": int(sum(bool(r["steady"]) for r in recs)),
           "n_quasistatic": int(sum(bool(r["quasistatic"]) for r in recs)),
           "n_dropped_by_gate": n_before - len(recs),
           "gate": "+".join(gate_names) or "ungated"}
    res["p_hat_max"] = float(np.nanmax(p))
    gd = np.array([r.get("gdot", np.nan) for r in recs], float)
    res["gdot_mean"] = float(np.nanmean(gd)) if np.isfinite(gd).any() else np.nan
    res.update(e1_phi_of_I(phi, ivals))
    res.update(_pressure_root(phi, p, None))
    res.update(_pressure_root(phi, p, 1.0))
    mu_cell = float(recs[0].get("mu_p", 1.0))
    res.update(e4_z_crossing(phi, z, mu_p=mu_cell))
    # the same crossing against the MEASURED isostatic count, reported beside it and
    # never instead of it -- see Z_ISO_MEASURED for why neither alone is the answer
    res.update(e4_z_crossing(
        phi, z, mu_p=mu_cell, tag="E4m",
        z_iso=z_iso_measured(str(recs[0].get("family", "lognormal")),
                             float(recs[0].get("cv", np.nan)), mu_cell)))

    # mu_s = mu_macro(I -> 0): an anchor on the rheology, independent of phi_c
    ok = np.isfinite(mu) & np.isfinite(ivals) & (ivals <= I_FIT_MAX) & (ivals > 0)
    if ok.sum() >= 3:
        A = np.vstack([np.ones_like(ivals[ok]), ivals[ok]]).T
        coef, *_ = np.linalg.lstsq(A, mu[ok], rcond=None)
        res["mu_s"] = float(coef[0])
    else:
        res["mu_s"] = np.nan

    # REFUSE a root outside the ladder it was fitted on.
    for tag in ("E2", "E3", "E4", "E4m"):
        v = res.get(f"phi_c_{tag}", np.nan)
        if np.isfinite(v) and not (res["phi_lo"] - 0.05 <= v <= res["phi_hi"]):
            res[f"phi_c_{tag}_refused"] = float(v)
            res[f"phi_c_{tag}"] = np.nan
    return res


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pattern", default="phase8_prod2*.jsonl",
                    help="glob under out/. NOTE the default is prod2, not prod*: "
                         "prod* would also match the protocol-v1 records")
    ap.add_argument("--out", default="phase8_verdict.txt")
    args = ap.parse_args()

    recs = load(args.pattern)
    errs = [r for r in recs if "error" in r]
    good = [r for r in recs if "error" not in r and np.isfinite(r.get("p_ss_hat", np.nan))]
    lines: list[str] = []

    def p(m: str = "") -> None:
        print(m)
        lines.append(m)

    p(f"phase8: {len(recs)} records from {args.pattern}; "
      f"{len(good)} usable, {len(errs)} errored")
    for r in errs:
        p(f"  ERRORED {r.get('key')}: {r.get('error')}")
    if not good:
        p("nothing to analyse")
        return 0

    # duplicate-state guard: two rungs landing on the same realised phi is the defect
    # that the 500-step compression chunk produced on 2026-09-08.
    by_g = defaultdict(list)
    for r in good:
        by_g[group_key(r)].append(r)
    n_dupe = 0
    for k, rs in by_g.items():
        phis = sorted(round(r["phi"], 6) for r in rs)
        for a, b in zip(phis, phis[1:]):
            if abs(a - b) < 1e-4:
                p(f"  🔴 DUPLICATE STATE in {k}: two rungs at phi = {a:.6f}, {b:.6f}")
                n_dupe += 1
    p(f"duplicate-state check: {n_dupe} offending pairs")
    p()

    protos = sorted({k[5] for k in by_g})
    if len(protos) > 1:
        p(f"🔴 REFUSING to report contrasts: {len(protos)} protocols present in this "
          f"selection -- {protos}. Re-run with a --pattern that picks one.")
        p("   (per-cell rows are still printed below, keyed by protocol)")

    def rows_under(gate_names: tuple[str, ...]) -> list[dict]:
        out = []
        for k in sorted(by_g):
            r = analyse_group(by_g[k], gate_names)
            fam, cv, mu, seed, n, proto, gmax, steps = k
            out.append({"family": fam, "cv": cv, "mu_p": mu, "seed": seed,
                        "n_requested": n, "protocol": proto, "gamma_max": gmax,
                        "shear_steps": steps, **r})
        return out

    rows = rows_under(())
    p(f"{'family':10s} {'cv':>5s} {'mu':>5s} {'seed':>4s} {'cells':>5s} "
      f"{'stdy':>4s} {'qs':>3s} {'phi_c_E2*':>9s} {'b_E2':>6s} "
      f"{'phi_c_E3':>9s} {'phi_c_E4':>9s} {'phiDenseE1':>10s} {'mu_s':>6s} "
      f"{'p_hat_max':>10s} {'gdot':>10s}")
    p("   (* E2 is PRIMARY; phiDenseE1 is the retracted E1, relabelled -- see docstring)")
    p("   (this table is UNGATED; `stdy`/`qs` are how many of `cells` pass each gate)")
    for r in rows:
        p(f"{r['family']:10s} {r['cv']:5.2f} {r['mu_p']:5.2f} {r['seed']:4d} "
          f"{r['n_cells']:5d} {r['n_steady']:4d} {r['n_quasistatic']:3d} "
          f"{r['phi_c_E2']:9.4f} {r.get('b_E2', np.nan):6.2f} "
          f"{r['phi_c_E3']:9.4f} {r['phi_c_E4']:9.4f} "
          f"{r['phi_dense_E1']:10.4f} {r['mu_s']:6.3f} {r['p_hat_max']:10.3e} "
          f"{r['gdot_mean']:10.3e}")

    # --- CV contrasts, per family/mu, averaged over seeds, UNDER EVERY GATE SET ------
    # 🔴 Never one alone. See the docstring: the gates were defined and not applied
    # until 2026-09-09, the headline falls ~30 % under any of them, and there is no
    # defensible single choice because the quasi-static gate removes precisely the rung
    # nearest the root. Printing all four is the honest report.
    p()
    p("=== CV contrasts (mean over seeds, +/- sd) -- ALL FOUR GATE SETS ===")
    if len(protos) > 1:
        p("   SKIPPED -- see the protocol refusal above")
        gate_rows = {}
    else:
        gate_rows = {label: rows_under(g) for label, g in GATE_SETS}

    for fam in sorted({r["family"] for r in rows}) if gate_rows else []:
        for mu in sorted({r["mu_p"] for r in rows if r["family"] == fam}):
            base = [r for r in rows if r["family"] == fam and r["mu_p"] == mu]
            if len({r["cv"] for r in base}) < 2:
                continue
            p(f"-- {fam}, mu_p = {mu}")

            # THE SHEAR RATE IS NOT CONTROLLED ALONG THIS LADDER. Print it, and flag it.
            gd = {}
            for cv in sorted({r["cv"] for r in base}):
                v = [r["gdot_mean"] for r in base
                     if r["cv"] == cv and np.isfinite(r["gdot_mean"])]
                if v:
                    gd[cv] = (float(np.mean(v)), float(np.min(v)), float(np.max(v)))
            if gd:
                p("   realised gdot: " + "  ".join(
                    f"cv{cv:.2f}={m:.3e}" for cv, (m, _, _) in sorted(gd.items())))
                span = max(m for m, _, _ in gd.values()) / min(
                    m for m, _, _ in gd.values())
                if span > GDOT_MATCH_TOL:
                    p(f"   🔴 SHEAR RATE VARIES {span:.2f}x ACROSS THIS CONTRAST "
                      f"(bar {GDOT_MATCH_TOL:.2f}x). gdot = gamma_max/(shear_steps*dt) "
                      f"and dt is set by the SMALLEST particle, so widening the size "
                      f"distribution raises the rate. Every dphi_c below is confounded "
                      f"with a rate contrast of this size; the rate-matched control is "
                      f"--pattern \"phase8_ratem*.jsonl\".")
                else:
                    p(f"   shear rate matched to {span:.2f}x across this contrast "
                      f"(bar {GDOT_MATCH_TOL:.2f}x) -- not a rate confound")

            for tag in ("E2", "E3", "E4", "E4m", "dense_E1"):
                col = "phi_dense_E1" if tag == "dense_E1" else f"phi_c_{tag}"
                p(f"   {tag}:")
                for label, _gates in GATE_SETS:
                    sub = [r for r in gate_rows[label]
                           if r["family"] == fam and r["mu_p"] == mu]
                    vals = {}
                    for cv in sorted({r["cv"] for r in sub}):
                        v = [r[col] for r in sub
                             if r["cv"] == cv and np.isfinite(r[col])]
                        if v:
                            vals[cv] = (float(np.mean(v)),
                                        float(np.std(v, ddof=1)) if len(v) > 1
                                        else float("nan"), len(v))
                    if len(vals) < 2:
                        p(f"      {label:12s} too few resolved CV levels "
                          f"({len(vals)}) -- no contrast")
                        continue
                    lo, hi = min(vals), max(vals)
                    d = vals[hi][0] - vals[lo][0]
                    sds = [x for x in (vals[hi][1], vals[lo][1]) if np.isfinite(x)]
                    sd = max(sds) if sds else float("nan")
                    sig = d / sd if np.isfinite(sd) and sd > 0 else float("nan")
                    slope = d / (hi - lo)
                    txt = " ".join(f"cv{cv:.2f}={m:.4f}+-{s:.4f}(n{n})"
                                   for cv, (m, s, n) in sorted(vals.items()))
                    p(f"      {label:12s} d(cv {lo:.2f}->{hi:.2f}) = {d:+.4f} "
                      f"({sig:+.1f} sd) slope {slope:+.4f}/CV   {txt}")

    # --- P3: the CV contrast at MATCHED phi -----------------------------------------
    # The ladders are centred per CV, so they share no exact phi. Where they OVERLAP the
    # nearest rungs are still a direct test of the prediction that raising CV at fixed
    # volume fraction collapses the undrained grain pressure -- and it needs no fit at
    # all, so no estimator systematic enters.
    #
    # 🔴 EVERY overlapping pair is reported, not the closest one. The first version took
    # the single nearest pair; at this ladder spacing TWO pairs tie to within rounding
    # (both 0.0023 apart in phi) and they give 43x and 5.8x, so which one a min() returns
    # is a floating-point accident. Reporting one of them would have been cherry-picking
    # by tie-break.
    p()
    p("=== P3: the CV contrast at MATCHED phi (no fit, no estimator) ===")
    if len(protos) > 1:
        p("   SKIPPED -- more than one protocol present")
    else:
        good2 = [r for r in good if np.isfinite(r.get("phi", np.nan))]
        cvs = sorted({round(float(r["cv"]), 3) for r in good2})
        pairs = []
        for i_, c1 in enumerate(cvs):
            for c2 in cvs[i_ + 1:]:
                A_ = [r for r in good2 if round(float(r["cv"]), 3) == c1]
                B_ = [r for r in good2 if round(float(r["cv"]), 3) == c2]
                seen = set()
                for ra in A_:
                    for rb in B_:
                        d = abs(ra["phi"] - rb["phi"])
                        if d > PHI_MATCH_TOL:
                            continue
                        key = (c1, c2, round(ra["phi"], 4), round(rb["phi"], 4))
                        if key in seen:
                            continue
                        seen.add(key)
                        pairs.append((c1, c2, ra["phi"], rb["phi"], d))
        if not pairs:
            p(f"   no cross-CV phi pair within {PHI_MATCH_TOL}")
        else:
            p(f"   {len(pairs)} overlapping pairs within {PHI_MATCH_TOL} in phi; "
              f"ALL are reported")
            ratios = []
            for c1, c2, p1, p2, d in sorted(pairs, key=lambda t: t[2]):
                A_ = [r for r in good2 if round(float(r["cv"]), 3) == c1
                      and abs(r["phi"] - p1) < PHI_GROUP_TOL]
                B_ = [r for r in good2 if round(float(r["cv"]), 3) == c2
                      and abs(r["phi"] - p2) < PHI_GROUP_TOL]
                pa = np.array([r["p_ss_hat"] for r in A_], float)
                pb = np.array([r["p_ss_hat"] for r in B_], float)
                za = np.array([r["z_nr_ss"] for r in A_], float)
                zb = np.array([r["z_nr_ss"] for r in B_], float)
                ra_ = np.array([r["rattler_frac_ss"] for r in A_], float)
                rb_ = np.array([r["rattler_frac_ss"] for r in B_], float)
                sd = max(pa.std(ddof=1) if len(pa) > 1 else 0.0,
                         pb.std(ddof=1) if len(pb) > 1 else 0.0)
                ratio = pa.mean() / pb.mean()
                ratios.append(ratio)
                nsd = abs(pa.mean() - pb.mean()) / sd if sd > 0 else float("inf")
                p(f"   phi {p1:.5f}(cv{c1:.2f}, n{len(pa)}) vs {p2:.5f}"
                  f"(cv{c2:.2f}, n{len(pb)}), mismatch {d:.5f}")
                p(f"      p_ss {pa.mean():.4e} -> {pb.mean():.4e}  "
                  f"ratio {ratio:6.1f}x  ({nsd:.0f} sd)   "
                  f"Z_nr {za.mean():.2f} -> {zb.mean():.2f}   "
                  f"rattlers {ra_.mean():.3f} -> {rb_.mean():.3f}")
                if p2 > p1:
                    p("      (the higher-CV cell is the DENSER of the two, which works "
                      "AGAINST the effect: conservative)")
                # ROBUSTNESS: analyse_rheology flags the cv = 0.5 cells as not having
                # settled to +/-10 % within gamma_max. Recompute the ratio from a LATE
                # window (the last 15 % of the strain) and print both, rather than
                # asserting that the flag does not matter.
                def _late(rs):
                    v = []
                    for rr in rs:
                        g = np.asarray(rr.get("gamma", []), float)
                        pp = np.asarray(rr.get("p_hat_series", []), float)
                        if g.size and g.size == pp.size:
                            v.append(float(np.mean(pp[g >= 0.85 * g.max()])))
                    return float(np.mean(v)) if v else float("nan")
                la, lb = _late(A_), _late(B_)
                if np.isfinite(la) and np.isfinite(lb) and lb != 0:
                    p(f"      robustness: on a LATE window (last 15 % of strain) the "
                      f"ratio is {la / lb:.1f}x against {ratio:.1f}x -- "
                      f"a {abs(la / lb / ratio - 1) * 100:.0f} % change")
            if len(ratios) > 1:
                p(f"   spread across pairs: {min(ratios):.1f}x to {max(ratios):.1f}x "
                  f"-- a factor {max(ratios) / min(ratios):.1f} between pairs that "
                  f"differ only in which rung they sit on. Quote the RANGE.")

    path = os.path.join(OUT_DIR, args.out)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    # 🔴 The CSV name is derived from --out, not hardcoded. Until 2026-09-09 every run
    # wrote `phase8_summary.csv` whatever `--out` said, so analysing the P5 uniform arm
    # silently clobbered the production arm's summary table with uniform-family rows --
    # a derived artifact overwritten by a different measurement, which is the
    # sync_back.sh defect of 2026-09-08 in a second place. Found by a numerics reviewer.
    stem = os.path.splitext(os.path.basename(args.out))[0]
    csv = os.path.join(OUT_DIR, f"{stem.replace('_verdict', '')}_summary.csv")
    import csv as _csv
    if rows:
        keys = sorted({k for r in rows for k in r})
        with open(csv, "w", encoding="utf-8", newline="") as fh:
            w = _csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(rows)
    print(f"\nwritten -> {path}\nwritten -> {csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
