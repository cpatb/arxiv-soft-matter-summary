"""Gate N1: the six corn-starch sweeps load, in SI, with the units pinned.

Pass conditions, fixed before the run:

  N1a  all six files load with the expected schema
  N1b  eta * gamma_dot agrees with the exported stress column to < 1% relative
       on every point where the exported column is positive
  N1c  every sweep is physically admissible: positive rates and viscosities,
       no NaN in the required columns, temperature within 1 degC of 25

N1b is a test of the unit interpretation. A misread of MPa as Pa shows up as a
factor of 1e6, so passing it pins the conversion.

Run:  python scripts/gate_n1.py
Exit code is 0 only if every gate passes.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dst_fit import check_stress_consistency, read_series  # noqa: E402

DATA_DIR = Path(
    r"C:\Users\compl\Desktop\Projects\Corn Starch Bulk Rheology AM"
)
EXPECTED_PHI = [0.40, 0.45, 0.47, 0.50, 0.52, 0.54]

RTOL_STRESS = 1e-2  # N1b
TEMP_SETPOINT_C = 25.0  # N1c
TEMP_TOL_C = 1.0


def main() -> int:
    paths = sorted(DATA_DIR.glob("*_e1.xls"))
    print(f"data dir: {DATA_DIR}")
    print(f"found {len(paths)} files\n")

    failures: list[str] = []

    # ---- N1a: load ------------------------------------------------------
    try:
        sweeps = read_series(paths)
    except Exception as exc:  # noqa: BLE001 - the gate reports, it does not raise
        print(f"N1a FAIL: {type(exc).__name__}: {exc}")
        return 1

    got_phi = [s.phi for s in sweeps]
    if got_phi != EXPECTED_PHI:
        failures.append(f"N1a: phi series {got_phi} != expected {EXPECTED_PHI}")
    else:
        print(f"N1a PASS: 6 files loaded, phi = {got_phi}")

    print()
    for s in sweeps:
        print("  " + s.summary())
    print()

    # ---- N1b: unit interpretation ---------------------------------------
    worst_rel = 0.0
    for s in sweeps:
        chk = check_stress_consistency(s, rtol=RTOL_STRESS)
        if not chk["ok"]:
            failures.append(f"N1b: {s.source.name}: {chk}")
            continue
        worst_rel = max(worst_rel, chk["max_rel_resid"])
        if chk["n_excluded"]:
            print(
                f"  note: {s.source.name} excluded {chk['n_excluded']} "
                f"non-positive stress points from the check"
            )
    if not any(f.startswith("N1b") for f in failures):
        print(
            f"N1b PASS: eta*gdot matches exported stress on all six, "
            f"worst {worst_rel:.2e} relative (tol {RTOL_STRESS:.0e}, "
            f"{RTOL_STRESS / worst_rel:.0f}x margin)"
        )

    # ---- N1c: physical admissibility ------------------------------------
    for s in sweeps:
        if not np.all(np.isfinite(s.shear_rate)) or not np.all(np.isfinite(s.viscosity)):
            failures.append(f"N1c: {s.source.name}: NaN in a required column")
        if not (np.all(s.shear_rate > 0) and np.all(s.viscosity > 0)):
            failures.append(f"N1c: {s.source.name}: non-positive rate or viscosity")
        if s.temperature is not None:
            drift = float(np.max(np.abs(s.temperature - TEMP_SETPOINT_C)))
            if drift > TEMP_TOL_C:
                failures.append(
                    f"N1c: {s.source.name}: temperature drifts {drift:.2f} degC "
                    f"from {TEMP_SETPOINT_C} (tol {TEMP_TOL_C})"
                )
    if not any(f.startswith("N1c") for f in failures):
        print("N1c PASS: all sweeps positive, finite, isothermal within 1 degC")

    # ---- observations, not gates ----------------------------------------
    print("\nobservations (not gated):")
    for s in sweeps:
        eta = s.viscosity
        i = int(np.argmin(eta))
        rise = eta[-1] / eta[i]
        onset = s.stress[i]
        ns = s.normal_stress
        ns_txt = (
            f", |N| max {np.max(np.abs(ns)):.3g} Pa" if ns is not None else ", no N column"
        )
        print(
            f"  phi={s.phi:.2f}  eta_min at gdot={s.shear_rate[i]:8.3f} 1/s "
            f"(sigma={onset:7.2f} Pa), thickening {rise:5.2f}x{ns_txt}"
        )

    print()
    if failures:
        print(f"GATE N1 FAIL ({len(failures)} problem(s)):")
        for f in failures:
            print("  " + f)
        return 1
    print("GATE N1 PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
