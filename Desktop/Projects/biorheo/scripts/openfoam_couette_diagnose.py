"""Why the L3 OpenFOAM anchor's first run failed at the low-rate end, measured.

`openfoam_couette_anchor.py` gate A5 failed on 2026-09-05: 6 of 7 cases hit the
2000-iteration cap, and at the lowest Carbopol rate OpenFOAM returned a viscosity 172 %
wrong with rc = 0 and no warning. This script measures WHICH knob controls that, so the
anchor's report can quote a table rather than an adjective. Two sweeps on
`carbopol_etd2020_0p1wt_30uL`:

  * the numerics knob: the relaxation factor on U and the iteration cap, at the lowest
    rate in the envelope (the hardest corner);
  * the physical knob: the shear rate, hence the local Bingham number
    B = tau_y / (k gd^n), at fixed relaxation.

  * a DISCRIMINATOR, because the rate sweep is confounded: moving gd moves B, the
    viscosity eta, and the wall speed U0 all at once, so "B controls the convergence"
    is not established by it. Fixing gd = 0.01 1/s (hence U0, the mesh and the boundary
    conditions) and choosing synthetic Herschel-Bulkley parameters that hit a target
    (B, eta) pair separates them: at fixed gd and n,

        k gd^n = eta gd / (1 + B),      tau_y = B eta gd / (1 + B),

    so B and eta can be set independently. A 2x2 with B in {2, 29} and eta in
    {0.5, 50} Pa s. PRE-REGISTERED (written before the run): report
    R_B = the factor err_nu changes across B at fixed eta, and R_eta = the factor across
    eta at fixed B. Call it "B controls" only if R_B >= 10 * R_eta on BOTH eta rows;
    "eta controls" only if R_eta >= 10 * R_B on both B columns; otherwise "not
    separated", which is a real verdict and is what goes in the report. The synthetic
    parameter sets are not registry records and the dictionary is built by the same
    arithmetic the export uses, stated inline.

Writes `data/openfoam/couette_diagnosis.json`, which `openfoam_couette_anchor.py` reads
to build the diagnosis table in `couette_anchor.md`. Nothing here is a gate; the gates
live in the anchor. Run this before the anchor if the JSON is missing or stale.

Run:  python scripts/openfoam_couette_diagnose.py     (~10 min of solver time)
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import openfoam_couette_anchor as A  # noqa: E402

from biorheo.registry.loader import load_all  # noqa: E402
from biorheo.solver.openfoam import export_transport_properties  # noqa: E402

REC = "carbopol_etd2020_0p1wt_30uL"
# (relaxation on U, iteration cap). 0.99 is omitted from the shipped sweep because it
# dies with SIGFPE; that fact is recorded as `sigfpe_at_relax` below, measured once.
KNOBS = [(0.3, 2000), (0.3, 20000), (0.7, 2000), (0.7, 20000),
         (0.9, 20000), (0.95, 20000), (0.99, 20000)]
RATES = [None, 0.003, 0.01, 0.03, 0.1]      # None -> the envelope minimum
RELAX_FOR_RATE_SWEEP = 0.7
RESID_CONTROL = 1e-12                        # the criterion A5 originally used

# the (B, eta) discriminator, at one fixed rate so U0 and the mesh cannot move
GD_2X2 = 0.01
N_2X2 = 0.579                                # the Carbopol exponent, held fixed
B_2X2 = (2.0, 29.0)
ETA_2X2 = (0.5, 50.0)                        # Pa s
RELAX_2X2, END_2X2 = 0.7, 2000               # a setting where the failure was large
SEPARATION = 10.0                            # pre-registered factor for a verdict


def _fvsolution(relax: float, resid: float | None) -> str:
    """Patch the anchor's fvSolution. Both substitutions are asserted to have BITTEN:
    `str.replace` on a literal from another module returns the string unchanged if that
    module's wording or constant moves, and a silently un-patched fvSolution would run
    every row at the same relaxation with no residualControl -- which is the exact shape
    of the claim this script is testing, so it would look correct. Raised by the code
    review 2026-09-05."""
    fv = A._FVSOLUTION.replace(f"U {A.RELAX_U:g};", f"U {relax:g};")
    if relax != A.RELAX_U and fv == A._FVSOLUTION:
        raise RuntimeError("relaxation substitution did not bite; "
                           "openfoam_couette_anchor._FVSOLUTION has changed")
    if resid is not None:
        marker = next((ln for ln in fv.splitlines() if "no residualControl" in ln), None)
        if marker is None:
            raise RuntimeError("no residualControl marker line in _FVSOLUTION")
        fv = fv.replace(marker, f"    residualControl {{ p {resid:g}; U {resid:g}; }}")
    return fv


def run(name: str, text: str, gd: float, eta: float, relax: float, end: int) -> dict:
    """One case with residualControl restored, so 'did it meet 1e-12' is measurable."""
    base_fv, base_end = A._FVSOLUTION, A.END_TIME
    try:
        A._FVSOLUTION = _fvsolution(relax, RESID_CONTROL)
        A.END_TIME = end
        row = A.one_case(name, text, gd=gd, eta_ref=eta)
    finally:
        A._FVSOLUTION, A.END_TIME = base_fv, base_end
    row["relax_U"] = relax
    row["end_time"] = end
    row["residual_control_met"] = row["iterations"] < end
    return row


def _hb_dict(tau_y: float, k: float, n: float, eta_max: float) -> str:
    """A HerschelBulkley transportProperties for synthetic parameters, using the same
    kinematic map as `biorheo.solver.openfoam.coefficients`: tau0 = tau_y/rho, k = k/rho,
    n = n, nu0 = eta_max/rho. Not a registry record -- these parameters are constructed
    to hit a (B, eta) target, not measured, and must never enter the registry."""
    rho = A.RHO
    return ("FoamFile" + chr(10) + "{" + chr(10) + "    version     2.0;" + chr(10) +
            "    format      ascii;" + chr(10) + "    class       dictionary;" + chr(10) +
            "    object      transportProperties;" + chr(10) + "}" + chr(10) +
            "// synthetic (B, eta) discriminator cell, openfoam_couette_diagnose.py"
            + chr(10) + chr(10) + "transportModel  HerschelBulkley;" + chr(10) + chr(10) +
            "HerschelBulkleyCoeffs" + chr(10) + "{" + chr(10) +
            f"    k       {k / rho:.17g};" + chr(10) +
            f"    n       {n:.17g};" + chr(10) +
            f"    tau0    {tau_y / rho:.17g};" + chr(10) +
            f"    nu0     {eta_max / rho:.17g};" + chr(10) + "}" + chr(10))


def two_by_two() -> dict:
    """Separate B from eta at fixed gd (hence fixed U0, mesh and boundary conditions)."""
    gd, n = GD_2X2, N_2X2
    cells = []
    for B in B_2X2:
        for eta in ETA_2X2:
            k = eta * gd ** (1.0 - n) / (1.0 + B)      # from k gd^n = eta gd / (1 + B)
            tau_y = B * k * gd ** n
            eta_chk = (tau_y + k * gd ** n) / gd       # must reproduce the target eta
            B_chk = tau_y / (k * gd ** n)
            assert abs(eta_chk - eta) / eta < 1e-12 and abs(B_chk - B) / B < 1e-12
            base_fv, base_end = A._FVSOLUTION, A.END_TIME
            try:
                A._FVSOLUTION = _fvsolution(RELAX_2X2, RESID_CONTROL)
                A.END_TIME = END_2X2
                r = A.one_case(f"disc_B{B:g}_eta{eta:g}", _hb_dict(tau_y, k, n, 10.0 * eta),
                               gd=gd, eta_ref=eta)
            finally:
                A._FVSOLUTION, A.END_TIME = base_fv, base_end
            r |= {"B": B, "eta": eta, "tau_y": tau_y, "k": k, "n": n,
                  "residual_control_met": r["iterations"] < END_2X2}
            cells.append(r)
            print(f"2x2 B={B:g} eta={eta:g}: err_nu {r['err_nu']:.3e} "
                  f"resid_met {r['residual_control_met']}")

    def e(B, eta):
        return next(c["err_nu"] for c in cells if c["B"] == B and c["eta"] == eta)

    lo_B, hi_B = min(B_2X2), max(B_2X2)
    lo_e, hi_e = min(ETA_2X2), max(ETA_2X2)
    R_B = {eta: e(hi_B, eta) / e(lo_B, eta) for eta in ETA_2X2}
    R_eta = {B: e(B, hi_e) / e(B, lo_e) for B in B_2X2}
    b_wins = all(abs(np.log(R_B[eta])) >= np.log(SEPARATION) + abs(np.log(R_eta[B]))
                 for eta in ETA_2X2 for B in B_2X2)
    e_wins = all(abs(np.log(R_eta[B])) >= np.log(SEPARATION) + abs(np.log(R_B[eta]))
                 for eta in ETA_2X2 for B in B_2X2)
    verdict = "B controls" if b_wins else ("eta controls" if e_wins else "not separated")
    print(f"2x2 verdict: {verdict}  R_B {R_B}  R_eta {R_eta}")
    return {"gammadot": gd, "n": n, "relax_U": RELAX_2X2, "end_time": END_2X2,
            "B_values": list(B_2X2), "eta_values": list(ETA_2X2),
            "separation_factor": SEPARATION, "cells": cells,
            "R_B_at_fixed_eta": {str(k_): v for k_, v in R_B.items()},
            "R_eta_at_fixed_B": {str(k_): v for k_, v in R_eta.items()},
            "verdict": verdict}


def main() -> int:
    if not A.has_openfoam():
        print("SKIP: no WSL OpenFOAM v2606")
        return 0
    rec = load_all()[REC]
    p = dict(zip(rec.model.param_names, [float(v) for v in rec.params_si]))
    gd_lo = float(rec.envelope.gammadot_min)
    eta_max = 10.0 * float(rec.eta_at(gd_lo))
    text = export_transport_properties(rec, A.RHO, eta_max)

    def bingham(gd: float) -> float:
        return p["tau_y"] / (p["k"] * gd ** p["n"])

    knob_rows, sigfpe = [], []
    for relax, end in KNOBS:
        try:
            r = run(f"diag_r{relax}_e{end}", text, gd_lo, float(rec.eta_at(gd_lo)), relax, end)
            r["B"] = bingham(gd_lo)
            knob_rows.append(r)
            print(f"relax {relax} end {end}: cauchy {r['cauchy']:.3e} "
                  f"err_nu {r['err_nu']:.3e} resid_met {r['residual_control_met']}")
        except (RuntimeError, subprocess.TimeoutExpired) as exc:
            # a timeout is a crash row, not a reason to abandon the sweep
            sigfpe.append({"relax_U": relax, "end_time": end, "error": str(exc)[:200]})
            print(f"relax {relax} end {end}: CRASH {str(exc)[:80]}")

    rate_rows = []
    for gd in RATES:
        gd = gd_lo if gd is None else gd
        r = run(f"diag_g{gd:g}", text, gd, float(rec.eta_at(gd)), RELAX_FOR_RATE_SWEEP, 20000)
        r["B"] = bingham(gd)
        rate_rows.append(r)
        print(f"gd {gd:g} (B {r['B']:.4g}): err_nu {r['err_nu']:.3e} "
              f"its {r['iterations']} resid_met {r['residual_control_met']}")

    grid = two_by_two()

    out = {"record": REC, "params": p, "rho": A.RHO, "H": A.H, "NY": A.NY,
           "discriminator": grid,
           "residual_control": RESID_CONTROL, "relax_for_rate_sweep": RELAX_FOR_RATE_SWEEP,
           "knob_sweep": knob_rows, "rate_sweep": rate_rows, "crashes": sigfpe,
           "envelope": [gd_lo, float(rec.envelope.gammadot_max)]}
    path = A.OUT / "couette_diagnosis.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2), encoding="ascii")
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
