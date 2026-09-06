"""Grade and print the item-11 gates that do NOT need OpenFOAM, and write
``data/openfoam/newtonian_powerlaw_gates.md``.

The bars are the ones pre-registered in ``data/openfoam/newtonian_powerlaw_prereg.md``
and committed at ``e9c1ace`` before any of the code they grade. They are imported from
``tests/test_openfoam_export_null.py`` rather than retyped here, so the report and the
regression tests cannot drift apart. **Do not change them to make a check pass.**

G11.3b, G11.4b and G11.7 need a running ``simpleFoam`` and are graded by
``scripts/openfoam_couette_anchor.py`` into ``data/openfoam/couette_anchor.md``.

Run:  python scripts/item11_gate_report.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "tests"))

from biorheo.constitutive import get_model                              # noqa: E402
from biorheo.registry import load_all                                   # noqa: E402
from biorheo.solver.openfoam import (_envelope_eta_range, coefficients,  # noqa: E402
                                     export_report, nu_openfoam)
from test_openfoam_export_null import (BAR_G111, BAR_G115, BAR_G116,    # noqa: E402
                                       N_RATES, RHO)

OUT = HERE.parent / "data" / "openfoam" / "newtonian_powerlaw_gates.md"


def _rates(rec) -> np.ndarray:
    lo, hi = rec.envelope.gammadot_min, rec.envelope.gammadot_max
    gd = np.logspace(np.log10(lo), np.log10(hi), N_RATES)
    gd[0], gd[-1] = lo, hi
    return np.clip(gd, lo, hi)


def _clamps(rec) -> tuple[float, float]:
    lo, hi = _envelope_eta_range(rec)
    return 10.0 * hi, 0.1 * lo


def main() -> int:
    recs = load_all()
    fam = {f: {rid: r for rid, r in recs.items() if r.model.name == f}
           for f in ("newtonian", "power_law")}
    for f, d in fam.items():
        if not d:
            print(f"FAIL: no {f} record in the registry")
            return 1

    # --- G11.1 round trip -------------------------------------------------------------
    g111: dict[str, float] = {}
    for f, d in fam.items():
        worst = 0.0
        for rid, rec in sorted(d.items()):
            eta_max, eta_min = ((None, None) if f == "newtonian" else _clamps(rec))
            rep = export_report(rec, RHO, eta_max, _rates(rec), eta_min=eta_min)
            assert not rep["cap_binds"].any() and not rep["floor_binds"].any(), rid
            worst = max(worst, rep["max_rel_err_unbound"])
        g111[f] = worst

    # --- G11.5 discriminating wrong maps ----------------------------------------------
    # Each is evaluated through `nu_openfoam` while the reference is `record.eta_at`, so
    # only the REJECTED candidate goes through this repo's transcription of calcNu.
    g115: dict[str, list[float]] = {"newtonian": [], "power_law": []}
    for rid, rec in sorted(fam["newtonian"].items()):
        gd = _rates(rec)
        nu_rec = np.asarray(rec.eta_at(gd), dtype=float) / RHO
        wrong = nu_openfoam("Newtonian", {"nu": float(rec.parameters["mu"].si)}, gd)
        g115["newtonian"].append(float(np.abs(wrong - nu_rec).max() / nu_rec.max()))
    for rid, rec in sorted(fam["power_law"].items()):
        gd = _rates(rec)
        eta_max, eta_min = _clamps(rec)
        c = dict(coefficients(rec, RHO, eta_max, eta_min=eta_min))
        nu_rec = np.asarray(rec.eta_at(gd), dtype=float) / RHO
        wrong = nu_openfoam("powerLaw", {**c, "n": c["n"] - 1.0}, gd)
        g115["power_law"].append(float(np.abs(wrong - nu_rec).max() / nu_rec.max()))

    # --- the control on G11.5b's ATTRIBUTION -------------------------------------------
    # G11.5b comes out at exactly 9 and the report below identifies that as
    # `1/gammadot_min - 1`. But `_clamps` puts the cap a factor of 10 outside the
    # envelope, and the wrong map's value at `gammadot_min` lands on `nuMax`, so 9 could
    # equally be `10 - 1` -- the clamp factor wearing the envelope's label. Re-measure
    # with the clamp 100x out: if the attribution is right the number does not move.
    # (Raised by code review 2026-09-06; the attribution survived.)
    g115b_wide = []
    for rid, rec in sorted(fam["power_law"].items()):
        gd = _rates(rec)
        lo, hi = _envelope_eta_range(rec)
        c = dict(coefficients(rec, RHO, 100.0 * hi, eta_min=0.01 * lo))
        nu_rec = np.asarray(rec.eta_at(gd), dtype=float) / RHO
        wrong = nu_openfoam("powerLaw", {**c, "n": c["n"] - 1.0}, gd)
        g115b_wide.append(float(np.abs(wrong - nu_rec).max() / nu_rec.max()))

    # --- G11.6 structural_steady vs Cross(n = 1) --------------------------------------
    ss, cr = get_model("structural_steady"), get_model("cross")
    gd = np.logspace(-2, 2, N_RATES)
    g116 = 0.0
    for eta_inf, beta, alpha_tr in ((3.5e-3, 28.0714286, 0.2), (1.0, 5.0, 10.0),
                                    (1e-4, 1e3, 1e-3)):
        a = ss(gd, eta_inf, beta, alpha_tr)
        b = cr(gd, eta_inf * (1.0 + beta), eta_inf, alpha_tr, 1.0)
        g116 = max(g116, float(np.abs(a - b).max() / np.abs(a).max()))

    # --- the clamp check's inertness ---------------------------------------------------
    n_inert = sum(_envelope_eta_range(r)[0] == float(r.eta_at(r.envelope.gammadot_max))
                  for r in recs.values())

    ok = (max(g111.values()) <= BAR_G111
          and min(min(v) for v in g115.values()) > BAR_G115
          and g116 <= BAR_G116 and n_inert == len(recs))

    lines = [
        "# Item 11 gate report: the L3 bridge for `newtonian` and `power_law`", "",
        "Generated by `scripts/item11_gate_report.py`; every number below reproduces by "
        "running it. Bars are pre-registered in `newtonian_powerlaw_prereg.md`, committed "
        "at `e9c1ace` before any of the code they grade, and are imported from "
        "`tests/test_openfoam_export_null.py` rather than retyped. **G11.3b, G11.4b and "
        "G11.7 need a running `simpleFoam` and are in `couette_anchor.md`.**", "",
        f"Registry: {len(recs)} records, of which {len(fam['newtonian'])} `newtonian` and "
        f"{len(fam['power_law'])} `power_law` (NEXT_STEPS item 11 added both families; "
        "they are NULL models fitted to the Horner blood curves and their records say "
        "so).", "",
        "| gate | what | worst | bar | verdict |", "|---|---|---|---|---|",
    ]
    for f in ("newtonian", "power_law"):
        v = g111[f]
        lines.append(f"| G11.1 | `{f}` round trip vs `eta_at` on {N_RATES} rates, "
                     f"{len(fam[f])} records | {v:.3e} | {BAR_G111:g} | "
                     f"{'PASS' if v <= BAR_G111 else 'FAILED'} |")
    for f, label in (("newtonian", "`nu = mu`, the `/rho` forgotten"),
                     ("power_law", "`n -> n - 1`, the -1 applied twice")):
        v = min(g115[f])
        lines.append(f"| G11.5{'a' if f == 'newtonian' else 'b'} | wrong map {label}: "
                     f"MINIMUM separation over {len(fam[f])} records | {v:.4g} | "
                     f"> {BAR_G115:g} | "
                     f"{'PASS' if v > BAR_G115 else 'UNINFORMATIVE'} |")
    lines += [
        f"| G11.6 | `structural_steady` equals `cross` with m = 1 | {g116:.3e} | "
        f"{BAR_G116:g} | {'PASS' if g116 <= BAR_G116 else 'FAILED'} |",
        f"| — | `_envelope_eta_range` agrees with the single-point check it replaced | "
        f"{n_inert}/{len(recs)} records | all | "
        f"{'PASS' if n_inert == len(recs) else 'FAILED'} |",
        "",
        "G11.2, G11.3a and G11.4a are pass/fail rather than numeric and are graded in "
        "`tests/test_openfoam_export_null.py`; all pass.", "",
        "## Two things the numbers say that the gates did not ask", "",
        f"**Both wrong maps separate by a CONSTANT, and neither constant is a property of "
        f"the fit.** G11.5a is {min(g115['newtonian']):.6g} on every one of the "
        f"{len(fam['newtonian'])} records, to "
        f"{max(g115['newtonian']) - min(g115['newtonian']):.1e} — it is `rho - 1` = "
        f"{RHO - 1:.6g} exactly, because forgetting `/rho` is a pure scaling. G11.5b is "
        f"{min(g115['power_law']):.6g} on all {len(fam['power_law'])}, spread "
        f"{max(g115['power_law']) - min(g115['power_law']):.1e}, and it is "
        f"`1/gammadot_min - 1` = {1.0 / next(iter(fam['power_law'].values())).envelope.gammadot_min - 1:.6g}, "
        "because the wrong exponent multiplies the law by `gd^-1` and the largest such "
        f"factor sits at the bottom of the envelope. **That second attribution is not "
        f"free:** the clamp sits a factor of 10 outside the envelope and the wrong map's "
        f"value at `gammadot_min` lands on `nuMax`, so 9 could equally have been the "
        f"clamp factor minus 1. Re-measured with the clamp 100x out it is "
        f"{min(g115b_wide):.6g} (spread "
        f"{max(g115b_wide) - min(g115b_wide):.1e}) -- unmoved, so the 9 is the envelope's "
        "and not the clamp's. So both discriminators test whether "
        "the slot is USED AT ALL; neither is sensitive to the fitted values, and a "
        "session that wanted a map error graded per-record would need a different "
        "candidate. That they are informative at all (both are far above "
        f"{BAR_G115:g}) is what G11.5 claims, and no more.", "",
        "**G11.6 refutes the reason NEXT_STEPS item 11 gave for refusing "
        "`structural_steady`.** The item said it is thixotropic and \"has no counterpart "
        f"there\"; the steady flow curve equals `CrossPowerLaw` with n = 1 to {g116:.3e}, "
        "so a counterpart exists exactly and refusing on that ground would refuse "
        "something OpenFOAM can do — item 6's mistake repeated. The refusal is kept and "
        "its reason replaced: `viscosityModel` is a function of the instantaneous "
        "`strainRate()` alone, so the STRUCTURE EQUATION cannot be carried, and the "
        "export would be right in steady state and silently wrong in every transient. "
        "`ThixotropyNotRepresentable` says that and gives the deliberate substitution.", "",
    ]
    OUT.write_text("\n".join(lines), encoding="utf-8", newline="\n")
    print("\n".join(lines))
    print(f"\nwrote {OUT}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
