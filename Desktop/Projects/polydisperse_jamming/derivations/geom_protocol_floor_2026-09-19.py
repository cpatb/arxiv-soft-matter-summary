"""geom_protocol_floor_2026-09-19.py -- re-grades the geometry batch's t-tests.

WHY THIS EXISTS
---------------
A physics-referee pass on the six `geom_*_2026-09-19` derivations returned one FATAL:
every significance in items A and C2 was computed from the SEED-PAIRED standard error
alone, and this repository documents in its own artifact that that is the wrong bar.

`mixture_model.py`, in the text it prints:

    "grading scale PHI_SCALE = 2.1e-3: the protocol's own run-to-run non-determinism,
     NOT the seed sem (~2e-04), which is 10x smaller than the reproducibility of the
     thing being modelled."

and `test_mixture.py` ENFORCES that bar for a structurally identical claim
(`assert abs(d0-d3) < 1.0*mm.PHI_SCALE`). The manuscript also already carries
`\\phiScale` = 2.1e-3 as a macro, so this is grading against a scale the project declares
and then not using it.

THE PROPAGATION IS DERIVED, BECAUSE IT CHANGES VERDICTS
-------------------------------------------------------
`PHI_SCALE` is the run-to-run scatter of ONE cell's `phi_J`. A statistic built from `L`
INDEPENDENT cells, averaged over `n` seeds, therefore has a protocol floor of

    se_proto = sqrt(L) * PHI_SCALE / sqrt(n)

and `L` is the number of cells that do NOT cancel. That last clause matters: where two
families share a cell (uniform and bidisperse are bit-identical at CV = 0), that term
cancels EXACTLY and `L` drops by two. Getting it wrong moves a verdict -- see the
zero-skewness pair in `paper/extract_numbers.py`, where L = 2 gives t = 6.01 (resolved)
and L = 4 gives t = 4.25 (not).

The bar actually used is `max(seed-paired s.e., se_proto)`: the protocol floor is a floor,
not a replacement, and where the seed scatter is larger it is the honest number.
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, ROOT)

T_CRIT_2DOF = 4.303
FAILURES: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"   [{'PASS' if ok else 'FAIL'}] {name}" + (f"  --  {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


def main() -> int:
    import mixture_model as mm
    import sweep_descriptive as sd
    sd.P = lambda *a, **k: None
    fits = sd.load_fits()
    S = mm.PHI_SCALE

    print("=" * 92)
    print("Re-grading the geometry batch against the protocol's own scale")
    print("=" * 92)
    print(f"\n   PHI_SCALE = {S:.4g}   (mixture_model.py, printed in its own artifact)")
    print(f"   t_crit at 2 dof = {T_CRIT_2DOF}")

    def ser(fam, mu, cv):
        return fits[(fits.family == fam) & (fits.mu_p == mu)
                    & (fits.cv_nom == cv)].set_index("seed").phi_J_lin

    # ---- ITEM A : separability --------------------------------------------------------
    print("\n=== ITEM A  SEPARABILITY, re-graded ===")
    print("   The statistic combines FOUR cells -- (mu, cv) = (0,0), (0,0.5), (0.3,0),")
    print("   (0.3,0.5) -- none of which is shared, so L = 4.")
    print(f"   {'family':10s} {'form':10s} {'value':>10s} {'se_seed':>9s} {'t_seed':>8s} "
          f"{'se_proto':>9s} {'t_proto':>8s}  verdict")
    verdicts = {}
    for fam in ("lognormal", "uniform"):
        a0, a5 = ser(fam, 0.0, 0.0), ser(fam, 0.0, 0.5)
        b0, b5 = ser(fam, 0.3, 0.0), ser(fam, 0.3, 0.5)
        cs = sorted(set(a0.index) & set(a5.index) & set(b0.index) & set(b5.index))
        n = len(cs)
        forms = {
            "additive": (np.array([(a5[s] - a0[s]) - (b5[s] - b0[s]) for s in cs]),
                         math.sqrt(4) * S / math.sqrt(n)),
            # a RATIO of packing fractions: the absolute floor S becomes a RELATIVE one,
            # S / phi, with phi ~ the monodisperse value. Stated because it is a choice.
            "multiplic": (np.array([(a5[s] / a0[s]) - (b5[s] / b0[s]) for s in cs]),
                          math.sqrt(4) * (S / float(a0.mean())) / math.sqrt(n)),
        }
        for nm, (arr, se_proto) in forms.items():
            se_seed = float(arr.std(ddof=1) / math.sqrt(n))
            t_seed = abs(float(arr.mean())) / se_seed
            t_proto = abs(float(arr.mean())) / max(se_seed, se_proto)
            res = t_proto > T_CRIT_2DOF
            verdicts[(fam, nm)] = res
            print(f"   {fam:10s} {nm:10s} {arr.mean():+10.5f} {se_seed:9.5f} {t_seed:8.2f} "
                  f"{se_proto:9.5f} {t_proto:8.2f}  "
                  f"{'REFUTED (resolved)' if res else 'NOT resolved'}")
    print("\n   🔴 CORRECTION TO ITEM A'S VERDICT. It reported BOTH forms REFUTED, on")
    print("      t = 27.70 (additive) and 65.43 (multiplicative) from the seed s.e. alone.")
    print("      On the protocol floor NEITHER is resolved on the lognormal family, so the")
    print("      verdict for both is UNDETERMINED, not REFUTED.")
    check("🔴 item A's additive verdict moves from REFUTED to UNDETERMINED",
          not verdicts[("lognormal", "additive")], "not resolved against the protocol floor")
    check("🔴 item A's MULTIPLICATIVE verdict also moves to UNDETERMINED -- a stronger "
          "correction than the referee proposed, because it computed no separate floor "
          "for a ratio", not verdicts[("lognormal", "multiplic")],
          "not resolved against the relative protocol floor")
    print("\n   ⚠ WHAT SURVIVES: the additive form fits BETTER than the multiplicative one")
    print("     (smaller residual, on every family). That ORDERING is a descriptive fact")
    print("     about two one-parameter forms and does not need a significance test. What")
    print("     does not survive is 'REFUTED', i.e. the claim that either is excluded.")

    # ---- ITEM C2 : the zero-skewness pair ---------------------------------------------
    print("\n=== ITEM C2  THE ZERO-SKEWNESS PAIR, re-graded ===")
    print("   Here uniform and bidisperse SHARE their CV = 0 cell bit-identically, so that")
    print("   term cancels exactly and L = 2, NOT 4. Using L = 4 would flip CV = 0.5.")
    print(f"   {'CV':>5s} {'diff':>10s} {'se_seed':>9s} {'t_seed':>8s} {'se_proto':>9s} "
          f"{'t_proto':>8s}  verdict")
    n_seed_only = n_proto = 0
    for cv in (0.1, 0.2, 0.3, 0.4, 0.5):
        d = {}
        for fam in ("uniform", "bidisperse"):
            a, b = ser(fam, 0.3, 0.0), ser(fam, 0.3, cv)
            d[fam] = {s: b[s] - a[s] for s in sorted(set(a.index) & set(b.index))}
        cs = sorted(set(d["uniform"]) & set(d["bidisperse"]))
        arr = np.array([d["uniform"][s] - d["bidisperse"][s] for s in cs])
        n = len(cs)
        se_seed = float(arr.std(ddof=1) / math.sqrt(n))
        se_proto = math.sqrt(2) * S / math.sqrt(n)
        t_seed = abs(float(arr.mean())) / se_seed
        t_proto = abs(float(arr.mean())) / max(se_seed, se_proto)
        n_seed_only += t_seed > T_CRIT_2DOF
        n_proto += t_proto > T_CRIT_2DOF
        print(f"   {cv:5.1f} {arr.mean():+10.5f} {se_seed:9.5f} {t_seed:8.2f} "
              f"{se_proto:9.5f} {t_proto:8.2f}  "
              f"{'RESOLVED' if t_proto > T_CRIT_2DOF else 'not resolved'}")
    print(f"\n   resolved on the seed s.e. alone: {n_seed_only} of 5")
    print(f"   resolved on the protocol floor  : {n_proto} of 5")
    check("🔴 the zero-skewness pair's robustness drops on the protocol floor",
          n_proto < n_seed_only, f"{n_seed_only} -> {n_proto} of 5")
    check("but it SURVIVES at the widest CV, so the conclusion (skewness alone does not "
          "control phi_J) still stands", n_proto >= 1, f"{n_proto} of 5")
    print("   ⚠ paper/main.tex and paper/extract_numbers.py were corrected to this bar on")
    print("     2026-09-19 and now print BOTH counts.")

    # ---- what this correction does NOT touch -------------------------------------------
    print("\n=== WHAT THIS CORRECTION DOES NOT TOUCH ===")
    print("   Items C1, B, D's z_iso derivation and E's Furnas identity are ALGEBRA. They")
    print("   cannot fail on data and carry no t-test, so the protocol floor is irrelevant")
    print("   to them. The referee re-derived all four independently and they held.")
    print("   Item B's verdict (UNDETERMINED) leans toward 'not resolved' already, and")
    print("   under-counting noise cannot manufacture that, so it is unaffected.")
    print("   ⚠ Item D's isostatic drift (t = -7.25) uses the same machinery but there is")
    print("     NO z_iso analogue of PHI_SCALE on record, so it CANNOT be re-graded here.")
    print("     Flagged as a risk, not corrected -- and the honest statement is that its")
    print("     significance is unaudited rather than wrong.")

    print("\n" + "=" * 92)
    if FAILURES:
        print(f"FAILED {len(FAILURES)} check(s):")
        for f in FAILURES:
            print(f"   - {f}")
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
