"""check_discrete_psd.py -- the gate `DiscretePSD` has to pass before any cell runs.

Two kinds of check, and the difference between them is the point.

  (1) BIT-IDENTITY. A two-class `DiscretePSD` must return the SAME radii array, element
      for element, as the `BinaryPSD` it is the same distribution as. Not "agrees to
      1e-12": `np.array_equal`. Only the two code paths being the same code can pass it,
      so it cannot be passed by a coincidence of two nearly-right numbers -- which is
      exactly the failure mode a tolerance-based gate has.

  (2) A DISCRIMINATING THREE-CLASS CASE. Bit-identity at K = 2 says nothing about K = 3,
      because the K = 2 path is degenerate in every way the new code could be wrong:
      one pair, one rounding remainder, one ordering. So the three-class checks are
      computed by hand, independently of the class, from the definitions:

          n_i  proportional to  y_i / d_i^3,    then scaled so <d> = d10.

      A check whose expected value is produced by the code under test certifies nothing
      (lessons.md, "two anchors that share one identity cannot check that identity"),
      so every expected number below is written out in closed form.

Run:  python check_discrete_psd.py      (exit 0 = all gates pass)
"""

from __future__ import annotations

import sys

import numpy as np

from psd import BinaryPSD, DiscretePSD, make_psd, parse_discrete_family

FAILURES: list[str] = []
LINES: list[str] = []


def say(s: str = "") -> None:
    LINES.append(s)
    print(s)


def check(name: str, ok: bool, detail: str) -> None:
    say(f"   [{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    if not ok:
        FAILURES.append(name)


def main() -> int:
    say("=== check_discrete_psd.py -- gates for DiscretePSD ===")
    say()

    # ---------------------------------------------------------------- G1
    say("=== G1 -- BIT-IDENTITY against BinaryPSD (same radii array) ===")
    say("   DiscretePSD(d_rel=[1, r], y=[x_s, 1-x_s])  vs  BinaryPSD(ratio=r, x_small=x_s)")
    cases = [(2.0, 0.05), (2.0, 0.30), (3.0, 0.25), (3.0, 0.035714285714285712),
             (4.0, 0.15), (5.0, 0.10), (3.0, 0.60), (2.0, 0.95)]
    n_tot = 0
    for r, xs in cases:
        bp = BinaryPSD(d10=1.0, ratio=r, x_small=xs)
        dp = DiscretePSD(d10=1.0, d_rel=(1.0, r), y=(xs, 1.0 - xs))
        rb = bp.radii_for_yade(1500, np.random.default_rng(20260914))
        rd = dp.radii_for_yade(1500, np.random.default_rng(20260914))
        same = bool(np.array_equal(rb, rd))
        n_tot += rb.size
        check(f"r={r} x_s={xs:.6f}", same,
              f"array_equal={same}, n={rb.size}, "
              f"max|diff|={float(np.max(np.abs(rb - rd))):.3e}, "
              f"unique d = {np.unique(rb).size}")
    say(f"   N compared before any filter: {n_tot}; after: {n_tot} (nothing is dropped)")
    say()

    # ---------------------------------------------------------------- G2
    say("=== G2 -- bit-identity of the DERIVED scalars at K = 2 ===")
    for r, xs in [(3.0, 0.25), (4.0, 0.15)]:
        bp = BinaryPSD(d10=1.0, ratio=r, x_small=xs)
        dp = DiscretePSD(d10=1.0, d_rel=(1.0, r), y=(xs, 1.0 - xs))
        for label, a_, b_ in (("n_small", bp.n_small, float(dp.number_fractions[0])),
                              ("d_small", bp.d_small, float(dp.diameters[0])),
                              ("d_large", bp.d_large, float(dp.diameters[1])),
                              ("cv", bp.cv, dp.cv),
                              ("d32", bp.d32, dp.d32),
                              ("d43", bp.d43, dp.d43)):
            check(f"r={r} x_s={xs} {label}", a_ == b_,
                  f"BinaryPSD {a_!r} vs DiscretePSD {b_!r} (exact equality {a_ == b_})")
    say()

    # ---------------------------------------------------------------- G3
    say("=== G3 -- a DISCRIMINATING three-class case, expectations in closed form ===")
    # (y_L, y_M, y_S) = (0.60, 0.25, 0.15) at d = (4, 2, 1). By hand:
    #   w_S = 0.15/1 = 0.15;  w_M = 0.25/8 = 0.03125;  w_L = 0.60/64 = 0.009375
    #   sum w = 0.190625
    w = np.array([0.15 / 1.0**3, 0.25 / 2.0**3, 0.60 / 4.0**3])
    nf_hand = w / w.sum()
    # number mean of the RELATIVE diameters, then the scale that puts it at d10 = 1
    mean_rel = float(np.sum(nf_hand * np.array([1.0, 2.0, 4.0])))
    scale_hand = 1.0 / mean_rel
    d_hand = np.array([1.0, 2.0, 4.0]) * scale_hand
    m1 = float(np.sum(nf_hand * d_hand))
    m2 = float(np.sum(nf_hand * d_hand**2))
    m3 = float(np.sum(nf_hand * d_hand**3))
    m4 = float(np.sum(nf_hand * d_hand**4))
    cv_hand = float(np.sqrt(m2 - m1 * m1) / m1)

    p = DiscretePSD(d10=1.0, d_rel=(1.0, 2.0, 4.0), y=(0.15, 0.25, 0.60))
    check("K=3 number fractions", bool(np.allclose(p.number_fractions, nf_hand, atol=0, rtol=1e-14)),
          f"code {np.array2string(p.number_fractions, precision=9)} vs hand "
          f"{np.array2string(nf_hand, precision=9)}")
    check("K=3 diameters", bool(np.allclose(p.diameters, d_hand, atol=0, rtol=1e-14)),
          f"code {np.array2string(p.diameters, precision=9)} vs hand "
          f"{np.array2string(d_hand, precision=9)}")
    check("K=3 d10 == 1", abs(m1 - 1.0) < 1e-14 and abs(p.moment(1) - 1.0) < 1e-14,
          f"number mean {p.moment(1):.15f} (hand {m1:.15f})")
    check("K=3 cv", abs(p.cv - cv_hand) < 1e-14, f"{p.cv:.12f} vs hand {cv_hand:.12f}")
    check("K=3 d32", abs(p.d32 - m3 / m2) < 1e-14, f"{p.d32:.12f} vs hand {m3/m2:.12f}")
    check("K=3 d43", abs(p.d43 - m4 / m3) < 1e-14, f"{p.d43:.12f} vs hand {m4/m3:.12f}")
    # volume fractions must come back out as they went in, ASCENDING in diameter
    check("K=3 volume fractions round-trip",
          bool(np.allclose(p.volume_fractions, [0.15, 0.25, 0.60], atol=0, rtol=1e-14)),
          f"{np.array2string(p.volume_fractions, precision=12)}")
    say()

    # ---------------------------------------------------------------- G4
    say("=== G4 -- the three-class case is NOT reachable from any two-class one ===")
    say("   (a gate that a K=2 implementation could also pass has discriminated nothing)")
    d3 = p.sample(1500, np.random.default_rng(1))
    check("K=3 sample has three distinct diameters", np.unique(d3).size == 3,
          f"unique diameters {np.unique(d3).size}, values "
          f"{np.array2string(np.unique(d3), precision=6)}")
    counts = p.class_counts(1500)
    check("K=3 counts sum to N", int(counts.sum()) == 1500,
          f"counts {counts.tolist()} sum {int(counts.sum())}")
    # realised volume fractions of the ACTUAL sample -- the composition is the control
    # variable, so it is reported, never assumed.
    vreal = np.array([float(np.sum(d3[d3 == dd] ** 3)) for dd in np.unique(d3)])
    vreal = vreal / vreal.sum()
    check("K=3 realised volume fractions within 0.005 of design",
          bool(np.max(np.abs(vreal - np.array([0.15, 0.25, 0.60]))) < 0.005),
          f"realised {np.array2string(vreal, precision=6)} vs design [0.15 0.25 0.6], "
          f"max dev {float(np.max(np.abs(vreal - np.array([0.15, 0.25, 0.60])))):.2e}")
    say()

    # ---------------------------------------------------------------- G5
    say("=== G5 -- make_psd hook and the parser's strictness ===")
    fam = "discrete:d1/2/4:y0.15/0.25/0.60"
    q = make_psd(fam, d10=1.0, cv=0.0)
    rq = q.radii_for_yade(1500, np.random.default_rng(7))
    rp = p.radii_for_yade(1500, np.random.default_rng(7))
    check("make_psd round-trips the family string bit-identically",
          bool(np.array_equal(rq, rp)), f"array_equal={np.array_equal(rq, rp)}")
    check("parse_discrete_family", parse_discrete_family(fam) == ((1.0, 2.0, 4.0), (0.15, 0.25, 0.6)),
          f"{parse_discrete_family(fam)}")
    check("a non-discrete family returns None (no interception)",
          parse_discrete_family("binary:r3.0:x0.25") is None, "None")
    for bad, why in (("discrete:d1/2/4:y0.5/0.5", "length mismatch"),
                     ("discrete:d1/2/4", "missing y block"),
                     ("discrete:d1/2/4:z0.1/0.2/0.7", "wrong prefix"),
                     ("discrete:dA/2/4:y0.1/0.2/0.7", "unparsable number")):
        try:
            parse_discrete_family(bad)
            check(f"malformed {bad!r} ({why}) RAISES", False, "it did not raise")
        except ValueError as exc:
            check(f"malformed {bad!r} ({why}) RAISES", True, f"{exc}")
    try:
        make_psd(fam, d10=1.0, cv=0.3)
        check("make_psd refuses a non-zero cv for a discrete family", False, "it did not raise")
    except ValueError as exc:
        check("make_psd refuses a non-zero cv for a discrete family", True, f"{exc}")
    say()

    # ---------------------------------------------------------------- G6
    say("=== G6 -- degenerate inputs behave, and the class refuses the inadmissible ===")
    mono = DiscretePSD(d10=1.0, d_rel=(1.0, 2.0, 4.0), y=(0.0, 0.0, 1.0))
    check("a one-class composition is monodisperse at d10",
          bool(np.array_equal(mono.sample(100, np.random.default_rng(3)), np.full(100, 1.0))),
          f"cv {mono.cv}, unique d {np.unique(mono.sample(100, np.random.default_rng(3)))}")
    check("a zeroed class drops out (K=3 spec, K=2 pack)",
          bool(np.array_equal(
              DiscretePSD(d10=1.0, d_rel=(1.0, 2.0, 4.0), y=(0.25, 0.0, 0.75)
                          ).radii_for_yade(1500, np.random.default_rng(11)),
              BinaryPSD(d10=1.0, ratio=4.0, x_small=0.25
                        ).radii_for_yade(1500, np.random.default_rng(11)))),
          "identical to the equivalent BinaryPSD")
    for bad_kw, why in (
            (dict(d_rel=(1.0, 2.0), y=(0.5,)), "length mismatch"),
            (dict(d_rel=(1.0, 1.0), y=(0.5, 0.5)), "repeated diameter"),
            (dict(d_rel=(1.0, -2.0), y=(0.5, 0.5)), "negative diameter"),
            (dict(d_rel=(1.0, 2.0), y=(0.0, 0.0)), "zero total volume"),
            (dict(d_rel=(1.0, 2.0), y=(-0.1, 1.1)), "negative volume fraction")):
        try:
            DiscretePSD(d10=1.0, **bad_kw)
            check(f"refuses {why}", False, "it did not raise")
        except ValueError as exc:
            check(f"refuses {why}", True, f"{exc}")
    say()

    say(f"=== {len(FAILURES)} FAILURES: {FAILURES} ===")
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
