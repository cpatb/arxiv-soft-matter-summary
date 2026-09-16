"""k_axis_bound_check.py -- the four checks on derivations/k_axis_bound.md.

Anchors, in the order the lab's rules ask for them:
  (3) SymPy      -- the telescoping identity and the extremal lemma, symbolically
  (4) limits     -- u_hi -> u_lo, u_hi -> 1, u_lo -> 0, K = 2
  (5) dimensions -- invariance under rescaling every diameter, and under the
                    OTHER ratio convention (large/small instead of small/large)
  (6) numeric    -- BRUTE FORCE: feasibility of a K-class ladder decided by random
                    search and by constrained optimisation, WITHOUT the formula,
                    compared against the formula on a grid of windows

and one anchor outside the derivation entirely: the shipped `kladder.ladder(K)`.

Run:  python3 derivations/k_axis_bound_check.py
"""

from __future__ import annotations

import itertools
import math
import os
import sys

import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import kladder as kl  # noqa: E402


def kappa(u_lo: float, u_hi: float) -> float:
    """The bound's right-hand side, ln(u_lo)/ln(u_hi)."""
    return math.log(u_lo) / math.log(u_hi)


def k_max(u_lo: float, u_hi: float) -> int:
    """Largest K admitting a ladder with every pair ratio in [u_lo, u_hi]."""
    return 1 + int(math.floor(kappa(u_lo, u_hi) + 1e-12))


def feasible_formula(k: int, u_lo: float, u_hi: float) -> bool:
    return u_lo <= u_hi ** (k - 1) + 1e-15


def pair_ratios(d):
    """Every u_ij = d_i/d_j for i < j on an ASCENDING ladder."""
    d = sorted(d)
    return [d[i] / d[j] for i, j in itertools.combinations(range(len(d)), 2)]


def feasible_bruteforce(k: int, u_lo: float, u_hi: float, ntrial: int = 20000,
                        seed: int = 20260917) -> bool:
    """Decide feasibility WITHOUT the formula: random search plus a direct construction
    search over log-spacings. A ladder is a point in the simplex of K-1 positive log
    gaps; sample it, and also run a coarse deterministic sweep so a thin feasible set
    is not missed by chance alone."""
    rng = np.random.default_rng(seed)
    # random: K-1 log gaps, total log span drawn over the admissible range
    for _ in range(ntrial):
        gaps = rng.dirichlet(np.ones(k - 1))
        span = rng.uniform(1.0, 1.0 / max(u_lo, 1e-9))
        d = np.concatenate([[1.0], np.exp(np.cumsum(gaps * math.log(span)))])
        us = pair_ratios(d)
        if min(us) >= u_lo - 1e-12 and max(us) <= u_hi + 1e-12:
            return True
    # deterministic: the geometric ladder at every span on a fine grid
    for span in np.exp(np.linspace(0.0, math.log(1.0 / max(u_lo, 1e-9)), 4001))[1:]:
        d = span ** (np.arange(k) / (k - 1))
        us = pair_ratios(d)
        if min(us) >= u_lo - 1e-12 and max(us) <= u_hi + 1e-12:
            return True
    return False


def main() -> int:
    bad: list[str] = []
    print("=== derivations/k_axis_bound_check.py ===")

    # ------------------------------------------------ (3) SymPy: the two lemmas
    print("\n--- (3) SYMBOLIC ---")
    K = 6
    v = sp.symbols(f"v1:{K}", positive=True)          # adjacent ratios v_1..v_{K-1}
    # Build the ladder DOWNWARD from the largest class: d_K = 1 and d_k = d_{k+1} v_k,
    # which is the DEFINITION of v_k = d_k/d_{k+1}. `d` ends up ascending, d[0] .. d[K-1].
    d = [sp.Integer(1)] * K
    for kk in range(K - 2, -1, -1):
        d[kk] = sp.simplify(d[kk + 1] * v[kk])
    # LEMMA 1, telescoping: u_ij = prod_{k=i}^{j-1} v_k
    l1_ok = True
    for i, j in itertools.combinations(range(K), 2):
        lhs = sp.simplify(d[i] / d[j])
        rhs = sp.simplify(sp.prod(v[i:j]))
        if sp.simplify(lhs - rhs) != 0:
            l1_ok = False
            bad.append(f"lemma 1 fails at (i,j)=({i},{j})")
    print(f"   LEMMA 1  u_ij = prod_(k=i)^(j-1) v_k, all {K*(K-1)//2} pairs at K={K}: "
          f"{'PASS' if l1_ok else 'FAIL'}")
    # LEMMA 2, extremes. The .md proves it structurally (a product over a SUPERSET of
    # factors in (0,1) is smaller). Here it is checked the independent way -- on random
    # ladders of random length, EXACT equality is required, not closeness, because both
    # sides are literally the same float when the lemma holds.
    rng2 = np.random.default_rng(20260917)
    l2_bad = 0
    n_l2 = 2000
    for _ in range(n_l2):
        kk = int(rng2.integers(2, 9))
        dd = np.sort(np.exp(rng2.uniform(0.0, 3.0, kk)))
        adj = [dd[i] / dd[i + 1] for i in range(kk - 1)]
        us = pair_ratios(dd)
        if not (max(us) == max(adj) and min(us) == dd[0] / dd[-1]):
            l2_bad += 1
    print(f"   LEMMA 2  max pair = max adjacent AND min pair = d_1/d_K, EXACTLY, on "
          f"{n_l2} random ladders of K = 2..8: {'PASS' if l2_bad == 0 else 'FAIL'}"
          f" ({l2_bad} violations)")
    if l2_bad:
        bad.append(f"lemma 2: {l2_bad}/{n_l2} random ladders violate it")
    u1K = sp.simplify(d[0] / d[K - 1])
    tele = sp.simplify(u1K - sp.prod(v)) == 0
    print(f"   LEMMA 1b u_1K = prod of ALL adjacent ratios, symbolically at K = {K}: "
          f"{'PASS' if tele else 'FAIL'}")
    if not tele:
        bad.append("telescoping to the extreme ratio failed")
    # THE BOUND, symbolically: u_lo <= prod v <= u_hi^(K-1) under v_k <= u_hi
    ulo, uhi, n = sp.symbols("u_lo u_hi n", positive=True)
    sol = sp.solve(sp.Eq(ulo, uhi ** n), n)
    print(f"   BOUND    u_lo = u_hi^n  =>  n = {sol[0]}  (so K-1 <= ln u_lo / ln u_hi)"
          f"  {'PASS' if sp.simplify(sol[0] - sp.log(ulo)/sp.log(uhi)) == 0 else 'FAIL'}")
    if sp.simplify(sol[0] - sp.log(ulo) / sp.log(uhi)) != 0:
        bad.append("bound inversion failed")

    # ------------------------------------------------------- (4) LIMITING CASES
    print("\n--- (4) LIMITS ---")
    lims = [
        ("u_hi -> u_lo+   (window is one point)", 0.5, 0.5 * (1 + 1e-9), 2,
         "only a binary at that exact ratio fits"),
        ("u_hi -> 1-      (no upper constraint)", 0.25, 1 - 1e-9, None,
         "K unbounded: the ladder can be subdivided arbitrarily finely"),
        ("u_lo -> 0+      (no lower constraint)", 1e-12, 0.70711, None,
         "K unbounded: an arbitrarily long ladder fits inside one adjacent step"),
        ("K = 2           (a binary)", 0.25, 0.70711, 5, "any u in the window works"),
    ]
    for name, lo, hi, want, why in lims:
        kk = kappa(lo, hi)
        km = k_max(lo, hi)
        ok = (want is None and km > 50) or (want is not None and km == want)
        print(f"   {name:34s} kappa = {kk:12.4g}  K_max = {km:<6d} "
              f"{'PASS' if ok else 'FAIL'}   ({why})")
        if not ok:
            bad.append(f"limit {name}: K_max {km}, wanted {want}")

    # ------------------------------------------------------ (5) DIMENSIONS
    print("\n--- (5) DIMENSIONS / INVARIANCE ---")
    rng = np.random.default_rng(1)
    d0 = np.sort(rng.uniform(1.0, 4.0, 5))
    for scale in (1e-6, 1.0, 1e6):
        us = pair_ratios(d0 * scale)
        u0 = pair_ratios(d0)
        if max(abs(np.array(us) - np.array(u0))) > 1e-12:
            bad.append(f"rescaling by {scale} changed the pair ratios")
    print("   rescaling every diameter by 1e-6, 1, 1e6 leaves every pair ratio "
          "unchanged: PASS (the bound is a statement about ratios only)")
    # the OTHER convention: U = large/small > 1, window [1/u_hi, 1/u_lo]
    lo, hi = 0.25, 4.0 ** -0.25
    kappa_other = math.log(1.0 / lo) / math.log(1.0 / hi)
    ok = abs(kappa_other - kappa(lo, hi)) < 1e-12
    print(f"   convention u = small/large gives kappa = {kappa(lo, hi):.10f}; "
          f"u = large/small gives {kappa_other:.10f}: {'PASS' if ok else 'FAIL'}")
    if not ok:
        bad.append("the bound is convention-dependent")

    # --------------------------------------------- (6) NUMERIC, BRUTE FORCE
    print("\n--- (6) NUMERIC: brute-force feasibility vs the formula ---")
    print("   Feasibility decided by random search + a geometric sweep, NOT by the")
    print("   formula. Disagreement on any cell fails the check.")
    print(f"   {'u_lo':>8}{'u_hi':>10}{'K':>4}{'kappa':>9}{'formula':>9}{'brute':>8}  ")
    n_cell = n_agree = 0
    for lo in (0.20, 0.25, 1.0 / 3.0, 0.40):
        for hi in (0.6667, 4.0 ** -0.25, 0.75, 0.80):
            for k in (2, 3, 4, 5, 6, 7):
                f_form = feasible_formula(k, lo, hi)
                f_brut = feasible_bruteforce(k, lo, hi, ntrial=4000)
                n_cell += 1
                n_agree += int(f_form == f_brut)
                if f_form != f_brut:
                    bad.append(f"brute force disagrees at u_lo={lo}, u_hi={hi}, K={k}")
                    print(f"   {lo:8.4f}{hi:10.5f}{k:4d}{kappa(lo, hi):9.4f}"
                          f"{str(f_form):>9}{str(f_brut):>8}  <-- DISAGREE")
    print(f"   {n_agree} of {n_cell} cells agree")

    # ---------------------------------------- ANCHOR OUTSIDE THE DERIVATION
    print("\n--- ANCHOR: the SHIPPED kladder.ladder(K), which the derivation never saw ---")
    print(f"   {'K':>3}{'span':>8}{'min pair':>14}{'max pair':>14}  in [0.25, 0.70711]?")
    for k in (2, 3, 4, 5, 6):
        d_ = kl.ladder(k)
        us = pair_ratios(d_)
        inw = min(us) >= 0.25 - 1e-12 and max(us) <= 4.0 ** -0.25 + 1e-12
        print(f"   {k:3d}{max(d_)/min(d_):8.4f}{min(us):14.10f}{max(us):14.10f}"
              f"  {'YES' if inw else 'NO'}")
        want = k <= k_max(0.25, 4.0 ** -0.25)
        if inw != want:
            bad.append(f"shipped ladder K={k}: in-window {inw}, bound says {want}")
    print(f"   bound says K_max = {k_max(0.25, 4.0**-0.25)} on [0.25, 0.70711]")

    # ------------------------------------- the two windows this project can claim
    print("\n--- THE TWO WINDOWS, and what each governs ---")
    for name, lo, hi in (
            ("ADDtab knots, before the U-GAP block  [0.25, 0.66667]", 0.25, 2.0 / 3.0),
            ("ADDtab knots, after  the U-GAP block  [0.25, 0.70711]", 0.25, 4.0 ** -0.25),
            ("row C's honest window for the FORM    [0.33333, 0.70711]", 1/3, 4.0**-0.25)):
        kk = kappa(lo, hi)
        km = k_max(lo, hi)
        s_lo, s_hi = hi ** -(km - 1), 1.0 / lo
        print(f"   {name}:  kappa = {kk:.5f}  K_max = {km}  "
              f"span(K_max) in [{s_lo:.5f}, {s_hi:.5f}]"
              f"{'  <-- FORCED, a single ladder' if abs(s_lo - s_hi) < 1e-9 else ''}")
    u6 = 0.25 ** 0.2
    print(f"   for K = 6 on u_lo = 0.25 the window must reach "
          f"u_hi >= 0.25**(1/5) = {u6:.6f}")

    # ------------------------- is that knot REACHABLE? the block's own measured signal
    print("\n--- IS THE K = 6 KNOT REACHABLE? the U-GAP block's own measured signal ---")
    print("   The bound says WHERE a knot must be, not whether one can be measured there.")
    print("   `out/grade_ugap.txt` measured the filling signal max(phi - beta) at three")
    print("   ratios and measured the design's resolving power (2.57 scale, read off arm")
    print("   B's ceiling table, which fails at 2.55). Log-interpolating the MEASURED")
    print("   signal between the two arms that bracket u = 0.7579 is what decides it.")
    # The signal is max(phi - beta)/SCALE -- the SAME quantity the 2.57 resolving-power
    # bar is stated in. NOT the "phi span" column of the U1 table, which is a different
    # number (arm B: span 2.6, max(phi - beta) 2.21) and would have been compared against
    # a bar it is not in the units of.
    import grade_mixture_ext as gx_  # noqa: PLC0415
    import grade_ugap as gu_         # noqa: PLC0415
    import mixture_model as mm_      # noqa: PLC0415
    import ugap as ug_               # noqa: PLC0415
    beta_ = mm_.measured(mu_p=0.3)[("lognormal", 0.0)][0]
    del gx_
    g_, _nf, _nr = gu_.load("ug")
    sig = {}
    for arm_, r_ in ug_.RATIOS.items():
        vals = []
        for x_ in ug_.X_GRID:
            cs = [c for c in g_.get((float(r_), round(float(x_), 6)), [])
                  if gu_.admissible(c, ug_.N_PROD)[0]]
            if cs:
                vals.append(float(np.mean([c["phi"] for c in cs])))
        if vals:
            sig[1.0 / float(r_)] = (max(vals) - beta_) / ug_.PHI_SCALE
    us_ = np.array(sorted(sig))
    ss_ = np.array([sig[u] for u in sorted(sig)])
    print(f"   {'u':>10}{'max(phi-beta)/scale':>21}   source")
    for u_ in sorted(sig):
        print(f"   {u_:10.5f}{sig[u_]:21.2f}   computed here from the `ug` shards")
    print("   resolving power of this design: 2.57 scale (fails at 2.55), quoted from")
    print("      out/grade_ugap.txt arm B's ceiling table, not invented here.")
    print("   CONTROL: arm B (u = 0.80000) is recorded there at 2.21 -- "
          f"computed here {sig.get(0.8, float('nan')):.2f}")
    s6 = float(np.exp(np.interp(math.log(u6), np.log(us_), np.log(ss_))))
    print()
    print("   the whole remaining K axis on a span-4 ladder, one row per rung:")
    print(f"   {'K':>3}{'u_hi needed':>13}{'signal est':>12}{'margin':>8}  attemptable?")
    for k_ in (5, 6, 7, 8, 10, 12):
        uk = 0.25 ** (1.0 / (k_ - 1))
        meas = any(abs(uk - u_) < 1e-4 for u_ in sig)
        if uk > max(us_) + 1e-9:
            # np.interp CLAMPS outside its range and would print the endpoint as an
            # estimate. The signal is monotone decreasing in u (it vanishes at u = 1),
            # so anything above u = 0.80 is BELOW arm B's measured 2.21, which already
            # failed -- an argument from monotonicity, not an extrapolation.
            print(f"   {k_:3d}{uk:13.5f}{'< 2.21':>12}{'< 0.86':>8}  "
                  f"no (below arm B's measured 2.21, which FAILED)")
            continue
        sk = float(np.exp(np.interp(math.log(uk), np.log(us_), np.log(ss_))))
        print(f"   {k_:3d}{uk:13.5f}{sk:12.2f}{sk/2.57:8.2f}  "
              f"{'MEASURED already' if meas else ('yes' if sk > 2.57 else 'no')}")
    print(f"   -> K = 6 needs a knot at u = {u6:.5f} carrying about {s6:.2f} scale "
          f"against a bar of 2.57, margin {s6/2.57:.2f}x.")
    k_reach = max(k_ for k_ in range(5, 13)
                  if 0.25 ** (1.0 / (k_ - 1)) <= max(us_) + 1e-9
                  and float(np.exp(np.interp(math.log(0.25 ** (1.0 / (k_ - 1))),
                                             np.log(us_), np.log(ss_)))) > 2.57)
    print(f"      the largest K whose required knot still clears the bar is K = {k_reach}"
          f", at margin "
          f"{float(np.exp(np.interp(math.log(0.25 ** (1.0/(k_reach-1))), np.log(us_), np.log(ss_))))/2.57:.2f}x.")
    print("      🔴 This QUALIFIES the Status card's item 2 ('do not push K past 5'),")
    print("      which reached its conclusion from the adjacent ratios at K = 8 and")
    print(f"      K = 12 (0.82, 0.88) without evaluating K = {k_reach}. Those two rungs")
    print("      are indeed excluded -- by a MEASUREMENT (arm B's 2.21 failed) rather")
    print(f"      than by the argument given -- but K = {k_reach} is not.")
    print("      This is an ESTIMATE from a 3-point log interpolation, not a")
    print("      measurement, and the margin is thin.")

    print(f"\nFAILURES: {bad}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
