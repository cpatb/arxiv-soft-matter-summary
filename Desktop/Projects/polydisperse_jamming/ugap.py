"""ugap.py -- the design of the U-GAP block: measure a(u), b(u) where the K-ladder reads them.

PRE-REGISTRATION: PREREG_UGAP_2026-09-16.md. Predict with predict_ugap.py (BLIND, before
any cell). Grade with grade_ugap.py.

WHY THIS BLOCK EXISTS
---------------------
`CLAIMS_KLADDER_2026-09-15.md` C1a(ii): the tabulated closure has NO measured knot above
u = 0.667, and 4 of K = 5's 10 pairs -- its most numerous pair type -- sit at
u = 4^(-1/4) = 0.70711, inside that gap. Two diagnostics decided that this matters and
what to do about it, both run before this file was written:

  * `out/diag_ugap.txt`: sweeping a(u*), b(u*) over the MONOTONE-admissible range moves
    the K = 5 cells' ADDtab by up to 9.59 scale -- 4.7x TAB's own in-sample rms of 2.02 --
    and carries R(5) - R(4) across zero on 12 of 12 rows. The K = 2/3 control moves
    0.00 exactly and K = 4 moves 3.49 (through the interpolant's shape at u = 0.630, not
    through a pair at u*). So the region is NOT degenerate and the caveat is load-bearing.
  * `out/diag_uident.txt`: the failure at r = 1.25 is NOT noise and NOT sample size. The
    profile bar is delta chi2 = 1 PER POINT, so duplicating a point set leaves the
    interval IDENTICAL (demonstrated: n = 8/16/32 all give w(a) = 0.0400). What decides
    the width is the CURVATURE, i.e. where in composition the block samples. And the
    closure is a MIN over branches with an escape at a = 1 or b = 1 that collapses it to
    the CONSTANT phi = beta; at r = 1.25 that escape scores chi2 = 7.09 against a bar of
    8, so every (a, b) that kills the binding branch fits. u = 0.8 constrains neither
    parameter NOR any combination of them.

THREE THINGS THAT FOLLOW, and they are the whole design
-------------------------------------------------------
1. THE RATIO IS 1.41421, NOT 1.25. K = 5 reads the closure at u = 0.70711. The noise-free
   identifiability ceiling is already YES there on the OLD composition grid
   (w(a) 0.040, w(b) 0.080) and NO only at r = 1.25. The block that would rescue u = 0.8
   is not the block the K-ladder needs.
2. THE COMPOSITION GRID MUST STRADDLE THE BRANCH CROSSOVER. a enters only the coarse
   branch and b only the fine one, and they swap at x* where (1 - a) x = (1 - beta)(1 - b)
   (1 - x) -- x* ~ 0.356 at these ratios. The OLD grid puts 6 points below and 2 above, so
   b is poorly placed; the MID grid (x in [0.30, 0.75]) is far worse and fails a at EVERY
   ratio. A grid spanning both sides is what the ceiling prefers.
3. MORE SEEDS BUY NOTHING under the registered rule, so this block does not spend on them.
   Three seeds, as every binary cell on disk has.

THE RISK THIS DESIGN IS SIZED AGAINST, and it is measured rather than assumed
-----------------------------------------------------------------------------
The ceiling is computed AT the shipped extrapolated a(u), which is the number the block
exists to replace. The r = 1.25 data is HELD OUT of the closure (it failed
identifiability, so `tabulated_closure` never put a knot there), which makes it a legal
test of that extrapolation -- and it fails it in sign on 8 of 8 compositions, mean
-1.03 scale: the real pack at u = 0.8 fills LESS than the PCHIP says, i.e. the true a is
CLOSER TO 1. That is the direction that kills identifiability. The robustness table in
`diag_uident.txt` section 6 puts the OLD grid's ceiling limit at a(u*) = 0.960 and the
WIDE grid's at 0.980, so the wider grid is chosen for a reason that was measured before
the choice, not after.

Tag `ug` / `ug3k`. It matches NONE of `mixr*`, `extr*`, `tern*`, `kl*`, so M0-M6, F1-F3,
X0-X5, T0-T6, NA1-NA4 and G0-G8 cannot be re-decided by anything in this block;
`grade_ugap.py` asserts that structurally. The runner is `phase12_mixture.py` UNCHANGED,
driven by its existing --ratio/--xgrid/--tag options: a new binary block needs no new
runner, and reusing the shipped one is what makes the bit-identity gate below possible.
"""

from __future__ import annotations

import numpy as np

from psd import BinaryPSD

# ---------------------------------------------------------------- protocol constants
P_HAT = 1e-3                 # the rung every binary, ternary and K-ladder cell has
N_PROD = 1500                # the N the tabulated closure is calibrated at. NOT optional:
                             # mixture_closure.binary_points drops any other N, so a
                             # knot measured at 3000 could never enter the closure.
N_BIG = 3000
MU_P = 0.3
SEEDS = (1, 2, 3)
MIN_PHASE_COUNT = 50         # inherited verbatim from phase12_mixture.py
PHI_SCALE = 2.1e-3           # the project's scale unit, unchanged

# The identifiability rule, quoted VERBATIM from grade_mixture_ext.tabulated_closure.
# NOT MOVABLE. A block designed to meet a bar must not be graded against a bar it moved.
IDENT_WIDTH = 0.25

# ---------------------------------------------------------------------------- the grid
U_STAR = 4.0 ** -0.25                    # 0.70711 -- K = 5's adjacent pair on span 4
R_STAR = 4.0 ** 0.25                     # 1.41421

# Arm A is the ratio the K-ladder reads. Arm B is the far end of the gap, and is the
# direct retest of the ratio that came back UNDETERMINED. Arm C is the CONTROL: r = 1.5
# is an EXISTING identifiable knot, so running the new composition grid there separates
# "the grid changed the answer" from "the ratio changed the answer". Without arm C a
# disagreement between the new a(0.7071) and the PCHIP could not be attributed.
RATIOS = {
    "A": R_STAR,        # u = 0.70711 -- the knot the K-ladder needs
    "B": 1.25,          # u = 0.80000 -- the ratio that failed
    "C": 1.5,           # u = 0.66667 -- the control: a knot that already exists
}

# The composition grid. Four points below the branch crossover (x* ~ 0.356), five above,
# so BOTH parameters get their own half of the curve. 0.05 and 0.20 are kept because
# r = 1.25 and r = 1.5 already have those exact cells on disk under tag `mixr`, at the
# same seeds, N, p_hat and mu_p, and `phase12_mixture.cell_key` carries no tag -- so
# those six cells must come back BIT-IDENTICAL. 0.88 is the largest composition
# affordable at MIN_PHASE_COUNT on ALL THREE ratios at N = 1500 (the binding one is
# r = 1.5, where x = 0.896 is the limit).
X_GRID = (0.05, 0.12, 0.20, 0.28, 0.36, 0.45, 0.55, 0.70, 0.88)

# The subset that must reproduce bit-identically against the `mixr` block on disk.
BITID = tuple((r, x) for r in (1.25, 1.5) for x in (0.05, 0.20))


def n_min_frac(ratio: float, x: float) -> float:
    """Minority NUMBER fraction of a binary at size ratio `ratio`, fines volume `x`."""
    p = BinaryPSD(d10=1.0, ratio=float(ratio), x_small=float(x))
    return float(min(p.n_small, 1.0 - p.n_small))


def crossover(beta: float, a: float, b: float) -> float:
    """The composition where the binding branch swaps, solved in closed form.

    The two branch totals are (1 - a) x and (1 - beta)(1 - b)(1 - x); they are equal at
    x* = c / (c + (1 - a)) with c = (1 - beta)(1 - b). Below x* the COARSE branch binds
    and the model sees only `a`; above it the FINE branch binds and it sees only `b`.
    Returns nan when neither is positive (a = b = 1, where phi = beta everywhere).
    """
    c = (1.0 - beta) * (1.0 - b)
    den = c + (1.0 - a)
    return float(c / den) if den > 0 else float("nan")


def phi_fast(a, b, x, beta: float):
    """Vectorised binary closure, from the branch algebra above.

    A SECOND implementation of `mixture_closure.packing_density_ab` for the binary case,
    not a wrapper -- `check_kernel` asserts it reproduces the shipped one before any
    number computed with it is printed. It exists because the registered profile grid is
    201 x 201 re-minimisations per parameter and the shipped pure-Python scan of it does
    not finish inside an agent call.
    """
    tot = np.minimum((1.0 - a) * x, (1.0 - beta) * (1.0 - b) * (1.0 - x))
    return beta / (1.0 - tot)


def check_kernel(beta: float, rng, ntrial: int = 300) -> float:
    """Worst |phi_fast - mixture_closure.packing_density_ab| over random (r, a, b, x)."""
    import mixture_closure as mc            # noqa: PLC0415 - avoids an import cycle

    worst = 0.0
    for _ in range(ntrial):
        a, b, x = rng.uniform(1e-4, 1 - 1e-4, 3)
        r = float(rng.uniform(1.05, 6.0))
        ref = mc.packing_density_ab(np.array([r, 1.0]), np.array([1.0 - x, x]),
                                    beta, float(a), float(b))
        worst = max(worst, abs(float(phi_fast(a, b, x, beta)) - ref))
    return worst


def profile_fast(pts, beta: float, which: int, ngrid: int = 201):
    """`mixture_closure.profile_interval`, same grid, same delta chi2 = n bar, vectorised."""
    x = np.array([p[1] for p in pts], dtype=float)
    y = np.array([p[2] for p in pts], dtype=float)
    grid = np.linspace(1e-4, 1.0 - 1e-4, ngrid)
    chi = (((phi_fast(grid[:, None, None], grid[None, :, None], x[None, None, :], beta)
             - y[None, None, :]) / PHI_SCALE) ** 2).sum(axis=2)
    prof = chi.min(axis=1) if which == 0 else chi.min(axis=0)
    lo = float(prof.min())
    ok = grid[prof <= lo + len(pts)]
    return (float(ok.min()), float(ok.max())) if ok.size else (0.0, 1.0)


def identifiable(pts, beta: float) -> tuple[bool, tuple, tuple, float, float]:
    """The registered rule, applied. Returns (ok, a_interval, b_interval, a_hat, b_hat).

    Quoted verbatim from grade_mixture_ext.tabulated_closure: both profile intervals
    narrower than IDENT_WIDTH *and* containing their own best fit. The second clause is
    not decoration -- it is what caught the r = 5 cusp.
    """
    import mixture_closure as mc            # noqa: PLC0415

    a_hat, b_hat, _ = mc.fit_ab_direct(pts, beta)
    ia = profile_fast(pts, beta, 0)
    ib = profile_fast(pts, beta, 1)
    ok = (ia[1] - ia[0] < IDENT_WIDTH and ib[1] - ib[0] < IDENT_WIDTH
          and ia[0] - 1e-9 <= a_hat <= ia[1] + 1e-9
          and ib[0] - 1e-9 <= b_hat <= ib[1] + 1e-9)
    return bool(ok), ia, ib, float(a_hat), float(b_hat)


def cells(n: int = N_PROD) -> list[dict]:
    """Every cell of the block, ADMISSIBILITY APPLIED HERE, judged at N_PROD always.

    Judging at N_PROD whatever `n` is asked for keeps the N = 1500 and N = 3000 arms on
    the IDENTICAL cell list, so the N comparison is a comparison of the same points
    (lessons.md 2026-09-08). A cell the coarse class cannot afford is REFUSED and
    recorded, never silently dropped.
    """
    out: list[dict] = []
    for arm, r in RATIOS.items():
        for x in X_GRID:
            rec = {"arm": arm, "ratio": float(r), "u": 1.0 / float(r), "x": float(x),
                   "n": int(n), "seeds": SEEDS,
                   "n_min_at_judge": n_min_frac(r, x) * N_PROD,
                   "label": f"ug{arm}_r{r:.5f}_x{x:.2f}",
                   "bitid": (round(float(r), 6), round(float(x), 6))
                            in {(round(a_, 6), round(b_, 6)) for a_, b_ in BITID}}
            if rec["n_min_at_judge"] < MIN_PHASE_COUNT:
                rec["why_refused"] = (f"min phase count {rec['n_min_at_judge']:.0f} < "
                                      f"{MIN_PHASE_COUNT} at N = {N_PROD}")
            out.append(rec)
    return out


def admitted(n: int = N_PROD) -> list[dict]:
    return [c for c in cells(n) if "why_refused" not in c]


def xgrid_arg() -> str:
    return ",".join(repr(float(x)) for x in X_GRID)


def summary() -> str:
    adm = admitted()
    ref = [c for c in cells() if "why_refused" in c]
    return (f"U-GAP block: {len(RATIOS)} ratios x {len(X_GRID)} compositions = "
            f"{len(cells())} cells, {len(adm)} admitted, {len(ref)} refused; "
            f"{len(adm) * len(SEEDS)} runs at N = {N_PROD}")


if __name__ == "__main__":
    print(summary())
    print(f"u* = {U_STAR:.6f}  r* = {R_STAR!r}")
    print(f"--xgrid {xgrid_arg()}")
    for c in cells():
        flag = c.get("why_refused", "BITID" if c["bitid"] else "")
        print(f"   {c['arm']}  r {c['ratio']:.5f}  u {c['u']:.5f}  x {c['x']:.2f}  "
              f"n_min {c['n_min_at_judge']:6.0f}  {flag}")
    print("beta-free check -- crossover at a = 0.94, b = 0.92, beta = 0.5947: "
          f"{crossover(0.5947, 0.94, 0.92):.4f}")
    print(f"np {np.__version__}")
