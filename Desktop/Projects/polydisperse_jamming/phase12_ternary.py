"""phase12_ternary.py -- is the pairwise-filling functional ADDITIVE?

PRE-REGISTRATION: PREREG_TERNARY_2026-09-14.md. Grade with grade_ternary.py.

THE QUESTION
------------
`grade_mixture.py` and `grade_mixture_ext.py` established that the control parameter for
the jamming shift is a PAIRWISE-FILLING FUNCTIONAL of the size distribution rather than
any descriptor of it, and calibrated the two interaction functions `a(u)` and `b(u)` on
BINARY packs at size ratios 1.25, 2, 3, 4, 5. The closure built from those measured
interactions then transfers to the continuous families WORSE than the two-exponent form
it was meant to improve on -- held-out chi2 108.1 against 34.9, under-predicting the
lognormal monotonically with CV to -6.9 scale at CV 0.5.

The linear packing model sums an interaction over every ORDERED PAIR of classes. A
binary pack has exactly one pair, so no binary campaign, at any ratio, can test whether
that SUM is the right way to combine them. The under-prediction on the continuous
families is precisely what a NON-ADDITIVE functional evaluated additively would look
like -- and it is also what a mis-calibrated `a(u)` at small `u` would look like. The two
are confounded in every measurement this project has.

A TERNARY separates them. Three classes at d = 1, 2, 4 give three pairs at
u = d_fine/d_coarse of 1/2 (1-2), 1/2 (2-4) and 1/4 (1-4). Both values sit inside the
range the binaries measure directly, so the interaction VALUES are not extrapolated: the
only new thing being asked of the model is whether three simultaneous pairs combine the
way it says they do. If the ternary is predicted well and the continuous families are
not, additivity holds and the residual is a calibration problem at small `u`. If the
ternary is under-predicted in the same direction and of a size the closure cannot absorb,
additivity is the defect.

WHY d = 1, 2, 4 AND NOT A WIDER LADDER
--------------------------------------
Two constraints bracket it from both sides. The classes must be far enough apart for the
filling physics to be ON -- at u = 1/2 the measured `a` is 0.785, so a pair at ratio 2
still interacts strongly -- and close enough together that the coarse class has enough
grains to be a phase rather than a sampling fluctuation. Number fraction goes as
y_i/d_i^3, so at d = 4 a 60 % volume fraction is only 5 % of the grains: 74 of 1500.
`MIN_PHASE_COUNT` = 50 is inherited verbatim from `phase12_mixture.py` and it is what
excludes (0.50, 0.30, 0.20), which would give n_L = 48.

d = 1, 2, 4 is also the only ladder on which the two pair ratios are BOTH inside the
binary calibration set AND distinct, which is what makes the additivity question sharp
rather than an extrapolation dressed up as one.

WHAT IS AND IS NOT NEW CODE
---------------------------
This file is `phase12_mixture.py` with the PSD swapped and the record fields renamed. The
protocol -- `build_and_compress` at p_hat = 1e-3, N = 1500, mu_p = 0.3, seeds 1/2/3, the
same `contact_graph`/`rattler_by_class` structure read off the live scene -- is byte-for-
byte the production one, so the monodisperse control must land on the same `beta` the
binary campaign measured. That control is gate T0 and it is what proves this block is on
the same protocol rather than merely near it.

The tag is `tern`, chosen NOT to match the `phase12_mixr*` or `phase12_extr*` globs, so
M0-M6, F1-F3 and X0-X5 are CLOSED and this campaign cannot re-decide any of them.
(lessons.md 2026-09-06: appending to a shared list is an edit to every gate that
iterates it.)

Usage:
    python3 phase12_ternary.py --shard 0 --nshards 12 --tag tern
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time

from measure import rattler_by_class
from packing import build_and_compress, contact_graph
from phase12_mixture import MIN_PHASE_COUNT
from psd import DiscretePSD

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

P_HAT = 1e-3          # the rung every family and binary record on disk already has
N_PROD = 1500
MU_P = 0.3
SEEDS = (1, 2, 3)

# The diameter ladder. Relative, exact integers: the family string carries them with
# full round-trip precision and the pair ratios are exactly 1/2, 1/2 and 1/4.
D_REL = (1.0, 2.0, 4.0)

# The compositions, as (y_L, y_M, y_S) -- VOLUME fractions of the d = 4, 2, 1 classes.
# Registered in PREREG_TERNARY_2026-09-14.md section 2. All six satisfy n_L >= 50 at
# N = 1500; (0.50, 0.30, 0.20) gives 48 and is REFUSED rather than run and dropped.
COMPOSITIONS = (
    (0.60, 0.25, 0.15),
    (0.55, 0.30, 0.15),
    (0.55, 0.35, 0.10),
    (0.65, 0.20, 0.15),
    (0.70, 0.20, 0.10),
    (0.60, 0.30, 0.10),
)

# T0, the protocol control: a monodisperse pack through the SAME DiscretePSD code path.
# It is not a spare cell -- if the discrete path does not reproduce beta, every number in
# the block is measuring the code rather than the composition.
MONO = (1.0, 0.0, 0.0)


def y_string(y_lms: tuple[float, float, float]) -> str:
    """(y_L, y_M, y_S) -> the family string, ASCENDING in diameter.

    The design tuple is written coarse-first because that is how the compositions were
    registered and how a bench reader thinks about them; `DiscretePSD` pairs its two
    sequences position by position. Reversing here, ONCE, in a named function, is the
    same defence `psd.radii_for_yade` uses for the factor of two.
    """
    y_l, y_m, y_s = y_lms
    return (f"discrete:d{D_REL[0]!r}/{D_REL[1]!r}/{D_REL[2]!r}"
            f":y{y_s!r}/{y_m!r}/{y_l!r}")


def make_ternary(y_lms: tuple[float, float, float]) -> DiscretePSD:
    y_l, y_m, y_s = y_lms
    return DiscretePSD(d10=1.0, d_rel=D_REL, y=(y_s, y_m, y_l))


def build_grid(n: int, with_mono: bool, comps=None):
    comps = list(COMPOSITIONS if comps is None else comps)
    if with_mono:
        comps = [MONO] + comps
    return [(c, seed, n) for c in comps for seed in SEEDS]


def cell_key(y_lms, seed: int, n: int, mu_p: float) -> str:
    y_l, y_m, y_s = y_lms
    return (f"tern_L{y_l:.4f}_M{y_m:.4f}_S{y_s:.4f}"
            f"_s{seed}_p{P_HAT:.0e}_N{n}_mu{mu_p:.2f}")


def load_done(path: str) -> set[str]:
    done: set[str] = set()
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                try:
                    done.add(json.loads(line)["key"])
                except Exception:      # noqa: BLE001 - a truncated last line is not fatal
                    continue
    return done


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="tern")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--n", type=int, default=N_PROD)
    ap.add_argument("--no-mono", action="store_true")
    ap.add_argument("--mu", type=float, default=MU_P)
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    grid = build_grid(args.n, not args.no_mono)
    suffix = f"_{args.tag}"
    if args.nshards > 1:
        suffix += f"_shard{args.shard:02d}of{args.nshards:02d}"
    path = os.path.join(OUT_DIR, f"phase12{suffix}.jsonl")
    done = load_done(path)
    mine = sum(1 for i in range(len(grid)) if i % args.nshards == args.shard)
    print(f"phase12 TERNARY: d {D_REL}, mu_p {args.mu}, {len(grid)} cells, {mine} on "
          f"shard {args.shard}/{args.nshards}, N={args.n}, {len(done)} done -> {path}",
          flush=True)

    t0 = time.time()
    nrun = 0
    for idx, (y_lms, seed, ncell) in enumerate(grid):
        if idx % args.nshards != args.shard:
            continue
        key = cell_key(y_lms, seed, ncell, args.mu)
        if key in done:
            continue
        fam = y_string(y_lms)
        psd = make_ternary(y_lms)
        try:
            st = build_and_compress(cv=0.0, mu_p=args.mu, p_hat=P_HAT, n=ncell, seed=seed,
                                    family=fam, truncate=None, verbose=True)
            rec = json.loads(st.to_json())
            # Size-resolved structure off the SAME live scene, exactly as the binary
            # runner does. The split is the geometric mean of the two EXTREME design
            # diameters, so for a ternary it puts the middle class on whichever side its
            # own diameter falls -- sqrt(1*4) = 2 is the middle class's own diameter, so
            # the split is recorded and the middle class is reported separately below
            # rather than being silently assigned. `by_class` is kept for continuity with
            # the binary block; `by_design_class` is the unambiguous one.
            from yade import O                    # noqa: PLC0415 - WSL only
            counts, neigh, _ = contact_graph(O)
            diam = {b.id: 2.0 * float(b.shape.radius) for b in O.bodies if b is not None}
            dd = psd.diameters
            rec["by_class"] = rattler_by_class(
                counts, neigh, diam, d_split=math.sqrt(float(dd[0]) * float(dd[-1])))
            # Per-DESIGN-class counts and rattler fractions: a ternary has three classes
            # and a two-way split cannot report them. Grains are matched to the design
            # diameter they were built from, which is exact (the PSD is discrete).
            per = {}
            for ci, dc in enumerate(dd):
                ids = [i for i, dv in diam.items() if abs(dv - float(dc)) < 1e-12]
                zs = [counts.get(i, 0) for i in ids]
                per[f"c{ci}"] = {
                    "d": float(dc),
                    "n": len(ids),
                    "z_mean": (float(sum(zs)) / len(zs)) if ids else float("nan"),
                    "rattler_frac": (sum(1 for z in zs if z < 4) / len(ids)) if ids else float("nan"),
                }
            rec["by_design_class"] = per
            rec["key"] = key
            rec["tern_y_L"] = y_lms[0]
            rec["tern_y_M"] = y_lms[1]
            rec["tern_y_S"] = y_lms[2]
            rec["d_rel"] = list(D_REL)
            rec["n_frac_design"] = [float(v) for v in psd.number_fractions]
            rec["d_design"] = [float(v) for v in psd.diameters]
            rec["cv_design"] = psd.cv
            rec["d32_design"] = psd.d32
            rec["d43_design"] = psd.d43
            rec["size_ratio_design"] = psd.size_ratio
            rec["n_design_counts"] = [int(v) for v in psd.class_counts(ncell)]
            rec["is_mono"] = bool(y_lms == MONO)
            rec["min_phase_count"] = MIN_PHASE_COUNT
            rec["campaign"] = "ternary-2026-09-14"
        except Exception as exc:            # noqa: BLE001 - record, never hide
            rec = {"key": key, "tern_y_L": y_lms[0], "tern_y_M": y_lms[1],
                   "tern_y_S": y_lms[2], "seed": seed, "p_hat_target": P_HAT,
                   "family": fam, "error": repr(exc),
                   "campaign": "ternary-2026-09-14"}
            print(f"CELL FAILED {key}: {exc!r}", flush=True)
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(rec, sort_keys=True) + "\n")
        nrun += 1
    print(f"phase12 ternary shard {args.shard}: ran {nrun} cells in "
          f"{time.time() - t0:.0f} s -> {path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
