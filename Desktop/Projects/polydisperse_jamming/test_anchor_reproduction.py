"""test_anchor_reproduction.py — the paper's frictionless anchor, re-run independently.

Added 2026-09-22. `\\phiJFrictionless` in `paper/numbers.tex` is built by
`extract_numbers.py` from the `mu_p = 0, CV = 0, lognormal, N = 1500` cells. The
`mustarcap` campaign re-ran exactly those cells, months after the originals and on the
current code, without being designed to — which makes them a **free reproducibility
measurement on a published number** (the 2026-09-19 lesson: design campaigns so they
re-run a few already-measured cells).

Two things are pinned here, and they are different questions:

  1. REPRODUCTION — the re-run agrees with the original to floating-point level.
  2. NON-TAKEOVER — the re-run does not silently become the record the paper reads.

(2) matters because `sweep_descriptive.load_fits` dedups with
`groupby([...])["max_iter"].idxmax()`, and `idxmax` returns the FIRST row at a tie. Both
campaigns ran at `max_iter = 1e6`, so **which one feeds the macro is decided by filename
sort order** — `capfix_lg0` before `mustarcap`. That is the "first in file order"
tie-break this repository has already flagged as arbitrary. The test does not object to
it; it pins that the macro has not MOVED, so a future campaign cannot change a published
number without a test going red.

Skipped with a printed reason if `mustarcap` has not landed.
"""

import glob
import json
import os

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

# The published value and its printed precision, from paper/numbers.tex.
PHI_J_FRICTIONLESS = 0.6324
PRINTED_PRECISION = 1e-4

# How far apart two independent runs of the SAME cell may be. This is not a physics
# tolerance: it is floating point. Measured max |relative| = 1.9e-14 over 39 cross-era
# converged pairs, so 1e-12 leaves ~50x headroom and would still catch any real change.
REPRO_RTOL = 1e-12


def _cells(family="lognormal", mu=0.0, cv=0.0, n_req=1500):
    """(seed, p_hat) -> {tag: record} for one (family, mu_p, CV) at N = n_req."""
    out = {}
    for fn in sorted(glob.glob(os.path.join(OUT, "phase2*.jsonl"))):
        tag = os.path.basename(fn).split("phase2_")[1].split("shard")[0].rstrip("_")
        for ln in open(fn, encoding="utf-8"):
            ln = ln.strip()
            if not ln:
                continue
            r = json.loads(ln)
            if r.get("failed") or r.get("n_requested") != n_req:
                continue
            if (r.get("family") or "lognormal") != family:
                continue
            if r.get("mu_p") != mu or abs(r.get("cv", 1.0) - cv) > 1e-12:
                continue
            out.setdefault((r["seed"], r["p_hat_target"]), {})[tag] = r
    return out


def test_the_papers_frictionless_anchor_cells_reproduce():
    """Every anchor cell re-run by `mustarcap` must reproduce its original.

    ⚠ NOT asserted bit-identical, and that is a measurement rather than a concession:
    cross-era records in this repository differ at ~1e-14 in `phi` while `iters` and
    `n_contacts` are EXACTLY equal — the trajectory is the same and only the final
    reduction moves. The offset is NOT a single constant factor (checked: it ranges
    -1.93e-14 to +1.59e-15 over 39 converged pairs), so it is floating-point level
    variation in the reduction, not a changed constant and not non-determinism in the
    dynamics.
    """
    cells = _cells()
    pairs = [(k, d) for k, d in cells.items()
             if "mustarcap" in d and any(t != "mustarcap" for t in d)]
    if not pairs:
        pytest.skip("mustarcap has not landed; no anchor cell has been re-run")
    n_cmp, worst = 0, 0.0
    for k, d in pairs:
        new = d["mustarcap"]
        old_tag = next(t for t in d if t != "mustarcap")
        old = d[old_tag]
        if not (old.get("converged") and new.get("converged")):
            continue            # a capped cell genuinely moved; not a reproduction test
        n_cmp += 1
        rel = abs(new["phi"] - old["phi"]) / old["phi"]
        worst = max(worst, rel)
        assert rel < REPRO_RTOL, (k, old_tag, old["phi"], new["phi"], rel)
        # the TRAJECTORY must be identical, which is the stronger statement
        for f in ("iters", "n_contacts"):
            if f in old and f in new:
                assert old[f] == new[f], (k, f, old[f], new[f])
    assert n_cmp >= 12, f"only {n_cmp} converged anchor pairs available"
    # the check must not be vacuous: an EXACT zero everywhere would mean the two tags
    # are the same records read twice, not an independent repeat.
    assert worst > 0.0, "every pair identical to the bit -- is this really two campaigns?"


def test_the_published_macro_has_not_moved():
    """NON-TAKEOVER: the value `extract_numbers.py` would print is still the published one.

    A new campaign that re-runs already-published cells can silently replace them in the
    loader's dedup and move a number in the paper. This asserts it has not.
    """
    import sys
    sys.path.insert(0, HERE)
    import sweep_descriptive as sd

    f = sd.load_fits()
    g0 = f[(f.mu_p == 0.0) & (f.cv_nom == 0.0) & (f.family == "lognormal")].phi_J_lin
    assert len(g0) == 3, f"expected 3 seeds feeding the anchor, got {len(g0)}"
    assert float(g0.mean()) == pytest.approx(PHI_J_FRICTIONLESS,
                                             abs=PRINTED_PRECISION / 2), float(g0.mean())


def test_the_reproduction_check_can_actually_fail():
    """DISCRIMINATING: the tolerance must reject a difference that matters.

    Without this, `REPRO_RTOL = 1e-12` could be any number at all and the test above
    would pass for a trivial reason.
    """
    cells = _cells()
    if not cells:
        pytest.skip("no anchor cells on disk")
    phi = next(iter(cells.values()))
    v = next(iter(phi.values()))["phi"]
    perturbed = v * (1 + 1e-9)          # a difference the PRINTED precision would hide
    assert abs(perturbed - v) / v > REPRO_RTOL, "the tolerance is too loose to be a test"


def test_capped_cells_are_excluded_from_the_reproduction_claim():
    """A previously-CAPPED cell genuinely moved, and must not count as a reproduction.

    This is the trap the whole `mustarcap` campaign exists around: capped cells read LOW
    and the re-run moves them. Including them would make the reproduction tolerance
    absurd, or — worse — make a real physical shift look like floating-point noise.
    """
    cells = _cells(mu=3.0, cv=0.5)
    moved = []
    for k, d in cells.items():
        if "mustarcap" not in d:
            continue
        for t, old in d.items():
            if t == "mustarcap" or old.get("converged"):
                continue
            moved.append(abs(d["mustarcap"]["phi"] - old["phi"]) / old["phi"])
    if not moved:
        pytest.skip("no previously-capped mu_p=3.0, CV=0.5 cell has been re-run yet")
    assert max(moved) > 100 * REPRO_RTOL, (
        f"capped cells moved by at most {max(moved):.2e}, which is the same order as the "
        f"reproduction tolerance -- the two categories are not separable and the "
        f"reproduction claim above would be meaningless")
