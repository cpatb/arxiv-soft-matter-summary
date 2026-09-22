"""grade_shapechi.py -- PREREG_SHAPECHI_2026-09-22.md.

Is SHAPE a WINDOW-knob (moves chi = phim/phi0) or a POSITION-knob (moves x = phi/phi0)?

  W1  WIRING (PRIMARY, gates everything). Three existing N = 1700 records re-run here
      must come back BIT-IDENTICAL.
  W2  HEALTH. 72 of 72 converge at 1e6.
  S1  PRIMARY. chi(alpha = 2) against chi(alpha = 1) at the matched rung, seed-paired,
      converged only, against the protocol floor. Three registered outcomes.
  S2  REPORTED: phi0(alpha) and phim(alpha) separately, so a Delta chi can be read as
      WHICH BRANCH moved.
  S3  REPORTED: the same verdict at all four rungs. chi is NOT rung-independent, so a
      verdict that flips between rungs is rung-dependent, not resolved.

Run:  python3 grade_shapechi.py     (writes out/grade_shapechi.txt)
"""

from __future__ import annotations

import glob
import json
import math
import os
import sys

import numpy as np

import mixture_model as mm

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
SCALE = mm.PHI_SCALE
TAG = "shapechi"
N_REQ = 1700
CV = 0.0
ALPHAS = (1.0, 2.0)
MU_0, MU_M = 0.0, 3.0
RUNGS = (1e-3, 3e-4, 1e-4, 3e-5)
QUOTE_RUNG = 3e-5
N_REGISTERED = 72
T_CRIT = {2: 12.706, 3: 4.303}

# W1's three registered reproductions, quoted in the prereg before launch.
W1 = {(1.0, 1e-3): 0.595694, (1.0, 3e-5): 0.575051, (2.0, 1e-3): 0.566148}
W1_SEED = 1

FAIL: list[str] = []
L: list[str] = []


def say(s: str = "") -> None:
    print(s)
    L.append(s)


def gate(name: str, ok: bool, detail: str = "") -> None:
    say(f"   [{'PASS' if ok else 'FAIL'}] {name}" + (f"  --  {detail}" if detail else ""))
    if not ok:
        FAIL.append(name)


def load(tag: str | None = TAG) -> dict:
    """(alpha, mu_p, seed, rung) -> record. Searches out/ AND out/shape/."""
    cells: dict = {}
    pats = [os.path.join(OUT, "phase2*.jsonl"),
            os.path.join(OUT, "shape", "*.jsonl")]
    for pat in pats:
        for fn in sorted(glob.glob(pat)):
            if tag is not None and tag not in os.path.basename(fn):
                continue
            for ln in open(fn, encoding="utf-8"):
                ln = ln.strip()
                if not ln:
                    continue
                r = json.loads(ln)
                if r.get("failed") or "phi" not in r:
                    continue
                if r.get("n_requested") != N_REQ or abs(r.get("cv", 1.0) - CV) > 1e-12:
                    continue
                k = (round(r.get("alpha", 1.0), 3), r["mu_p"], r["seed"],
                     r["p_hat_target"])
                prev = cells.get(k)
                if prev is None or r.get("max_iter", 0) > prev.get("max_iter", 0):
                    cells[k] = r
    return cells


def chi_at(cells: dict, alpha: float, rung: float):
    """chi = phi(MU_M)/phi(MU_0) at one alpha and rung. SEED-PAIRED.

    Paired because `packing.py` seeds the radius/orientation draw from the seed
    regardless of `mu_p`. The floor propagates SCALE through the ratio and does NOT
    cancel under pairing.
    """
    a, b = {}, {}
    for s in (1, 2, 3):
        ra, rb = cells.get((alpha, MU_0, s, rung)), cells.get((alpha, MU_M, s, rung))
        if ra is None or rb is None:
            continue
        if not (ra.get("converged") and rb.get("converged")):
            continue
        a[s], b[s] = ra["phi"], rb["phi"]
    ss = sorted(a)
    if not ss:
        return None, None, 0, None, None
    ch = np.array([b[s] / a[s] for s in ss])
    m = float(ch.mean())
    p0 = float(np.mean([a[s] for s in ss]))
    pm = float(np.mean([b[s] for s in ss]))
    se_seed = float(ch.std(ddof=1) / math.sqrt(len(ss))) if len(ss) > 1 else float("inf")
    se_floor = SCALE * math.sqrt(1 + m * m) / (p0 * math.sqrt(len(ss)))
    return m, max(se_seed, se_floor), len(ss), p0, pm


def main() -> int:
    say("=" * 96)
    say("grade_shapechi.py -- PREREG_SHAPECHI_2026-09-22.md")
    say("=" * 96)
    say("   chi = phim/phi0 sets the SHEAR-JAMMING WINDOW (width 1-chi).")
    say("   Friction moves chi by construction; the size DISTRIBUTION does not")
    say("   (|t| = 0.40). Shape is the remaining candidate window-knob, and the shape")
    say("   campaign ran at ONE friction so it could not form chi at all.")

    cells = load()
    say(f"\n   `{TAG}` records on disk: {len(cells)}  (registered {N_REGISTERED})")
    if not cells:
        say("   🔴 NOTHING DELIVERED. Not graded.")
        return 1
    caps = sorted(set(str(r.get("max_iter")) for r in cells.values()))
    say(f"   caps present: {caps}")

    # ---------------------------------------------------------------- W1
    say("\n=== W1  WIRING (PRIMARY): the three pre-existing records must reproduce ===")
    ref = load(tag=None)          # every N=1700 CV=0 record, any tag
    n_cmp = n_bad = 0
    for (al, rung), want in sorted(W1.items()):
        got = cells.get((al, 0.3, W1_SEED, rung))
        if got is None:
            say(f"      alpha={al}, p={rung:.0e}: NOT YET RUN")
            continue
        n_cmp += 1
        old = ref.get((al, 0.3, W1_SEED, rung))
        d = abs(got["phi"] - want)
        traj = ""
        if old is not None and old is not got:
            same = all(old.get(f) == got.get(f) for f in ("iters", "n_contacts")
                       if f in old and f in got)
            traj = f"   trajectory identical: {same}"
        ok = d == 0.0
        n_bad += (not ok)
        say(f"      alpha={al}, p={rung:.0e}: got {got['phi']:.9f} vs registered "
            f"{want:.6f}  |d| = {d:.2e}{traj}")
    gate("W1: the three registered records reproduce EXACTLY",
         n_cmp > 0 and n_bad == 0, f"{n_cmp - n_bad}/{n_cmp} exact")
    if n_cmp and n_bad:
        say("      ⚠ If the only difference is ~1e-14 in `phi` while `iters` and")
        say("      `n_contacts` match, this is the cross-code-era effect this session")
        say("      characterised, NOT a wiring fault -- but it is reported, not excused.")

    # ---------------------------------------------------------------- W2
    say("\n=== W2  HEALTH ===")
    nconv = sum(1 for r in cells.values() if r.get("converged"))
    gate("W2: every delivered record converged", nconv == len(cells),
         f"{nconv}/{len(cells)}")
    bad = [k for k, r in cells.items() if not r.get("converged")]
    for k in bad[:8]:
        say(f"      🔴 capped: alpha={k[0]}, mu={k[1]}, seed={k[2]}, p={k[3]:.0e}")
    if bad:
        say("      (reported, NOT re-raised -- the prereg says so)")

    # ---------------------------------------------------------------- S2
    say("\n=== S2  THE TWO BRANCHES SEPARATELY (reported, not gated) ===")
    say(f"   at the quoted rung p_hat = {QUOTE_RUNG:.0e}")
    say(f"   {'alpha':>7}{'n':>4}{'phi_0':>10}{'phi_m':>10}{'chi':>9}{'se':>9}"
        f"{'window 1-chi':>14}")
    q: dict = {}
    for al in ALPHAS:
        m, se, n, p0, pm = chi_at(cells, al, QUOTE_RUNG)
        if m is None:
            say(f"   {al:7.1f}   -- no converged pair --")
            continue
        q[al] = (m, se, n)
        say(f"   {al:7.1f}{n:4d}{p0:10.5f}{pm:10.5f}{m:9.4f}{se:9.4f}{1-m:14.4f}")
    say("")
    say("   ⚠ alpha is NON-MONOTONE on the records already on disk (phi = 0.5957 /")
    say("   0.6289 / 0.5661 at alpha = 1.0 / 1.5 / 2.0, mu_p = 0.3, p = 1e-3), so these")
    say("   two rungs show a DIRECTION and not a shape. No curvature claim.")

    # ---------------------------------------------------------------- S1
    say("\n=== S1  PRIMARY: does SHAPE move the window? ===")
    if len(q) < 2:
        say("   🔴 NOT DECIDABLE: fewer than two alpha values have a converged chi.")
        FAIL.append("S1: both alphas measurable")
    else:
        (m1, s1, n1), (m2, s2, n2) = q[ALPHAS[0]], q[ALPHAS[1]]
        d = m2 - m1
        se = math.hypot(s1, s2)
        t = abs(d) / se if se > 0 else float("inf")
        tc = T_CRIT.get(min(n1, n2), 4.303)
        say(f"   chi(alpha=2) - chi(alpha=1) = {d:+.4f} +- {se:.4f}   |t| = {t:.2f} "
            f"vs t_crit {tc} (n = {min(n1, n2)})")
        say("")
        if t > tc and d < 0:
            say("   🔑 VERDICT (W) WINDOW-KNOB: shape WIDENS the shear-jamming window.")
            say("   Elongation acts like friction -- it lowers phi_m more than phi_0.")
            say(f"   The window goes {1-m1:.4f} -> {1-m2:.4f}, a factor {(1-m2)/(1-m1):.2f}.")
        elif t > tc and d > 0:
            say("   🔴 VERDICT (N) shape NARROWS the window -- the rotational-constraint")
            say("   argument is REFUTED, and this was registered as the third outcome.")
        else:
            say("   🔑 VERDICT (P) POSITION-KNOB: no resolved move in chi. Shape joins the")
            say("   size distribution -- it moves phi_0 and phi_m together and changes")
            say("   only WHERE the suspension sits, not whether it can shear-jam.")
            say("   ⚠ This leaves FRICTION as the only measured window-knob, which is the")
            say("   more useful answer and the one that most simplifies the programme.")
            say("   ⚠ 'Not resolved' at n = 3 is a BOUND, never a proof of identity.")
        say("")
        say("   ⚠ A CLUMP MOVES THREE THINGS AT ONCE -- elongation, surface corrugation")
        say("   and contact multiplicity. This verdict is about NON-SPHERICITY, not about")
        say("   elongation; CLAIMS_SHAPE_2026-09-18 records that no clump-only run")
        say("   separates them, and the LevelSet control is not run here.")

    # ---------------------------------------------------------------- S3
    say("\n=== S3  IS THE VERDICT THE SAME AT EVERY RUNG? ===")
    say("   chi is NOT rung-independent (drift 0.016-0.025 across the ladder, larger")
    say("   than phi_0's own), so S1 is a fixed-rung comparison by construction.")
    say("")
    say(f"   {'rung':>10}{'chi(a=1)':>11}{'chi(a=2)':>11}{'d(chi)':>10}{'|t|':>8}  verdict")
    signs = []
    for p in RUNGS:
        a1 = chi_at(cells, ALPHAS[0], p)
        a2 = chi_at(cells, ALPHAS[1], p)
        if a1[0] is None or a2[0] is None:
            say(f"   {p:10.0e}{'--':>11}{'--':>11}")
            continue
        d = a2[0] - a1[0]
        se = math.hypot(a1[1], a2[1])
        t = abs(d) / se if se > 0 else float("inf")
        tc = T_CRIT.get(min(a1[2], a2[2]), 4.303)
        v = "resolved" if t > tc else "not resolved"
        signs.append((p, d, t > tc))
        say(f"   {p:10.0e}{a1[0]:11.4f}{a2[0]:11.4f}{d:+10.4f}{t:8.2f}  {v}")
    res = [s for s in signs if s[2]]
    if len(res) > 1 and len({np.sign(s[1]) for s in res}) > 1:
        say("")
        say("   🔴 THE SIGN FLIPS BETWEEN RUNGS on resolved comparisons -- the ladder")
        say("   systematic is larger than the effect and S1's verdict is RUNG-DEPENDENT.")
        FAIL.append("S3: the verdict does not flip sign between rungs")
    elif signs:
        say("")
        say(f"   the sign of d(chi) is {'consistent' if len({np.sign(s[1]) for s in signs}) == 1 else 'MIXED'} "
            f"across the {len(signs)} rungs with data "
            f"({len(res)} resolved).")

    say("\n" + "=" * 96)
    if FAIL:
        say(f"FAILED {len(FAIL)} gate(s):")
        for f in FAIL:
            say(f"   - {f}")
    else:
        say("ALL GATES PASSED")
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "grade_shapechi.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
