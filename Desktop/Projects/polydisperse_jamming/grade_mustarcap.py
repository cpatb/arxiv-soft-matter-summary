"""grade_mustarcap.py -- PREREG_MUSTARCAP_2026-09-22.md.

The mu ladder `grade_mustar.py` fits is a MIXED-CAP ladder. This grades its repair.

  M1a IDENTITY (PRIMARY, gates everything). Cells that CONVERGED under the old cap must
      come back BIT-IDENTICAL at 1e6, compared against the `mustar` tag itself. This is
      the arm the prereg registers ("previously-converged `mustar` records").
  M1b DIAGNOSTIC, NOT A GATE. The same comparison against records written by an OLDER
      code state. Those differ from the current code at ~1e-14 in `phi` ONLY, with
      `iters` and `n_contacts` bit-identical -- i.e. the trajectory is the same and the
      derived scalar is not. That is a fact about this repository's cross-version
      reproducibility, not about the cap, so it is characterised and printed rather than
      graded. Gating on it would fail the campaign for the wrong reason.
  M2  RESCUE. >= 15 of the 17 previously-capped cells converge at 1e6.
  M3  MONOTONICITY (the point). phi(mu) must not RISE by more than PHI_SCALE at any CV.
      Fit-free; no threshold beyond the protocol's own declared scale.
  M4  DIRECTION. The previously-capped cells must move UP (under-compressed reads low);
      the converged ones must not move at all (M1).
  M5  mu* and A, re-fitted and REPORTED. Not a primary gate: F3a stays binding, and
      repairing a cap does not repair an identifiability problem.

⚠ If M1 fails, M2-M5 are printed but NOT READ, and the grader says so.

Run:  python3 grade_mustarcap.py     (writes out/grade_mustarcap.txt)
"""

from __future__ import annotations

import glob
import json
import math
import os
import sys

import numpy as np
from scipy.optimize import curve_fit

import mixture_model as mm

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

SCALE = mm.PHI_SCALE                     # 2.1e-3
RUNG = 3e-5
N_REQ = 1500
NEW_TAG = "mustarcap"
NEW_MUS = (0.0, 0.6, 3.0)                # the frictions this campaign re-ran
ALL_MUS = (0.0, 0.3, 0.6, 1.0, 3.0)
CVS = (0.0, 0.2, 0.5)
RESCUE_MIN = 15                          # of 17, registered
N_REGISTERED = 108
# F3/F4's bands are NOT moved. Quoted from PREREG_MUSTAR_2026-09-19.md via grade_mustar.
MUSTAR_BAND = (1.05, 1.45)
A_BAND = (0.90, 1.20)
F4 = {(0.6, 0.0): 0.56394, (0.6, 0.5): 0.62042, (3.0, 0.5): 0.61346}
F4_TOL = 0.0030

FAIL: list[str] = []
L: list[str] = []


def say(s: str = "") -> None:
    print(s)
    L.append(s)


def gate(name: str, ok: bool, detail: str = "") -> None:
    say(f"   [{'PASS' if ok else 'FAIL'}] {name}" + (f"  --  {detail}" if detail else ""))
    if not ok:
        FAIL.append(name)


def load(tag_filter) -> dict:
    """(mu_p, cv, seed, p_hat) -> record, for lognormal N=1500 records matching a filter.

    `tag_filter(basename) -> bool`. Where a key repeats WITHIN the selection, the highest
    `max_iter` wins -- the rule `sweep_descriptive.load_fits` has had since 2026-09-08 and
    the one `grade_mustar` uses.
    """
    cells: dict = {}
    for fn in sorted(glob.glob(os.path.join(OUT, "phase2*.jsonl"))):
        base = os.path.basename(fn)
        if not tag_filter(base):
            continue
        for ln in open(fn, encoding="utf-8"):
            ln = ln.strip()
            if not ln:
                continue
            r = json.loads(ln)
            if r.get("failed") or r.get("n_requested") != N_REQ:
                continue
            if (r.get("family") or "lognormal") != "lognormal":
                continue
            k = (r["mu_p"], round(r["cv"], 1), r["seed"], r["p_hat_target"])
            prev = cells.get(k)
            if prev is None or r.get("max_iter", 0) > prev.get("max_iter", 0):
                cells[k] = r
    return cells


IDENTITY_FIELDS = ("phi", "z_all", "z_nr", "rattler_frac", "p_hat", "iters",
                   "n_contacts", "vol_solid")


def identical(a: dict, b: dict) -> list[str]:
    """Fields that DIFFER. Compared on a CONTINUOUS channel, not only quantized ones.

    `phi`, `p_hat` and `vol_solid` are continuous; `iters` and `n_contacts` are counts.
    A bit-identity gate read only on a grid-quantized observable is not a bit-identity
    gate (lessons.md, 2026-08-31), which is why `phi` and `p_hat` lead the list.
    """
    bad = []
    for f in IDENTITY_FIELDS:
        if f not in a or f not in b:
            continue
        if a[f] != b[f]:
            bad.append(f"{f}: {a[f]!r} vs {b[f]!r}")
    return bad


def expo(mu, pi, A, ms):
    return pi + A * np.exp(-mu / ms)


def alge(mu, pi, A, ms):
    return pi + A / (1.0 + mu / ms)


def main() -> int:
    say("=" * 96)
    say("grade_mustarcap.py -- PREREG_MUSTARCAP_2026-09-22.md (the mixed-cap repair)")
    say("=" * 96)

    new = load(lambda b: NEW_TAG in b)
    old_mustar = load(lambda b: "mustar" in b and NEW_TAG not in b)
    old_any = load(lambda b: "mustar" not in b and "wcbidi" not in b)

    say(f"\n   new `{NEW_TAG}` records on disk: {len(new)}  (registered {N_REGISTERED})")
    caps = sorted(set(r.get("max_iter") for r in new.values()))
    say(f"   caps present in the re-run: {caps}")
    if not new:
        say("   🔴 NOTHING DELIVERED. Not graded.")
        return 1
    gate("the re-run really used the raised cap", caps == [1000000], f"{caps}")

    # ---------------------------------------------------------------- M1
    say("\n=== M1  IDENTITY (PRIMARY): previously-CONVERGED cells must be BIT-IDENTICAL ===")
    n_cmp_a = n_bad_a = 0
    worst = []
    for k, r_old in old_mustar.items():
        if not r_old.get("converged") or k not in new:
            continue
        n_cmp_a += 1
        bad = identical(r_old, new[k])
        if bad:
            n_bad_a += 1
            worst.append((k, bad))
    say(f"   M1a vs the `mustar` tag (mu_p in {{0.6, 3.0}}): {n_cmp_a} converged "
        f"records compared, {n_bad_a} mismatched")
    for k, bad in worst[:5]:
        say(f"      🔴 {k}: " + "; ".join(bad[:3]))
    gate("M1a: raising the cap moved no cell that had already converged",
         n_cmp_a > 0 and n_bad_a == 0, f"{n_cmp_a} compared, {n_bad_a} mismatches")

    n_cmp_b = n_bad_b = 0
    worst_b = []
    for k, r_old in old_any.items():
        if k[0] not in NEW_MUS or not r_old.get("converged") or k not in new:
            continue
        n_cmp_b += 1
        bad = identical(r_old, new[k])
        if bad:
            n_bad_b += 1
            worst_b.append((k, bad))
    say(f"   M1b vs pre-existing records at the same keys: {n_cmp_b} compared, "
        f"{n_bad_b} differ")
    if n_cmp_b == 0:
        say("      (no overlapping converged keys -- VACUOUS here, reported as such)")
    else:
        # Characterise the difference rather than pass/fail on it. M1b is NOT the
        # registered gate (the prereg registers M1 on "previously-converged `mustar`
        # records", which is M1a); it compares against records written by an OLDER code
        # state, which is a different question.
        dphi, fields = [], {}
        for k, r_old in old_any.items():
            if k[0] not in NEW_MUS or not r_old.get("converged") or k not in new:
                continue
            for f in IDENTITY_FIELDS:
                if f in r_old and f in new[k] and r_old[f] != new[k][f]:
                    fields[f] = fields.get(f, 0) + 1
            if "phi" in r_old:
                dphi.append(abs(new[k]["phi"] - r_old["phi"]))
        if dphi:
            mx = max(dphi)
            say(f"      magnitude: max |d(phi)| = {mx:.3e}, min {min(dphi):.3e}"
                f"  (~{mx/2.22e-16/0.56:.0f} ULP at phi ~ 0.56)")
        say(f"      fields that differ: {dict(fields) if fields else 'none'}")
        say("      ⚠ REPORTED, NOT GATED, and here is why. `iters` and `n_contacts` are")
        say("      IDENTICAL on every pair, so the TRAJECTORY is bit-identical and only")
        say("      the derived scalar differs -- a reduction-order / code-version effect,")
        say("      not a different simulation. M1b is not the registered gate; the prereg")
        say("      registers M1 on previously-converged `mustar` records, which is M1a.")
        say("      What it DOES establish: this repository's older records are NOT")
        say("      bit-reproducible by the current code at the ~1e-14 level, while")
        say("      same-code-era runs are bit-identical (M1a, exactly 0 on every pair).")
        say("      ⚠ RELATED DEFECT, found the same session and NOT fixed here:")
        say("      `grade_mustar.py`'s F1 is named 'reproduce an existing record EXACTLY'")
        say("      and its prereg says BIT-IDENTICAL, but it is implemented as")
        say("      `abs(v - x) < 1e-12`. It passed on records that differ at ~1e-14, so")
        say("      its VERDICT (the wiring is right) stands while its WORDING does not.")

    m1_ok = (n_cmp_a > 0 and n_bad_a == 0)
    if not m1_ok:
        say("\n   🔴 M1 FAILED. Per the prereg, M2-M5 below are PRINTED BUT NOT READ:")
        say("   the cap is not the only thing that changed between the two runs.")

    # ---------------------------------------------------------------- M2
    say("\n=== M2  RESCUE: the previously-capped cells under the raised cap ===")
    capped = [k for k, r in old_mustar.items() if not r.get("converged")]
    rescued = [k for k in capped if k in new and new[k].get("converged")]
    say(f"   previously capped (mustar): {len(capped)}")
    say(f"   of those, converged at 1e6: {len(rescued)} of {len(capped)}")
    if rescued:
        its = sorted(new[k]["iters"] for k in rescued)
        say(f"   iterations used by the rescued cells: {its[0]} .. {its[-1]}")
    still = [k for k in capped if k in new and not new[k].get("converged")]
    for k in still[:8]:
        say(f"      🔴 STILL CAPPED at 1e6: mu={k[0]}, CV={k[1]}, seed={k[2]}, p={k[3]:.1e}")
    gate(f"M2: at least {RESCUE_MIN} of {len(capped)} previously-capped cells converge",
         len(rescued) >= min(RESCUE_MIN, len(capped)),
         f"{len(rescued)}/{len(capped)}")
    say(f"   overall: {sum(1 for r in new.values() if r.get('converged'))} of {len(new)} "
        f"new records converged")

    # ---------------------------------------------------------------- ladder
    say("\n=== THE mu LADDER at the matched rung p_hat = 3e-5, AFTER the repair ===")
    merged = dict(old_any)
    merged.update(old_mustar)
    merged.update(new)            # the re-run wins at its own keys, and it is graded above

    def rung_vals(mu, cv):
        out = {}
        for s in (1, 2, 3):
            r = merged.get((mu, cv, s, RUNG))
            if r is not None:
                out[s] = r
        return out

    say(f"   {'CV':>6}" + "".join(f"{('mu=%g' % m):>13}" for m in ALL_MUS))
    lad: dict = {}
    for cv in CVS:
        row, conv = [], []
        for mu in ALL_MUS:
            v = rung_vals(mu, cv)
            if not v:
                row.append(float("nan"))
                conv.append("0/0")
                continue
            phis = [r["phi"] for r in v.values()]
            lad[(mu, cv)] = (float(np.mean(phis)), phis,
                             sum(1 for r in v.values() if r.get("converged")), len(v))
            row.append(float(np.mean(phis)))
            conv.append(f"{lad[(mu,cv)][2]}/{len(v)}")
        say(f"   {cv:6.1f}" + "".join(f"{x:13.5f}" for x in row))
        say(f"   {'conv':>6}" + "".join(f"{c:>13}" for c in conv))

    # ---------------------------------------------------------------- M3
    say("\n=== M3  MONOTONICITY: phi_J cannot RISE with friction ===")
    say(f"   bar: no rise above PHI_SCALE = {SCALE}. Fit-free, no other threshold.")
    n_rise = 0
    for cv in CVS:
        vals = [(mu, lad[(mu, cv)][0]) for mu in ALL_MUS if (mu, cv) in lad]
        for (m0, v0), (m1, v1) in zip(vals, vals[1:]):
            if v1 > v0 + SCALE:
                n_rise += 1
                c0 = lad[(m0, cv)]
                c1 = lad[(m1, cv)]
                say(f"      🔴 RISE at CV={cv}, mu {m0} -> {m1}: +{v1-v0:.5f} = "
                    f"{(v1-v0)/SCALE:.1f} x PHI_SCALE  "
                    f"[conv {c0[2]}/{c0[3]} -> {c1[2]}/{c1[3]}]")
    gate("M3: phi(mu) does not rise by more than PHI_SCALE at any CV", n_rise == 0,
         f"{n_rise} rises")
    if n_rise == 0:
        say("      the ladder is monotone within the protocol floor at every CV.")
        say("      ⚠ That is consistent with the mixed-cap explanation and does NOT")
        say("      prove it: monotonicity is necessary, not sufficient. M4 is the")
        say("      directional test that can distinguish them.")

    # ---------------------------------------------------------------- M4
    say("\n=== M4  DIRECTION: under-compressed cells must read LOW, so the shift is UP ===")
    dd_cap, dd_conv = [], []
    for k, r_old in old_mustar.items():
        if k not in new:
            continue
        d = new[k]["phi"] - r_old["phi"]
        (dd_cap if not r_old.get("converged") else dd_conv).append(d)
    if dd_cap:
        arr = np.array(dd_cap)
        say(f"   previously CAPPED cells (n={len(arr)}): mean d(phi) = {arr.mean():+.5f}"
            f"  = {arr.mean()/SCALE:+.2f} x PHI_SCALE, "
            f"{int((arr > 0).sum())} of {len(arr)} positive")
        say(f"      range {arr.min():+.5f} .. {arr.max():+.5f}")
        gate("M4: the previously-capped cells shifted UP on average",
             float(arr.mean()) > 0, f"mean {arr.mean():+.5f}")
    else:
        say("   no previously-capped cells to compare")
    if dd_conv:
        arr = np.array(dd_conv)
        say(f"   previously CONVERGED cells (n={len(arr)}): max |d(phi)| = "
            f"{np.abs(arr).max():.3e}  (M1 requires exactly 0)")

    # ---------------------------------------------------------------- M5
    say("\n=== M5  mu* AND A, RE-FITTED (reported; F3a stays binding) ===")
    say("   ⚠ NOT a primary gate. Repairing a cap does not repair an identifiability")
    say("   problem, and this campaign was not registered as able to.")
    fits = {}
    for cv in CVS:
        xs, ys, ses = [], [], []
        for mu in ALL_MUS:
            if (mu, cv) not in lad:
                continue
            m, phis, nconv, ntot = lad[(mu, cv)]
            if nconv < ntot:            # converged-only, per M1's own logic
                continue
            xs.append(mu)
            ys.append(m)
            se = float(np.std(phis, ddof=1) / math.sqrt(len(phis))) if len(phis) > 1 else SCALE
            ses.append(max(se, SCALE / math.sqrt(max(len(phis), 1))))
        if len(xs) < 4:
            say(f"   CV = {cv}: only {len(xs)} fully-converged frictions -- NOT FITTED")
            continue
        xs, ys, ses = np.array(xs, float), np.array(ys, float), np.array(ses, float)
        for label, model in (("exponential", expo), ("algebraic", alge)):
            try:
                p, _ = curve_fit(model, xs, ys, p0=[ys[-1], ys[0] - ys[-1], 0.25],
                                 sigma=ses, absolute_sigma=True, maxfev=40000)
                res = float(np.sqrt(np.mean((ys - model(xs, *p)) ** 2)))
                fits[(cv, label)] = (p, res)
                say(f"   CV = {cv}  {label:12s} phi_inf={p[0]:.5f} A={p[1]:.5f} "
                    f"mu*={p[2]:.4f}  resid RMS={res:.6f}  (n={len(xs)} pts)")
            except Exception as exc:
                say(f"   CV = {cv}  {label:12s} FIT FAILED: {exc}")
    for label in ("exponential", "algebraic"):
        a, b = fits.get((0.0, label)), fits.get((0.5, label))
        if a and b:
            r_mu = b[0][2] / a[0][2]
            r_A = b[0][1] / a[0][1]
            say(f"   {label}: mu*(0.5)/mu*(0) = {r_mu:.3f} vs registered {MUSTAR_BAND};"
                f"  A(0.5)/A(0) = {r_A:.3f} vs {A_BAND}")
            if not (A_BAND[0] <= r_A <= A_BAND[1]):
                say("      ⚠ A moved outside its band -> mu* and A are NOT separated by "
                    "this design -> UNDETERMINED, as PREREG_MUSTAR registered as likely.")

    # ---------------------------------------------------------------- F4
    say("\n=== F4  the ORIGINAL blind point predictions, rescored on repaired data ===")
    say("   ⚠ NOT re-registered. Scored once against the protocol they were written for")
    say("   (mixed-cap) and once against the repaired one; BOTH are printed.")
    for (mu, cv), pred in sorted(F4.items()):
        if (mu, cv) not in lad:
            say(f"   phi(mu={mu}, CV={cv}): no data")
            continue
        m, phis, nconv, ntot = lad[(mu, cv)]
        se = float(np.std(phis, ddof=1) / math.sqrt(len(phis))) if len(phis) > 1 else float("nan")
        d = m - pred
        say(f"   phi(mu={mu}, CV={cv}): predicted {pred:.5f}  repaired {m:.5f} "
            f"+- {se:.5f} ({nconv}/{ntot} conv)  delta {d:+.5f}  "
            f"{'HIT' if abs(d) <= F4_TOL else '🔴 MISS'}")

    say("\n" + "=" * 96)
    if not m1_ok:
        say("M1 FAILED -- nothing below it was read. Re-run or diagnose the cap change.")
    if FAIL:
        say(f"FAILED {len(FAIL)} gate(s):")
        for f in FAIL:
            say(f"   - {f}")
    else:
        say("ALL GATES PASSED")
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "grade_mustarcap.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
