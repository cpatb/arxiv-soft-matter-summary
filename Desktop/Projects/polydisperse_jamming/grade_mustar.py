"""grade_mustar.py -- PREREG_MUSTAR_2026-09-19.md.

Does polydispersity move the FRICTION SATURATION SCALE mu*, or only the amplitude?

    phi(mu) = phi_inf + A exp(-mu/mu*)

  (M1) SCALE      -- mu* grows with CV: the approach to the frictional limit is genuinely
                     slower for a wider distribution.
  (M2) AMPLITUDE  -- mu* is CV-free and A grows with CV: there is simply more to mobilise.

Gates, quoted from the prereg (Amendment 1 moved the CHANNEL to the matched lowest rung
p_hat = 3e-5, fit-free, BEFORE launch, because phi_J_lin is a fit to noise at mu_p >= 1.0):

  F0  -- every (mu_p, CV) group reports SPAN/SD; below 2.0 is NOT READABLE on the
         extrapolated channel.
  F1  -- WIRING: the re-run (mu_p = 3.0, CV = 0) cells must be BIT-IDENTICAL to the existing
         ones. Registered extrapolated value to reproduce: phi_J_lin = 0.55742.
  F2  -- 72 of 72 converge.
  F3  -- PRIMARY. REGISTERED: mu*(0.5)/mu*(0) in [1.05, 1.45] AND A(0.5)/A(0) in
         [0.90, 1.20]. ⚠ "If A moves as much as mu* does, the two mechanisms are NOT
         separated by this design and the verdict is UNDETERMINED regardless of what mu*
         does." Errors by seed bootstrap; per-seed fits printed so one bad seed is visible.
  F4  -- blind point predictions, matched-rung channel, each within +-0.0030:
             phi(mu=0.6, CV=0.0) = 0.56394
             phi(mu=0.6, CV=0.5) = 0.62042
             phi(mu=3.0, CV=0.5) = 0.61346
  F5  -- the FORM. With 5 points and 3 parameters the fit gains 2 residual dof. REGISTERED:
         if the residual RMS exceeds the mean seed s.e. (~0.0015) the exponential form is
         INADEQUATE and F3's ratio is quoted as form-dependent. An algebraic alternative
         phi_inf + A/(1 + mu/mu*) is named in advance and BOTH are reported.
  F6  -- mu*(CV) is monotone across the three CV levels.
"""

from __future__ import annotations

import collections
import glob
import json
import math
import os
import sys

import numpy as np
from scipy.optimize import curve_fit

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

RUNGS = (1e-3, 3e-4, 1e-4, 3e-5)
MUS = (0.0, 0.3, 0.6, 1.0, 3.0)
CVS = (0.0, 0.2, 0.5)
MUSTAR_BAND = (1.05, 1.45)
A_BAND = (0.90, 1.20)
F4 = {(0.6, 0.0): 0.56394, (0.6, 0.5): 0.62042, (3.0, 0.5): 0.61346}
F4_TOL = 0.0030
FAIL: list[str] = []


def gate(name: str, ok: bool, detail: str = "") -> None:
    print(f"   [{'PASS' if ok else 'FAIL'}] {name}" + (f"  --  {detail}" if detail else ""))
    if not ok:
        FAIL.append(name)


def load_all():
    """Every lognormal N=1500 record, keyed (mu_p, cv, seed) -> rung -> record.
    New `mustar` records win over old ones at the same key, and that is stated: they are
    the campaign under grade."""
    cells: dict = collections.defaultdict(dict)
    n_new = 0
    for fn in sorted(glob.glob(os.path.join(OUT, "phase2*.jsonl"))):
        is_new = "mustar" in os.path.basename(fn)
        if "wcbidi" in os.path.basename(fn):
            continue
        for ln in open(fn, encoding="utf-8"):
            ln = ln.strip()
            if not ln:
                continue
            r = json.loads(ln)
            if r.get("failed") or r.get("n_requested") != 1500:
                continue
            if (r.get("family") or "lognormal") != "lognormal":
                continue
            n_new += is_new
            k = (r["mu_p"], round(r["cv"], 1), r["seed"])
            p = r["p_hat_target"]
            r["_new"] = is_new
            prev = cells[k].get(p)
            if prev is None or is_new or r.get("max_iter", 0) > prev.get("max_iter", 0):
                cells[k][p] = r
    return cells, n_new


def rung(cells, mu, cv):
    return {s: cells[(mu, cv, s)][3e-05]["phi"] for s in (1, 2, 3)
            if (mu, cv, s) in cells and 3e-05 in cells[(mu, cv, s)]}


def expo(mu, pi, A, ms):
    return pi + A * np.exp(-mu / ms)


def alge(mu, pi, A, ms):
    return pi + A / (1.0 + mu / ms)


def fit(model, xs, ys, ses):
    p0 = [ys[-1], ys[0] - ys[-1], 0.25]
    p, _ = curve_fit(model, xs, ys, p0=p0, sigma=ses, absolute_sigma=True, maxfev=40000)
    res = ys - model(xs, *p)
    return p, float(np.sqrt(np.mean(res ** 2)))


def main() -> int:
    print("=" * 92)
    print("grade_mustar.py -- PREREG_MUSTAR_2026-09-19.md")
    print("=" * 92)
    cells, n_new = load_all()
    print(f"\n   new `mustar` records on disk: {n_new}  (registered 72)")
    if n_new == 0:
        print("   🔴 NOTHING DELIVERED. Not graded.")
        return 1

    # ---------------- F1 : wiring -------------------------------------------------------
    print("\n=== F1  WIRING: the re-run (mu_p = 3.0, CV = 0) cells ===")
    old3 = collections.defaultdict(dict)
    for fn in sorted(glob.glob(os.path.join(OUT, "phase2*.jsonl"))):
        if "mustar" in os.path.basename(fn) or "wcbidi" in os.path.basename(fn):
            continue
        for ln in open(fn, encoding="utf-8"):
            ln = ln.strip()
            if not ln:
                continue
            r = json.loads(ln)
            if (r.get("failed") or r.get("mu_p") != 3.0 or abs(r.get("cv", 1)) > 1e-12
                    or r.get("n_requested") != 1500):
                continue
            old3[r["seed"]].setdefault(r["p_hat_target"], []).append(r["phi"])
    newv = rung(cells, 3.0, 0.0)
    hit = tot = 0
    for s, v in sorted(newv.items()):
        olds = old3.get(s, {}).get(3e-05, [])
        if not olds:
            continue
        tot += 1
        m = any(abs(v - x) < 1e-12 for x in olds)
        hit += m
        print(f"      seed {s}: new {v:.9f} vs existing {[round(x, 9) for x in olds]}  "
              f"{'MATCH' if m else '🔴 no match'}")
    gate("F1: the re-run mu_p = 3.0, CV = 0 cells reproduce an existing record exactly",
         tot > 0 and hit == tot, f"{hit}/{tot}")

    # ---------------- F2 : health -------------------------------------------------------
    print("\n=== F2  LADDER HEALTH (registered 72 of 72) ===")
    conv = tot2 = 0
    for k, d in cells.items():
        if k[0] not in (0.6, 3.0) or k[1] not in CVS:
            continue
        for p, r in d.items():
            if not r.get("_new"):
                continue
            tot2 += 1
            conv += bool(r.get("converged"))
    print(f"   converged {conv} of {tot2} new records")
    gate("F2: every new record converged", conv == tot2 and tot2 == 72, f"{conv}/{tot2}")

    # ---------------- F0 : readability --------------------------------------------------
    print("\n=== F0  SPAN/SD per (mu_p, CV) group ===")
    print(f"   {'mu_p':>5s} {'CV':>5s} {'span':>9s} {'sd':>9s} {'SPAN/SD':>8s}")
    for mu in MUS:
        for cv in CVS:
            sp, bot = [], []
            for s in (1, 2, 3):
                d = cells.get((mu, cv, s))
                if d and all(p in d for p in RUNGS):
                    ys = [d[p]["phi"] for p in RUNGS]
                    sp.append(ys[0] - ys[-1])
                    bot.append(ys[-1])
            if len(bot) < 3:
                continue
            sd_ = float(np.std(bot, ddof=1))
            print(f"   {mu:5.1f} {cv:5.1f} {np.mean(sp):9.5f} {sd_:9.5f} "
                  f"{np.mean(sp)/sd_ if sd_ else float('inf'):8.1f}")

    # ---------------- the ladders -------------------------------------------------------
    print("\n=== THE mu LADDERS, matched rung p_hat = 3e-5 (Amendment 1's channel) ===")
    lad = {}
    print(f"   {'CV':>5s} " + "".join(f"{m:>12}" for m in MUS))
    for cv in CVS:
        row = f"   {cv:5.1f} "
        pts = {}
        for mu in MUS:
            v = rung(cells, mu, cv)
            if len(v) == 3:
                a = np.array(list(v.values()))
                pts[mu] = (float(a.mean()), float(a.std(ddof=1) / math.sqrt(3)))
                row += f"{pts[mu][0]:12.5f}"
            else:
                row += f"{'--':>12s}"
        lad[cv] = pts
        print(row)

    # ---------------- F3 / F5 : the fits ------------------------------------------------
    print("\n=== F3 / F5  FITS: exponential and the pre-named algebraic alternative ===")
    fits = {}
    for nm, model in (("exponential", expo), ("algebraic", alge)):
        print(f"   --- {nm} ---")
        print(f"   {'CV':>5s} {'n pts':>6s} {'phi_inf':>10s} {'A':>10s} {'mu*':>9s} "
              f"{'resid RMS':>11s}")
        fits[nm] = {}
        for cv in CVS:
            pts = lad[cv]
            if len(pts) < 4:
                print(f"   {cv:5.1f} {len(pts):6d}  too few points")
                continue
            xs = np.array(sorted(pts))
            ys = np.array([pts[m][0] for m in sorted(pts)])
            ses = np.array([max(pts[m][1], 1e-6) for m in sorted(pts)])
            try:
                p, rms = fit(model, xs, ys, ses)
            except Exception as exc:               # noqa: BLE001
                print(f"   {cv:5.1f} fit failed: {exc!r}")
                continue
            fits[nm][cv] = (p, rms, len(xs))
            print(f"   {cv:5.1f} {len(xs):6d} {p[0]:10.5f} {p[1]:10.5f} {p[2]:9.4f} "
                  f"{rms:11.6f}")
    mean_se = float(np.mean([pts[m][1] for cv in CVS for m in lad[cv] for pts in [lad[cv]]]))
    print(f"\n   mean seed s.e. over all points: {mean_se:.6f}  (F5's bar)")
    for nm in fits:
        bad = [cv for cv, (p, rms, n) in fits[nm].items() if rms > mean_se]
        print(f"   {nm}: residual RMS exceeds the seed s.e. at CV {bad or 'nowhere'}")
    ex = fits.get("exponential", {})
    if 0.0 in ex and 0.5 in ex:
        r_mu = ex[0.5][0][2] / ex[0.0][0][2]
        r_A = ex[0.5][0][1] / ex[0.0][0][1]
        print(f"\n   🔑 mu*(0.5)/mu*(0) = {r_mu:.3f}   registered {MUSTAR_BAND}")
        print(f"      A(0.5)/A(0)     = {r_A:.3f}   registered {A_BAND}")
        a_ok = A_BAND[0] <= r_A <= A_BAND[1]
        mu_ok = MUSTAR_BAND[0] <= r_mu <= MUSTAR_BAND[1]
        if not a_ok:
            print("      ⚠ A moved outside its registered band, so the two mechanisms are "
                  "NOT separated by this design -> UNDETERMINED, as registered.")
        gate("F3a: the amplitude A is stable, so mu* and A ARE separated by this design",
             a_ok, f"A ratio {r_A:.3f} vs {A_BAND}")
        gate("F3b: mu*(0.5)/mu*(0) lands in the registered band, i.e. (M1) SCALE",
             mu_ok and a_ok, f"mu* ratio {r_mu:.3f} vs {MUSTAR_BAND}"
             + ("" if a_ok else "  [moot: F3a failed]"))
    if len(ex) >= 3:
        ms = [ex[cv][0][2] for cv in sorted(ex)]
        mono = all(ms[i] <= ms[i + 1] for i in range(len(ms) - 1)) or \
               all(ms[i] >= ms[i + 1] for i in range(len(ms) - 1))
        print(f"\n   F6 mu*(CV) = {[round(x, 4) for x in ms]} at CV {sorted(ex)}")
        gate("F6: mu*(CV) is monotone across the three CV levels", mono, f"{ms}")

    # ---------------- F4 : blind predictions --------------------------------------------
    print("\n=== F4  BLIND POINT PREDICTIONS (matched rung, +-0.0030) ===")
    for (mu, cv), pred in sorted(F4.items()):
        v = rung(cells, mu, cv)
        if len(v) < 2:
            print(f"   mu={mu} CV={cv}: not delivered")
            continue
        a = np.array(list(v.values()))
        m, se = float(a.mean()), float(a.std(ddof=1) / math.sqrt(len(a)))
        d = m - pred
        print(f"   phi(mu={mu}, CV={cv}): predicted {pred:.5f}  measured {m:.5f} +- {se:.5f}"
              f"  delta {d:+.5f}  {'HIT' if abs(d) <= F4_TOL else '🔴 MISS'}")
        gate(f"F4: phi(mu={mu}, CV={cv}) within +-{F4_TOL} of its blind prediction",
             abs(d) <= F4_TOL, f"{d:+.5f}")

    print("\n" + "=" * 92)
    if FAIL:
        print(f"FAILED {len(FAIL)} gate(s):")
        for f in FAIL:
            print(f"   - {f}")
        return 1
    print("ALL GATES PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
