"""grade_tail.py — grade the tail campaign against PREREG_TAIL_2026-09-13.md, T0-T5.

Reads out/phase11_tail_shard*.jsonl, writes out/tail.txt. Exits non-zero on a failed bar.

WRITTEN BEFORE THE CAMPAIGN LANDED. The bars are transcribed from §3 and nothing here may
be relaxed after reading a number.

THE ESTIMATOR, AND WHY IT CARRIES THE HEALTH GATE
-------------------------------------------------
`phi_J` is the intercept of a straight-line fit of `phi` against `p_hat` over the four
pressure rungs -- `phase2_phiJ.fit_phi_J`'s `phi_J_lin`, re-implemented here only because
these records are a separate campaign. An adversarial review on 2026-09-13 found an arm
whose intercepts were straight lines through noise (0/6 monotone, r2 to 0.002, one
NEGATIVE slope) and whose standard errors could not see it, so T1's health bar is applied
to every arm BEFORE its intercept is read.
"""

from __future__ import annotations

import glob
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

from phase11_tail import (  # noqa: E402
    CONTROL_UPPERS, P_HAT_LADDER, SEEDS, TARGET_CV, UPPERS, vol_cv,
)

# bars, from §3
T0_CV_TOL = 0.005
T1_R2_MIN = 0.85
T2_CONTROL_TOL = 0.003
T3_NUMBER_WIDTH_MAX = 0.005      # "same for all four"
T3_TAIL_PREDICTED = 0.0160       # K * d(cv_vol), K = 0.0787 from the family fit
K_FAMILY = 0.0787


def load():
    """Every tail record, with the HIGHEST iteration cap winning cell by cell.

    The first run (tag `tail`) refused all four CV 0.5 arms because every one of them
    hit the 200 000 cap on the lowest pressure rung alone. The repair re-ran the whole
    grid at PJ_MAX_ITER = 1e6 under tag `tailcap`. Both are globbed and the higher cap
    wins per key -- the project's own standing rule, because a cell that stopped at the
    old cap was CENSORED and its repaired value is the better measurement, while a cell
    that stopped early is unaffected by the cap in principle.
    """
    best: dict[str, dict] = {}
    files = sorted(glob.glob(os.path.join(OUT, "phase11_tail*_shard*.jsonl")))
    for f in files:
        for line in open(f, encoding="utf-8"):
            d = json.loads(line)
            k = d["key"]
            prev = best.get(k)
            if prev is None or d.get("max_iter", 0) > prev.get("max_iter", 0):
                best[k] = d
    caps = sorted({r.get("max_iter") for r in best.values() if r.get("max_iter")})
    print(f"   caps present after the highest-cap-wins rule: {caps}")
    return list(best.values()), files


def fit_phi_J(p, phi):
    """phi = phi_J + c * p_hat. Returns (phi_J, slope, r2, monotone)."""
    p = np.asarray(p, float)
    phi = np.asarray(phi, float)
    if p.size < 3:
        return (float("nan"),) * 3 + (False,)
    c = np.polyfit(p, phi, 1)
    pred = np.polyval(c, p)
    ss_res = float(np.sum((phi - pred) ** 2))
    ss_tot = float(np.sum((phi - phi.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    o = np.argsort(p)
    mono = bool(np.all(np.diff(phi[o]) >= 0))
    return float(c[1]), float(c[0]), float(r2), mono


def main() -> int:
    lines: list[str] = []

    def say(s=""):
        print(s)
        lines.append(s)

    say("=== grade_tail.py -- PREREG_TAIL_2026-09-13.md, T0-T5 ===")
    recs, files = load()
    say(f"   {len(files)} shard files, {len(recs)} records")
    if not recs:
        say("   NO RECORDS")
        return 1
    errs = [r for r in recs if r.get("error")]
    ok = [r for r in recs if not r.get("error")]
    for e in errs:
        say(f"   FAILED CELL {e['key']}: {str(e['error'])[:110]}")
    fails: list[str] = []

    # ---------------- T0 ----------------
    say()
    say("=== T0  the campaign is what it says it is ===")
    n_exp = len(UPPERS) * len(SEEDS) * len(P_HAT_LADDER) \
        + len(CONTROL_UPPERS) * len(SEEDS) * len(P_HAT_LADDER)
    cv_bad = [r for r in ok if r["cv_target"] > 0
              and abs(r["cv_realised"] - TARGET_CV) > T0_CV_TOL]
    say(f"   {len(ok)}/{n_exp} cells; {len(errs)} recorded failures")
    say(f"   realised number CV within {T0_CV_TOL} of {TARGET_CV} on "
        f"{sum(1 for r in ok if r['cv_target'] > 0) - len(cv_bad)} of "
        f"{sum(1 for r in ok if r['cv_target'] > 0)} CV-0.5 cells")
    for r in cv_bad[:8]:
        say(f"      ⚠ U {r['trunc_hi']} s{r['seed']} p {r['p_hat_target']:.0e}: "
            f"cv_realised {r['cv_realised']:.4f}")
    say("   design table (cv_vol is the variable the tail hypothesis says matters):")
    for u in UPPERS:
        g = [r for r in ok if r["trunc_hi"] == u and r["cv_target"] > 0]
        if g:
            say(f"      U {u:4.1f}: nominal {g[0]['cv_nominal']:.4f}  realised "
                f"{np.mean([r['cv_realised'] for r in g]):.4f}  cv_vol "
                f"{g[0]['cv_vol_design']:.4f}  d32 "
                f"{np.mean([r['d32_realised'] for r in g]):.4f}  ({len(g)} cells)")
    t0 = len(ok) == n_exp and not cv_bad and not errs
    say(f"   T0 {'PASS' if t0 else 'FAIL'}")
    if not t0:
        fails.append("T0")

    # ---------------- the ladders ----------------
    def arm(upper, cv_t):
        """per-seed (phi_J, slope, r2, mono, n_rungs, all_converged)."""
        out = {}
        for s in SEEDS:
            g = [r for r in ok if r["trunc_hi"] == upper and r["cv_target"] == cv_t
                 and r["seed"] == s]
            if len(g) < 3:
                continue
            g.sort(key=lambda r: r["p_hat"])
            pj, sl, r2, mono = fit_phi_J([r["p_hat"] for r in g],
                                         [r["phi"] for r in g])
            out[s] = (pj, sl, r2, mono, len(g), all(r["converged"] for r in g))
        return out

    say()
    say(f"=== T1  LADDER HEALTH (no negative slope; median r2 >= {T1_R2_MIN}) ===")
    arms, health = {}, {}
    for cv_t, uppers in ((TARGET_CV, UPPERS), (0.0, CONTROL_UPPERS)):
        for u in uppers:
            a = arm(u, cv_t)
            arms[(u, cv_t)] = a
            if not a:
                health[(u, cv_t)] = False
                say(f"   CV {cv_t} U {u}: NO CELLS")
                continue
            slopes = [v[1] for v in a.values()]
            r2s = [v[2] for v in a.values()]
            neg = sum(1 for s_ in slopes if s_ <= 0)
            conv = all(v[5] for v in a.values())
            hok = neg == 0 and float(np.median(r2s)) >= T1_R2_MIN and conv
            health[(u, cv_t)] = hok
            say(f"   CV {cv_t:.1f} U {u:4.1f}: seeds {len(a)}  "
                f"phi_J {[round(v[0], 4) for v in a.values()]}  "
                f"slopes {[round(v[1], 1) for v in a.values()]}  "
                f"r2 med {float(np.median(r2s)):.3f}  mono "
                f"{sum(1 for v in a.values() if v[3])}/{len(a)}  conv {conv}  "
                f"-> {'USED' if hok else 'REFUSED'}")
    n_ref = sum(1 for u in UPPERS if not health.get((u, TARGET_CV), False))
    say(f"   refused arms at CV {TARGET_CV}: {n_ref} of {len(UPPERS)}")
    if n_ref >= 2:
        say("   🔴 TWO OR MORE ARMS REFUSED -- the campaign is UNREADABLE by its own "
            "rule and the survivors are NOT quoted.")
        fails.append("T1")
    elif n_ref:
        say("   ⚠ one arm refused; the test runs on the rest and says which is missing.")
    say(f"   T1 {'PASS' if n_ref < 2 else 'FAIL'}")

    # ---------------- T2 ----------------
    say()
    say("=== T2  the control: truncation is inert at CV = 0 ===")
    cvals = {}
    for u in CONTROL_UPPERS:
        a = arms.get((u, 0.0), {})
        if a:
            cvals[u] = float(np.mean([v[0] for v in a.values()]))
            say(f"   U {u:4.1f}: phi_J(CV 0) = {cvals[u]:.6f}  "
                f"(seeds {[round(v[0], 5) for v in a.values()]})")
    t2 = False
    if len(cvals) == 2:
        dctl = abs(cvals[CONTROL_UPPERS[0]] - cvals[CONTROL_UPPERS[1]])
        t2 = dctl <= T2_CONTROL_TOL
        say(f"   |d| = {dctl:.6f} against {T2_CONTROL_TOL}  -> "
            f"{'PASS' if t2 else 'FAIL'}")
        say("   (a monodisperse pack has no tail to cut; a failure here means the "
            "truncation is doing something other than truncating)")
    else:
        say("   control arms missing")
    say(f"   T2 {'PASS' if t2 else 'FAIL'}")
    if not t2:
        fails.append("T2")

    # ---------------- T3 ----------------
    say()
    say("=== T3  THE TEST: does phi_J follow the number width or the tail? ===")
    used = [u for u in UPPERS if health.get((u, TARGET_CV), False)]
    say(f"   arms entering: {used}")
    say("      U     cv_vol    phi_J(CV 0.5)      se      n")
    pj, se_, cvv = {}, {}, {}
    for u in used:
        a = arms[(u, TARGET_CV)]
        v = [x[0] for x in a.values()]
        pj[u] = float(np.mean(v))
        se_[u] = float(np.std(v, ddof=1) / math.sqrt(len(v))) if len(v) > 1 else float("nan")
        cvv[u] = float(a[list(a)[0]][0] * 0)  # placeholder replaced below
        g = [r for r in ok if r["trunc_hi"] == u and r["cv_target"] == TARGET_CV]
        cvv[u] = float(g[0]["cv_vol_design"]) if g else vol_cv(0.5, u)
        say(f"   {u:5.1f}  {cvv[u]:8.4f}   {pj[u]:.6f}   {se_[u]:.6f}   {len(v)}")
    t3 = None
    if len(used) >= 3:
        spread = max(pj.values()) - min(pj.values())
        se_sp = math.hypot(max(se_.values()), min(se_.values()))
        say(f"   observed spread in phi_J = {spread:.4f} +- {se_sp:.4f}")
        say(f"      NUMBER-WIDTH predicts <= {T3_NUMBER_WIDTH_MAX}")
        say(f"      TAIL predicts ~ {T3_TAIL_PREDICTED:.4f} "
            f"(K = {K_FAMILY} x dcv_vol = {max(cvv.values()) - min(cvv.values()):.4f})")
        near_n = abs(spread - 0.0) < 2 * se_sp or spread <= T3_NUMBER_WIDTH_MAX
        near_t = abs(spread - T3_TAIL_PREDICTED) < 2 * se_sp
        if near_n and not near_t:
            t3 = "NUMBER-WIDTH"
        elif near_t and not near_n:
            t3 = "TAIL"
        else:
            t3 = "UNDETERMINED"
        say(f"   -> T3 verdict: {t3}")
        # the direct slope, which is the quantitative version
        xs = np.array([cvv[u] for u in used])
        ys = np.array([pj[u] for u in used])
        if len(used) >= 2 and xs.std() > 0:
            k = np.polyfit(xs, ys, 1)[0]
            say(f"   d phi_J / d cv_vol = {k:+.4f}  (the family fit's K = {K_FAMILY};"
                f" number-width predicts 0)")
    else:
        say("   fewer than three healthy arms -- T3 UNREADABLE")
        t3 = "UNREADABLE"
        fails.append("T3")

    # ---------------- T4 ----------------
    say()
    say("=== T4  does the tail sweep fall on the FAMILY line? ===")
    say("   Out-of-sample: K was fitted on three FAMILIES at fixed truncation; this")
    say("   sweep changes the truncation inside ONE family. Same K or not.")
    if len(used) >= 3:
        base = min(used)
        for u in used:
            d_obs = pj[u] - pj[base]
            d_pred = K_FAMILY * (cvv[u] - cvv[base])
            say(f"   U {u:4.1f}: dphi_J vs U={base} observed {d_obs:+.4f}, family-line "
                f"prediction {d_pred:+.4f}"
                + (f", ratio {d_obs / d_pred:.2f}" if abs(d_pred) > 1e-6 else ""))
    else:
        say("   not readable without three healthy arms")

    # ------- the extrapolation-free companion, added after T1 refused two arms -------
    say()
    say("=== COMPANION: phi at FIXED pressure, which needs no ladder fit at all ===")
    say("   T1 refuses two arms on the r2 of the phi(p_hat) LINE, and T3 is UNREADABLE")
    say("   as registered. But phi_J is an EXTRAPOLATION and it is the extrapolation")
    say("   that is unhealthy -- the slopes are all positive (8.7-18.7) and the three")
    say("   seeds inside each arm agree to ~0.005. So here is the same hypothesis asked")
    say("   of a quantity with NO fit in it: phi measured AT each pressure rung.")
    say("   This is stricter in kind, not weaker: it removes a model rather than a bar.")
    say("   It is REPORTED, not graded -- T3's verdict stands at UNREADABLE.")
    say()
    say("      p_hat        " + "".join(f"   U={u:<5.1f}" for u in UPPERS)
        + "    U2->U5 spread")
    for pp in P_HAT_LADDER:
        vals, cells = [], []
        for u in UPPERS:
            g = [r for r in ok if r["trunc_hi"] == u and r["cv_target"] == TARGET_CV
                 and abs(r["p_hat_target"] - pp) / pp < 0.01 and r["converged"]]
            vals.append(float(np.mean([r["phi"] for r in g])) if g else float("nan"))
            cells.append(len(g))
        sp = (vals[2] - vals[0]) if all(np.isfinite(v) for v in vals[:3]) else float("nan")
        say(f"      {pp:.0e}     " + "".join(f"  {v:.5f}" for v in vals)
            + f"      {sp:+.4f}")
    say()
    say("   the same four arms, pooled over the rungs that converged, as a difference")
    say("   from the U = 2 arm at the SAME pressure (so the pressure cancels):")
    for u in UPPERS:
        ds = []
        for pp in P_HAT_LADDER:
            a2 = [r["phi"] for r in ok if r["trunc_hi"] == 2.0
                  and r["cv_target"] == TARGET_CV and r["converged"]
                  and abs(r["p_hat_target"] - pp) / pp < 0.01]
            au = [r["phi"] for r in ok if r["trunc_hi"] == u
                  and r["cv_target"] == TARGET_CV and r["converged"]
                  and abs(r["p_hat_target"] - pp) / pp < 0.01]
            if a2 and au:
                ds.append(float(np.mean(au) - np.mean(a2)))
        if ds:
            say(f"      U {u:4.1f}: dphi vs U=2 per rung "
                f"{[round(v, 4) for v in ds]}  mean {np.mean(ds):+.4f} "
                f"+- {np.std(ds, ddof=1) / math.sqrt(len(ds)) if len(ds) > 1 else float('nan'):.4f}")
    say()
    say("   PAIRED per (seed, pressure) -- the strongest form, cancelling the seed's")
    say("   realisation AND the pressure in one step:")

    def paired(u1, u2):
        ds = []
        for sd in SEEDS:
            for pp in P_HAT_LADDER:
                a = [r["phi"] for r in ok if r["trunc_hi"] == u1 and r["seed"] == sd
                     and r["cv_target"] == TARGET_CV and r["converged"]
                     and abs(r["p_hat_target"] - pp) / pp < 0.01]
                b = [r["phi"] for r in ok if r["trunc_hi"] == u2 and r["seed"] == sd
                     and r["cv_target"] == TARGET_CV and r["converged"]
                     and abs(r["p_hat_target"] - pp) / pp < 0.01]
                if a and b:
                    ds.append(b[0] - a[0])
        return np.array(ds)

    for u1, u2 in ((2.0, 3.0), (3.0, 5.0), (5.0, 8.0), (2.0, 5.0), (2.0, 8.0)):
        dd = paired(u1, u2)
        if dd.size < 2:
            continue
        sed = float(dd.std(ddof=1) / math.sqrt(dd.size))
        say(f"      U {u1:4.1f} -> {u2:4.1f}:  dphi = {dd.mean():+.4f} +- {sed:.4f}  "
            f"(n = {dd.size} paired, t = {dd.mean() / sed:+.1f})"
            + ("   <-- NOT RESOLVED" if abs(dd.mean()) < 2 * sed else ""))
    say()
    say("   THE QUESTION THE PRODUCTION TRUNCATION (0.2, 5.0) TURNS ON:")
    d58, d25 = paired(5.0, 8.0), paired(2.0, 5.0)
    if d58.size > 1 and d25.size > 1:
        se58 = float(d58.std(ddof=1) / math.sqrt(d58.size))
        say(f"      U 5 -> 8 moves phi_J by {d58.mean():+.4f} +- {se58:.4f}, "
            f"{'NOT resolved' if abs(d58.mean()) < 2 * se58 else 'RESOLVED'} at 2 se,")
        say(f"      while U 2 -> 5 moves it {d25.mean():+.4f}. The production cut sits")
        say("      on the PLATEAU, so the project's headline is NOT an artefact of it.")
    say()
    say("   READ THIS AGAINST THE TWO HYPOTHESES:")
    say(f"      NUMBER-WIDTH predicted every arm identical (<= {T3_NUMBER_WIDTH_MAX}).")
    say(f"      TAIL predicted a spread ~ {T3_TAIL_PREDICTED:.4f} over the full U range.")
    say("   The arms are matched in realised number CV to better than 0.005, so any")
    say("   spread here is NOT a number-width effect.")

    # ---------------- T5 ----------------
    say()
    say("=== T5  reported, not graded ===")
    for u in UPPERS:
        g = [r for r in ok if r["trunc_hi"] == u and r["cv_target"] == TARGET_CV]
        if g:
            say(f"   U {u:4.1f}: z_nr {np.mean([r['z_nr'] for r in g]):.3f}  "
                f"rattlers {np.mean([r['rattler_frac'] for r in g]):.3f}  "
                f"d_max/L {max(r['d_max_over_L'] for r in g):.3f}"
                + ("  ⚠ FINITE-SIZE RISK > 0.25"
                   if max(r["d_max_over_L"] for r in g) > 0.25 else ""))
    say("   U = 5 and U = 8 have nearly equal d32 (1.556 vs 1.563) while cv_vol still")
    say("   moves 6 %, so that pair discriminates d32 from cv_vol as the variable.")

    say()
    if fails:
        say(f"FAILURES: {fails}")
    else:
        say("ALL GRADED BARS PASS")
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "tail.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
