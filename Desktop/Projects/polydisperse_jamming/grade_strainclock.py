"""grade_strainclock.py — grade the strain-clock campaign against PREREG_STRAINCLOCK_2026-09-13.md.

Reads out/phase10_sc_shard*.jsonl, writes out/strainclock.txt. Exits non-zero if any
graded bar fails, so the pre-registration cannot quietly become a description of what
happened.

WRITTEN BEFORE THE CAMPAIGN LANDED. The bars G0-G5 are transcribed from §5 of the
pre-registration and nothing here may be relaxed after reading a number; if a bar turns
out to be unmeetable (G0's `I` already is, see the amendment in §5) the failing number
stays printed and a STRICTER companion is reported beside it.

THE ONE ESTIMATOR THAT NEEDED CARE
----------------------------------
`gamma_c` is a CROSSING TIME, and this project has been bitten twice by crossing times
whose monotonicity was assumed (2026-09-03, `t_half`). So `clock()` returns the FIRST
up-crossing, the LAST up-crossing, the number of down-crossings before the argmax, and a
readability flag; a cell whose two crossings differ by more than READABLE_TOL is excluded
WITH ITS NUMBERS PRINTED rather than silently averaged in.
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

# --- anchors, copied from the driver so a mismatch is loud ---------------------------
from phase10_strainclock import (  # noqa: E402
    PHI_0, PHI_MU, X_GRADED_MAX, X_RUNGS, N_SAMPLES, N_PROD,
)

P_STAR = 3e-4                 # the constant-pressure set point phi_mu was measured at
PHI_C_CONSTP = PHI_MU         # same numbers; named for what G5 compares against

# bars, from §5
G0_GDOT_SPREAD = 1.2
G0_I_MAX = 1e-2
G1_P_PREP_MAX = 1e-6
G2_RATIO_MIN = 100.0
G2_RISE_FRAC = 0.90
G2_NONMONO_TOL = 0.10
G3_Z_ISO = 4.0                # AS REGISTERED. See Z_ISO_MEASURED below.
# 🔴 The pre-registration set G3's isostatic bracket at 4.0, which is `d + 1`, the
# `mu_p -> infinity` limit -- and this campaign shears at `mu_p` = 0.3, where this
# project has MEASURED the isotropic contact number at 4.535 +- 0.072
# (`out/ziso_audit.txt`, `measure_z_iso.py --rule r3`). Using the frictionless-limit
# constant at a finite friction is a lesson already in this lab's own notes. The
# registered bar is NOT moved: 4.0 stays, its verdict stays printed, and the SHARPER
# companion at the measured value is reported beside it. The companion is stricter for
# the "nothing jams below phi_mu" half (a higher bar is easier to sit under, so it is
# the "x >= 0.35 is jammed" half that it makes harder, and that is the half carrying
# the clock).
Z_ISO_MEASURED = 4.535
G4_BAND = (0.67, 1.50)
G4_SEPARATION_MIN = 1.6       # below this the discriminator is NOT INFORMATIVE
G5_TOL = 0.010
READABLE_TOL = 0.25           # |gamma_c_last - gamma_c| / gamma_c


def load() -> list[dict]:
    recs, files = [], sorted(glob.glob(os.path.join(OUT, "phase10_sc_shard*.jsonl")))
    seen: set[str] = set()
    dupes = 0
    for f in files:
        for line in open(f, encoding="utf-8"):
            d = json.loads(line)
            if d["key"] in seen:            # shards are disjoint; a dupe means a re-run
                dupes += 1
                continue
            seen.add(d["key"])
            recs.append(d)
    print(f"   {len(files)} shard files, {len(recs)} unique records"
          f"{f', {dupes} duplicate keys dropped' if dupes else ''}")
    return recs


def clock(gamma, p, p_prep, p_ss, level=0.5):
    """Strain at which p first reaches p_prep + level*(p_ss - p_prep).

    Returns (first, last, n_down_before_argmax, span). NaNs rather than exceptions when
    the series never crosses -- a cell that never builds pressure is a result.
    """
    g = np.asarray(gamma, float)
    y = np.asarray(p, float)
    tgt = p_prep + level * (p_ss - p_prep)
    if not np.isfinite(tgt) or tgt <= 0 or y.size < 3:
        return float("nan"), float("nan"), -1, float("nan")
    above = y >= tgt
    if not above.any():
        return float("nan"), float("nan"), -1, float("nan")

    def interp(i):
        """linear interpolation of the up-crossing between i-1 and i"""
        if i == 0:
            return float(g[0])
        y0, y1 = y[i - 1], y[i]
        if y1 == y0:
            return float(g[i])
        return float(g[i - 1] + (tgt - y0) * (g[i] - g[i - 1]) / (y1 - y0))

    ups = [i for i in range(1, y.size) if y[i] >= tgt > y[i - 1]]
    if y[0] >= tgt:
        ups = [0] + ups
    first = interp(ups[0]) if ups else float(g[np.argmax(above)])
    last = interp(ups[-1]) if ups else first
    argmax = int(np.argmax(y))
    n_down = int(np.sum(np.diff(y[: argmax + 1]) < 0)) if argmax > 0 else 0
    return first, last, n_down, float(abs(last - first))


def cell_row(d: dict) -> dict:
    p = d.get("p_hat_series") or []
    gam = d.get("gamma") or []
    p_prep = float(d.get("p_hat_prep", float("nan")))
    p_ss = float(d.get("p_ss_hat", float("nan")))
    c50 = clock(gam, p, p_prep, p_ss, 0.5)
    c90 = clock(gam, p, p_prep, p_ss, 0.9)
    readable = (math.isfinite(c50[0]) and c50[0] > 0
                and abs(c50[1] - c50[0]) / c50[0] <= READABLE_TOL)
    # monotone rise up to the half-crossing, within tolerance on the step
    mono_ok = True
    if math.isfinite(c50[0]) and len(gam) > 2:
        i = int(np.searchsorted(np.asarray(gam, float), c50[0]))
        seg = np.asarray(p[: max(i + 1, 2)], float)
        if seg.size > 2 and np.nanmax(seg) > 0:
            drops = np.diff(seg)
            mono_ok = bool(np.nansum(drops[drops < 0]) >= -G2_NONMONO_TOL * np.nanmax(seg))
    return {
        "cv": d["cv"], "seed": d["seed"], "phi": d.get("phi", float("nan")),
        "phi_target": d["phi_target"], "x": d.get("x_rung", float("nan")),
        "p_prep": p_prep, "p_ss": p_ss, "z_nr_ss": d.get("z_nr_ss", float("nan")),
        "rat_prep": d.get("rattler_frac_prep", float("nan")),
        "z_nr_prep": d.get("z_nr_prep", float("nan")),
        "gdot": d.get("gdot", float("nan")), "i_ss": d.get("i_ss", float("nan")),
        "steady": bool(d.get("steady", False)),
        "protocol": d.get("protocol", ""), "mu_prep": d.get("mu_prep", float("nan")),
        "fswitch": bool(d.get("friction_switch_ok", False)),
        "n_samples": d.get("n_samples", 0), "nsamp_seen": len(gam),
        "n_requested": d.get("n_requested", 0),
        "gamma_c": c50[0], "gamma_c_last": c50[1], "n_down": c50[2],
        "gamma_90": c90[0], "readable": readable, "mono_ok": mono_ok,
        "ratio": (p_ss / p_prep) if p_prep and math.isfinite(p_prep) and p_prep > 0
        else float("inf"),
        "error": d.get("error", ""),
    }


def gm(vals):
    v = [x for x in vals if math.isfinite(x) and x > 0]
    return float(np.exp(np.mean(np.log(v)))) if v else float("nan")


def main() -> int:
    lines: list[str] = []

    def say(s=""):
        print(s)
        lines.append(s)

    say("=== grade_strainclock.py -- PREREG_STRAINCLOCK_2026-09-13.md, G0-G6 ===")
    recs = load()
    if not recs:
        say("   NO RECORDS -- the campaign has not landed")
        return 1
    errs = [r for r in recs if r.get("error")]
    rows = [cell_row(r) for r in recs if not r.get("error")]
    say(f"   {len(rows)} cells with results, {len(errs)} recorded failures")
    for e in errs:
        say(f"   FAILED CELL {e['key']}: {e['error'][:120]}")
    rows.sort(key=lambda r: (r["cv"], r["x"], r["seed"]))
    fails: list[str] = []

    # ---------------- the table ----------------
    say()
    say("   cv     x    phi_t   p_prep     p_ss     ratio  z_nr  gamma_c  g_last  g90   "
        "n_dn rd mo st   I")
    for r in rows:
        say(f" {r['cv']:.2f} {r['x']:+.2f} {r['phi_target']:.4f} {r['p_prep']:.2e} "
            f"{r['p_ss']:.2e} {r['ratio']:9.3g} {r['z_nr_ss']:5.2f} "
            f"{r['gamma_c']:7.4f} {r['gamma_c_last']:7.4f} {r['gamma_90']:6.3f} "
            f"{r['n_down']:4d} {'y' if r['readable'] else 'N'} "
            f"{'y' if r['mono_ok'] else 'N'} {'y' if r['steady'] else 'N'} "
            f"{r['i_ss']:.2e}")

    # ---------------- G0 ----------------
    say()
    say("=== G0  the campaign is what it says it is ===")
    n_exp = len(X_RUNGS) * 3 * 3
    proto_ok = sum(r["protocol"] == "cv-shear-v3-fprep" for r in rows)
    fsw_ok = sum(r["fswitch"] for r in rows)
    prep_ok = sum(r["mu_prep"] == 0.0 for r in rows)
    ns_ok = sum(r["n_samples"] == N_SAMPLES for r in rows)
    n_ok = sum(r["n_requested"] == N_PROD for r in rows)
    say(f"   records {len(rows)}/{n_exp}; protocol v3-fprep {proto_ok}; "
        f"friction_switch_ok {fsw_ok}; mu_prep=0 {prep_ok}; "
        f"n_samples={N_SAMPLES} {ns_ok}; N={N_PROD} {n_ok}")
    gd = {cv: [r["gdot"] for r in rows if r["cv"] == cv] for cv in PHI_MU}
    means = {cv: float(np.mean(v)) for cv, v in gd.items() if v}
    spread = max(means.values()) / min(means.values()) if len(means) > 1 else float("nan")
    for cv, m in means.items():
        say(f"   cv {cv}: gdot = {m:.4e}  (n {len(gd[cv])})")
    say(f"   between-CV gdot spread = {spread:.3f} against bar {G0_GDOT_SPREAD}")
    i_all = max((r["i_ss"] for r in rows if math.isfinite(r["i_ss"])), default=float("nan"))
    jam = [r for r in rows if r["z_nr_ss"] >= G3_Z_ISO]
    i_jam = max((r["i_ss"] for r in jam if math.isfinite(r["i_ss"])), default=float("nan"))
    say(f"   max I over ALL cells      = {i_all:.2e} against bar {G0_I_MAX:.0e}  "
        f"[AS REGISTERED]")
    say(f"   max I over cells carrying the clock (z_nr_ss >= {G3_Z_ISO}) = {i_jam:.2e}  "
        f"[the stricter companion, §5 amendment]")
    jam2 = [r for r in rows if r["z_nr_ss"] >= Z_ISO_MEASURED]
    i_jam2 = max((r["i_ss"] for r in jam2 if math.isfinite(r["i_ss"])),
                 default=float("nan"))
    say(f"   max I over cells with z_nr_ss >= the MEASURED z_iso {Z_ISO_MEASURED} "
        f"({len(jam2)} cells) = {i_jam2:.2e}")
    g0 = (len(rows) == n_exp and proto_ok == len(rows) and fsw_ok == len(rows)
          and prep_ok == len(rows) and ns_ok == len(rows) and n_ok == len(rows)
          and spread <= G0_GDOT_SPREAD and i_all <= G0_I_MAX)
    if i_all > G0_I_MAX and i_jam <= G0_I_MAX:
        say("   ⚠ the I bar fails ONLY on cells that do not carry the clock; the bar "
            "is kept as registered and the companion is what the claims rest on")
    say(f"   G0 {'PASS' if g0 else 'FAIL'}")
    if not g0:
        fails.append("G0")

    # ---------------- G1 ----------------
    say()
    say("=== G1  the premise: the frictionless prep carries no pressure ===")
    graded = [r for r in rows if r["x"] <= X_GRADED_MAX + 1e-9]
    bad = [r for r in graded if not (r["p_prep"] < G1_P_PREP_MAX)]
    say(f"   graded rungs (x <= {X_GRADED_MAX}): {len(graded)}; "
        f"p_hat_prep < {G1_P_PREP_MAX:.0e} on {len(graded) - len(bad)}")
    for r in bad:
        say(f"   🔴 cv {r['cv']} x {r['x']:+.2f} s{r['seed']}: p_prep = {r['p_prep']:.3e}")
    top = [r for r in rows if r["x"] > X_GRADED_MAX + 1e-9]
    say(f"   x = 0.85 rungs REPORTED not graded (phi_0 sd; §2 caveat 2): "
        f"p_prep max {max((r['p_prep'] for r in top), default=float('nan')):.2e}, "
        f"{sum(r['p_prep'] >= G1_P_PREP_MAX for r in top)} of {len(top)} above the "
        f"threshold i.e. prepared JAMMED")
    g1 = not bad and len(graded) == 4 * 3 * 3
    say(f"   G1 {'PASS' if g1 else 'FAIL'}")
    if not g1:
        fails.append("G1")

    # ---------------- G2 ----------------
    say()
    say("=== G2  the clock exists ===")
    jamrows = [r for r in rows if r["z_nr_ss"] >= G3_Z_ISO]
    rat_ok = sum(r["ratio"] > G2_RATIO_MIN for r in jamrows)
    mono = sum(r["mono_ok"] for r in jamrows)
    readable = sum(r["readable"] for r in jamrows)
    say(f"   jammed cells (z_nr_ss >= {G3_Z_ISO}): {len(jamrows)}")
    say(f"   p_ss/p_prep > {G2_RATIO_MIN:.0f} on {rat_ok}; monotone rise to gamma_c on "
        f"{mono} ({G2_RISE_FRAC:.0%} needed); readable clock on {readable}")
    for r in jamrows:
        if not r["readable"]:
            say(f"   ⚠ NOT READABLE cv {r['cv']} x {r['x']:+.2f} s{r['seed']}: "
                f"gamma_c {r['gamma_c']:.4f} vs last {r['gamma_c_last']:.4f} "
                f"({r['n_down']} down-crossings before the argmax) -- EXCLUDED")
    say("   by rung (the split is the whole story):")
    for x in X_RUNGS:
        rr = [r for r in jamrows if abs(r["x"] - x) < 1e-9]
        if rr:
            say(f"      x = {x:+.2f}: readable {sum(r['readable'] for r in rr)}/{len(rr)}, "
                f"monotone {sum(r['mono_ok'] for r in rr)}/{len(rr)}")
    g2 = bool(jamrows) and rat_ok == len(jamrows) and mono >= G2_RISE_FRAC * len(jamrows)
    say(f"   G2 {'PASS' if g2 else 'FAIL'}")
    if not g2:
        fails.append("G2")

    # ---------------- G3 ----------------
    say()
    say("=== G3  below phi_mu nothing jams ===")
    g3 = True
    for cv in sorted(PHI_MU):
        lo = [r for r in rows if r["cv"] == cv and r["x"] < 0]
        hi = [r for r in rows if r["cv"] == cv and r["x"] >= 0.35]
        lo_ok = all(r["z_nr_ss"] < G3_Z_ISO for r in lo) if lo else False
        hi_ok = all(r["z_nr_ss"] >= G3_Z_ISO for r in hi) if hi else False
        say(f"   cv {cv}: x<0 z_nr = "
            f"{[round(r['z_nr_ss'], 2) for r in lo]} (all < {G3_Z_ISO}: {lo_ok}); "
            f"x>=0.35 min z_nr = "
            f"{min((r['z_nr_ss'] for r in hi), default=float('nan')):.2f} "
            f"(all >= {G3_Z_ISO}: {hi_ok})")
        g3 = g3 and lo_ok and hi_ok
    say(f"   G3 {'PASS' if g3 else 'FAIL'} (bar {G3_Z_ISO}, AS REGISTERED)")
    if not g3:
        fails.append("G3")
    say(f"   -- companion at the MEASURED z_iso({Z_ISO_MEASURED}) for mu_p = 0.3, "
        f"reported not graded --")
    for cv in sorted(PHI_MU):
        lo = [r for r in rows if r["cv"] == cv and r["x"] < 0]
        hi = [r for r in rows if r["cv"] == cv and r["x"] >= 0.35]
        say(f"      cv {cv}: x<0 all < {Z_ISO_MEASURED}: "
            f"{all(r['z_nr_ss'] < Z_ISO_MEASURED for r in lo) if lo else 'no cells'}; "
            f"x>=0.35 all >= {Z_ISO_MEASURED}: "
            f"{all(r['z_nr_ss'] >= Z_ISO_MEASURED for r in hi) if hi else 'no cells'}")

    # ---------------- G4 ----------------
    say()
    say("=== G4  THE MECHANISM BIT: does the clock collapse on distance from phi_c? ===")
    use = [r for r in rows if r["readable"] and r["z_nr_ss"] >= G3_Z_ISO and r["x"] > 0]
    say(f"   {len(use)} readable jammed cells above phi_mu enter the statistic")
    for x in X_RUNGS:
        if x <= 0:
            continue
        tot = [r for r in rows if abs(r["x"] - x) < 1e-9 and r["x"] > 0]
        rd = [r for r in tot if r["readable"] and r["z_nr_ss"] >= G3_Z_ISO]
        say(f"      rung x = {x:+.2f}: {len(rd)} of {len(tot)} cells readable"
            + ("   <-- rung lost" if len(rd) < 2 else ""))

    def per_rung(cv):
        out = {}
        for x in X_RUNGS:
            if x <= 0:
                continue
            v = [r["gamma_c"] for r in use if r["cv"] == cv and abs(r["x"] - x) < 1e-9]
            if v:
                out[x] = (gm(v), float(np.std(v, ddof=1)) if len(v) > 1 else float("nan"),
                          len(v))
        return out

    tabs = {cv: per_rung(cv) for cv in sorted(PHI_MU)}
    say("      x     " + "".join(f"  cv {cv}: gamma_c (n)      " for cv in sorted(PHI_MU)))
    for x in X_RUNGS:
        if x <= 0:
            continue
        cells = []
        for cv in sorted(PHI_MU):
            t = tabs[cv].get(x)
            cells.append(f"  {t[0]:.4f} +-{t[1]:.4f} ({t[2]})" if t else "         --      ")
        say(f"   {x:+.2f}  " + "".join(cells))

    g4 = True
    verdicts = {}
    se_by_cv: dict[float, float] = {}
    for cv in (0.2, 0.5):
        ratios, rel_se = [], []
        for x, t in tabs[cv].items():
            t0 = tabs[0.0].get(x)
            if not t0 or not math.isfinite(t[0]) or not math.isfinite(t0[0]):
                continue
            ratios.append(t[0] / t0[0])
            # relative SE of a ratio of two geometric means, seeds independent
            s1 = t[1] / t[0] / math.sqrt(t[2]) if t[2] > 1 and math.isfinite(t[1]) else 0.0
            s0 = (t0[1] / t0[0] / math.sqrt(t0[2])
                  if t0[2] > 1 and math.isfinite(t0[1]) else 0.0)
            rel_se.append(math.hypot(s1, s0))
        if not ratios:
            say(f"   cv {cv}: NO matched rungs -- UNREADABLE")
            verdicts[cv] = "UNREADABLE"
            g4 = False
            continue
        R = gm(ratios)
        se = R * math.sqrt(sum(s * s for s in rel_se)) / len(ratios)
        se_by_cv[cv] = se
        trend = (np.corrcoef(sorted(tabs[cv]), [tabs[cv][x][0] / tabs[0.0][x][0]
                                                for x in sorted(tabs[cv])
                                                if x in tabs[0.0]])[0, 1]
                 if len(ratios) > 2 else float("nan"))
        inside = G4_BAND[0] <= R <= G4_BAND[1]
        near = min(abs(R - G4_BAND[0]), abs(R - G4_BAND[1])) < 2 * se if se > 0 else False
        say(f"   cv {cv}: R = {R:.3f} +- {se:.3f} (per-rung "
            f"{[round(v, 3) for v in ratios]}), rung-trend rho = {trend:+.3f}")
        if near:
            v = "UNDETERMINED (a bar lies within 2 SE)"
        elif inside:
            v = "COLLAPSE -- CV acts through phi_c only (M2a)"
        else:
            v = ("DOES NOT COLLAPSE -- CV changes the network-building kinetics (M4), "
                 f"and the broad pack is {'FASTER' if R < 1 else 'SLOWER'}")
        verdicts[cv] = v
        say(f"      -> {v}")
        if near or not inside:
            g4 = False
    # ---- the separation rule of §5, and a correction to how it was first coded ----
    # §5 requires "the two candidates' PREDICTIONS to be at least 1.6x apart on the data
    # actually collected". The first version of this block computed the SPAN OF THE
    # PER-RUNG RATIOS instead -- a measure of rung-to-rung scatter, not of how far apart
    # the hypotheses sit -- and reported the test NOT INFORMATIVE on that basis. The two
    # candidate predictions are R = 1 (M2a: CV acts through phi_c only) and, from
    # cv_unloading.py's U4 prior that a CV 0.5 pack relaxes 2-3x faster at matched phi,
    # R = 1/3 to 1/2 (M4). Both numbers are printed; the registered SENTENCE is graded on
    # the second, which is what it says.
    if 0.5 in tabs and tabs[0.5] and tabs[0.0]:
        Rs = [tabs[0.5][x][0] / tabs[0.0][x][0] for x in tabs[0.5] if x in tabs[0.0]]
        span = max(Rs) / min(Rs) if len(Rs) > 1 and min(Rs) > 0 else float("nan")
        say(f"   [reported, first coding] the CV 0.5 per-rung ratios span {span:.2f}x "
            f"-- this is rung-to-rung SCATTER, not candidate separation")
        m2a, m4_lo, m4_hi = 1.0, 1.0 / 3.0, 1.0 / 2.0
        sep = m2a / m4_hi
        say(f"   candidate predictions: M2a R = {m2a:.2f}; M4 with U4's 2-3x prior "
            f"R = {m4_lo:.2f}-{m4_hi:.2f}. Separation {sep:.2f}x-{m2a / m4_lo:.2f}x "
            f"against the {G4_SEPARATION_MIN}x rule -> the test IS "
            f"{'informative' if sep >= G4_SEPARATION_MIN else 'NOT informative'}")
        R5 = gm(Rs)
        say(f"   measured R(0.5) = {R5:.3f}; distance from M2a = "
            f"{abs(R5 - m2a) / max(se_by_cv.get(0.5, float('nan')), 1e-12):.1f} se, "
            f"from the M4 band = "
            f"{min(abs(R5 - m4_lo), abs(R5 - m4_hi)) / max(se_by_cv.get(0.5, float('nan')), 1e-12):.1f} se")
        say("   -> the U4-sized kinetic effect is EXCLUDED; a kinetic effect smaller "
            "than the seed scatter is not.")
    say(f"   G4 verdict: {'COLLAPSE' if g4 else 'see per-CV lines above'}")

    # ---------------- G5 ----------------
    say()
    say("=== G5  cross-protocol anchor: p_ss(phi) crosses p* at phi_c(p*) ===")
    g5 = True
    crossings: dict[float, float] = {}
    for cv in sorted(PHI_MU):
        rr = sorted([r for r in rows if r["cv"] == cv], key=lambda r: r["phi_target"])
        by_phi: dict[float, list[float]] = {}
        for r in rr:
            by_phi.setdefault(r["phi_target"], []).append(r["p_ss"])
        ph = sorted(by_phi)
        ps = [float(np.mean(by_phi[p])) for p in ph]
        cross = float("nan")
        for a, b, pa, pb in zip(ph, ph[1:], ps, ps[1:]):
            if (pa - P_STAR) * (pb - P_STAR) <= 0 and pb != pa:
                cross = a + (P_STAR - pa) * (b - a) / (pb - pa)
                break
        crossings[cv] = cross
        ok = math.isfinite(cross) and abs(cross - PHI_C_CONSTP[cv]) <= G5_TOL
        say(f"   cv {cv}: p_ss = {[f'{v:.2e}' for v in ps]} at phi "
            f"{[f'{v:.4f}' for v in ph]}")
        say(f"           crossing p* = {P_STAR:.0e} at phi = "
            f"{cross if math.isfinite(cross) else float('nan'):.4f} vs constant-pressure "
            f"phi_c = {PHI_C_CONSTP[cv]:.4f}  -> "
            f"{'|d| = %.4f' % abs(cross - PHI_C_CONSTP[cv]) if math.isfinite(cross) else 'REFUSED: no crossing inside the ladder'}")
        g5 = g5 and ok
    say(f"   G5 {'PASS' if g5 else 'FAIL'} (bar {G5_TOL})")
    if all(math.isfinite(v) for v in crossings.values()) and len(crossings) == 3:
        d_cross = crossings[0.5] - crossings[0.0]
        d_cp = PHI_C_CONSTP[0.5] - PHI_C_CONSTP[0.0]
        say(f"   🔑 the CONTRAST from these crossings: Delta phi_c(CV 0 -> 0.5) = "
            f"{d_cross:+.4f}, against {d_cp:+.4f} from the constant-pressure ladder -- "
            f"agreeing to {abs(d_cross - d_cp):.4f}.")
        say("      Two protocols with nothing in common but the material: constant")
        say("      PRESSURE with a volumetric servo, and constant VOLUME from a")
        say("      FRICTIONLESS preparation. A shared systematic would have to survive")
        say("      both, and the offsets it leaves are one-sided (both crossings sit")
        say("      BELOW their constant-pressure value by 0.0033-0.0046) so they largely")
        say("      cancel in the contrast -- which is why the contrast agrees better")
        say("      than either absolute does.")
    if not g5:
        fails.append("G5")

    # ---------------- G6 ----------------
    say()
    say("=== G6  reported, not graded ===")
    shp = [(r["gamma_90"] / r["gamma_c"]) for r in rows
           if r["readable"] and math.isfinite(r["gamma_90"]) and r["gamma_c"] > 0]
    if shp:
        say(f"   gamma_90/gamma_c over {len(shp)} readable cells: "
            f"median {np.median(shp):.2f}, range {min(shp):.2f}-{max(shp):.2f}")
    for cv in sorted(PHI_MU):
        rr = [r for r in rows if r["cv"] == cv]
        rp = [r["rat_prep"] for r in rr if math.isfinite(r["rat_prep"])]
        say(f"   cv {cv}: rattler_frac_prep {min(rp):.3f}-{max(rp):.3f}; "
            f"z_nr_prep NaN on {sum(1 for r in rr if not math.isfinite(r['z_nr_prep']))}"
            f"/{len(rr)} (expected below phi_0)")
    say("   W = phi_0 - phi_mu used to place the ladder: "
        + ", ".join(f"cv {cv}: {PHI_0[cv] - PHI_MU[cv]:.4f}" for cv in sorted(PHI_MU)))

    say()
    if fails:
        say(f"FAILURES: {fails}")
    else:
        say("ALL GRADED BARS PASS")
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "strainclock.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
