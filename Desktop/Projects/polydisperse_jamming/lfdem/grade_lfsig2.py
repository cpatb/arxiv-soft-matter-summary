"""grade_lfsig2.py -- grade PREREG_LFSIG2_2026-09-27.md.

    python lfdem/grade_lfsig2.py            ->  out/grade_lfsig2.txt       (the registered campaign, out/lfsig2/)
    python lfdem/grade_lfsig2.py --dry-lfsig ->  out/grade_lfsig2_dry.txt   (the SAME readouts on the old lfsig
                                                 seed-1 cells: a reproduce-what-is-known check of the grader;
                                                 arms have one seed there, so no paired statistics)

Per cell (out/lfsig2/<label>/): META (LABEL, CV, PHI, SIG, SEED, rc), conf.json (cv_real, a32), data_*.dat.
Label: c<CV*100>_x90_s<sigma>_r<seed>.  Every number printed here is computed from those files.

Written 2026-09-27 ~21:10 CDT, after the queue was launched (20:56) and before any cell finished (the first
cells need >= 25 min). The readout DEFINITIONS are the pre-registration's; this file implements them.
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import os
import re
import sys

import numpy as np
from scipy.optimize import curve_fit

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "out")
W0, W1 = 2.0, 4.0
SIGMAS = (1.0, 2.0, 3.0, 4.5, 7.0, 10.0, 15.0, 30.0, 100.0)
SIGMAS_DRY = (1.0, 3.0, 10.0, 30.0, 100.0)      # the lfsig ladder points >= 1 (0.1, 0.3 are below every readout)
CVS = (0.2, 0.35, 0.5)
SEEDS = (1, 2, 3)
MAX_OVERLAP = 0.05
R_GEO = 0.735 / 0.957                            # dry chain, sigma*(0.5)/sigma*(0.2), d10 (predict_cvmap.txt)
R_INERT = 1.0
R_GEO_35 = math.exp(math.log(R_GEO) * (0.35**2 - 0.2**2) / (0.5**2 - 0.2**2))   # CV^2 interpolation, printed only
DELTA_GEO = 0.815 / 0.938                        # dry chain delta(0.5)/delta(0.2) (predict_cvmap.txt)
SE_MAX = 0.25 * math.log(1.0 / R_GEO)            # G0: a quarter of the candidate separation
T_CRIT = 4.303                                   # df 2, 95 %
X_SENS = {0.2: 0.05 / 0.05, 0.5: 0.09 / 0.05}    # lfsig: d ln sigma_half / d x_m (per unit x_m), stated confound
X_MISMATCH = 0.025                               # LF-relative x of the CV 0.5 arm above the CV 0.2 arm (lfsig P2)

lines: list[str] = []


def say(s: str = "") -> None:
    print(s)
    lines.append(s)


def wc_form(s, ss, d):
    return np.exp(-(ss / s) ** d)


def load_cell(d: str, legacy: bool) -> dict | None:
    meta = open(os.path.join(d, "META"), encoding="utf-8").read()
    kv = dict(re.findall(r"(\w+)=(\S+)", meta))
    files = glob.glob(os.path.join(d, "data_*.dat"))
    if not files:
        return None
    a = np.loadtxt(files[0], comments="#", ndmin=2)
    conf = json.load(open(os.path.join(d, "conf.json"), encoding="utf-8"))
    rc = int(kv["rc"]) if "rc" in kv else None
    if rc is None:   # runner killed: read LF_DEM's own completion line (lfsig Amendment 4 rule)
        logp = os.path.join(d, "run.log")
        done = os.path.exists(logp) and "Job done ok" in open(logp, encoding="utf-8", errors="replace").read()[-400:]
        rc = 0 if done else 124
    strain = a[:, 1]
    w = (strain >= W0) & (strain <= W1)
    c = {"label": kv["LABEL"], "cv": float(kv["CV"]), "phi": float(kv["PHI"]), "sig": float(kv["SIG"]),
         "seed": int(kv.get("SEED", 1)), "rc": rc, "smax": float(strain.max()), "cv_real": conf["cv_real"],
         "a32": conf["a32"], "n_win": int(w.sum()), "raw": a}
    if legacy and "_x90_" not in c["label"]:
        return None
    if w.sum() >= 8:
        blocks = np.array_split(np.where(w)[0], 4)
        for key, col in (("eta", 3), ("zc", 19), ("etac", 4), ("etah", 6), ("etar", 7)):
            bm = np.array([a[b, col].mean() for b in blocks])
            c[key], c[key + "_se"] = float(bm.mean()), float(bm.std(ddof=1) / 2.0)
        c["fcont"] = c["etac"] / c["eta"]
        c["overlap"] = float(-a[w, 15].min())
    else:
        for key in ("eta", "zc", "etac", "etah", "etar", "fcont", "overlap"):
            c[key] = float("nan")
    c["complete"] = c["smax"] >= W1 - 0.05
    c["healthy"] = rc in (0, 124) and c["complete"] and c["overlap"] < MAX_OVERLAP
    return c


def interp_log(s: np.ndarray, y: np.ndarray, thr: float):
    """Stress at which y first crosses thr, log-linear in sigma between the bracketing rungs (NaN if none)."""
    for i in range(1, len(y)):
        if y[i - 1] < thr <= y[i]:
            t = (thr - y[i - 1]) / (y[i] - y[i - 1])
            return float(math.exp(math.log(s[i - 1]) + t * (math.log(s[i]) - math.log(s[i - 1]))))
    return float("nan")


def interp_at(s: np.ndarray, y: np.ndarray, s0: float) -> float:
    """y log-linearly interpolated in sigma at s0."""
    if not (s[0] <= s0 <= s[-1]):
        return float("nan")
    return float(np.interp(math.log(s0), np.log(s), y))


def arm_readouts(cells: dict, cv: float, seed: int, sigmas) -> dict:
    r = {"cv": cv, "seed": seed, "complete": all((cv, seed, s) in cells for s in sigmas)}
    if not r["complete"]:
        r["missing"] = [s for s in sigmas if (cv, seed, s) not in cells]
        return r
    s = np.array(sigmas, dtype=float)
    z = np.array([cells[(cv, seed, x)]["zc"] for x in sigmas])
    eta = np.array([cells[(cv, seed, x)]["eta"] for x in sigmas])
    fc = np.array([cells[(cv, seed, x)]["fcont"] for x in sigmas])
    f = z / z[-1]
    r["z"], r["f"], r["eta"], r["fcont"] = z, f, eta, fc
    r["a32"] = cells[(cv, seed, sigmas[0])]["a32"]
    r["cv_real"] = cells[(cv, seed, sigmas[0])]["cv_real"]
    # R1: WC activation form
    try:
        p, cov = curve_fit(wc_form, s, f, p0=(5.0, 1.5), maxfev=20000)
        resid = float(np.abs(wc_form(s, *p) - f).max())
        r["sstar"], r["delta"], r["resid"] = float(p[0]), float(p[1]), resid
        r["sstar_se"], r["delta_se"] = [float(x) for x in np.sqrt(np.diag(cov))]
        r["r1_ok"] = resid < 0.05
    except Exception as e:  # noqa: BLE001
        r["sstar"] = r["delta"] = r["resid"] = float("nan"); r["r1_ok"] = False; r["fit_err"] = str(e)
    # R2: midpoint (requires monotone non-decreasing f)
    r["mono"] = bool(np.all(np.diff(z) >= 0))
    r["shalf"] = interp_log(s, f, 0.5) if r["mono"] else float("nan")
    # R3: absolute Z_c = 1
    r["sz1"] = interp_log(s, z, 1.0) if r["mono"] else float("nan")
    # R4: contact share at the midpoint
    r["fcont_half"] = interp_at(s, fc, r["shalf"]) if np.isfinite(r["shalf"]) else float("nan")
    # R5: thickening ratio
    r["thick"] = float(eta[-1] / eta[0])
    return r


def paired(arms: dict, key: str, cv_a: float, cv_b: float, seeds) -> dict:
    """Paired per-seed ln(readout(cv_b)/readout(cv_a)) with mean, sd, se, t (n = seeds present in both)."""
    vals = []
    for sd in seeds:
        a, b = arms.get((cv_a, sd)), arms.get((cv_b, sd))
        if a and b and a.get("complete") and b.get("complete") and np.isfinite(a.get(key, np.nan)) and np.isfinite(b.get(key, np.nan)):
            vals.append(math.log(b[key] / a[key]))
    v = np.array(vals)
    out = {"n": len(v), "vals": v}
    if len(v) >= 2:
        out["mean"] = float(v.mean()); out["sd"] = float(v.std(ddof=1)); out["se"] = out["sd"] / math.sqrt(len(v))
        out["t"] = out["mean"] / out["se"] if out["se"] > 0 else float("inf")
    elif len(v) == 1:
        out["mean"] = float(v[0]); out["sd"] = out["se"] = out["t"] = float("nan")
    return out


def fmt(x, w=7, p=3):
    return f"{x:{w}.{p}f}" if isinstance(x, (float, int, np.floating)) and np.isfinite(x) else f"{'nan':>{w}s}"


def wiring_gate(cells2: dict) -> None:
    """P4: the ten seed-1 cells lfsig already ran must reproduce its data rows bit-identically for strain <= 4."""
    say("\n==== P4 wiring: seed-1 cells vs lfsig (bit-identity of every data row with strain <= 4) ====")
    n_ok = n_tot = n_diff = 0
    for cv, tag in ((0.2, "c20"), (0.5, "c50")):
        for s in (1, 3, 10, 30, 100):
            n_tot += 1
            old_d = os.path.join(OUT, "lfsig", f"{tag}_x90_s{s:g}")
            new = cells2.get((cv, 1, float(s)))
            fo = glob.glob(os.path.join(old_d, "data_*.dat"))
            if new is None or not fo:
                say(f"  {tag}_x90_s{s:g}: NOT GRADED (new cell {'present' if new else 'missing'}, old file {'present' if fo else 'missing'})")
                continue
            old = np.loadtxt(fo[0], comments="#", ndmin=2)
            new_a = new["raw"]
            mo, mn = old[old[:, 1] <= 4.0 + 1e-9], new_a[new_a[:, 1] <= 4.0 + 1e-9]
            n = min(len(mo), len(mn))
            same = n > 0 and mo.shape[1] == mn.shape[1] and np.array_equal(mo[:n], mn[:n])
            n_ok += int(same); n_diff += int(not same)
            say(f"  {tag}_x90_s{s:g}: rows old {len(mo)} new {len(mn)} compared {n}: {'IDENTICAL' if same else 'DIFFER'}"
                + ("" if same else f" (max |d eta| {np.abs(mo[:n, 3] - mn[:n, 3]).max():.3e})"))
    verdict = "PASS" if n_ok == n_tot else ("FAIL" if n_diff else "INCOMPLETE")
    say(f"  P4 verdict: {verdict} ({n_ok}/{n_tot} identical, {n_diff} differ) -- any DIFFER line above is a wiring fault")


SIGMAS_X = (3.0, 4.5, 7.0, 10.0, 100.0)
X1_PRED = -0.047          # Amendment 2: predicted ln[R2(0.5x)/R2(0.5)] from the lfsig x-sensitivity


def x_arm(cells: dict, arms: dict, con: dict) -> None:
    """Amendment 2: the x-matched CV 0.5 arm (out/lfsig2x/, labels c50x_x85_s<sigma>_r<seed>, phi 0.55643).
    X1 = paired ln[R2(0.5x)/R2(0.5)] (the confound's own size); X2 = R on the x-matched pair, 5-rung R1 and R2.
    The registered 0.2 and 0.5 arms are refitted on the SAME five rungs for X2 (like for like)."""
    say("\n==== Amendment 2: x-matched CV 0.5 arm (phi 0.55643, five rungs) ====")
    dirs = [x for x in sorted(glob.glob(os.path.join(OUT, "lfsig2x", "c50x*"))) if ".killed_" not in x]
    xc = [c for c in (load_cell(d, False) for d in dirs) if c]
    if not xc:
        say("  no x-arm cells yet")
        return
    say(f"  {'label':20s} rc {'strain':>6s} {'eta':>8s} {'Zc':>6s} {'fcont':>6s} {'ovlp':>7s} ok")
    for c in sorted(xc, key=lambda c: (c["seed"], c["sig"])):
        say(f"  {c['label']:20s} {c['rc']:2d} {c['smax']:6.2f} {fmt(c['eta'], 8)} {fmt(c['zc'], 6)} {fmt(c['fcont'], 6)} {fmt(c['overlap'], 7, 4)} {int(c['healthy'])}")
    xcells = {(0.5, c["seed"], c["sig"]): c for c in xc if c["healthy"]}
    xarms = {sd: arm_readouts(xcells, 0.5, sd, SIGMAS_X) for sd in SEEDS}
    # the registered arms on the same five rungs
    sub = {k: v for k, v in cells.items() if k[2] in SIGMAS_X}
    a5 = {(cv, sd): arm_readouts(sub, cv, sd, SIGMAS_X) for cv in (0.2, 0.5) for sd in SEEDS}
    say(f"  {'arm':8s} {'seed':>4s} {'sigma*5':>8s} {'delta5':>7s} {'resid':>6s} {'s_half':>7s}")
    for sd in SEEDS:
        for name, r in (("0.2 (5)", a5[(0.2, sd)]), ("0.5 (5)", a5[(0.5, sd)]), ("0.5x", xarms[sd])):
            if r["complete"]:
                say(f"  {name:8s} {sd:4d} {fmt(r['sstar'], 8)} {fmt(r['delta'], 7)} {fmt(r['resid'], 6)} {fmt(r['shalf'])}")
            else:
                say(f"  {name:8s} {sd:4d}  INCOMPLETE (missing {r['missing']})")
    # X1: paired ln R2(0.5x)/R2(0.5), registered 0.5 arm's 9-rung sigma_half against the x arm's 5-rung one
    # (both interpolate between the same bracketing rungs 3..10 when the crossing lies there, so like for like)
    v = []
    for sd in SEEDS:
        a, b = arms.get((0.5, sd)), xarms[sd]
        if a and a.get("complete") and b["complete"] and np.isfinite(a["shalf"]) and np.isfinite(b["shalf"]):
            v.append(math.log(b["shalf"] / a["shalf"]))
    if len(v) == 3:
        v = np.array(v); m, se = float(v.mean()), float(v.std(ddof=1) / math.sqrt(3))
        say(f"  X1: ln[R2(0.5x)/R2(0.5)] = {m:+.4f} se {se:.4f} (ratio {math.exp(m):.3f}, band [{math.exp(m - se):.3f}, {math.exp(m + se):.3f}]); "
            f"predicted {X1_PRED:+.3f}; per seed " + " ".join(f"{x:+.3f}" for x in v))
        p = con.get((0.2, 0.5, "sstar"), {"n": 0})
        if p["n"] == 3 and p["se"] < SE_MAX:
            mR, seR = p["mean"], p["se"]
            # margin from the raw ln R to the nearest candidate boundary (ln R_cand +- 2 se)
            margins = []
            for rc_ in (R_GEO, R_INERT):
                lo, hi = math.log(rc_) - 2 * seR, math.log(rc_) + 2 * seR
                margins.append(min(abs(mR - lo), abs(mR - hi)))
            margin = min(margins)
            say(f"  X1 vs P1 margin: |X1| {abs(m):.4f} vs the raw ln R's distance to the nearest candidate boundary {margin:.4f} -> "
                f"confound {'CLOSED' if abs(m) < margin else 'NOT closed: P1 verdict is CONDITIONAL on the x-mismatch'}")
        else:
            say("  X1 vs P1 margin: P1 not graded / not informative, nothing to close")
    else:
        say(f"  X1: NOT GRADED ({len(v)}/3 seeds)")
    # X2: R on the x-matched pair, five rungs
    for key, name in (("sstar", "5-rung R1 sigma*"), ("shalf", "5-rung R2 sigma_half")):
        vals = [math.log(xarms[sd][key] / a5[(0.2, sd)][key]) for sd in SEEDS
                if xarms[sd]["complete"] and a5[(0.2, sd)]["complete"] and np.isfinite(xarms[sd][key]) and np.isfinite(a5[(0.2, sd)][key])]
        reg = [math.log(a5[(0.5, sd)][key] / a5[(0.2, sd)][key]) for sd in SEEDS
               if a5[(0.5, sd)]["complete"] and a5[(0.2, sd)]["complete"] and np.isfinite(a5[(0.5, sd)][key]) and np.isfinite(a5[(0.2, sd)][key])]
        if len(vals) >= 2 and len(reg) >= 2:
            vx, vr = np.array(vals), np.array(reg)
            say(f"  X2 {name:22s}: R(0.5x/0.2) {math.exp(vx.mean()):.3f} band [{math.exp(vx.mean() - vx.std(ddof=1) / math.sqrt(len(vx))):.3f}, "
                f"{math.exp(vx.mean() + vx.std(ddof=1) / math.sqrt(len(vx))):.3f}] (n {len(vx)})  vs registered pair on the same rungs "
                f"R(0.5/0.2) {math.exp(vr.mean()):.3f} (n {len(vr)})")
        else:
            say(f"  X2 {name:22s}: NOT GRADED ({len(vals)} x pairs, {len(reg)} registered pairs)")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-lfsig", action="store_true", help="run the readouts on out/lfsig (seed-1 legacy cells)")
    args = ap.parse_args()
    legacy = args.dry_lfsig
    sub = "lfsig" if legacy else "lfsig2"
    sigmas = SIGMAS_DRY if legacy else SIGMAS
    say(f"grade_lfsig2 -- PREREG_LFSIG2_2026-09-27.md   [{'DRY RUN on lfsig seed-1 cells' if legacy else 'registered campaign'}]")
    dirs = [x for x in sorted(glob.glob(os.path.join(OUT, sub, "c*"))) if ".killed_" not in x]
    cells_all = [c for c in (load_cell(d, legacy) for d in dirs) if c]
    n_reg = len(CVS) * len(SEEDS) * len(SIGMAS)
    say(f"cells with data: {len(cells_all)}" + ("" if legacy else f" of {n_reg} registered"))
    say(f"  {'label':18s} {'CVreal':>6s} {'a32':>6s} rc {'strain':>6s} {'n_win':>5s} {'eta (se)':>18s} {'Zc (se)':>16s} {'fcont':>6s} {'ovlp':>7s} ok")
    for c in sorted(cells_all, key=lambda c: (c["cv"], c["seed"], c["sig"])):
        say(f"  {c['label']:18s} {c['cv_real']:6.3f} {c['a32']:6.3f} {c['rc']:2d} {c['smax']:6.2f} {c['n_win']:5d} "
            f"{fmt(c['eta'], 9, 3)} ({fmt(c['eta_se'], 5, 3)}) {fmt(c['zc'], 7, 3)} ({fmt(c['zc_se'], 5, 3)}) "
            f"{fmt(c['fcont'], 6, 3)} {fmt(c['overlap'], 7, 4)} {int(c['healthy'])}")
    cells = {(c["cv"], c["seed"], c["sig"]): c for c in cells_all if c["healthy"]}
    unhealthy = [c["label"] for c in cells_all if not c["healthy"]]
    say(f"healthy (registered rule: rc 0/124, strain >= {W1 - 0.05:.2f}, overlap < {MAX_OVERLAP}): {len(cells)}; excluded: {unhealthy or 'none'}")

    seeds = (1,) if legacy else SEEDS
    cvs = (0.2, 0.5) if legacy else CVS
    arms = {(cv, sd): arm_readouts(cells, cv, sd, sigmas) for cv in cvs for sd in seeds}
    say("\n==== per-arm readouts (R1 WC form sigma*, delta; R2 sigma_half; R3 sigma at Z_c = 1; R4 contact share at sigma_half; R5 eta(100)/eta(1)) ====")
    say(f"  {'CV':>4s} {'seed':>4s} {'sigma*':>7s} {'(se)':>6s} {'delta':>6s} {'(se)':>6s} {'resid':>6s} R1ok {'mono':>4s} {'s_half':>7s} {'s_Z1':>7s} {'fc_half':>7s} {'thick':>6s}  f(sigma) at the ladder")
    for (cv, sd), r in sorted(arms.items()):
        if not r["complete"]:
            say(f"  {cv:4.2f} {sd:4d}  INCOMPLETE -- missing sigma {r['missing']}")
            continue
        say(f"  {cv:4.2f} {sd:4d} {fmt(r['sstar'])} {fmt(r['sstar_se'], 6)} {fmt(r['delta'], 6)} {fmt(r['delta_se'], 6)} {fmt(r['resid'], 6)} "
            f"{int(r['r1_ok'])}    {int(r['mono']):4d} {fmt(r['shalf'])} {fmt(r['sz1'])} {fmt(r['fcont_half'])} {fmt(r['thick'], 6, 2)}  "
            + " ".join(f"{v:.3f}" for v in r["f"]))

    # ---- paired contrasts -------------------------------------------------------------------------------
    say("\n==== paired per-seed contrasts, ln(readout(CV_b)/readout(CV_a)); mean, sd, se = sd/sqrt(n), t (df n-1) ====")
    con = {}
    for cv_a, cv_b in ((0.2, 0.5), (0.2, 0.35), (0.35, 0.5)):
        if cv_a not in cvs or cv_b not in cvs:
            continue
        for key, name in (("sstar", "R1 sigma*_WC"), ("delta", "R1 delta_WC"), ("shalf", "R2 sigma_half"),
                          ("sz1", "R3 sigma_Z1"), ("fcont_half", "R4 fcont(sigma_half)"), ("thick", "R5 thickening")):
            p = paired(arms, key, cv_a, cv_b, seeds)
            con[(cv_a, cv_b, key)] = p
            if p["n"] == 0:
                say(f"  {cv_a:.2f}->{cv_b:.2f} {name:22s}: no complete pair")
                continue
            vals = " ".join(f"{math.exp(v):.3f}" for v in p["vals"])
            if p["n"] >= 2:
                say(f"  {cv_a:.2f}->{cv_b:.2f} {name:22s}: ratio {math.exp(p['mean']):.3f}  ln mean {p['mean']:+.4f} sd {p['sd']:.4f} se {p['se']:.4f} "
                    f"t {p['t']:+.2f} (n {p['n']}; band [{math.exp(p['mean'] - p['se']):.3f}, {math.exp(p['mean'] + p['se']):.3f}])  per seed: {vals}")
            else:
                say(f"  {cv_a:.2f}->{cv_b:.2f} {name:22s}: ratio {math.exp(p['mean']):.3f}  (n 1, no band)  per seed: {vals}")

    # ---- convention factor ---------------------------------------------------------------------------
    a32 = {cv: np.mean([arms[(cv, sd)]["a32"] for sd in seeds if arms[(cv, sd)].get("complete")]) for cv in cvs}
    if all(np.isfinite(a32.get(cv, np.nan)) for cv in (0.2, 0.5)):
        fac = a32[0.5] / a32[0.2]   # FIRST power (prereg Erratum 1): F_R scales with the reduced radius, so sigma_0 = F_R0/a
        say(f"\n  convention: a32 (number-mean radius = 1) CV 0.2 {a32[0.2]:.4f}, CV 0.5 {a32[0.5]:.4f}; R_d32 = R_d10 x {fac:.3f} "
            f"(first power, sampled truncated a32; predict_cvmap used the analytic (1+CV^2)^2 = 1.558 for 0.5 vs 0)")

    # ---- P1 ----------------------------------------------------------------------------------------------
    say("\n==== P1: R = sigma*_WC(0.5)/sigma*_WC(0.2); H_geo 0.768, H_inert 1.000 ====")
    p = con.get((0.2, 0.5, "sstar"), {"n": 0})
    if p["n"] < 3:
        say(f"  P1: NOT GRADED ({p['n']}/3 seed pairs complete)")
    else:
        se = p["se"]; m = p["mean"]
        say(f"  G0 resolving power: se(ln R) = {se:.4f} vs bar {SE_MAX:.4f} -> {'INFORMATIVE' if se < SE_MAX else 'NOT INFORMATIVE'}")
        if se < SE_MAX:
            geo = abs(m - math.log(R_GEO)) <= 2 * se
            ine = abs(m - math.log(R_INERT)) <= 2 * se
            verdict = "BOTH" if geo and ine else "GEO" if geo else "INERT" if ine else "NEITHER"
            say(f"  R = {math.exp(m):.3f}, 2 se band [{math.exp(m - 2 * se):.3f}, {math.exp(m + 2 * se):.3f}]; "
                f"|ln R - ln 0.768| = {abs(m - math.log(R_GEO)):.4f} ({abs(m - math.log(R_GEO)) / se:.1f} se); "
                f"|ln R - ln 1| = {abs(m):.4f} ({abs(m) / se:.1f} se)")
            say(f"  P1 verdict: {verdict}")
            # stated confound: the LF-relative x mismatch
            dl = (X_SENS[0.5] - 0.0) * X_MISMATCH  # the CV 0.5 arm sits ~0.025 higher in x_LF: its sigma_half is ~ +9 %/0.05 x that
            say(f"  stated confound (x_LF mismatch {X_MISMATCH}, lfsig sensitivity {X_SENS[0.5] * 0.05:.2f}/0.05 x_m at CV 0.5): "
                f"R bracket [{math.exp(m - dl):.3f}, {math.exp(m):.3f}] (verdict is on the raw R)"
                + ("" if (abs(m - dl - math.log(R_GEO)) <= 2 * se) == geo and (abs(m - dl) <= 2 * se) == ine
                   else "  ** the bracket straddles a candidate boundary **"))
        else:
            say(f"  R = {math.exp(m):.3f} (se {se:.4f}); no candidate graded")

    # ---- P2 ----------------------------------------------------------------------------------------------
    say("\n==== P2: delta_WC falls with CV (3/3 seeds paired, mean - 2 se < 0); dry-DEM reference 0.869 ====")
    p = con.get((0.2, 0.5, "delta"), {"n": 0})
    if p["n"] < 3:
        say(f"  P2: NOT GRADED ({p['n']}/3)")
    else:
        neg = int(np.sum(p["vals"] < 0)); m, se = p["mean"], p["se"]
        ok = neg == 3 and (m - 2 * se) < 0
        say(f"  seeds with ln ratio < 0: {neg}/3; mean - 2 se = {m - 2 * se:+.4f}; ratio {math.exp(m):.3f} band [{math.exp(m - 2 * se):.3f}, {math.exp(m + 2 * se):.3f}]")
        say(f"  P2 verdict: {'PASS' if ok else 'FAIL'}; dry-DEM 0.869 {'CONSISTENT' if abs(m - math.log(DELTA_GEO)) <= 2 * se else 'NOT consistent'} "
            f"(|ln diff| {abs(m - math.log(DELTA_GEO)):.4f} = {abs(m - math.log(DELTA_GEO)) / se if se > 0 else float('inf'):.1f} se)")

    # ---- P3 ----------------------------------------------------------------------------------------------
    say("\n==== P3: contact share at sigma_half rises with CV, 0.2 -> 0.5, 3/3 seeds paired ====")
    p = con.get((0.2, 0.5, "fcont_half"), {"n": 0})
    if p["n"] < 3:
        say(f"  P3: NOT GRADED ({p['n']}/3)")
    else:
        pos = int(np.sum(p["vals"] > 0))
        say(f"  seeds with fcont(0.5) > fcont(0.2): {pos}/3; per seed ratios " + " ".join(f"{math.exp(v):.3f}" for v in p["vals"]))
        say(f"  P3 verdict: {'PASS' if pos == 3 else 'FAIL'}")

    # ---- P4 ----------------------------------------------------------------------------------------------
    if not legacy:
        wiring_gate(cells)

    # ---- P5 ----------------------------------------------------------------------------------------------
    if 0.35 in cvs:
        say("\n==== P5 (descriptive): sigma*_WC(0.35) between the 0.2 and 0.5 values, >= 2/3 seeds ====")
        k = 0; n = 0
        for sd in seeds:
            a, b, c = arms.get((0.2, sd)), arms.get((0.35, sd)), arms.get((0.5, sd))
            if all(x and x.get("complete") and np.isfinite(x.get("sstar", np.nan)) for x in (a, b, c)):
                n += 1
                between = min(a["sstar"], c["sstar"]) <= b["sstar"] <= max(a["sstar"], c["sstar"])
                k += int(between)
                say(f"  seed {sd}: sigma* {a['sstar']:.3f} / {b['sstar']:.3f} / {c['sstar']:.3f}  between: {between}")
        p35 = con.get((0.2, 0.35, "sstar"), {"n": 0})
        if p35.get("n", 0) >= 1:
            say(f"  R(0.35/0.2) = {math.exp(p35['mean']):.3f}  (H_geo CV^2-interpolated {R_GEO_35:.3f}, printed only)")
        say(f"  P5: {'PASS' if n == 3 and k >= 2 else 'NOT GRADED' if n < 3 else 'FAIL'} ({k}/{n} seeds between)")

    # ---- P6 ----------------------------------------------------------------------------------------------
    say("\n==== P6 (for the record): the lfsig post-hoc readouts on 3-seed data ====")
    for key, name, old in (("shalf", "R2 self-normalised sigma_half", "0.956"), ("sz1", "R3 absolute Z_c = 1", "1.107")):
        p = con.get((0.2, 0.5, key), {"n": 0})
        if p["n"] >= 2:
            say(f"  {name}: ratio {math.exp(p['mean']):.3f} band [{math.exp(p['mean'] - p['se']):.3f}, {math.exp(p['mean'] + p['se']):.3f}] (lfsig single-seed: {old})")
        elif p["n"] == 1:
            say(f"  {name}: ratio {math.exp(p['mean']):.3f} (n 1) (lfsig single-seed: {old})")
        else:
            say(f"  {name}: no pair")

    if not legacy:
        x_arm(cells, arms, con)

    outp = os.path.join(OUT, "grade_lfsig2_dry.txt" if legacy else "grade_lfsig2.txt")
    with open(outp, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwrote {outp}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
