"""grade_mixture_ext.py -- grade the EXTENSION block, PREREG_MIXTURE_EXT_2026-09-14.md.

Separate from `grade_mixture.py` / `grade_mixture_mu0.py` for the same reason those two
are separate: M0-M6 and F1-F3 are CLOSED gates and this script must not be able to
re-decide them. It loads ONLY `phase12_extr*_shard*.jsonl` (ratios 4 and 5) plus the
ratio <= 3 binaries through `mixture_closure.binary_points`, which is the calibration
set the predictions were written from.

    python grade_mixture_ext.py            grade X0-X5 and print the registered analysis
    python grade_mixture_ext.py --predict  print the blind table only (no data loaded)
"""

from __future__ import annotations

import collections
import glob
import json
import os
import sys

import numpy as np
from scipy.interpolate import PchipInterpolator

import mixture_closure as mc
import mixture_model as mm
from measure import RATTLER_THRESHOLDS
from phase12_mixture import MIN_PHASE_COUNT, N_PROD

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
SCALE = mm.PHI_SCALE
MU = 0.3
GLOB = "phase12_extr*_shard*.jsonl"          # NOT phase12_mixr*: see module docstring

# --- registered bars, PREREG_MIXTURE_EXT_2026-09-14.md section 4 --------------------
X0_BAR = 2.0                 # scale; monodisperse endpoints vs beta
X1_RMS_BAR = 4.0             # scale; 1.4 x the in-sample rms 2.83
X1_MAX_BAR = 6.7             # scale; the in-sample maximum |residual|
X2_POINTS = ((4.0, 0.20), (4.0, 0.25), (5.0, 0.15))
X2_SEP_BAR = 1.6             # rms ratio loser/winner for "informative"
X3_STEP_BAR = 3.0            # scale; every 0.05 step at r = 4 and the r4-over-r3 max
X4_SLACK = 1.0               # scale
X5_RAT_L_BAR = 0.10
GRID4 = (0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 1.0)
GRID5 = (0.0, 0.05, 0.10, 0.15, 1.0)
R3_MAX_MEAS = None           # filled from the ratio-3 data at run time, printed


def load():
    g: dict[tuple[float, float], list[dict]] = collections.defaultdict(list)
    files = sorted(glob.glob(os.path.join(OUT, GLOB)))
    nraw = 0
    for f in files:
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                nraw += 1
                r = json.loads(line)
                key = (float(r["mix_ratio"]), round(float(r["x_small_vol"]), 6))
                g[key].append(r)
    return g, len(files), nraw


def admissible(c: dict) -> tuple[bool, str]:
    if c.get("error"):
        return False, "error"
    if not c.get("converged"):
        return False, "not converged"
    if c.get("n_requested") != N_PROD:
        return False, f"N {c.get('n_requested')}"
    if abs(float(c.get("mu_p", -1)) - MU) > 1e-9:
        return False, f"mu {c.get('mu_p')}"
    bc = c.get("by_class", {})
    x = float(c["x_small_vol"])
    if 0.0 < x < 1.0 and min(bc.get("n_S", 0), bc.get("n_L", 0)) < MIN_PHASE_COUNT:
        return False, f"phase count {bc.get('n_S')}/{bc.get('n_L')}"
    return True, ""


def candidates(beta: float, pts):
    """(PA, PB) for J, R3, EXT exactly as registered in section 3."""
    paJ, pbJ, _ = mc.fit_on_binaries(pts, beta)
    per = {}
    for r in sorted({p[0] for p in pts}):
        sub = [p for p in pts if abs(p[0] - r) < 1e-9]
        per[r] = mc.fit_on_binaries(sub, beta)[:2]
    u = np.array([1.0 / r for r in per])
    pa = np.array([per[r][0] for r in per])
    pb = np.array([per[r][1] for r in per])
    sa, ia = np.polyfit(u, pa, 1)
    sb, ib = np.polyfit(u, pb, 1)
    out = {"J": lambda r: (paJ, pbJ),
           "R3": lambda r, p=per[3.0]: p,
           "EXT": lambda r: (ia + sa / r, ib + sb / r)}
    return out, per, (ia, sa, ib, sb)


def ext_points(beta: float):
    """Seed-mean (ratio, x, phi) for every admissible interior cell of the ext block,
    appended to the ratio <= 3 calibration points. What section 5 fits."""
    g, _, _ = load()
    pts = mc.binary_points(mu_p=MU)
    out = list(pts)
    for (r, x), cs in sorted(g.items()):
        if not (0.0 < x < 1.0):
            continue
        ok = [c for c in cs if admissible(c)[0]]
        if ok:
            out.append((r, x, float(np.mean([c["phi"] for c in ok]))))
    return out


def tabulated_closure(beta: float, allpts):
    """The nonparametric closure of section 5: direct (a, b) per ratio, kept only where
    IDENTIFIABLE, joined by shape-preserving interpolation to the two proved limits.

    Identifiable = both profile intervals narrower than 0.25 (the E-section's rule,
    verbatim) AND containing their own best fit (added after the r = 5 profile returned
    an interval that excluded its Nelder-Mead optimum: three coarse-branch points cannot
    place b, and the grid found a cusp where the model's crossing passes through the
    last point). Returns ({u: (a, b)}, a(u), b(u)).
    """
    tab: dict[float, tuple[float, float]] = {}
    for r in sorted({p[0] for p in allpts}):
        sub = [p for p in allpts if abs(p[0] - r) < 1e-9 and 0.0 < p[1] < 1.0]
        a_, b_, _ = mc.fit_ab_direct(sub, beta)
        alo, ahi = mc.profile_interval(sub, beta, 0)
        blo, bhi = mc.profile_interval(sub, beta, 1)
        if (ahi - alo < 0.25 and bhi - blo < 0.25
                and alo - 1e-9 <= a_ <= ahi + 1e-9 and blo - 1e-9 <= b_ <= bhi + 1e-9):
            tab[1.0 / r] = (a_, b_)
    us = [0.0, 1.0] + sorted(tab)
    aa = [0.0, 1.0] + [tab[u][0] for u in sorted(tab)]
    bb = [0.0, 1.0] + [tab[u][1] for u in sorted(tab)]
    order = np.argsort(us)
    uu = np.asarray(us)[order]
    fa = PchipInterpolator(uu, np.asarray(aa)[order])
    fb = PchipInterpolator(uu, np.asarray(bb)[order])
    return tab, fa, fb


def main() -> int:
    predict_only = "--predict" in sys.argv
    lines: list[str] = []
    failures: list[str] = []

    def say(s: str = "") -> None:
        lines.append(s)
        print(s)

    meas3 = mm.measured(mu_p=MU)
    beta = meas3[("lognormal", 0.0)][0]
    pts = mc.binary_points(mu_p=MU)
    cand, per, (ia, sa, ib, sb) = candidates(beta, pts)
    say("=== grade_mixture_ext.py -- PREREG_MIXTURE_EXT_2026-09-14.md ===")
    say(f"   beta = {beta:.4f} (measured, mu_p {MU}); calibration set {len(pts)} binary "
        f"points at ratios {sorted({p[0] for p in pts})}")
    say("   per-ratio (PA, PB): " + ", ".join(f"r{r:.2f} ({a:.3f}, {b:.3f})"
                                             for r, (a, b) in per.items()))
    say(f"   J = {cand['J'](0)[0]:.3f}, {cand['J'](0)[1]:.3f};  R3 = {cand['R3'](0)[0]:.3f}, "
        f"{cand['R3'](0)[1]:.3f};  EXT: PA(u) = {ia:.3f} + {sa:.3f} u, PB(u) = {ib:.3f} + {sb:.3f} u")
    say()

    say("=== blind predictions (section 3) ===")
    say(f"   {'r':>4}{'X_s':>6}{'J':>9}{'R3':>9}{'EXT':>9}{'ideal':>9}{'J-EXT/scale':>13}")
    pred: dict[tuple[float, float], dict[str, float]] = {}
    for r, grid in ((4.0, GRID4), (5.0, GRID5)):
        for x in grid:
            d = {k: mm.binary_density(r, x, beta, *f(r)) for k, f in cand.items()}
            d["ideal"] = mm.furnas_ideal(x, beta)
            pred[(r, x)] = d
            say(f"   {r:4.1f}{x:6.2f}{d['J']:9.4f}{d['R3']:9.4f}{d['EXT']:9.4f}"
                f"{d['ideal']:9.4f}{(d['J'] - d['EXT']) / SCALE:13.1f}")
    say()
    if predict_only:
        return 0

    g, nfiles, nraw = load()
    say(f"=== data: {nraw} records from {nfiles} shard files -> {len(g)} compositions ===")
    adm: dict[tuple[float, float], list[dict]] = {}
    for key in sorted(g):
        ok = []
        for c in g[key]:
            a, why = admissible(c)
            if a:
                ok.append(c)
            else:
                say(f"   REFUSED {c.get('key')}: {why}")
        adm[key] = ok
    nadm = sum(len(v) for v in adm.values())
    say(f"   admitted {nadm} of {nraw}")
    # the ratio-3 block, for X3/X4 (seed means from the calibration set, N = 1500)
    r3 = {round(p[1], 6): p[2] for p in pts if abs(p[0] - 3.0) < 1e-9}
    say(f"   ratio-3 seed means available at X_s = {sorted(r3)}")
    say()

    def mean_sd(key):
        v = np.array([c["phi"] for c in adm[key]], dtype=float)
        return float(v.mean()), (float(v.std(ddof=1)) if v.size > 1 else float("nan")), v.size

    # ---------------------------------------------------------------- X0
    say("=== X0 -- monodisperse controls ===")
    x0_ok = True
    for key in ((4.0, 0.0), (4.0, 1.0), (5.0, 0.0), (5.0, 1.0)):
        m, sd, n = mean_sd(key)
        dev = (m - beta) / SCALE
        say(f"   r{key[0]:.0f} X_s {key[1]:.0f}: phi {m:.4f} +- {sd:.4f} (n {n})  "
            f"vs beta {beta:.4f}: {dev:+.1f} scale")
        x0_ok &= abs(dev) <= X0_BAR
    say(f"   -> X0 {'PASS' if x0_ok else 'FAIL'} (bar {X0_BAR} scale)")
    if not x0_ok:
        failures.append("X0")
    say()

    # ---------------------------------------------------------------- X1
    say("=== X1 -- blind transfer of the shipped closure J ===")
    say(f"   {'r':>4}{'X_s':>6}{'phi':>9}{'sd':>8}{'n':>3}{'J':>9}{'res/scale':>11}"
        f"{'R3':>9}{'res':>7}{'EXT':>9}{'res':>7}")
    resJ, resR3, resE = {}, {}, {}
    meas: dict[tuple[float, float], float] = {}
    for r, grid in ((4.0, GRID4), (5.0, GRID5)):
        for x in grid:
            if x in (0.0, 1.0):
                continue
            key = (r, x)
            m, sd, n = mean_sd(key)
            meas[key] = m
            resJ[key] = (m - pred[key]["J"]) / SCALE
            resR3[key] = (m - pred[key]["R3"]) / SCALE
            resE[key] = (m - pred[key]["EXT"]) / SCALE
            say(f"   {r:4.1f}{x:6.2f}{m:9.4f}{sd:8.4f}{n:3d}{pred[key]['J']:9.4f}"
                f"{resJ[key]:11.1f}{pred[key]['R3']:9.4f}{resR3[key]:7.1f}"
                f"{pred[key]['EXT']:9.4f}{resE[key]:7.1f}")
    rj = np.array(list(resJ.values()))
    rms_j, max_j = float(np.sqrt((rj ** 2).mean())), float(np.abs(rj).max())
    say(f"   J: rms {rms_j:.2f} (bar {X1_RMS_BAR}), max |res| {max_j:.2f} (bar {X1_MAX_BAR}), "
        f"mean signed {rj.mean():+.2f}, {int((rj > 0).sum())} of {rj.size} positive")
    for name, rr in (("R3", resR3), ("EXT", resE)):
        a = np.array(list(rr.values()))
        say(f"   {name}: rms {np.sqrt((a ** 2).mean()):.2f}, max |res| {np.abs(a).max():.2f}, "
            f"mean signed {a.mean():+.2f}  (reported, not graded here)")
    x1_ok = rms_j <= X1_RMS_BAR and max_j <= X1_MAX_BAR
    say(f"   -> X1 {'PASS' if x1_ok else 'FAIL'}")
    if not x1_ok:
        failures.append("X1")
    say()

    # ---------------------------------------------------------------- X2
    say("=== X2 -- which way does the form bend? (J vs EXT on the 3 separated points) ===")
    rj2 = np.array([resJ[k] for k in X2_POINTS])
    re2 = np.array([resE[k] for k in X2_POINTS])
    rms_j2, rms_e2 = float(np.sqrt((rj2 ** 2).mean())), float(np.sqrt((re2 ** 2).mean()))
    seps = [abs(pred[k]["J"] - pred[k]["EXT"]) / SCALE for k in X2_POINTS]
    say(f"   points {X2_POINTS}; J-EXT separations {', '.join(f'{s:.1f}' for s in seps)} scale")
    say(f"   J residuals {', '.join(f'{v:+.1f}' for v in rj2)} -> rms {rms_j2:.2f}")
    say(f"   EXT residuals {', '.join(f'{v:+.1f}' for v in re2)} -> rms {rms_e2:.2f}")
    winner = "EXT" if rms_e2 < rms_j2 else "J"
    ratio = max(rms_j2, rms_e2) / max(min(rms_j2, rms_e2), 1e-12)
    informative = ratio >= X2_SEP_BAR
    say(f"   winner {winner} by rms ratio {ratio:.2f} (bar {X2_SEP_BAR}); registered prediction: EXT")
    say(f"   -> X2 {'INFORMATIVE, ' + ('PASS' if winner == 'EXT' else 'FAIL') if informative else 'NOT INFORMATIVE'}")
    if informative and winner != "EXT":
        failures.append("X2")
    say()

    # ---------------------------------------------------------------- X3
    say("=== X3 -- the maximum keeps rising with ratio ===")
    m4 = {x: meas[(4.0, x)] for x in GRID4 if 0.0 < x < 1.0}
    m4[0.0] = mean_sd((4.0, 0.0))[0]
    x4max = max(m4, key=m4.get)
    r3max_x = max(r3, key=r3.get)
    d_max = (m4[x4max] - r3[r3max_x]) / SCALE
    say(f"   r4 max {m4[x4max]:.4f} at X_s {x4max:.2f}; r3 max {r3[r3max_x]:.4f} at "
        f"X_s {r3max_x:.2f}; difference {d_max:+.1f} scale (bar >= {X3_STEP_BAR})")
    steps_ok = True
    xs4 = [x for x in GRID4 if x < 1.0]
    for x0, x1 in zip(xs4[:-1], xs4[1:]):
        st = (m4[x1] - m4[x0]) / SCALE
        pj = (pred[(4.0, x1)]["J"] - pred[(4.0, x0)]["J"]) / SCALE
        # se of the DIFFERENCE of two seed means (lessons.md 2026-09-02): reported so a
        # step that misses the bar can be read as resolved or not; the bar is unchanged
        _, s0, n0 = mean_sd((4.0, x0))
        _, s1, n1 = mean_sd((4.0, x1))
        se = float(np.sqrt(s0 ** 2 / n0 + s1 ** 2 / n1)) / SCALE
        say(f"   step {x0:.2f}->{x1:.2f}: {st:+.1f} +- {se:.1f} scale (J predicted {pj:+.1f})")
        steps_ok &= st >= X3_STEP_BAR
    x3_ok = d_max >= X3_STEP_BAR and steps_ok
    say(f"   -> X3 {'PASS' if x3_ok else 'FAIL'}")
    if not x3_ok:
        failures.append("X3")
    say()

    # ---------------------------------------------------------------- X4
    say("=== X4 -- ratio ordering at fixed composition (model-free) ===")
    n_ok = 0
    for x in (0.05, 0.10, 0.15):
        a5, a4, a3 = meas[(5.0, x)], meas[(4.0, x)], r3[round(x, 6)]
        ok54 = a5 >= a4 - X4_SLACK * SCALE
        ok43 = a4 >= a3 - X4_SLACK * SCALE
        n_ok += int(ok54) + int(ok43)
        _, s5, n5 = mean_sd((5.0, x))
        _, s4, n4 = mean_sd((4.0, x))
        se54 = float(np.sqrt(s5 ** 2 / n5 + s4 ** 2 / n4)) / SCALE
        pj54 = (pred[(5.0, x)]["J"] - pred[(4.0, x)]["J"]) / SCALE
        say(f"   X_s {x:.2f}: r5 {a5:.4f}  r4 {a4:.4f}  r3 {a3:.4f}   "
            f"r5-r4 {(a5 - a4) / SCALE:+.1f} (se {se54:.1f}, J predicts {pj54:+.1f})  "
            f"r4-r3 {(a4 - a3) / SCALE:+.1f} scale  "
            f"{'ok' if ok54 else 'VIOLATED'}/{'ok' if ok43 else 'VIOLATED'}")
    x4_ok = n_ok == 6
    say(f"   -> X4 {'PASS' if x4_ok else 'FAIL'} ({n_ok} of 6)")
    if not x4_ok:
        failures.append("X4")
    say()

    # ---------------------------------------------------------------- X5
    say("=== X5 -- size-resolved structure (exploratory; k = 2 primary, all k printed) ===")
    say(f"   {'r':>4}{'X_s':>6}{'n_S':>6}{'n_L':>5}{'rat_S':>7}{'rat_L':>7}{'ratvol':>8}"
        f"{'pooled':>8}{'z_S':>6}{'z_L':>6}{'z_LL':>6}{'z_LS':>6}{'z_SL':>6}{'z_SS':>6}"
        f"{'sd zLL':>7}")
    bc_mean: dict[tuple[float, float], dict[str, float]] = {}
    for r, grid in ((4.0, GRID4), (5.0, GRID5)):
        for x in grid:
            key = (r, x)
            cs = adm[key]
            keys = ["n_S", "n_L", "rat_S_k2", "rat_L_k2", "ratvol_k2", "z_S", "z_L",
                    "z_LL", "z_LS", "z_SL", "z_SS"]
            d = {k: float(np.nanmean([c["by_class"][k] for c in cs])) for k in keys}
            d["pooled"] = float(np.mean([c["rattler"]["rattler_frac_k2"] for c in cs]))
            zll = np.array([c["by_class"]["z_LL"] for c in cs], dtype=float)
            d["sd_zLL"] = float(np.nanstd(zll, ddof=1)) if zll.size > 1 else float("nan")
            for k in RATTLER_THRESHOLDS:
                d[f"rat_S_k{k}"] = float(np.nanmean([c["by_class"][f"rat_S_k{k}"] for c in cs]))
                d[f"rat_L_k{k}"] = float(np.nanmean([c["by_class"][f"rat_L_k{k}"] for c in cs]))
            bc_mean[key] = d
            say(f"   {r:4.1f}{x:6.2f}{d['n_S']:6.0f}{d['n_L']:5.0f}{d['rat_S_k2']:7.3f}"
                f"{d['rat_L_k2']:7.3f}{d['ratvol_k2']:8.3f}{d['pooled']:8.3f}{d['z_S']:6.2f}"
                f"{d['z_L']:6.2f}{d['z_LL']:6.2f}{d['z_LS']:6.2f}{d['z_SL']:6.2f}{d['z_SS']:6.2f}"
                f"{d['sd_zLL']:7.2f}")
    say("   rat_S / rat_L at k = 3 and k = 4:")
    for key in sorted(bc_mean):
        if key[1] in (0.0, 1.0):
            continue
        d = bc_mean[key]
        say(f"   r{key[0]:.0f} X_s {key[1]:.2f}: k3 {d['rat_S_k3']:.3f}/{d['rat_L_k3']:.3f}  "
            f"k4 {d['rat_S_k4']:.3f}/{d['rat_L_k4']:.3f}")
    inner = [k for k in bc_mean if 0.0 < k[1] < 1.0]
    a_ok = all(bc_mean[k]["rat_S_k2"] > bc_mean[k]["rat_L_k2"] for k in inner)
    b_ok = all(bc_mean[k]["rat_L_k2"] <= X5_RAT_L_BAR for k in inner)
    z4 = [bc_mean[(4.0, x)]["z_LL"] for x in GRID4 if 0.0 < x < 1.0]
    c_mono = all(z4[i + 1] <= z4[i] for i in range(len(z4) - 1))
    sd_ref = bc_mean[(4.0, 0.05)]["sd_zLL"]
    c_drop = (z4[0] - z4[-1]) > sd_ref
    pooled4 = {x: bc_mean[(4.0, x)]["pooled"] for x in GRID4 if 0.0 < x < 1.0}
    x_pool = max(pooled4, key=pooled4.get)
    d_ok = x_pool < x4max
    say(f"   (a) rat_S > rat_L at every interior X_s, both ratios: {'PASS' if a_ok else 'FAIL'}")
    say(f"   (b) rat_L <= {X5_RAT_L_BAR} everywhere: {'PASS' if b_ok else 'FAIL'} "
        f"(max {max(bc_mean[k]['rat_L_k2'] for k in inner):.3f})")
    say(f"   (c) z_LL non-increasing at r4: {'PASS' if c_mono else 'FAIL'}; drop 0.05->0.30 = "
        f"{z4[0] - z4[-1]:+.2f} vs seed sd {sd_ref:.2f}: {'PASS' if c_drop else 'FAIL'}")
    say(f"   (d) pooled rattler_frac_k2 peaks at X_s {x_pool:.2f}, density max at {x4max:.2f}: "
        f"{'PASS' if d_ok else 'FAIL'}")
    say("   X5 is exploratory: recorded, not in FAILURES.")
    say()

    # ---------------------------------------------------------------- registered analysis
    say("=== registered analysis (section 5): direct per-ratio (a, b) at six u ===")
    allpts = list(pts) + [(r, x, meas[(r, x)]) for (r, x) in meas]
    say(f"   {'r':>5}{'u':>7}{'n':>3}{'a':>8}{'b':>8}{'chi2':>8}{'a[lo,hi]':>18}{'b[lo,hi]':>18}")
    for r in sorted({p[0] for p in allpts}):
        sub = [p for p in allpts if abs(p[0] - r) < 1e-9 and 0.0 < p[1] < 1.0]
        try:
            a_, b_, chi = mc.fit_ab_direct(sub, beta)
            alo, ahi = mc.profile_interval(sub, beta, 0)
            blo, bhi = mc.profile_interval(sub, beta, 1)
            # how many points the fit puts on the FINE branch: b enters the model only
            # there, so with zero fine-branch points b is bounded above and not measured
            nfine = sum(1 for _, x, _ in sub
                        if beta / (1 - (1 - a_) * x) > beta / (1 - (1 - beta) * (1 - b_) * (1 - x)))
            flags = []
            if not (alo - 1e-9 <= a_ <= ahi + 1e-9):
                flags.append("a: best fit OUTSIDE its own interval (cusp)")
            if not (blo - 1e-9 <= b_ <= bhi + 1e-9):
                flags.append("b: best fit OUTSIDE its own interval (cusp)")
            if nfine == 0:
                flags.append("b: 0 fine-branch points -> upper bound only")
            say(f"   {r:5.2f}{1 / r:7.3f}{len(sub):3d}{a_:8.3f}{b_:8.3f}{chi:8.1f}"
                f"   [{alo:.3f}, {ahi:.3f}]   [{blo:.3f}, {bhi:.3f}]   fine-branch pts {nfine}"
                + (("   ** " + "; ".join(flags)) if flags else ""))
        except Exception as exc:  # noqa: BLE001 - print, never hide
            say(f"   {r:5.2f}{1 / r:7.3f}{len(sub):3d}   fit refused: {exc!r}"[:150])
    say("   shipped forms at the same u:  a = sqrt(1-(1-u)^PA), b = 1-(1-u)^PB with J:")
    tab_points, fa, fb = tabulated_closure(beta, allpts)
    paJ, pbJ = cand["J"](0)
    for r in (1.25, 1.5, 2.0, 3.0, 4.0, 5.0):
        u = 1.0 / r
        say(f"   r{r:.2f}: a {np.sqrt(1 - (1 - u) ** paJ):.3f}  b {1 - (1 - u) ** pbJ:.3f}")
    say()

    # ---- registered analysis, second item: re-transfer of the TABULATED closure -------
    say("=== registered analysis (section 5): tabulated a(u), b(u) -> the three continuous families ===")
    us = [0.0, 1.0] + sorted(tab_points)
    aa = [0.0, 1.0] + [tab_points[u][0] for u in sorted(tab_points)]
    bb = [0.0, 1.0] + [tab_points[u][1] for u in sorted(tab_points)]
    order = np.argsort(us)
    aa = np.asarray(aa)[order]
    bb = np.asarray(bb)[order]
    say("   identifiable u (both intervals < 0.25 wide and containing the best fit): "
        + ", ".join(f"{u:.3f}" for u in sorted(tab_points)) + "; plus the proved limits u = 0, 1")
    say(f"   monotone in u:  a {bool(np.all(np.diff(aa) >= -1e-9))}   b {bool(np.all(np.diff(bb) >= -1e-9))}")
    fam_meas = mm.measured()
    pts23 = [p for p in pts if p[0] in (2.0, 3.0)]
    pa23, pb23, _ = mc.fit_on_binaries(pts23, beta)
    chi = {"textbook": 0.0, "refit23": 0.0, "J": 0.0, "tab": 0.0}
    npt = 0
    say(f"   {'family':<11}{'CV':>5}{'phi_meas':>10}{'textbook':>10}{'refit23':>9}{'J':>9}{'tab':>9}{'tab/scale':>10}")
    for fam in mm.FAMILIES:
        for cv in mc.CVS:
            if (fam, cv) not in fam_meas:
                continue
            d, y = mm.psd_classes(fam, cv)
            phi = fam_meas[(fam, cv)][0]
            v = {"textbook": mm.packing_density(d, y, beta, 1.02, 1.50),
                 "refit23": mm.packing_density(d, y, beta, pa23, pb23),
                 "J": mm.packing_density(d, y, beta, paJ, pbJ),
                 "tab": mc.packing_density_tabulated(d, y, beta, fa, fb)}
            for k in chi:
                chi[k] += ((v[k] - phi) / SCALE) ** 2
            npt += 1
            say(f"   {fam:<11}{cv:5.1f}{phi:10.4f}{v['textbook']:10.4f}{v['refit23']:9.4f}"
                f"{v['J']:9.4f}{v['tab']:9.4f}{(v['tab'] - phi) / SCALE:+10.1f}")
    say(f"   chi2 over {npt} held-out family points: textbook {chi['textbook']:.1f}, "
        f"two-ratio refit {chi['refit23']:.1f}, J {chi['J']:.1f}, tabulated {chi['tab']:.1f}")
    say("   (record, mixture_closure.txt: textbook 437.7, refit 63.0, tabulated-at-four-u 375.9)")
    say()

    say(f"FAILURES: {failures}")
    with open(os.path.join(OUT, "grade_mixture_ext.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
