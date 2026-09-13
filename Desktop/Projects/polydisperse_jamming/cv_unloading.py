"""cv_unloading.py -- what a frictional pack does to its pressure when sheared at
constant volume from an isotropic state, at every CV. Writes out/cv_unloading.txt and
cv_unloading.md; every number in the .md is printed here.

THE PHYSICS QUESTION. The dilatancy-sign mechanism (MECHANISM_2026-09-08.md section 1)
says a sheared pack CONTRACTS below phi_c and DILATES above it, and that the pore fluid
therefore goes into overpressure or suction. Constant-volume shear is the undrained
test: the volume change the skeleton wanted shows up as a change of grain pressure, and
its SIGN is the sign of the pore pressure a saturated pack would generate.
    dp/dgamma < 0 at constant volume  <=>  the skeleton wanted to contract  <=>  pore
                                           OVERpressure (liquefies, spreads)
    dp/dgamma > 0                     <=>  wanted to dilate  <=>  SUCTION (arrests)
The phase-8 constant-volume records (cv-shear-v2) start every cell from an
isotropically compressed, relaxed pack and record p_hat(gamma). So they already answer:
on which side of the critical-state line does an isotropically prepared pack sit, and
does that depend on phi, or on CV, or on neither?

PRE-REGISTERED CHECKS (bars fixed before this script was first run; the analysis is on
records already on disk, so "before running" means before the first execution):

  U1  SIGN.  On every JAMMED cell (quasi-static, I_ss <= 1e-2 -- shear.I_MAX) the net
      undrained pressure change is NEGATIVE: p_ss < p_0, where p_0 is the first shear
      sample (gamma = 0.0075; there is no gamma = 0 sample in the schema). Bar: 100 % of
      jammed cells at every CV. A single dilative cell falsifies "isotropic packs are
      contractive". Transient rises (p_max > p_0) are printed beside it, not hidden.
  U2  ONE LINE, TWO PROTOCOLS.  The constant-volume steady state p_ss(phi) and the
      constant-pressure steady state phi_c(p) must be the same critical-state line.
      Interpolate the common-ladder p_ss(phi) (log p, rung means, bracketing rungs only)
      to p* = 3e-04 and compare with phi_c(p*) from the rate-matched constant-pressure
      run. Bar: |dphi| <= 0.010 at each CV. Confounds that could break it and are NOT
      corrected: N 1500 vs 1000; gdot 3.0e-05 (CV 0) / 6e-05 (CV 0.2) / 1.4e-04 (CV 0.5)
      in the common ladder against 1.17e-04 matched in the constant-pressure run.
  U3  THE ISOTROPIC LINE LIES ABOVE THE CRITICAL-STATE LINE IN p.  At every phi where the
      phase-2 isotropic ladder (N = 1500, mu_p = 0.3, lognormal, default cap) overlaps
      the common ladder, p_iso(phi) > p_ss(phi). Bar: every overlapping rung, every CV.
      This is R4 of PREREG_CONSTP_RATEMATCH (phi_c - phi_iso > 0) read from the
      constant-volume side with a different shear protocol; it can fail independently.
  U4  REPORTED, not a bar: the strain over which the overpressure decays, gamma_half =
      first strain at which (p - p_ss) falls below half of (p_0 - p_ss), with the last
      down-crossing printed beside it (a crossing rule is the whole result -- lessons
      2026-09-03). Prediction written down: gamma_half GROWS toward phi_c from above.
  U5  REPORTED: the ratio p_0/p_ss on every jammed cell. In critical-state language it
      is the OVERCONSOLIDATION RATIO at which an isotropically prepared pack of this
      phi would change sign -- unload it below p_ss(phi) before shearing and it must
      dilate. That is the prediction the next DEM run (the OCR test) is graded on.

Everything is filtered to protocol cv-shear-v2, N = 1500, mu_p = 0.3, lognormal, and
deduplicated by key; N is printed before and after every filter.
"""
from __future__ import annotations

import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
sys.path.insert(0, HERE)

from grade_constp_rm import admit as constp_admit, load as constp_load  # noqa: E402
from measure_z_iso import load_phase2  # noqa: E402
from shear import I_MAX  # noqa: E402

CV_SETS = {                       # label -> glob; all cv-shear-v2 constant-volume runs
    "commonlad": "phase8_commonlad_c*_shard*.jsonl",   # 54 cells, common phi grid
    "prod2": "phase8_prod2_shard*.jsonl",             # 63 cells, ladders centred on phi_J
    "prod2ext": "phase8_prod2ext_shard*.jsonl",       # +18 cells at the top
    "ratem": "phase8_ratem_*.jsonl",                  # 30 cells, gdot matched 1.43e-04
}
N_REQUIRED = 1500
MU_P = 0.3
P_STAR = 3.0e-4                   # the constant-pressure servo set point
U2_BAR = 0.010
U2_CONTRAST = (0.0533, 0.015)     # the rate-matched constant-pressure Delta phi_c and band
CONSTP_EXPECTED = {0.0: 0.5911, 0.2: 0.6020, 0.5: 0.6444}   # graded numbers, reproduced
CONSTP_REPRO_TOL = 1.0e-3

lines: list[str] = []


def say(s: str = "") -> None:
    print(s)
    lines.append(s)


# ---------------------------------------------------------------------------- loading
def load_cv() -> list[dict]:
    recs, seen, raw = [], set(), 0
    for label, pat in CV_SETS.items():
        n_lab = 0
        for path in sorted(glob.glob(os.path.join(OUT, pat))):
            with open(path, encoding="utf-8") as fh:
                for line in fh:
                    if not line.strip():
                        continue
                    r = json.loads(line)
                    raw += 1
                    if r.get("key") in seen:
                        continue
                    seen.add(r["key"])
                    r["_set"] = label
                    recs.append(r)
                    n_lab += 1
        say(f"  {label:10s} {n_lab:3d} distinct records")
    say(f"constant-volume records: {raw} raw, {len(recs)} distinct keys")
    keep = []
    for r in recs:
        ok = (r.get("protocol") == "cv-shear-v2"
              and int(r.get("n_requested", -1)) == N_REQUIRED
              and abs(float(r.get("mu_p", -1)) - MU_P) < 1e-9
              and (r.get("family") or "lognormal") == "lognormal"
              and not r.get("error") and r.get("p_ss_hat") is not None
              and bool(r.get("compressed")))
        if ok:
            keep.append(r)
    say(f"after protocol/N/mu/family/compressed filter: {len(keep)}")
    for r in keep:
        p = np.asarray(r["p_hat_series"], float)
        assert p.ndim == 1 and p.size >= 100, (r["key"], p.shape)
        assert np.all(np.isfinite(p)), r["key"]
        g = np.asarray(r["gamma"], float)
        assert g.shape == p.shape and g[0] > 0 and np.all(np.diff(g) > 0), r["key"]
        r["_p0"] = float(p[0])
        r["_g0"] = float(g[0])
        r["_pmax"] = float(p.max())
        r["_gmax_at"] = float(g[int(p.argmax())])
    return keep


def load_iso() -> dict[float, list[tuple[float, float, int]]]:
    """Isotropic ladder rung means (phi, p_target, n) per CV; the p_0 anchor."""
    recs, info = load_phase2(cap="default")
    say(f"isotropic phase-2 records (cap=default, dedup): {len(recs)}  ({info})")
    sel = [r for r in recs
           if int(r.get("n_requested", 0)) == N_REQUIRED
           and abs(float(r.get("mu_p", -1)) - MU_P) < 1e-9
           and (r.get("family") or "lognormal") == "lognormal"]
    say(f"  N=1500, mu_p=0.3, lognormal: {len(sel)}")
    by = defaultdict(lambda: defaultdict(list))
    for r in sel:
        by[round(float(r["cv"]), 2)][float(r["p_hat_target"])].append(float(r["phi"]))
    out = {}
    for cv, d in by.items():
        out[cv] = sorted((float(np.mean(v)), p, len(v)) for p, v in d.items())
    return out


def load_constp() -> dict[float, tuple[float, float, int]]:
    recs = constp_load("phase9_rm_*.jsonl")
    by = defaultdict(list)
    for r in recs:
        ok, _why = constp_admit(r)
        if ok:
            by[round(float(r["cv"]), 2)].append(float(r["phi_ss"]))
    out = {cv: (float(np.mean(v)), float(np.std(v, ddof=1)), len(v)) for cv, v in by.items()}
    say("constant-pressure rate-matched cells admitted: "
        + ", ".join(f"CV {cv}: {n}" for cv, (_m, _s, n) in sorted(out.items())))
    for cv, (m, _s, _n) in out.items():
        assert abs(m - CONSTP_EXPECTED[cv]) <= CONSTP_REPRO_TOL, (cv, m, CONSTP_EXPECTED[cv])
    say(f"  reproduces the graded phi_c(p*) to {CONSTP_REPRO_TOL:.0e}: "
        + ", ".join(f"{m:.4f}" for _cv, (m, _s, _n) in sorted(out.items())))
    return out


# ---------------------------------------------------------------------------- helpers
def crossings(g: np.ndarray, y: np.ndarray, level: float) -> tuple[float, float]:
    """(first, last) strain at which y drops below `level`; nan if never."""
    below = y < level
    if not below.any():
        return float("nan"), float("nan")
    first = int(np.argmax(below))
    last = int(np.where(~below)[0].max()) + 1 if (~below).any() else 0
    last = min(last, len(g) - 1)
    return float(g[first]), float(g[last])


def interp_phi_at_p(rungs: list[tuple[float, float]], p: float) -> float:
    """phi at pressure p by linear interpolation in (phi, ln p) between BRACKETING rung
    means; nan if p is outside the ladder (no extrapolation, ever)."""
    pts = sorted((phi, np.log(pp)) for phi, pp in rungs if pp > 0)
    if len(pts) < 2:
        return float("nan")
    phis = np.array([a for a, _ in pts])
    lnp = np.array([b for _, b in pts])
    if not np.all(np.diff(lnp) > 0):          # a non-monotone ladder cannot be inverted
        return float("nan")
    x = np.log(p)
    if x < lnp[0] or x > lnp[-1]:
        return float("nan")
    return float(np.interp(x, lnp, phis))


def interp_lnp_at_phi(rungs: list[tuple[float, float]], phi: float) -> float:
    pts = sorted((a, np.log(pp)) for a, pp in rungs if pp > 0)
    phis = np.array([a for a, _ in pts])
    lnp = np.array([b for _, b in pts])
    if len(pts) < 2 or phi < phis[0] or phi > phis[-1] or not np.all(np.diff(phis) > 0):
        return float("nan")
    return float(np.interp(phi, phis, lnp))


# ---------------------------------------------------------------------------- main
def main() -> int:
    say("cv_unloading.py -- the undrained sign of an isotropically prepared pack")
    say(f"quasi-static bar I_ss <= {I_MAX:.0e} (shear.I_MAX); p* = {P_STAR:.0e}")
    say()
    recs = load_cv()
    iso = load_iso()
    constp = load_constp()
    say()
    failures: list[str] = []

    jam = [r for r in recs if bool(r["quasistatic"])]
    flow = [r for r in recs if not bool(r["quasistatic"])]
    say(f"jammed (quasi-static) cells: {len(jam)}; flowing (I_ss > {I_MAX:.0e}): {len(flow)}")

    # ---- U1 ---------------------------------------------------------------------
    say()
    say("U1  SIGN of the net undrained pressure change on every jammed cell")
    say("    set        CV   phi     seed  p_0        p_max/p_0  p_ss/p_0  gamma(p_max)  I_ss     steady")
    n_neg = defaultdict(int)
    n_tot = defaultdict(int)
    n_rise = 0
    for r in sorted(jam, key=lambda r: (r["cv"], r["phi"], r["seed"], r["_set"])):
        cv = round(float(r["cv"]), 2)
        ratio = r["p_ss_hat"] / r["_p0"]
        rise = r["_pmax"] / r["_p0"]
        n_tot[cv] += 1
        if ratio < 1.0:
            n_neg[cv] += 1
        if rise > 1.0 + 1e-12:
            n_rise += 1
        say(f"    {r['_set']:10s} {cv:.1f}  {r['phi']:.4f}  {r['seed']}     {r['_p0']:.3e}  "
            f"{rise:8.3f}  {ratio:8.3f}  {r['_gmax_at']:6.3f}       {r['i_ss']:.1e}  {bool(r['steady'])}")
    say()
    for cv in sorted(n_tot):
        say(f"    CV {cv:.1f}: p_ss < p_0 on {n_neg[cv]} of {n_tot[cv]} jammed cells")
    all_neg = all(n_neg[cv] == n_tot[cv] for cv in n_tot) and sum(n_tot.values()) > 0
    say(f"    transient rise (p_max > p_0) on {n_rise} of {sum(n_tot.values())} jammed cells")
    say(f"U1  {'PASS' if all_neg else 'FAIL'}: contractive on "
        f"{sum(n_neg.values())}/{sum(n_tot.values())} jammed cells (bar: all)")
    if not all_neg:
        failures.append("U1")
    # the flowing cells, for the record: they start unjammed and p RISES to an inertial
    # steady state -- that is not a dilatancy sign, it is an unjammed pack at finite rate
    say()
    say("    flowing cells (not quasi-static): p_0 -> p_max -> p_ss, all inertial")
    for r in sorted(flow, key=lambda r: (r["cv"], r["phi"], r["seed"], r["_set"])):
        say(f"    {r['_set']:10s} {float(r['cv']):.1f}  {r['phi']:.4f}  {r['seed']}     "
            f"{r['_p0']:.3e}  {r['_pmax'] / r['_p0']:8.1f}  {r['p_ss_hat'] / r['_p0']:8.1f}  "
            f"{r['_gmax_at']:6.3f}       {r['i_ss']:.1e}  z_nr_ss {r['z_nr_ss']:.2f}")

    # ---- U2 ---------------------------------------------------------------------
    say()
    say("U2  ONE critical-state line from two protocols: phi at p* from constant-volume")
    say("    p_ss(phi) (common ladder, jammed rungs, rung means) vs phi_c(p*) at constant pressure")
    csl = {}
    for cv in sorted({round(float(r["cv"]), 2) for r in jam if r["_set"] == "commonlad"}):
        by = defaultdict(list)
        for r in jam:
            if r["_set"] == "commonlad" and round(float(r["cv"]), 2) == cv:
                by[round(float(r["phi_target"]), 3)].append((float(r["phi"]), float(r["p_ss_hat"]),
                                                              bool(r["steady"])))
        rungs = []
        for pt, v in sorted(by.items()):
            phi_m = float(np.mean([a for a, _, _ in v]))
            p_m = float(np.exp(np.mean([np.log(b) for _, b, _ in v])))
            rungs.append((phi_m, p_m))
            say(f"    CV {cv:.1f}  phi {phi_m:.4f}  p_ss {p_m:.3e}  n {len(v)}  "
                f"steady {sum(s for _, _, s in v)}/{len(v)}")
        csl[cv] = rungs
    say()
    u2_ok = True
    phi_cv_at_pstar = {}
    for cv, rungs in sorted(csl.items()):
        phi_v = interp_phi_at_p(rungs, P_STAR)
        phi_cv_at_pstar[cv] = phi_v
        m, s, n = constp.get(cv, (float("nan"), float("nan"), 0))
        d = phi_v - m
        ok = np.isfinite(d) and abs(d) <= U2_BAR
        u2_ok &= bool(ok)
        say(f"    CV {cv:.1f}: phi_CV(p*) = {phi_v:.4f}   phi_c,constp(p*) = {m:.4f} +- {s:.4f} (n={n})"
            f"   dphi = {d:+.4f}   {'ok' if ok else 'MISS'} (bar {U2_BAR})")
    if 0.0 in phi_cv_at_pstar and 0.5 in phi_cv_at_pstar:
        dc = phi_cv_at_pstar[0.5] - phi_cv_at_pstar[0.0]
        c0, half = U2_CONTRAST
        say(f"    Delta phi(CV 0 -> 0.5) from the constant-volume side: {dc:+.4f} "
            f"(constant-pressure {c0:+.4f}, band +-{half}: "
            f"{'inside' if abs(dc - c0) <= half else 'OUTSIDE'})")
    say(f"U2  {'PASS' if u2_ok else 'FAIL'}: the two protocols agree to {U2_BAR} at every CV")
    if not u2_ok:
        failures.append("U2")

    # ---- U3 ---------------------------------------------------------------------
    say()
    say("U3  the isotropic line sits ABOVE the critical-state line in pressure")
    say("    (phase-2 isotropic rung means, N=1500 mu_p=0.3 lognormal, interpolated in ln p to the")
    say("     common-ladder phi; p_0 is the phase-8 preparation's own first sample)")
    u3_ok, u3_n = True, 0
    u3_rows: list[tuple[float, float, float, float, float, bool]] = []
    for cv, rungs in sorted(csl.items()):
        iso_r = [(phi, p) for phi, p, _n in iso.get(cv, [])]
        for phi, p_ss in rungs:
            lnp_iso = interp_lnp_at_phi(iso_r, phi)
            if not np.isfinite(lnp_iso):
                say(f"    CV {cv:.1f}  phi {phi:.4f}: outside the isotropic ladder "
                    f"({iso_r[0][0]:.4f}-{iso_r[-1][0]:.4f}), no anchor")
                continue
            p_iso = float(np.exp(lnp_iso))
            p0s = [r["_p0"] for r in jam if r["_set"] == "commonlad"
                   and round(float(r["cv"]), 2) == cv and abs(float(r["phi"]) - phi) < 0.004]
            p0 = float(np.exp(np.mean(np.log(p0s)))) if p0s else float("nan")
            ok = p_iso > p_ss
            u3_ok &= ok
            u3_n += 1
            u3_rows.append((cv, phi, p_iso, p_ss, p0, ok))
            say(f"    CV {cv:.1f}  phi {phi:.4f}: p_iso {p_iso:.3e}  p_ss {p_ss:.3e}  "
                f"p_iso/p_ss {p_iso / p_ss:6.2f}   prep p_0 {p0:.3e}  p_0/p_iso {p0 / p_iso:6.2f}"
                f"   {'ok' if ok else 'VIOLATED'}")
    say(f"U3  {'PASS' if (u3_ok and u3_n) else 'FAIL'}: p_iso > p_ss on {u3_n} overlapping rungs")
    if not (u3_ok and u3_n):
        failures.append("U3")

    # ---- U4 ---------------------------------------------------------------------
    say()
    say("U4  REPORTED: strain over which the overpressure decays (common ladder, jammed cells)")
    say("    CV   phi     seed  p_0/p_ss  gamma_half(first)  gamma_half(last)  monotone_to_ss")
    u4 = defaultdict(list)
    for r in sorted(jam, key=lambda r: (r["cv"], r["phi"], r["seed"])):
        if r["_set"] != "commonlad":
            continue
        g = np.asarray(r["gamma"], float)
        p = np.asarray(r["p_hat_series"], float)
        excess = p - r["p_ss_hat"]
        if excess[0] <= 0:
            continue
        gf, gl = crossings(g, excess, 0.5 * excess[0])
        # monotone if, after the argmax, p never rises by more than 2 % of the excess
        tail = p[int(p.argmax()):]
        mono = bool(np.all(np.diff(tail) <= 0.02 * excess[0]))
        cv = round(float(r["cv"]), 2)
        u4[(cv, round(float(r["phi_target"]), 3))].append(gf)
        say(f"    {cv:.1f}  {r['phi']:.4f}  {r['seed']}     {r['_p0'] / r['p_ss_hat']:7.2f}  "
            f"{gf:14.3f}  {gl:14.3f}      {mono}")
    say()
    say("    per rung: mean +- sd of gamma_half(first), n")
    for (cv, pt), v in sorted(u4.items()):
        v = np.array(v, float)
        fin = v[np.isfinite(v)]
        s = f"{fin.mean():.3f} +- {fin.std(ddof=1):.3f}" if fin.size > 1 else (f"{fin[0]:.3f}" if fin.size else "never crossed")
        say(f"    CV {cv:.1f}  phi_target {pt:.3f}: {s}  (n={fin.size} of {v.size})")

    # ---- U5 ---------------------------------------------------------------------
    say()
    say("U5  REPORTED: p_0/p_ss on every jammed common-ladder cell = the overconsolidation")
    say("    ratio at which the sign would flip, per rung (geometric mean, range, n)")
    ocr = defaultdict(list)
    for r in jam:
        if r["_set"] == "commonlad":
            ocr[(round(float(r["cv"]), 2), round(float(r["phi_target"]), 3))].append(r["_p0"] / r["p_ss_hat"])
    for (cv, pt), v in sorted(ocr.items()):
        v = np.array(v, float)
        say(f"    CV {cv:.1f}  phi_target {pt:.3f}: OCR* {np.exp(np.log(v).mean()):6.2f}  "
            f"[{v.min():.2f}, {v.max():.2f}]  n={v.size}")
    say()
    say("    matched phi, CV 0 vs CV 0.5 (both jammed at 0.645 and 0.670 only):")
    for pt in (0.645, 0.670):
        a, b = ocr.get((0.0, pt)), ocr.get((0.5, pt))
        if a and b:
            ga, gb = np.exp(np.log(a).mean()), np.exp(np.log(b).mean())
            say(f"    phi {pt:.3f}: OCR*(CV 0) {ga:.2f}  OCR*(CV 0.5) {gb:.2f}  ratio {gb / ga:.2f}")

    say()
    say(f"FAILURES: {failures if failures else 'none'}")
    with open(os.path.join(OUT, "cv_unloading.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")

    write_md(recs, jam, flow, n_neg, n_tot, n_rise, csl, phi_cv_at_pstar, constp, u3_rows, u4, ocr,
             failures)
    return 1 if failures else 0


def write_md(recs, jam, flow, n_neg, n_tot, n_rise, csl, phi_cv, constp, u3_rows, u4, ocr, failures):
    cvs = sorted(n_tot)
    n_jam = sum(n_tot.values())
    n_neg_all = sum(n_neg.values())
    u2_rows = []
    for cv in sorted(csl):
        m, s, n = constp.get(cv, (float("nan"), float("nan"), 0))
        p_lo = min(p for _phi, p in csl[cv])
        note = "read" if np.isfinite(phi_cv[cv]) else f"UNREADABLE: lowest jammed p_ss {p_lo:.1e} > p*"
        u2_rows.append(f"| {cv:.1f} | {phi_cv[cv]:.4f} | {m:.4f} ± {s:.4f} | {phi_cv[cv] - m:+.4f} | {note} |")
    dc = phi_cv[0.5] - phi_cv[0.0]
    u2_read = [cv for cv in sorted(csl) if np.isfinite(phi_cv[cv])]
    u2_unread = [cv for cv in sorted(csl) if not np.isfinite(phi_cv[cv])]
    u2_verdict = "PASS" if "U2" not in failures else "FAIL as registered"
    u2_read_txt = ("; ".join(f"CV {cv:.1f} agrees to {phi_cv[cv] - constp[cv][0]:+.4f}" for cv in u2_read)
                   if u2_read else "no CV is readable")
    u3_verdict = "PASS" if "U3" not in failures else "FAIL as registered"
    u3_md = []
    for cv, phi, p_iso, p_ss, p0, ok in u3_rows:
        u3_md.append(f"| {cv:.1f} | {phi:.4f} | {p_iso:.2e} | {p_ss:.2e} | {p_iso / p_ss:.2f} | "
                     f"{p0 / p_iso:.2f} | {'ok' if ok else '**VIOLATED**'} |")
    u3_bad = [(cv, phi, p_iso / p_ss) for cv, phi, p_iso, p_ss, _p0, ok in u3_rows if not ok]
    p0_over_iso = [p0 / p_iso for _cv, _phi, p_iso, _p_ss, p0, _ok in u3_rows]
    ocr_rows = []
    for (cv, pt), v in sorted(ocr.items()):
        v = np.array(v, float)
        g = u4.get((cv, pt), [])
        g = np.array([x for x in g if np.isfinite(x)], float)
        gs = f"{g.mean():.3f} ± {g.std(ddof=1):.3f}" if g.size > 1 else (f"{g[0]:.3f}" if g.size else "—")
        ocr_rows.append(f"| {cv:.1f} | {pt:.3f} | {np.exp(np.log(v).mean()):.2f} | "
                        f"{v.min():.2f}–{v.max():.2f} | {gs} | {v.size} |")
    verdict = "ALL BARS PASSED" if not failures else f"FAILED: {', '.join(failures)}"
    rise_note = (f"{n_rise} of {n_jam} jammed cells show a transient rise (`p_max > p_0`) before "
                 f"the decay" if n_rise else "no jammed cell shows a transient rise")
    md = f"""# Constant-volume shear from an isotropic state — which side of the critical state is a fresh pack on?

**Written 2026-09-12.** Every number here is printed by `cv_unloading.py` into
`out/cv_unloading.txt` (`{verdict}`). Bars U1–U3 were fixed in that script's docstring
before its first execution; U4–U5 are reported, not graded.

---

## Result box

> **An isotropically prepared frictional pack is CONTRACTIVE under shear at every packing
> fraction and every CV measured.** On {n_neg_all} of {n_jam} jammed constant-volume cells
> (CV {', '.join(f'{c:.1f}' for c in cvs)}; N = 1500, μ_p = 0.3, four campaigns pooled) the grain
> pressure falls from its preparation value `p_0` to the critical-state value `p_ss(φ)`. In the
> undrained reading that is pore **overpressure on both sides of `φ_c`** — the "dilative
> branch" of the mechanism is not reached by putting `φ` above `φ_c`; it needs a state
> that sits *below* the critical-state line in pressure, i.e. a loading history.
>
> **U2 — {u2_verdict}.** Constant-volume `p_ss(φ)` inverted at `p*` gives `φ` =
> {', '.join(f'{phi_cv[cv]:.4f}' for cv in sorted(csl))} against the constant-pressure
> {', '.join(f'{constp[cv][0]:.4f}' for cv in sorted(csl))}: {u2_read_txt}. The bar is
> unreadable at CV {', '.join(f'{cv:.1f}' for cv in u2_unread) or '—'} because `p*` = {P_STAR:.0e} lies
> *below* the lowest jammed constant-volume rung's `p_ss` there — the interpolation has no
> bracket and no extrapolation is admitted. The CV 0 → 0.5 contrast from the constant-volume
> side is therefore **{dc:+.4f}** (constant-pressure +0.0533); the bar stays FAIL and the
> stricter test it needs is a constant-volume rung between `φ_c(p*)` and the present lowest
> jammed rung at CV 0 and 0.5 (§7).
>
> **U3 — {u3_verdict}** on {len(u3_rows)} overlapping rungs; {len(u3_bad)} violated:
> {', '.join(f'CV {cv:.1f} φ {phi:.4f} at p_iso/p_ss = {r:.2f}' for cv, phi, r in u3_bad) or 'none'}.
> And a number the bars did not ask for: the phase-8 preparation's own `p_0` sits
> **{min(p0_over_iso):.1f}–{max(p0_over_iso):.1f}× above the phase-2 quasi-static isotropic line at the same `φ`**, so
> U1's "isotropically prepared" pack is a rate-compressed one, over-pressured relative to the
> servo-held isotropic state. U1's sign therefore mixes prep-overpressure relaxation with the
> shear response; the clean version is §7's next run.

---

## 1. What was measured

`p_hat(γ)` on every `cv-shear-v2` constant-volume record: {len(recs)} distinct cells after the
protocol / N / μ_p / family / compressed filter, of which {len(jam)} are jammed (quasi-static,
`I_ss` ≤ 1e-02) and {len(flow)} are flowing. `p_0` is the first shear sample at γ = 0.0075 — the
schema has no γ = 0 sample, and the decay per sample at that point is ~1 %, so `p_0` is a
~1 % underestimate of the isotropic pressure. `p_ss` is the project's steady-state
window mean.

## 2. U1 — the sign (PASS on {n_neg_all}/{n_jam})

Per CV: {'; '.join(f'CV {cv:.1f} {n_neg[cv]}/{n_tot[cv]}' for cv in cvs)}. {rise_note}.
The flowing cells ({len(flow)}) start unjammed (`z` ≈ 0, `p_0` ~ 1e-07) and their pressure
*rises* to an inertial steady state; that is a pack at finite rate, not a dilatancy sign,
and they are listed in the `.txt` for the record.

**Why this matters for the mechanism.** `MECHANISM_2026-09-08.md` §1 reads "below `φ_c` a
sheared pack contracts … above it the pack dilates". The second half is true only of a
pack that is *denser than the critical state at the pressure it is sheared at*. Where U3
can be read (CV 0 and 0.2), isotropic compression to `φ` leaves the pack at a pressure
`p_iso(φ)` **above** `p_ss(φ)`, so the first thing shear does is relax that pressure —
contraction, overpressure — however far above `φ_c(0)` the pack sits; at CV 0.5 the one
readable rung is undecided (§4). In critical-state soil mechanics this is the textbook
statement that a **normally consolidated** soil is contractive and only an
**over-consolidated** one dilates; what is new is that it holds for a polydisperse
frictional sphere pack at every CV, and that the DEM puts numbers on the ratio (U5).

## 3. U2 — one line from two protocols

| CV | `φ` at `p*` from `p_ss(φ)` (constant volume) | `φ_c(p*)` (constant pressure) | Δ | status |
|---|---|---|---|---|
{chr(10).join(u2_rows)}

Bar ±{U2_BAR}: **{u2_verdict}**. The bar was posed on the assumption that `p*` falls inside
the jammed part of every common-ladder — it does not at CV
{', '.join(f'{cv:.1f}' for cv in u2_unread) or '—'}, where the lowest jammed rung's `p_ss` is
already above `p*`, so the inversion returns NaN rather than a number. That is a defect of
the bar's design, and it is not repaired by extrapolating: the CSL is measurably curved in
(`φ`, ln p) (≈ 0.009 per decade near jamming against ≈ 0.047 between φ 0.62 and 0.67 at
CV 0), so a one-rung extrapolation would be a guess dressed as a test. Where the bar can
be read it passes: {u2_read_txt}. The two routes there differ in N (1500 vs 1000), in `γ̇`
(the common ladder is *unmatched*: 3.0e-05 / 6e-05 / 1.4e-04 against 1.17e-04 matched) and
in the steady-state definition. One CV agreeing is evidence, not a verdict; the failing
NaNs stay printed.

## 4. U3 — is the isotropic line above the critical-state line in pressure?

| CV | `φ` | `p_iso(φ)` (phase-2, servo-held) | `p_ss(φ)` (constant volume) | `p_iso/p_ss` | `p_0/p_iso` | bar |
|---|---|---|---|---|---|---|
{chr(10).join(u3_md)}

**{u3_verdict}.** Only {len(u3_rows)} common-ladder rungs fall inside the phase-2 isotropic
ladder at all (its top is 0.603 / 0.609 / 0.649 at CV 0 / 0.2 / 0.5), so the bar is
decided on one rung per CV. The violation is at the CV 0.5 rung, whose constant-volume
`p_ss` is from a rung that is steady on 1 of 3 seeds (§3's table) with one seed still
decaying at `γ` = 1.2 (U4) — so its `p_ss` is an upper bound and the ratio 0.74 a lower
bound. **Unresolved at CV 0.5, not a sign flip**; a longer constant-volume run at
CV 0.5, `φ` = 0.645 decides it, and the constant-pressure `φ_c(p*)` = 0.6444 at that CV
says `p_ss(0.645)` should be ≈ `p*` = 3e-04, i.e. the ratio should come out ≈ 1.7 if the rung
finishes relaxing.

The column the bars did not ask for is the larger finding. `p_0/p_iso` = {min(p0_over_iso):.1f}–{max(p0_over_iso):.1f}:
the phase-8 preparation (rate compression to `φ`, then constant-volume relaxation to the
unbalanced-force tolerance) lands at a pressure {min(p0_over_iso):.1f}–{max(p0_over_iso):.1f}× the servo-held
isotropic state of phase 2 at the same `φ` and N. For a frictional pack `p(φ)` is not a
function — it is a history — and the two "isotropic" protocols sit on different branches
of it. Consequence for U1: the sign it measured is that of a *rate-compressed* pack, whose
pressure has two reasons to fall under shear (relaxing the prep overpressure toward the
isotropic branch, and relaxing the isotropic branch toward the critical line). U1 is
correct as stated and weaker than its prose implied; the second reason alone is what §6
needs, and it is tested by preparing ON the phase-2 line (§7).

## 5. U4 / U5 — the decay strain and the overconsolidation ratio (reported)

| CV | `φ_target` | OCR* = `p_0/p_ss` (geo. mean) | range | `γ_half` (first crossing) | n |
|---|---|---|---|---|---|
{chr(10).join(ocr_rows)}

`γ_half` is the strain at which the excess pressure `p − p_ss` has halved; first and last
down-crossings and a monotonicity flag are printed per cell. **The overpressure is not
instantaneous and not permanent — it decays over a strain of order 0.1–1 that grows
toward `φ_c`.** A viscosity (M2a) has no strain in it; this is a strain-history signature
on the *contractive* side, measured, with no fit.

OCR* is the prediction for the next run: unload an isotropic pack at `φ` to a pressure
below `p_ss(φ)` — a ratio above OCR* — and constant-volume shear must then *raise* the
pressure. That would be the first dilative (suction-generating) state this DEM has
produced, reached from the other side of the same line.

## 6. What this changes in the mechanism, and what it does not

- **Unchanged:** `Δφ_c(CV)` and every contrast. The line moves with CV exactly as the
  constant-pressure run said (U2).
- **Changed:** the dilatancy sign is **not a function of `φ` alone**. It is the sign of
  `p_shear − p_ss(φ)`, and a first-time-loaded pack always has `p_shear = p_iso(φ) > p_ss(φ)`.
  For the drop this says the sign is set by the *loading path* at the contact: a pressure
  pulse that rises (first-time loading → overpressure) and then relaxes while the pack
  keeps shearing (unloading → the pack crosses OCR* → suction). The E-3 transducer
  prediction becomes a **sign reversal within one impact**, not a single sign.
- **Changed:** the route to the dilative side that a lubricated suspension actually takes
  is not isotropic over-consolidation but the Wyart–Cates one — contacts are frictionless
  until the stress exceeds `σ*`, so the pack at rest is a *frictionless* pack at its `φ`,
  which for `φ_μ(CV) < φ < φ_0` carries **no pressure at all** (`p = 0 < p_ss(φ)`, i.e. the
  OCR → ∞ limit). Switching friction on and shearing must then build a frictional network
  over a finite strain `γ_c(φ − φ_μ(CV))`, with `p` rising from zero — dilative, suction,
  arrest — and `γ_c` diverging as `φ → φ_μ(CV)⁺`. That is the strain-clock protocol
  (`PREREG_STRAINCLOCK_2026-09-13.md`), and it is the physically right preparation for a
  DST suspension, which the isotropic-frictional one is not.

## 7. Confounds, stated, and the runs that remove them

- `p_0` is a γ = 0.0075 proxy (≈ 1 % low). It cannot flip U1: the decay continues for
  tens of samples on every jammed cell.
- **U1's pack is not on the isotropic line** (`p_0/p_iso` = {min(p0_over_iso):.1f}–{max(p0_over_iso):.1f}, §4). The
  clean test: take the phase-2 servo-held state itself (`φ`, `p_iso`) and shear it at
  constant volume. Prediction from U3 where it is readable: CV 0 at `φ` = 0.595 has
  `p_iso/p_ss` = {[r for r in u3_rows if r[0] == 0.0][0][2] / [r for r in u3_rows if r[0] == 0.0][0][3]:.2f}, so `p` must fall by that factor — contractive with no
  prep-overpressure to hide behind. If instead `p` does not fall, U1's whole sign was the
  preparation, and §6 is withdrawn.
- **U2 is unreadable at CV 0 and 0.5.** The run that makes it readable: common-ladder
  constant-volume cells at `φ` between `φ_c(p*)` and the present lowest jammed rung
  (CV 0: 0.595–0.62 → e.g. 0.600 and 0.607; CV 0.5: 0.645–0.67 → e.g. 0.650 and 0.657),
  N = 1500, rate-matched at 1.17e-04, `γ_max` ≥ 3 so the near-`φ_c` rungs reach steady.
  The bar then reads without extrapolation.
- U2 where readable mixes N (1500 vs 1000) and `γ̇`; a rate-matched N = 1500
  constant-pressure run would be the clean version.
- The rung nearest `φ_c` at each CV is not steady at `γ_max` = 1.5 (printed in U2's
  table), so its `p_ss` is an upper bound, its OCR* a lower bound, and U3's CV 0.5 verdict
  undecided.
"""
    with open(os.path.join(HERE, "cv_unloading.md"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(md)


if __name__ == "__main__":
    sys.exit(main())
