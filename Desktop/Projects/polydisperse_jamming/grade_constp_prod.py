"""grade_constp_prod.py — grade PREREG_CONSTP_PROD_2026-09-12.md, P1-P5.

Written before the production cells landed. Every bar is quoted from that file's §3-§4.

    python grade_constp_prod.py          -> out/constp_prod.txt

Exit 0 only if P1 is RESOLVED and passes.
"""

from __future__ import annotations

import glob
import json
import os
from collections import defaultdict

import numpy as np

from grade_constp import stability

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

PATTERN = "phase9_prod*.jsonl"
CVS = (0.0, 0.2, 0.5)
STARTS = (0.55, 0.64)

# --- bars, QUOTED from PREREG_CONSTP_PROD_2026-09-12.md ----------------------------
ADMIT_Z_MIN = 2.0            # §3.1
ADMIT_P_FLOOR_FRAC = 0.1     # §3.1
ADMIT_P_TOL = 0.05           # §3.2
MIN_SEEDS = 2                # §3.2 / §5: fewer than this and the block is UNRESOLVED
P1_REF = +0.0525             # the common ladder's midpoint
P1_BAND = 0.015              # §4 P1, the same band PREREG_COMMON_LADDER §3 used
P2_SHAKEOUT = 0.5935         # phi_ss(CV 0) at N = 600
P2_BAR = 0.005
P4_BAR = 0.005
# isotropic phi at p = 3e-04, mu_p = 0.3, N = 1500, from out/phase2_summary.csv via
# measure_z_iso.load_phase2; recomputed here rather than trusted (see below).
PHI_ISO_FALLBACK = {0.0: 0.583314}


def phi_iso_at_p(p_target: float = 3.0e-4) -> dict:
    """Isotropic phi at the servo's set pressure, per CV, from the project's own records.

    Recomputed rather than hardcoded so that P3 cannot drift away from the ladder it
    compares against. Uses the DEFAULT-cap selection, which is what every published
    number rests on (PREREG_ZISO_2026-09-12.md §8).
    """
    try:
        from measure_z_iso import load_phase2
    except Exception:
        return dict(PHI_ISO_FALLBACK)
    recs, _info = load_phase2(cap="default")
    by = defaultdict(list)
    for r in recs:
        if int(r.get("n_requested", 0)) != 1500:
            continue
        if abs(float(r.get("mu_p", -1)) - 0.3) > 1e-9:
            continue
        if abs(float(r.get("p_hat_target", -1)) - p_target) > 1e-12:
            continue
        if str(r.get("family", "lognormal")) != "lognormal":
            continue
        by[round(float(r["cv"]), 3)].append(float(r["phi"]))
    return {cv: float(np.mean(v)) for cv, v in by.items() if v}


def load() -> list[dict]:
    recs, seen = [], set()
    for path in sorted(glob.glob(os.path.join(OUT, PATTERN))):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                r = json.loads(line)
                if r.get("key") in seen:
                    continue
                seen.add(r["key"])
                recs.append(r)
    return recs


def admit(r: dict) -> tuple[bool, str]:
    if r.get("error"):
        return False, "errored"
    d = stability(r)
    if not (np.isfinite(d["z_win_min"]) and d["z_win_min"] >= ADMIT_Z_MIN):
        return False, f"UNJAMMED (z_win_min {d['z_win_min']:.2f})"
    if d["p_win_min"] < ADMIT_P_FLOOR_FRAC * float(r["p_hat_target"]):
        return False, f"p visited {d['p_win_min']:.1e}, below 0.1 p*"
    if abs(r["p_ss_hat"] / r["p_hat_target"] - 1) > ADMIT_P_TOL:
        return False, f"servo off by {r['p_ss_hat'] / r['p_hat_target'] - 1:+.1%}"
    if not r.get("steady"):
        return False, (f"not steady (p {r['steady_split']:.1%}, "
                       f"phi {r['phi_ss_split']:.1%})")
    if not r.get("quasistatic"):
        return False, f"not quasi-static (I {r['i_ss']:.1e})"
    return True, ""


def main() -> int:
    lines: list[str] = []

    def p(m: str = "") -> None:
        print(m)
        lines.append(m)

    recs = load()
    p("=== grade_constp_prod.py -- PREREG_CONSTP_PROD_2026-09-12.md, P1-P5 ===")
    p(f"   {len(recs)} records on disk (18 expected)")
    if not recs:
        p("   🔴 no records -- nothing to grade")
        return 1

    ok_cells, refused = [], []
    for r in recs:
        good, why = admit(r)
        (ok_cells if good else refused).append((r, why))

    p("")
    p(f"{'cv':>5s}{'phi0':>6s}{'seed':>5s}{'phi_ss':>10s}{'p_err':>9s}{'z_wmin':>8s}"
      f"{'alpha':>9s}{'dec':>6s}{'adm':>5s}  why refused")
    for r, why in sorted(recs and [(r, w) for r, w in ok_cells + refused],
                         key=lambda t: (t[0].get("cv", 0), t[0].get("phi_target", 0),
                                        t[0].get("seed", 0))):
        if r.get("error"):
            p(f"{r.get('cv', -1):5.2f}{r.get('phi_target', -1):6.2f}"
              f"{r.get('seed', -1):5d}{'':>10s}{'':>9s}{'':>8s}{'':>9s}{'':>6s}"
              f"{'no':>5s}  {r['error'][:60]}")
            continue
        d = stability(r)
        adm = "yes" if not why else "no"
        p(f"{r['cv']:5.2f}{r['phi_target']:6.2f}{r['seed']:5d}{r['phi_ss']:10.6f}"
          f"{r['p_ss_hat'] / r['p_hat_target'] - 1:+9.1%}{d['z_win_min']:8.2f}"
          f"{r['alpha']:9.3g}{r['alpha_decades']:6.2f}{adm:>5s}  {why}")
    p("")
    p(f"   admitted {len(ok_cells)} of {len(recs)}; refused {len(refused)}")

    # per (cv, start) blocks
    blocks: dict[tuple, list[dict]] = defaultdict(list)
    for r, _w in ok_cells:
        blocks[(round(float(r["cv"]), 3), round(float(r["phi_target"]), 3))].append(r)

    p("")
    p(f"=== the blocks (>= {MIN_SEEDS} admitted seeds to be usable) ===")
    means: dict[tuple, float] = {}
    for k in sorted(blocks):
        v = np.array([r["phi_ss"] for r in blocks[k]], float)
        usable = v.size >= MIN_SEEDS
        sd = float(v.std(ddof=1)) if v.size > 1 else float("nan")
        p(f"   cv {k[0]:.1f} phi0 {k[1]:.2f}: n = {v.size}, phi_ss = {v.mean():.6f}"
          f" +- {sd:.6f}   {'usable' if usable else 'UNRESOLVED (too few seeds)'}")
        if usable:
            means[k] = float(v.mean())
    for cv in CVS:
        for s in STARTS:
            if (round(cv, 3), round(s, 3)) not in blocks:
                p(f"   cv {cv:.1f} phi0 {s:.2f}: NO admitted cell")

    # per-CV mean, pooling the two preparations (legitimate only if P4 holds)
    per_cv: dict[float, float] = {}
    for cv in CVS:
        vals = [means[k] for k in means if abs(k[0] - cv) < 1e-9]
        if vals:
            per_cv[cv] = float(np.mean(vals))

    ok: dict[str, bool] = {}

    # ---------------------------------------------------------------- P4 first
    p("")
    p(f"=== P4 start-independence at every CV (bar {P4_BAR}) ===")
    worst4, n4 = 0.0, 0
    for cv in CVS:
        a = means.get((round(cv, 3), 0.55))
        b = means.get((round(cv, 3), 0.64))
        if a is None or b is None:
            p(f"   cv {cv:.1f}: UNREADABLE (one preparation has no usable block)")
            continue
        n4 += 1
        gap = abs(a - b)
        worst4 = max(worst4, gap)
        p(f"   cv {cv:.1f}: {a:.6f} vs {b:.6f}   gap {gap:.6f}"
          f"{'' if gap <= P4_BAR else '   🔴 OVER BAR'}")
    ok["P4"] = n4 == len(CVS) and worst4 <= P4_BAR
    p(f"   worst gap {worst4:.6f} on {n4} of {len(CVS)} CVs -> "
      f"{'PASS -- pooling the two preparations is legitimate' if ok['P4'] else 'NOT PASSED -- quote the preparations SEPARATELY below'}")

    # ---------------------------------------------------------------- P1
    p("")
    p(f"=== P1 🔑 the contrast, from a route with no phi_J in it "
      f"(|d - {P1_REF:+.4f}| <= {P1_BAND}) ===")
    if 0.0 in per_cv and 0.5 in per_cv:
        d = per_cv[0.5] - per_cv[0.0]
        p("   phi_ss: " + "  ".join(f"cv{cv:.1f}={per_cv[cv]:.6f}"
                                     for cv in sorted(per_cv)))
        p(f"   delta phi_c (CV 0 -> 0.5) = {d:+.6f}")
        p(f"   common ladder: E4 +0.0554, E4m +0.0497, midpoint {P1_REF:+.4f}")
        ok["P1"] = abs(d - P1_REF) <= P1_BAND
        p(f"   |d - midpoint| = {abs(d - P1_REF):.6f} (bar {P1_BAND}) -> "
          f"{'PASS -- two structurally different routes agree' if ok['P1'] else '🔴 FAIL -- the two routes DISAGREE; this puts the headline in doubt, not just the new route'}")
        # the same contrast, per preparation, since a pooled mean can hide a split
        for s in STARTS:
            a = means.get((0.0, round(s, 3)))
            b = means.get((0.5, round(s, 3)))
            if a is not None and b is not None:
                p(f"   (from phi0 = {s:.2f} alone: {b - a:+.6f})")
    else:
        ok["P1"] = False
        p("   🔴 UNRESOLVED: CV 0 or CV 0.5 has no usable block (§5)")

    # ---------------------------------------------------------------- P2
    p("")
    p(f"=== P2 CV 0 reproduces the shakeout across N 600 -> 1500 (bar {P2_BAR}) ===")
    if 0.0 in per_cv:
        dd = abs(per_cv[0.0] - P2_SHAKEOUT)
        ok["P2"] = dd <= P2_BAR
        p(f"   phi_ss(CV 0) = {per_cv[0.0]:.6f} vs shakeout {P2_SHAKEOUT:.4f} "
          f"(N = 600)   |d| = {dd:.6f} -> "
          f"{'PASS' if ok['P2'] else 'FAIL -- a finite-size effect ON THIS PROTOCOL, and the shakeout S7 corroboration of E4m rests on N = 600'}")
    else:
        ok["P2"] = False
        p("   UNRESOLVED")

    # ---------------------------------------------------------------- P3
    p("")
    p("=== P3 steady shear stays denser than isotropic compression, at every CV ===")
    iso = phi_iso_at_p()
    p("   isotropic phi at p = 3e-04 (recomputed): "
      + "  ".join(f"cv{cv:.1f}={iso[cv]:.6f}" for cv in sorted(iso)))
    n_pos, n_tot = 0, 0
    for cv in sorted(per_cv):
        if cv not in iso:
            p(f"   cv {cv:.1f}: no isotropic value at this CV and pressure")
            continue
        n_tot += 1
        off = per_cv[cv] - iso[cv]
        n_pos += off > 0
        p(f"   cv {cv:.1f}: phi_ss - phi_iso = {off:+.6f}"
          f"{'' if off > 0 else '   🔴 SIGN FLIP'}")
    ok["P3"] = n_tot > 0 and n_pos == n_tot
    p(f"   positive on {n_pos} of {n_tot} -> "
      f"{'PASS' if ok['P3'] else 'NOT PASSED -- a sign flip is a NEW statement about polydispersity and the critical state, not a bug'}")

    # ---------------------------------------------------------------- P5
    p("")
    p("=== P5 alpha: is the branch asymmetry bigger than the seed scatter? ===")
    n_sign = sum(1 for r, _w in ok_cells if r["alpha_sign_ok"])
    p(f"   fitted rate negative on {n_sign} of {len(ok_cells)} admitted cells")
    for cv in CVS:
        row = []
        for s in STARTS:
            v = [r["alpha"] for r in blocks.get((round(cv, 3), round(s, 3)), [])
                 if np.isfinite(r["alpha"])]
            row.append((s, np.array(v, float)))
        if all(a.size >= 2 for _s, a in row):
            (s1, a1), (s2, a2) = row
            within = float(np.mean([a1.std(ddof=1), a2.std(ddof=1)]))
            between = abs(float(a1.mean() - a2.mean()))
            p(f"   cv {cv:.1f}: alpha({s1:.2f}) = {a1.mean():.3g} +- "
              f"{a1.std(ddof=1):.3g}, alpha({s2:.2f}) = {a2.mean():.3g} +- "
              f"{a2.std(ddof=1):.3g};  between/within = "
              f"{between / within if within > 0 else float('inf'):.2f}")
        else:
            p(f"   cv {cv:.1f}: fewer than 2 admitted seeds on a branch -- "
              f"the asymmetry is not testable here")
    p("   (reported, not gated: §4 P5)")

    p("")
    p("=== SUMMARY ===")
    for g in ("P1", "P2", "P3", "P4"):
        p(f"   {g}: {'PASS' if ok.get(g) else 'NOT PASSED'}")
    p("   P5: reported")
    path = os.path.join(OUT, "constp_prod.txt")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwritten -> {path}")
    return 0 if ok.get("P1") else 1


if __name__ == "__main__":
    raise SystemExit(main())
