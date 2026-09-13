"""grade_constp_rm.py — grade PREREG_CONSTP_RATEMATCH_2026-09-12.md, R1-R5.

    python grade_constp_rm.py          -> out/constp_ratematch.txt

Every bar is quoted from that file's §3. Admission is PREREG_CONSTP_PROD §3 with the
relative z floor of §2 here. Exit 0 only if R1, R2 and R3 pass (R4 is a sign prediction
that is reported either way; R5 grades a quantity retracted after the pre-registration --
it is graded AS REGISTERED and its verdict is printed, but it cannot decide anything).
"""

from __future__ import annotations

import glob
import json
import os
from collections import defaultdict

import numpy as np

from grade_constp import stability
from grade_constp_prod import (ADMIT_P_FLOOR_FRAC, ADMIT_P_TOL, N_REQUIRED, phi_iso_at_p)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
RM_PATTERN = "phase9_rm_*.jsonl"
PROD_PATTERN = "phase9_prod*.jsonl"
CVS = (0.0, 0.2, 0.5)
STARTS = (0.55, 0.64)
MIN_SEEDS = 2

# --- §2 admission: PROD §3 with the relative z floor ---------------------------------
ADMIT_Z_REL_MIN = 0.5
# --- §3 bars -------------------------------------------------------------------------
R1_SPREAD_MAX = 1.20          # GDOT_MATCH_TOL
R1_I_MAX = 1.0e-2
R1_PROD_SPREAD = 4.70         # the unmatched production spread, for the record
R2_CENTRE, R2_HALF = 0.0525, 0.015
R2_UNMATCHED = 0.0412         # production Delta phi_c(0 -> 0.5), which R2 must move UP from
R3_STEPS_CV02 = 200000
R5_ALPHA_PROD = {0.55: (6.21, 1.14), 0.64: (3.49, 0.68)}   # CV 0, by phi0, from the prereg
R5_N_SD = 2.0
R5_ASYM_MIN = 2.0
P_STAR = 3.0e-4


def load(pattern: str) -> list[dict]:
    recs, seen = [], set()
    for path in sorted(glob.glob(os.path.join(OUT, pattern))):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                r = json.loads(line)
                if r.get("key") in seen:
                    continue
                seen.add(r["key"])
                if int(r.get("n_requested", -1)) != N_REQUIRED:
                    continue
                recs.append(r)
    return recs


def admit(r: dict) -> tuple[bool, str]:
    if r.get("error"):
        return False, "errored"
    d = stability(r)
    rel = d["z_win_min"] / float(r["z_ss"]) if r.get("z_ss") else float("nan")
    if not (np.isfinite(rel) and rel >= ADMIT_Z_REL_MIN):
        return False, f"UNJAMMED (z_win_min/z_ss {rel:.2f})"
    if d["p_win_min"] < ADMIT_P_FLOOR_FRAC * float(r["p_hat_target"]):
        return False, f"p visited {d['p_win_min']:.1e}, below 0.1 p*"
    if abs(r["p_ss_hat"] / r["p_hat_target"] - 1) > ADMIT_P_TOL:
        return False, f"servo off by {r['p_ss_hat'] / r['p_hat_target'] - 1:+.1%}"
    if not r.get("steady"):
        return False, f"not steady (p {r['steady_split']:.1%}, phi {r['phi_ss_split']:.1%})"
    if not r.get("quasistatic"):
        return False, f"not quasi-static (I {r['i_ss']:.1e})"
    return True, ""


def cell(r: dict) -> tuple:
    return (round(float(r["cv"]), 3), round(float(r["phi_target"]), 3), int(r["seed"]))


def main() -> int:
    lines: list[str] = []

    def p(m: str = "") -> None:
        print(m)
        lines.append(m)

    rm = load(RM_PATTERN)
    prod = load(PROD_PATTERN)
    p("=== grade_constp_rm.py -- PREREG_CONSTP_RATEMATCH_2026-09-12.md, R1-R5 ===")
    p(f"   {len(rm)} rate-matched records at N = {N_REQUIRED} (18 expected); "
      f"{len(prod)} production records for R3/R5")
    if not rm:
        p("   no records -- nothing to grade")
        return 1

    ok_cells, refused = [], []
    for r in rm:
        good, why = admit(r)
        (ok_cells if good else refused).append((r, why))
    p("")
    p(f"{'cv':>5s}{'phi0':>6s}{'seed':>5s}{'steps':>8s}{'gdot':>11s}{'I':>9s}{'phi_ss':>10s}"
      f"{'p_err':>8s}{'z_rel':>7s}{'alpha':>8s}{'adm':>5s}  why refused")
    for r, why in sorted(ok_cells + refused, key=lambda t: cell(t[0])):
        d = stability(r)
        rel = d["z_win_min"] / float(r["z_ss"]) if r.get("z_ss") else float("nan")
        p(f"{r['cv']:5.2f}{r['phi_target']:6.2f}{r['seed']:5d}{r['shear_steps']:8d}"
          f"{r['gdot']:11.4e}{r['i_ss']:9.2e}{r['phi_ss']:10.6f}"
          f"{r['p_ss_hat'] / r['p_hat_target'] - 1:+8.1%}{rel:7.2f}{r['alpha']:8.2f}"
          f"{('yes' if not why else 'no'):>5s}  {why}")
    p(f"   admitted {len(ok_cells)} of {len(rm)}; refused {len(refused)}")

    ok: dict[str, bool] = {}

    # ---------------------------------------------------------------- R1
    p("")
    p("=== R1  the rate is matched ===")
    gd = defaultdict(list)
    for r in rm:
        gd[round(float(r["cv"]), 3)].append(float(r["gdot"]))
    means = {cv: float(np.mean(v)) for cv, v in gd.items()}
    for cv in sorted(means):
        p(f"   cv {cv:.1f}: gdot = {means[cv]:.4e}  (n {len(gd[cv])})")
    spread = max(means.values()) / min(means.values())
    i_max = max(float(r["i_ss"]) for r in rm)
    p(f"   between-CV spread max/min = {spread:.3f} against bar {R1_SPREAD_MAX} "
      f"(production: {R1_PROD_SPREAD})")
    p(f"   max I over all cells = {i_max:.2e} against bar {R1_I_MAX:.0e}")
    ok["R1"] = spread <= R1_SPREAD_MAX and i_max <= R1_I_MAX
    p(f"   R1 {'PASS' if ok['R1'] else 'FAIL'}")

    # ---------------------------------------------------------------- blocks
    blocks: dict[tuple, list[float]] = defaultdict(list)
    for r, _ in ok_cells:
        blocks[(round(float(r["cv"]), 3), round(float(r["phi_target"]), 3))].append(
            float(r["phi_ss"]))
    p("")
    p(f"=== blocks (>= {MIN_SEEDS} admitted seeds to be usable) ===")
    bmean: dict[tuple, float] = {}
    for cv in CVS:
        for s in STARTS:
            k = (round(cv, 3), round(s, 3))
            v = np.array(blocks.get(k, []), float)
            if v.size == 0:
                p(f"   cv {cv:.1f} phi0 {s:.2f}: NO admitted cell")
                continue
            sd = float(v.std(ddof=1)) if v.size > 1 else float("nan")
            usable = v.size >= MIN_SEEDS
            p(f"   cv {cv:.1f} phi0 {s:.2f}: n = {v.size}, phi_ss = {v.mean():.6f} +- {sd:.6f}"
              f"   {'usable' if usable else 'UNRESOLVED (too few seeds)'}")
            if usable:
                bmean[k] = float(v.mean())
    per_cv: dict[float, float] = {}
    per_cv_sd: dict[float, float] = {}
    for cv in CVS:
        vals = [bmean[k] for k in bmean if abs(k[0] - cv) < 1e-9]
        allv = [x for k in blocks if abs(k[0] - cv) < 1e-9 for x in blocks[k]]
        if len(vals) == 2:     # both preparations, as P4 of the production run licensed
            per_cv[cv] = float(np.mean(vals))
            per_cv_sd[cv] = float(np.std(allv, ddof=1))
            p(f"   phi_c(CV {cv:.1f}) = {per_cv[cv]:.6f}  (pooled over both starts, "
              f"sd over {len(allv)} cells {per_cv_sd[cv]:.6f})")
        else:
            p(f"   phi_c(CV {cv:.1f}): UNRESOLVED ({len(vals)} usable block(s))")

    # ---------------------------------------------------------------- R2
    p("")
    p("=== R2  Delta phi_c(CV 0 -> 0.5) ===")
    if 0.0 in per_cv and 0.5 in per_cv:
        dphi = per_cv[0.5] - per_cv[0.0]
        lo, hi = R2_CENTRE - R2_HALF, R2_CENTRE + R2_HALF
        p(f"   Delta phi_c = {dphi:+.4f}; bar [{lo:+.4f}, {hi:+.4f}]; unmatched production "
          f"{R2_UNMATCHED:+.4f}")
        in_band = lo <= dphi <= hi
        up = dphi > R2_UNMATCHED
        p(f"   in band: {in_band}; moved UP from the unmatched value: {up} "
          f"(by {dphi - R2_UNMATCHED:+.4f})")
        if 0.2 in per_cv:
            d02 = per_cv[0.2] - per_cv[0.0]
            p(f"   Delta phi_c(0 -> 0.2) = {d02:+.4f}; (0.2 -> 0.5) = {per_cv[0.5] - per_cv[0.2]:+.4f}"
              f" -- per unit CV: {d02 / 0.2:+.3f} then {(per_cv[0.5] - per_cv[0.2]) / 0.3:+.3f}"
              f" (REPORTED: the move is not linear in CV)")
        ok["R2"] = in_band and up
        if not up:
            p("   the direction disagrees with the constant-volume control's -- a finding, reported")
    else:
        ok["R2"] = False
        p("   UNRESOLVED -- a CV level has fewer than two usable blocks")
    p(f"   R2 {'PASS' if ok['R2'] else 'FAIL'}")

    # ---------------------------------------------------------------- R3
    p("")
    p("=== R3  bit-identity of CV 0.2 against production (phi_ss, p_ss_hat, alpha) ===")
    prod_by = {cell(r): r for r in prod if int(r["shear_steps"]) == R3_STEPS_CV02}
    n_pairs, n_ident = 0, 0
    for r in sorted(rm, key=cell):
        if abs(float(r["cv"]) - 0.2) > 1e-9:
            continue
        q = prod_by.get(cell(r))
        if q is None:
            p(f"   {cell(r)}: no production twin")
            continue
        n_pairs += 1
        same = all(float(r[k]) == float(q[k]) for k in ("phi_ss", "p_ss_hat", "alpha"))
        n_ident += same
        p(f"   cv 0.2 phi0 {r['phi_target']:.2f} s{r['seed']}: steps {r['shear_steps']} vs "
          f"{q['shear_steps']}; phi_ss {r['phi_ss']:.9g} vs {q['phi_ss']:.9g}; p_ss "
          f"{r['p_ss_hat']:.9g} vs {q['p_ss_hat']:.9g}; alpha {r['alpha']:.9g} vs "
          f"{q['alpha']:.9g}  -> {'IDENTICAL' if same else 'DIFFER'}")
    ok["R3"] = n_pairs == 6 and n_ident == n_pairs
    p(f"   {n_ident} of {n_pairs} pairs bit-identical on all three continuous channels")
    p(f"   R3 {'PASS' if ok['R3'] else 'FAIL'}")

    # ---------------------------------------------------------------- R4
    p("")
    p("=== R4  sign of phi_ss - phi_iso(p*) at CV 0.5, on admissible cells ===")
    phi_iso = phi_iso_at_p(P_STAR)
    r4 = None
    for cv in CVS:
        if cv in per_cv and cv in phi_iso:
            off = per_cv[cv] - phi_iso[cv]
            p(f"   cv {cv:.1f}: phi_c {per_cv[cv]:.4f} - phi_iso {phi_iso[cv]:.4f} = {off:+.4f}")
            if cv == 0.5:
                r4 = off
        else:
            p(f"   cv {cv:.1f}: phi_iso {'missing' if cv not in phi_iso else 'ok'}, "
              f"phi_c {'missing' if cv not in per_cv else 'ok'}")
    if r4 is None:
        p("   R4 UNRESOLVED")
        ok["R4"] = False
    else:
        ok["R4"] = r4 < 0
        p(f"   predicted negative at CV 0.5: {'PASS' if ok['R4'] else 'FAIL'} ({r4:+.4f})")
        p("   (phi_iso is N = 1500, these are N = 1000: the magnitude mixes finite size in)")

    # ---------------------------------------------------------------- R5
    p("")
    p("=== R5  alpha(CV 0) at 1.96x the production rate -- GRADED AS REGISTERED ===")
    p("   NOTE: alpha was RETRACTED after this pre-registration (derivations/alpha_from_dem.md):")
    p("   the fitted rate is the servo clamp on the phi0 = 0.64 side and the servo compressing")
    p("   an unjammed pack on the phi0 = 0.55 side, and the prereg's branch labels are swapped")
    p("   (phi0 = 0.55 is BELOW phi_c and contracts). The bar is graded on the numbers as")
    p("   registered; whatever it says is about the servo, not about a material coefficient.")
    r5_ok = True
    branch_means = {}
    for s in STARTS:
        v = np.array([float(r["alpha"]) for r, _ in ok_cells
                      if abs(float(r["cv"])) < 1e-9 and abs(float(r["phi_target"]) - s) < 1e-9
                      and np.isfinite(r["alpha"])], float)
        if v.size < MIN_SEEDS:
            p(f"   phi0 {s:.2f}: too few cells ({v.size})")
            r5_ok = False
            continue
        m0, sd0 = R5_ALPHA_PROD[s]
        within = abs(v.mean() - m0) <= R5_N_SD * sd0
        branch_means[s] = (float(v.mean()), float(v.std(ddof=1)))
        p(f"   phi0 {s:.2f} (prereg label {'dilating' if s == 0.55 else 'contracting'}, "
          f"actual {'CONTRACTING' if s == 0.55 else 'DILATING'}): alpha = {v.mean():.2f} +- "
          f"{v.std(ddof=1):.2f} (n {v.size}) vs production {m0:.2f} +- {sd0:.2f}; "
          f"|diff| = {abs(v.mean() - m0):.2f} against {R5_N_SD * sd0:.2f} -> "
          f"{'within' if within else 'OUTSIDE'}")
        r5_ok &= within
    if len(branch_means) == 2:
        a, b = branch_means[0.55], branch_means[0.64]
        between = abs(a[0] - b[0])
        within_sd = 0.5 * (a[1] + b[1])
        asym = between / within_sd if within_sd > 0 else float("inf")
        p(f"   asymmetry between/within = {between:.2f}/{within_sd:.2f} = {asym:.2f} against "
          f"{R5_ASYM_MIN} (production 2.98, 3.62)")
        r5_ok &= asym >= R5_ASYM_MIN
    ok["R5"] = r5_ok
    p(f"   R5 {'PASS' if ok['R5'] else 'FAIL'} as registered -- and decides nothing (retracted quantity)")

    # ---------------------------------------------------------------- reported
    p("")
    p("=== REPORTED (not a bar): rate sensitivity a = -dphi_ss/dI from the two runs ===")
    p("   two-point differences, ALL production cells (refused ones included, since the")
    p("   refusal was the rate itself); the CV 0.5 value is the one the sign flip rested on")
    prod_cv = defaultdict(list)
    for r in prod:
        if not r.get("error"):
            prod_cv[round(float(r["cv"]), 3)].append((float(r["i_ss"]), float(r["phi_ss"])))
    rm_cv = defaultdict(list)
    for r in rm:
        rm_cv[round(float(r["cv"]), 3)].append((float(r["i_ss"]), float(r["phi_ss"])))
    for cv in CVS:
        if cv not in prod_cv or cv not in rm_cv:
            continue
        pi = np.array(prod_cv[cv]); ri = np.array(rm_cv[cv])
        dI = pi[:, 0].mean() - ri[:, 0].mean()
        dphi = pi[:, 1].mean() - ri[:, 1].mean()
        # error on the difference of two means, seed scatter pooled
        se = float(np.hypot(pi[:, 1].std(ddof=1) / np.sqrt(len(pi)),
                            ri[:, 1].std(ddof=1) / np.sqrt(len(ri))))
        if abs(dI) < 1e-4:
            p(f"   cv {cv:.1f}: I unchanged ({ri[:, 0].mean():.2e} -> {pi[:, 0].mean():.2e}), "
              f"no lever (dphi = {dphi:+.4f} +- {se:.4f})")
            continue
        a = -dphi / dI
        p(f"   cv {cv:.1f}: I {ri[:, 0].mean():.2e} -> {pi[:, 0].mean():.2e}, phi_ss "
          f"{ri[:, 1].mean():.4f} -> {pi[:, 1].mean():.4f} (dphi {dphi:+.4f} +- {se:.4f}); "
          f"a = {a:.2f} +- {se / abs(dI):.2f}")

    # ---------------------------------------------------------------- verdict
    p("")
    p("=== verdict ===")
    for k in ("R1", "R2", "R3", "R4", "R5"):
        p(f"   {k}: {'PASS' if ok[k] else 'FAIL'}")
    decisive = ok["R1"] and ok["R2"] and ok["R3"]
    p(f"   rate-matched Delta phi_c(0 -> 0.5) is {'a MEASUREMENT' if decisive else 'NOT quotable'}")

    with open(os.path.join(OUT, "constp_ratematch.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0 if decisive else 1


if __name__ == "__main__":
    raise SystemExit(main())
