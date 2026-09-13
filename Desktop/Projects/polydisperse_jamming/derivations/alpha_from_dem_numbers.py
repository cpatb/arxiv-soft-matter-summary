"""alpha_from_dem_numbers.py — the dilatancy coefficient: RETRACTED on 2026-09-12, the
day it was written. This script prints why, from the same records the claim was fitted on.

    python derivations/alpha_from_dem_numbers.py   -> derivations/alpha_from_dem_numbers.txt

What the first version of this file claimed (kept printed below as "the 2026-09-12 claim",
never deleted): alpha = 6.62 +- 1.25 on the "dilating" branch and 3.31 +- 0.65 on the
"contracting" one, a 2.00x asymmetry, propagated into Pi, ID, M1's last corner and an
upgraded E-3. Three defects, found by reading the phi(gamma) series the fit ran on:

  D1  THE LABELS WERE SWAPPED. phi0 = 0.55 sits BELOW phi_c = 0.591 (CV 0), so that pack
      must COMPACT -- it is the contracting branch. The file called it "dilating".
  D2  THE ESTIMATOR READS THE SERVO, NOT THE GRAINS. A pack compressed to phi0 = 0.64
      starts at p/p* ~ 50 (CV 0); the servo edot_v = K gdot (p/p* - 1) is then clamped at
      0.5 gdot, and d ln(phi)/dgamma = -0.5 takes 0.64 to 0.589 in a strain of 0.165 --
      which is, to the sample, the strain in which the "dilation" completes. The dense
      branch's whole transit is the box controller relaxing an over-compressed pack to
      p*. The loose branch starts UNJAMMED (p = 0) and its first ~0.3 strain is the servo
      compressing a pack that offers no resistance. The straight line ln|phi - phi_ss|
      vs gamma over the first 40 % of strain is then fitted through a window that is
      mostly servo on one branch and mostly noise about phi_ss on the other -- which is
      why alpha_fit(loose) > alpha_fit(dense) while the FIT-FREE half-strain says the
      dense pack reaches phi_c ~3x FASTER.
  D3  AT CV 0.5 THERE IS NO DILATING BRANCH AT ALL: phi_ss = 0.641-0.652 and the "dense"
      start is 0.64, so |phi0 - phi_ss| < 0.012 and the fit has no signal.

Checks (each can fail; the script exits nonzero if the retraction's own evidence fails):
  L1  every phi0 = 0.55 cell has phi0 < phi_ss (CONTRACTING); every phi0 = 0.64 cell at
      CV <= 0.2 has phi0 > phi_ss (DILATING) -- the 2026-09-12 labels were backwards on
      every cell, not merely on some
  L2  on the dilating cells the servo sits on its clamp for >= 80 % of the transit to
      within 10 % of the initial gap, and the transit strain matches the clamp-only
      prediction ln(phi0/phi_10)/clamp within a factor 1.5
  L3  the fitted alpha ORDERS the branches opposite to the fit-free half-strain
  L4  at CV 0.5 every phi0 = 0.64 cell carries |phi0 - phi_ss| < 0.012
  A1  alpha is dimensionless (still true; the .md cites it)
  A2  Pi is linear in alpha (still true; the .md cites it)

Nothing here propagates a measured alpha anywhere. poro_impact.ALPHA_DIL stays at 1.0 and
the ALPHA_DIL_MEASURED_* constants added to poro_impact.py on 2026-09-12 are removed. The
protocol that CAN measure alpha is stated at the end of the output.
"""

from __future__ import annotations

import glob
import json
import math
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "out")
sys.path.insert(0, ROOT)

import shear as SH  # noqa: E402  -- SERVO_GAIN, SERVO_EDOT_CLAMP, ALPHA_FIT_FRAC

PROD_PATTERN = "phase9_prod*.jsonl"
RM_PATTERN = "phase9_rm_*.jsonl"
N_REQUIRED = 1000
# the 2026-09-12 claim, verbatim, so the numbers being retracted stay printed
CLAIM_2026_09_12 = {"dilating": (6.616, 1.249, 6), "contracting": (3.313, 0.647, 6),
                    "ratio": 2.00, "labels": {0.55: "dilating", 0.64: "contracting"}}
GAP_FRAC_TRANSIT = 0.10        # "transit" = strain to bring |phi - phi_ss| under 10 % of its start
CLAMP_TOL = 0.999              # a sample is "at the clamp" when |edot_v| >= CLAMP_TOL * clamp * gdot
MIN_GAP_FOR_A_BRANCH = 0.012   # below this the fit has no signal (L4)

_fails: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        _fails.append(name)


def load(pattern: str) -> list[dict]:
    recs, seen = [], set()
    for path in sorted(glob.glob(os.path.join(OUT, pattern))):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                t = line.strip()
                if not t:
                    continue
                r = json.loads(t)
                if r.get("error") or r.get("key") in seen:
                    continue
                seen.add(r["key"])
                if int(r.get("n_requested", -1)) != N_REQUIRED:
                    continue
                recs.append(r)
    return recs


def admitted(r: dict) -> bool:
    """PREREG_CONSTP_PROD §3, with §6's relative z floor (PREREG_CONSTP_RATEMATCH §2)."""
    g = np.asarray(r.get("gamma", []), float)
    z = np.asarray(r.get("z_series", []), float)
    pp = np.asarray(r.get("p_hat_series", []), float)
    if g.size == 0 or z.size != g.size or pp.size != g.size:
        return False
    w = g >= 0.5 * float(g.max())
    if not (r.get("z_ss") and float(z[w].min()) / float(r["z_ss"]) >= 0.5):
        return False
    if float(pp[w].min()) < 0.1 * float(r["p_hat_target"]):
        return False
    if abs(r["p_ss_hat"] / r["p_hat_target"] - 1) > 0.05:
        return False
    return bool(r.get("steady") and r.get("quasistatic"))


def transit(r: dict) -> dict:
    """Fit-free description of one cell's approach to phi_ss.

    gap0        phi(0) - phi_ss; its SIGN is the branch
    g_half      first strain at which |gap| < 0.5 |gap0|
    g_10        first strain at which |gap| < GAP_FRAC_TRANSIT |gap0|  (the "transit")
    clamp_frac  fraction of samples in [0, g_10) with the servo on its clamp
    underp_frac fraction of samples in [0, g_10) with p < 0.5 p* (servo compressing an
                unjammed pack)
    g_servo     clamp-only prediction of g_10 on the dilating side: with
                d ln(phi)/dgamma = -clamp, the strain from phi0 to phi_10 is
                ln(phi0/phi_10)/clamp
    First crossings are used; the series is checked to be monotone toward phi_ss up to
    g_10 (mono_ok), which is the 2026-09-03 lesson on crossing times.
    """
    g = np.asarray(r["gamma"], float)
    phi = np.asarray(r["phi_series"], float)
    e = np.asarray(r["edot_v_series"], float)
    p = np.asarray(r["p_hat_series"], float)
    phi_ss = float(r["phi_ss"])
    gdot = float(r["gdot"])
    p_star = float(r["p_hat_target"])
    gap = phi - phi_ss
    gap0 = float(gap[0])
    d = np.abs(gap)
    i_half = int(np.argmax(d < 0.5 * abs(gap0)))
    i_10 = int(np.argmax(d < GAP_FRAC_TRANSIT * abs(gap0)))
    ok_half = bool(d[i_half] < 0.5 * abs(gap0))
    ok_10 = bool(d[i_10] < GAP_FRAC_TRANSIT * abs(gap0))
    upto = slice(0, max(i_10, 1))
    rel = np.abs(e[upto]) / gdot
    clamp_frac = float(np.mean(rel >= CLAMP_TOL * SH.SERVO_EDOT_CLAMP))
    underp_frac = float(np.mean(p[upto] < 0.5 * p_star))
    # monotone toward phi_ss up to the transit: |gap| non-increasing within one sample's noise
    mono_ok = bool(np.all(np.diff(d[: i_10 + 1]) <= 0.002)) if ok_10 else False
    phi_10 = phi_ss + GAP_FRAC_TRANSIT * gap0
    g_servo = (math.log(phi[0] / phi_10) / SH.SERVO_EDOT_CLAMP) if gap0 > 0 else float("nan")
    return dict(gap0=gap0, g_half=float(g[i_half]) if ok_half else float("nan"),
                g_10=float(g[i_10]) if ok_10 else float("nan"), clamp_frac=clamp_frac,
                underp_frac=underp_frac, g_servo=g_servo, p0_over_pstar=float(p[0]) / p_star,
                mono_ok=mono_ok, branch="dilating" if gap0 > 0 else "contracting")


def main() -> int:
    print("=" * 96)
    print("alpha_from_dem_numbers.py -- the dilatancy coefficient: RETRACTED 2026-09-12")
    print("=" * 96)
    prod = load(PROD_PATTERN)
    rm = load(RM_PATTERN)
    allr = prod + rm
    gains = sorted({float(r["servo_gain"]) for r in allr})
    print(f"servo: K = {gains} (from the records; shear.SERVO_GAIN default is {SH.SERVO_GAIN}), "
          f"clamp = {SH.SERVO_EDOT_CLAMP} gdot,")
    print(f"       alpha fit window = first {SH.ALPHA_FIT_FRAC:.0%} of strain; "
          f"transit = |gap| < {GAP_FRAC_TRANSIT:.0%} of gap0")
    recs = [r for r in allr if admitted(r)]
    # the CV 0.2 rate-matched cells are bit-identical re-runs of the production ones (same
    # shear_steps; PREREG_CONSTP_RATEMATCH R3) -- keep one per physical cell for the stats
    seen_cell, uniq = set(), []
    for r in recs:
        cell = (round(float(r["cv"]), 2), round(float(r["phi_target"]), 2), int(r["seed"]),
                int(r["shear_steps"]))
        if cell in seen_cell:
            continue
        seen_cell.add(cell)
        uniq.append(r)
    print(f"records: {len(prod)} production + {len(rm)} rate-matched = {len(allr)}, "
          f"{len(recs)} admitted (N = {N_REQUIRED}), {len(uniq)} distinct "
          f"(cv, phi0, seed, steps) cells after dropping the bit-identical CV 0.2 twins")
    recs = uniq

    c = CLAIM_2026_09_12
    lab_dil = [k for k, v in c["labels"].items() if v == "dilating"][0]
    lab_con = [k for k, v in c["labels"].items() if v == "contracting"][0]
    print("\n=== the 2026-09-12 claim, verbatim (kept; this is what is retracted) ===")
    print(f"  alpha (dilating)    = {c['dilating'][0]:.2f} +- {c['dilating'][1]:.2f}  "
          f"(n = {c['dilating'][2]})   <- fitted on the phi0 = {lab_dil} cells")
    print(f"  alpha (contracting) = {c['contracting'][0]:.2f} +- {c['contracting'][1]:.2f}  "
          f"(n = {c['contracting'][2]})   <- fitted on the phi0 = {lab_con} cells")
    print(f"  ratio {c['ratio']:.2f}x; propagated into Pi, ID, M1's corner and E-3")

    print("\n=== per admitted cell: what the series actually do ===")
    print(f"  {'cv':>4s}{'seed':>5s}{'phi0':>6s}{'steps':>8s}{'phi_ss':>9s}{'gap0':>9s}  "
          f"{'branch':<12s}{'p0/p*':>7s}{'g_half':>8s}{'g_10':>7s}{'g_servo':>8s}"
          f"{'clamp%':>8s}{'p<p*/2%':>9s}{'mono':>6s}{'alpha_fit':>10s}")
    rows = []
    for r in sorted(recs, key=lambda r: (r["cv"], r["phi_target"], r["seed"], r["shear_steps"])):
        t = transit(r)
        t.update(cv=round(float(r["cv"]), 2), seed=int(r["seed"]),
                 phi0=round(float(r["phi_target"]), 2), steps=int(r["shear_steps"]),
                 phi_ss=float(r["phi_ss"]), alpha_fit=r.get("alpha"))
        rows.append(t)
        a = t["alpha_fit"]
        a = a if (a is not None and np.isfinite(a)) else float("nan")
        print(f"  {t['cv']:4.1f}{t['seed']:5d}{t['phi0']:6.2f}{t['steps']:8d}{t['phi_ss']:9.4f}"
              f"{t['gap0']:+9.4f}  {t['branch']:<12s}{t['p0_over_pstar']:7.1f}"
              f"{t['g_half']:8.3f}{t['g_10']:7.3f}{t['g_servo']:8.3f}"
              f"{100 * t['clamp_frac']:8.0f}{100 * t['underp_frac']:9.0f}"
              f"{('yes' if t['mono_ok'] else 'no'):>6s}{a:10.2f}")

    # ---------------------------------------------------------------- L1
    print("\n=== L1  the labels: which way does each start actually move? ===")
    lo = [t for t in rows if t["phi0"] == 0.55]
    hi_low_cv = [t for t in rows if t["phi0"] == 0.64 and t["cv"] <= 0.2]
    n_lo_contract = sum(1 for t in lo if t["branch"] == "contracting")
    n_hi_dilate = sum(1 for t in hi_low_cv if t["branch"] == "dilating")
    check("L1a every phi0 = 0.55 cell is CONTRACTING (phi0 < phi_ss)",
          n_lo_contract == len(lo) and len(lo) > 0, f"{n_lo_contract} of {len(lo)}")
    check("L1b every phi0 = 0.64 cell at CV <= 0.2 is DILATING (phi0 > phi_ss)",
          n_hi_dilate == len(hi_low_cv) and len(hi_low_cv) > 0,
          f"{n_hi_dilate} of {len(hi_low_cv)}")
    print(f"  => the 2026-09-12 file labelled phi0 = 0.55 '{c['labels'][0.55]}' and "
          f"phi0 = 0.64 '{c['labels'][0.64]}': backwards on every cell.")

    # ---------------------------------------------------------------- L2
    print("\n=== L2  the dilating transit is the servo clamp ===")
    dil = [t for t in rows if t["branch"] == "dilating" and abs(t["gap0"]) >= MIN_GAP_FOR_A_BRANCH]
    clamp_fracs = np.array([t["clamp_frac"] for t in dil])
    ratios = np.array([t["g_10"] / t["g_servo"] for t in dil])
    print(f"  dilating cells with a readable gap (|gap0| >= {MIN_GAP_FOR_A_BRANCH}): {len(dil)}")
    print(f"  fraction of the transit at the clamp: min {clamp_fracs.min():.2f}, "
          f"median {np.median(clamp_fracs):.2f}, max {clamp_fracs.max():.2f}")
    print(f"  measured transit / clamp-only prediction: min {ratios.min():.2f}, "
          f"median {np.median(ratios):.2f}, max {ratios.max():.2f}")
    check("L2a the servo is at its clamp for >= 80 % of every dilating transit",
          bool(len(dil)) and clamp_fracs.min() >= 0.80, f"min {clamp_fracs.min():.2f}")
    check("L2b the transit strain matches the clamp-only prediction within 1.5x on every cell",
          bool(len(dil)) and ratios.max() <= 1.5 and ratios.min() >= 1 / 1.5,
          f"range {ratios.min():.2f}-{ratios.max():.2f}")
    print("  => the pack is not dilating under shear; the box is expanding at its maximum rate")
    print(f"     to relax p/p* ~ {min(t['p0_over_pstar'] for t in dil):.0f}-"
          f"{max(t['p0_over_pstar'] for t in dil):.0f} to 1. alpha(dilating) is a readout of "
          f"SERVO_EDOT_CLAMP = {SH.SERVO_EDOT_CLAMP}.")
    con = [t for t in rows if t["branch"] == "contracting" and abs(t["gap0"]) >= MIN_GAP_FOR_A_BRANCH]
    underp = np.array([t["underp_frac"] for t in con])
    print(f"  contracting cells ({len(con)}): fraction of the transit spent UNDER-pressured "
          f"(p < p*/2, the servo")
    print(f"     compressing an unjammed pack): min {underp.min():.2f}, median "
          f"{np.median(underp):.2f}, max {underp.max():.2f} -- that part is not dilatancy either")

    # ---------------------------------------------------------------- L3
    print("\n=== L3  the fitted alpha orders the branches against the fit-free strain ===")
    by = defaultdict(list)
    for t in rows:
        if abs(t["gap0"]) >= MIN_GAP_FOR_A_BRANCH and t["alpha_fit"] is not None \
                and np.isfinite(t["alpha_fit"]) and np.isfinite(t["g_half"]):
            by[t["branch"]].append((t["alpha_fit"], t["g_half"]))
    a_dil = np.array([a for a, _ in by["dilating"]])
    a_con = np.array([a for a, _ in by["contracting"]])
    gh_dil = np.array([g for _, g in by["dilating"]])
    gh_con = np.array([g for _, g in by["contracting"]])
    print(f"  fitted alpha:    dilating {a_dil.mean():.2f} +- {a_dil.std(ddof=1):.2f} (n {a_dil.size}), "
          f"contracting {a_con.mean():.2f} +- {a_con.std(ddof=1):.2f} (n {a_con.size})")
    print(f"  fit-free g_half: dilating {gh_dil.mean():.3f} +- {gh_dil.std(ddof=1):.3f}, "
          f"contracting {gh_con.mean():.3f} +- {gh_con.std(ddof=1):.3f}")
    print(f"  a larger alpha is a FASTER approach, so the fit says the contracting pack is "
          f"{a_con.mean() / a_dil.mean():.2f}x faster; the half-strain says the dilating pack "
          f"is {gh_con.mean() / gh_dil.mean():.2f}x faster")
    check("L3 the two orderings contradict each other",
          (a_con.mean() > a_dil.mean()) and (gh_dil.mean() < gh_con.mean()))
    print("  => the estimator is not measuring the transit at all. On the dilating cells the")
    print("     transit is over in ~0.15 strain and the remaining ~1.05 of the fit window is")
    print("     noise about phi_ss; on the contracting cells the window straddles the servo")
    print("     compression and the material approach. Neither slope is alpha.")

    # ---------------------------------------------------------------- L4
    print("\n=== L4  CV 0.5 has no dilating branch ===")
    c50 = [t for t in rows if t["cv"] == 0.5 and t["phi0"] == 0.64]
    gaps = np.array([abs(t["gap0"]) for t in c50])
    check(f"L4 every CV 0.5 phi0 = 0.64 cell has |gap0| < {MIN_GAP_FOR_A_BRANCH}",
          bool(len(c50)) and gaps.max() < MIN_GAP_FOR_A_BRANCH,
          f"max |gap0| = {gaps.max():.4f} over {len(c50)} cells")
    print(f"  => phi_c(CV 0.5) = {min(t['phi_ss'] for t in c50):.3f}-"
          f"{max(t['phi_ss'] for t in c50):.3f} sits ABOVE the 'dense' start of 0.64. The "
          f"alpha values there")
    print("     are fits to a gap of 0.001-0.012, i.e. to seed scatter.")

    # ---------------------------------------------------------------- A1, A2 (still true)
    print("\n=== A1, A2  the algebra that still stands (and is now moot) ===")
    import pint
    u = pint.UnitRegistry()
    eps_v = (1.0 * u.dimensionless) * (1.0 / u.s) * (0.02 * u.dimensionless)
    check("A1 edot_v = alpha gdot dphi has dimensions of 1/time for dimensionless alpha",
          eps_v.dimensionality == (1 / u.s).dimensionality)
    mu_w, dphi, ell, rho, U0, D0, kappa = 1e-3, 0.02, 2.5e-3, 1440.0, 3.0, 2.5e-3, 1e-15
    pi_of = lambda a: a * mu_w * dphi * ell ** 2 / (rho * kappa * U0 * D0)  # noqa: E731
    check("A2 Pi(2 alpha) = 2 Pi(alpha) exactly",
          math.isclose(pi_of(2.0), 2 * pi_of(1.0), rel_tol=1e-12))
    print("  So a measured alpha WOULD propagate as a factor -- there just is not one.")

    # ---------------------------------------------------------------- what replaces it
    print("\n=== WITHDRAWN, STANDING, AND THE PROTOCOL THAT MEASURES alpha ===")
    gh_c = {cv: [t["g_half"] for t in con if t["cv"] == cv] for cv in (0.0, 0.2, 0.5)}
    for line in (
        "WITHDRAWN: alpha = 6.62 / 3.31, the 2.00x asymmetry, the 3.3-6.6x factor on Pi and",
        "  1/ID, the conditional 12-24x closure of M1's corner, and the E-3 'magnitude'",
        "  upgrade (dense-side suction ~2x loose-side overpressure). E-3 is a SIGN test again.",
        "STANDS: the sign of the approach (every admitted cell moves toward phi_ss) -- a",
        "  property of a stable steady state, not evidence about alpha.",
        "THE PROTOCOL: the servo must be much faster than the grains and the pack must start",
        "  ON the pressure it is sheared at. Consolidate isotropically to p = p* -- a",
        "  frictional pack for the contracting side (phi_iso(p*) < phi_c), a FRICTIONLESS",
        "  consolidation with friction switched on at the first shear step for the dilating",
        "  side (phi_iso,0(p*) > phi_c) -- and shear at p*. Then edot_v is the grains' and",
        "  phi(gamma) from phi_iso(p*) to phi_c(p*) is the Reynolds transient. Stop the fit",
        "  where |phi - phi_ss| reaches the seed scatter, not at a fixed fraction of strain.",
        "FIT-FREE NUMBERS THAT SURVIVE (contracting side, servo compression INCLUDED, so an",
        "  UPPER bound on the material strain), g_half from phi0 = 0.55:",
        f"  CV 0: {', '.join(f'{x:.3f}' for x in gh_c[0.0])};  CV 0.2: "
        f"{', '.join(f'{x:.3f}' for x in gh_c[0.2])};  CV 0.5: "
        f"{', '.join(f'{x:.3f}' for x in gh_c[0.5])}",
        "  -- from a common phi0, so the CV 0.5 pack also has ~2x the gap to close.",
    ):
        print("  " + line)

    print("\n" + "=" * 96)
    if _fails:
        print(f"{len(_fails)} CHECK(S) FAILED: {_fails}")
        return 1
    print("ALL CHECKS PASSED -- the retraction's evidence holds")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
