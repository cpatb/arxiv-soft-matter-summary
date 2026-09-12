"""measure_z_iso.py — what IS the isostatic contact number at the production friction?

WHY THIS EXISTS
---------------
`analyse_phase8.e4_z_crossing` (estimator E4) returns the packing fraction at which the
non-rattler contact number crosses `Z_ISO`. The module sets

    Z_ISO_FRICTIONLESS = 6.0    # 2d, exact for frictionless spheres
    Z_ISO_FRICTIONAL   = 4.0    # d+1, exact for INFINITELY frictional spheres

and its own comment records the anchors that pin them: "6.37-6.92 at mu_p = 0 and
3.98-4.07 at mu_p = 3". Neither anchor is at the production friction. **The production
arm runs at mu_p = 0.3**, where d+1 = 4 is a limit the material is nowhere near, and
where no theory fixes the count at all -- it lies somewhere between 4 and 6 and has to
be measured.

This measures it, from the project's own isotropic pressure ladders (`phase2_*.jsonl`),
by extrapolating `z_nr` to zero confining pressure along the standard `z - z_c ~ p^1/2`
approach to jamming. That is the same construction `phase2_phiJ.py` uses for `phi_J`, on
the same records, so it introduces no new assumption.

WHAT IT FINDS, AND WHY IT MATTERS TWICE
---------------------------------------
1. `z_iso(mu_p = 0.3)` is about 4.54, not 4.0. E4's target is out by ~0.54 contacts, and
   `dz/dphi` along the shear ladder is ~34 per unit phi, so E4's absolute `phi_c` is out
   by ~0.016 -- a third of the CV contrast it is used to measure.
2. Worse for a CONTRAST: `z_iso` itself drifts with CV (lognormal 4.61 at CV = 0 down to
   4.38 at CV = 0.5), so a single constant puts a CV-DEPENDENT bias into E4's own CV
   contrast, in the direction that inflates it.

AND THE CONSEQUENCE FOR THE WHOLE PROTOCOL
------------------------------------------
E4 was the candidate replacement for E2 after E2's root was found to follow the bottom
of its own fitting window (PREREG_ROOT_PINNING_2026-09-09.md). But E4 at mu_p = 0.3 can
only be calibrated by reading `z_iso` off the ISOTROPIC ladder -- the same measurement
the shear ladder is centred on. So at the production friction this protocol has no
estimator of `phi_c` that is independent of the isotropic `phi_J`, and the claim that
the two "agree" cannot be made from these runs at all. E4 IS theory-anchored at
`mu_p = 0`, where 2d = 6 is exact, which is exactly the arm the P4 anchor tests.

    python measure_z_iso.py                      # --cap default --rule r3
    python measure_z_iso.py --cap all --rule broken --no-dedup   # the stale behaviour

AUDITED AND REPAIRED 2026-09-12 (PREREG_ZISO_2026-09-12.md, out/ziso_audit.txt)
------------------------------------------------------------------------------
Two defects, neither of which changed the conclusion by much, and both of which had to
be measured rather than argued:

  D1  this script read the SHARED GLOB `phase2_*.jsonl`, which the `capfix` campaign
      began writing into on 2026-09-09, with no cap filter and no dedup by cell key.
      Re-running it unchanged today mixes two iteration caps and double-weights the low
      rungs: mu_p = 0.3 reads 4.549 instead of 4.535 and the mu_p = 1 arm goes from 1
      cell to 2. Hence `--cap`.
  D2  `N_LOW_RUNGS = 5` was documented as RUNGS and implemented as RECORDS. With three
      seeds per rung that is TWO distinct pressures, on ladders that offer four -- a
      two-point fit with zero residual degrees of freedom, wearing a five-point label.
      Hence `--rule`, and the rung-correct default.

What it cost: z_iso(mu_p = 0.3) moves 0.024 contacts across the three rules (4.535 /
4.533 / 4.510), i.e. 0.0011 in phi, and E4m's CV contrast moves +0.0003. The code was
wrong and the number it produced was right to within a twentieth of its own error bar.

What it found: the CV drift in z_iso is a LOGNORMAL effect (-0.213) and essentially
absent on the zero-skew families (uniform -0.056, bidisperse -0.012) -- the project's
own CV/skewness separation, from an observable with no phi_J in it. And uniform's
apparent -0.193 under the shipped rule rested entirely on a two-pressure ladder.

The repaired-cap table is NOT adopted: it covers CV 0.4/0.5 and nothing at CV 0 or 0.2,
so a half-adopted table would put the iteration cap on the CV axis under test.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
from collections import defaultdict

import numpy as np

import analyse_phase8 as A

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "out")

N_REQ = 1500            # the production system size
MIN_RUNGS = 3           # fewest DISTINCT pressures an extrapolation may be fitted on

# 🔴 2026-09-12, D2 of PREREG_ZISO_2026-09-12.md. This constant was named
# `N_LOW_RUNGS = 5` and documented as "how many lowest-pressure rungs the extrapolation
# uses", while the code kept the 5 lowest *RECORDS*. Every ladder carries 3 seeds per
# rung, so 5 records is TWO distinct pressures: the intercept the whole E4m correction
# rests on was a straight line through two pressures with zero residual degrees of
# freedom, on ladders that have four pressures available. The rule is now explicit and
# selectable, the number counts rungs, and the shipped behaviour is kept as "broken" so
# that the old artifact can still be reproduced.
RULES = {
    "broken": None,     # the shipped rule: 5 lowest RECORDS (~2 distinct pressures)
    "r2": 2,            # every record at the 2 lowest distinct rungs
    "r3": 3,            # every record at the 3 lowest distinct rungs
}
N_LOW_RECORDS_BROKEN = 5
# The rule SELECTED by PREREG_ZISO_2026-09-12.md §4, using the one anchor outside this
# code: z_iso(mu_p = 0) = 2d = 6 exactly. R3 lands 0.132 from it against R2's 0.148,
# with an sd of 0.027 against 0.060 -- with only four rungs the two-rung fit is
# NOISE-limited, not curvature-limited, which is the opposite of what §5 predicted.
# `--rule broken` is kept so the historical artifact can be reproduced ON PURPOSE.
DEFAULT_RULE = "r3"

# 🔴 D1. `phase2_*.jsonl` is a SHARED GLOB that two campaigns write into: 537 records at
# the default iteration cap and 180 from `capfix` at PJ_MAX_ITER = 1e6. The shipped
# loader filtered on `converged` alone, so a ladder was fitted through both protocols at
# once with its low rungs double-weighted, and no dedup by cell key was done at all.
DEFAULT_MAX_ITER = 200000


def _cap_of(r: dict) -> int:
    """The iteration cap this record ran under. Absent means the default."""
    mi = r.get("max_iter")
    return DEFAULT_MAX_ITER if mi is None else int(mi)


def load_phase2(cap: str = "default", dedup: bool = True) -> tuple[list[dict], dict]:
    """Converged isotropic records, selected by iteration cap and deduplicated by key.

    `cap` is one of "default" (the cap every published number was measured under),
    "repaired" (the `capfix` campaign only) or "all" (what the shipped script did, kept
    so that the stale artifact can be reproduced on purpose).

    Returns (records, info). `info` carries the counts and the duplicate-state
    conflicts, which are printed by the caller -- a selection that silently drops or
    double-counts cells is the defect this function exists to make visible.
    """
    if cap not in ("default", "repaired", "all"):
        raise ValueError(f"cap must be default|repaired|all, got {cap!r}")
    kept: dict[str, dict] = {}
    recs: list[dict] = []
    n_lines = n_err = n_unconv = n_cap = n_dupe = 0
    conflicts: list[tuple[str, float, float]] = []
    files = sorted(glob.glob(os.path.join(OUT_DIR, "phase2_*.jsonl")))
    for path in files:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                n_lines += 1
                if "error" in r:
                    n_err += 1
                    continue
                if not r.get("converged"):
                    n_unconv += 1
                    continue
                is_default = _cap_of(r) == DEFAULT_MAX_ITER
                if (cap == "default" and not is_default) or \
                        (cap == "repaired" and is_default):
                    n_cap += 1
                    continue
                if not dedup:
                    recs.append(r)
                    continue
                k = r.get("key")
                if k is None:
                    recs.append(r)
                    continue
                if k in kept:
                    n_dupe += 1
                    prev = float(kept[k].get("phi", np.nan))
                    cur = float(r.get("phi", np.nan))
                    if not np.isclose(prev, cur, rtol=0, atol=1e-12, equal_nan=True):
                        conflicts.append((k, prev, cur))
                    continue
                kept[k] = r
    if dedup:
        recs = list(kept.values())
    info = {"files": len(files), "lines": n_lines, "errors": n_err,
            "unconverged": n_unconv, "wrong_cap": n_cap, "dupes": n_dupe,
            "conflicts": conflicts, "n": len(recs), "cap": cap, "dedup": dedup}
    return recs, info


def z_iso_fit(rs: list[dict], rule: str = "r2") -> dict:
    """Extrapolate z_nr to p -> 0 along z = z_iso + a sqrt(p), on the lowest rungs.

    `rule` selects which points enter the fit:
      "broken" -- the 5 lowest RECORDS, i.e. the shipped behaviour (see RULES).
      "r2"/"r3" -- every record at the 2 or 3 lowest DISTINCT rungs.

    Returns the intercept plus what it was fitted on, because "how many rungs" is the
    thing that was wrong here and a caller must be able to see it.
    """
    if rule not in RULES:
        raise ValueError(f"rule must be one of {sorted(RULES)}, got {rule!r}")
    bad = {"z_iso": float("nan"), "n_pts": 0, "n_rungs": 0, "p_min": float("nan"),
           "p_max_fit": float("nan"), "rule": rule}
    p_ = np.array([r["p_hat"] for r in rs], float)
    z = np.array([r["z_nr"] for r in rs], float)
    tgt = np.array([r.get("p_hat_target", np.nan) for r in rs], float)
    ok = np.isfinite(p_) & np.isfinite(z) & (p_ > 0)
    if ok.sum() < MIN_RUNGS:
        return bad
    x, y, t = np.sqrt(p_[ok]), z[ok], tgt[ok]
    if rule == "broken":
        keep = x <= np.sort(x)[min(N_LOW_RECORDS_BROKEN - 1, len(x) - 1)]
    else:
        # group by the TARGET pressure, which is the rung label; the realised p_hat
        # scatters within a rung and sorting on it mixes neighbouring rungs.
        label = t if np.all(np.isfinite(t)) else p_[ok]
        rungs = np.unique(label)
        if len(rungs) < MIN_RUNGS:
            return bad
        keep = np.isin(label, rungs[:RULES[rule]])
    n_rungs = int(len(np.unique((t if np.all(np.isfinite(t)) else p_[ok])[keep])))
    if keep.sum() < 2 or n_rungs < 2:
        return bad
    coef = np.polyfit(x[keep], y[keep], 1)
    return {"z_iso": float(coef[1]), "n_pts": int(keep.sum()), "n_rungs": n_rungs,
            "p_min": float(p_[ok].min()), "p_max_fit": float((x[keep] ** 2).max()),
            "rule": rule}


def main(cap: str = "default", rule: str = DEFAULT_RULE, dedup: bool = True,
         out: str | None = None, quiet: bool = False) -> dict:
    """Measure z_iso under one (cap, rule) and return the per-arm and per-CV tables."""
    lines: list[str] = []

    def p(m: str = "") -> None:
        if not quiet:
            print(m)
        lines.append(m)

    recs, info = load_phase2(cap=cap, dedup=dedup)
    by = defaultdict(list)
    for r in recs:
        by[(r.get("family", "lognormal"), round(float(r["cv"]), 3),
            round(float(r["mu_p"]), 3), int(r["n_requested"]))].append(r)

    p(f"cap = {cap!r}   rule = {rule!r}   dedup = {dedup}")
    p(f"selection: {info['lines']} lines in {info['files']} files -> {info['n']} cells "
      f"({info['errors']} errors, {info['unconverged']} unconverged, "
      f"{info['wrong_cap']} wrong cap, {info['dupes']} duplicate keys dropped)")
    if info["conflicts"]:
        w = max(info["conflicts"], key=lambda c: abs(c[2] - c[1]))
        p(f"   🔴 {len(info['conflicts'])} duplicated keys return a DIFFERENT phi; "
          f"worst |dphi| = {abs(w[2] - w[1]):.2e} ({w[0]})")
    how = ("5 lowest RECORDS -- the shipped rule" if rule == "broken"
           else f"every record at the {RULES[rule]} lowest distinct rungs")
    p(f"extrapolating z_nr to p -> 0 under rule {rule!r} ({how}), N = {N_REQ}")
    p(f"\n{'family':11s}{'cv':>5s}{'mu_p':>6s}{'pts':>5s}{'rungs':>6s}{'p_min':>10s}"
      f"{'p_fit_max':>11s}{'z_iso':>8s}")
    per_mu = defaultdict(list)
    per_cv = defaultdict(list)
    for k in sorted(by):
        if k[3] != N_REQ:
            continue
        f = z_iso_fit(by[k], rule=rule)
        zi = f["z_iso"]
        if not np.isfinite(zi):
            continue
        per_mu[k[2]].append(zi)
        if abs(k[2] - 0.3) < 1e-9:
            per_cv[(k[0], k[1])].append(zi)
        p(f"{k[0]:11s}{k[1]:5.2f}{k[2]:6.2f}{f['n_pts']:5d}{f['n_rungs']:6d}"
          f"{f['p_min']:10.2e}{f['p_max_fit']:11.2e}{zi:8.3f}")

    p("\n=== z_iso by friction ===")
    p(f"   {'mu_p':>6s} {'n cells':>8s} {'z_iso':>8s} {'sd':>7s}   "
      f"{'in analyse_phase8':>18s}")
    for mu in sorted(per_mu):
        v = np.array(per_mu[mu], float)
        used = (A.Z_ISO_FRICTIONLESS if mu <= A.Z_ISO_MU_BAR else A.Z_ISO_FRICTIONAL)
        p(f"   {mu:6.2f} {len(v):8d} {v.mean():8.3f} "
          f"{v.std(ddof=1) if len(v) > 1 else float('nan'):7.3f}   {used:18.1f}")

    prod = np.array(per_mu.get(0.3, []), float)
    if prod.size:
        used = A.Z_ISO_FRICTIONAL
        p(f"\n🔴 AT THE PRODUCTION FRICTION mu_p = 0.3 the measured isostatic count is "
          f"{prod.mean():.3f} +- {prod.std(ddof=1):.3f} ({prod.size} cells), while "
          f"E4 uses {used:.1f}.")
        p(f"   The gap is {prod.mean() - used:+.3f} contacts. d+1 = 4 is the "
          f"mu_p -> infinity limit; nothing pins the count at finite friction, so 4.0 "
          f"is not a convention here, it is a wrong number.")

        # what that costs E4, using the shear ladder's own dz/dphi
        sh = [r for r in A.load("phase8_prod2*.jsonl,phase8_downlad*.jsonl")
              if "error" not in r and abs(float(r.get("mu_p", -1)) - 0.3) < 1e-9]
        byc = defaultdict(list)
        for r in sh:
            byc[round(float(r["cv"]), 3)].append(r)
        p(f"\n   what a {prod.mean() - used:+.3f} target error costs E4, using each "
          f"CV ladder's own dz/dphi:")
        p(f"   {'cv':>5s} {'dz/dphi':>9s} {'dphi_c':>9s}")
        for cv in sorted(byc):
            rs = sorted(byc[cv], key=lambda r: r["phi"])
            x = np.array([r["phi"] for r in rs], float)
            y = np.array([r["z_nr_ss"] for r in rs], float)
            if len(x) < 3:
                continue
            slope = float(np.polyfit(x, y, 1)[0])
            p(f"   {cv:5.2f} {slope:9.1f} {(prod.mean() - used) / slope:9.4f}")

    if per_cv:
        p("\n=== and z_iso is not constant in CV, which biases the CONTRAST ===")
        for fam in sorted({k[0] for k in per_cv}):
            row = [(cv, float(np.mean(v))) for (f_, cv), v in sorted(per_cv.items())
                   if f_ == fam]
            if len(row) < 2:
                continue
            txt = "  ".join(f"cv{cv:.2f}={z:.3f}" for cv, z in row)
            p(f"   {fam:11s} {txt}   drift {row[-1][1] - row[0][1]:+.3f}")
        p("   A single Z_ISO therefore puts a CV-dependent offset into E4's own CV "
          "contrast. E4 needs a per-cell target, and the only place to read one is the "
          "isotropic ladder -- which is the measurement the shear ladder is already "
          "centred on.")

    # ⚠ The historical artifact out/z_iso_measured.txt is what G2 of
    # PREREG_ZISO_2026-09-12.md reproduces against, so this never overwrites it unless
    # asked for by name.
    path = out or os.path.join(OUT_DIR, f"z_iso_{cap}_{rule}.txt")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    if not quiet:
        print(f"\nwritten -> {path}")
    return {"cap": cap, "rule": rule, "info": info, "path": path,
            "per_mu": {mu: list(v) for mu, v in per_mu.items()},
            "per_cv": {k: list(v) for k, v in per_cv.items()},
            "text": "\n".join(lines)}


def _cli() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cap", default="default", choices=["default", "repaired", "all"],
                    help="which iteration cap to select (default: the published one)")
    ap.add_argument("--rule", default=DEFAULT_RULE, choices=sorted(RULES),
                    help="which points enter the p -> 0 extrapolation")
    ap.add_argument("--no-dedup", action="store_true",
                    help="reproduce the shipped behaviour, which did not deduplicate")
    ap.add_argument("--out", default=None, help="output path")
    a = ap.parse_args()
    main(cap=a.cap, rule=a.rule, dedup=not a.no_dedup, out=a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(_cli())
