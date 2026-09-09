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

    python measure_z_iso.py
"""

from __future__ import annotations

import glob
import json
import os
from collections import defaultdict

import numpy as np

import analyse_phase8 as A

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "out")

N_REQ = 1500            # the production system size
N_LOW_RUNGS = 5         # how many lowest-pressure rungs the extrapolation uses
MIN_RUNGS = 3


def _load_phase2() -> list[dict]:
    recs = []
    for path in sorted(glob.glob(os.path.join(OUT_DIR, "phase2_*.jsonl"))):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if "error" in r or not r.get("converged"):
                    continue
                recs.append(r)
    return recs


def _z_iso(rs: list[dict]) -> tuple[float, int, float]:
    """Extrapolate z_nr to p -> 0 along z = z_iso + a sqrt(p), on the lowest rungs."""
    p_ = np.array([r["p_hat"] for r in rs], float)
    z = np.array([r["z_nr"] for r in rs], float)
    ok = np.isfinite(p_) & np.isfinite(z) & (p_ > 0)
    if ok.sum() < MIN_RUNGS:
        return float("nan"), int(ok.sum()), float("nan")
    x, y = np.sqrt(p_[ok]), z[ok]
    keep = x <= np.sort(x)[min(N_LOW_RUNGS - 1, len(x) - 1)]
    coef = np.polyfit(x[keep], y[keep], 1)
    return float(coef[1]), int(keep.sum()), float(p_[ok].min())


def main() -> int:
    lines: list[str] = []

    def p(m: str = "") -> None:
        print(m)
        lines.append(m)

    recs = _load_phase2()
    by = defaultdict(list)
    for r in recs:
        by[(r.get("family", "lognormal"), round(float(r["cv"]), 3),
            round(float(r["mu_p"]), 3), int(r["n_requested"]))].append(r)

    p(f"{len(recs)} converged isotropic cells; extrapolating z_nr to p -> 0 on the "
      f"{N_LOW_RUNGS} lowest rungs of each ladder, N = {N_REQ}")
    p(f"\n{'family':11s}{'cv':>5s}{'mu_p':>6s}{'rungs':>6s}{'p_min':>10s}{'z_iso':>8s}")
    per_mu = defaultdict(list)
    per_cv = defaultdict(list)
    for k in sorted(by):
        if k[3] != N_REQ:
            continue
        zi, n, pmin = _z_iso(by[k])
        if not np.isfinite(zi):
            continue
        per_mu[k[2]].append(zi)
        if abs(k[2] - 0.3) < 1e-9:
            per_cv[(k[0], k[1])].append(zi)
        p(f"{k[0]:11s}{k[1]:5.2f}{k[2]:6.2f}{n:6d}{pmin:10.2e}{zi:8.3f}")

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

    path = os.path.join(OUT_DIR, "z_iso_measured.txt")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwritten -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
