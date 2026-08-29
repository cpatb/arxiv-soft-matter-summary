#!/usr/bin/env python3
"""
compare_domain_extent.py -- does widening the radial domain change the answer?

The manuscript's weakest sentence is that the tau_c = 1 Pa cleaning contour lies
OUTSIDE the R = 25 mm computation (tau_w has fallen only to 1.957 Pa at the outer
boundary), so it can only be bounded, r > 25 mm, or reached by extrapolating the
fitted closure -- which gives 39 mm from the full band and 35 mm from an outer-half
refit. This script runs the check that settles it, in two parts:

  1. VALIDATION. On the overlap 0 < r < 25 mm the two cases must agree. If they do
     not, the wider case is not a superset of the narrower one and nothing below
     means anything -- the outflow boundary was influencing the interior.
  2. THE ANSWER. Read the tau_c = 1 Pa root directly off the R = 50 mm field and
     compare it with the two extrapolations the paper currently quotes.

Usage: python compare_domain_extent.py [narrow_case] [wide_case]
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RHO = 1000.0
D = 1.0e-3
# The paper's two extrapolated values, for the comparison this script exists to make.
EXTRAP_FULL_BAND_MM = 39.14
EXTRAP_OUTER_HALF_MM = 34.75


def load(case: str):
    path = os.path.join(HERE, case, "tau_w_r.csv")
    if not os.path.exists(path):
        raise SystemExit(f"missing {path}; run fit_closure.py on {case} first")
    d = np.genfromtxt(path, delimiter=",", names=True)
    r, tau = d["r_mm"], d["tau_w_Pa"]
    order = np.argsort(r)
    return r[order], tau[order]


def info(case: str) -> dict:
    out = {}
    for line in open(os.path.join(HERE, case, "case_info.txt")):
        kv = line.strip().split("=")
        if len(kv) == 2:
            try:
                out[kv[0]] = float(kv[1])
            except ValueError:
                out[kv[0]] = kv[1]
    return out


def root_at(r, tau, level):
    """Radius where tau_w first falls to `level`, on the branch past the peak.

    Returns None if the level is never reached inside the domain, which is exactly
    the situation the narrow case is in and the reason this script exists.
    """
    ipk = int(np.argmax(tau))
    ro, to = r[ipk:], tau[ipk:]
    if to.min() > level:
        return None
    crossings = np.where(np.diff(np.sign(to - level)) != 0)[0]
    if len(crossings) != 1:
        raise SystemExit(f"expected one crossing of {level} Pa, found {len(crossings)}")
    i = crossings[0]
    # linear interpolation in r between the bracketing samples
    t0, t1, r0, r1 = to[i], to[i + 1], ro[i], ro[i + 1]
    return float(r0 + (level - t0) * (r1 - r0) / (t1 - t0))


def main() -> None:
    narrow = sys.argv[1] if len(sys.argv) > 1 else "case_Re10000_H6"
    wide = sys.argv[2] if len(sys.argv) > 2 else "case_Re10000_H6_R50"

    rn, tn = load(narrow)
    rw, tw = load(wide)
    inf_n, inf_w = info(narrow), info(wide)
    print(f"narrow: {narrow}  Rdom={inf_n['Rdom']*1e3:.1f} mm  {len(rn)} samples, "
          f"r to {rn.max():.2f} mm, tau_w(end) = {tn[-1]:.3f} Pa")
    print(f"wide  : {wide}  Rdom={inf_w['Rdom']*1e3:.1f} mm  {len(rw)} samples, "
          f"r to {rw.max():.2f} mm, tau_w(end) = {tw[-1]:.3f} Pa")
    for key in ("Re_d", "H_over_d", "V0"):
        assert abs(inf_n[key] - inf_w[key]) < 1e-9, f"{key} differs between cases"
    print("Re_d, H/d and V0 match between the two cases.\n")

    # ---- 1. validation on the overlap -------------------------------------
    # Compare on the narrow case's own sample radii, interpolating the wide case,
    # over the band the paper actually reads numbers from.
    lo, hi = 2.0 * D * 1e3, rn.max()
    m = (rn >= lo) & (rn <= hi)
    tw_at_rn = np.interp(rn[m], rw, tw)
    rel = np.abs(tw_at_rn - tn[m]) / tn[m]
    print(f"[1] overlap {lo:.1f}-{hi:.2f} mm, {m.sum()} points")
    print(f"    relative difference: median {np.median(rel)*100:.3f}%  "
          f"max {rel.max()*100:.3f}% at r = {rn[m][np.argmax(rel)]:.2f} mm")
    verdict = "AGREE" if rel.max() < 0.05 else "DISAGREE"
    print(f"    -> the wide case {verdict} with the narrow one on the overlap "
          f"(bar: max relative difference < 5%)\n")

    # ---- 2. the cleaning radii, computed rather than extrapolated ----------
    print("[2] cleaning radius, tau_w(r) = tau_c, read off each computed field")
    print(f"    {'tau_c':>6s}  {'narrow (R=25)':>16s}  {'wide (R=50)':>14s}")
    for tc in (1.0, 3.0, 5.0, 10.0):
        a, b = root_at(rn, tn, tc), root_at(rw, tw, tc)
        sa = f"{a:.2f} mm" if a else "not in domain"
        sb = f"{b:.2f} mm" if b else "not in domain"
        print(f"    {tc:6.1f}  {sa:>16s}  {sb:>14s}")

    r1 = root_at(rw, tw, 1.0)
    print()
    if r1 is None:
        print(f"    tau_c = 1 Pa is STILL outside the domain at R = {inf_w['Rdom']*1e3:.0f} mm "
              f"(tau_w falls only to {tw[-1]:.3f} Pa). A wider domain is needed.")
    else:
        print(f"    tau_c = 1 Pa contour COMPUTED at r = {r1:.2f} mm")
        print(f"      paper's full-band extrapolation : {EXTRAP_FULL_BAND_MM:.2f} mm "
              f"({100*(EXTRAP_FULL_BAND_MM/r1 - 1):+.1f}%)")
        print(f"      paper's outer-half extrapolation: {EXTRAP_OUTER_HALF_MM:.2f} mm "
              f"({100*(EXTRAP_OUTER_HALF_MM/r1 - 1):+.1f}%)")
        print(f"      over-provisioning vs r_max = 15 mm: {r1/15.0:.2f}x "
              f"(paper quotes an extrapolated 2.3-2.6x)")


if __name__ == "__main__":
    main()
