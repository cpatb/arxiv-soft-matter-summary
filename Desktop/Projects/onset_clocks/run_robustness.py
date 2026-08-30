"""Sweep the withdrawn calibration and report which conclusions survive it.

Run:  python run_robustness.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from data import is_arrested, load_drops  # noqa: E402
from robustness import (  # noqa: E402
    ETA_S_RANGE,
    PHI_0_RANGE,
    PHI_M_RANGE,
    SIGMA_STAR_RANGE,
    Params,
    all_above_ceiling,
    can_any_drop_jam,
    corners,
    f_range,
    marginal_jamming_offset,
)

OUT = HERE / "out"
FITTED = Params(eta_s=4.736e-3, phi_0=0.70, phi_m=0.5176, sigma_star=44.2327)
# The clock is judged flat if its dynamic range cannot span the response's.
RESPONSE_RANGE = 6.35


def main() -> int:
    drops, _ = load_drops()
    arrested = [d for d in drops if is_arrested(d)]
    OUT.mkdir(exist_ok=True)

    print("=" * 96)
    print("ROBUSTNESS TO A WITHDRAWN CALIBRATION")
    print("=" * 96)
    print("\nsi_fit was withdrawn 2026-08-26; E:\rheology was retracted 2026-07-29.")
    print("So the question is not 'what are the parameters' but 'over what range")
    print("of them does each conclusion hold'.\n")
    print(f"  eta_s   {ETA_S_RANGE[0]*1e3:.3f} to {ETA_S_RANGE[1]*1e3:.3f} mPa.s "
          f"(water -> the withdrawn fit)")
    print(f"  phi_0   {PHI_0_RANGE[0]} to {PHI_0_RANGE[1]}")
    print(f"  phi_m   {PHI_M_RANGE[0]} to {PHI_M_RANGE[1]}")
    print(f"  sigma*  {SIGMA_STAR_RANGE[0]:.2f} to {SIGMA_STAR_RANGE[1]:.1f} Pa "
          f"(fitted x0.1 to x10)")

    box = list(corners(n_per_axis=4))
    print(f"\n{len(box)} valid corners of a 4^4 grid "
          f"(corners with phi_0 <= phi_m are skipped: the model thins there).\n")

    res = {"n_corners": len(box), "fitted": FITTED.__dict__}

    # C1
    jammable = [q for q in box if can_any_drop_jam(drops, q)]
    print(f"C1  'no packing in this dataset can shear-jam'")
    print(f"      holds on {len(box) - len(jammable)} of {len(box)} corners "
          f"({100*(1-len(jammable)/len(box)):.0f}%)")
    if jammable:
        pm = [q.phi_m for q in jammable]
        print(f"      it fails only where phi_m <= {max(pm):.4f}, i.e. at or below "
              f"the packings that arrest")
    res["C1_holds_frac"] = 1 - len(jammable) / len(box)

    # C2
    ranges = np.array([f_range(drops, q) for q in box])
    flat = ranges < RESPONSE_RANGE
    print(f"\nC2  'the stress-activation clock cannot span the response' "
          f"(f range < {RESPONSE_RANGE})")
    print(f"      holds on {flat.sum()} of {len(box)} corners "
          f"({100*flat.mean():.0f}%)")
    print(f"      f range over the box: {ranges.min():.2f}x to {ranges.max():.2f}x")
    worst = box[int(np.argmax(ranges))]
    print(f"      widest at sigma* = {worst.sigma_star:.1f} Pa "
          f"(the top of the swept range)")
    res["C2_holds_frac"] = float(flat.mean())
    res["f_range_min"] = float(ranges.min())
    res["f_range_max"] = float(ranges.max())

    # C3
    ok3, mins = zip(*[all_above_ceiling(arrested, q) for q in box])
    ok3 = np.array(ok3)
    print(f"\nC3  'every arrested drop demands more viscosity than f = 1 allows'")
    print(f"      holds on {ok3.sum()} of {len(box)} corners "
          f"({100*ok3.mean():.0f}%)")
    fails = [box[i] for i in range(len(box)) if not ok3[i]]
    if fails:
        print(f"      fails only where phi_m <= {max(q.phi_m for q in fails):.4f}, "
              f"i.e. where the model jams outright and has no finite ceiling")
    res["C3_holds_frac"] = float(ok3.mean())

    # C4
    offs = np.array([marginal_jamming_offset(arrested, q) for q in box])
    offs = offs[np.isfinite(offs)]
    print(f"\nC4  'drops arrest at marginal jamming: |phi_J - phi| stays small'")
    print(f"      worst offset over the box: {offs.max():.4f}  "
          f"(median over corners {np.median(offs):.4f})")
    print(f"      it depends only on eta_s, and only as eta_s^(1/2)")
    res["C4_worst_offset"] = float(offs.max())

    (OUT / "robustness.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"\nwrote {OUT / 'robustness.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
