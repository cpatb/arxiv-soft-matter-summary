"""grid_sign_check.py -- does the FINE mesh give higher or lower disc-band speed?

WHY THIS EXISTS
---------------
WP4b item 3 asks for a caveat saying the fine mesh gives *higher* speeds than the coarse,
so the paper's "24-28 % of the reference" is biased LOW and errs toward overstating
shelter. That direction matters -- it decides whether the study's central claim is
conservative or flattering -- and it is **not in the validation record**:
`validate_cfd.test_grid_independence` stores `abs(ua - ub) / ua`, so
`cfd_validation.json` keeps only a magnitude (`worst_dU` = 0.571, i.e. 57.1 % relative,
10.3 % of `u_ref`). The per-station signed values are printed to the console and never
saved.

Rather than repeat the brief's unverified 1.297 / 1.021, this recovers the sign from the
two field files that are on disk, using the paper's own definition of the reported
quantity: the horizontal speed at 2 m averaged over the playing surface.

    python grid_sign_check.py

⚠ SCOPE. This is heading 090 only, because `flow_fields_coarse090.npz` is the only coarse
field stored. A single heading fixes the sign for that heading; it does not establish that
the sign is the same at 210 and 270, and the script says so rather than generalising.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "P5_Synthesis"))
from conventions import (disc_altitude_index, horizontal_speed,   # noqa: E402
                         pitch_mask)
FINE = HERE / "flow_fields.npz"
COARSE = HERE / "flow_fields_coarse090.npz"
Z_TARGET = 2.0


def field_mean_speed(npz, key_u: str, z_target: float):
    """Horizontal speed at the z-plane nearest `z_target`, meaned over the pitch mask."""
    z = npz["z"]
    # Use the CANONICAL altitude picker, not a bare argmin. `disc_altitude_index` carries a
    # grid-derived tolerance (half the smallest z-spacing) and RAISES if the nearest sampled
    # level is not actually within it -- which a bare argmin cannot do, because argmin always
    # returns something. A units audit on 2026-09-03 flagged the bare argmin here as the same
    # shape of defect this file's own comment already records fixing for `pitch_mask`: a
    # convenience copy of a convention is a convention that can drift. Both stored fields
    # happen to sample z = 2.0 exactly, so nothing was wrong -- but only by luck, and a
    # regenerated field with different spacing would have been accepted silently.
    zk, z_off, within = disc_altitude_index(z, target=z_target)
    if not within:
        raise ValueError(
            "nearest sampled altitude is %.4f m, i.e. %.4f m off the %.2f m target, which "
            "is outside the grid-derived tolerance. Refusing to report a '2 m' speed from "
            "a level that is not 2 m." % (float(z[zk]), z_off, z_target))
    # `horizontal_speed` too, for the same reason as `pitch_mask` -- it is what
    # make_results_tex.py forms the Q table's |U_2m| column with.
    mag = horizontal_speed(npz[key_u][:, :, zk, :])
    # IMPORT the canonical mask; do not reproduce it. The first version of this script
    # wrote |x| <= 45, |y| <= 25 "so it does not depend on the synthesis package", and
    # that was simply wrong -- the pitch is PITCH_X = (-31, 47), PITCH_Y = (-10, 15), an
    # asymmetric rectangle less than a third the area. Averaging over the wrong region
    # gave 1.1557 m/s against the paper's own 1.21 for the same heading. A convenience
    # copy of a convention is a convention that can drift, which is the whole reason
    # `conventions.pitch_mask` exists.
    x, y = npz["x"], npz["y"]
    fx, fy = pitch_mask(x, y)
    sub = mag[np.ix_(fx, fy)]
    return float(sub.mean()), float(z[zk]), int(fx.sum()), int(fy.sum())


def main() -> int:
    for p in (FINE, COARSE):
        if not p.exists():
            print("missing:", p)
            return 1
    f = np.load(FINE)
    c = np.load(COARSE)

    fkeys = [k for k in f.files if k.endswith("__U") and "090" in k]
    ckeys = [k for k in c.files if k.endswith("__U")]
    if not fkeys or not ckeys:
        print("no 090 velocity key found. fine:", f.files, "coarse:", c.files)
        return 1

    uf, zf, nx, ny = field_mean_speed(f, fkeys[0], Z_TARGET)
    uc, zc, _, _ = field_mean_speed(c, ckeys[0], Z_TARGET)

    print("heading 090, horizontal speed at z = %.2f m (fine) / %.2f m (coarse),"
          % (zf, zc))
    print("meaned over the %d x %d pitch mask:" % (nx, ny))
    print("  fine mesh   |U| = %.4f m/s   (%s)" % (uf, fkeys[0]))
    print("  coarse mesh |U| = %.4f m/s   (%s)" % (uc, ckeys[0]))
    print("  fine - coarse   = %+.4f m/s   (%+.1f %% of coarse)"
          % (uf - uc, 100 * (uf - uc) / uc))
    print()
    if uf > uc:
        print("  => the FINE mesh is FASTER. Refining raises the disc-band speed, so the")
        print("     reported shelter fraction is biased LOW and the study ERRS TOWARD")
        print("     OVERSTATING shelter. That is the conservative direction for the")
        print("     'deeply sheltered' claim and the unfavourable one for any number")
        print("     quoted as a shelter fraction.")
    else:
        print("  => the FINE mesh is SLOWER. Refining lowers the disc-band speed, so the")
        print("     reported shelter fraction is biased HIGH and the study errs toward")
        print("     UNDERSTATING shelter -- the flattering direction for 'deeply")
        print("     sheltered', and it must be stated as such.")
    print()
    print("⚠ One heading only (090); `flow_fields_coarse090.npz` is the only coarse field")
    print("  stored. This fixes the sign at 090 and says nothing about 210 or 270.")

    out = HERE / "grid_sign_check.json"
    out.write_text(json.dumps({
        "heading": 90, "z_target_m": Z_TARGET, "z_used_fine_m": zf,
        "z_used_coarse_m": zc, "U_fine_ms": uf, "U_coarse_ms": uc,
        "fine_minus_coarse_ms": uf - uc,
        "fine_minus_coarse_pct_of_coarse": 100 * (uf - uc) / uc,
        "fine_is_faster": bool(uf > uc),
        "caveat": "heading 090 only; the only coarse field stored",
    }, indent=2), encoding="utf-8")
    print("\nwrote", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
