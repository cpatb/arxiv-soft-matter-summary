"""Single definitions for the two quantities that were drifting between files.

Both of these were audit findings (AUDIT_2026-08-29.md, F3 and F4).

1. DISC-ALTITUDE SPEED. The manuscript's headline sheltering number is a fraction:
   the mean 2 m speed over the playing surface, divided by the ABL reference speed at
   10 m. Those two must be the same kind of quantity. `U_ref` is the Richards-Hoxey
   log-law inlet speed and is horizontal *by construction* -- the inlet BC sets v = w = 0
   (`P3_CFD/generateCases.py`, `abl_reference`). The field sampled at 2 m, however, is a
   full 3-vector, and near the resolved buildings it carries a genuine vertical
   component that the reference cannot contain. Taking `norm` over all three components
   and dividing by a horizontal reference compares unlike with unlike, and does so
   exactly where the paper looks: in the building wakes.

   The convention is therefore: **any speed that is normalised by `U_ref`, or reported
   as a fraction of it, is the HORIZONTAL magnitude sqrt(u^2 + v^2).**

   This is a diagnostic convention only. Nothing is thrown away physically: the Monte
   Carlo flies the full 3-component field through `cfd_field.py`, and the disc feels w
   as it always did. Only the reported scalar changes.

   Measured effect on the box means at z = 2 m over the pitch (2026-08-29):

       heading   3-component   horizontal   shift
       090         1.2123        1.2057     +0.55%
       210         1.4195        1.4034     +1.15%
       270         1.3956        1.3882     +0.53%

   so the published "24-28% of the reference" is unchanged either way. Locally it is
   not negligible: |w|/|U_h| reaches 3.5 in a wake cell, with a p95 of 0.23-0.41, so
   the per-cell Q map does carry the difference.

2. THE REFERENCE SPEED ITSELF. `U_ref` was a hardcoded 5.0 in three separate files,
   while the value that actually drove each CFD case is written into that case's
   `meta_json`. Nothing read it back, so the match held only by manual bookkeeping --
   the same shape of defect as the `--exclude-throws hammer` mismatch, where the
   shipped artifact and the committed script disagreed about what had been run.
   `u_ref_from_field` reads it out of the artifact and fails loudly on disagreement.
"""

from __future__ import annotations

import json

import numpy as np


def horizontal_speed(U):
    """sqrt(u^2 + v^2) from a (..., 3) velocity array.

    The convention for every speed that is compared to, or normalised by, U_ref.
    See this module's docstring for why the vertical component is excluded.
    """
    U = np.asarray(U)
    if U.shape[-1] != 3:
        raise ValueError(f"expected a (..., 3) velocity array, got {U.shape}")
    return np.hypot(U[..., 0], U[..., 1])


def u_ref_from_field(F, headings=None, expected=None, rtol=1e-9):
    """Read the ABL reference speed out of the flow field's own metadata.

    `F` is an open ``np.load`` of ``flow_fields.npz``. Every heading present must
    agree, since the manuscript quotes a single reference speed for all of them.

    Raises if the headings disagree with each other, or if `expected` is given and
    does not match -- the point of the function is to make a silent mismatch
    impossible, so it must never fall back to a default.
    """
    if "meta_json" not in F.files:
        raise KeyError(
            "flow_fields.npz carries no meta_json, so u_ref cannot be established "
            "from the artifact. Regenerate the field with the current P3 scripts "
            "rather than passing a hardcoded reference speed.")

    meta = json.loads(str(F["meta_json"]))
    if headings is None:
        headings = [int(h) for h in F["headings"]]

    found = {}
    for h in headings:
        for key in (f"dir_{int(h):03d}", str(int(h)), int(h)):
            if isinstance(meta, dict) and key in meta:
                entry = meta[key]
                if isinstance(entry, dict) and "u_ref" in entry:
                    found[int(h)] = float(entry["u_ref"])
                break

    if not found:
        raise KeyError(
            f"no u_ref found in meta_json for headings {headings}; "
            f"meta_json top-level keys are {list(meta)[:8]}")

    values = set(found.values())
    if len(values) > 1:
        raise ValueError(
            f"the CFD cases were driven at different reference speeds {found}; "
            "a single sheltering fraction cannot be quoted across them")

    u_ref = found.popitem()[1]
    if expected is not None and not np.isclose(u_ref, expected, rtol=rtol):
        raise ValueError(
            f"u_ref mismatch: the field was driven at {u_ref} m/s but this script "
            f"was told {expected} m/s. Fix the caller, do not edit this check.")
    return u_ref


# The playing surface, from the OSM outline (see runDownstream.sh). Defined once here
# because it was previously duplicated as inline literals in make_results_tex.py, as
# named constants in vertical_structure.py, and in rounded form in paper_figures.py.
PITCH_X = (-31.0, 47.0)
PITCH_Y = (-10.0, 15.0)


def pitch_mask(x, y):
    """Boolean index arrays selecting the playing surface from a field grid."""
    fx = (np.asarray(x) >= PITCH_X[0]) & (np.asarray(x) <= PITCH_X[1])
    fy = (np.asarray(y) >= PITCH_Y[0]) & (np.asarray(y) <= PITCH_Y[1])
    return fx, fy


def disc_altitude_index(z, target=2.0):
    """Index of the sampled altitude nearest `target`, with a grid-derived tolerance.

    The previous form hardcoded ``atol=0.51``, which is not tied to the grid: a re-run
    with different z-spacing could silently accept a "2 m" sample half a metre off, in
    the region where the log law is steepest. The tolerance is half the smallest
    spacing, so it can only ever accept the genuinely nearest level.
    """
    z = np.asarray(z, dtype=float)
    zk = int(np.argmin(np.abs(z - target)))
    if z.size > 1:
        atol = 0.5 * float(np.min(np.diff(np.sort(z))))
    else:
        atol = 0.0
    off = abs(float(z[zk]) - target)
    return zk, off, (off <= atol)
