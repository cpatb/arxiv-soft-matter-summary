"""Release-parameter distributions spanning the range of real Ultimate throws.

WHY THIS EXISTS
---------------
`trajectories.RELEASE_DISTRIBUTIONS` holds exactly one backhand, one forehand
and one hammer, all at a 15 m/s release. In this model that backhand carries
27.1 +- 4.7 m in still air (measured, see `calibrate_catalogue.py`), which is a
mid-range pass. Sean Earl Field is 78 m long. A study that only ever throws
27 m never puts a disc into the part of the field -- or the part of the height
band -- where the interesting flow is, and it cannot say anything about how
throw LENGTH changes wind sensitivity.

This module adds the short and the long ends of the range. It does NOT touch
`RELEASE_DISTRIBUTIONS`: every production artifact was generated from those and
must stay reproducible.

HOW THE NEW ENTRIES WERE SET, AND WHAT THAT DOES AND DOES NOT LICENCE
---------------------------------------------------------------------
I did not have measured release kinematics for a short pass or a huck, and I am
not going to invent a release speed and present it as observed. So the new
entries are defined by a TARGET STILL-AIR RANGE -- which is a model output that
can be measured and reported -- and the release speed is solved for:

    `calibrate_catalogue.py` bisects `speed_mean` until the median still-air
    range of 400 draws hits the target, holding every other release parameter
    at the published backhand's value except where noted below.

So the CLAIM these entries support is "a throw that carries R metres in this
model", not "a throw released at V m/s as measured on a player". The targets:

    pass_15m     15 m   a give-and-go / dump, thrown flat and low
    backhand_27m 27 m   the existing `right_backhand`, RE-EXPORTED UNCHANGED
    huck_50m     50 m   a long attacking throw
    pull_65m     65 m   a pull, the longest throw in the game

15 / 50 / 65 m are round numbers chosen to span the game, not measurements.
The 27 m entry is a pointer to the published distribution, so the catalogue
contains at least one arm whose release parameters are not of my choosing and
against which the calibrated ones can be read.

`elevation_mean` is raised with range because a 65 m throw released at 8 deg
does not stay airborne long enough to get there at any release speed the model
tolerates; the calibration script reports the elevation used with each solved
speed, and `CALIBRATION` below is written by that script, not by hand.
"""

from __future__ import annotations

import copy

from trajectories import RELEASE_DISTRIBUTIONS

# --- calibrated release speeds ------------------------------------------------
# WRITTEN BY calibrate_catalogue.py -- do not hand-edit. Each entry records the
# target range, the solved release speed, and the still-air range actually
# achieved at that speed, so the calibration can be checked without re-running.
CALIBRATION = {}          # filled in by _load_calibration()

_CATALOGUE_SPEC = {
    # name          target_range_m  elevation_deg  spin_rad_s  base
    "pass_15m":     dict(target_range_m=15.0, elevation_mean=4.0,  spin_mean=85.0,
                         base="right_backhand"),
    "huck_50m":     dict(target_range_m=50.0, elevation_mean=14.0, spin_mean=115.0,
                         base="right_backhand"),
    "pull_65m":     dict(target_range_m=65.0, elevation_mean=22.0, spin_mean=125.0,
                         base="right_backhand"),
}


def _load_calibration():
    """Read the solved release speeds written by calibrate_catalogue.py."""
    import json
    from pathlib import Path
    p = Path(__file__).with_name("throw_catalogue.json")
    if not p.exists():
        return {}
    return json.loads(p.read_text())


def build_catalogue():
    """Return {name: release-distribution dict}, published entries included.

    Raises if a calibrated entry has no solved release speed on disk, rather
    than silently falling back to the published 15 m/s -- a catalogue whose
    `huck_50m` quietly throws 27 m would corrupt every downstream comparison.
    """
    cal = _load_calibration()
    CALIBRATION.clear()
    CALIBRATION.update(cal)

    cat = {}
    # The published distributions, unchanged and under their own names.
    for name, d in RELEASE_DISTRIBUTIONS.items():
        cat[name] = copy.deepcopy(d)
    cat["backhand_27m"] = copy.deepcopy(RELEASE_DISTRIBUTIONS["right_backhand"])

    for name, spec in _CATALOGUE_SPEC.items():
        if name not in cal:
            # Absent, not defaulted. A calibrated entry that quietly fell back
            # to the published 15 m/s would throw 27 m under the name
            # `huck_50m` and corrupt every comparison drawn from it, so it is
            # left out of the catalogue entirely and `get_throw` explains why.
            continue
        d = copy.deepcopy(RELEASE_DISTRIBUTIONS[spec["base"]])
        d["speed_mean"] = float(cal[name]["speed_mean"])
        d["elevation_mean"] = float(spec["elevation_mean"])
        d["spin_mean"] = float(spec["spin_mean"])
        # Release scatter is held at the published backhand's values so that
        # differences between catalogue entries are differences in the THROW,
        # not in how tightly it was modelled.
        cat[name] = d
    return cat


def target_range(name):
    """Nominal still-air range a catalogue entry was calibrated to, or None."""
    spec = _CATALOGUE_SPEC.get(name)
    return None if spec is None else spec["target_range_m"]


def get_throw(name):
    """Return one catalogue entry, with a useful error if it is uncalibrated."""
    cat = build_catalogue()
    if name in cat:
        return cat[name]
    if name in _CATALOGUE_SPEC:
        raise RuntimeError(
            f"throw '{name}' has no calibrated release speed on disk. Run "
            f"`python calibrate_catalogue.py` in P1_AnalyticalModel; it writes "
            f"throw_catalogue.json.")
    raise KeyError(f"unknown throw {name!r}; have {sorted(cat)}")
