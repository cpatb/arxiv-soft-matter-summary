"""
Canonical Ultimate throws with sensible default release parameters.

These return ready-to-integrate initial-condition tuples. Use with
disc_6dof.simulate(). Numbers are practitioner-typical, not extreme.
"""

from __future__ import annotations
import numpy as np
from disc_6dof import build_initial_conditions


def right_handed_backhand(
    launch_pos=(0, 0, 1.5),
    speed=15.0,
    heading_deg=0.0,
    elevation_deg=8.0,
    release_alpha_deg=8.0,    # nose-up release
    roll_deg=-5.0,            # slight hyzer
    spin_rad_s=100.0,
):
    return build_initial_conditions(
        launch_pos, speed, heading_deg, elevation_deg,
        release_alpha_deg, roll_deg, spin_rad_s,
        handedness="right_backhand",
    )


def right_handed_forehand(
    launch_pos=(0, 0, 1.5),
    speed=15.0,
    heading_deg=0.0,
    elevation_deg=6.0,
    release_alpha_deg=6.0,
    roll_deg=+5.0,            # mild anhyzer typical
    spin_rad_s=90.0,
):
    return build_initial_conditions(
        launch_pos, speed, heading_deg, elevation_deg,
        release_alpha_deg, roll_deg, spin_rad_s,
        handedness="right_forehand",
    )


def hammer(
    launch_pos=(0, 0, 2.0),
    speed=11.0,
    heading_deg=0.0,
    elevation_deg=40.0,       # overhead launch
    release_alpha_deg=40.0,   # disc plane near-vertical at release
    roll_deg=+80.0,           # disc tilted almost on its side
    spin_rad_s=110.0,
):
    return build_initial_conditions(
        launch_pos, speed, heading_deg, elevation_deg,
        release_alpha_deg, roll_deg, spin_rad_s,
        handedness="right_forehand",
    )


# Convenient release-parameter distributions for Monte Carlo
RELEASE_DISTRIBUTIONS = {
    "right_backhand": dict(
        speed_mean=15.0, speed_std=1.5,
        elevation_mean=8.0, elevation_std=3.0,
        release_alpha_mean=8.0, release_alpha_std=4.0,
        roll_mean=-5.0, roll_std=5.0,
        spin_mean=100.0, spin_std=10.0,
        handedness="right_backhand",
    ),
    "right_forehand": dict(
        speed_mean=15.0, speed_std=1.5,
        elevation_mean=6.0, elevation_std=3.0,
        release_alpha_mean=6.0, release_alpha_std=4.0,
        roll_mean=5.0, roll_std=5.0,
        spin_mean=90.0, spin_std=10.0,
        handedness="right_forehand",
    ),
    "hammer": dict(
        speed_mean=11.0, speed_std=1.5,
        elevation_mean=40.0, elevation_std=5.0,
        release_alpha_mean=40.0, release_alpha_std=8.0,
        roll_mean=80.0, roll_std=10.0,
        spin_mean=110.0, spin_std=15.0,
        handedness="right_forehand",
    ),
}


def sample_release(
    throw_type: str,
    heading_deg: float,
    launch_pos=(0, 0, 1.5),
    rng: np.random.Generator | None = None,
):
    """Draw a single release-parameter realization from the named distribution."""
    rng = rng or np.random.default_rng()
    d = RELEASE_DISTRIBUTIONS[throw_type]
    return build_initial_conditions(
        launch_pos=launch_pos,
        speed=rng.normal(d["speed_mean"], d["speed_std"]),
        heading_deg=heading_deg,
        elevation_deg=rng.normal(d["elevation_mean"], d["elevation_std"]),
        release_alpha_deg=rng.normal(d["release_alpha_mean"], d["release_alpha_std"]),
        roll_deg=rng.normal(d["roll_mean"], d["roll_std"]),
        spin_rad_s=max(20.0, rng.normal(d["spin_mean"], d["spin_std"])),
        handedness=d["handedness"],
    )
