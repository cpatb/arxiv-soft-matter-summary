"""Acceptance gates for the windbreak design study.

WRITTEN BEFORE ANY WINDBREAK CASE HAD FINISHED SOLVING. The three designs were
generated and submitted first; this file was committed while they were still in
`simpleFoam`, so the thresholds could not have been chosen to fit the answer.
Any later change to a threshold must be recorded as a change, with its reason,
not edited in silence.

THE PROBLEM
-----------
Under a westerly (`theta_wind` = 270) the disc-band turbulence over the pitch
is strongly non-uniform. Measured on the production field at z = 2 m,
y-averaged over the playing surface:

    sigma_u = 0.39 m/s at x = -30      (west end)
    sigma_u = 0.77 m/s at x = +40      (east end)      ratio 2.03

Campus Towers' footprint starts at x = +56 on the pitch's y-band, so the
eastern third of the field sits in the windward stagnation and horseshoe-vortex
system of a 45 m block. The middle of the pitch is its quietest part. The other
two headings do not do this: east/west ratio 1.17 at 90 deg and 1.03 at 210.

WHAT WOULD COUNT AS A FIX, AND WHAT WOULD COUNT AS CHEATING
-----------------------------------------------------------
The trap is specific and this project has already named it. A barrier that
simply removes wind would score well on any mean-wind criterion and could still
make the field WORSE to throw on, because `FRISBEE_GENERALIZATION.md`'s central
claim is that buildings cut mean wind roughly in half while RAISING turbulence
intensity from 12% to 25-38% -- and a disc is disturbed by the fluctuation, not
by the mean, which a thrower can aim off. A solid fence is the obvious way to
fail this way: it would trade a stagnation zone for a shear layer sitting in
the 1-6 m band the disc flies in.

So the gates below reward reducing sigma_u, penalise raising it anywhere,
penalise killing the mean wind, penalise putting a new shear layer in the disc
band, and finally require that discs actually fly better -- measured by flying
them, not inferred from the field.

    G-W1  THE GOAL.        Mean sigma_u over the EAST THIRD of the pitch
                           (x in [21, 47]) within the disc band z in [1, 6]
                           falls by >= 25% against the baseline, at 270 deg.

    G-W2  DO NO HARM, SPACE. Mean sigma_u over NO third of the pitch rises by
                           more than 10% against the baseline.

    G-W3  DO NO HARM, MEAN. Pitch-mean horizontal |U| in the disc band stays
                           at or above 0.5x the baseline. A dead-calm pitch is
                           a different failure, not a success.

    G-W4  NO NEW SHEAR.    The pitch-maximum of |d sigma_u / dz| inside
                           z in [1, 6] does not exceed 1.5x the baseline. This
                           is the gate a solid fence is expected to fail.

    G-W5  IT MUST FLY.     The measured landing dispersion of the standard
                           throw, released from the affected part of the field,
                           improves. Flown through the new field with the same
                           seeds, paired against the baseline.

    G-W6  OTHER HEADINGS.  Whichever design passes G-W1..G-W5 is re-solved at
                           90 and 210 deg and must not fail G-W2 or G-W3 there.
                           A barrier that fixes one wind direction and ruins
                           the other two is not a fix. This gate needs two more
                           CFD runs and is scored separately.

A design passes only if G-W1 through G-W5 all pass. Partial passes are reported
as partial, not rounded up.
"""

from __future__ import annotations

PITCH_X = (-31.0, 47.0)
PITCH_Y = (-10.0, 15.0)
DISC_BAND = (1.0, 6.0)
EAST_THIRD = (21.0, 47.0)
THIRDS = {"west": (-31.0, -5.0), "middle": (-5.0, 21.0), "east": (21.0, 47.0)}

GATES = {
    "G-W1": dict(what="mean sigma_u over the east third, disc band",
                 rule="<= 0.75 x baseline", threshold=0.75, direction="below"),
    "G-W2": dict(what="mean sigma_u over each third, disc band",
                 rule="<= 1.10 x baseline in every third", threshold=1.10,
                 direction="below"),
    "G-W3": dict(what="pitch-mean horizontal |U|, disc band",
                 rule=">= 0.50 x baseline", threshold=0.50, direction="above"),
    "G-W4": dict(what="pitch-max |d sigma_u/dz| inside the disc band",
                 rule="<= 1.50 x baseline", threshold=1.50, direction="below"),
    "G-W5": dict(what="paired landing dispersion of the standard throw",
                 rule="< 1.00 x baseline", threshold=1.00, direction="below"),
    "G-W6": dict(what="G-W2 and G-W3 at theta_wind 90 and 210",
                 rule="both pass at both headings", threshold=None,
                 direction=None),
}

DESIGNS = {
    "W1": dict(name="shelterbelt", x=21.0, height=5.0, k_r=3.00,
               idea="a porous fence on the east-third boundary, sheltering the "
                    "affected zone from upwind in the classic shelterbelt way"),
    "W2": dict(name="tower deflector", x=52.0, height=10.0, k_r=3.00,
               idea="a porous screen 4 m in front of Campus Towers' windward "
                    "face, breaking the stagnation up before it forms rather "
                    "than shielding the field from it"),
    "W3": dict(name="low open spoiler", x=30.0, height=3.0, k_r=1.78,
               idea="a low, more open barrier meant to lift the shear layer "
                    "above the disc band instead of blocking the flow"),
}


# ---------------------------------------------------------------------------
# W4, added after W1-W3 had all failed G-W1, to test WHY rather than to keep
# guessing. Not a fourth guess: a measurement-driven hypothesis test.
#
# The upstream influence of Campus Towers, measured on the baseline field, is
# 56 m -- the disc-band sigma_u begins rising 10% above the pitch minimum at
# x = 0 and the tower face is at x = +56. Against H = 44.8 m that is **1.25 H**,
# which is the textbook upstream influence length of a bluff body. The three
# barriers tried are 3, 5 and 10 m tall, i.e. **0.07, 0.11 and 0.22 H**: each
# intercepts only the bottom fifth or less of a disturbance whose length scale
# is set by a 45 m block.
#
# So the hypothesis is that a barrier must be a substantial fraction of H to
# move this zone at all, and W4 tests it at 25 m = 0.56 H, in the same place as
# W2. If W4 also fails G-W1, the conclusion is that the zone is not fixable by
# a barrier at any plausible scale and the honest recommendation is to use the
# field differently. If it passes, the design question becomes an architectural
# one about how tall a screen anybody would build.
#
# The prediction is recorded before the case has solved.
DESIGNS["W4"] = dict(name="tall tower screen", x=52.0, height=25.0, k_r=3.00,
                     idea="the same position as W2 at 0.56 H instead of 0.22 H, "
                          "testing whether the barrier must scale with the "
                          "building rather than with the disc band")
UPSTREAM_INFLUENCE = dict(onset_x_m=0.0, tower_face_x_m=56.0,
                          length_m=56.0, tower_height_m=44.8,
                          length_over_H=1.25,
                          barrier_H_ratio={"W3": 0.07, "W1": 0.11,
                                           "W2": 0.22, "W4": 0.56})


# ---------------------------------------------------------------------------
# W4 FALSIFIED THE HEIGHT HYPOTHESIS. At 25 m (0.56 H) it gives east-third
# sigma_u = 0.722 m/s against W2's 0.725 at 10 m (0.22 H) -- tripling the
# barrier height changed nothing, and the two profiles lie on top of each other
# across the whole pitch. The prediction recorded above was wrong.
#
# W5 tests the one remaining knob: RESISTANCE. Every barrier so far has been
# porous at k_r = 1.78-3.0 (optical porosity ~0.5-0.6). W5 is near-solid,
# k_r = 20 (beta ~ 0.22), at the east-third boundary where a barrier has the
# best geometric chance of shielding the target.
#
# This is deliberately the design the shelterbelt literature and this project's
# own headline both say should FAIL, and G-W4 exists to catch it: a solid fence
# should trade the stagnation zone for a shear layer inside the 1-6 m disc band.
# Running it is not an attempt to make something work; it is closing the last
# door, so that "no barrier fixes this" is a tested statement rather than three
# porous designs and an assumption.
#
# Two outcomes, both useful. If W5 also fails G-W1, the zone is not barrier-
# fixable at any resistance or height and the honest recommendation is to use
# the field differently. If W5 passes G-W1 but fails G-W4, that is the
# literature's own trade-off measured on this site.
DESIGNS["W5"] = dict(name="near-solid fence", x=21.0, height=6.0, k_r=20.00,
                     idea="the last knob: resistance rather than height or "
                          "position. beta ~ 0.22, the design the shelterbelt "
                          "literature says should swap a stagnation zone for a "
                          "shear layer in the disc band")
