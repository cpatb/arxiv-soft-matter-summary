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
