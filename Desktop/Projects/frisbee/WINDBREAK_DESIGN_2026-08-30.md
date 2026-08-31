# Can a barrier fix the turbulent end of Sean Earl Field?

**2026-08-30. Code in `SeanEarl_DiscFlight/P6_Theory` and `P3_CFD`. Gates
pre-registered in `windbreak_gates.py`, committed before any design finished
solving (`2ac4f47`).**

---

## 1. The zone is real, and it is a westerly phenomenon

Disc-band `sigma_u` at z = 2 m, y-averaged over the playing surface, along the
pitch:

| x (m) | −30 | −18 | −6 | +6 | +18 | +30 | +42 |
|---|---|---|---|---|---|---|---|
| **270° westerly** | 0.39 | 0.41 | 0.40 | 0.48 | 0.52 | 0.71 | **0.77** |
| 90° easterly | 0.30 | 0.27 | 0.28 | 0.32 | 0.33 | 0.34 | 0.30 |
| 210° SSW | 0.66 | 0.54 | 0.50 | 0.37 | 0.34 | 0.61 | 0.61 |

Under a **westerly** it rises monotonically west to east by a factor **2.03**,
and the rise is concentrated past x ≈ +18. Campus Towers' footprint starts at
x = +56 on the pitch's y-band, so the eastern third of the field sits in the
windward stagnation and horseshoe-vortex system of a 45 m block, spilling
roughly 16 m back onto the playing surface.

The other two headings do not do this: east/west ratio 1.17 at 90° (the tower is
upwind, so the field is in its wake, which is calmer at 2 m) and a
minimum-in-the-middle shape at 210°.

`P6_Theory/figs_windbreak/fig_zone_270.png` shows it in plan and in a centreline
section. The section is the informative one: the disc band (1–6 m) is dark and
quiet over the western and middle pitch and brightens steeply toward the tower
face, while a much stronger layer sits above z ≈ 10 m across the whole site.

---

## 2. What would count as a fix

The trap is specific and this project already named it.
`FRISBEE_GENERALIZATION.md`'s central claim is that buildings cut mean wind
roughly in half while raising turbulence intensity from 12 % to 25–38 % — and a
disc is disturbed by the fluctuation, not by the mean, which a thrower can aim
off. **A barrier that simply removes wind would score well on any mean-wind
criterion and could still make the field worse to throw on.** A solid fence is
the obvious way to fail like that: it trades a stagnation zone for a shear layer
sitting in exactly the 1–6 m band the disc flies in.

So barriers are modelled as **porous** DarcyForchheimer cell zones — the
standard shelterbelt representation — with the resistance a swept parameter:

```
pressure loss    dp = 0.5 rho k_r U^2
Forchheimer      dp/L = 0.5 rho f |U| U        so   f = k_r / t
```

and `k_r ≈ (1 − β²)/β²` for a thin screen of optical porosity β, so β = 0.5
gives k_r = 3.0 and β = 0.6 gives 1.78. That mapping is a standard orifice-plate
estimate, quoted as such, not measured here.

### The gates, fixed before any result existed

| gate | what | rule |
|---|---|---|
| **G-W1** | mean σ_u over the east third, disc band | ≤ **0.75×** baseline |
| **G-W2** | mean σ_u over *each* third, disc band | ≤ 1.10× baseline everywhere |
| **G-W3** | pitch-mean horizontal \|U\|, disc band | ≥ 0.50× baseline |
| **G-W4** | pitch-max \|∂σ_u/∂z\| inside the disc band | ≤ 1.50× baseline |
| **G-W5** | flown landing dispersion, paired | < 1.00× baseline |
| **G-W6** | G-W2 and G-W3 at 90° and 210° | both pass at both |

G-W4 exists specifically to catch a barrier that swaps a stagnation zone for a
shear layer. G-W5 requires discs to actually fly better, measured by flying them
from the two launch points in and just upwind of the affected zone — a barrier
evaluated on throws from the west end, which land before they reach it, would
score a null it did not earn.

All at **production mesh resolution**, not the cheap grid: the grid study found
the coarse meshes span only 3–6 cells across Campus Towers' 24 m narrow side
against a COST 732 minimum of 10, and this study is entirely about that
building's near field.

### Three architectures

| tag | idea | x | H | k_r |
|---|---|---|---|---|
| **W1** | shelterbelt — shelter the zone from upwind | +21 m | 5 m | 3.00 |
| **W2** | tower deflector — break the stagnation before it forms | +52 m | 10 m | 3.00 |
| **W3** | low open spoiler — lift the shear layer above the disc band | +30 m | 3 m | 1.78 |

---

## 3. W1, the classic shelterbelt: 4 of 5 gates, and it fails the one that matters

| gate | verdict | measured |
|---|---|---|
| G-W1 | **FAIL** | east-third σ_u 0.748 → 0.720 m/s, **×0.962** against a ×0.75 target |
| G-W2 | PASS | worst third is the **west**, at ×1.080 |
| G-W3 | PASS | pitch \|U\| 1.304 → 1.234 m/s, ×0.947 |
| G-W4 | PASS | max \|∂σ_u/∂z\| 0.1437 → 0.1510 s⁻¹, ×1.051 |
| G-W5 | PASS | flown σ 8.34 → 8.07 m and 5.70 → 5.64 m, ×0.978 |

A 3.8 % reduction where 25 % was asked for.

**The failure is informative, and the profile says why.** In
`figs_windbreak/fig_profile_270.png` the fence puts a local notch in σ_u right
at x = +21, the flow recovers immediately downstream, and by x = +40 the two
curves have rejoined. Meanwhile σ_u is **raised** upstream — west third 0.438 →
0.473, and the middle bulges from 0.53 to 0.59.

That is the project's own thesis reproduced in a controlled experiment: **the
barrier removed mean wind and added turbulence, and it did so where the field
had been quietest.** G-W2 caught the 8 % rise, just inside its 10 % bar.

The reason a shelterbelt cannot work here is geometric. **A shelterbelt shelters
from upwind. The turbulence in the east third is not advected in from the west —
it is generated locally by Campus Towers, downwind of the fence.** Putting a
screen at x = +21 shields the zone from air that was not the problem, and adds a
new disturbance over the part of the pitch that was fine.

W2 attacks the source instead of shielding from it. That is the discriminating
comparison, and it is the reason both were in the sweep from the start rather
than W2 being added after W1 disappointed.

---

## 4. All three architectures land in the same place

| design | position | H | H/H_tower | k_r | **G-W1** | worst third (G-W2) | shear (G-W4) | flown (G-W5) |
|---|---|---|---|---|---|---|---|---|
| W1 shelterbelt | +21 m | 5 m | 0.11 | 3.00 | ×0.962 | **×1.080** | ×1.051 | ×0.978 |
| W2 tower deflector | +52 m | 10 m | 0.22 | 3.00 | ×0.969 | ×0.996 | ×0.965 | ×0.978 |
| W3 low open spoiler | +30 m | 3 m | 0.07 | 1.78 | ×0.954 | ×1.066 | ×1.040 | ×0.999 |
| | | | | *target* | **×0.75** | ≤ ×1.10 | ≤ ×1.50 | < ×1.00 |

**Every one fails G-W1, and they all fail it by the same amount.** A 3–5 %
reduction where 25 % was asked for. Position (+21, +30, +52 m), height (3, 5,
10 m) and porosity (k_r 1.78 vs 3.00) barely move the answer: the spread between
the three designs is smaller than the distance any of them is from the goal.

The four-curve profile figure shows the mechanism. Each design puts a **local
notch** in σ_u at its own fence, the flow recovers within about 15 m, and every
curve has rejoined the baseline before x = +35 where the peak actually is. Even
W2, sitting 4 m off the tower's windward face, does not reduce the peak.

W1 and W3 additionally **raise** σ_u over the west and middle thirds (+8.0 % and
+6.6 %) — the project's own thesis reproduced in a controlled experiment: the
barrier removed mean wind and added turbulence, over the part of the pitch that
had been quietest. W2 is the only one that harms nothing, and it is also the
only one placed where the flow is already disturbed.

## 5. The height hypothesis, pre-registered and then falsified

All three failing at once is a pattern, so it got a hypothesis rather than a
fourth guess. Measured on the baseline field: disc-band σ_u begins rising 10 %
above the pitch minimum at x = 0 and the tower face is at x = +56, so the
**upstream influence is 56 m = 1.25 H** against H = 44.8 m — the textbook value
for a bluff body. The three barriers are 0.07, 0.11 and 0.22 H. Each intercepts
at most the bottom fifth of a disturbance whose length scale is set by a 45 m
block.

Prediction: a barrier must be a substantial fraction of H to matter. **W4** put
the same screen as W2 at **25 m = 0.56 H**. The reasoning and the prediction were
committed before the case solved.

| | H | H/H_tower | east-third σ_u |
|---|---|---|---|
| baseline | — | — | 0.748 m/s |
| W2 | 10 m | 0.22 | 0.725 |
| **W4** | **25 m** | **0.56** | **0.722** |

**The prediction is wrong.** Tripling the barrier height moved east-third σ_u by
0.4 %, and W4's profile lies on top of W2's across the entire pitch.

The revised explanation is more restrictive than the one it replaces, and it
follows from what a porous screen can and cannot do. The disturbance over the
pitch is the **upstream influence of the tower's own blockage**: a 45 m solid
block spanning the field width decelerates and lifts the approaching flow, and
that is a far-field consequence of total obstruction. A porous screen placed in
front of it — at any height — barely changes the total blockage the flow sees,
so it cannot undo the field that blockage creates. A barrier placed further
upwind (W1, W3) is not in the causal path at all; it only adds its own local
disturbance to a part of the pitch that was fine.

Being wrong here was more informative than being right would have been. Had W4
worked, the answer would have been "build a taller screen". Because it did not,
the answer is that the mechanism is not interception.

## 6. The last knob: resistance

One thing is still untested. Every barrier so far is porous at k_r = 1.78–3.00
(β ≈ 0.5–0.6). **W5** is near-solid, k_r = 20 (β ≈ 0.22), at the east-third
boundary where a barrier has the best geometric chance of shielding the target.

This is deliberately the design the shelterbelt literature and this project's
own headline both say should fail, and G-W4 exists to catch it. Running it is
not an attempt to make something work; it closes the last door, so that "no
barrier fixes this" is a tested statement rather than three porous designs and
an assumption.

*(W5 result to follow; `P6_Theory/out/windbreak.json` carries every gate value.)*

---

## 5. Limitations that travel with any verdict here

- **The porous zone is 4 m thick**, two production cells, not a physically thin
  screen. The resistance is set from the thickness so a thicker zone with a
  proportionally smaller `f` is the same barrier, but the near-field detail
  inside a fence thinner than the mesh is not represented.
- **Steady RANS.** The quantity being optimised is `sigma_u` inferred from a
  modelled `k`, not a resolved fluctuation. A barrier that changes the character
  of the unsteadiness rather than its energy would not show up.
- **Neutral stratification only.**
- **One geometry per architecture.** Height, distance and porosity were each
  fixed at one value per design rather than swept, so this is a comparison of
  three architectures and not an optimisation within any of them. The gate
  values in `out/windbreak.json` are the starting point for that sweep.
- **G-W6 is not yet run.** Whichever design passes must be re-solved at 90° and
  210°; a barrier that fixes the westerly and ruins the other two is not a fix,
  and nothing here has tested that.
