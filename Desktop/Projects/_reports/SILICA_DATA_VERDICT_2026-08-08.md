# Is the silica data still worth using?

Written 2026-08-08, after building the rebound/modulus package
(`E:\CPB\rebound_modulus_2026-08-08\`) and the literature-framing package
(`E:\CPB\literature_framing_2026-08-08\`) over all six sessions.

---

## The short answer

**Keep it. It supports a real paper. But the paper is not the one the six
sessions were designed for, and the publishable core is smaller than six
sessions.**

You are right that the data is not reducible to a comparison with a single
flow curve. That door is closed and nothing in the existing files reopens it.
What is not closed is everything that never needed the flow curve, and that
turns out to be most of what is interesting.

The reframe in one line: **stop trying to put a number on the φ axis, and use
the data for orderings, timescales and mechanism instead.**

---

## 1. What actually dies

These need a known, matched material and cannot be rescued from what is on
disk:

- Any statement of the form "at φ = 0.44, σ\* = X Pa".
- Any comparison of your fitted γ̇_c to a literature or lab flow curve. §8
  already records this: `silica_full_rheo_800_nm_colloidal.csv` is not on your
  φ scale, by an unquantified amount.
- Absolute φ_c, φ_m or φ_0. The brief shows a +0.078 φ-axis offset would
  reconcile your Wyart-Cates fit, and that the most likely single cause is
  ρ_p for porous Stöber silica, on which your own files carry an unresolved
  1.9-2.6 g/cm³ spread. A 15% error in ρ_p moves every φ you have.
- The advisor's item 3 as originally posed, "how does the predicted γ̇_c
  compare to the steady-state flow curve". That moves to TPM, as §2b already
  concluded.

Note this list was already mostly known. Your batch worry does not add a new
category of loss; it removes the last hope of the φ axis being quantitative.

## 2. What survives, and it is more than it looks

Everything below is **ordinal or kinematic**, and none of it contains ρ, φ, η,
σ\* or a maximum packing fraction.

**Within a session, everything is apples to apples.** One session is one
suspension, one day, one syringe. The velocity ladder inside a session is
unimpeachable regardless of what the suspension actually was. That is where
the front-speed vs shear-rate slopes, the arrest ladders and the rebound
series all live.

**The salt axis is much stronger than the φ axis**, and this is the point most
worth internalising. Salt is an *intervention on a fixed stock*: 0, 10 and
100 mM come from the suspension combined on 22JUN26 and were run on 23-24JUN.
Different φ points are different preparations. So a salt contrast asks "what
did adding LiCl do to *this* material", which is answerable without knowing
what the material was. A φ contrast asks "how do these two materials differ",
which is not.

Concretely, these now exist and hold:

| result | where | needs φ? |
|---|---|---|
| t\*_arrest falls monotonically with salt, ρ = −0.82, p < 1e-4, velocity- threshold- and calibration-matched | framing pkg | no |
| Early spreading is universal at 0.64 × Wagner across all six sessions | framing pkg | no |
| Jamming is established before the impact force peaks (t\* < 0.315) | framing pkg | no |
| Four decades of t\* separate jam formation from final relaxation | framing pkg | no |
| Salt makes the drop 1.57× bouncier and shortens contact 2.3× at matched u₀ | rebound pkg | no |
| Elastic wave speed c, density-free by construction | rebound pkg | no |
| u_front = 6.85 u₀, matching the independent 7.06 | rebound pkg | no |

That is a coherent paper: **drop impact resolves the formation and the
un-jamming of a shear-jammed solid in one record, four decades of time apart,
and screening the electrostatic repulsion moves both.** No flow curve appears
in that sentence.

## 3. The 50 mM session is probably a different material, and that is a finding

Your non-monotonic jamming is real, and it localises to one session.

**Four independent observables put 50 mM out of order:**

| observable | 0 mM | 10 mM | 50 mM | 100 mM | verdict |
|---|---|---|---|---|---|
| t\*_arrest, band B (new) | 0.266 | 0.228 | **0.179** | 0.187 | 50 below 100 |
| d/d₀ at jam, front anchor (§6) | 0.751 | 0.667 | **0.544** | 0.612 | 50 below 100 |
| γ̇_c own-fits (§6) | 632 | 479 | **471** | 195 | 10 and 50 tied |
| restitution (new, unmatched u₀) | n/a | 0.120 | **0.169** | 0.132 | 50 above 100 |

**And 50 mM is the only salt point from a different stock date.** 0 and
100 mM are 23JUN26, 10 mM is 24JUN26, both within two days of the 22JUN26
"combined various silica suspensions into a new stock" log entry. 50 mM is
13JUL26, three weeks later, and is also the only salt point in the JUL
campaign (37 px/mm, threshold 120).

So the honest reading is not "the salt effect is non-monotone". It is **"three
salt points from one stock are monotone in every observable, and a fourth
point from a different stock and campaign is not."** Drop 50 mM to a
robustness appendix and the ladder is clean everywhere.

This is the single most useful thing to come out of the batch question: it
converts an embarrassment into a provenance argument, and it costs nothing.

## 4. Do you need new synthesis?

**Not for the paper in §2. Yes for the paper you originally wanted.**

New, characterised particles buy exactly one thing: the right to put numbers
on the φ axis and to compare to a flow curve. If the thesis is "drop impact
resolves jamming and un-jamming spatiotemporally, and screening moves both",
you do not need it. If the thesis is "the impact-derived critical shear rate
equals/differs from the rheometric one by X", you need it and there is no way
around it.

The TPM program in §2b is already the right answer to the second. Nothing here
changes that, and nothing here should delay it.

**What would materially strengthen the existing data, cheaply, without new
synthesis:**

1. **One repeat cell.** Re-run φ = 0.44, 0 mM at band B on a fresh day from the
   current stock. If t\*_arrest and e reproduce, batch drift is bounded and
   every ordinal claim hardens. This is what §3's interleaved repeatability
   session is for, and the new kinematic observables make it cheaper to read
   out than the old ones.
2. **Retain and rheometer the actual stocks you still have.** Even one flow
   curve on the *same vial* that produced a session converts that session from
   ordinal to quantitative. You do not need all six.
3. **A degassed-vs-not pair** (see the rebound package §6). Half a session.

---

## 5. On φ_c in a rheometer versus in a free surface

**Your intuition is correct, it has a literature home, and your own brief
already says the mechanism is required.** Taking those in turn.

**It is correct and it is Brown & Jaeger's.** In a dense suspension the
thickened state must dilate, and dilation has to be resisted by *something*.
Whatever provides the weakest confining stress sets the ceiling on the
thickened stress. In a rheometer that is the instrument: the gap, the normal
force, the plate stiffness. At a free surface it is surface tension, and the
particle-scale scale is σ ~ Γ/a. So σ\* is not a material constant; it is a
material-plus-boundary constant. A cone-plate number is a statement about a
suspension *in a cone-plate*.

**Your brief already needs this mechanism for an unrelated reason.** §"the
Wyart-Cates gap" states that frustrated dilatancy against the free surface, or
confinement, is required *in any case*, because steady-state Wyart-Cates is a
local constitutive law with **no front solution at any φ** and a propagating
solidification front is exactly what you measure. So the free surface is not
an optional extra in your system; without it your central observable has no
mechanism.

**Where the numbers land, and why drop impact is well placed.** Two
surface-tension stresses differ by three orders of magnitude:

| scale | expression | value (800 nm silica, water) |
|---|---|---|
| particle-scale meniscus | Γ/a | ≈ 90 kPa |
| whole-drop Laplace | 2Γ/R | ≈ 120 Pa |
| **impact stress (brief)** | Wagner, on axis / at rim | **7.5 / 42 kPa** |

**The impact stress sits between them.** That is a genuinely favourable
position: drop impact loads the suspension hard enough to matter against the
particle-scale meniscus but not so hard as to overwhelm it, whereas a
rheometer's confining stress is set by hardware and can be put anywhere. It is
a real argument that drop impact probes a branch a rheometer cannot.

### Is it a paper?

**Yes, as a mechanism paper, and the framing is stronger than the one you have
been trying to defend.** "The σ\* from a confined rheometer does not govern a
free-surface flow, and here is the free-surface scale measured directly" turns
your missing flow curve from a hole into the thesis.

**One caution, and it is important.** To claim the rheometric value does not
transfer, you need both measurements *on the same suspension*. Otherwise
"doesn't transfer" is indistinguishable from "different material", which is
exactly the trap the current dataset is already in. So this framing does not
escape the need for matched rheology. It changes what the rheology is *for*:
not an anchor to agree with, but the confined-geometry arm of a two-geometry
comparison.

### The discriminating experiment, and it is cheap

Your brief already names it: **vary surface tension at fixed φ and salt, and
re-measure the arrest threshold.** A nonionic surfactant takes Γ from 72 to
30 mN/m. If the arrest moves with Γ, the free surface is load-bearing and the
mechanism is frustrated dilatancy. If it does not, the mechanism is bulk.
One afternoon, no rheometer, no new synthesis.

Two additions from the new packages:

- **Read it out on t\*_arrest, not on u₀\*.** t\*_arrest is dimensionless,
  needs no φ, has an IQR of 0.005-0.038 within a cell, and orders cleanly. It
  is a far better readout than the fitted onset velocity.
- **A drop-size ladder separates the two surface-tension scales.** Γ/a does
  not depend on drop radius; 2Γ/R does. One nozzle swap at fixed φ and Γ tells
  you which scale is load-bearing. This is a different experiment from the
  surfactant one and answers a different question, and the open-problems slate
  already ranks a drop-size ladder highly for unrelated reasons.

⚠ The surfactant experiment has one confound worth planning around: most
surfactants adsorb onto silica and can change the interparticle interaction as
well as Γ. Pair it with the drop-size ladder, which changes the confining
stress geometrically and touches no chemistry.

---

## Bottom line

- **Keep the data.** The ordinal, kinematic and timescale results are sound and
  they are enough for a paper.
- **Drop 50 mM from the headline.** It is a different stock and four
  observables say so. Three salt points from one stock are monotone in
  everything.
- **Do not spend more effort trying to put numbers on the φ axis** from these
  sessions. That is unrecoverable and it is not needed for the story that is
  actually supported.
- **The free-surface-versus-rheometer idea is your strongest new framing**, it
  is grounded in Brown & Jaeger, your own brief already requires the
  mechanism, and it has a one-afternoon discriminating experiment. It does
  still need one well-characterised suspension measured both ways, so it
  argues *for* the TPM program rather than around it.
