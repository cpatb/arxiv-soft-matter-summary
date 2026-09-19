# Size, shape, size variation, material properties: what the data on disk can decide

**2026-09-19.** A working session's worth of analytic and geometric argument aimed at one
question: when and how does an impacting dense-suspension drop cross from the ball-like
shear-jammed mode to the inertial-dilatant one, and does particle size decide it?

Numbers here are quoted from the file that printed them, with its path. Nothing is typed
from memory. One of the five planned axes was carried through. The other four were not
started, and saying which is the point of §3.

---

## 1. The one item carried through: does size decide the regime?

`derivations/onset_size_degeneracy_2026-09-19.md`. **Verdict: UNDETERMINED.**

Two headlines were written and both were withdrawn by referee passes.

The first said the two regime criteria carry the same power of particle size, so no size
ladder could separate them. It rested on the claim that no admissible impact-stress gauge
carries a particle size. Three of six do, and the one that breaks it is the Reynolds-Darcy
suction, which this laboratory's own `impact_stress/dilatancy.py` names as the mechanism.

The second said the degenerate pair was the two grain-scale stresses. That failed on a
distinction the first version had blurred. The onset stress `σ*` is the threshold both
routes must cross, not one route's driving stress, so the ratio in question is the dilatant
branch's activation number and says nothing about the frictional route. The discriminator
between the two mechanisms is driving stress against driving stress, and its size exponent
is +2, about a hundredfold per size decade. Both versions asserted the opposite a few lines
below printing that exponent. It is also, exactly, [M26]'s own `ID` times a geometric
factor, so the page had added no discriminator at all.

What the third version reports is a negative with its reason enumerated. The dilatant
stress's size exponent is one of {−2, −1, 0, +1} across five conventions, each sourced to a
file: a span of three decades per size decade. The drainage length enters squared, and
`SIZE_THRESHOLD_2026-09-08.md` §1 already prints both candidate readings without choosing
between them, a fact reproduced here to 0.12 % from the ratio the two readings must have.
On top of that, `dilatant_pressure`'s own docstring instructs the caller to check the result
against a validity ceiling, and neither earlier version did. The Darcy formula is being used
at 1.6x to 459x above that ceiling at every parameter point in play. And at the lab's own
800 nm particle in its own zero-salt chemistry, the double layer filling the pore throat
adds more size dependence than the whole effect under argument.

**What would change the verdict** is three measurements, none of them algebra: whether the
suction relaxes over the drop scale or the grain scale, whether the pore fluid cavitates
during impact, and what `φ_c` is for this material. The first is worth the most, because it
alone is four decades of `ID` and it needs no new simulation.

**What survives** is a sign. The margin `Δφ = φ − φ_c` enters the dilatant stress linearly
and `ID` as its reciprocal, and widening the size distribution at fixed `φ` raises `φ_c` and
so pushes the drop toward the inertial side. That is the direction the author reported from
the bench on 2026-08-31. Nothing quantitative about it survives. Read locally from each
family's own ladder rather than across the endpoints, the three families' slopes differ by
1.17x rather than 4.47x; the bidisperse family's local slope above CV = 0.3 is negative; and
the one `φ_c` on record that names the model's own quantity is flagged a placeholder.

---

---

## 1b. The item that did work, found in the last hour

`derivations/wc_endpoints_move_together_2026-09-19.md`. **Verdict: SUPPORTED as a bound.**

The brief's material-properties item asked what polydispersity does to each Wyart-Cates
endpoint. Both endpoints turn out to be on disk already, at matched CV, on the full six-rung
ladder with three seeds, for the lognormal and uniform families. Their CV dependences are
not resolved apart: lognormal +0.05803 ± 0.00268 against +0.06207 ± 0.00255 (|t| = 1.09),
uniform +0.02385 ± 0.00162 against +0.02420 ± 0.00080 (|t| = 0.20), both against t_crit =
4.303 at 2 dof. The window the model interpolates across, `W = φ_0 − φ_m`, is flat over the
whole ladder: range 0.00646 (1.75 s.e.) for lognormal, 0.00191 (1.06 s.e.) for uniform.

So polydispersity translates the Wyart-Cates φ-axis rather than deforming it. Both ends move
up together, the gap between them does not change, and `σ*` is untouched because it is a
contact-force property. The prediction that follows is clean and needs none of the contested
conventions of §1: widening the size distribution at fixed `φ` moves the suspension down in
`φ/φ_m` without moving the stress at which thickening begins, and the packing fraction at
which DST appears rises by 0.124 per unit CV for a lognormal family, 0.048 for a uniform
one. Same sign as the bench observation, reached by a different route.

⚠ Three things keep this a bound rather than a result. It was **not pre-registered**, which
the brief asked for. At three seeds the design could only have resolved a difference larger
than 26-32 % of the shift, so "move together" means "to within a quarter to a third", never
"identically". And the family that is missing is bidisperse, the non-monotone one, which has
no frictionless block: the family least like the others is the one not tested.

---

## 2. A correction the manuscript should carry, and it is not an error

`paper/main.tex` quotes the family spread in `φ_J` as 4.47x and brackets it correctly
(2.89x peak-defined, 5.51x inverse-variance). Read as a statement about the CV 0 to 0.5
endpoint, that is right. Read as "the families differ by 4.5x at a given CV", it is not.
Locally at CV ≈ 0.2, where the regime question sits, the spread is 1.17x, computed from the
same committed fits. The three families are nearly indistinguishable exactly where the
impact question asks about them, and separate only toward the top of the ladder. One
sentence in the Discussion would close a gap a referee will otherwise find.

---

## 3. The five axes: which the data on disk can test, and which cannot

| axis | can the existing data test it? | what it would take |
|---|---|---|
| **absolute size** | No. §1: the dilatant stress's size exponent is unpinned across 3 decades per size decade, and the one anchor (`α = 1.9(1)`) is measured on a sterically stabilised system where salt does not move `σ*`, while the lab's is charge stabilised. | the three measurements in §1, plus `α` in a second chemistry. Of the five candidates searched and DOI-resolved this session, only Guy et al. 2015 has a particle-size series of `σ*`. ⚠ The other four were checked at abstract level only, so that is CANNOT TELL rather than a hard negative. |
| **size variation (CV)** | Partly. `φ_J(CV)` is measured for three families on a committed ladder, and the sign of its effect on the regime follows. The magnitude does not, because it needs `φ_c`. | `φ_c` for the lab silica, and the contact-weighted PSD moment (§4), which is free. |
| **shape** | Partly, and not started this session. The 83-cell clump campaign measures `φ_J(α)` and `Z_nr(α)`, and the constraint-counting literature is now API-resolved: van Hecke 2010 for `d+1 ≤ Z ≤ 2d`, Donev 2004 for the near-sphere limit (spheres `Z ≈ 6.4`; oblate 9.82 measured against an isostatic 10; and their own observation that the contact number is *singular at the sphere value*, rising sharply for small deviations while the degree-of-freedom count jumps discontinuously). That last point is the one a two-rung α ladder is worst placed to see. | the LevelSet control, already item 1 of the 09-16 list. Until it exists, every α > 1 number is "elongation + interlocking + contact multiplicity". |
| **material properties (friction, stiffness)** | **Yes, and it was done: see §1b.** Both Wyart-Cates endpoints are on disk at matched CV for two of three families, and their CV dependences are not resolved apart. Stiffness was not touched. | more seeds (the bound is 26-32 % of the effect at n = 3), a frictionless bidisperse block (six cells), and the `μ_p = 1.0` ladder to check whether 0.3 stands in for the frictional limit. |
| **chemical potential (salt, evaporation, osmotic pressure)** | Not started, and §1 found the axis is not clean. Salt moves the permeability through `λ_D/r_throat` (0.484 at 800 nm) and also moves `σ*`. A salt ladder is therefore not a discriminator between the mechanisms; what it does cleanly is remove the extra size exponent the double layer contributes. | the derivation is writable with symbolic parameters today. The Paper 1 numbers it would need are behind the raw-data hook. |

---

## 4. One free measurement this project can make today

`κ` reads the Sauter diameter by construction: Carman-Kozeny's length is defined by the
specific surface. Which moment of the size distribution `σ*` reads is a modelling choice
nobody has pinned, and it matters, because the two candidates differ by the full
`d₃₂/d₁₀ = (1+CV²)²` along a CV ladder.

The onset condition is per contact, so the object that answers it is the contact-weighted
moment of the PSD, and the DEM contact lists are already on disk. This replaces a guess with
a measurement and costs no cells. It is item 3 of the Status card's next-action list.

⚠ It also runs into an older open question: whether the lab's "800 nm" is a number mean or a
`d₄₃`. If it is a `d₄₃`, `d₃₂` falls as CV rises and the moment channel reverses sign. That
is Q1 on the Paper 1 card, and it remains the single measurement that decides the most here.

---

## 5. What this session's failures are worth

Two headlines withdrawn in one day, by two referees, on the same underlying mistake. Each
version chose a convention where the literature offers several, and reported the consequence
as a result. The second time is the instructive one, because the repair repeated the error
it was repairing.

The check script could not catch either, for a reason worth writing down. The symbolic pass
carried the disputed factor as a free symbol, so it could not see the dependence. The
numeric pass finite-differenced only the one knob where that factor is inert. A finite
difference checks the knobs you vary and no others.

The negative that came out of it is a better deliverable than either headline would have
been. It names three measurements, one of which is cheap, and it removes a question from the
project's list rather than adding a claim to it.
