# Size, shape, size variation, material properties: what the data on disk can decide

**2026-09-19.** A working session's worth of analytic and geometric argument aimed at one
question: when and how does an impacting dense-suspension drop cross from the ball-like
shear-jammed mode to the inertial-dilatant one, and does particle size decide it?

Numbers here are quoted from the file that printed them, with its path. Nothing is typed
from memory. One of the five planned axes was carried through. The other four were not
started, and saying which is the point of §3.

---

## 1. The one item carried through: does size decide the regime?

`derivations/onset_size_degeneracy_2026-09-19.md`. **Three headlines, three withdrawals.**

> 🔴 **The verdict below is itself withdrawn.** Referee pass 3 found that the file this page
> names as its own companion settles two of the three knobs the negative rests on.
> `impact_stress/PORE_PHYSICS.md` §3 is titled "The drainage length is a modelling choice"
> and then makes the choice for exactly this question: for whether frustrated dilatancy can
> raise the particle pressure in the spreading film, the contact length is correct,
> `L = 0.49 D₀`, a computed Wagner `r_c` rather than an assumption, and the length implied by
> Mobaseri et al.'s Eq. 3. The `ℓ = 10a` row is not a rival convention: `size_threshold.py`
> says of that very line that the group is size-blind by construction and predicts no
> threshold. And §9.1 already ran the capacity-versus-drainage comparison: for the lab silica
> the uncapped `p_dil = 6.2e7 Pa` against a capacity ceiling of 94 Pa, which binds by six
> decades. So the binding ceiling is the capacity one, three decades below the cavitation
> ceiling §1 used, and above it the exponent is +1 for a double-layer contact and 0 for a
> Hertzian one: determined per material, not one of five.
>
> Two smaller failures in the same page. The "reproduced to 0.12 %" check is a tautology,
> because `size_threshold.py:140` computes the 10a row as the D₀ row times `(D₀/10a)²`. And
> the bidisperse "refuted" rests on −0.0379 ± 0.0187, t = 2.03 against t_crit = 4.303, which
> is not resolved. Both verified against the files before acting.

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

## 2. 🔴 A retraction: the "correction" this section proposed is withdrawn

This section said `paper/main.tex`'s 4.47x family spread is an endpoint quantity that falls
to 1.17x locally at CV ≈ 0.2, and that one sentence in the Discussion would close a gap.
**Both halves fail, and nothing in `paper/` should change.**

The spread of the 0 to CV chord slope, at every CV, with the max-min difference's own t at
the design's 2 dof:

| CV | spread | max−min | s.e. | t |
|---|---|---|---|---|
| 0.1 | 2.74 | 0.0242 | 0.0080 | 3.04 |
| **0.2** | **1.17** | **0.0100** | **0.0108** | **0.93** |
| 0.3 | 1.42 | 0.0289 | 0.0075 | 3.84 |
| 0.4 | 2.38 | 0.0655 | 0.0064 | 10.23 |
| 0.5 | 4.47 | 0.0964 | 0.0071 | 13.52 |

1.17x is the **argmin** of its own sweep: I quoted the single CV at which the families are
closest, while §1's own `CV*` sits at 0.30 to 0.35 where the spread is 1.42 and rising. And
at that CV the difference is not resolved (t = 0.93 against t_crit = 4.303), so it
establishes nothing. An unresolved number cannot demote one that passes three gates in
`main.tex` (ANOVA F = 144.7, Welch t = 13.5 against a selection-calibrated 5.87), and
"artefact" was the wrong word for a result the manuscript already brackets as 2.89 to 5.51.

Found by referee pass 3 and verified against the committed fits before retracting.

---

## 3. The five axes: which the data on disk can test, and which cannot

| axis | can the existing data test it? | what it would take |
|---|---|---|
| **absolute size** | **More than §1 claimed, less than a result.** After the §1 retraction: `PORE_PHYSICS` §3 pins the drainage length and §9.1 shows the capacity ceiling binds by six decades for the lab silica, so the exponent there is +1 (double-layer contact) or 0 (Hertzian), per material. What is open is narrower: `φ_c`, which contact law the material has, and `λ_D`. | `φ_c` for this material; a contact-law determination; `λ_D`, which `PORE_PHYSICS` §5.2 says spans a decade and is not measured. Plus `α` in a second chemistry: of the five candidates DOI-resolved this session only Guy et al. 2015 has a particle-size series of `σ*`, and ⚠ the other four were checked at abstract level only, so that is CANNOT TELL rather than a hard negative. |
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

Three headlines withdrawn in one day by three referees. The first two shared a mistake: each
version chose a convention where the literature offers several, and reported the consequence
as a result. The second is instructive because the repair repeated the error it was
repairing.

The third failure is different and worse. Having been caught twice choosing a convention, I
wrote a page whose whole claim was that the conventions are unpinned, in bold, while the
file I named as the page's own companion has a section titled "The drainage length is a
modelling choice" that makes the choice, and a worked section that runs the comparison I
said nobody had run. I cited that file six times and did not read the two sections that
answer the question. Over-claiming a negative is still over-claiming: "nobody has pinned
this" is a literature claim and it needs the same reading a positive one does.

The check script could not catch either, for a reason worth writing down. The symbolic pass
carried the disputed factor as a free symbol, so it could not see the dependence. The
numeric pass finite-differenced only the one knob where that factor is inert. A finite
difference checks the knobs you vary and no others.

The negative that came out of it is a better deliverable than either headline would have
been. It names three measurements, one of which is cheap, and it removes a question from the
project's list rather than adding a claim to it.
