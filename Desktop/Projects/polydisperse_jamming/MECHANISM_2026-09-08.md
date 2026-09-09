# The dilatancy-sign mechanism — what it is, what it replaces, and how to kill it

**2026-09-08.** Derivation and checks: `derivations/dilatancy_sign_mechanism.md`
(17 checks, `ALL CHECKS PASSED`). Measurement: `phase8_criticalstate.py` →
`analyse_phase8.py`. Predictions P1–P7 were committed (`25f1517`) **before** the sweep
returned a single `phi_c`.

---

## 1. The claim, in one paragraph

Mean particle size and size *width* act on an impacting dense suspension through **two
different terms of the same expression, and only one of them is the one the project has
been sweeping**. Permeability — the `d₃₂²` channel — sets *how much* pore pressure the
skeleton's volume change generates. The critical-state fraction `φ_c` — a packing
geometry that polydispersity moves — sets *which way it points*. Below `φ_c` a sheared
pack contracts, the pore fluid goes into **overpressure**, effective stress **falls**,
and the drop liquefies and spreads. Above it the pack dilates, the fluid goes into
**suction**, effective stress **rises**, and the drop arrests. Raising CV at fixed `φ`
raises `φ_c`, so **there is a critical width CV\* at which the drop crosses from one
branch to the other** — and no smooth permeability factor can produce a crossing.

That is a threshold, not a factor, which is why it survives an argument that killed the
factor.

---

## 2. What it replaces

The Status card's mechanism table has three candidates. This changes the standing of
two of them and it is worth being precise about how.

| | as it stood | after this derivation |
|---|---|---|
| **M1** permeability | *"DEAD over most of the range"* — `ID` moves 2.44× over CV 0→0.5 against a boundary 3–4 decades away | **The exclusion swept one of the two CV-dependent factors in `ID` and held the other fixed.** `Δφ = φ − φ_c` was pinned at 0.02 while only `κ` moved. `d ln ID/dCV` has two terms (derivation Eq. 11); the permeability one is `O(CV)` and **vanishes at CV → 0**, the critical-state one is `O(1)` and **has a pole**. M1 remains the weaker channel — that part of the old verdict stands — but the *reason* it is weak is new, and the thing that is strong was never computed. |
| **M2** distance to jamming | *"real and insufficient"* — `φ_J` rises +0.062 but at `φ` = 0.44 that buys only 0.65–0.77× in viscosity | **Graded against the wrong observable.** `η ∝ (φ_J − φ)⁻²` is a smooth factor and 0.7× is indeed nothing. But the same packing shift enters the pore-pressure problem as a **sign**, through `φ − φ_c`. A shift of +0.062 against a margin of 0.02 does not multiply the answer by 0.7, it changes the answer's sign. **M2 is revived as M2′ = the dilatancy-sign mechanism**, and the number that was called insufficient is the number that does the work. ⚠ **Read §7a before quoting this row.** For *its own* hypothesis the old evaluation used the right observable; what happened is that a **second** mechanism was found that consumes the same input number, and §7a says plainly that this row revives **both** and that which one acts is open. |
| **M3** force-network structure | *"survives; polydispersity REMOVES grains, it does not restructure"* (rattlers 0.040 → 0.234) | **Untouched, and now has a role rather than being the residual.** Rattlers are why `φ_c` moves: grains that carry no force do not resist the volume change, so a broader pack reaches its critical state at a higher `φ`. M3 is the microscopic account of M2′'s macroscopic parameter. The constant-volume runs measure both on the same cells. |

**The methodological point, and it is one this repo has made before.** A candidate was
excluded by sweeping the factor that was easy to sweep. `κ(CV)` is a closed-form
identity — `d₃₂/d₁₀ = (1+CV²)²` — and `φ_c(CV)` is a DEM campaign, so the arithmetic got
done and the measurement did not. The exclusion was not wrong about `κ`; it was silent
about the term next to it.

---

## 3. The regime map

At fixed `φ`, two axes and three regimes (derivation §3):

```
                    |  φ < φ_c(CV)  CONTRACTIVE     |  φ > φ_c(CV)  DILATANT
  ------------------+-------------------------------+--------------------------------
   Π ≫ 1 UNDRAINED  |  R2  LIQUEFIED                |  R1  SOLID-LIKE / JAMMED
   (small d₃₂)      |  pore OVERpressure, σ' falls  |  pore SUCTION, σ' rises
                    |  spreads; LOOKS inertial      |  arrests; can rebound
  ------------------+-------------------------------+--------------------------------
   Π ≪ 1 DRAINED    |  R3  RE-ENTRANT LIQUID, both columns: the fluid keeps up and
   (large d₃₂)      |  there is no pore-pressure feedback either way ([M26] Regime 3)
```

- **R1** is the lab's low-CV silica.
- **R2** is what raising CV alone should produce, at fixed chemistry, fixed `φ`, fixed
  mean size. **This is the author's bench observation.**
- **R3** is [M26]'s cornstarch at high `U₀`.

**R2 and R3 look the same from outside and have opposite pore pressure.** That is the
discriminator (E-3 below), and it is the reason the claim is worth making: the paper the
lab is answering to has no R2 in it.

---

## 4. Why silica can show both regimes without reaching `ID` ≈ 1

The premise of the question — *both* solid-like jamming *and* an inertial-looking
thickened state, in **silica** — sat awkwardly with `ID ∝ a²` putting 800 nm silica 156×
below 10 µm cornstarch and 3–4 decades below the nominal boundary. Two independent
resolutions now exist, and they are not alternatives; either alone suffices.

1. **R2 is not R3.** The high-CV state does not need `ID` ≳ 1. It is on the other side
   of the *sign* boundary, not the drainage one, and reaching it costs no permeability
   at all.
2. **The `a²` in `ID` is a property of a length choice, not of the physics.** Deriving
   `ID` rather than quoting it shows a drainage length `ℓ` entering **squared**
   (derivation Eq. 8). [M26] set `ℓ ~ D₀` explicitly and flag it with "~". At `ℓ = D₀`
   the lab's silica has `Π` = 1.39e4 (deeply undrained). At `ℓ = 10 d₃₂` = 8 µm it has
   `Π` = 0.142, i.e. `ID` = 7.0 — **on the inertial side** — and every factor of particle
   size cancels exactly, at every size:

   | `a` [µm] | `ID` (ℓ = D₀) | `ID` (ℓ = 10 d₃₂) |
   |---|---|---|
   | 0.8 | 7.19e-05 | 7.03 |
   | 2.5 | 7.03e-04 | 7.03 |
   | 8.0 | 7.19e-03 | 7.03 |

   *(printed by `derivations/dilatancy_sign_mechanism_numbers.py`)*

**Nothing here chooses between the two length readings, and they differ by 10⁵.** That
is now the largest single uncertainty in the lab's use of `ID`, larger than the
boundary-calibration landmine already on the card (`ID` = 6.6e-3 at [M26]'s own
Regime-2/3 transition), and it is separately testable — see E-4.

---

## 5. The simulation: what is being measured and why it needed new code

`phase2_phiJ.py` measures `φ_J` by isotropic compression. **`φ_c` cannot be obtained
from it at any resolution**: isotropic compression asks how dense a pack can get, and
the critical state asks what density sustained *shear* drives it to. They are different
numbers with different friction dependences.

`shear.py` + `phase8_criticalstate.py` run **constant-volume simple shear**, which is
the undrained test of soil mechanics: hold the total volume, shear, and the volume
change the skeleton wanted appears as a change in grain pressure instead. `p_ss(φ)` at
constant volume is — up to the pore-fluid bookkeeping — the pore-pressure sign the
experiment sees. `velGrad = (0, γ̇, 0; 0,0,0; 0,0,0)` has trace exactly zero, so `φ`
during the shear is an input known to machine precision, not a controller output
(measured drift: `0.0e+00`).

Production grid: lognormal, `μ_p` = 0.3, CV ∈ {0, 0.2, 0.5}, 7 volume-fraction rungs at
0.015 spacing centred on the published isotropic `φ_J`, 3 seeds, N = 1500, `γ_max` = 2.0,
300 000 steps — 63 cells. Four `φ_c` estimators are reported side by side
(`analyse_phase8.py`): the `φ(I → 0)` extrapolation (primary; the standard granular
dilatancy law, **da Cruz, Emam, Prochnow, Roux & Chevoir, PRE 72, 021309 (2005),
doi:10.1103/PhysRevE.72.021309**, and **Forterre & Pouliquen, Annu. Rev. Fluid Mech. 40,
1 (2008), doi:10.1146/annurev.fluid.40.111406.102142** — ⚠ **their numerical coefficients
are 2D disks and are NOT quoted here**, only the form of the law), a free-exponent
pressure root, a pinned-exponent pressure root, and the coordination crossing at
`Z_iso` = 4. A contrast is quoted only if all four agree in sign; their spread is the
systematic and it is printed.

**Two defects were found while building it**, both by checks that could fail:

- **The shear-stress sign.** The first version returned `μ_macro` = −0.362. A negative
  macroscopic friction is impossible, so the flip was caught on the very first cell. The
  correct value, +0.362 at `μ_p` = 0.3, is also a plausible steady-state friction, which
  is a second reason to believe the protocol.
- **🔴 Stage-B compression quantized by its own chunk.** Compression ran in fixed
  500-step blocks and stopped at the first block past the target. One block is a linear
  strain of 2.5e-2, i.e. a jump of ~0.04 in `φ` — **nearly three times the 0.015 ladder
  spacing the whole measurement rests on.** Two rungs 0.015 apart landed on the *same
  state* and returned bit-identical `p_hat` (1.6786e-04), `Z` and rattler fraction. This
  is the "detector versus its own quantum" failure again, in a new place: the *stopping
  rule* was quantized. Fixed with a signed, three-scale targeting loop; `analyse_phase8`
  now carries a duplicate-state guard and fits on the **realised** `φ`, never the target.

Pre-registered predictions P1–P7 are in `derivations/dilatancy_sign_mechanism.md` §5.
Results and grading: §8 below.

---

## 6. The experiment

Five items, in dependency order. E-0 is blocking; E-2 is the one that decides it.

### E-0 (blocking, and already on the card) — which moment was matched?

Everything downstream is conditioned on the two existing bench samples being matched on
the **number** mean. If they were matched on a volume-weighted or DLS `z`-average, `d₃₂`
moved much more than assumed and the permeability channel has to be re-graded before the
packing one can be believed. **Raw PSDs, not a summary statistic.**

### E-1 (primary) — a CV ladder, run TWICE

Blend monodisperse Stöber batches so the PSD is *constructed*, not fitted. Two ladders
over the same CV values:

| ladder | held fixed | what it isolates |
|---|---|---|
| **A** | number mean `d₁₀` | both channels move (this is the bench observation as it stands) |
| **B** | **Sauter mean `d₃₂`** | `κ` is CONSTANT by construction ⇒ **any regime change is the packing channel, with no arithmetic required** |

Ladder B is the M1 null done properly: it does not *argue* that permeability is too weak,
it *removes* it. Blending makes it available — `d₃₂` is a computable functional of the
blend, so `d₁₀` is simply adjusted downward as CV rises to hold `d₃₂`.

**Prediction:** both ladders switch regime, at CV\* values that differ by less than the
`d₃₂` change between them.
**Falsified if:** ladder B shows no switch.

### E-2 (the decisive test) — the φ-reversal, with no free parameter

Take the high-CV sample that spreads. **Raise `φ` and the solid-like branch must come
back.** The mechanism predicts *how much*: the required increment is exactly the DEM's
measured `Δφ_c(CV)`, because the boundary is `φ = φ_c(CV)`.

```
   φ_needed(CV_high) − φ_needed(CV_low)  =  Δφ_c(CV_low → CV_high)     [DEM, no fit]
```

No other candidate mechanism predicts a number here. A viscosity-shift account predicts
a smooth weakening with no return; a permeability account predicts the increment is
independent of CV.

⚠ **This is why the evaporation problem is now load-bearing rather than merely annoying.**
With `dφ_c/dCV` ≈ 0.12, resolving CV\* to ±0.02 needs `φ` held to ±0.0025. The Status
card records the drop drying on the nozzle at **3.4× that tolerance before impact**. The
reversal test cannot be run until that is closed.

### E-3 (the discriminator) — the SIGN of the pore pressure

R2 and R3 both spread. They differ in the sign of the pore pressure. This is also the
**only** clean separation from the rival mechanism M2a (§7a), so it is load-bearing
rather than confirmatory.

**A flush-mounted transducer** in the impact substrate: suction below atmospheric on R1,
overpressure on R2, neither on R3. Bandwidth has to beat the contact time `D₀/U₀` ≈ 0.8 ms.

⚠ **A cheaper version was proposed and is WITHDRAWN.** The first draft suggested looking
for cavitation bubbles in footage already recorded, on the grounds that the dilatant
suction exceeds one atmosphere. It does not: with the correct ceiling (derivation §7,
rebuilt after the referee pass) the suction is bounded at **26.8 Pa**, four decades below
atmospheric, so **nothing cavitates and there is nothing to see.** The free measurement
was an artefact of a withdrawn cap.

**Falsified if** the high-CV spreading state shows pore suction.

**And a magnitude caveat that came out of the same correction.** Equation (17) of the
derivation gives the dilatant strengthening as a fraction `(a/6λ_D)(Δφ/φ)` of the
particle pressure: **0.20 at nominal zero salt, 0.61 at 1 mM, 6.06 at 100 mM.** So the
pore-pressure route is *marginal at the lab's usual low-salt condition* and becomes
available with screening. That is a prediction on an axis Paper 1 already sweeps, and it
is the mechanism's most exposed flank.

### E-4 — the size ladder, which settles the drainage length

Three particle sizes at fixed CV, `φ`, `U₀` (0.8 / 2.5 / 8 µm). **Three hypotheses, three
different exponents, one decade of `a`:**

| reading | prediction | over the decade |
|---|---|---|
| macroscopic Darcy drainage, `ℓ = D₀` (this is what `ID ∝ a²` assumes) | `∝ a²` | **100×** |
| grain-scale drainage, `ℓ ∝ a` | `∝ a⁰` | **1×** |
| **undrained, network-modulus limited — what this analysis predicts** | **`∝ a¹`** | **10×** |

The third row is new (derivation §6a, from `diag_kappa_sweep.py`): in the undrained limit
the pressure follows `K_p`, not `κ`, and `K_p = P a/(6λ_D)` for a charge-stabilised
contact is linear in the particle size. ⚠ That exponent is specific to colloidal silica —
a Hertzian contact carries no size at all — so it must not be carried to cornstarch.

The Paper 1 card already wants exactly this ladder for the relaxation-clock question
("one material, three particle sizes, same chemistry, same `φ/φ_M`" — *"a null result is
publishable on its own"*). **The same samples answer both questions**, which is the
cheapest thing on this page.

### E-5 (cheap, from data already recorded) — the arrest boundary in φ

The mechanism says arrest requires `φ > φ_c` at the contact. The existing per-drop table
has an arrest classification against `φ`. Two things to look for: whether the boundary is
a **threshold** rather than a gradual weakening, and whether the drops nearest it show
the largest run-to-run scatter (a threshold crossed by evaporation drift would).

---

## 6a. The second simulation leg: the transient, and the strain clock

`poro_impact.py` solves the 1D consolidation problem (Terzaghi + a Reynolds source +
`φ` evolution) that the closed form only takes the quasi-steady limit of. It exists
because the closed form cannot say two things the bench observation needs.

1. **The O(1) prefactor.** [M26]'s Eq. 6 is written with "~". Solving the boundary-value
   problem gives it: for an impermeable substrate and a drained free surface,
   `p̂_peak = Π/2` exactly. Pinned by `test_poro.py`.
2. **Whether the balance is reached at all.** It is not, and by six decades. At the lab's
   `Λ` ≈ 2.5 × 10⁻⁷ the diffusive steady state of Eq. (6) would need a strain of order
   `1/Λ` ≈ 4 × 10⁶. **Darcy plays no part during the impact**; the response is undrained
   from the first instant.
3. **What the undrained pack then does — and this reversed a conclusion.** As it dilates
   plastically the pore pressure falls, so the effective stress rises, so the skeleton
   compresses elastically by exactly what it dilated. Net volume change zero, `φ − φ_c`
   **locked**, and the suction **grows linearly with strain without bound**. This is
   **dilatancy hardening**, the standard reason a dense saturated sand is strong under
   fast loading, and it predicts a body that **stiffens as it spreads**.

🔴 **A referee found the first version of this module missing the elastic half of its own
`φ`-evolution law**, which produced a spurious self-limiting decay. On that decay this
page previously claimed a **"strain clock"** and offered it as one of three separations
from the rival M2a. **That claim is withdrawn.** The clock is real only in the *drained*
limit, which is not where this system is. What bounds the undrained growth is a closure
the model does not carry — the critical state is a *line*, so `φ_c` falls as `σ'` rises —
and the network-modulus ceiling stands in for it.

**Why no test caught it:** the one test that exercised the evolution law re-integrated
the *same incomplete equation* in a separate scalar solver. Two anchors sharing one
identity cannot check that identity. Replaced by a closed form — in the undrained limit
`φ − φ_c` must be constant to machine precision and `p(t) = −SαΔφ₀ t` exactly — plus a
negative case that the two laws disagree, and a test that the clock appears **only** at
large `Λ`. **15 tests.**

The other anchors are outside the code under test: the published Terzaghi consolidation
series (max error 2.4e-7, observed grid order > 1.7) and the hand-derived `Π/2` with a
negative case that rejects a prefactor of 1.

🔴 **The Terzaghi anchor also caught a real defect on its first run**: the analytic series
had been written with the coordinate measured from the *sealed* face while the solver
drains at the other end, i.e. the mirror image, and it disagreed by 0.61 in a field of
order 1. An anchor that reproduced the code's own convention would not have found it.

---

## 7. What would kill the mechanism outright

| observation | verdict |
|---|---|
| `φ_c` flat or falling with CV in the DEM (P1) | dead — the sign channel cannot move |
| ladder B (matched `d₃₂`) shows no regime switch | dead — it *was* permeability |
| the high-CV spreading state shows pore **suction** (E-3) | dead — the sign is the mechanism |
| raising `φ` at high CV does not restore the solid branch (E-2) | dead |
| raising `φ` restores it, but by an increment unrelated to `Δφ_c` | mechanism survives, the quantitative version does not |

## 7a. 🔴 The rival that the DEM cannot separate

A referee named a mechanism this page had not, and it is the most serious open issue on
it.

> **M2a, plain suspension rheology.** Raising CV raises `φ_m` (≈ `φ_c`). At fixed `φ`
> that moves the suspension further from jamming, `η ∝ (φ_m − φ)⁻²` falls, and the drop
> spreads. **No pore fluid, no Terzaghi, no dilatancy sign.**

M2a uses the *same measured quantity* (`dφ_c/dCV`), predicts the *same direction*, and
puts the crossing in *roughly the same place*. **Every DEM prediction P1–P6 is consumed
identically by both.** So the phase-8 campaign, however it comes out, establishes that
the CV observation is a packing-geometry effect — and does not establish *which* packing
mechanism.

Three separations, in order of cost:

| separation | M2′ (this mechanism) | M2a | where |
|---|---|---|---|
| **pore-pressure sign** | overpressure when spreading, suction when arrested | **makes no prediction** (a scalar viscosity has no signed-pressure coupling) | E-3 |
| **strain history** ⚠ *restated* | undrained **locking**: the suction grows with strain, so the body **stiffens as it spreads** | none — a viscosity is a state function with no strain history | `poro_impact.py` + high-speed footage |
| **drainage number** | the effect weakens as `Λ` rises; locking is an undrained phenomenon | `Λ` does not appear | E-4's size ladder moves `Λ` by 10² |
| **salt** | strengthens with screening, `∝ a/6λ_D` | `φ_m` is geometric; salt enters only via `σ*` | Paper 1's existing salt series |

⚠ The **strain-history** row previously read "a strain clock; the effect decays as
`φ → φ_c`". That was wrong — see §6a — and the corrected version is a *better*
separation, not a weaker one. ⚠ And the "cannot separate" claim is itself contingent on
`φ_m ≈ φ_c`; if those two shifts were measured to diverge, P1–P6 would partially
discriminate without E-3.

Note that §2's revival of M2 as "M2′" was, strictly, a revival of *both*. The honest
statement until E-3 is run: **the width acts through the packing, and whether it acts
through the pore fluid or through the viscosity is open.**

---

## 8. Results and grading of P1–P7

*(To be filled from `out/phase8_verdict.txt`. Nothing may be written into this section
before the sweep returns, and the bars in §5 of the derivation do not move.)*

---

## 9. Referee

*(To be filled from the `physics-referee` pass.)*

---

## 10. Position in the literature

- **[M26]** Mobaseri, Gordillo, Burton, Yoon, Lee, Kumar, Driscoll & Cheng, *PRL* **136**,
  148201 (2026), doi:10.1103/fyx7-jb1d, arXiv:2601.12642 — the source of `ID`. This work
  *derives* their Eq. 6 rather than adopting it, and in doing so exposes the drainage
  length and the `Δφ` pole that their fixed-`φ_c` presentation does not show. **Driscoll
  is a co-author and is local**, which makes both the boundary calibration and the
  drainage length settleable by asking.
- **Desmond & Weeks**, *PRE* **90**, 022204 (2014), **doi:10.1103/PhysRevE.90.022204** —
  ✅ **resolved 2026-09-08**; the card carried it as `[UNVERIFIED]`. Frictionless RCP vs
  polydispersity; this project extends it to frictional shear.
- **Critical-state soil mechanics** (Schofield & Wroth 1968) — the framework for `φ_c`
  and for the liquefaction/dilatancy dichotomy. `[UNVERIFIED — not resolved to a DOI in
  this session; it is a book and the citation should be checked before it is printed.]`
- **What is new**, as far as the project's own literature pass can tell: not
  "polydispersity weakens thickening", which is published and must not be claimed; but
  **(i)** that the width acts through the *sign* of the dilatancy rather than through its
  magnitude, **(ii)** the existence of a critical width CV\* with a parameter-free
  prediction for how `φ` moves it, and **(iii)** that `ID`'s particle-size scaling is a
  property of an unmeasured drainage length. `[UNVERIFIED — a literature pass against
  (i)–(iii) has not been run and must be before any of this is written up.]`
