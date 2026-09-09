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

No other candidate mechanism predicts a number here. A permeability account predicts the
increment is independent of CV. ⚠ A *viscosity* account (M2a, §7a) predicts a smooth
weakening — but note it also predicts a return on raising `φ`, so E-2 tests the
**quantitative** prediction (the size of the increment) and not the qualitative one.

**Now that `Δφ_c` is measured, the requirement can be stated properly — and it is looser
than I first wrote.**

```
   predicted increment = Δφ_c(CV 0 → 0.5) = 0.0623 ± 0.003   (measured, 20.6 sd)
```

The experiment must resolve the **difference of two thresholds** against that 0.062. If
each threshold is located by stepping `φ` in increments of `δ`, each is known to about
`δ/2` and the difference to `δ/√2`. Testing 0.062 at 3σ therefore needs
`δ/√2 ≲ 0.021`, i.e.

> **`φ` steps of ≲ 0.03 and drift held below ≈ 0.015** — comfortably coarser than a
> concentration series already has to be.

🔴 **CORRECTION to what this section said earlier today.** It read "resolving CV\* to
±0.02 needs `φ` held to ±0.0025 … the reversal test cannot be run until [evaporation] is
closed". That tolerance belongs to a *different* measurement — **locating CV\*** on the
CV axis, which does need `φ` to ±0.0025 because `dφ_c/dCV` ≈ 0.12 converts a CV
resolution into a much tighter `φ` one. **E-2 does not locate CV\***; it measures a
*shift* of 0.062, and a shift that large needs an order of magnitude less `φ` control.
E-2 is **not blocked by evaporation**; locating `CV*` is.

⚠ Evaporation still matters here — a drift comparable to the step size blurs each
threshold — but the bar is ≈ 0.015, not 0.0025.

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

### E-6 🆕 (added 2026-09-08) — the nozzle test, and it needs no new material either

`drop_size_test.py`. The drainage length is the derivation's largest stated uncertainty
and E-4 proposed a particle-size ladder to settle it. **Drop diameter tests the same
length with a different nozzle, and the readings differ in SIGN.**

| reading | `γ̇` | pressure length | boundary velocity |
|---|---|---|---|
| (M) macroscopic, = [M26] Eq. 6 | `U₀/D₀` | `D₀` | **`U₀* ∝ D₀`** |
| (B) band drainage | `U₀/D₀` | `ξ` | **`U₀* ∝ 1/D₀`** |
| (L) fully localised | `U₀/ξ` | `ξ` | **independent of `D₀`** |

Pinned to the same observation so all three agree at `D₀` = 4 mm, halving the drop to
2 mm gives **3.6 / 14.4 / 7.2 m/s** — a 4× separation between the extremes, decided by
the sign, so it depends on neither the band thickness nor the `ID` prefactor nor where
the boundary sits.

**Run it on cornstarch, not silica** (silica's boundary is 278 m/s, unreachable), and
score the **regime**, not the spread factor — `We` and `Re` move with `D₀` by themselves,
so the control is a Newtonian liquid of matched viscosity and density on the same ladder.

🔑 **This is also the second hardware-free separator from the rival M2a**, which has no
boundary in `(U₀, D₀)` at all. Until now the only clean separator was E-3's
pore-pressure sign, which needs a transducer.

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

## 6b. 🔴 The largest quantitative gap: which side of `φ_c` is the lab actually on?

`phi_axis_offset.py` → `out/phi_axis_offset.txt`.

The mechanism turns on the **sign** of `φ − φ_c`. So the first question in applying it to
the lab's own experiment is which side the lab sits on — and on the numbers now measured,
the answer is the wrong one:

```
   lab working point                          φ   = 0.4400
   MEASURED critical state (μ_p = 0.3, CV = 0) φ_c = 0.5695 ± 0.0030
   φ − φ_c                                        = −0.1295
```

Taken at face value the lab is **0.13 below its critical state**, i.e. permanently
contractive and never solid-like — which contradicts the observation the project exists
to explain.

**Friction cannot rescue it.** `φ_J` saturates with `μ_p` on the project's own A4 ladder:

| `μ_p` | 0 | 0.3 | 1 | 3 |
|---|---|---|---|---|
| `φ_J` | 0.6332 | 0.5745 | 0.5590 | **0.5574** |

The step from `μ_p` = 1 to 3 is −0.0016, inside the seed scatter. So there is a **floor**
at ≈ 0.557, and

> **No amount of interparticle friction can bring a hard-sphere critical state down to
> φ = 0.44. The gap friction cannot close is 0.117.**

**One of three things must therefore be true**, and they are distinguishable:

1. **The local `φ` during impact far exceeds the bulk 0.44.** Reaching the measured
   `φ_c` needs **1.29×** densification, a **23 % volumetric compaction** at the contact.
   That is the number the Status card's next-action #2 ("measure the local `φ` during
   impact") has to beat.
2. **The material is not a frictional hard-sphere pack.** 800 nm silica at low salt is
   charge-stabilised — the particles never touch — and a repulsive network is rigid at
   far lower solid fractions than a hard-sphere pack. **This DEM has no repulsion, so its
   `φ_c` is the wrong quantity for that material and the floor does not apply.**
3. **The experiment's `φ` axis is systematically low.** The Paper 1 card already records
   three irreconcilable numbers here (fronts at 0.42–0.45, first solid patches at 0.47 in
   Shah *et al.*, a fitted `φ_M` of 0.5176). What is new is that a controlled simulation
   now contributes a **floor** rather than a fourth free number.

**What this does and does not touch.**

- **Untouched: every CV claim.** They are all *differences at fixed protocol* —
  `dφ_c/dCV` = +0.125, the 43× collapse at matched `φ` — and a common offset in the `φ`
  axis cancels out of all of them.
- **Touched: any statement about which side the lab is on**, including the critical width
  `CV*` and the regime map's assignment of the lab to R1. Those are **conditional** on
  resolution (1), (2) or (3) and should be read that way.

🔑 **And it makes one cheap thing urgent.** Option (2) is decidable *in this same DEM* by
adding a repulsive contact law and re-running one ladder — a bigger change than a sweep,
far smaller than a new experiment, and it would say whether a charge-stabilised pack has
a critical state anywhere near 0.44.

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
| **a boundary in `(U₀, D₀)`** 🆕 | there is one, and its slope has a sign per drainage reading | **there is none at all** — `φ_m` and `φ` are material properties | **E-6, a nozzle change** |
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

45 cells, protocol `cv-shear-v2`, lognormal, `μ_p` = 0.3, N = 1500, CV ∈ {0, 0.2, 0.5},
5 volume-fraction rungs at 0.020 spacing, 3 seeds, `γ_max` = 1.5, 200 000 steps.
**0 failed, 0 duplicate states, all cells quasi-static** (`I` = 2.3e-4 … 1.7e-2, median
~1e-3). Source: `out/phase8_verdict.txt`, `out/phase8_summary.csv`. Bars unchanged from
`derivations/dilatancy_sign_mechanism.md` §5.

| # | prediction | bar | result | verdict |
|---|---|---|---|---|
| **P1** | `φ_c` rises with CV | Δ > 0 and > 3× seed scatter | **+0.0578 (E1, 46 sd)**, **+0.0623 (E2, 21 sd)**; all resolved estimators agree in sign | ✅ **PASS**, by >6× the bar |
| **P2** | critical-state slope ≈ isotropic `φ_J` slope | within 2× of +0.1246/CV | **E1 0.1157, E2 0.1247** — within **8 %**, and E2 within 0.1 % | ✅ **PASS** |
| **P3** | at fixed `φ`, raising CV drives the pack toward/past critical | `p_ss` falls and crosses the unjammed bar | at matched `φ`, `p_ss` falls **5.4×–43×** (15–63 sd); `Z_nr` 6.1 → 4.4, rattlers **0.004 → 0.448** | ✅ **PASS** on the fall; ⚠ *partial* on "crosses" — see below |
| **P4** | frictionless control `φ_c(μ_p=0)` ≈ `φ_J`(μ_p=0) = 0.6332 | within 0.01 | **0.6347 ± 0.0010** (E4, the isostatic crossing at Z = 6) — off by **0.0015** | ✅ **PASS**, 7× inside the bar |
| **P5** | uniform vs lognormal at matched CV give different slopes | ratio > 1.5 | **NOT RUN** | ⬜ **not tested** |
| P6, P7 | bench | — | for the author | — |

### The headline measurement, and it needs no fit at all

The estimators disagree on the *absolute* `φ_c` by 0.074 (E1 ≈ 0.643, E2 ≈ 0.570 at
CV = 0) — **a systematic larger than the effect being measured**, which is why this
project's convention is to quote contrasts and never absolutes. The contrast is robust:
every resolved estimator gives the same sign and a slope within 8 %.

But the cleanest result involves no estimator at all. Where the CV = 0 and CV = 0.5
ladders overlap, two rungs nearly coincide in `φ`:

| `φ` | CV | `p_ss` | `Z_nr` | rattlers |
|---|---|---|---|---|
| 0.63453 | 0.00 | 1.082e-2 | 6.12 | 0.004 |
| 0.63687 | 0.50 | **2.506e-4** | **4.43** | **0.448** |
| 0.65453 | 0.00 | 1.636e-2 | 6.49 | 0.002 |
| 0.65686 | 0.50 | **3.054e-3** | **5.15** | **0.337** |

> **At matched packing fraction, raising CV from 0 to 0.5 collapses the undrained grain
> pressure by 5.4× to 43×.** In the mechanism's terms that is the pore-suction capacity
> collapsing at fixed `φ` — the thing the whole argument turns on — measured directly.

Two honesty notes. **Both pairs are reported**, not the closest: they tie to within
rounding on `φ` mismatch (0.0023) and differ 8× in the answer, so a `min()` tie-break
would have been cherry-picking — the ratio grows as you approach `φ_c`, which is expected
and is why the **range** is quoted. And in both pairs the higher-CV cell sits at the
*denser* `φ`, which works **against** the effect, so the ratios are conservative.

⚠ **P3's "crosses the unjammed bar" is only partial.** `Z_nr` = 4.43 at the closer pair
is *at* the frictional isostatic value, i.e. marginal, not clearly unjammed. The fall is
decisive; the crossing is not demonstrated by these cells.

### ✅ P4 ran, and the protocol now has an external anchor

15 cells at `μ_p` = 0, ladder centred on the published `φ_J`(μ_p=0) with rungs on both
sides. `grade_p4_anchor.py` → `out/p4_anchor_verdict.txt`.

| estimator | mean | \|diff\| from 0.6332 | verdict |
|---|---|---|---|
| **E4** — where `Z_nr` crosses the **isostatic** count (6, frictionless) | **0.6347 ± 0.0010** | **0.0015** | ✅ **PASS** |
| E2 — free-exponent pressure root | 0.5839 (1 of 3 seeds; 2 refused) | 0.0492 | FAIL |
| E3 — same, `b` pinned to 1 | 0.6128 | 0.0203 | FAIL |

This matters more than the other predictions because `φ_J`(μ_p=0) was measured months
earlier by a completely different route (an isotropic pressure ladder), and the physics
that makes the two coincide — frictionless packings have no shear-jamming branch below
the isotropic point — is not something this protocol can arrange for itself.

⚠ **The estimators disagree and the pre-registration did not name one.** That is a real
weakness in how P4 was written, so the choice is argued on **definitions**, not on which
one passes: for a *frictionless* pack, jamming **is** the isostatic point, so E4 measures
exactly the anchored quantity with no fitted form at all. E2's exponent ran to its
**bound** on 2 of 3 seeds (`b` = 8.00, 7.33, 8.00) — a fit on a bound is not a
measurement, and `b ≥ 7` says the transition is far sharper than any low-order power law,
which is *correct* for frictionless jamming and is the fit reporting its own inadequacy.
E3 pins `b` = 1 while the data says `b ≥ 7`.

⚠ **On the frictional production arm the ordering reverses** — there `b_E2` = 1.5–3.6,
nowhere near the bound, and E4 mostly finds no crossing, which is why E2 is primary
there. **The right estimator depends on the friction**, and one rule would have been
silently wrong in one of the two arms.

🔴 **A bug found on the way in.** The isostatic count was hardcoded at 4. It is `2d` = 6
for frictionless spheres and `d+1` = 4 for frictional ones — right for the `μ_p` = 0.3
arm and silently wrong for this one, where E4 would have hunted for a crossing at 4 that
happens far below jamming. Now chosen per cell from `μ_p`, with a test pinning both.

**A free internal check:** `μ_s` = **0.105 ± 0.003** at `μ_p` = 0 against
**0.359 ± 0.010** at `μ_p` = 0.3, a factor **3.42**. Removing particle friction cuts the
macroscopic friction by that factor and does *not* take it to zero — a frictionless pack
still resists shear geometrically. Both the direction and the non-zero floor are what the
`μ(I)` reading must show if it is measuring friction at all. ⚠ Both numbers come from
this protocol, so it is an internal consistency check, not an anchor.

### What is still not tested

**P5 is running at the time of writing** — the uniform, zero-skew PSD at matched CV,
which is the discriminator that separates a CV effect from a `d₃₂` effect. The isotropic
analogue is already 2.7×.

### One estimator was disqualified mid-analysis

E3 (pressure root with the exponent pinned to 1) first returned Δ`φ_c` = +0.0612 with a
seed sd of **0.0000** and an apparent significance of **161 sd** — and it was returning
its own **fit bound** (`lo − 1e-6`, i.e. the ladder's bottom rung) on 6 of 9 groups. A fit
on a bound is not a measurement; `analyse_phase8.py` now detects and NaNs it. Had it been
believed it would have been the most significant number on the page.

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
