# The dilatancy-sign mechanism: how mean size and polydispersity set the impact regime

**Written 2026-09-08.** Checks: `derivations/dilatancy_sign_mechanism_check.py`
(`python derivations/dilatancy_sign_mechanism_check.py` → `ALL CHECKS PASSED`,
17 checks). Every number in this prose is printed by
`derivations/dilatancy_sign_mechanism_numbers.py` into
`dilatancy_sign_mechanism_numbers.txt`; none is typed by hand. Referee: see §9.

---

## Result box

> A saturated dense-suspension drop has **two independent grain-scale control axes**,
> and the size distribution moves one of each.
>
> **Axis 1 — drainage (set by the MEAN size).** How much pore pressure the skeleton's
> volume change generates, relative to the impact stress:
>
> ```
> Π  =  α μ_w |φ − φ_c| ℓ² / (ρ κ U₀ D₀) ,      κ = (1−φ)³ d₃₂² / (180 φ²)
> ```
>
> With the drainage length ℓ = D₀ this is **exactly 1/ID of [M26] Eq. 6** — its
> *functional form* derived here from Reynolds + Darcy rather than quoted (check C1.3).
> ⚠ The **coefficient is imported, not derived**: `α` occupies the same slot as [M26]'s
> `A`, which is their fitted granular-friction closure (their Eq. 5) and carries no
> Reynolds–Darcy content. Π ≫ 1 is undrained, Π ≪ 1 is drained.
>
> **Axis 2 — sign (set by the WIDTH).** Which way the pore pressure goes:
>
> ```
> sign(Δp_f)  =  −sign(φ − φ_c(CV))
> φ > φ_c :  dilates → suction  → effective stress RISES → solid-like, jammed
> φ < φ_c :  contracts → overpressure → effective stress FALLS → liquefied, spreads
> ```
>
> **The two CV channels are not comparable in kind.** Permeability enters ID as
> `(1+CV²)⁴`, whose logarithmic derivative **vanishes at CV → 0**; the critical state
> enters as `1/(φ − φ_c(CV))`, whose logarithmic derivative **does not** (checks C5.3,
> C5.4). So for any non-zero `dφ_c/dCV`, the leading-order CV dependence of the regime
> is entirely the packing one — **and it has a pole at CV\*, where φ_c(CV\*) = φ**,
> which no smooth permeability factor can imitate.

---

## 1. Frame

**Quantity sought.** Which impact regime a saturated dense-suspension drop lands in —
the non-inertial solid-like/jammed branch or the spreading, inertia-dominated one — as
a function of the particle-size distribution at fixed solid fraction φ.

**The observation being explained** (author's bench note, 2026-08-31, verbal, recorded
in the Status card): *at similar particle size, raising CV moves the impacting drop from
purely dilatant to spreading and inertial-dominated for some time.* Both regimes are
seen in **silica**, i.e. within one chemistry and one particle shape.

**Governing physics, four ingredients, all standard.**

| | ingredient | statement | source |
|---|---|---|---|
| G1 | Reynolds dilatancy | a dense skeleton sheared at rate `γ̇` changes volume at `ε̇_v = γ̇ tan ψ`, with `tan ψ = α (φ − φ_c)`, `α ~ 1` | [M26] Eq. 3 and refs therein; Reynolds 1885 |
| G2 | critical state | under sustained shear a granular material approaches `φ_c` from either side, independent of its initial density | critical-state soil mechanics (Schofield & Wroth 1968) `[UNVERIFIED — not resolved to a DOI in this session]` |
| G3 | Darcy | `u_D = (κ/μ_w) ∇p_f`, `κ` from Carman–Kozeny with the **Sauter mean** `d₃₂` | [M26] Eq. 3; `psd.carman_kozeny_kappa` |
| G4 | Terzaghi | `σ_total = σ' + p_f` (compression positive); grain contacts carry `σ'` only | standard |

**Assumptions, and the small parameter for each.**

- **A1** Saturated, incompressible pore fluid. Fails if air is entrained or the fluid
  cavitates — see §7, where it *does* fail for silica.
- **A2** One drainage length `ℓ` describes the pore-pressure field. This is the
  substitution `r_c ~ h ~ D₀` that [M26] make explicitly and flag with "~". §6 makes it
  the sharpest open question in the whole argument.
- **A3** Carman–Kozeny applies, with `a = d₃₂`. CK's length is the surface-area
  equivalent diameter, which for a polydisperse pack *is* the Sauter mean — this is an
  identity, not an approximation, but CK itself is marginal here: `λ_D/r_throat` =
  0.10–0.48 for these particles (`impact_stress/PORE_PHYSICS.md` §5.2).
- **A4** `φ_c` is a property of the packing geometry (hence of the PSD and the
  interparticle friction), not of the fluid. This is what makes a hydrodynamics-free
  DEM the right instrument for `dφ_c/dCV`.
- **A5** Small strain: `κ`, `φ` and `φ_c` are evaluated at the instantaneous state and
  not evolved.

**Conventions.** Compression-positive for stresses. `d` is a DIAMETER throughout.
`φ` includes rattlers. `α` is [M26]'s dilatancy coefficient; their `A` (Eq. 5) is the
same slot in the algebra and the check treats them as one symbol (C1.3).

**Dimensionless groups at the actual experiment** (800 nm silica, φ = 0.44, U₀ = 3 m/s,
D₀ = 2.5 mm, water):

| group | value | meaning |
|---|---|---|
| `Re = ρU₀D₀/μ_w` | 1.1e4 | inertia dominates the *fluid* |
| `We = ρU₀²D₀/γ` | 4.5e2 | surface tension is not the impact-scale resistance |
| `St = ρ_p a²U₀/(9μ_w D₀)` | 1.7e-4 | grains follow the fluid at the drop scale |
| `Π (= 1/ID)` at ℓ = D₀ | 1.4e4 | **deeply undrained** on the macroscopic length |
| `Π` at ℓ = 10 d₃₂ | 1.4e-1 | **drained** on the grain scale, but only by 7× — §6 |

The last two lines are the whole of §6: the same suspension is on opposite sides of the
drained/undrained boundary depending on a length nobody has measured.

---

## 2. Derivation, every step

**Step 1 — what the skeleton demands.** By G1, a region of size `ℓ` sheared at `γ̇`
changes its volume at a fractional rate

```
(1)   ε̇_v = γ̇ tan ψ = α γ̇ (φ − φ_c)          [1/s]
```

so the volume the pore space must gain or lose, per unit area of that region per unit
time, is

```
(2)   q = ε̇_v ℓ = α γ̇ (φ − φ_c) ℓ            [m/s]
```

**Step 2 — what the fluid can supply.** By G3, with the pressure difference `Δp_f`
developing over the same length `ℓ`,

```
(3)   u_D = (κ/μ_w) (Δp_f/ℓ)                   [m/s]
```

**Step 3 — the balance.** Incompressibility (A1) forces `u_D = q`, hence

```
(4)   Δp_f = α γ̇ (φ − φ_c) μ_w ℓ² / κ         [Pa]
```

*(SymPy: check C1.1.)* The sign convention: `φ > φ_c` means the skeleton dilates, the
fluid must be drawn **in**, and the pore pressure must fall — `Δp_f` in (4) is the
magnitude of a **suction** in that case, and of an **overpressure** when `φ < φ_c`.

**Step 4 — effective stress.** By G4, with the impact supplying a total normal stress
`σ_n`,

```
(5)   σ'_n = σ_n − p_f = σ_n − (p_0 + Δp_f^signed) ,   Δp_f^signed = −sign(φ − φ_c)|Δp_f|
```

so a **dilating** skeleton carries *more* than the applied stress and a **contracting**
one carries *less*. Frictional strength is `τ_max = μ σ'_n`, so the first is a
stiffening/solidifying branch and the second a liquefying one. This is the liquefaction/
dilatancy dichotomy of soil mechanics, transplanted to impact.

**Step 5 — nondimensionalise.** With `γ̇ ~ U₀/D₀` and the impact stress scale `ρU₀²`,

```
(6)   Π ≡ |Δp_f| / (ρU₀²) = α μ_w |φ − φ_c| ℓ² / (ρ κ U₀ D₀)
```

*(SymPy: check C1.2.)*

**Step 6 — the identification with the published number.** Putting `ℓ = D₀` and writing
`α → A`:

```
(7)   Π(ℓ = D₀) = A μ_w (φ − φ_c) D₀ / (ρ κ U₀) = 1 / ID          [M26] Eq. 6
```

*(SymPy: check C1.3 — exact, symbolically.)* So [M26]'s inertial-dilatant number is
recovered here from Reynolds + Darcy without being assumed, and the derivation exposes
what their "~" hides: **ID contains a choice of drainage length, and it enters squared.**

```
(8)   Π(ℓ) / Π(D₀) = (ℓ/D₀)²        ⇔        ID(ℓ) = ID(D₀) · (D₀/ℓ)²
```

*(check C1.4.)*

**Step 7 — where the size distribution enters.** Two places, and only two:

```
(9)    κ  = (1−φ)³ d₃₂² / (180 φ²) ,      d₃₂/d₁₀ = (1 + CV²)²   (lognormal, fixed d₁₀)
(10)   φ_c = φ_c(CV)   — measured, not assumed
```

Taking logs of ID and differentiating:

```
(11)   d ln ID / dCV  =  8CV/(1+CV²)   +   φ_c'(CV)/(φ − φ_c(CV))   −  d ln A/dCV
                         └ permeability ┘   └──── critical state ────┘
```

*(checks C5.1, C5.2.)* The third term is small: `A` depends on CV only through
`tan β = tan β₀ + α Δφ` with `β₀ = 44°`, so `d ln A/dCV = O(Δφ')` against the second
term's `O(Δφ'/Δφ)`, i.e. smaller by the factor `Δφ ≈ 0.02`.

**Step 8 — the two channels are different in kind.**

```
(12)   lim(CV→0) [permeability channel] = 0
(13)   lim(CV→0) [critical-state channel] = φ_c'(0)/(φ − φ_c(0)) ≠ 0
```

*(checks C5.3, C5.4.)* The permeability channel is `O(CV)`; the critical-state channel
is `O(1)`. **At small CV the packing channel wins for any non-zero `φ_c'`, whatever its
size** — this is a statement about orders, not about numbers, and it is the part of the
argument that does not depend on a measurement.

**Step 9 — the pole, and the critical width.** The critical-state channel diverges where
`φ_c(CV) = φ`. Linearising `φ_c(CV) ≈ φ_c(0) + s·CV`,

```
(14)   CV*  =  (φ − φ_c(0)) / s
```

Below `CV*` the drop is dilatant (solid-like branch). **Above `CV*` it is contractive,
and [M26] Eq. 6 does not describe it at all** — `ID` as coded returns `inf` for
`Δφ ≤ 0`, which is a placeholder for "no dilatant branch", not a physical regime label.
The contractive branch is a *third* regime with its own physics (positive pore pressure,
loss of effective stress), and it is the natural candidate for what the author sees at
high CV.

---

## 3. The regime map

Two axes, three regimes, at fixed φ:

```
                      |  φ < φ_c(CV)  (contractive)   |  φ > φ_c(CV)  (dilatant)
  --------------------+-------------------------------+---------------------------
   Π ≫ 1  undrained   |  R2  LIQUEFIED                |  R1  SOLID-LIKE / JAMMED
   (small d₃₂)        |  pore OVERpressure, σ' falls  |  pore SUCTION, σ' rises
                      |  spreads, looks "inertial"    |  arrests, can rebound
  --------------------+-------------------------------+---------------------------
   Π ≪ 1  drained     |  R3  RE-ENTRANT LIQUID (both columns): the pore fluid keeps
   (large d₃₂)        |  up, no pore-pressure feedback, [M26]'s Regime 3
```

- **R1** is the lab's low-CV silica: undrained *and* dilatant.
- **R2** is the prediction. It is reached from R1 **by raising CV alone**, at fixed
  chemistry, fixed φ and fixed mean size — which is exactly the author's experiment.
- **R3** is [M26]'s cornstarch at high `U₀`.

**R2 and R3 look alike macroscopically and have opposite pore pressure.** That is the
discriminating measurement (§5).

---

## 4. What this says that the project's earlier answer did not

The Status card records `M1` (permeability) as *"quantitatively excluded across most of
the range"*, on the grounds that CV 0 → 0.5 multiplies `ID` by 2.44× against a boundary
three to four decades away. That calculation held `Δφ` **fixed at 0.02** and swept only
`κ`. Equation (11) says `Δφ` is the *other* CV-dependent factor in the same expression,
and (14) says it has a zero inside the experimental range.

Using the isotropic `φ_J` slope as a placeholder for `s` (+0.0623 over CV 0 → 0.5, from
`out/phase2_summary.csv`, lognormal, `μ_p` = 0.3, N = 1500):

| CV | `Δφ = φ − φ_c` | critical-state channel ÷ permeability channel |
|---|---|---|
| 0.05 | +0.0138 | **22.7×** |
| 0.10 | +0.0075 | **20.9×** |
| 0.15 | +0.0013 | **81.0×** |
| 0.16 | 0 | **pole** |
| ≥ 0.20 | negative | past the pole — contractive branch, Eq. 6 does not apply |

✅ **MEASURED 2026-09-08, and the substitution holds.** `s` above was the ISOTROPIC
jamming slope standing in for the critical-state one — a different quantity, and an
assumption rather than a measurement. `phase8_criticalstate.py` has now measured the
critical-state slope directly, by constant-volume shear, over 45 cells:

| estimator | Δ`φ_c` (CV 0 → 0.5) | slope `s` per unit CV | vs the isotropic 0.1246 |
|---|---|---|---|
| E1, `φ(I → 0)` | +0.0578 (46 sd) | **0.1157** | 0.93× |
| E2, free-exponent pressure root | +0.0623 (21 sd) | **0.1247** | 1.0006× |

So the table above now stands on a measured slope and **P2 passes**. ⚠ A third estimator
appeared to agree to 0.1 % with a seed sd of *zero* and was disqualified — it was sitting
on its own fit bound (`MECHANISM_2026-09-08.md` §8). ⚠ `CV*` ≈ 0.16 still carries its
*other* placeholder, `φ_c(0)` = 0.42, which is unmeasured; the slope is no longer one.

What was *never* provisional is Step 8: the ordering of the two channels at small CV
holds for any `s ≠ 0`.

---

## 5. Predictions, pre-registered

Registered before `phase8_criticalstate.py` returned any `φ_c`. Bars are stated so they
can fail.

| # | prediction | bar | status |
|---|---|---|---|
| **P1** | `φ_c` rises with CV for the lognormal family at `μ_p` = 0.3 | Δ`φ_c`(CV 0→0.5) > 0 and > 3× seed scatter | ✅ **PASS** +0.0578 (46 sd) / +0.0623 (21 sd) |
| **P2** | the critical-state slope `s` agrees with the isotropic `φ_J` slope | within a factor of 2 of +0.1246 per unit CV | ✅ **PASS** — 0.1157 and 0.1247, within 8 % |
| **P3** | there exists a `φ` between `φ_c(CV=0)` and `φ_c(CV=0.5)` where raising CV alone takes the pack from dilatant to contractive | `p_ss` at fixed φ falls with CV and crosses the unjammed bar | ✅ **PASS on the fall** (5.4–43× at matched φ, 15–63 sd); ⚠ partial on the crossing |
| **P4** | frictionless control: `φ_c(μ_p=0, CV=0)` coincides with the isotropic `φ_J`(μ_p=0) = 0.6332, because frictionless shear jamming and isotropic jamming coincide | agree within 0.01 | ⬜ **NOT RUN** — and it is the protocol's only external anchor |
| **P5** | the CV effect is **not** explained by `d₃₂`: a uniform (zero-skew) PSD and a lognormal at matched CV have `d₃₂` within 9 % but must give different `s` | `s_lognormal/s_uniform > 1.5` | ⬜ **NOT RUN** — the isotropic analogue is already 2.7× (+0.0623 vs +0.0228) |
| **P6** | **experiment**: at fixed φ and fixed `d₁₀`, the impact regime switches at a critical CV\*, and raising φ moves CV\* **up** | a monotone CV\*(φ) with the predicted sign | *for the bench* |
| **P7** | **experiment, the discriminator**: in the high-CV spreading state the pore pressure is **positive** (overpressure), where [M26]'s R3 has none and R1 has suction | sign of a pore-pressure measurement, or the absence of cavitation bubbles | *for the bench* |

**What would falsify the mechanism outright:** `φ_c` flat or falling in CV (P1), or a
high-CV spreading drop showing pore *suction* (P7).

### 5a. ⚠ A rival that P1–P6 cannot separate, added after the referee pass

The referee named a mechanism this derivation had not, and it is degenerate with the one
here on almost every axis.

> **M2a, plain suspension rheology.** Raising CV raises the frictional jamming fraction
> `φ_m` (≈ `φ_c` for frictional spheres). At fixed `φ` that moves the suspension
> *further from jamming*, so `η ∝ (φ_m − φ)^-2` falls and the drop spreads. **No pore
> fluid, no Terzaghi, no dilatancy sign — and it uses the same measured `dφ_c/dCV` that
> the mechanism here needs.**

M2a therefore predicts the *same direction*, driven by the *same number*, and — because
`φ_m ≈ φ_c` — with the crossing in *roughly the same place*. P1 and P2 are measurements
of `φ_c(CV)` and are consumed identically by both. **P1–P6 passing would not establish
the poroelastic mechanism over M2a**, and the first version of this document did not say
so.

What separates them:

⚠ **The "cannot separate" claim is itself contingent on `φ_m ≈ φ_c`.** If a measurement
were to show `dφ_m/dCV` and `dφ_c/dCV` diverging, P1–P6 could partially discriminate the
two without E-3. As things stand they are taken to be the same shift.

| | this mechanism (M2′) | M2a, rheology |
|---|---|---|
| pore-pressure sign (**P7**) | overpressure on the spreading branch, suction on the arrested one | **makes no prediction** — a scalar viscosity divergence has no signed-pressure coupling. (Not "predicts zero": it is silent, which is still a separation, but a weaker kind.) |
| shape of the crossing | a **threshold** in `φ − φ_c`: the *sign* of a term flips | a smooth divergence/weakening; `η` is continuous through `φ_m` from below |
| transient | ⚠ **restated 2026-09-08** — the strain-clock version of this row was wrong (§9, round 2). The undrained pack **locks** and the suction grows *with strain*, so the mechanism predicts a body that **stiffens as it spreads** | none; a viscosity is a state function and has no strain history |
| drainage number | the effect switches off as `Λ` rises: locking is an undrained phenomenon, and the drained limit really does relax to `φ_c` | `Λ` does not appear at all |
| salt | strength `∝ a/6λ_D` (Eq. 17), so it **strengthens with screening** | `φ_m` is geometric; salt enters only through `σ*` |

**P7 is therefore load-bearing rather than confirmatory**, and the transient and salt
rows are two further separations that cost nothing extra to look for. The honest position
until P7 is run: the CV observation is explained by *a* packing-geometry channel, and
which of the two it is remains open.

---

## 6. The open question this derivation creates, and it is the sharp one

Equation (8) says `ID ∝ ℓ⁻²`. [M26] take `ℓ = D₀`. If instead the dilating zone is a
shear band of thickness `ξ = n·d₃₂` with `n = O(10)`, then

```
(15)   ID_local = ID(D₀) · (D₀/ξ)² = ρ D₀ U₀ /(α μ_w Δφ) · (1−φ)³/(180 φ² n²)
```

— and **every factor of `d₃₂` cancels.** The `a²` scaling that puts 800 nm silica 156×
below 10 µm cornstarch is *entirely* a consequence of choosing the drop diameter as the
drainage length. On the grain-scale choice, mean particle size does not set the regime at
all.

This is not a detail: it decides whether the author's silica *can* reach the
inertial-dilatant regime, which is the premise of the whole question. The two readings
are cleanly separable:

> **A size ladder at fixed CV, φ and U₀ separates them by 10²** over one decade of `a`:
> macroscopic drainage predicts `ID ∝ a²`, grain-scale drainage predicts `ID ∝ a⁰`.

Nothing in this derivation chooses between them, and the numbers in §1's table
(`Π` = 1.39e4 vs 0.142, a factor 9.77e4 = (D₀/ℓ)²) show the choice spans the entire
boundary: on the grain-scale reading `ID` = 7.0 and this drop sits on the INERTIAL side.

### 6a. A THIRD candidate exponent, and it is the one this analysis actually predicts

`diag_kappa_sweep.py` (added 2026-09-08; `out/kappa_sweep.txt`) integrates the transient
at the lab's own physical numbers and sweeps `κ` over ten decades:

| `κ/κ₀` | `Λ` | `p(γ=1)`/ceiling |
|---|---|---|
| 1 → 10⁵ | 2.5e-7 → 2.5e-2 | **0.4400** — flat to 1.000002× |
| 10⁶ | 0.25 | 0.377 |
| 10⁸ | 25 | 0.0057 |

So M1's absence is not just an algebraic property of Eq. (16): **the dynamics are flat
over five decades of permeability**, and they only stop being flat once `Λ` ≳ 0.25, which
is six decades away from where this system is. That is about three decades more headroom
than the static argument alone establishes.

But the same run shows the pressure is **not** size-independent, because `K_p` is:

| `a` [µm] | `κ` [m²] | `p(γ=1)` [Pa] | ratio |
|---|---|---|---|
| 0.8 | 3.23e-15 | 11.8 | 1.000 |
| 2.5 | 3.15e-14 | 36.9 | 3.125 |
| 8.0 | 3.23e-13 | 118 | 10.000 |

**`p ∝ a¹`, exactly** — because in the undrained limit `p` follows `K_p`, and for a
double-layer contact `K_p = P a/(6λ_D)` is linear in the particle size. So the size
ladder E-4 proposes now separates **three** hypotheses, not two:

| reading | scaling | over one decade of `a` |
|---|---|---|
| macroscopic Darcy drainage (`ℓ = D₀`) | `∝ a²` | 100× |
| grain-scale drainage (`ℓ ∝ a`) | `∝ a⁰` | 1× |
| **undrained, network-modulus limited (this analysis)** | **`∝ a¹`** | **10×** |

⚠ The `a¹` is specific to a **charge-stabilised** contact. For a Hertzian contact
`K_p ∝ E^(2/3)P^(1/3)` carries no particle size at all and the exponent would be 0, so
this prediction is a statement about colloidal silica and **must not** be carried across
to cornstarch or glass beads.

---

## 7. The ceiling that bounds the dilatant branch — and it is not cavitation

⚠ **This section was rebuilt on 2026-09-08 after a referee pass.** Its first version
capped equation (4) at cavitation of the pore water (≈1 atm) and reported a 459×
overshoot. That cap is the one `impact_stress/inertial_dilatant.py`'s own `ID_capped`
docstring calls **"largely withdrawn after review"**, and which `PORE_PHYSICS.md` §4.4
(2026-08-10, a month before this derivation) replaces. Building on a withdrawn cap is
the failure mode this repo already has a lesson about; the section below uses the
replacement.

Equation (4) is unbounded as `κ → 0`, and for silica it goes far past what the material
can do. At the lab's numbers (`α`=1, `Δφ`=0.02, `ℓ`=`D₀`=2.5 mm, `d₃₂`=800 nm,
`γ̇`=`U₀`/`D₀`) it gives **4.65 × 10⁷ Pa**.

⚠ **That number is a steady state the impact never reaches, and saying "Darcy builds up
to it and then gets clipped" is the wrong causal story.** Corrected after the referee's
second pass: the drainage number here is `Λ` ≈ 2.5 × 10⁻⁷, so reaching Eq. (4)'s value
would take a strain of order `1/Λ` ≈ 4 × 10⁶. Integrating the transient
(`poro_impact.py`) at the lab's own `κ`, `K_p` and `ℓ` shows the pressure never gets
within six decades of it. **The response is undrained from the first instant; Darcy
plays no part during the impact at all**, and the pressure is set by the local
dilatancy–elasticity coupling and bounded by the ceiling below. The conclusions of this
section are unaffected — the magnitude and the `κ`-independence are what they are — but
the mechanism is "there is no time to drain", not "it drains and overshoots".

The binding limit is not what the *fluid* can survive but what the *packing has left to
give*. A pack at `φ` can only dilate until it reaches `φ_c`, a volumetric strain of
`Δφ/φ`, so with `K_p` the drained network modulus (`PORE_PHYSICS.md` Eq. 4.5):

```
(16)   p_dil^max = K_p Δφ / φ
```

| quantity, 800 nm silica at φ = 0.44, Δφ = 0.02 | value |
|---|---|
| uncapped Darcy, Eq. (4), ℓ = D₀ | 4.65e7 Pa |
| **network-modulus ceiling, Eq. (16)** | **26.8 Pa** (bracket 8.9 → 268 Pa over λ_D 30 → 1 nm) |
| overshoot of the uncapped expression | **1.7 × 10⁶**, i.e. six decades |
| impact stress `ρU₀²` | 1.30e4 Pa |
| withdrawn cavitation cap | 1.01e5 Pa — **not used** |

**Three consequences, and the first two are good for the argument.**

1. **Where the ceiling binds, permeability leaves the problem entirely.** Equation (16)
   contains no `κ`. For this system it binds by six decades, so **M1 is not merely
   sub-dominant there, it is absent** — a stronger statement than the order-counting of
   Step 8, and one that does not depend on `dφ_c/dCV` at all. What survives as a PSD
   channel is `Δφ(CV)` and whatever moment of the distribution `K_p` reads; for a
   double-layer contact `K_p ∝ a`, so that is a *first* moment, not `d₃₂`, and it is
   unmeasured.
2. **The sign structure and the pole are untouched.** Equation (16) is `∝ Δφ` exactly as
   equation (4) is, so it vanishes and changes sign at the same place. Capping changes
   the *magnitude* of the dilatant branch, not *which side of it you are on*.
3. **The magnitude is small, and that is a real constraint.** 26.8 Pa against a 1.3e4 Pa
   impact stress. The meaningful comparison is not against the impact stress but against
   the **particle pressure** that already presses the contacts together, and there
   `P_conf` cancels (`PORE_PHYSICS.md` §9.4):

   ```
   (17)   p_dil^max / P_conf = (a / 6λ_D) · (Δφ/φ)
   ```

   | added salt | λ_D | ratio at φ = 0.44, Δφ = 0.02 | frustrated dilation available as strengthening? |
   |---|---|---|---|
   | nominal 0 | 30 nm | 0.20 | no |
   | 1 mM | 10 nm | 0.61 | no |
   | 100 mM | 1 nm | 6.06 | **yes** |

   ⚠ **State this at full strength, because "marginal" undersells it** (referee, round
   2). At nominal salt the frustrated-dilation channel supplies a **20 % boost to a
   confining pressure that is already there**. Nothing in this document shows that is
   enough to hold a drop **arrested and rebounding** for a whole contact time — that
   would need the frictional network essentially fully mobilised, which is a stricter
   condition than "some extra confining stress", and the Wyart–Cates activation curve
   that would decide it is not computed anywhere here (and its natural calibration,
   `si_fit`'s `σ*`, is withdrawn). **So the honest position is that at low salt the
   poroelastic channel is not shown to be sufficient on its own, and M2a and M3 may be
   doing real work there too.** What the section does establish is the *direction*: the
   route should be **strongest at high salt**, which is a prediction on an axis Paper 1
   already sweeps. Equation (17) also survives the withdrawal of `si_fit`'s `σ*`,
   because `P_conf` cancels out of it — only the absolute Pa figures in the table above
   depend on that withdrawn number.

**And one thing the first version got backwards.** With the correct ceiling at ~27 Pa,
the suction never approaches one atmosphere, so **there is no cavitation and the
"look for cavitation bubbles" discriminator is dead.** It is withdrawn from
`MECHANISM_2026-09-08.md` §6 E-3; the pressure-transducer version survives.

---

## 8. Checks that ran

```
$ python derivations/dilatancy_sign_mechanism_check.py
ALL CHECKS PASSED          (21 checks)
```

- **C1 SymPy** (4): (4), (6), the identification with [M26] Eq. 6, and the `ℓ²` scaling.
- **C2 limits** (7): `μ_w→0`, `κ→∞`, `κ→0`, `U₀→∞`, `U₀→0`, `Δφ→0`, `ℓ→0`. Each recovered.
- **C3 dimensions** (3): `pint`; `Π` dimensionless, `κ` an area, `Δp_f` a pressure.
- **C4 numeric** (2): the closed form against `impact_stress.inertial_dilatant.ID` at
  5 points spanning the range, plus a discriminating negative case (halving `Δφ` must
  move `ID` by > 1.5×; it moves it 2.03×).
- **C5 channels** (4): the two logarithmic derivatives and their `CV → 0` limits.
- **C6 the `A` channel** (1), added after the referee pass: Step 7 *asserted* that the
  third term of (11) is small and nothing checked it. Differentiating `[M26]`'s own `A`
  numerically gives `|d ln A/dCV|` = 0.213 against a critical-state term of 6.2–95, i.e.
  **0.2 %–3.4 %** of it. Assertion confirmed, and now it is a check rather than a claim.
- **C7 external anchor** (2), added after the referee pass: the same
  Reynolds–Darcy–network-modulus chain applied to **Jérôme, Vandenberghe & Forterre's
  170 µm glass beads, the one system in this comparison where the pore suction was
  MEASURED** (−10 to −100 kPa). Equation (16) predicts 4.24e5 Pa — within a decade and
  above, which is the expected direction for an affine upper bound. Its negative case:
  the same formula on 800 nm silica gives 26.8 Pa, a factor 373 outside the bead range,
  so C7 is sensitive to its inputs.

⚠ **What C4 does and does not establish.** It agrees to *exactly* zero, because it is the
same algebra typed twice — it anchors the **transcription** of Eq. 6 and says nothing
about the physics. C1.3 *derives* Eq. 6 from G1 + G3 rather than restating it, which is
worth more, but it is still internal consistency with G1–G4 rather than evidence that
G1–G4 describe reality. **C7 is the only check here that compares against data**, and it
was absent from the first version.

⚠ **C2's seven limits are weaker than their number suggests.** `Π` is a monomial with
`μ_w`, `Δφ`, `ℓ²` in the numerator and `κ`, `U₀`, `D₀` in the denominator, so *every*
limit in the list is guaranteed by that structure alone and none of them can catch a
wrong exponent. They confirm the slots, not the powers; C1's exact symbolic identity
already did that.

**Unchecked.** (i) `φ_c(CV)` itself — that is a measurement, in flight. (ii) The
one-length approximation A2; no check here can validate it, and §6 says why it matters.
(iii) Whether a *drop* under impact reaches the critical state at all within the impact
time — the critical state is a large-strain concept and the impact delivers `γ ~ O(1)`;
the DEM runs to `γ = 2` and the steadiness gate reports whether that is enough. (iv) The
coefficient `α`/`A` itself: G1–G4 fix the functional form of (4) and (6), but `A` is
[M26]'s **fitted** friction closure (Eq. 5, exponential in `tan β` with `β₀` = 44° taken
from cornstarch's repose angle) and has no Reynolds–Darcy content. It is **imported, not
derived**, and the Result Box says so.

---

## 9. Referee — findings and what changed

`physics-referee` pass, 2026-09-08, adversarial, seeing only this document, its checks,
`inertial_dilatant.py` and `PORE_PHYSICS.md`. It ran the check script and re-derived the
claims independently.

| grade | finding | response |
|---|---|---|
| **FATAL** | §7 capped the suction at **cavitation**, which `inertial_dilatant.py`'s own `ID_capped` docstring calls "largely withdrawn after review" and `PORE_PHYSICS.md` §4.4 replaces with the network-modulus ceiling. The overshoot is 1.7e6×, not 459×, and the withdrawn cap removes all particle-size dependence — the very defect it was withdrawn for. | **§7 rebuilt** on Eq. (4.5). Confirmed the referee's arithmetic: 26.8 Pa, not 101 kPa. The correction turned out to *strengthen* the argument (permeability leaves the problem entirely where the ceiling binds) and to *weaken* the magnitude (the branch is marginal at low salt, Eq. 17). The cavitation-bubble discriminator is **withdrawn**. |
| **MATERIAL** | No check compares against **data**. C4 agrees to exactly zero by construction; C1.3 is internal consistency. `PORE_PHYSICS.md` §9.3 already validates this machinery against Jérôme *et al.*'s measured suctions and the derivation never invoked it. | **C7 added**, with a discriminating negative case. 21 checks now. |
| **MATERIAL** | An alternative mechanism — plain suspension rheology through `φ_m(CV)` — is never named, uses the *same* measured quantity, predicts the *same* direction, and would not be excluded by P1–P6. | **§5a added.** P7 is re-labelled load-bearing; two further separations (the transient and the salt scaling) are named. |
| MINOR | The Result Box overstates "derived rather than quoted": `A`/`α` is imported. | Result Box amended; §8 "Unchecked" item (iv). |
| MINOR | Step 7's "the `A` term is small" was asserted and unchecked. | **C6 added.** Confirmed: 0.2 %–3.4 %. |
| MINOR | C2's seven limits are guaranteed by the monomial structure and cannot catch a wrong exponent. | Said so in §8 rather than counting them as seven independent checks. |

### Round 2, after the above were acted on

| grade | finding | response |
|---|---|---|
| **FATAL** | `poro_impact.py` eq. (iv) omitted the **elastic** half of the volumetric strain its own eq. (ii) requires. Restored, `φ − φ_c` does not decay at all: the pack **locks** and the suction grows linearly without bound (checked to strain 24 at 10.6× the ceiling). The referee derived (iv) by hand and re-integrated. **The "strain clock" — one of only three separations from M2a — was backwards.** | **Corrected**, with the flag `elastic_feedback` defaulting to the right law. The strain-clock claim is **withdrawn** from the module, §5a, `MECHANISM_2026-09-08.md` and figure M2. Replaced by the correct statement: **undrained locking / dilatancy hardening**, which is a *stronger* separation from M2a (a viscosity has no strain history), plus a new `Λ` row. What bounds the growth is a closure this model does not carry — `φ_c(σ')`, the critical-state *line* — and that is now said. |
| — | **None of the 14 tests caught it**: the one that exercised (iv) re-integrated the *same* incomplete equation in a separate scalar solver. Two anchors sharing one identity. | Replaced by a **closed-form** anchor: in the undrained limit with `α` = 1, `φ − φ_c` must be constant to machine precision and `p(t) = −SαΔφ₀ t` exactly. Both now asserted, plus a discriminating negative case that the two laws disagree, plus a test that the clock exists **only** at large `Λ`. 15 tests. |
| MATERIAL | §7's causal story ("Darcy builds up to 46.5 MPa then gets clipped") is not what the transient shows: at `Λ` ≈ 2.5e-7 there is no time for diffusion at all, and nothing approaches 46.5 MPa. | §7 amended. The conclusions (magnitude, `κ`-independence) are unaffected; the referee's own `κ` sweep over ten decades found `p_peak/ceiling` flat to four significant figures across the whole experimentally reachable range, i.e. **"M1 absent" holds with about three decades more headroom than the static argument alone gives.** |
| MATERIAL | Eq. (17)'s 0.20–0.61 is not shown sufficient for **arrest and rebound** specifically; "marginal" undersells the gap. | §7 restated at full strength. |
| MATERIAL | `MECHANISM_2026-09-08.md` §2's "M2 was graded against the wrong observable" overclaims and contradicts its own §7a two sections later — the pattern where a table row gets quoted into a Status card on its own. | §2's row now carries a forward pointer to §7a. |
| MINOR | "M2a predicts no coherent signature" is more accurately "M2a makes no prediction"; and the claim that P1–P6 cannot separate the two is itself **contingent on `φ_m ≈ φ_c`**. | Both said in §5a. |
| survives | the sign and pole structure of Eqs. (16)/(17); the `Π/2` prefactor, re-derived by hand; the Terzaghi anchor. | — |

**Survived the attack unqualified:** Eq. (4) and its assumptions (A); the Terzaghi sign
algebra (C); the order-counting in Eq. (11), which the referee re-derived numerically
(D); the pole's physical status as distinct from its provisional location (E); the
drainage-length cancellation (F), which the referee found independently corroborated by
`PORE_PHYSICS.md` §3 Eq. 3.1, derived a month earlier; and the `φ_c` incommensurability
caveat, which was already disclosed. No sign errors and no prose/artifact number
mismatches were found.

---

## 10. Open questions

1. **The drainage length (§6).** The single largest lever in the argument and the one
   nobody has measured. A size ladder decides it.
2. **`φ_c` vs `φ_J`.** P2 tests whether the isotropic slope can stand in for the
   critical-state one. If it cannot, §4's table must be rebuilt and nothing from
   `phase2` transfers to this mechanism.
3. **Where is φ during impact?** The mechanism's sign turns on `φ − φ_c` at the
   *contact*, not on the bulk 0.44. This is already the Status card's next-action #2 and
   it is now load-bearing rather than merely desirable.
4. **`φ_c` for real silica.** A frictional-sphere DEM gives `φ_c` near 0.58; [M26] use
   0.37 for cornstarch and this repo guesses 0.42 for silica. Those cannot all be on one
   axis. **The mechanism is quoted on the SHIFT `dφ_c/dCV`, never on the absolute level**,
   and the shift is what a sphere DEM transfers to a sphere experiment.
5. Does cavitation (§7) change the *sign* structure, or only cap the magnitude?
