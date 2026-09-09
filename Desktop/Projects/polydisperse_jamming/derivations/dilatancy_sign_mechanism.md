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
> With the drainage length ℓ = D₀ this is **exactly 1/ID of [M26] Eq. 6** — derived
> here from Reynolds + Darcy rather than quoted (check C1.3). Π ≫ 1 is undrained,
> Π ≪ 1 is drained.
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

⚠ **`s` here is the ISOTROPIC jamming slope, not the critical-state one.** They are
different quantities and the substitution is an assumption, not a measurement. Closing it
is what `phase8_criticalstate.py` exists for, and until it returns, **every number in
this table is provisional and the `CV* = 0.16` must not be quoted.** What is *not*
provisional is Step 8: the ordering of the two channels at small CV holds for any
`s ≠ 0`.

---

## 5. Predictions, pre-registered

Registered before `phase8_criticalstate.py` returned any `φ_c`. Bars are stated so they
can fail.

| # | prediction | bar | status |
|---|---|---|---|
| **P1** | `φ_c` rises with CV for the lognormal family at `μ_p` = 0.3 | Δ`φ_c`(CV 0→0.5) > 0 and > 3× seed scatter | *pending* |
| **P2** | the critical-state slope `s` agrees with the isotropic `φ_J` slope | within a factor of 2 of +0.1246 per unit CV | *pending* — **a failure here does not damage the mechanism, it damages §4's placeholder** |
| **P3** | there exists a `φ` between `φ_c(CV=0)` and `φ_c(CV=0.5)` where raising CV alone takes the pack from dilatant to contractive | `p_ss` at fixed φ falls with CV and crosses the unjammed bar | *pending* |
| **P4** | frictionless control: `φ_c(μ_p=0, CV=0)` coincides with the isotropic `φ_J`(μ_p=0) = 0.6332, because frictionless shear jamming and isotropic jamming coincide | agree within 0.01 | *pending* |
| **P5** | the CV effect is **not** explained by `d₃₂`: a uniform (zero-skew) PSD and a lognormal at matched CV have `d₃₂` within 9 % but must give different `s` | `s_lognormal/s_uniform > 1.5` | *pending* — the isotropic analogue is already 2.7× (+0.0623 vs +0.0228) |
| **P6** | **experiment**: at fixed φ and fixed `d₁₀`, the impact regime switches at a critical CV\*, and raising φ moves CV\* **up** | a monotone CV\*(φ) with the predicted sign | *for the bench* |
| **P7** | **experiment, the discriminator**: in the high-CV spreading state the pore pressure is **positive** (overpressure), where [M26]'s R3 has none and R1 has suction | sign of a pore-pressure measurement, or the absence of cavitation bubbles | *for the bench* |

**What would falsify the mechanism outright:** `φ_c` flat or falling in CV (P1), or a
high-CV spreading drop showing pore *suction* (P7).

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

---

## 7. Two ceilings that bound the dilatant branch

Equation (4) is unbounded as `κ → 0`, and for silica it goes past what the material can
sustain. At the lab's numbers (`α`=1, `Δφ`=0.02, `ℓ`=`D₀`=2.5 mm, `d₃₂`=800 nm,
`γ̇`=`U₀`/`D₀`) check C3.3 gives

```
   |Δp_f| = 4.7e4 kPa = 47 MPa
```

against two ceilings already implemented in `impact_stress.inertial_dilatant`:

- **cavitation** of the pore water at ≈ 1 atm of suction (101 kPa)
- **capillary entry** of air through a pore throat, `2γ/r_throat` = 2.3 MPa for 800 nm

So the dilatant branch **saturates at cavitation, by a factor of ~460**, and `ID` as an
uncapped ratio overstates the dilatant force there. `ID_capped` exists in that module for
exactly this reason. Consequence for this derivation: **R1's strength is set by the
cavitation ceiling, not by (4)** — which also means R1 is only weakly dependent on `κ`,
i.e. only weakly dependent on the mean size, once saturated. The *contractive* branch R2
has no such ceiling (an overpressure is bounded only by the total stress), so the R1↔R2
asymmetry is real and is another thing to look for at the bench.

---

## 8. Checks that ran

```
$ python derivations/dilatancy_sign_mechanism_check.py
ALL CHECKS PASSED          (17 checks)
```

- **C1 SymPy** (4): (4), (6), the identification with [M26] Eq. 6, and the `ℓ²` scaling.
- **C2 limits** (7): `μ_w→0`, `κ→∞`, `κ→0`, `U₀→∞`, `U₀→0`, `Δφ→0`, `ℓ→0`. Each recovered.
- **C3 dimensions** (3): `pint`; `Π` dimensionless, `κ` an area, `Δp_f` a pressure.
- **C4 numeric** (2): the closed form against `impact_stress.inertial_dilatant.ID` at
  5 points spanning the range, plus a discriminating negative case (halving `Δφ` must
  move `ID` by > 1.5×; it moves it 2.03×).
- **C5 channels** (4): the two logarithmic derivatives and their `CV → 0` limits.

⚠ **What C4 does and does not establish.** It agrees to *exactly* zero, because it is the
same algebra typed twice — it anchors the **transcription** of Eq. 6 and says nothing
about the physics. The check that carries physical weight is **C1.3**, which *derives*
Eq. 6 from G1 + G3 rather than restating it, and in doing so identifies their `A` with
the dilatancy coefficient `α` and their implicit length with `D₀`.

**Unchecked.** (i) `φ_c(CV)` itself — that is a measurement, in flight. (ii) The
one-length approximation A2; no check here can validate it, and §6 says why it matters.
(iii) Whether a *drop* under impact reaches the critical state at all within the impact
time — the critical state is a large-strain concept and the impact delivers `γ ~ O(1)`;
the DEM runs to `γ = 2` and the steadiness gate reports whether that is enough.

---

## 9. Referee

`physics-referee` pass: see §9 of `MECHANISM_2026-09-08.md` for the findings and what
was changed in response.

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
