# Does particle size decide the regime? — three versions, three withdrawals

> # 🔴 v3's VERDICT IS WITHDRAWN TOO. Referee pass 3, the same day.
>
> v3 (below) reported **UNDETERMINED**: "the dilatant stress's size exponent is one of
> {−2, −1, 0, +1} across five conventions", and in bold, "**NOTHING ON DISK PINS `ℓ`**".
> **That is over-claimed, and the file this page names as its own companion settles two of
> the three knobs.** Verified against the files before acting:
>
> 1. **`impact_stress/PORE_PHYSICS.md` §3 pins the drainage length for exactly this
>    question.** Its title is "The drainage length is a modelling choice" and its closing
>    paragraph makes the choice: *"For whether frustrated dilatancy can raise the particle
>    pressure in the spreading film, the contact length is correct: the fluid must be drawn
>    from somewhere in the film to fill the dilating region, and the film's own radius bounds
>    that distance."* `L = 0.49 D₀`, and it is a **computed** Wagner `r_c` at peak force from
>    `stress.py`, not an assumption — *"and it is the length implied by Mobaseri et al.'s
>    Eq. 3."* §3's row 2 (`ℓ = 10a`) is not a rival convention either: `size_threshold.py:138`
>    says of that very line, *"this group is SIZE-BLIND by construction and predicts no
>    threshold"* — a declared null.
> 2. **§9.1 and §9.2 already ran the capacity-vs-drainage comparison.** Lab silica, uncapped
>    `p_dil = 6.2e7 Pa` against `p_dil^max = K_p Δφ/φ = 94 Pa`: ***"The ceiling binds by six
>    decades."*** Cornstarch: 1.3e5 against 4.0e5, *"only by a factor of 3."* So the ceiling
>    that binds for the lab silica is the **capacity** one, three decades below the
>    `min(P_atm, 2γ/r_throat)` §2 used, and above it the exponent is **+1** for an EDL
>    contact and **0** for a Hertzian one — i.e. **determined, per material**, not one of
>    five. ⚠ And applying `p = +1` material-independently is the EDL-vs-Hertz error
>    `SIZE_THRESHOLD` §3 already recorded against its own v1. Committed again.
> 3. **§3's "reproduced to 0.12 %" check is a tautology.** `size_threshold.py:140` literally
>    computes `ID(ℓ=10a) = ID(ℓ=D₀)·(D₀/10a)²`. The check verified a multiplication, not that
>    `ℓ` enters squared — and its two inputs were typed literals in a page whose first line
>    says nothing is typed from memory.
> 4. **§5's bidisperse "REFUTED" rests on an unresolved number.** Paired per seed, the local
>    0.3 → 0.5 slope is **−0.0379 ± 0.0187, t = 2.03 against `t_crit` = 4.303**. `main.tex`
>    says in its own Statistics section that a flat 2 s.e. bar is a test at p = 0.184. The
>    page states −0.0379 bare, three times, and marks it REFUTED.
> 5. **§5's "1.17× locally" is the argmin of its own sweep.** The spread of the 0 → CV chord
>    slope runs **2.74 / 1.17 / 1.42 / 2.38 / 4.47** at CV = 0.1 … 0.5. The page's own `CV*`
>    is 0.30–0.35, where it is 1.42 and rising. The local max−min difference is
>    **+0.0100 ± 0.0108, t = 0.93** — it establishes nothing, and it was used to call an
>    "endpoint artefact" a result that passes three gates in `main.tex`. **That claim is
>    retracted and nothing in `paper/` should change on account of it.**
>
> ### What is actually established, after three passes
>
> * `σ_in/P_f = A·ID` exactly, and the force-vs-stress distinction does not change the
>   exponent: `F_d = A π r_c² · P_f`, so the `r_c²` cancels. [M26]'s number **is** the
>   A-vs-B discriminator and this page never had a new one.
> * For the lab silica the dilatant stress is **capacity-limited**, and the capacity ceiling
>   carries the particle size through the contact law: `K_p ∝ a` for a double layer, nothing
>   for Hertz. That is a statement per material, and it is in `PORE_PHYSICS` §9.
> * What is genuinely open is narrower than v3 claimed: **`φ_c`** for this material (the only
>   value naming the model's own quantity is a placeholder), **which contact law it has**,
>   and **`λ_D`**, which `PORE_PHYSICS` §5.2 says "spans a decade" and "is not measured".
> * ⚠ Unexamined by all three versions: `poroelastic_time_ratio` returns **τ_χ/τ_imp = 1.058**
>   for the lab silica, against 1.5e-3 and 0.021 for cornstarch. The lab's own particle sits
>   on the incompressible-pore-fluid boundary that every row of §3 assumes.
>
> **Everything below is left as the record of v3, with its errors marked here rather than
> edited out. Do not quote §3's "3 decades per size decade" or §5's 1.17× and −0.0379.**

---

# v3, as written and now withdrawn: UNDETERMINED, and this page says exactly why

**2026-09-19.** Every number below is printed by `onset_size_degeneracy_check.py` into
`onset_size_degeneracy_check.txt` (rc 0, ALL CHECKS PASSED); nothing is typed from memory.
Companions: `../SIZE_THRESHOLD_2026-09-08.md` (which this extends),
`../../impact_stress/inertial_dilatant.py`, `../../impact_stress/dilatancy.py`,
`../../impact_stress/PORE_PHYSICS.md`.

> ## VERDICT: **UNDETERMINED.** Two headlines were written today and both were withdrawn.
>
> | version | claim | withdrawn because |
> |---|---|---|
> | **v1** | `ID` and `σ_imp/σ*` carry the same power of `d`, so no size ladder can separate the frictional and inertial-dilatant regimes | it rested on "no admissible impact-stress gauge carries a particle size". **Three of six do.** |
> | **v2** | the degenerate pair is `(σ*, P_f)`, the two grain-scale stresses, so a size ladder cannot separate frictional jamming from frustrated dilatancy | **(a)** `σ*` is the **threshold both routes must cross**, not one route's driving stress, so `P_f/σ*` is the dilatant branch's *activation number* — a statement about **one** mechanism. **(b)** `P_f`'s own size exponent is not −2; the conventions that fix it are unpinned and span **−2, −1, 0, +1**. |
>
> **The second withdrawal is the one worth keeping: the repair committed the same error as
> the thing it repaired.** v1 chose a convention and called the question settled; v2 chose a
> different convention and called the question settled. Both are one row of the table in §3.
>
> ### What is actually true
>
> 1. **A size ladder is NOT blind to the mechanism.** The A-vs-B discriminator is
>    driving-stress against driving-stress, `σ_in/P_f`, and on [M26]'s own convention its
>    exponent is **+2.00** — ~100× per size decade. Both v1 and v2 asserted the opposite,
>    four lines below printing that exponent. And `σ_in/P_f` **is** `A·ID`: [M26]'s number
>    already is the discriminator, and this page adds no new one.
> 2. **The dilatant stress's size exponent is undetermined over 3 decades per size
>    decade.** Three conventions are unpinned on disk — the drainage length `ℓ`, the
>    ceiling regime, and whether capacity or drainage binds. §3.
> 3. 🔴 **The Darcy branch is not valid at any parameter point either version used.**
>    `dilatant_pressure`'s own docstring says *"Compare it against `p_dilatant_ceiling`
>    before believing it."* Neither version did. `P_f/ceiling` = **459 / 77 / 23 / 1.6** at
>    the four points. Above the ceiling the exponent is 0 (cavitation) or −1 (capillary
>    entry), not −2. §2.
> 4. **The double layer adds more size dependence at 800 nm than the whole effect being
>    argued about**, and v2's salt escape clause fails on v2's own model once the condition
>    is stated on the *exponent* rather than the *value*. §4.
> 5. **The CV lever's sign is the only part that transfers**, and even that is a statement
>    about the chord: the bidisperse family's local slope above CV = 0.3 is **−0.0379**. §5.
>
> ### What would settle it — three measurements, none of them algebra
>
> | # | question | how |
> |---|---|---|
> | 1 | does the suction relax over the **drop scale** or the **grain scale**? `ℓ` enters squared, so this alone is 4 decades of `ID` | a pore-pressure probe, or the `U₀`-dependence of the force plateau |
> | 2 | does the pore fluid **cavitate** during impact? | §2 says both materials are 1.6–459× above 1 atm of suction on the Darcy formula |
> | 3 | what is **`φ_c` for the lab silica**? | the one value on record that names the model's own `φ_c` is flagged a placeholder (§5.3) |
>
> **This page can decide** which conventions the answer depends on, and by how much.
> **It cannot decide** the answer.

---

## 1. The algebra — correct, and an identity given the definitions

With `κ = C d²`, `C = (1−φ)³/(180φ²) > 0` (checked numerically on 20001 points):

```
(1)   σ*    = B d^(−α)                          α = 1.9 ± 0.1, MEASURED (§6)
(2)   P_f   ~ (μ/κ) Δφ γ̇ ℓ²                     the Reynolds–Darcy suction
(3)   σ_in  = ρ U₀²                             the inertial impact stress
(4)   ID    = ρ κ U₀ / (A μ Δφ D₀)              [M26] Eq. 6
```

Three ratios, and **what each one asks** — the distinction v2 lost:

| ratio | exponent in `d` (at `ℓ ~ D₀`) | the question it asks |
|---|---|---|
| `σ_in/σ*` | `α` = +1.90 | is the **imposed** stress above the frictional onset? — *activation, route A* |
| `P_f/σ*` | `α−2` = −0.10 | is the **suction** above the frictional onset? — *activation, route B* |
| **`σ_in/P_f`** | **+2.00** | **does inertia beat the suction? — the A-vs-B discriminator** |

> 🔴 **v2 read the middle row as "the two mechanisms are size-degenerate."** It is not.
> `σ*` is the threshold **both** routes cross; the middle row is route B's activation
> number and says nothing about route A. The discriminator is the last row, and its
> exponent is +2 — a size ladder separates the mechanisms at ~100× per size decade.
>
> And `σ_in/P_f = A·ID` exactly (checked symbolically). So the discriminator is [M26]'s own
> number, and the only thing this page has added to it is a list of the conventions it
> depends on.

All four size exponents and all three ratio exponents are reproduced symbolically, and the
`α = 2` exact cancellation holds. **These are identities given (1)–(4), not results.** The
physics is entirely in which of (1)–(4) is the right form — which is §3.

---

## 2. 🔴 The Darcy branch is not valid where either version used it

`dilatant_pressure` returns a positive magnitude and its docstring instructs the caller to
compare it with `p_dilatant_ceiling = min(P_atm, 2γ/r_throat)`. Neither v1 nor v2 did.

| point | `P_f` [Pa] | ceiling | `p_cap` | `P_f`/ceiling |
|---|---|---|---|---|
| silica, 3.0 m/s | 4.651e+07 | 1.013e+05 | 2.323e+06 | **459.1** |
| silica, 0.5 m/s | 7.751e+06 | 1.013e+05 | 2.323e+06 | **76.5** |
| cornstarch, 7.2 m/s | 2.343e+06 | 1.013e+05 | 1.858e+05 | **23.1** |
| cornstarch, 0.5 m/s | 1.627e+05 | 1.013e+05 | 1.858e+05 | **1.6** |

**All four are above the ceiling.** Above it the suction is set by cavitation of the pore
fluid (exponent **0**) or by the interface entering a pore (`2γ/r_throat`, exponent **−1**),
not by Darcy (−2). The cornstarch/0.5 m/s point sits at 1.6×, i.e. the boundary runs
straight through the parameter range in question.

⚠ The counter-argument is real and is on record: `ID_capped`'s docstring lists four
objections to applying the cap. That makes the cap *contestable*; it does not make the
regime question skippable. Both versions of this page claimed a size exponent for a
material their own numbers place outside the formula's stated validity.

Related: v2 wrote "the suction is 37× to 3600× larger than `ρU₀²`". That quantity is
algebraically `1/(A·ID)` — the very ratio §6.1(A3) declares carries an undetermined
**O(10)–O(100)** prefactor. At cornstarch it is 36.7×, *inside* that range, so its direction
is not determined there.

---

## 3. The convention scan — the result

`P_f ~ (μ/κ)Δφ γ̇ ℓ²` with `κ ~ d²`, so the exponent is `2·(d ln ℓ/d ln d) − 2`, and
**nothing on disk pins `ℓ`.** `SIZE_THRESHOLD` §1 tabulates **both** readings and this page
reproduces the factor between them exactly: `ID(ℓ = D₀) = 7.19e−05` and `ID(ℓ = 10a) = 7.03`,
against a predicted `7.19e−05 × (D₀/10a)² = 7.021` — **0.12 %**, confirming `ℓ` enters
squared and that both readings are on record.

| convention for the dilatant stress | `p` | vs `σ*`: `d^` | vs `σ_in`: `d^` |
|---|---|---|---|
| Darcy, `ℓ ~ D₀` (drop-scale; [M26]'s choice, **v2's pick**) | −2 | −0.10 | +2.00 |
| Darcy, `ℓ ~ 10a` (grain-scale; `SIZE_THRESHOLD` §1) | 0 | +1.90 | 0.00 |
| cavitation-limited (§2 says it binds) | 0 | +1.90 | 0.00 |
| capillary-entry limited, `2γ/r_throat` (§2 says it binds too) | −1 | +0.90 | +1.00 |
| **capacity** criterion, `p_dil^max/P_conf ~ (a/6λ_D)(Δφ/φ)` (`PORE_PHYSICS` 9.1) | **+1** | +2.90 | −1.00 |

**Span: 3 decades per size decade.** The last row is the one `SIZE_THRESHOLD` §3 already
names — *"capacity as a¹, drainage as a²"* — and that section withdrew the window the two
bracket as unquotable, for the same reason.

Note the second column also flips sign across the table: `P_f/σ*` is −0.10, +1.90, +0.90 or
+2.90. **v2's degeneracy exists on exactly one of five rows**, and the row where the two
grain-scale stresses track each other is also the row §2 says is out of validity.

⚠ Two of the five rows are not independent constructions: `μU₀D₀/κ` (which v2 listed
separately) **is** `P_f` with the `Δφ` factor removed. v2's "three of six gauges carry a
size" double-counted; the honest count is five distinct constructions, of which three carry
a non-zero exponent.

---

## 4. The double layer — a condition on the exponent, which is stricter than v2's

`κ_eff/κ = (1−x)^n`, `x = λ_D/r_throat`, `n = 2` (immobile shell) or `4` (Poiseuille
throat); both readings are **estimates made here**. The throat convention
`r_throat = 0.155R = 0.0775d` is re-derived by hand (centroid-to-centre `2R/√3 = 1.1547R`,
inscribed sphere `0.1547R`) and matches `PORE_PHYSICS` §5.2 and `inertial_dilatant.r_throat`.

⚠ **v2 required the *value* correction to be small (`x ≤ 0.1` → 2–4 %). The claim is about
an *exponent*.** `d ln(1−x)^n/d ln d = n x/(1−x)`, which at `x = 0.1` is **0.222** (n = 2) or
**0.444** (n = 4) — 2.2× to 4.4× the residual `|α−2| = 0.10` the condition exists to protect.

Requiring the **extra exponent ≤ 0.10**:

| `λ_D` [nm] | added salt [mM] | `d_min` (n = 2) [µm] | `d_min` (n = 4) [µm] |
|---|---|---|---|
| 30 | 0 | 8.13 | 15.87 |
| 95 | 0 | 25.74 | 50.26 |
| 10 | 1 | **2.71** | 5.29 |
| 3 | 10 | 0.81 | 1.59 |
| 1 | 100 | 0.27 | 0.53 |

> 🔴 **v2's escape clause — "or the lab's own 800 nm particle at 1–10 mM added salt" —
> fails on v2's own model.** At 1 mM the particle needs `d > 2.71 µm`. At 10 mM it clears
> only on the mildest reading (0.81 µm against 0.80 µm), i.e. by 1 %.
>
> And `λ_D = 30 nm` is itself the mild choice: `PORE_PHYSICS` §5.2 calls it a compromise
> against its own **primary figure of 95 nm**, says the value "spans a decade" and is "not
> measured". At 95 nm the condition excludes the 10 µm cornstarch as well.

What survives, and it is robust across both readings and both `λ_D`: **at 800 nm and zero
added salt the double layer dominates the residual exponent**, so any size statement about
the lab's own particle at its own chemistry is contaminated before the conventions of §3
are even considered.

---

## 5. The CV lever — the sign survives; nothing quantitative does

### 5.1 Read each family from its own ladder

`φ_J(CV)` per family, mean over 3 seeds, `μ_p = 0.3` (`out/phase2_fits.csv`):

| family | 0.0 | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 |
|---|---|---|---|---|---|---|
| lognormal | 0.57453 | 0.57811 | 0.58805 | 0.60399 | 0.61970 | 0.63661 |
| uniform | 0.57514 | 0.57896 | 0.58665 | 0.59591 | 0.59744 | 0.59934 |
| bidisperse | 0.57514 | 0.57653 | 0.58717 | 0.59661 | 0.59409 | 0.58902 |

| family | chord 0→0.5 | local 0→0.2 | local 0.3→0.5 | local/chord |
|---|---|---|---|---|
| lognormal | 0.1241 | 0.0676 | **+0.1631** | 0.544 |
| uniform | 0.0484 | 0.0576 | +0.0171 | **1.189** |
| bidisperse | 0.0278 | 0.0602 | **−0.0379** | **2.167** |

> 🔴 **v2 rescaled every family's chord by the lognormal local/chord ratio (0.544)** — i.e.
> assumed family-independent convexity. The two non-lognormal families are **concave** there
> (ratio > 1), so the correction had the **wrong sign** for two of three. Corrected `CV*` at
> `Δφ₀ = 0.02`: lognormal 0.296 (v2: 0.296), uniform **0.348** (v2: 0.759), bidisperse
> **0.332** (v2: 1.323). All three land at CV ≈ 0.30–0.35, not the spread v2 reported.
>
> ⚠ **And the 4.47× family spread is an endpoint artefact where the lever is asked about.**
> Over the chord it is 4.47× (the manuscript's headline, correctly, for the chord);
> **locally at CV ≈ 0.2 it is 1.17×.** The three families are nearly indistinguishable
> exactly where `CV*` sits. This does not contradict the manuscript — which quotes the
> endpoint definition and brackets it (2.89× peak-defined, 5.51× inverse-variance) — but it
> does mean §5's lever cannot be given a family-dependent magnitude from that number.
>
> 🔴 **The sign is not family-universal either.** The bidisperse family is non-monotone
> (`main.tex`: it peaks near CV = 0.3) and its local slope above 0.3 is **−0.0379**. v2's
> "`dφ_c/dCV > 0` on every family" is a statement about the **chord** only, and v2's quoted
> bidisperse `CV* = 1.262` lies past a turning point where no crossing exists.

### 5.2 The channel is lognormal-only

`phase8`'s critical-state ladder contains **`family = "lognormal"` and nothing else**
(checked on the shards). Any `s = dφ_c/dCV` taken from it is a one-family number. ⚠ v2
quoted it as "measured by this project" without the label, and additionally hardcoded three
`φ_c` values that are **not in `out/phase8_verdict.txt`** (its CV = 0.2 entry was out by
0.0010 = **10 seed s.d.**) and a "0.089–0.125" bracket matching none of that file's four
gate sets. All of that is now read from the files.

### 5.3 🔴 The model's `φ_c` is unmeasured for the lab silica

v2 listed four values and read "3 of 4 close the dilatant branch" as a sign ambiguity in one
quantity. **They are not one quantity**, and `PORE_PHYSICS` §9.5 and `inertial_dilatant.py`
both already say so:

| value | what it is | admissible here? |
|---|---|---|
| **0.4200** | drop-impact solidification fronts | the model's `φ_c` — and the module calls it a **PLACEHOLDER** |
| 0.4700 | Shah et al., first localised solid | an *observation of solidity* under their loading, not a dilatancy onset |
| 0.4900 | Shah et al., bulk | same, different criterion |
| 0.5176 | repo Wyart–Cates `φ_M` | 🔴 a **viscosity-divergence** fraction, and part of the **withdrawn `si_fit`** — the same page's result box says `si_fit` is withdrawn and then counted it as a vote |

The honest statement needs no vote: **the one value that names the model's own `φ_c` is
flagged a placeholder**, so the lever has no measured input for this material.

### 5.4 What survives

`Δφ` enters the dilatant stress linearly and `ID` as `1/Δφ`, and CV raises `φ_c`, so
widening the size distribution at fixed `φ` moves the drop toward the inertial side. **That
sign matches the author's bench observation** (2026-08-31: *"at similar particle size,
raising CV moves the impacting drop from purely dilatant to spreading and inertial-dominated
for some time"*), and the chord slope is positive on all three families. Everything
quantitative — the magnitude, the family dependence, the CV at which the branch closes —
needs a measured `φ_c`, which does not exist.

---

## 6. Inputs and their provenance

### 6.1 Assumptions

- **(A1)** `κ = (1−φ)³d²/(180φ²)`, Carman–Kozeny, `d` the **Sauter** diameter — re-derived
  from `κ = ε³/(k_K S_v²(1−ε)²)` with `k_K = 5`, `S_v = 6/d₃₂`, `ε = 1−φ`; the `180 = 5·36`
  and the Sauter identification are exact.
- **(A2)** `σ* = F*/(6πd²)`; the `6π` is a convention. **Only the exponent is used.**
- **(A3)** `P_f` uses [M26]'s own substitutions `r_c ~ D₀`, `|dh/dt| ~ U₀` — order-unity
  stand-ins for quantities that vary by factors of several during one impact, so `P_f` and
  `ID` share an undetermined **O(10)–O(100)** prefactor. Being a *constant*, it drops out of
  every log-derivative here; it does **not** drop out of §2's ratios.
- **(A4)** Wyart–Cates: `σ*` is independent of `φ`.
- **(A5)** `A` ([M26] Eq. 5) carries no particle size. ⚠ It **does** carry `Δφ` through
  `tanβ = tanβ₀ + tanΨ`, `tanΨ ∝ Δφ`, so `d ln ID/d ln Δφ` is not exactly −1. v1 and v2
  asserted −1. This was invisible to the symbolic check (which carries `A` as a free symbol)
  and to the numeric check (which finite-differenced only in `d`, the one knob where `A` is
  inert). **A finite difference checks only the knobs you vary.**

### 6.2 `α = 1.9 ± 0.1` — the one external anchor, and it is an open problem

Guy, Hermes & Poon, *PRL* **115**, 088304 (2015), doi:10.1103/PhysRevLett.115.088304,
arXiv:1506.09041, full text fetched: `σ̂* = σ*d³/k_BT = f*d/k_BT`, and "`σ* ≃ Bd^α` with
`B = 2×10⁻¹¹ J·m^(−1.1)` and **α = −1.9(1)**", PMMA 268–3770 nm. The units of `B` are
dimensionally consistent with `σ* = B d^(−1.9)` (Pa·m^1.9 = J·m^−1.1) — an independent check
on how the exponent was read out. Independently on record in `dst_fit`
(`src/dst_fit/io_literature.py:41`, `scripts/gate_g12.py:112`).

> The naive geometric expectation is **α = 1**: Derjaguin gives `F ∝ R` generically,
> independent of whether the flat-plate law is electrostatic, van der Waals or a steric
> brush. Guy et al. resolve this only by positing a size-independent `f*` (3 k_BT/nm) and
> write that its microscopic origin "is at present unclear".

Four caveats: **(i)** an erratum, *PRL* **116**, 059901 (2016), was **not read** (paywalled,
two attempts) — `[UNVERIFIED against the erratum]`; **(ii)** Guy's system is **sterically**
stabilised and states salt does not move its `σ*`, while the lab's silica is **charge**
stabilised — the transfer rests on Mewis & Wagner (2011), which is **`[UNVERIFIED]`** (not
resolved to a page or equation) and is the single citation carrying `α` to this page's own
material; **(iii)** only one paper in the searched literature has a size series at all (Lin
2015, Royer 2016, Comtet 2017, Cwalina & Wagner 2014 all resolved, none reports one — an
abstract-level check, so CANNOT TELL rather than a hard negative); **(iv)** `α` is applied
here to 10 µm, 2.65× outside the measured range.

### 6.3 Prior art

The surviving algebraic content is Wyart–Cates activation applied to a Reynolds–Darcy
suction in the Brown–Jaeger frustrated-dilatancy picture with a Carman–Kozeny closure — **a
rediscovery in a new combination, not a new law.**

---

## 7. The checks

`python onset_size_degeneracy_check.py` → rc **0**, ALL CHECKS PASSED.

- **§S1 symbolic** — every size exponent and all three ratio exponents; the `α = 2`
  cancellation; and `σ_in/P_f = A·ID` exactly. Labelled identities, not results.
- **§S2 dimensions (pint)** — `κ` an area; `σ*`, `P_f`, `ρU₀²` in Pa; `ID` dimensionless;
  plus the `B`-units check above.
- **§S3 limits** — `d → 0` at α = 1.9 and at α = 2; `Δφ → 0⁺` for `P_f` and `ID`; `μ → 0`;
  `φ → 1⁻`. ⚠ v2's "finite nonzero constant" check tested only `d ∉ free_symbols`, which
  `oo` and `0` also pass; it now excludes both.
- **§S4 shipped code** — size exponents +2 and −2 measured by finite difference on
  `ID()` and `dilatant_pressure()`, **and the ceiling comparison of §2**.
- **§S5** the drainage-length factor reproduced to 0.12 %, and the convention scan.
- **§S6** the double-layer exponent condition.
- **§S7** per-family ladders and slopes read from `out/phase2_fits.csv`; the phase8 family
  set read from the shards.

⚠ **Two honest limits on these checks.** (a) The §S1 exponent checks differentiate
expressions defined one line above — they verify the algebra, not the physics. (b) The
convention table in §S5 is a hand-typed list, so the check over it ("the candidates span 3
decades") can only fail if the list is typed wrong; its defence is that every row is sourced
to a file named beside it. **The page's one external anchor, `α`, is named in §6.2 and is
not asserted by any check.**

---

## 8. Referee

Two passes, `physics-referee`, fresh context each, against stated acceptance criteria.

- **Pass 1: 4 FATAL, 7 MATERIAL, 9 MINOR.** Withdrew v1. Found the gauge enumeration
  incomplete, the double layer excluded, the `s` bracket unlabelled by family, and the `Δφ`
  sign filed as a magnitude caveat. It also **cleared** one worry I had flagged: re-reading
  the phase8 pairs at *matched pressure* gives a CV 0 → 0.5 shift of +0.0488, inside the
  `p → 0` root's own gate spread, so the critical-state-vs-jamming object difference does
  not change the sign or materially the size.
- **Pass 2: 3 FATAL, 6 MATERIAL, 4 MINOR.** Withdrew v2. §1's threshold-vs-driving-stress
  point, §3's further conventions, §5.1's rescaling, §2's ceiling, §4's exponent condition,
  §5.2's hardcoded numbers and §5.3's incommensurable `φ_c` are all its findings, verified
  here against the files before being acted on. It also confirmed, by hand: the Terzaghi
  sign convention is coherent (`σ' = σ − u` with `u < 0` gives total **plus** suction, and
  `dilatant_pressure` returns a positive magnitude, so like-signed stresses are compared);
  `r_throat = 0.155R`; and that `Δφ_J` and 4.47× reproduce `paper/numbers.tex` exactly.

⚠ **This third version has not been refereed.** Its verdict is a negative, which is the
cheaper thing to get right, but the enumeration in §3 is the claim most likely to be
incomplete — that is exactly the failure mode that killed v1.
