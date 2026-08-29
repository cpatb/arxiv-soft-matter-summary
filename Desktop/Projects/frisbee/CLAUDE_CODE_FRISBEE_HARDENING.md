# CLAUDE_CODE_FRISBEE_HARDENING.md

**Target repo:** `C:\Users\compl\Desktop\Projects\frisbee\`
**Generated:** 2026-07-26, from a five-panel expert screen of the current state
**Status of the work being reviewed:** `FRISBEE_STATUS.md` (pipeline), `FRISBEE_GENERALIZATION.md` (the generalization memo), `Project1_Frisbee_LoyolaWind_PAPER.REVIEWED.tex` (22-page preprint, "POST-READY")

---

## 0. How to use this document

This is a self-contained execution brief. You do not need the panel reports; everything
load-bearing is restated here, including the physics.

Five expert panels reviewed the project independently: an urban-wind-engineering CFD
referee, a turbulence/dispersion theorist, a flight-dynamics and sports-engineering
auditor, a UQ / V&V / design-of-experiments referee, and a research-strategy reviewer.
They were given the same artifacts and did not see each other's findings.

**Work strictly in package order.** The packages are ordered by (information gained) /
(cost), and three of them are gates: if a gate fails, stop and report rather than
continuing. WP0 through WP2 cost almost no compute and can invalidate everything
downstream, which is exactly why they come first.

**House rule inherited from this repo's audit convention:** conservative-fix /
flag-don't-edit. Fix code defects and provenance defects directly. Do **not** silently
change physics claims in the manuscript — mark them, produce the corrected number, and
leave the author the final call. Every retraction in WP4 is written as a specific
edit to make, but flag each one in the commit message rather than burying it.

**Do not post, submit, or push the preprint anywhere until WP4 is complete.** Two
panels independently found statements in the current `.tex` that contradict its own
abstract.

---

## 1. The headline: what the panels agreed on

Ranked by how much each should change what happens next. Items 1–3 are unanimous or
near-unanimous across independent panels that could not see each other's work — treat
that convergence as strong evidence.

### C1. The gust generator is in the wrong asymptotic limit. (Found independently by 4 of 5 panels.)

`P4_Coupling/synthetic_turbulence.py` sets the gust correlation time to the RANS
turbulence decay time `T_L = k/ε`, evaluated once at the launch point. Measured at this
site, `T_L = 10.5–18.0 s`. Flight time is `T_f ≈ 2.5 s`.

The disc therefore experiences a **frozen random offset**, not turbulence. Measured
autocorrelation at 4 s lag is 0.74; the within-flight standard deviation of the gust
signal is 0.28 against an ensemble standard deviation of 1.0, i.e. ~90% of the
"fluctuation" is a per-throw constant. Because the wind couples through `V_rel = v − U`
and Galilean invariance is exact in this code (verified: residual 2.4e-7 m), **a
spatially-uniform quasi-constant `u'` is literally a randomly drawn mean wind.**

This is self-defeating, because the paper's central conceptual claim is exactly that
*"mean wind is a bias the thrower can aim off; turbulence is variance they cannot."*
The code implements turbulence as a randomly-drawn mean wind.

`FRISBEE_GENERALIZATION.md` §2 contains the correct argument and the code contradicts
it. The disc supplies its own advection (`V_d/U ≈ 8–10`), so the field is frozen in
*space* and the relevant decorrelation is over **path arclength**, not wall-clock time:

```
τ_c  =  L_int / V_d  ≈  6 m / 14 m/s  ≈  0.43 s        (not k/ε ≈ 10–18 s)
```

A factor of ~30 in the one number the entire averaging argument turns on.

**Measured consequence** (60 throws/cell, σ_u = 0.40 m/s): σ_landing = 2.62 m at
T_L = 12 s vs 1.01 m at τ_c = 0.36 s. **The current pipeline over-predicts
turbulence-induced dispersion by 2.3–2.6×**, which is `√N_eddy` for N ≈ 4–6 — exactly
what the memo's own scaling predicts. Two panels found this from opposite directions
(one from the code, one from the theory) and got the same factor. That coincidence is
the strongest single piece of evidence in this review that the underlying physics is
right and the implementation is wrong.

Every landing-scatter, success-rate and Q number in `FRISBEE_STATUS.md` §Stage 5–6
inherits this. The heading *ranking* may survive (the bias is common-mode); the absolute
values do not.

Secondary defects in the same file:
- `T_L` is fixed at the launch point (`monte_carlo.py:85`, `ref_pos=launch`) and never
  updated along the trajectory — so the height dependence of the integral scale, which
  is half of the "throw low wins twice" argument, is **not in the simulation at all**.
- σ is isotropic (`synthetic_turbulence.py:163`). Real ABL turbulence has
  σ_w ≈ 0.5 σ_u. Since `w'` drives α and hence the pitching moment, isotropy inflates
  the moment channel by ~2× in variance.
- `T_L_MAX = 20.0 s` silently clips in low-ε cells and the realised per-cell `T_L` is
  never recorded, so no one can reconstruct what the gust model actually did.
- The module docstring claims `T_L` is "O(1 s) ... comparable to the flight time." Its
  own artifact (`vertical_structure.json`) says 10.5–18.0 s. Stale and misleading.

### C2. The "shelter is not protection" headline is confounded, and the confound is in the two directions the claim rests on.

`FRISBEE_GENERALIZATION.md` §1 divides all three production cases by a single control
(`abl_control.npz`, sampled from `dir_090_empty`). **I verified the case metadata
directly.** The picture is more specific than the panels stated, and this correction
matters:

| case | `z0_inlet` | `z0_ground` | `u*` | `k_inlet` | matched to the control? |
|---|---|---|---|---|---|
| `dir_090_empty` (the control) | 0.001 | 0.03 | 0.2226 | 0.1651 | — |
| `dir_090` | 0.001 | 0.03 | 0.2226 | 0.1651 | **yes, exactly** |
| `dir_210` | **0.5** | 0.03 | **0.6733** | **1.5113** | **no** |
| `dir_270` | **0.5** | 0.03 | **0.6733** | **1.5113** | **no** |

So — correcting a claim one panel made about ground roughness, which is in fact matched —
the real confound is the **inlet** roughness:

```
σ_u,inlet ratio (urban / control) = sqrt(1.5113 / 0.1651) = 3.03
```

The 210° and 270° production cases are fed inflow carrying **3.03× the fluctuation** of
the control before they meet a single building. Divide the memo's σ_u ratios by that:

| θ | memo's σ_u ratio | vs. *matched* approach flow |
|---|---|---|
| 90° | 0.91 | 0.91 (comparison is roughness-clean) |
| 210° | **1.51** | **0.50** |
| 270° | **1.57** | **0.52** |

**The claim "buildings raised σ_u by half again in two of three directions" comes
exclusively from the two confounded directions.** Against matched approach flow the
buildings *reduce* σ_u in all three. The headline as written is dead.

A second, independent defect: the control has no disc-band mesh refinement. All cases
use `first_dz = 4.0 m` background; production cases add a `fieldBox` at refinement
level 3 (→ 2.0 × 0.50 m cells) over the pitch, and the empty case appears not to. The
control's profile (U = 0.834, 1.587, 2.256, 2.810 m/s at z = 0.5, 1, 1.5, 2 m) is
**linear from the no-slip wall to the first cell centre**, not a log law.
Richards–Hoxey with the fitted `u*, z0` gives U(2 m) = 4.073 m/s, so **the control
under-reads mean speed by 31% at 2 m** and recovers to <1% error only at z ≈ 6 m.

Consequences, and one of them is in the author's favour:
- Mean-wind sheltering is **0.30–0.35×, not 0.44–0.50×** — i.e. *stronger* than reported.
- Control TI at 2 m is ~0.083, not 0.120, so the 90° TI amplification is ×3.0, not ×2.05.
- σ_u in the control is essentially unaffected (k is height-constant in Richards–Hoxey),
  which is why these two defects are separable and act on different rows of the table.

**What survives C2 and should become the paper's spine:** the site is deeply and
uniformly sheltered in the mean; Campus Towers (44.8 m, footprint starting 9 m off the
eastern sideline) and not the lakefront IBL controls disc-altitude flow; and shelter
removes mean wind *more efficiently than it removes fluctuation*, so TI rises. That last
one is real, counterintuitive, and defensible — it is just weaker than "buildings add
turbulence," and on present evidence it is demonstrated cleanly in **one** direction.

### C3. Q is not a usable estimator, and the ranking it produces is not resolvable.

Read `P5_Synthesis/throwability.py`: `Q = clip(1 − |U|/U_ref) · clip(1 − √(2k/3)/|U|) · mean(success_rate)`.

Failures, each independently confirmed:

1. **Non-monotone in wind speed, with indeterminate sign.** At fixed k,
   `∂ln w_f/∂ln|U| = −0.37` but `∂ln t_f/∂ln|U| = +0.65`, because `I = σ_u/|U|` *falls*
   as U rises. Net from the two penalty factors is **+0.28: Q increases with mean wind.**
2. **The calm-spot pathology.** As `|U| → 0` at finite k, `I → ∞`, `t_f → 0`, `Q → 0`.
   The field contains such cells (`\UfieldMin` = 0.02, 0.07, 0.12 m/s). **Q scores its
   calmest locations as unplayable** — a sign error in intent for a metric whose purpose
   is to find sheltered throwing lanes. The code already prints a saturation warning; it
   fires because the physics of the metric is wrong, not because the site is unusual.
3. **Top of scale wrong by ~7×.** The docstring says "Q = 1 is calm-air-equivalent."
   In calm air `w_f = t_f = 1` but `S̄ = P(hit) ≈ 0.09–0.15`, so **calm air scores
   Q ≈ 0.13.** A reader sees Q = 0.048 and concludes 5% playable; the correct reading is
   ~35% of calm-air-equivalent.
4. **Double counting.** `S̄` is already a functional of U and k (the trajectory is
   advected by U; the gust amplitude is exactly `√(2k/3) = I·|U|`). Multiplying by
   `w_f(U)·t_f(U,k)` counts the same two fields again with an arbitrary functional form.
   **Any monotone reweighting changes the ranking**, and no principle fixes the weight.
5. **`σ_u/U` is the wrong normaliser for a projectile.** A disc has never heard of `U`;
   it responds to `σ_u/V_d ≈ 0.03–0.06`, which never saturates. This is not just a
   metric bug — it is the paper's most valuable idea, buried in a broken estimator
   (see C5).

**And the ranking is not resolvable.** The between-heading spread in field-mean 2 m
speed is 0.21 m/s (1.21 / 1.42 / 1.40). The fine-vs-coarse mesh difference in the *same
quantity for a single heading* is 1.297 → 1.021 = **0.28 m/s**. Discretisation error for
one heading exceeds the entire physical signal across all three, before any other error
source is counted. Smallest meaningful ΔQ ≈ 0.02 (charitable) to 0.04 (Roache
F_s = 3); largest observed pairwise ΔQ is **0.0127**. No heading comparison in this
study is resolvable, and with `z0_inlet` differing by direction it is *confounded* as
well as unresolved — an unresolved comparison can be published with error bars; a
confounded one cannot be published at all.

### C4. The flight-dynamics claims that must be withdrawn.

The 6-DOF core came through the audit well: the five prior Hummel corrections are
correctly implemented, the frame reconciliation is right, Euler's equations check
against an analytic gyroscope to 0.09%, the advance-ratio definition matches the tunnel
convention, and Galilean invariance is exact. The inertia ratio even passes an untested
free check — `I_s/I_p = 1.926` predicts a nutation period of 0.0796 s at ω = 41 rad/s
against Lorenz's measured ~0.09 s wobble, 12% agreement from a purely inertial quantity.

But four downstream claims do not survive:

**(a) "No side force ⇒ crosswind drift is a lower bound" is false, and the stated
mechanism is backwards.** The paper says the model has "no yaw stiffness." An
axisymmetric spinning disc has no yaw stiffness *by symmetry* — Potts & Masters (2015)
and Sung et al. (2018) both state this explicitly — so adding `C_S` would not create
any. The only admissible side force is Robins–Magnus, `F ∝ ω_z b3 × V_par`, whose
**direction flips with handedness**. Measured with `C_Yr = 0.02–0.05`: the correction is
+0.744 m for a right backhand in a +5 m/s crosswind and −0.722 m in a −5 m/s crosswind,
and swaps sign between backhand and forehand. **For half the (handedness × crosswind)
combinations the missing term would *reduce* drift.** Not a bound in any sense. Its
dominant effect is on the *still-air lateral landing point* (−3.7 to −8.9 m over a 29 m
throw), i.e. it corrupts **aim**, which is correctable, not **scatter**, which is not.
Hummel (2003) reports her own best-coefficient simulation of flight f2302 curving 1.7 m
*left* where the measured flight curved 1 m *right*, and attributes the residual to
exactly this omission. This model inherits that 2.7 m sign error.

**(b) The gyroscopic-immunity hypothesis has no threshold.** The precession relation
gives a *persistent* bank change per gust `δφ = δM·τ_g/(I_s ω)` — persistent because a
gyroscope has no restoring term — so the attitude-mediated deflection goes as **1/ω,
monotonically, with no critical rate**. But there is a second channel spin cannot touch:
the gust rotates the relative wind by `β = u'/V`, giving `δa ≈ g u'/V`, which is
**spin-independent**. Their ratio is `t_r q̄Ad C_Mα / (2 I_s ω)` — the gust amplitude and
duration cancel — so the crossover is a pure property of the disc and the throw:

```
ω* = C_Mα · ½ρV²Ad · t_r / (2 I_s)
```

≈ 160 rad/s analytically, 67–97 rad/s measured. Human release spins are 46.5 rad/s
(Hummel, 57% effort), 59–81 (Sasakawa, unskilled→skilled forehand), ~94 (Winograd,
elite backhand). **Throwers already straddle ω*; there is no untapped high-spin
regime.** Measured full Monte Carlo: total landing scatter goes 2.30 → 3.00 m as spin
goes 50 → 400 rad/s. **Scatter increases with spin.** The narrow grain of truth in the
coaching folklore: *cross-range* scatter does fall (2.01 → 1.30 m, ≈ ω^−0.21) — paid for
in down-range scatter, because a higher-spin disc fades less and integrates more wind.

**(c) The hammer is a model failure, not a short-range throw.** Instrumented as
parameterised in `PAPER.tex:511–513` (V₀ ~ N(11, 1.5), α_rel ~ N(+40°, 8°)):
**100.0% of integration steps lie outside the coefficient table's valid range**
(`alpha_valid_deg = (−10, 30)`), with α spanning −78°…−40°. Every coefficient over the
entire flight is pinned at the α = −10° edge: `C_L = −0.225, C_D = 0.15, C_M = −0.038`.
A real inverted disc at |α| = 40° has |C_L| ≈ 1.2–1.7 and C_D ≈ 0.7–0.9. The model gives
it about **one sixth** of the lift and drag it should have. Separately, the release speed
is roughly half a real one: every measured Ultimate release is 20–22 m/s, and at
V₀ = 11 m/s the *vacuum ballistic* range from 1.8 m at 20–25° is 11.4 m against the
model's 11.9 m — **the reported "still-air range of 12.1 m" is `v²/g` and contains
essentially no disc aerodynamics.** Author decision #4 in `FRISBEE_STATUS.md`
("re-score at ~12 m, its own characteristic range") **should be rejected**: it would
produce a confident-looking number from clamped coefficients and an unmeasured release
speed. Note also that no published measurement of Ultimate hammer release speed, spin,
or distance exists anywhere — 25–40 m is folklore, not data.

**(d) Release wobble is identically zero for all 97,200 trajectories.**
`disc_6dof.py:466` sets `omega_B0 = [0, 0, ±spin]`. Real throws always wobble; Hummel's
own published IC for flight f2302 includes a **measured φ̇ = −14.94 rad/s**. Measured
sensitivity: a 20 rad/s transverse rate moves the landing point **3.5 m downrange and
1.1 m across** — comparable to the *entire* turbulence-induced scatter. This directly
undermines the H2b verdict: the model reports forehand beating backhand at every
heading, but gives both throw types identical spin (N(100,10)) and zero wobble, so it
cannot represent the one measured difference between them (Winograd & Engsberg 2012:
backhands carry 35% less precession — and found **no accuracy difference between throw
types** in real players). The current result comes entirely from a 2° release-pitch
difference.

### C5. The genuinely generalizable result is in the project but not in the paper.

The transferable claim is **not** "shelter is not protection." It is the normalisation:

> For a projectile the governing turbulence measure is **σ_u/V_d**, not the turbulence
> intensity σ_u/U that underpins pedestrian-comfort criteria (Lawson, NEN 8100, COST
> 732). A sheltered site whose *turbulence intensity* triples can be **easier** to throw
> on, because TI rises through a collapsing denominator while the physically relevant
> σ_u/V_d falls.

This is sharper, survives C2, and is demonstrated by this site rather than contradicted
by it. It also comes with its own diagnostic and its own falsification.

Two corollaries the panels derived and the project should own:

**The corrected dispersion law.** The memo's proposed `σ_landing ~ σ_u T_f/√N_eddy` is
the **tracer** result (Taylor 1921) — appropriate to a fluid particle, not to a disc.
The disc's ballistic length `ℓ_b = 2m/(ρ C_D S) = 51 m` gives a Stokes number
`St = ℓ_b/L ≈ 8.5`: **high Stokes, the disc does not follow eddies.** It responds
through *force*, so displacement is a **double** time-integral, not a single one. The
dominant channel is lift modulation about trim (`L ≈ mg`), so gravity and the glide
ratio enter and two extra powers of `V_d` appear. Deriving the response about trim and
integrating along a straight path through a frozen exponentially-correlated field:

```
σ_landing = (κ · g · σ_u · T_f²) / (2 V_d) · sqrt(Φ(N)),        N ≡ L_f / L

Φ(N) = 8 [ 1/(3N) − 1/(2N²) − e^(−N)/N³ + (1 − e^(−N))/N⁴ ]

κ = sqrt( 4 + (C_Lα/C_L)² (σ_w/σ_u)² )     ≈ 4.96 isotropic, 3.03 at σ_w/σ_u = 0.5
```

Limits: `Φ → 1` as `N → 0` (coherent gust — **saturates**, the proposed form diverges);
`Φ → 8/(3N)` as `N → ∞`, recovering `1/√N` with prefactor **1.633**. The asymptotes
cross at `N* = 8/3`. At the site's `N ≈ 4` the bare `1/√N` form is **22% wrong** — this
is the crossover, not an asymptote.

Ratio of correct to proposed amplitude is `(κ/2)·Fr⁻¹` where `Fr⁻¹ = gL_f/V_d² = 2.45`
here, i.e. **~6–8×**, and it is not a constant — it is a variable of the problem. At
fixed `L_f` the correct law goes as `V_d⁻³` where the proposed goes as `V_d⁻¹`; at
fixed `V_d`, as `L_f^1.7` vs `L_f^0.5`. **These are trivially distinguishable in the
existing 6-DOF model**, which is why WP2 is cheap and decisive.

Corrected numbers give σ_L = 3.4 / 5.5 / 6.1 m for 90/210/270 — the right order to
explain the observed MC scatter (S̄ ≈ 0.09 at r = 3 m implies σ_land ≈ 6.9 m), where the
proposed law gives 0.4–0.8 m. **The proposed law under-predicts by an order of
magnitude and would have looked wrong the instant it was compared to the pipeline's own
output. It never was.**

**The frozenness criterion should be `V_d/σ_u`, not `T_L/T_f`.** Frozenness asks whether
an eddy evolves while being *traversed*. Scale by scale, `τ_turnover/τ_cross = V_d/u(ℓ)`,
worst at the largest scale, so a single number bounds it: `V_d/σ_u = 17–36` here — frozen
to 3–6%, **better than a hot-wire in a typical wind tunnel.** The memo has been
*under-claiming* by using a criterion with a margin of 4–7 instead of one with a margin
of 20–36. This also matters defensively: `T_L = k/ε` is the Eulerian turnover, and a
referee substituting the standard Lagrangian value `T_L^Lag = (4/3C₀)k/ε ≈ 2.3–6 s` gets
`T_L/T_f = 0.9–2.4` and the "frozen" claim looks marginal. `V_d/σ_u` survives that
substitution.

### C6. Priority risk that must be closed before writing anything.

None of the CFD/theory panels flagged this; the strategy panel did, and it is the single
largest threat to the novelty claim:

- **Aircraft gust-load theory already owns the sweeping correction.** CS-25 Appendix G
  and the von Kármán continuous-turbulence model map spatial to temporal frequency via
  true airspeed. Etkin (1981) and Hoblit, *Gust Loads on Aircraft*, are exactly this
  problem. A referee from aeroelasticity will say "this is Etkin with a frisbee."
- **Curcin & Harris (1986)** and the windborne-debris literature (Kordi & Kopp 2012) may
  already contain the Monte-Carlo-through-synthetic-turbulence result.
- **CAP 437 / CAA Paper 2008/03** already define an operational turbulence criterion for
  a vehicle in the built environment (helideck σ_w threshold). So "nobody has written
  down a projectile criterion" as it appears in `FRISBEE_GENERALIZATION.md` §1 is **too
  strong and will be challenged.**

The reframe that survives all three, and it is stronger than the original: **CAP 437
specifies an *absolute* σ_w threshold, not a TI threshold — precisely because a
helicopter supplies its own velocity scale. The one existing operational standard
already got the denominator right, empirically, and nobody wrote down why.** You are not
competing with CAP 437; you are deriving it and generalising it to a class of objects it
was never written for. What remains unclaimed is (i) the **closed-form landing-variance
law for a *lifting* body** — aircraft work computes loads and response, not landing-point
variance; ballistics computes landing variance for non-lifting bodies; the lifting
projectile falls between them — and (ii) the demonstration that a pedestrian criterion
and a projectile criterion **rank the same site oppositely**.

---

## 2. Where the panels disagreed, and the adjudication

| Question | Positions | Adjudication |
|---|---|---|
| What first? | CFD panel: matched control. Theory panel: correlation time + validate the law. Flight-dynamics: fix the gust model. UQ: free artifact checks. | **The gust-correlation fix and the "frozenness criterion" fix are the same bug.** Do the free artifact checks (WP0, 30 min), then the gust fix (WP1), then validate the law (WP2). The CFD control (WP3) is cheap and important but it only affects the site paper, which is being demoted anyway. |
| LES of the easterly? | Status memo ranks it #2; review dossier recommends it; CFD panel rejects as *next*; theory panel reframes it; strategy panel kills it. | **Not now.** LES against a mismatched control reproduces the same wrong ratio at 100× the cost, and it needs a matched-z₀ precursor anyway — so WP3 is a prerequisite, not an alternative. If it is ever run, the deliverable is `L_11`, `L_22` and `σ_w/σ_u` in the disc band — **not** "gust structure." Those two numbers move `√Φ` by ×3 and `κ` by ×1.6; peak gusts move nothing. |
| The idealised H/D shelter sweep (memo's #1)? | Memo: highest value, 12 cases. CFD panel: right instrument, wrong parameterisation. Theory: premature. Strategy: gate it. | **Gated to WP7, and restructured.** `H/D` is not a case parameter — `x/H` is a *sampling coordinate*, so one run of a block in a 20H domain gives the whole curve. The parameter that decides whether an optimum exists is **porosity**, then approach-flow roughness. Do not run it until the corrected law predicts a non-monotone response. |
| Third mesh? | Status memo ranks it #1. CFD panel: right idea, wrong mesh. UQ: needed for a real GCI. Strategy: kill. | **Reformulated and deferred.** Do not add another `fieldBox` level — you have already shown that refining where the answer is *sampled* doesn't move it. The turbulence is *produced* in the separating shear layers off Campus Towers and Mertz, whose surfaces sit at 4–8 m cells (3–6 cells across a 23.5 m building, against an AIJ guideline of ≥10). Refine the **source**, not the sample. And with two meshes there is no observed order `p` and therefore no GCI — the "57.1%" is a max-difference, not an error estimate. |
| 13-direction sweep, second site, second reference speed, building-height sensitivity? | Status memo §8 items 3–6. | **All rejected or redirected.** 13 directions multiplies every defect by 13. A second site yields n = 2 — the idealised sweep *is* the generalisation. A second reference speed needs **no CFD**: the flow is Re-independent so the field rescales linearly; scale the existing field and re-run the MC. Height sensitivity perturbs mostly far-field type-defaults, while the two buildings that set the answer both carry real `building:levels` tags. |
| Is the site study still publishable? | CFD panel: yes, with corrections. Strategy: demote to a demonstration. | **Demote.** With C2 applied the surviving site findings are "a 45 m building 9 m off the sideline dominates the flow at 2 m" and "shelter reduces mean speed to 0.3×" — which is the *null hypothesis* of urban wind engineering, not a result. It earns its place as the demonstration section of a theory paper, because it supplies a real σ_u and L at a real site and is the case where TI and σ_u/V_d rank the site oppositely. |

---

## 3. Work packages

### WP0 — Free artifact checks — **ALREADY RUN, 2026-07-26. Results below.**

Five checks against artifacts already on disk. WP0.1–WP0.3 were executed against
`landing_stats.npz` during the review; **read these results before starting WP1**, because
one of them found a defect no panel had identified and it changes what WP0.5 and WP5 must
do. WP0.4 and WP0.5 remain to be done.

---

#### WP0.1 RESULT — blocking is a real differential, but it does **not** explain the ranking

```
θ      blocked (excl hammer)    S̄ (excl hammer)
 90          15.463%                8.287%
210          14.389%               10.806%
270          15.704%                9.241%
spread        1.315 pp              2.519 pp
```

The spread is **1.315 pp, above the 1.0 pp threshold** — so per-heading blocking is a real
differential effect and must be reported rather than left in `landing_stats.npz` unread.
But correcting for it does **not** collapse the ranking:

```
S_cond = S / (1 - blocked):   90 -> 9.803%,  210 -> 12.622%,  270 -> 10.962%
raw spread 2.519 pp  ->  blocking-corrected 2.819 pp   (it slightly WIDENS)
```

**The UQ panel's leading hypothesis is therefore partially falsified.** Report the
blocked-fraction differential, but do not attribute the ranking to it. The ordering is
also not monotone in blocking (90° has less blocking than 270° but a lower S̄).

#### WP0.2 RESULT — cleanly closed, no action needed

`n_valid / n_throws = 1.0000` for every cell in all three headings; zero NaN
`success_rate` cells. **Nothing is being silently dropped**, and the "disc scored as
landed at `t_max` altitude" failure mode is not occurring in this dataset. Still add the
`landings[:,2] ≈ 0` assertion as a regression guard, but the concern is not live.

#### WP0.3 RESULT — **NEW DEFECT, found here, not on any panel's list**

**`landing_mean` / `landing_cov` / `ellipse_axes` are not stored in the same reference
frame as the target that `success_rate` is measured against.**

Evidence — launch cell [0,0], θ = 210°, right backhand:

```
bearing   stored |Δ| from "target"    success rate
    0°          1.24 m                   11.3%
   60°         26.21 m                   11.3%
  120°         45.71 m                   10.0%
  180°         51.61 m                   10.7%
  240°         44.66 m                   14.0%
  300°         26.34 m                    6.0%
```

A landing centroid 51.6 m from a 3 m target with σ = 6.0 m has a hit probability of
essentially **zero**, yet 10.7% is reported. The stored offsets are almost exactly the
distance between a **fixed** target at `(R, 0)` and the **bearing-rotated** target at
`(R cos b, R sin b)` — predicted 0 / 25.0 / 43.3 / 50.0 / 43.3 / 25.0 m against observed
1.24 / 26.21 / 45.71 / 51.61 / 44.66 / 26.34 m. Correcting for the rotation collapses
them to a true bias of **1.2 / 2.7 / 7.3 / 8.3 / 4.3 / 4.7 m**, which is physically
sensible and consistent with the reported success rates.

Note this is precisely the failure mode the earlier `success_rate` fix was meant to cure:
that fix corrected `success_rate` (which now uses the rotated target) but **left
`landing_mean` in the old frame.** Averaging over launch cells suggests a second,
compounding offset — the residual does not collapse the same way across all 12 cells —
so **determine the actual convention from `monte_carlo.py` before trusting any
correction.** Do not simply apply the rotation fix above; diagnose it.

**Consequences:**
- Any downstream use of `landing_mean` or `ellipse_axes` is wrong for 5 of 6 bearings.
  Check whether anything in P5 or the figures touches them.
- **`landing_cov` is safe.** A covariance is translation-invariant, so σ_land is
  unaffected by the frame error. This matters: WP5's proposed replacement metric
  `F = σ_land(wind)/σ_land(still)` is built on the covariance and is therefore **not
  contaminated** by this defect. It is another reason to prefer F over anything
  bias-based.
- Add this to the WP0.5 / provenance list: the schema must state, per field, which target
  frame it uses.

#### WP0.3 RESULT (cont.) — the scatter, and what does *not* explain the ranking

From the covariance (frame-safe):

```
θ      σ_land        blocked      S̄
 90     5.46 m       15.46%     8.29%
210     5.58 m       14.39%    10.81%
270     5.94 m       15.70%     9.24%
```

σ_land ≈ 5.5–5.9 m at a 25 m range — confirming the UQ panel's independent estimate of
~6.9 m to the right order, and confirming that this is **not a model of a club Ultimate
player**, who completes 25 m backhands well above 90% (risk-register item 4).

Critically, **the rankings do not line up**:

```
rank by S̄ (best first) : 210, 270, 90
rank by low σ_land      :  90, 210, 270
rank by low blocking    : 210,  90, 270
rank by low bias        : 270, 210,  90
```

No single mechanism orders the headings the way S̄ does. Combined with the mean-flow
argument (a 0.21 m/s speed spread shifts the centroid ~0.27 m, changing hit probability
by ≲1%, against an observed 27% relative S̄ spread), **the heading ranking is not
currently explicable by any identified physical or numerical mechanism.** That is a
stronger reason to withdraw it (WP4.3) than "it is unresolved" — and the frame defect
above must be diagnosed before the question can even be posed properly.

---

Remaining WP0 items:

**WP0.1 — Is the heading ranking a meshing artefact?** *(RUN — see result above; answer
is "partly, but it does not explain the ranking." The reasoning below is retained because
the `solid`-mask comparison was not run and is still worth doing.)*

There is a diagnostic that is cheaper than anything else in this document. The mean-flow
displacement is ~1.3 m per (m/s) of wind (linearised: τ = m/(ρA C_D V) ≈ 1.5 s, so
displacement = w[T − τ(1−e^(−T/τ))] ≈ 1.3w for T = 2.5 s). The 0.21 m/s between-heading
speed spread therefore shifts the landing centroid by **0.27 m**. Against σ_land ≈ 6.9 m
that changes hit probability by **≲1%**. **The observed S̄ spread is 30% relative
(8.3% → 10.8%).** A 1% mechanism cannot produce a 30% effect. Something else is
producing the ranking.

Leading candidate: `sample_disc_altitude.py` builds the `solid` mask independently per
heading from a different snappyHexMesh result via `(~valid).any(axis=2)`, and
`CFDField.is_solid` then dilates it by one cell in each direction.

```
# from landing_stats.npz — blocked_fraction is already saved per cell
print blocked_fraction, averaged by heading
# from flow_fields.npz
np.array_equal on dir_090__solid, dir_210__solid, dir_270__solid
```

**Acceptance / decision:** if per-heading blocked fractions differ by more than ~1 pp,
the S̄ ranking is substantially a per-mesh ray-casting artefact and the WP4 withdrawal
of the ranking becomes **mandatory** rather than precautionary. Report the numbers
either way.

**WP0.2 — `n_valid / n_throws` by heading.** *(RUN — clean, see result above. Add the
assertion as a regression guard anyway.)* `success_rate` removes non-landing throws
from the denominator while keeping blocked throws in it as failures — two failure modes,
two different treatments, neither documented. Worse, `valid = np.isfinite(landings[:,0])`
tests only finiteness: if `simulate` returns the state at `t_max = 8 s` without a ground
event, the disc is scored as having "landed" at its current altitude. **Add an assertion
that `landings[:,2] ≈ 0`** and report `n_valid/n_throws` per heading.

**WP0.3 — Landing scatter and success rate use different populations.** `mean`, `cov`,
`ellipse_axes` are computed from `landings[valid]`, which *includes* blocked throws;
`success_rate` excludes them. Make them consistent and note which changed. *(RUN — and it
turned out to be worse than a population mismatch: they use different target **frames**.
See the WP0.3 result above. Diagnose the convention in `monte_carlo.py` before fixing.)*

**WP0.4 — Record what the gust model actually did.** Add realised per-cell `T_L`,
whether `T_L_MAX` clipped, and `N_MODES`, `GUST_DT`, `T_L_MIN`, `t_max`, `MC_RTOL`,
`MC_ATOL`, the git SHA and library versions to `landing_stats.npz`. None is currently
recorded.

**WP0.5 — Fix the results table that contradicts itself.** `tab:Qbyheading` prints
`S̄` = 5.5 / 7.2 / 6.2%, the mean over **all three** throws including the hammer's
structural zeros — but Q is computed with the hammer **excluded**, from 8.3 / 10.8 /
9.2%. Check: 0.727 × 0.744 × 0.083 = 0.0449 ≈ 0.0443 ✓ versus 0.727 × 0.744 × 0.055 =
0.0298 ✗. **A reader who multiplies the printed row gets Q wrong by a factor of 1.5 and
concludes the pipeline is broken.** Emit the hammer-excluded `S̄`, or add a `\SRexcluded`
macro and a footnote.

---

### WP1 — Fix the gust model (≈3–4 h code + ~1 h compute) — GATE G1

**Actions, in `P4_Coupling/synthetic_turbulence.py` and `monte_carlo.py`:**

1. **Correlate along path arclength, not wall-clock time.** Replace the temporal
   correlation `T_L = k/ε` with a sweeping decorrelation: make `ξ` a function of the
   trajectory arclength `s` with decorrelation length `L_int` (equivalently
   `τ_c = L_int/V_d ≈ 0.4 s`). Keep the Shinozuka–Deodatis spectral representation —
   it is analytic in its argument and costs the integrator nothing, which is why it was
   adopted in the first place. Do **not** revert toward white noise; the original
   white-noise generator made the ODE non-deterministic and is correctly gone.
2. **Update `σ(z)` and `L(z)` along the path**, not fixed at the launch point
   (`monte_carlo.py:85`, `ref_pos=launch`). Until this is done, the "throw low wins
   twice" argument is not represented in the simulation at all.
3. **Add anisotropy:** `σ_v/σ_u ≈ 0.75`, `σ_w/σ_u ≈ 0.5`. Expose it as a parameter, since
   it enters `κ` and the panels flag it as one of the two dominant amplitude
   uncertainties.
4. **Record the realised `τ_c`, `L_int` and any clipping** per cell (WP0.4).
5. Fix the stale module docstring claiming `T_L` is "O(1 s)".
6. `max_step = 0.05` in `disc_6dof.py:319` will start to bind: it is 1.5× the nutation
   period at ω = 100 (0.0326 s) and only ~7 steps per correlation time at τ_c = 0.4 s.
   Check whether it binds after the change and reduce if so.

**Then re-run a 5k-throw subset** (~3 min at 27 throws/s) and compare to the old numbers.

**GATE G1 — acceptance:** turbulence-induced landing dispersion must **fall by a factor
of 2.0–2.6** (≈ `√N_eddy` for N ≈ 4–6). Two independent panels predict 2.3–2.6× from
different arguments.

- **If it falls in range:** proceed to WP2. The two panels agree and the physics is
  understood.
- **If the change is < 20%:** the defect was harmless, the existing Q numbers stand, and
  something in the panels' reasoning is wrong. **Stop and report** — do not proceed to
  WP2 on a theory whose central prediction just failed.
- **If it falls by ≫ 3×:** the disc is not behaving as a linear filter. Stop and report.

---

### WP2 — Validate the corrected dispersion law (≈1–2 days, zero CFD) — GATE G2

This is the actual scientific content and it costs no compute. Run **in homogeneous
synthetic turbulence** with `σ_u`, `L`, `V_d`, `T_f` swept independently — not in the
site field, where they co-vary and nothing is identifiable.

Use paired common-random-number arms throughout (still air / mean only / mean + gusts),
which the code already supports: `CorrelatedGusts.__init__` draws exactly
`3 × N_MODES = 96` uniforms regardless of `T_L` or field, so the arms genuinely share
release draws. **Add a test asserting identical RNG state advance across arms** — this
invariant is currently undefended and will break silently if `N_MODES` ever becomes
adaptive.

**Also required before the sweeps: drop the heading index from the seed key.**
`monte_carlo.py:247` uses `spawn_key = (ix, iy, i_h, it, ib)`, so **heading 90 and
heading 210 draw independent release parameters** — the single comparison the study
exists to make is the one comparison that gets no variance reduction. Change to
`(ix, iy, it, ib)` and pass heading as a non-random parameter. One line; an estimated
5–20× reduction in Var(ΔS̄), since ~98% of per-throw variance is the shared release draw.
The same applies to the forehand-vs-backhand contrast (`it` is also in the key).

**Predictions to test.** Each has a stated falsification.

| # | Test | Predicted | Falsified if |
|---|---|---|---|
| **P1** | The WP1 fix itself | σ_L ratio 0.67–0.70 | outside 0.55–0.85 |
| **P2** | **Range sweep** `L_f ∈ {12,18,25,32,40}` m at fixed release speed — *the exponent test, highest discriminating power* | `σ_L ∝ L_f^1.7` (local exponent 1.80 at N=2 → 1.54 at N=20) | fitted exponent < 1.2 or > 2.1. The proposed law predicts **0.5** — the two hypotheses differ by >3× across this range, so no realistic MC noise can hide it |
| **P3** | **Component decomposition:** gust in `u'` only / `w'` only / all three | `w'`-only carries **84%** of gust-induced variance, `u'`-only 16% | `u'`-only ≥ `w'`-only. This is the sharpest single test of the lift-modulation mechanism and costs one extra flag |
| **P4** | **Φ(N) measured directly:** `τ_c` is a free parameter — sweep it over two decades (0.05–5 s) at fixed σ_u, ~200 throws each | the closed-form `Φ`; the `−1/2` slope appears only for `N ≳ 20` | >15% deviation, or the `−1/2` slope reached at N < 8. Settles the `N_eddy` definitional ambiguity empirically |
| **P5** | **Spin sweep** `ω_z ∈ [25, 400]` rad/s, paired, decomposed into force-only / moment-only channels | moment channel `≈ ω^−1`, force channel `≈ ω^0`, crossover at ~10² rad/s, **total varies < 15%** | total `∝ ω^−1`, or variation > 25%, or the exponent *steepens* with ω (a genuine threshold) |
| **P6** | **Linearity:** scale `k` by {0.25, 1, 4} | `d ln σ_L / d ln σ_u = 1` exactly | outside 1 ± 0.05 — which would mean α excursions are leaving `alpha_valid_deg = (−10, 25)` and the linear-response framework is invalid at this turbulence level |
| **P7** | **Release-scatter sweep** `s ∈ {0.5, 1.0, 1.5}` on the release covariance — *the identifiability test* | wind signal `(S₁−S₀)/S₀ ∝ 1/s²` | if confirmed, **Q is not a site property** and must be reported at a stated thrower skill or not at all |

**P5 has a cheap deterministic form worth running first:** a paired impulse-response test
— a single Gaussian gust pulse (1 m/s, τ = 0.15 s, at t = 0.6 s), 5 spins × 3 gust
directions × 3 channel variants = **45 throws, ~2 minutes**, no Monte Carlo noise. It
isolates the mechanism rather than the outcome and returns `ω*` directly.

**Note the definitional trap in `N_eddy`.** For *lateral* deflection driven by the
*transverse* gust component sampled along a line, isotropy gives `L_22 = L_11/2`, so N
doubles. Disc-band L (6.2 m) vs 2 m L (4.4 m) vs `L_22` gives N = 4.0, 5.6, 11.3 and
`√Φ` = 0.67, 0.60, 0.46 — a **45% spread from definitional choices alone.** Separately,
the code's `L_int = σ_u k/ε = √(2/3) k^1.5/ε` is **exactly 4.97× larger** than the
standard RANS mixing length `C_μ^0.75 k^1.5/ε`, always. Pick one convention, define it
in the paper, report `Φ` rather than `N^−1/2`, and **report the law's sensitivity across
all three conventions.** If the regime claim is not robust across them — and it is not —
the honest statement is "N ∈ [1, 22]; here is the measurement that would pin it."

**GATE G2 — the project gate:** does the closed-form law reproduce the measured σ_L
within ~20% across the sweeps?

- **Yes:** you have the paper. Proceed.
- **No:** the theory is wrong and there is no theory paper. **Stop and reassess the
  entire program** before spending anything on WP3 or later.

---

### WP3 — Matched-roughness empty-fetch control (≈1.5 h compute) — GATE G3

**Read §1/C2 first — the specifics differ from what one panel reported, and I verified
the case metadata directly.**

`generateCases.py` already supports what is needed: `--empty`, `--z0-ground`, `--tag`,
and `z0_inlet` is derived from heading via `inlet_z0(theta)` with
`Z0_WATER = 0.001` / `Z0_URBAN = 0.5`. So:

1. **The urban control does not exist and is the one that matters.** Generate an empty
   case at an urban heading so it picks up `z0_inlet = 0.5`:
   `python3 generateCases.py --headings 210 --empty --resolution production --tag _empty`
   This is the matched counterfactual for **both** 210° and 270°, which is exactly where
   the dead claim came from.
2. **Check whether `dir_090_empty` got the `fieldBox` refinement.** Its profile is linear
   from the wall to 2 m, which says it did not. If so, regenerate it with the production
   refinement present, or the 2 m denominator stays a wall interpolation.
3. **Re-derive the `counterfactual_2m` table** in `vertical_structure.py` against the
   **matching** control per direction (not one control for all three), and **band-average
   over z ∈ [1, 6] m** rather than reading a single height.
4. This run simultaneously closes a validation gap the paper already flags: the existing
   `\ABLverdict{PASS}` was measured on a smooth-wall case (`z0_inlet = z0_ground = 0.001`)
   over the band `z ∈ [6, 100] m` — which **excludes the entire disc band** `(1, 6) m`
   where every reported number lives. The one variant with realistic ground roughness
   (`cfd_validation_abl.json`) is **MARGINAL with 107% k drift**. Evaluate homogeneity
   inside the disc band this time.

**GATE G3 — acceptance:**
- `|U(z) − Richards–Hoxey| < 5%` and `|k − u*²/√C_μ| < 10%` over `z ∈ [1,6] m` at the
  site's 304 m fetch, **for each control**.
- Then recompute the σ_u and TI ratios and **publish them whatever they say.**
- **Predicted outcome:** σ_u ratios ≈ 0.50–0.52 for 210°/270° and < 0.91 for 90° — i.e.
  buildings *reduce* σ_u in all three, and "shelter is not protection" as written must be
  formally retracted from the paper, the memo, and any notes.

**Also extend `validate_cfd.py::test_grid_independence` to report `k`, `σ_u`, `TI`,
`T_L` and `L`** at the same five stations from the two meshes you already have. It
currently compares **velocity only**, and the entire generalization rests on `k` — the
load-bearing quantity has no discretisation error estimate whatsoever. Zero compute,
~1 h scripting. **Acceptance: if `Δσ_u/σ_u > 30%`, no quantitative σ_u ratio may be
published — only signs.**

---

### WP4 — Retractions, corrections and manuscript hygiene (zero compute, ~1 day)

**Nothing gets posted until this is done.** Flag each item in the commit message per the
house rule; make the mechanical fixes directly, leave the physics judgement calls marked
for the author.

**Statements to withdraw or restate:**

1. **`PAPER.tex` line 259, contribution #3** — verified — still reads *"A case study
   identifying the lakefront internal boundary layer as the dominant local control on
   Ultimate-disc flight at the chosen site."* **This directly contradicts the abstract
   and the conclusions in the same PDF**, and it is the premise the study falsified.
   First thing a referee circles; highest-priority single edit in the document.
2. **§1 "Feasibility check" (line 117) and §2 "Project plan" (line 159)**, including the
   14-week work breakdown and the risks-and-mitigations table. The review dossier waved
   these through with "a status box states the paper sections are authoritative." It
   will not survive. **A document that opens with its own Gantt chart is a project
   report, not a paper.** Excise before anything is posted anywhere, including arXiv.
3. **The heading ranking, `Q_min`/`Q_max`, `\QbestWorstRatio`, and 3-significant-figure
   Q.** Replace with: *"all three computed directions give 2 m mean speeds of 24–28% of
   the 10 m reference and Q within mutual uncertainty; the site is uniformly sheltered
   and, at the resolution of this study, direction-insensitive."* This is cleaner,
   defensible, and consistent with the study's own `\UfieldFracMin`/`\UfieldFracMax`.
   Add the caveat that the fine mesh gives **higher** speeds than the coarse (1.297 vs
   1.021 m/s), so if the trend continues the "24–28%" is **biased low** and errs toward
   *overstating* shelter.
4. **"Shelter is not protection"** — retract as written (C2). The corrected, sharper
   version to replace it with: *"shelter is not protection **when it raises σ_u** — and
   TI is a misleading diagnostic for a projectile, because dividing by a collapsing mean
   wind manufactures a factor of 2–3 where the physical penalty can be zero."* Note the
   concrete example that makes it vivid: at 90°, TI goes 0.12 → 0.25 (the "factor of 2–3"
   headline) but the corrected law gives σ_L = 3.40 m against the control's 3.34 m — a
   **2% change**. The TI rise there is almost entirely the denominator.
5. **"No side force ⇒ crosswind drift is a lower bound"** — appears **twice**, at
   `PAPER.tex:339` (footnote) and again at `PAPER.tex:1039` (limitations), plus
   `P1_AnalyticalModel/README.md:111–116`. Fix all three. Replace with, in substance: *"The model omits
   the Robins–Magnus side force, measured by Potts & Crowther (2002) at C_S ≈ 0.02
   (AdvR < 0.5) rising to 0.04–0.05 (AdvR = 0.69), over an operating range of
   AdvR ≈ 0.5–0.9. Its leading-order effect is a handedness-dependent lateral bias of
   order 1–4 m over a 25 m throw; its effect on sensitivity to crosswind is of either
   sign depending on handedness and crosswind direction. Lateral predictions are
   therefore **uncertain in sign, not bounded below.** Hummel (2003) reports a 2.7 m
   curvature sign error in the same coefficient set and attributes it to this
   omission."*
6. **The hammer.** Reject author decision #4. Either drop the hammer entirely, or extend
   the coefficient table past stall (WP6.1) *and* re-parameterise the release at
   18–21 m/s first. State plainly that the "12.1 m still-air range" was a clamp artefact
   and essentially `v²/g`, not a model result. The H1b "NOT TESTABLE" verdict is right,
   but for the wrong reason — the correct reason is that the model has no valid
   aerodynamic data for the hammer at all.
7. **The gyroscopic-immunity hypothesis** (`FRISBEE_GENERALIZATION.md` item 3). Rewrite
   as the null result with the derivation in C4(b): no threshold, `ω*` inside the human
   range, total scatter increases with spin, cross-range scatter falls as ≈ ω^−0.21.
   Note the sport-general prediction that falls out: `ω* ∝ V²t_r/I_s`, so the knee moves
   *up* for faster throws and *down* for rim-weighted discs.
8. **The shear null result** (`FRISBEE_GENERALIZATION.md` item 7) — right answer, wrong
   argument. The stated evidence compares the sheltered pitch's 0.02–0.06 s⁻¹ against the
   control's "up to 1.5 s⁻¹", but 1.5 s⁻¹ is the control's value at z = 0.5 m, *below*
   the disc band, and the control's near-ground profile is a wall interpolation anyway
   (C2). Replace with the scaling the memo doesn't do: shear enters through the velocity
   difference **across the disc diameter**, `ΔU = (∂U/∂z)·d = 0.06 × 0.269 = 0.016 m/s`,
   i.e. `ΔU/V_d = 1.3e-3` against `σ_u/V_d = 3e-2`. **The gust-induced rolling moment is
   ~20× the shear-induced one.** Note the honest caveat: this holds in the *sheltered*
   case; in open fetch the two are comparable. And there is a second shear channel — the
   change in mean wind over the trajectory's height excursion, ~0.18 m/s against
   σ_u = 0.5 m/s — comparable, but **deterministic and therefore learnable**, exactly
   like mean wind.
9. **`Project1_AnalyticalTreatment.md` — do not post it.** It documents a *different
   model* from the one that ran: eq (26) boxes the precession rate as `M/(I_p Ω)` and
   flags it as "footnote-worthy," but the measured spin-axis tilt rate is `M/(I_s Ω)` —
   the document is **2.08× too fast**, and eq (46), which it calls "the core physics
   result," carries the same factor. It also still contains **three of the five errors
   the README says were fixed**: the `(α − α₀)` convention in eqs (11) and (14), the
   `d/(2V)` rate non-dimensionalisation in eqs (14)–(15), and a turn/fade transition of
   **6.7°** in eq (28) — described as "this single number... organizes the entire
   turn-and-fade pattern" — where the correct value is 8.8°/9.0°. Either rewrite it to
   match the shipped code or exclude it. **If it is the paper's supplement, the paper's
   stated model is not the model that produced the results.**

**Provenance and reproducibility fixes:**

10. **Add `\ConvVerdict`.** `cfd_validation.json` records `"test": "convergence"`,
    `"verdict": "FAIL"` for both meshes. `results_generated.tex` emits `\ConvFine{13.7}`
    and `\ConvCoarse{7.4}` but **no verdict macro** — while `\ABLverdict` and
    `\GridVerdict` both exist. The one FAIL that gets no macro is the one the prose then
    softens to "characterised." Emit it.
11. **Footnote or remove `\ABLverdict{PASS}`** — it was measured on a smooth-wall case
    over a band that excludes the disc band (see WP3.4). As it stands the macro name
    promises validation of the setup and the artifact validates something else.
12. **Reconcile the two contradictory grid artifacts.** `cfd_validation.json` reports
    worst 57.1%; `cfd_validation_grid.json` reports **306.96%** at station (0, −20)
    (0.353 vs 3.070 for the same station). `\GridDiff` takes the smaller. Add a run
    manifest (timestamps, case IDs, averaging window) and report both the relative and
    the `u_ref`-normalised measure — the 8.7× discrepancy also reveals that the relative
    metric is unstable because its denominator is near zero.
13. **Bind the artifacts.** `landing_stats.npz` records `field_source` as a path string
    and `compute_Q` never checks it. Re-run the CFD to the same path, re-run only P5, and
    Q silently mixes new `w_f`/`t_f` with stale `S̄`, with no error and no warning.
    Content-hash both artifacts and assert.
14. **Guard the silent still-air path.** `_get_field(None, h)` returns
    `StubField(np.zeros(3))`, so a `landing_stats.npz` produced with no CFD at all is
    indistinguishable downstream. The existing `is_stub` guard protects the fields file,
    not the landings file.
15. **`P4_Coupling/README.md` and `P3_CFD/README.md` are stale and will break a
    reproducer.** P4 advertises 10⁴ throws/cell, a `landing_stats.py` that doesn't exist,
    and an output schema the code never writes (`launch_grid[...,3]`, `headings_deg`,
    `throw_types`, `bearings_deg` vs the actual `launch_x`/`launch_y`/`launch_z`/
    `headings`/`throws`/`bearings`) — anyone following it gets a `KeyError`. Its "open
    task 1" still says the gust model is white noise and that this is "fine for
    first-pass Monte Carlo," the exact opposite of what the code now does. P3 still
    describes 13 runs, v2406, a 500×500×80 m domain, 5–8 M cells, and an AIJ-B benchmark
    "using identical mesh and BC stack" that was never run. **A referee handed the repo
    reads that last one as a claim.**
16. **Silent NaN absorption in Q.** `np.nanmean` over `success_rate` quietly drops cells
    with `n_valid < 2`. If a cell is NaN for one heading and not another, the headings are
    averaged over **different supports** and the difference is not like-for-like. Count
    and report NaN cells per heading.
17. **Check `alpha_landing`.** `alpha_lands[i] = res.diag[-1].get("alpha_deg", nan)` runs
    with `compute_diag=False` — best case all-NaN saved as if meaningful, worst case an
    `IndexError`. Verify before anyone cites it.
18. **Fix the release-convention docstring contradiction.** `disc_6dof.py:377–378` says
    `release_alpha_deg` is *"angle between disc plane and ground at release (NOT angle of
    attack — this is the body orientation)"*; `:410–413` says *"Convention:
    release_alpha_deg = angle of attack at release"* and implements the latter. Given that fix #1 was a convention bug of exactly this kind,
    a self-contradicting docstring on the release convention is how it comes back.
19. **Add a determinism test:** `--jobs 1` vs `--jobs 4` on a 2×2 grid, assert
    array-equal. The seeding design supports worker-count independence; nothing verifies
    it.
20. Also fix: `MC_RTOL = 1e-4` cites a convergence study that validated `rtol = 1e-5`;
    `DEFAULT_COEFFS` is a mutated module-level singleton so the reported "~7% clamped" is
    a nondeterministic per-process figure; the dead `TabulatedCoefficients.CR()` hard-codes
    `0.00171` inconsistently with `CR_tabulated` and is exactly the shadow path that would
    reintroduce fix #4. (`runDownstream.sh` **does** exist — one panel reported it missing,
    but that was a staging artefact of the review, not a repo defect. No action.)

**Do NOT do in WP4:** do not relax any validation criterion after seeing a result. The
project has already earned credit for not doing so once, and that credibility is worth
more than any individual number.

---

### WP5 — Replace Q (≈4 h) — after WP1

**Primary replacement — a paired dispersion-inflation factor:**

```
F(x, θ) = σ_land(wind) / σ_land(still air)        on CRN-paired throws
```

Dimensionless, ≥ 0, monotone, **no saturation, no clipping, no arbitrary `U_ref`**.
Cancels the thrower model, the disc-coefficient bias and the target-radius choice to
first order — precisely the common-mode error terms. Sampling distribution: ratio of two
paired sample dispersions; bootstrap over throw pairs. **And it makes the hammer testable
again**, since F is well-defined where S = 0 — recovering the 33% of the design that a
fixed 25 m target destroyed. It also matches the project's own theory directly: F is the
measurable form of the σ_L law from WP2.

**Secondary — a paired skill score for delivery:** `R = S_wind/S_still`, or the paired
difference with an exact McNemar CI and a Katz log-CI on the ratio. *"The wind at this
location removes (1−R)×100% of successful deliveries"* is a sentence a coach and a
referee can both check.

**If a single scalar index is wanted**, use the derived group rather than the asserted
one. The memo's `Π = σ_u T_f/(r_target √N_eddy)` fails limit checks (it diverges where
the physics saturates, is unchanged as `g → 0`, and goes as `V_d⁻¹` where the correct
law goes as `V_d⁻³`). The corrected form:

```
Π_play = σ_L / r_target = (κ/2)(σ_u/V_d)(g L_f/V_d²)(L_f/r_target) sqrt(Φ(N))

P_hit  = 1 − exp( −1 / (2 Π_play²) )              (unbiased; Rice-distributed with
                                                   non-centrality b/σ_L when mean wind
                                                   leaves a bias b)
```

Reporting `P_hit` rather than a raw ratio is what makes it the true analogue of Lawson,
which is an **exceedance-probability** criterion. Threshold bands: `Π_play < 0.47` (90%
delivery), `< 0.60` (75%), `< 0.85` (50%), `> 1.32` (unplayable, <25%).

Note the Buckingham-Pi reality check: the problem has ~15 dimensional variables and 3
dimensions, so **12 groups** — there is no unique Π, and the memo's index is one member
of that family with the other nine set silently to unity. The minimal working set is
five: `{σ_u/V_d, Fr⁻¹ = gL_f/V_d², N_eddy, κ, r/L_f}`. Say so in the paper; it is a
strength, not a weakness, to show you know which nine you fixed.

**If `(1 − I)` must be retained for continuity with the draft**, use `I/(1+I)` or
`exp(−I)` (both already recommended in the code's own docstring) and reference the wind
penalty to the disc's airspeed `|U|/V_release` rather than to an arbitrary `U_ref`.

---

### WP6 — Flight-dynamics hardening (≈6 h total, all cheap)

Ordered by value. Independent of WP1–WP5; can be interleaved.

**WP6.1 — Extend the coefficient table past stall to |α| ≤ 90° (~2 h).** Potts & Masters
(2015), *Procedia Engineering* 112:400–405, is **open access** (shura.shu.ac.uk/11092)
and carries digitised Potts-2005 curves through stall at 45° for both the Frisbee section
and a flat plate. Beyond that, blend to flat-plate theory (`C_L = 2 sinα cosα`,
`C_D = 2 sin²α`). *Expected:* the 25 m backhand barely moves — it clamps only in its final
second — while the hammer changes qualitatively. **Acceptance: if a canonical backhand's
landing point moves by < 0.5 m, clamping never mattered for the main throws.** That is a
*good* result to be able to state: it lets the paper defend the backhand/forehand
conclusions while dropping the hammer cleanly.

**WP6.2 — Add `C_S = C_Sr·AdvR` and re-run Hummel's flight f2302 (~1 h).** Calibrate so
`C_S(AdvR = 0.69) = 0.045`. f2302's full 12-state IC is published in Hummel Table 3-3:
`x,y,z = −0.90, −0.63, −0.91 m; ẋ,ẏ,ż = 13.42, −0.41, 0.001 m/s; φ,θ,ψ = −0.07, 0.21,
5.03 rad; φ̇,θ̇,ψ̇ = −14.94, −1.48, 54.25 rad/s.` The measured trajectory curved **1 m
right**; Hummel's own best simulation curved **1.7 m left**. *Acceptance:* the curvature
sign error shrinks materially. If it doesn't, the residual is in `C_Rr` or the pitching
moment — itself worth reporting, as the first quantitative test of Hummel's own
conjecture. Note this IC also carries the **measured release wobble** (φ̇ = −14.94 rad/s)
that WP6.3 is about, and the current validation uses only the §2.4.2 / §2.4.3 published
*summary* cases, not this fully-specified tracked flight.

**WP6.3 — Add release wobble to the Monte Carlo (~1 h).** Draw `ω_x, ω_y ~ N(0, σ_w)` at
release, calibrated so backhand carries ~35% less than forehand (Winograd & Engsberg
2012), anchored on Hummel's f2302. *Expected:* this may **reverse** the "forehand beats
backhand at every heading" result, which currently comes entirely from a 2° release-pitch
difference. *Acceptance:* if 10–20 rad/s of transverse rate moves the BH/FH gap by less
than the current +0.018 to +0.032, the existing verdict stands.

**WP6.4 — Re-source `C_Nr` against measured spin decay (~30 min).** Current
`τ_spin = I_s/|C_Nr| = 68.9 s`. Lorenz (2003, instrumented Frisbee) measured an
e-folding of **~14 s**; Weizman et al. (2020, onboard triaxial gyro) measured steady
aerodynamic spin-down torques of **0.030–0.078 N·m** against the model's 1.4e-3 N·m.
Over a 4 s flight the model loses ~6% of spin where Lorenz implies ~25%. *Acceptance:*
if the still-air 25 m validation degrades, the discrepancy may live in Lorenz's disc
(204 g, 15% heavy, visibly wobbling — and Weizman found wobble increases spin-down)
rather than in `C_Nr`. Also reportable. **Note this bears on WP2/P5**: the high-ω end of
any spin sweep currently holds spin unrealistically.

**WP6.5 — Free inertia check (~10 min).** The inertia ratio predicts a nutation period
of `2π/((I_s/I_p)ω)` = 0.0796 s at ω = 41 rad/s against Lorenz's measured ~0.09 s
wobble. 12% agreement from a purely inertial quantity, and it has never been run. Cheap
independent validation the model probably passes.

---

### WP7 — Gated / deferred (do not start without the gate)

**WP7.1 — AIJ Case B cube/high-rise benchmark (~4 h).** Only if any site content is going
somewhere refereed. Every quantitative statement now rests on **modelled k in a separated
wake**, which is the best-documented RANS failure mode in this literature — and the
benchmark designed to expose it was not run. *Acceptance:* mean U within 20% and k within
a factor 2 over `0 < z/H < 1.5` at 1.25H downstream. **A failure is still publishable —
report the k bias factor and carry it as a stated systematic on every σ_u in the paper.**
That converts "we used RANS" into "we know how wrong our k is, by this factor."

**WP7.2 — Third mesh, refining the *source* (~8–14 h, one direction).** Not another
`fieldBox` level. Raise refinement on the Campus Towers and Mertz **surfaces** and in a
wake box behind Campus Towers, targeting ≥10 cells per building side (currently 3–6
against a 23.5 m building). With three meshes you finally get an observed order `p` and a
GCI instead of a max-difference. *Acceptance:* `p ∈ [0.5, 2.5]` **and** GCI_fine < 20% of
local speed **and** < 5% of `u_ref` in the disc band at all 5 stations, reporting **k
alongside U**. If `p` falls outside that range the solution is not in the asymptotic
range — report the two-mesh band with `F_s = 3` and state that no direction ranking is
possible on this mesh family.

**WP7.3 — The generalizing parametric sweep (~6–12 h). GATE G4: run only if the
corrected law predicts a non-monotone σ_u response to porosity.** If sheltering
monotonically helps, there is no optimum and the paper is "shelter more" — don't spend
the compute. Design, if the gate opens:

| axis | levels | note |
|---|---|---|
| porosity β | 0, 0.20, 0.35, 0.50 | **the primary axis** — the shelterbelt literature is unambiguous that solid barriers buy a quiet zone with a violent bleed layer while porous ones don't. Sean Earl is the β = 0 endpoint |
| approach `z0/H` | 1e-4, 1e-3, 1e-2 | open water → suburban → urban at H = 50 m. **Must be in the sweep** — it is the exact variable the control error was confounded with |
| `W/H` | 1 (primary) + 4 (2 extra runs) | 3-D wake vs quasi-2-D fence. Sean Earl is W/H ≈ 0.8, strongly 3-D; a 2-D fence sweep would not represent it |
| matched empty controls | one per `z0` | **non-negotiable** given C2 |

4×3 + 2 + 3 = **17 runs**. Hold fixed: `z_obs/H`, blockage < 3% (a numerical artefact
parameter to be controlled, not studied), `Re_H > 1e4`, domain in H units, 20H
downstream. **Do not sweep `H/D` as a case parameter** — `x/H` is a sampling coordinate,
so one run in a 20H domain gives the whole curve and 11 of the memo's 12 cases are
wasted.

Plot everything against `x/H ∈ [1, 20]`, one curve family per β, one panel per `z0`:
mean deficit `U/U₀`, **`σ_u/σ_u,₀`** (the quantity that decides the paper), `TI/TI₀`, and
`Π_play/Π_play,₀`. Normalising by the matched control is also what makes RANS defensible
here — a systematic k bias largely cancels between case and control. **Master figure:**
`min over x of Π/Π₀` as a surface over `(β, z0)`.

**Set the falsification criterion before running:** an optimum exists iff
`min_x Π/Π₀ < 0.9` for at least one `(β, z0)`. If `Π/Π₀` rises monotonically from 1,
"optimal shelter" is dead — **and that is equally publishable**, provided you said so
beforehand.

**WP7.4 — LES.** Deferred indefinitely. If ever run: the deliverable is `L_11`, `L_22`
and `σ_w/σ_u` in the disc band, *not* gust structure. `√Φ` varies by ×3 over the plausible
range of L and `κ` by ×1.6 over the plausible range of `σ_w/σ_u`; peak gusts move
nothing. Realistic cost is 3–8 days (~5e6 cells, ~12.5 GB — tight in 15 GB; ~7e4
timesteps at CFL<1), not the "~2 weeks" in the status memo. It also needs a matched-z₀
synthetic-inflow precursor, so **WP3 is a prerequisite**.

---

## 4. Prior-art check — do this before writing anything (≈2–3 h)

**GATE G5.** The novelty claim in `FRISBEE_GENERALIZATION.md` §1 ("Nobody has written
down the equivalent criterion for a projectile") is too strong and a referee will
challenge it. Retrieve and read:

1. **Curcin, V. & Harris, R.I. (1986)** — the earliest candidate for Monte-Carlo of a
   body through synthetic turbulence. **If it contains the scaling law, Paper 1's
   novelty collapses to the lifting correction only.** Close this first; it is cheap and
   it changes the framing, not the work.
2. **Etkin, B. (1981)** and **Hoblit, F.M., *Gust Loads on Aircraft*** — the sweeping
   correction and the von Kármán continuous-turbulence model. Cite in the
   **introduction**, state that the correction is standard in aerospace and absent from
   the sports / debris / urban-UAV literature, and confine the novelty claim to the
   closed-form **lifting-body landing variance**. **If you cannot articulate what `Φ(N)`
   adds beyond a von Kármán response integral, there is no theory paper.**
3. **CAP 437** and **CAA Paper 2008/03** — the helideck σ_w criterion. Reframe per C6:
   you are deriving it, not competing with it.
4. **Kordi & Kopp (2012)** and the windborne-debris literature — landing dispersion for
   non-lifting bodies.
5. **Giljarhus, Gooding & Njærheim (2022), *Sports Engineering* 25:26** — open access,
   CFD coefficients + 6-DOF + measured release conditions + video-tracked trajectories
   with per-quantity uncertainties. **The closest published analogue of this entire
   exercise, and currently uncited.** It establishes what "good" means here (their model
   under-predicted range by 13–19 m, <20% error), and they set the side force to exactly
   zero and say so — useful precedent *and* useful contrast.
6. **Crowther & Potts (2007), *Sports Engineering* 10(1):3–21** — an independent 6-DOF
   disc model with **published, directly comparable results**: optimum launch pitch ≈ 10°
   for range and ≈ 20° for duration, straightest flight at initial roll −6°, trajectory
   similarity governed by static-margin-to-advance-ratio. Four cheap regression targets,
   none currently used. **This is the model this paper is closest to, in the journal it
   is being submitted to. Omitting it is the most exposed citation gap in the project.**
7. **Sasakawa & Sakurai (2008), *Sports Biomechanics* 7(3):311–321** — 17 players,
   release speed statistically indistinguishable between skilled and unskilled (21.7 vs
   20.7 m/s) but **spin 81 vs 59 rad/s and distance 51.4 vs 29.5 m**. *This is a direct,
   measured, human-throw test of exactly the spin sensitivity the gyroscopic hypothesis
   rests on:* hold V ≈ 21 m/s, vary ω from 59 to 81 rad/s, and the model must produce
   ~30 m → ~51 m. **Highest-value single validation test available, and it is free.**

Also worth having: Potts & Crowther (2002) AIAA 2002-3150 and Potts (2005) PhD thesis
(the primary source — this model currently uses *Hummel's transcription* of it, not the
source; both closed-access, acquisition is the gating step); Hubbard & Hummel (2000) for
two more fully-specified measured flights and the `C_Mα̇` lag term this model omits;
Weizman et al. (2020) for spin-down torques; Winograd & Engsberg (2012) for the
backhand/forehand precession difference and the finding of **no accuracy difference
between throw types**; Kamaruddin, Potts & Crowther (2018) for geometry sensitivity;
Stilley & Carstens (1972) for an independent coefficient set that would give the
model its first real uncertainty bracket (currently every variant traces to one
transcription).

**Warning:** the "C_S ≈ 0.8 at AdvR = 1.04" figure circulating via Lorenz's *Spinning
Flight* is almost certainly a typo for 0.08. Do not cite it without the primary.

**Gaps that genuinely do not exist — state these as such rather than searching
indefinitely:** no radar/Trackman/motion-capture trajectory dataset for Ultimate exists;
no measured hammer release parameters exist anywhere; no wind-tunnel yawing moment for a
*spinning* disc exists (the rig spins about the axis in question), so every published
`C_Nr` comes from flight-data parameter estimation, CFD, or IMU spin-decay; and disc
"flight numbers" have never been validated against measurement.

---

## 5. Publication plan (context for why the packages are ordered this way)

Not one paper. The panels converged on a decomposition:

**Paper 1 — write first.** *The wrong denominator: turbulence criteria for projectiles in
the built environment.* Core claim: a projectile supplies its own advection velocity, so
sensitivity scales with `σ_u/V_d`, not `σ_u/U`; landing dispersion follows a closed-form
path-averaging law; gust models built on an Eulerian correlation time overpredict by
`√N_eddy`. Content is WP2's homogeneous-turbulence validation — **no site CFD at all**,
with Sean Earl as a one-paragraph worked example. Venue: JWEIA primary (where the debris
and urban-UAV people making this error actually publish), *Physics of Fluids* fallback.
**Not JFM** — this is a gust-response calculation, not fluid mechanics, and a JFM referee
will say so.

Put the WP1/WP2 coincidence in the paper explicitly: two independent analyses, one from
the code and one from the theory, predicted the same 2.3–2.6× correction. That is the
strongest evidence the physics is right.

**Paper 2 — write in parallel, it is cheap and it is the most durable object here.** *A
validated open-source 6-DOF flying-disc model, with corrections to the Hummel (2003)
coefficient set.* No published errata to Hummel exist and Giljarhus et al. (2022), the
current state of the art, does not correct her coefficients. Venue: JOSS + a short
technical note in *Sports Engineering* (Giljarhus is the obvious referee and a friendly
one — the corrections improve his stack). **Blocking requirements: LICENSE, git history,
`pyproject.toml`, the 26-test suite actually present and runnable, coefficient-table
redistribution rights checked.** None of that exists today. This establishes priority on
the five corrections before someone else finds them, and it costs a week.

**Paper 3 — fold in or preprint only.** The Sean Earl site study is not publishable
standalone after C2. Either demote it to Paper 1's demonstration section — where it earns
its place by supplying a real σ_u and L at a real site, and by being the case where TI
and `σ_u/V_d` rank the site oppositely — or post it as a clearly-labelled case study and
move on.

**Paper 4 — conditional on G4.** *Optimal shelter for a projectile differs from optimal
shelter for erosion control.* The shelterbelt literature is fully mature on the
mean/turbulence tradeoff (Wilson 1987; Cornelis & Gabriels 2005) and you cannot claim it.
What is unclaimed is that **the objective function is different**: erosion control
minimises U, playability minimises σ_u, and those give different optimal porosity and
standoff.

**Abstract sentence to aim at:**

> For a projectile the governing turbulence measure is σ_u/V_d, not the turbulence
> intensity σ_u/U that underpins pedestrian-comfort criteria — so a sheltered site whose
> turbulence intensity triples is nonetheless easier to throw on — and the resulting
> landing dispersion obeys a closed-form path-averaging law that gust models built on an
> Eulerian correlation time overpredict by √N_eddy.

---

## 6. Risk register — flag these to the author, don't try to solve them in code

1. **Prior art in aircraft gust-load theory** (C6). Largest single threat to novelty.
   Mitigated by WP4's prior-art gate, not by more compute.
2. **The corrected result may be undramatic, and partly is.** `V_d/σ_u = 17–36` means
   `σ_u/V_d ≈ 0.03–0.05` — **the disc is insensitive.** The project's founding premise
   ("outdoor disc sports are dominated by wind") is not supported at the site that was
   chosen, because a sheltered site was chosen. The *theory* survives this; a law is a
   law. The *site study* does not. **Do not try to rescue the site study by finding
   something dramatic in it.**
3. **Nothing is measured.** Zero field data across the project. A solo preprint on
   frisbee physics with a CFD that fails grid independence and no measurement is a
   citation risk regardless of theory quality. Two cheap mitigations: one afternoon with
   a sonic anemometer at 2 m on the pitch (validates σ_u, the quantity the whole argument
   now rests on — the cheapest credibility purchase available); and better, **a disc with
   an IMU/GPS logger thrown ~200 times in measured wind.** Nobody has published measured
   disc landing dispersion against measured σ_u. That dataset would make every one of
   these papers referee-proof and is worth more than any additional compute in this
   document.
4. **σ_land ≈ 6.9 m at 25 m range is not a model of an Ultimate player** (implied by
   S̄ ≈ 0.09 at r = 3 m; club players complete 25 m backhands well above 90%). Since
   `p ∝ σ⁻²`, a 2× error in σ_land is a 4× error in Q — and it *suppresses the signal*:
   against σ_land = 6.9 m the synthetic turbulence accounts for **~0.8% of landing
   variance.** The variance decomposition the design was built to support has a
   predetermined answer, and it is an artefact of the nuisance model. Calibrate against
   *any* real data — even 50 throws with a tape measure.
5. **A wind-engineering co-author would halve the referee risk** on Papers 1 and 4.
6. **The hypothesis tests are not tests.** `check_hypotheses.py` decides verdicts by
   argmax with a 1e-9 tolerance. **H2a is decided by a difference of 1.8e-4** (CV 0.25221
   vs 0.25203 — 270° loses by 0.07% and is stamped CONTRADICTED). H2b's gap is ~2.1σ
   against MC noise alone and ~0.3σ against physics uncertainty. H1a takes `np.gradient`
   of a 4-point profile then the maximum — a noise amplifier by construction — against an
   undocumented 0.5 threshold. **None of the five verdicts is currently a test.** Replace
   with pre-registered effect sizes and paired CIs, or report them as descriptive.

---

## 7. The honest counterfactual, for the author to consider

The strategy panel's closing argument, recorded here because it is a real option and not
obviously wrong:

**The disc may be the wrong object.** `σ_u/V_d ≈ 0.03–0.05` puts it in the *insensitive*
corner of its own parameter space. The regime where this physics bites is
`σ_u/V_d → O(1)` — a slow vehicle: a 5–15 m/s delivery quadrotor, and in the limit a
hovering one, where `V_d → 0` and the criterion diverges. The urban-air-mobility
community has produced a wave of CFD risk maps since 2022, **has no turbulence
criterion**, is openly looking for one, and every current regulator wind limit is a mean
or gust speed with no turbulence content at all.

The high-impact version is the same theory aimed at urban UAV operations, with the
frisbee as the **validation case** — the disc has published wind-tunnel coefficients, a
validated 6-DOF model, and a landing point measurable with a tape, making it ideal ground
truth for a law you then apply where measurement is expensive and the stakes are
regulatory. Same theory, same code, same effort, and a plausible path from "a paper" to
"a criterion that gets adopted."

The low-risk version is Paper 2 alone: ship the validated disc model as a proper
open-source package with its errata. One week, and it will outlive everything else here.

What is *not* recommended is spending another month of compute making the Sean Earl case
study rigorous. It has already given what it had to give — a real σ_u, a real L, and the
discovery that TI misleads.

---

## 8. Execution checklist

```
[x] WP0  free artifact checks                    DONE 2026-07-26 -- see results in §3
      [x] WP0.1 blocking: 1.32 pp differential, does NOT explain the ranking
      [x] WP0.2 n_valid = 100% everywhere, concern closed
      [x] WP0.3 NEW DEFECT: landing_mean/ellipse_axes in the wrong target frame
      [ ] WP0.4 record realised T_L, RNG/version provenance
      [ ] WP0.5 fix tab:Qbyheading S-bar column
      [ ] compare the three per-heading `solid` masks
[ ] WP1  gust sweeping timescale                  4 h + 3 min run
      [ ] GATE G1: dispersion falls 2.0-2.6x?  -- if <20% change, STOP AND REPORT
[ ] WP2  validate the dispersion law              1-2 days, no CFD
      [ ] drop i_h from spawn_key first (one line)
      [ ] P5 impulse test first (45 throws, 2 min)
      [ ] GATE G2: law within ~20%?  -- if no, STOP AND REASSESS THE PROGRAM
[ ] WP4a prior-art gate (§4)                      2-3 h
      [ ] GATE G5: what survives Etkin / Curcin & Harris / CAP 437?
[ ] WP3  matched-roughness controls                1.5 h compute
      [ ] GATE G3: retract "shelter is not protection" if confirmed
      [ ] grid-check k as well as U (no compute)
[ ] WP4b retractions + provenance (§WP4)           1 day, no compute
      [ ] NOTHING IS POSTED UNTIL THIS IS DONE
[ ] WP5  replace Q with F and P_hit(Pi_play)       4 h
[ ] WP6  flight-dynamics hardening                 6 h
[ ] WP7  gated: AIJ-B / third mesh / porosity sweep / LES
      [ ] GATE G4: does the corrected law predict non-monotonicity?
```

**Total to the end of WP5: roughly one working week, of which under three hours is
compute.** The current 97,200-throw sweep cost one hour and addressed ~2% of the output
variance; the CFD and the thrower model, which set ~85% of it, have never had a single
sample allocated to them. Whatever else changes, invert that ratio.

**Report back after every gate.** Gates G1, G2 and G5 can each end the program in its
current form, and that is a legitimate outcome — the point of putting them first is that
they cost hours instead of weeks.
