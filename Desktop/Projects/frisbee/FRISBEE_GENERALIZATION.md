# Generalizing the Sean Earl result: what makes this more than a field report

**Status:** analysis memo, 2026-07-26. Numbers here trace to
`P3_CFD/vertical_structure.json`, `P3_CFD/vertical_fields.npz` and the paired
release-height sweep in `P4_Coupling/height_sweep.json`.

---

## 1. The finding that generalizes

> ~~**Shelter is not protection.** Buildings near a playing field remove *mean*
> wind and add *turbulence*. For a thrown object, mean wind is a bias the
> thrower can aim off; turbulence is variance they cannot. A "sheltered" pitch
> can therefore be systematically *harder* to throw on than an exposed one,
> while feeling calmer to a person standing on it.~~
>
> ⚠ **Withdrawn 2026-09-04.** The "add *turbulence*" half is not supported by these
> cases: it rests on two headings whose inlet roughness differs from the control's, and
> against matched approach flow the buildings *reduce* σ_u everywhere (see the box
> below). What replaces it is narrower but survives:
>
> > **Turbulence intensity is the wrong diagnostic for a thrown object.** Buildings
> > remove mean wind, which inflates TI whether or not the fluctuation the object
> > actually responds to has grown. At the one roughness-clean heading here, TI rises
> > 2.08× while σ_u *falls* to 0.91×. For a thrown object the group that matters is
> > σ_u/V_d, and a "sheltered" pitch cannot be graded by how calm it feels — or by TI.
>
> Whether shelter can be *anti*-protective in the original sense remains open and is
> worth one matched control run.

> ⚠ **THE σ_u COLUMN IS CONFOUNDED — corrected 2026-09-04.** The control is *not* an
> "identical inflow" run for two of the three headings. Verified directly from
> `case_info.json` by `check_control_matching.py`:
>
> | case | `z0_inlet` | `z0_ground` | `u*` | `k_inlet_ref` | matched to control? |
> |---|---|---|---|---|---|
> | `dir_090_empty` (control) | 0.001 | 0.03 | 0.2226 | 0.1651 | — |
> | `dir_090` | 0.001 | 0.03 | 0.2226 | 0.1651 | **yes** |
> | `dir_210` | **0.5** | 0.03 | **0.6733** | **1.5113** | **no** |
> | `dir_270` | **0.5** | 0.03 | **0.6733** | **1.5113** | **no** |
>
> Ground roughness is matched; **inlet** roughness is not. The 210° and 270° cases meet
> inflow carrying `sqrt(1.5113/0.1651) = 3.03×` the control's fluctuation *before
> reaching a single building*. Dividing that out:
>
> | θ_wind | memo σ_u ratio | approach flow | buildings' own effect |
> |---|---|---|---|
> | 90° | 0.91 | 1.00 | **0.91** (clean) |
> | 210° | 1.51 | 3.03 | **0.50** |
> | 270° | 1.57 | 3.03 | **0.52** |
>
> **Against matched approach flow the buildings *reduce* σ_u in all three directions.**
> The "raised by half again" claim came entirely from the two unmatched headings, and the
> one roughness-clean heading (90°) is the one where the ratio was already below 1.
> **This claim cannot be made from these cases.** Settling it needs one more run: an
> empty-domain control at `z0_inlet = 0.5`. Raised in
> `CLAUDE_CODE_FRISBEE_HARDENING.md` §C2 and left standing by the 2026-09-03 pass.

Measured at 2 m over the playing surface, against an empty-fetch ABL control
(`--control abl_control.npz`) — **roughness-matched at 90° only**, see the box above:

| θ_wind | mean speed vs control | σ_u vs control | TI pitch | TI control |
|---|---|---|---|---|
| 90° easterly | **0.44×** | 0.91× | 0.25 | 0.12 |
| 210° SSW | **0.50×** | ~~1.51×~~ → 0.50 | 0.36 | 0.12 |
| 270° westerly | **0.50×** | ~~1.57×~~ → 0.52 | 0.38 | 0.12 |

The buildings cut mean wind roughly in half, and turbulence intensity went from 12% to
25–38% — a factor of 2–3. **What survives is the diagnostic claim, not the shelter
claim.** At 90°, the only roughness-clean pair, σ_u *falls* to 0.91 while TI *rises* by
2.08×, and the rise is arithmetically the denominator: `σ-ratio / U-ratio = 0.91/0.44 =
2.07` against a measured TI ratio of `0.25/0.12 = 2.08`. Dividing by a collapsing mean
manufactures a factor of two where the disturbance the disc responds to has slightly
decreased.

So the transferable claim is: **turbulence intensity cannot tell you whether shelter
helped or hurt a thrown object; the group that can is σ_u/V_d.** Whether shelter ever
raises σ_u outright is a well-posed question these cases do not answer. This is *not*
what urban wind engineering optimises for either: the COST 732 / Lawson / NEN 8100
tradition is built around **pedestrian** comfort and safety, which is a
mean-and-gust-exceedance criterion for a person standing still.

> ⚠ **NARROWED 2026-09-02 by GATE G5** (`GATE_G5_PRIOR_ART_2026-09-02.md`). This
> paragraph used to end *"Nobody has written down the equivalent criterion for a
> projectile."* That is refutable from a single search. The urban-air-mobility
> literature has built exactly such criteria since about 2020 — site-resolved,
> CFD-derived thresholds on wind speed, shear and turbulence intensity for flying
> bodies (Watkins *et al.*, *Int. J. Micro Air Vehicles* **2**(2):55–68, 2010;
> *Drones* **9**(5):362, 2025, `10.3390/drones9050362`; and two 2026 urban-air-mobility
> reviews). **Do not claim the gap; claim the distinction.**
>
> The distinction is real and it is sharper than the original claim. Every one of
> those criteria is a **threshold on a controlled vehicle** — a go/no-go on whether
> a drone can hold station or a pilot can hold a path. A thrown disc is **unguided
> after release**, so there is no control loop to threshold: the quantity of
> interest is the **variance of the terminal position**, and the organising group is
> `σ_u/V_d` (fluctuation against the body's own airspeed) rather than the turbulence
> intensity `σ_u/U` a station-keeping criterion uses. That is the same argument §2
> below makes independently from `V_d/U ≈ 8–10`.

## 2. Why the disc cares — the dimensionless framing

This is what turns a site report into a general result. Three scales matter:

| symbol | meaning | value at Sean Earl, disc band 1–6 m |
|---|---|---|
| `T_f` | flight time | ≈ 2.5 s (25 m backhand) |
| `L_f` | flight length | 25 m |
| `T_L = k/ε` | eddy turnover time | **10.5–18.0 s** |
| `L = σ_u·T_L` | integral length scale | **5.7–6.7 m** |
| `V_d` | disc airspeed | 10–14 m/s |
| `U` | mean wind | 1.2–1.5 m/s |

(Band averages over 1–6 m from `vertical_structure.json`; at 2 m specifically,
`T_L` = 9.5–14.4 s and `L` = 4.4–5.1 m.)

Two ratios follow, and they place this problem in a specific regime:

- **`T_L / T_f = 4.2–7.2`** — the eddies barely evolve during a throw. The disc
  flies through a *frozen* turbulent field.
- **`N_eddy = L_f / L = 3.7–4.4`** — but it crosses about four independent
  eddies on the way.

Note `V_d / U ≈ 8–10`: the disc, not the wind, does the sweeping. This is
Taylor's frozen-turbulence hypothesis with the *projectile* supplying the
advection velocity. The practical consequence is partial averaging: scatter
scales not as the quasi-steady limit `σ_u·T_f` but roughly as

```
σ_landing  ~  σ_u · T_f / sqrt(N_eddy),      N_eddy = L_f / L ≈ 4
```

**This immediately generalizes across sports** by asking where each object sits
in `(T_L/T_f, L_f/L)` space:

| object | `V_d` | `L_f` | regime | expected sensitivity |
|---|---|---|---|---|
| badminton shuttle | slow, decelerating | ~10 m | few eddies, low speed | **worst** |
| ultimate disc | 10–14 m/s | 25 m | ~4 eddies, frozen field | high |
| soccer ball (free kick) | 25–30 m/s | 25 m | ~4 eddies, faster | moderate |
| golf drive | 60–70 m/s | 230 m | many eddies, strong averaging | low per-eddy, but long lever |
| delivery UAV | 5–15 m/s | 100s m | slow, long | **worst**, and it must *hover* |

The framework is object-agnostic: site CFD → integral scales at operating
altitude → object flight scales → dimensionless ratio → dispersion estimate.

## 3. Why throwing low wins twice

The vertical profiles say the disc band is the best place to be, for two
independent reasons that compound:

1. **σ_u increases with height.** Relative to 2 m, σ_u is 1.2× at 6 m, 1.3× at
   10 m, 1.6× at 20 m — in *all three* directions. More energy aloft.
2. **The integral scale grows with height** — `L` ≈ 4.4–5.1 m at 2 m rising to
   ~8–9 m at 6 m and ~10 m at 13 m, close to the classic surface-layer `L ~ z`.
   Larger eddies mean *fewer independent samples* per flight, so less averaging.

So throwing higher buys both a stronger disturbance and less cancellation of
it. The ground is not merely calmer; it is calmer *and* more forgiving.

Above that sits the real hazard. The θ = 90° centreline slice
(`fig6_centreline_slice.png`) shows the Campus Towers separating shear layer as
a band of TI ≈ 0.8–1.0 spanning **z ≈ 15–50 m directly over the pitch**, peaking
at TI = 0.86 at the 45 m roof height. **The disc band sits underneath it.** A
disc lofted into that layer — a high hammer, a skied huck, a blade — leaves the
quiet ground-hugging layer and enters flow whose fluctuation is three times
stronger. That is a concrete, physically-grounded piece of coaching advice, and
it is exactly the kind of statement a mean-wind-only analysis cannot make.

## 4. Deeper questions worth probing

Ordered by (impact × tractability). Items 1 and 2 would, in my judgement, turn
this from a case study into a paper with legs.

### 1. Is there an *optimal* amount of shelter? ★ highest value
Exposed: large correctable mean, modest TI. Deeply sheltered: small mean, large
TI. Both extremes are bad; the interior of that trade may have an optimum.
**Probe:** idealised parametric CFD — a single block of height `H` at distance
`D` upwind of a pitch, sweeping `H/D` over ~12 cases. Sean Earl becomes the
validation case study; the sweep is the general result. This is the paper that
gets cited by stadium designers, and it is *cheap* — idealised geometry meshes
in minutes.

### 2. A playability criterion to sit alongside pedestrian comfort ★
Formalise Q as a dimensionless index, e.g. `Π = σ_u·T_f / (r_target·sqrt(N_eddy))`,
and propose threshold bands the way Lawson did for pedestrians. Wind engineering
has no projectile-sport standard. Writing the first defensible one is a
genuinely high-impact move, and it makes every future site study comparable.

### 3. ~~Does gyroscopic stiffness buy gust immunity?~~ — ANSWERED, and the answer is no

> ⚠ **Rewritten 2026-09-03 (WP4b item 7). This was still posed as an open hypothesis
> after it had been measured and refuted.** The probe below was run on 2026-08-30
> (`THEORY_WIND_SPEED_2026-08-30.md`, §"Spin is not the explanation") and the memo was
> never updated.
>
> ⚠ **The hardening brief's own restatement of this item is also wrong and must not be
> copied.** It says to report "no threshold, `ω*` inside the human range, total scatter
> increases with spin, cross-range scatter falls as ≈ `ω^−0.21`". There is no `ω*` in the
> measurement and no negative exponent anywhere in it. Use the numbers below.

**The original hypothesis.** A disc is a gyroscope; a gust's rolling moment precesses into
a *heading* change a quarter cycle later, so higher spin ⇒ stiffer ⇒ smaller deflection per
gust, and there should be a critical spin rate above which turbulence is materially
suppressed. If it held, it would explain why coaches teach "more spin" in wind, which is
currently folklore.

**The measurement.** Release spin swept 50 → 200 rad s⁻¹ at two headings and two wind
speeds, everything else held and the release draws **paired**:

| wind from | `Uref` | fitted `σ_gust ~ spin^p` | 200 vs 50 rad s⁻¹ |
|---|---|---|---|
| 90° | 5 | +0.02 [−0.01, +0.06] | 1.03× |
| 90° | 10 | +0.03 [+0.02, +0.05] | 1.05× |
| 270° | 5 | +0.07 [+0.01, +0.15] | 1.12× |
| 270° | 10 | +0.07 [−0.03, +0.17] | 1.10× |

**A fourfold change in spin moves the gust dispersion by 3–12%, and in the wrong
direction** — more spin gives very slightly *more* scatter, not less. There is no
threshold, no knee, and no immunity. **The folklore is not reproduced.**

**The caveat that has to travel with the null.** Spin enters this model only through the
rate-damping and gyroscopic terms of `disc_6dof`. If those under-represent the real
coupling, the null is a property of the model rather than of a disc. What can be said
without qualification is that spin is not a free parameter that reconciles the ballistic
dispersion law with the 6-DOF model it is tested against.

⚠ **Do not upgrade this to "spin does not matter in wind."** Sasakawa & Sakurai (*Sports
Biomechanics* **7**(3):311–321, 2008) measured 17 players at statistically
indistinguishable release speeds (21.7 vs 20.7 m s⁻¹) and found distance 51.4 vs 29.5 m
with spin 81 vs 59 rad s⁻¹. Spin matters enormously for *range*. This null is about
**gust-induced scatter at fixed flight**, which is a different question, and running that
measured comparison through this model is the outstanding validation test (GATE G5, §4).

### 4. Do throw types re-rank under turbulence?
The hypothesis check already found forehand beating backhand at *every* heading
(+0.018 to +0.032), which contradicted the draft's westerly-specific prediction.
Is that a mean-flow effect or a turbulence-robustness effect? The paired
gust/no-gust machinery built for the height sweep answers this directly at no
new cost.

### 5. Porous windbreaks: can you engineer playability?
Agricultural shelterbelt literature is mature: **porous** barriers reduce mean
wind with far less shed turbulence than solid ones. Importing that to sports-
field design is a real cross-disciplinary contribution — and it is the
actionable answer to question 1 if the optimum turns out to be "shelter, but
not with a solid tower".

### 6. Is the optimal release height a flow property or a dynamics property?
If the optimum tracks the shear-layer base, it is site-specific. If it tracks
the disc's own flight envelope, it is universal. Distinguishing them requires
one more site — the cheapest real generalisation test available.

### 7. Null result worth reporting: shear is not the hazard

**The conclusion is right. The argument was not, and it is replaced here** (rewritten
2026-09-03, WP4b item 8).

> ⚠ **The retired argument.** It compared the sheltered pitch's 0.02–0.06 s⁻¹ against the
> control's "up to 1.5 s⁻¹". That 1.5 s⁻¹ is the control's value at z = 0.5 m — *below*
> the disc band — and the control's near-ground profile is a wall-function interpolation
> anyway. Comparing a disc-band number against a below-band interpolated one is not a
> comparison. It also quietly dropped the sign: the 270° shear over the disc band is
> **−0.018 s⁻¹**, i.e. the profile is *inverted* there, which a range written as
> "0.02–0.06" hides.

> ⚠ **The second retired argument (retired 2026-09-04, one day after it was written).**
> The replacement below was itself wrong, and a physics referee caught it. It compared
> `σ_u` against `ΔU = (∂U/∂z)·d` and called the quotient a ratio of *rolling moments*:
>
> | θ_wind | `∂U/∂z` (s⁻¹) | `ΔU = (∂U/∂z)·d` (m/s) | `σ_u` (m/s) | `σ_u / ΔU` |
> |---|---|---|---|---|
> | 90° | +0.024 | 0.0065 | 0.34 | ~~53×~~ |
> | 210° | +0.057 | 0.0153 | 0.56 | ~~37×~~ |
> | 270° | −0.018 | 0.0048 | 0.59 | ~~122×~~ |
>
> **A rolling moment comes from a velocity difference across the span.** `ΔU` is such a
> difference; `σ_u` is not — it is the RMS of the velocity *at a point*. A spatially
> uniform gust changes the disc's lift and drag and exerts **no roll couple at all**. The
> table therefore divided a couple by a point velocity and labelled the result a ratio of
> moments. The two quantities are not the same kind of object.
>
> Note the failure mode: this table was written *to fix* an earlier like-for-like problem,
> and reproduced it one level down. Getting both numbers from the same file and the same
> height band made them look comparable without making them comparable.

**The argument that works.** Put both channels on the same footing — the velocity
difference each imposes **across the disc diameter** `d = 0.269 m`
(`P1_AnalyticalModel/coefficients.py:516`, Potts & Crowther Table 3-1):

- shear: `ΔU_shear = |∂U/∂z|·d` — a **mean** gradient, smooth by construction, linear in `d`
- gust: `ΔU_gust` — a **turbulent** velocity increment, whose scaling depends on the regime

> ⚠ **A SECOND RETRACTION, 2026-09-04, and it goes back the other way.** The version of
> this section written on 2026-09-03 used the smooth-flow increment `σ_u·(d/L_int)` and
> defended it here as *"the standard eddy-gradient estimate"* with *"an unmeasured
> order-one prefactor"*. **The unmeasured thing was the exponent, not the prefactor.**
> That form holds only for a body smaller than the Kolmogorov scale. It gave
> ~~**1.7 to 4.9×**~~ and every number in the retired table below.
>
> | θ_wind | `ΔU_gust` (m/s), smooth | ratio, smooth |
> |---|---|---|
> | 90° | ~~0.014752~~ | ~~2.3×~~ |
> | 210° | ~~0.026428~~ | ~~1.7×~~ |
> | 270° | ~~0.023688~~ | ~~4.9×~~ |
>
> **The regime is decided by a computed number, not by assertion.** With
> `ε ~ σ_u³/L_int` the Kolmogorov scale is `η = (ν³/ε)^(1/4) = 0.58–0.85 mm` against
> `d = 269 mm`, so **`d/η` = 315 to 468**. The disc is firmly in the **inertial subrange**
> and the applicable estimate is the second-order longitudinal structure function,
> `ΔU_gust = √C₂·σ_u·(d/L_int)^(1/3)` with `C₂ ≈ 2.0` (Sreenivasan 1995; Pope, *Turbulent
> Flows*, ch. 6). The smooth limit understates it by `(L/d)^(2/3)` ≈ **10.8 to 12.1×**.
>
> Note the direction. The **original** 37–122× — derived wrongly, by dividing a couple by
> a point RMS — was *nearer the truth* than the 2026-09-03 "correction". Being closer to
> the answer is not the same as being right, and it is not a defence of the first table
> either. Both retractions stand.

`shear_vs_gust_moment.py` now prints `d/η` beside the answer so the regime choice can be
checked rather than trusted, and prints all three ratios so the size of each retraction is
visible. All inputs are disc-band quantities from `vertical_generated.tex`; regenerate
with `python shear_vs_gust_moment.py`.

| θ_wind | `∂U/∂z` (s⁻¹) | `σ_u` (m/s) | `L_int` (m) | `η` (mm) | `d/η` | `ΔU_shear` (m/s) | `ΔU_gust` (m/s) | ratio |
|---|---|---|---|---|---|---|---|---|
| 90° | +0.024 | 0.34 | 6.2 | 0.854 | 315 | 0.006456 | 0.168959 | **26×** |
| 210° | +0.057 | 0.56 | 5.7 | 0.575 | 468 | 0.015333 | 0.286195 | **19×** |
| 270° | −0.018 | 0.59 | 6.7 | 0.576 | 467 | 0.004842 | 0.285710 | **59×** |

**The gust-induced rolling moment exceeds the shear-induced one by 19 to 59×.**

**What that does to the null.** It **survives with room to spare**: more than 10× at every
heading, so the order-one prefactor in the structure-function estimate cannot reverse it.
"Shear is not the span-wise hazard at this site" is again a defensible statement — but see
caveat 2, which was the real finding of the 2026-09-03 pass and is untouched by either
retraction.

⚠ **Three caveats.**

1. **This is the *sheltered* case.** In open fetch `∂U/∂z` at 2 m is an order of magnitude
   larger, which leaves gust ahead by only **1.9 to 5.9×** there — a single-figure margin
   the order-one prefactor *can* reach at the weakest heading. The null is comfortable at
   this site and marginal away from it.
2. **There is a second shear channel and it is not in the table.** The mean wind changes
   over the trajectory's own height excursion, which is *metres* rather than the 0.27 m of
   the disc; in the inertial-subrange scaling that channel beats the span-wise one by
   roughly `(Δz/d)`, not `(Δz/d)^(1/3)`, because it is a **mean** gradient and stays
   linear. It is **deterministic and therefore learnable**, exactly like mean wind, so it
   remains a bias a thrower can aim off rather than a variance they cannot. **This is the
   part of the 2026-09-03 pass that was right and that neither retraction touches.**
3. ⚠ The hardening brief states this ratio as "~20×", using `σ_u = 0.36`. That 0.36 is the
   **turbulence intensity** at 210°, not `σ_u`; the actual `σ_u` there is 0.56 m/s. The
   brief's ~20× is *also* not moment-consistent, so it lands near the right answer for two
   wrong reasons at once. Use the table above, not the brief.

**What to report.** That **at this sheltered site the span-wise gust exceeds the span-wise
mean shear by an order of magnitude and more**, that the ordering narrows to single figures
in open fetch, and that the **trajectory-scale** mean shear — absent from every version of
this table — is the larger deterministic term and the one worth modelling. The history of
this section is itself worth a line in the paper's methods: two retractions, in opposite
directions, on a quantity that was never measured but only estimated.

## 5. What would make the current paper materially stronger

- **LES for one direction.** The gust field is currently *synthesised* from k by
  a spectral model, not resolved. Every statement about gust structure inherits
  that assumption. ~3–7 days for a nested domain (memory-bound at ~5M cells on
  15 GB). The easterly case is where RANS is least trustworthy and where the
  shear layer matters most.
- **One afternoon of field measurement.** A single sonic anemometer at 2 m on
  the pitch, in a known synoptic wind, would validate σ_u — the quantity the
  whole argument now rests on — and is the cheapest credibility purchase
  available.
- **Fix grid independence** or keep scoping claims to what survives both meshes.
  Currently a stated failure (57.1% worst disc-band difference).
