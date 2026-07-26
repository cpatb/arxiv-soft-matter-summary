# Generalizing the Sean Earl result: what makes this more than a field report

**Status:** analysis memo, 2026-07-26. Numbers here trace to
`P3_CFD/vertical_structure.json`, `P3_CFD/vertical_fields.npz` and the paired
release-height sweep in `P4_Coupling/height_sweep.json`.

---

## 1. The finding that generalizes

> **Shelter is not protection.** Buildings near a playing field remove *mean*
> wind and add *turbulence*. For a thrown object, mean wind is a bias the
> thrower can aim off; turbulence is variance they cannot. A "sheltered" pitch
> can therefore be systematically *harder* to throw on than an exposed one,
> while feeling calmer to a person standing on it.

Measured at 2 m over the playing surface, against the identical inflow ABL run
over empty fetch (`--control abl_control.npz`, the no-buildings
counterfactual):

| θ_wind | mean speed vs control | σ_u vs control | TI pitch | TI control |
|---|---|---|---|---|
| 90° easterly | **0.44×** | 0.91× | 0.25 | 0.12 |
| 210° SSW | **0.50×** | **1.51×** | 0.36 | 0.12 |
| 270° westerly | **0.50×** | **1.57×** | 0.38 | 0.12 |

The buildings cut mean wind roughly in half and, in two of three directions,
raised the absolute velocity fluctuation by half again. Turbulence intensity
went from 12% to 25–38% — a factor of 2–3.

This is the counterintuitive part and it is the paper's most transferable
claim. It is also *not* what urban wind engineering optimises for: the COST 732
/ Lawson / NEN 8100 tradition is built around **pedestrian** comfort and safety,
which is a mean-and-gust-exceedance criterion for a person standing still.
Nobody has written down the equivalent criterion for a *projectile*.

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

### 3. Does gyroscopic stiffness buy gust immunity?
A disc is a gyroscope; a gust's rolling moment precesses into a *heading* change
90° later. Higher spin ⇒ stiffer ⇒ smaller deflection per gust. There should be
a critical spin rate above which turbulence effects are materially suppressed.
**Probe:** cheap — sweep release spin in the existing 6-DOF model at fixed
release velocity. If it holds, it explains *why* coaches teach "more spin" in
wind, which is currently folklore. Novel, and entirely within reach today.

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
I expected vertical shear to be the disc-specific villain (a gyroscope in shear
should misbehave). The data says otherwise: the sheltered pitch has *lower*
near-ground shear (0.02–0.06 s⁻¹ over the disc band) than the open-fetch control
(up to 1.5 s⁻¹ near the ground, a log profile). Shelter flattens the profile.
The hazard is σ_u, not `∂U/∂z`. Reporting this stops the next group chasing it.

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
