# Polydispersity, particle size and the impact regime — where this stands

**2026-09-08, revised 2026-09-09.** One page, for a co-author or a fresh session.
🔴 **Sections 2(b) and 2(e) were revised and 2(e) withdrawn on 2026-09-09** — read
`NEXT_STEPS_2026-09-10.md` §0 and `PREREG_ROOT_PINNING_2026-09-09.md` §6–7 first.
Everything below is either
measured in this repo or derived and checked in it; the numbers are quoted from the
artifacts named beside them and are not typed from memory.

Full documents: `MECHANISM_2026-09-08.md` (the mechanism and the experiments),
`derivations/dilatancy_sign_mechanism.md` (the derivation + 21 checks),
`SIZE_THRESHOLD_2026-09-08.md` (the size hypothesis, v1 withdrawn),
`FIGURES_MECHANISM.md` (M1–M4).

---

## 1. The question

At similar particle size, raising the width of the size distribution moves an impacting
dense-suspension drop from solid-like/dilatant to spreading and inertia-dominated. Why?

## 2. What is now established

**(a) The size distribution enters the governing number in TWO places, and the project
had swept only one.** Deriving Mobaseri *et al.*'s inertial-dilatant number from Reynolds
dilatancy + Darcy + Terzaghi rather than quoting it gives

```
   ID = ρ κ(d₃₂(CV)) U₀ / (A μ_w (φ − φ_c(CV)) D₀)
```

The permeability channel is `(1+CV²)⁴`, whose log-derivative **vanishes** as CV → 0; the
critical-state channel is `1/(φ − φ_c(CV))`, whose log-derivative does not and which has
a **pole**. The 2026-08-31 exclusion of the permeability mechanism held `φ_c` fixed and
swept only `κ` — it measured the weaker channel. *(21 checks, `ALL CHECKS PASSED`.)*

**(b) `φ_c` rises with CV — but as a ladder-relative root, not a measured constant.**
🔴 **REVISED 2026-09-09**; read `PREREG_ROOT_PINNING_2026-09-09.md` §6–7. 63 cells:

| gate | Δ`φ_c` (CV 0 → 0.5) | slope | sd |
|---|---|---|---|
| ungated | **+0.0607** | 0.1214 | 42.6 |
| quasi-static | **+0.0446** | 0.0892 | 6.0 |
| steady | **+0.0471** | 0.0941 | 4.6 |
| both | **+0.0479** | 0.0957 | 6.8 |

The gates were **counted and never applied** until 2026-09-09; all four are now printed
side by side and none alone. **Quote +0.045 to +0.061 with the gate named.**

🔴 **And the "vs the isotropic slope 0.1246" comparison is WITHDRAWN as corroboration.**
The estimator's root follows the bottom of its own fitting window (response 1.15 when the
bottom moves, 0.06 when the top moves; on a ladder that brackets the root it refuses on
3/3 seeds with its exponent at the bound), and every ladder is CENTRED on those very
`φ_J` values — so the two agreeing to 2.6 % is one measurement read twice. What survives
is the sign and rough size, and a fit-free bracket: at CV = 0 the critical state lies
**between 0.5545 and 0.5945** (`p_ss` 7.65e-07 → 4.73e-04, `z_nr` 3.30 → 4.58).

⚠ A separate confound, found and controlled the same day: the realised shear rate is
not constant along the CV ladder (4.78×, because `dt` is set by the smallest particle). A
rate-matched control moves Δ`φ_c` by at most +0.0107, and in the **conservative**
direction. `grade_rate_matched.py`.

**(c) The headline needs no fit at all.** Where the CV = 0 and CV = 0.5 ladders overlap
in `φ`, at *matched packing fraction*:

| `φ` | CV | `p̂_ss` | `Z_nr` | rattlers |
|---|---|---|---|---|
| 0.63453 | 0.00 | 1.08e-2 | 6.12 | 0.004 |
| 0.63687 | 0.50 | **2.51e-4** | **4.43** | **0.448** |

**Raising CV collapses the undrained grain pressure by 43×** (63 sd). The second
overlapping pair gives 5.4×; **quote the range 5.4–43×**, because the two pairs tie on
`φ`-mismatch and the ratio grows as `φ_c` is approached. In both, the higher-CV cell is
the *denser*, so the ratios are conservative. Robustness on a late strain window: 41.0×
and 5.6×, a 4–5 % change.

**(d) The protocol has an external anchor.** Frictionless control: `φ_c(μ_p=0)` =
**0.6347 ± 0.0010** against `φ_J(μ_p=0)` = **0.6332**, measured months earlier by an
isotropic pressure ladder. Off by **0.0015**, seven times inside the pre-registered bar.

**(e)** 🔴 **WITHDRAWN 2026-09-09.** "`φ_c` and `φ_J` coincide to within 0.005 at
`μ_p` = 0.3" is forced by the ladder placement, not measured — see (b). It stands only at
`μ_p` = 0, where the estimator is E4 against the exact `Z` = 2d = 6 target and the
crossing is a bracketed interpolation, i.e. exactly the arm (d) tests.

**(f) `μ_s` = 0.359 ± 0.010 at `μ_p` = 0.3 and 0.105 ± 0.003 at `μ_p` = 0**, a factor
3.42: removing particle friction cuts macroscopic friction by that much and does not take
it to zero.

## 3. 🔴 What is NOT established, in order of how much it matters

**(i) The mechanism is degenerate with a simpler rival.** *Plain suspension rheology* —
CV raises `φ_m ≈ φ_c`, so `η ∝ (φ_m − φ)⁻²` falls and the drop spreads — consumes the
**same measured number**, predicts the **same direction**, and is consumed identically by
every DEM prediction here. **No simulation can separate them.** Three separators, all
experimental: the pore-pressure **sign**; the **strain history** (undrained locking
predicts a body that stiffens as it spreads, a viscosity has no strain history); and a
**boundary in (U₀, D₀)**, which a viscosity does not have at all.

**(ii) The lab's silica is on the WRONG side of the measured `φ_c`.** `φ` = 0.44 against
`φ_c` = 0.5695, i.e. 0.13 below — permanently contractive, never solid-like. Friction
cannot fix it: `φ_J` saturates at 0.5574, so **the gap friction cannot close is 0.117**.
Three possible resolutions, all live: a 23 % local compaction at the contact; the
material not being a frictional hard-sphere pack (the double layer *can* close the gap —
`φ_c^eff` spans 0.234–0.560 over plausible salt and stress, with 0.44 inside); or the
experiment's `φ` axis being low, which the Paper 1 card already flags.
⚠ **Every CV claim in §2 is a difference at fixed protocol and survives this**; every
statement about *which side* the lab is on does not.

**(iii) The drainage length is unmeasured and it moves `ID` by 10⁵.** `ℓ = D₀` gives
`Π` = 1.4e4 (undrained); `ℓ = 10 d₃₂` gives 0.14 (drained). `Π` = 1 sits at
**ℓ\* = 21.2 µm ≈ 26 grain diameters**, inside shear-band range — the system is *on* the
boundary.

**(iv) The dilatancy coefficient `α` is assumed constant in CV and is untested.** Testing
it needs constant-*pressure* shear, which this protocol does not do.

**(v) P5 RAN, and it does not discriminate.** 🔴 2026-09-09, `p5_degeneracy.py`: the
zero-skew family gives 2.823× less Δ`φ_c` than the lognormal one at matched `d₃₂` — but
the same shape effect is present in `φ_J`, measured with **no shear at all** (2.735×,
3.2 % apart). The effect is a property of the packing geometry, which both accounts read.
Original note: the zero-skew family at matched CV, which separates a CV effect
from a `d₃₂` effect.

## 4. The experiments, cheapest first

| | what | why it is decisive | cost |
|---|---|---|---|
| **E-0** | raw size distributions of the two bench samples | decides whether a pure size threshold can work at all (it needs CV ≈ 1.0–1.5); and whether "800 nm" is a number mean — if it is a volume mean, `d₃₂` *falls* with CV and the sign reverses | samples already exist |
| **E-6** | **drop-size ladder on cornstarch**, 2–6 mm | the three drainage readings move `U₀*` in three *directions* (∝D₀, ∝1/D₀, flat): 3.6 / 14.4 / 7.2 m/s at 2 mm. Decided by the **sign** | a different nozzle |
| **E-2** | raise `φ` at high CV until the solid branch returns | the increment should equal the measured Δ`φ_c` = **0.045 to 0.061** (the range is the four gate sets) — **no free parameter** | needs `φ` steps ≲0.03 and drift <0.015 — **not** blocked by evaporation (see below) |
| **E-3** | pore-pressure **sign** (flush transducer, >1 kHz) | the only *clean* separator from the rival | hardware |
| **E-4** | particle-size ladder 0.8/2.5/8 µm | three exponents: `a²` (100×), `a¹` (10×), `a⁰` (1×) over a decade | new particles |

⚠ **Withdrawn:** "look for cavitation in existing footage". The correct ceiling is 26.8 Pa,
four decades below atmospheric — nothing cavitates.

🔴 **And a correction made after `Δφ_c` was measured**: E-2 was written up as blocked by
the evaporation problem, on a tolerance of ±0.0025 in `φ`. That tolerance belongs to
*locating* CV\*, where `dφ_c/dCV` ≈ 0.12 converts a CV resolution into a much tighter `φ`
one. E-2 measures a **shift of 0.062**, needing an order of magnitude less control —
`φ` steps ≲0.03, drift <0.015. **E-2 is not blocked; locating CV\* is.**

### What the EXISTING data can already test — checked, `experiment_coverage.py`

Read-only from the per-drop table (**control variables only**; the project records that
the `arrested` flag once went stale, so no derived flag is used here):

| axis | range | levels | verdict |
|---|---|---|---|
| `φ` | 0.420 – 0.450, span **0.030** | 3 | 🔴 **2.1× too small for E-2** (needs 0.062) — and 3 levels cannot locate a threshold anyway |
| `d₀` | 2.10 – 2.61 mm, ratio **1.24** | ~continuous | 🔴 **nozzle scatter, not a designed axis** — E-6 needs deliberate nozzles |
| `U₀` | 1.15 – 5.13 m/s, ratio 4.5 | 104 | a real axis |
| **salt** | **0 / 10 / 50 / 100 mM** | **4** | 🔑 **the one designed ladder the lab already has** |

🔑 **So the cheapest test may be a re-analysis, not a new experiment.** The salt series
already exists, and the double-layer estimate has `φ_c^eff` moving 0.23 → 0.56 across
it — i.e. the *sign* of `φ − φ_c` should flip inside data already taken. ⚠ Doing that
needs a regime classification, and establishing which criterion is current must come
first.

## 5. What would kill it

`φ_c` flat in CV (it is not: 4.6–42.6 sd depending on the gate); the high-CV spreading state showing pore *suction*;
raising `φ` failing to restore the solid branch; a matched-`d₃₂` ladder still switching
regime; or no threshold anywhere over a decade of particle size.

## 6. Corrections made today, because they bear on how much to trust the rest

Six defects were found by checks that could fail, three of them in claims this same
session had already written down:

1. The **primary `φ_c` estimator was extrapolating the wrong way** — `γ̇` is fixed along a
   ladder, so `I → 0` is the *densest* rung. Retracted; the "0.074 systematic between
   estimators" was the two ends of the ladder, not a systematic.
2. The **suction ceiling was a withdrawn one** (cavitation, 101 kPa) where the repo's own
   replacement gives **26.8 Pa** — a 1.7e6 overshoot, not 459×.
3. The **transient model dropped a term**, reversing a "strain clock" into undrained
   locking. No test caught it because the test re-integrated the same equation.
4. The **size-threshold page's headline was an operating-point artefact** — one `a*`
   computed at silica's parameters and used to sort cornstarch. Withdrawn and rebuilt.
5. Compression **quantized by its own stopping rule**; two ladder rungs landed on the
   same state.
6. `sync_back.sh` **clobbered fifteen committed analysis outputs** with stale copies.
