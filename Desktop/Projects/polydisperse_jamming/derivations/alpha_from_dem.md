# The dilatancy coefficient — RETRACTED

**Written 2026-09-12; retracted 2026-09-12, the same day.** Every number in this page is
printed by `derivations/alpha_from_dem_numbers.py` into `alpha_from_dem_numbers.txt`
(`ALL CHECKS PASSED -- the retraction's evidence holds`). The original claim is kept in
§1 verbatim, beside the reason it is wrong.

---

## Result box

> **Withdrawn:** `α` (dilating) = 6.62 ± 1.25, `α` (contracting) = 3.31 ± 0.65, the
> 2.00× asymmetry, the 3.3–6.6× factor on `Π` and `1/ID`, the conditional 12–24× closure of
> M1's last corner, and the E-3 "magnitude" upgrade. **E-3 is a sign test again.**
>
> **What was actually measured:** on the dense branch, `SERVO_EDOT_CLAMP` = 0.5. On the
> loose branch, a window that is 60–100 % servo compression of an unjammed pack.
>
> **What stands:** `poro_impact.ALPHA_DIL` = 1.0, an *import* from [M26], unmeasured. The
> algebra (`α` dimensionless, `Π` linear in `α`, `Λ` free of `α`) is correct and now moot.

---

## 1. The claim, verbatim

```
    alpha (dilating,    phi > phi_c)  =  6.62 +- 1.25     (n = 6)   <- fitted on phi0 = 0.55
    alpha (contracting, phi < phi_c)  =  3.31 +- 0.65     (n = 6)   <- fitted on phi0 = 0.64
    ratio 2.00x
```

with `α` read as the slope of `ln|φ − φ_ss|` against `γ` over the first 40 % of strain
(`ALPHA_FIT_FRAC` = 0.4) on the 18 production cells of `PREREG_CONSTP_PROD_2026-09-12.md`.

## 2. Three defects

**D1 — the labels were swapped.** `φ_c` = 0.591 (CV 0) and 0.602 (CV 0.2). A pack started
at `φ₀` = 0.55 sits *below* `φ_c` and must **compact**: it is the contracting branch. The
script called it "dilating". Check L1: 12 of 12 `φ₀` = 0.55 cells have `φ₀ < φ_ss`; 9 of 9
`φ₀` = 0.64 cells at CV ≤ 0.2 have `φ₀ > φ_ss`. Backwards on every cell, not on some.

**D2 — the estimator reads the servo, not the grains.** The constant-pressure servo is
`ε̇_v = K γ̇ (p/p* − 1)`, `K` = 0.15 (from the records), clamped at `|ε̇_v| ≤ 0.5 γ̇`.

- *Dense side.* A pack compressed to `φ₀` = 0.64 starts at `p/p*` = 45–59, so the servo
  is on its clamp from the first step and `d ln φ/dγ = −0.5`. The strain that takes
  `φ₀` to within 10 % of the gap is then `ln(φ₀/φ₁₀)/0.5` = 0.10–0.15 by that arithmetic
  alone. Measured transit: 0.120–0.165, i.e. **1.12–1.22× the clamp-only prediction**,
  with the servo on its clamp for **86–90 %** of the transit on every cell (check L2).
  The "dilation" is the box expanding at its maximum permitted rate to relax an
  over-compressed pack to `p*`. `α`(dense) is a readout of `SERVO_EDOT_CLAMP`.
- *Loose side.* A pack at `φ₀` = 0.55 starts **unjammed** (`p` = 0 exactly, `p₀/p*` = 0.0
  on all 12 cells). The servo compresses it at `K γ̇` against no resistance until it
  jams; **60–100 %** (median 89 %) of the transit is spent at `p < p*/2`. That part is
  the controller, not dilatancy.
- *The fit window.* 40 % of `γ_max` = 3 is a strain of 1.2. On the dense side the transit
  is over by 0.165 and the remaining ~1.0 of the window is seed noise about `φ_ss`; on the
  loose side the window straddles the servo compression and whatever material approach
  follows it. Neither slope is `α`. Check L3 makes this concrete: the fitted `α` says the
  contracting pack approaches **1.63× faster**; the fit-free half-strain (`γ_half` =
  0.093 ± 0.010 dilating vs 0.370 ± 0.127 contracting) says the dilating pack is
  **3.96× faster**. The estimator and the raw series disagree on the *ordering*.

**D3 — at CV 0.5 there is no dilating branch.** `φ_c`(CV 0.5) = 0.641–0.652 sits *above*
the "dense" start of 0.64, so `|φ₀ − φ_ss|` = 0.0007–0.0116 (check L4) and the fitted
`α` there (0.40–1.44) is a fit to seed scatter.

## 3. What this withdraws downstream

| Where | What was written | Status |
|---|---|---|
| this page §2–4 (old) | `Π` and `1/ID` corrected by 3.3–6.6×; M1's corner closes 12–24×; E-3 becomes a 2:1 magnitude test | **withdrawn** |
| `poro_impact.py` | `ALPHA_DIL_MEASURED_DILATING/CONTRACTING/SD` | **removed**; comment records the retraction |
| `test_poro.py` | `test_measured_alpha_differs_from_the_default` | replaced by `test_no_measured_alpha_is_shipped` |
| `MECHANISM_2026-09-08.md` §6a / E-3 | magnitude asymmetry | to be struck; E-3 is a sign test |

Nothing in the constant-pressure campaign's **primary** results moves: `φ_c(CV)` is a
steady-state quantity read after the transient, and the transient's shape was never part
of its admission.

## 4. What measures `α` properly

`α` is the coefficient of the *grains'* volumetric response to shear strain. Reading it
needs two things the campaign did not have:

1. **The servo must be much faster than the material**, so `ε̇_v` is set by the grains.
   Here the clamp (0.5 γ̇) is *slower* than the transit it was supposed to observe.
2. **The pack must start at the pressure it is sheared at**, so the transient is the
   Reynolds one and not a pressure relaxation. Consolidate isotropically to `p*`, then
   shear at `p*`:
   - *contracting side*: a frictional pack, `φ_iso(p*) < φ_c(p*)` (measured: the
     isotropic jamming line lies below the critical-state line at every CV — that is the
     constant-volume unloading result written up in `cv_unloading.md`);
   - *dilating side*: a **frictionless** consolidation to `φ_iso,μ=0(p*) > φ_c(p*)`, with
     friction switched on at the first shear step (`utils.setContactFriction` in Yade).
3. **Stop the fit where `|φ − φ_ss|` reaches the seed scatter**, not at a fixed fraction
   of the strain.

That protocol is the strain-clock campaign (`PREREG_STRAINCLOCK_2026-09-13.md`), which
measures `γ_c(φ, CV)` for shear jamming and, on the same runs, the Reynolds transient
`φ(γ)` from `φ_iso(p*)` — so `α` comes out of it as a by-product with the servo removed
from the question.

## 5. Fit-free numbers that survive

`γ_half` on the contracting side, from a common `φ₀` = 0.55, *including* the servo
compression (so an **upper bound** on the material strain to jam from that state):

```
    CV 0.0:  0.285 0.285 0.285 0.285 0.285 0.270
    CV 0.2:  0.330 0.330 0.345 0.345 0.345 0.345
    CV 0.5:  0.555 0.600 0.570
```

The CV 0.5 pack also has ~2× the gap to close (`φ_c` is 0.05 higher), so this is not a
rate statement about `α`. It is recorded because it is the one transient number the
campaign can quote without a model.

## 6. Lesson

A transport-like coefficient fitted through a window where a *controller* is saturated is
a readout of the controller's limit. The checks that catch it are cheap and fit-free:
print the fraction of the window at the clamp; compare the transit strain with the
clamp-only prediction; and check the branch *ordering* against a fit-free landmark. The
label swap would have been caught by printing `sign(φ₀ − φ_ss)` beside each `α`.
