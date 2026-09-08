# Item 21 — a temperature axis, without inventing one

Written by `scripts/temperature_shift.py` from `data/temperature_shift_prereg.md`. biorheo 0.0.1. Every number on this page is printed by that script.

## What this is, and what it is not

**No record in this registry carries a temperature axis, and this item does not add one.** Every record is measured at a single temperature. Nothing here fits a temperature dependence, calibrates an activation energy, or produces a number at a temperature nobody measured. What it does is carry a shift the CALLER supplies, exactly, and refuse to guess — the same stance this package already takes on needle gauges.

The assumption, which is falsifiable by the client's own data:

```
eta(gd; T)  =  a * eta_ref(a * gd)         a > 1 is the colder, more viscous fluid
```

Its whole content is that ONE number describes the entire curve's move. Two measured curves either collapse under a single `a` or they do not, and `collapse_shift` returns the residual so that is checked rather than asserted (G21.9).

## Gates

| gate | statement | result |
|---|---|---|
| G21.1 | the derived parameter map reproduces `a*eta(a*gd)` on all six families | **PASS**, 24 cases, worst relative error 6.66e-16 |
| G21.2 | `a = 1` is the identity, bitwise | **PASS**, 6/6 |
| G21.3 | the shift composes | **PASS**, worst 2.22e-16 |
| G21.4 | the discriminating half — the map is not the identity | **PASS**, smallest move at `a = 2` is 14.3 % |
| G21.5 | the envelope moves and the sidecar is dropped | **PASS** |
| G21.6 | `tau_y` invariant, viscous parameter not | **PASS** |
| G21.7 | Newtonian force scales exactly as `a` | **PASS**, worst 1.84e-12 |
| G21.8 | `structural_steady` refused, for the RIGHT reason | **PASS** |
| G21.9 | the collapse residual discriminates | **FAILED**, 1.442e-05 against a 1e-10 bar |
| **G21.9′** | that residual CONVERGES at the interpolation order | **PASS**, observed order 1.832 over an 8x refinement |

**G21.9 FAILED and its bar has not been moved.** It required the clean residual to fall below 1e-10; it is 1.442e-05. The bar assumed exact recovery from a routine that compares against a LINEAR INTERPOLATION of the reference curve in `(log gd, log eta)`, so the residual is bounded below by the interpolation error of whatever grid it is handed — **no correct implementation can beat it**, which makes this a mis-specified bar rather than a failed one. Its number stays on this page.

**G21.9′ is the stricter test beside it**, and it is stricter in kind rather than looser in size: if the residual really is interpolation error and nothing else, it must CONVERGE at the interpolation order under grid refinement. A wrong parameter map, a broken objective or a genuine shape mismatch all leave a residual that does not fall with `h` at all, so this cannot be passed by being merely small.

| rates | grid spacing (decades) | rms log residual | recovered `a` relative error |
|---|---|---|---|
| 200 | 0.02513 | 1.4417e-05 | 4.16e-06 |
| 400 | 0.01253 | 4.8729e-06 | 1.43e-06 |
| 800 | 0.00626 | 1.5899e-06 | 4.70e-07 |
| 1600 | 0.00313 | 3.1681e-07 | 9.34e-08 |

Observed order **1.832**, against 2 for linear interpolation. The discriminating half of G21.9 — that a shape-mismatched pair does NOT give a small residual — passed on its own and is reported below.

## G21.6 — the yield stress falls out invariant, and that was derived first

The pre-registration's §2 works the algebra for all six families **before any code existed**, and records the non-obvious consequence: under this shift the yield stress is UNCHANGED in both yield-stress families. For Casson the `sqrt(a)` outside cancels the `sqrt(1/a)` inside; for Herschel–Bulkley the `a` cancels the `1/a` in `tau_y/(a gd)`.

That is the physically right answer — a yield stress is a stress scale set by structure, not a viscosity — but here it is a **consequence** of the superposition assumption rather than a choice, and the opposite result would have been a reason to distrust the whole construction. The discriminating half is in the same gate: the viscous parameter beside it must move, or the map would be reporting invariance by doing nothing.

| family | a | `tau_y` invariant | viscous parameter | its ratio |
|---|---|---|---|---|
| `casson` | 0.25 | yes | `eta_c` | 0.25 |
| `casson` | 0.5 | yes | `eta_c` | 0.5 |
| `casson` | 2 | yes | `eta_c` | 2 |
| `casson` | 4 | yes | `eta_c` | 4 |
| `herschel_bulkley` | 0.25 | yes | `k` | 0.53396 |
| `herschel_bulkley` | 0.5 | yes | `k` | 0.73072 |
| `herschel_bulkley` | 2 | yes | `k` | 1.3685 |
| `herschel_bulkley` | 4 | yes | `k` | 1.8728 |

## G21.9 — the gate that can embarrass the method

Time–temperature superposition is an **assumption**. A routine that reported a shift factor without a residual would give a confident number for two curves that do not collapse at all, because a best-fitting `a` always exists. **This gate's first half FAILED on a mis-specified bar (see above); its second half — the half that makes the first mean anything — passed on its own.**

- On a synthetic pair generated BY the shift (`a = 2.5`): recovered **a = 2.4999896** (relative error 4.2e-06), rms log residual **1.44e-05**, on 184 of 200 rates.
- On two curves matched in magnitude at mid-range but differing in SHAPE: rms log residual **0.943**, 6.5e+04x the clean case. The returned `a = 0.4714` is meaningless rather than approximate, and the returned dict says so in words.

`n_used` is reported because a factor fitted on three surviving rates is not a collapse either: shifting the rate axis takes part of the reference curve outside its own measured range, and scoring the fit there would be scoring an extrapolation.

## Where `a` comes from

Two routes, both requiring the caller to supply a measurement or a constant, and there is no third:

- `shift_from_two_viscosities(eta_cold, eta_hot)` — the ratio of two viscosities measured at the same shear rate. This is what a client can actually send.
- `arrhenius_shift(E_a, T, T_ref)` — with `R` from `scipy.constants`, not typed from memory. **No default activation energy is shipped for any fluid.** As an illustration of the arithmetic only: `E_a = 20000` J/mol from 298.15 K to 278.15 K gives `a = 1.7862`. **That is not a property of anything in this registry.**

## What a shifted record is

An ordinary `FluidRecord` — it exports, evaluates and reports like any other — with three differences, all stated on it:

1. the parameters are mapped;
2. the envelope is divided by `a` (checked exactly at 4 factors, and evaluating outside it raises); a shift that left the envelope alone would silently widen the record's claim;
3. **the identifiability sidecar is dropped, not mapped.** A profile-likelihood interval belongs to a fit against measured data, and there is no data at the shifted temperature. Carrying the old interval across would be a confidence statement about measurements nobody made. A shifted record therefore reports as a POINT ESTIMATE in `biorheo.apps.report`, which is the honest description of it.
