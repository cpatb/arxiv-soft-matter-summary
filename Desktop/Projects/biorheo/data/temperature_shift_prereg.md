# Item 21 pre-registration — a temperature axis, without inventing one

**Committed before any code.** Bars are not moved. A bar that turns out to be the wrong
bar is recorded with its number and a *stricter* test is added beside it. A
pre-registration is never rewritten; corrections go in an addendum at the foot.

## 0. Why, and what is deliberately NOT being claimed

`../../biofluids_venture/OFFER.md` §3 lists **"No temperature dependence"** as the single
largest missing capability for the first product, and §6 ranks it second among the things
that block revenue. A cold-chain question — *"what does this formulation cost me at 5 °C
instead of 25 °C?"* — cannot be answered by anything in this package today. Report 01 of
the market research records a vendor figure of 43 cP at 5 °C against 24 cP at 15 °C for a
biologic, so the effect is roughly a factor of two over the range a fridge spans.

**No record in this registry carries a temperature axis, and this item does not add one.**
Every record is measured at a single temperature. So the item cannot and will not fit a
temperature dependence, calibrate an activation energy, or produce a number at a
temperature nobody measured. What it can do is carry a shift the CALLER supplies, exactly,
and refuse to invent one — the same stance the package already takes on needle gauges.

## 1. The model, stated before it is coded

Take the fluid to be **thermorheologically simple over the range in question**: changing
the temperature rescales the viscosity and the rate axis together and does not change the
shape of the flow curve. Then, with a shift factor `a` (dimensionless, `a > 0`),

```
eta(gd; T)  =  a * eta_ref(a * gd)
```

`a > 1` is a colder, more viscous fluid. This is the standard time–temperature
superposition of a viscosity curve, and its only content is that ONE number describes the
whole curve's move.

**This is an assumption, not a derivation, and it is falsifiable by the client's own
data**: two measured curves at two temperatures either do collapse onto one another under
a single `a` or they do not. The module will provide the collapse residual so that the
assumption is checked rather than asserted.

## 2. The parameter transformations, derived here before any code

For each family, substituting `a * gd` and multiplying by `a` gives a map on the fitted
parameters. These are algebra, done here, and G21.1 checks them numerically.

| family | `eta_ref(gd)` | transformation under `eta -> a eta(a gd)` |
|---|---|---|
| `newtonian` | `mu` | `mu -> a mu` |
| `power_law` | `k gd^(n-1)` | `a k (a gd)^(n-1) = a^n k gd^(n-1)`, so **`k -> a^n k`**, `n` fixed |
| `herschel_bulkley` | `tau_y/gd + k gd^(n-1)` | `a[tau_y/(a gd) + k (a gd)^(n-1)] = tau_y/gd + a^n k gd^(n-1)`, so **`tau_y` UNCHANGED**, `k -> a^n k`, `n` fixed |
| `casson` | `(sqrt(tau_y/gd) + sqrt(eta_c))^2` | `a(sqrt(tau_y/(a gd)) + sqrt(eta_c))^2 = (sqrt(tau_y/gd) + sqrt(a eta_c))^2`, so **`tau_y` UNCHANGED**, `eta_c -> a eta_c` |
| `cross` | `eta_inf + (eta_0-eta_inf)/(1+(lam gd)^m)` | `eta_0 -> a eta_0`, `eta_inf -> a eta_inf`, **`lam -> a lam`**, `m` fixed |
| `carreau_yasuda` | `eta_inf + (eta_0-eta_inf)(1+(lam gd)^a_cy)^((n-1)/a_cy)` | `eta_0 -> a eta_0`, `eta_inf -> a eta_inf`, **`lam -> a lam`**, `a_cy`, `n` fixed |

**The yield stress comes out INVARIANT in both yield-stress families, and that was not
imposed.** In the Casson case the `sqrt(a)` outside cancels the `sqrt(1/a)` inside
(`sqrt(a) sqrt(tau_y/(a gd)) = sqrt(tau_y/gd)`), and in Herschel–Bulkley the `a` cancels
the `1/a` in `tau_y/(a gd)`. That is the physically right answer — a yield stress is a
stress scale set by structure, not a viscosity — but here it is a *consequence* of the
superposition assumption rather than a choice, and it is worth stating because the opposite
would have been a reason to distrust the whole construction.

`structural_steady` is REFUSED. Its steady curve could be shifted, but the structure
equation carries a relaxation time that this one-parameter shift does not describe, which
is the same reason the OpenFOAM export refuses it (item 11).

## 3. The envelope moves, and that is the trap

A record's envelope `[gd_min, gd_max]` is where its curve was MEASURED. The shifted curve
`a eta_ref(a gd)` needs `a gd` inside that envelope, so it is valid only on

```
[gd_min / a, gd_max / a]
```

**A shift that did not move the envelope would silently widen the record's claim**, which
is the failure mode this package already has a rule about (item 3: "the record's refusal to
extrapolate does not survive the export"). The shifted record carries the shifted envelope,
and G21.5 checks it.

## 4. Where `a` comes from — and it is never from here

Two routes, both requiring the caller to supply a measurement or a constant:

- **`shift_from_two_viscosities(eta_cold, eta_hot)`** — the ratio of two viscosities
  measured at the same shear rate. This is what a client can actually send.
- **`arrhenius_shift(E_a, T, T_ref)` = `exp(E_a/R (1/T - 1/T_ref))`** — requires the
  caller's own activation energy. `R` is taken from `scipy.constants`, not typed from
  memory. **No default `E_a` is shipped for any fluid.**

There is no third route. The module will not guess a shift factor from a fluid name, a
temperature, or a table.

## 5. Gates

| gate | statement | bar | prediction |
|---|---|---|---|
| **G21.1** | the parameter map reproduces the shifted curve | `max |a*eta_ref(a*gd) / eta_shift(gd) - 1| < 1e-12` over 200 log-spaced rates, for **all six families**, at `a` in {0.25, 0.5, 2, 4} | PASS |
| **G21.2** | `a = 1` is the identity | every transformed parameter **bit-identical** to the original, all six families | PASS |
| **G21.3** | the shift composes | `shift(shift(p, a), b)` equals `shift(p, a*b)` to `1e-12`, all six families | PASS |
| **G21.4** | the discriminating half | at `a = 2` the viscosity at the envelope midpoint must move by more than 10 %, all six families — a map that returned its input would pass G21.1–G21.3 | PASS |
| **G21.5** | the envelope moves with the curve | shifted envelope is exactly `[gd_min/a, gd_max/a]`, and evaluating the shifted record outside it raises | PASS |
| **G21.6** | the yield stress is invariant | `tau_y` bit-identical for `casson` and `herschel_bulkley` at every `a`; and, as the DISCRIMINATING half, `k` and `eta_c` are NOT | PASS |
| **G21.7** | end-to-end closed form | for a Newtonian fluid the plunger force scales exactly as `a` (`F = 8 mu L Q A_s / (pi R^4)`), checked against the force solver at `a` in {0.5, 2, 4} to `1e-10` | PASS |
| **G21.8** | `structural_steady` is refused by a named error | raises, and the message says the STRUCTURE EQUATION is what cannot be carried — not "no counterpart exists" | PASS |
| **G21.9** | the collapse residual is reported, not assumed | given two curves, the routine returns the best `a` AND the residual at it; on a synthetic pair generated BY the shift the residual is < 1e-10, and on a pair that differs in SHAPE it is not | PASS both halves |

### 5.1 G21.9 is the one that can embarrass the method

Time–temperature superposition is an assumption. G21.9's second half exists so that the
package cannot report a shift factor without also reporting whether the two curves it came
from actually collapse. A routine that always returned a small residual would pass the
first half and be worthless; the shape-mismatch case is what makes the first half mean
something.

## 6. What this item does NOT do

- It does not add a record, fit anything, or touch `dst_fit`.
- It does not put a number on any real fluid at any real temperature. Every worked example
  will use a caller-supplied `a` and say so.
- It does not claim thermorheological simplicity for blood, for Carbopol, or for any
  protein formulation. That is the caller's assumption to make and G21.9's residual is how
  they check it.
- It does not touch the OpenFOAM export. A shifted record exports like any other record
  because it *is* a record; nothing new is needed there and nothing is claimed about it.

## 7. Artifacts

- `src/biorheo/constitutive/shift.py` — `shift_factor` helpers, `shift_parameters`,
  `shift_record`, `collapse_shift`, `ThermalShiftNotRepresentable`.
- `scripts/temperature_shift.py` -> `data/temperature_shift.md` / `.json`.
- `tests/test_shift.py`.
