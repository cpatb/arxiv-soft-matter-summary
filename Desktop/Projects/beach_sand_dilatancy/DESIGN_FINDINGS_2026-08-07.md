# Design findings, before any sand is collected

**2026-08-07.** Output of `dilatancy_design.py` (22 tests in `test_dilatancy_design.py`).
`PROTOCOL.md` §2 says the regime estimate is "worth doing before buying anything". This is
that estimate, made computable. Three results, one of which changes the rig specification
and one of which changes what the project should claim to be about.

---

## 0. A notation collision, fixed first

`PROTOCOL.md` uses `φ` for two different quantities: the **porosity** in the Kozeny-Carman
line (`k ~ 10⁻³ d², φ ≈ 0.4`) and the **packing (solid) fraction** everywhere else
(`φ₀ = 0.560 vs 0.604`). A natural sand pack has porosity ≈ 0.40 and solid fraction ≈ 0.60,
so reading one for the other changes `k` by 7.6× and moves every derived speed with it. It
cost exactly that error while writing this module.

`dilatancy_design.py` never uses the bare name `φ`; it takes `porosity` explicitly, and a
regression test asserts the two are not interchangeable. Recommend the same fix in
`PROTOCOL.md`: write `n` for porosity throughout.

With that settled, the module reproduces the §2 table: `k = 4.0 × 10⁻¹¹ m²` and `Δp ≈ 101 kPa`
for 0.2 mm sand, `9.9 × 10⁻¹⁰ m²` and `4 kPa` for 1.0 mm, at `V = 2 m/s`, `D = 20 mm`,
`ψ = 0.1`. (The table's `ρ_s V² = 8 kPa` implies `ρ_s ≈ 2000`; quartz at 2650 gives 10.6 kPa,
so the ratios here are 9.6 and 0.4 rather than 12 and 0.5. Same conclusion.)

---

## 1. Q1 is well posed, and comfortably so

The design question: how well-sorted may the sand be before the three candidate diameters
stop making distinguishable predictions for `a_max ∝ d^-0.8`?

| C_u | d₁₀/d₅₀ | d_Sauter/d₅₀ | a_max ratio using d₁₀ | using Sauter |
|---|---|---|---|---|
| 1.2 | 0.859 | 0.993 | 1.13× | 1.01× |
| 1.5 | 0.713 | 0.966 | 1.31× | 1.03× |
| **2.0** | 0.561 | 0.903 | **1.59×** | 1.08× |
| 3.0 | 0.400 | 0.774 | 2.08× | 1.23× |
| 4.0 | 0.314 | 0.665 | 2.52× | 1.39× |

At `C_u = 2`, a well-sorted beach sand, the `d₁₀` and `d₅₀` forms of `a_max` differ by 59%,
far above any plausible scatter. The minimum sorting needed to discriminate at 15% scatter
is `C_u ≈ 1.23`, which essentially every real sand exceeds.

The Sauter mean sits much closer to `d₅₀` than `d₁₀` does, so the informative contrast is
**Hazen (`d₁₀`) against Kozeny-Carman (Sauter)**, not either against the median. Design
implication: measure the full grading curve, not just `d₅₀` and `d₁₀`, because the Sauter
mean needs the whole distribution.

**Verdict: Q1 survives contact with the numbers.** It was already the best question in the
protocol and it remains so.

---

## 2. The rig as specified is in the wrong regime, and the fix is free

The hardening ratio is `Δp/(ρ_s V²) = η ψ D / (k ρ_s V)`. At fixed sand it depends on the
intruder **only through `D/V`**. So matching the field phenomenon means matching `D/V`, not
matching `V`.

- A footstep: `D ≈ 0.10 m`, `V ≈ 1 m/s`, so `D/V = 0.1 s`.
- The proposed rig: 20 mm sphere at 2 m/s, so `D/V = 0.01 s`.

**That is a factor of 10 off.** The proposed rig sits an order of magnitude further into the
drained regime than the barefoot observation that motivated the project.

Two ways to fix it, and one is nearly free:

| fix | requirement |
|---|---|
| bigger intruder at the same speed | `D = 0.20 m` at 2 m/s, a 200 mm sphere |
| **same intruder, slower** | **20 mm sphere at 0.20 m/s, a 2 mm drop height** |

The second costs nothing and is the obvious move: keep the sphere, drop it from millimetres
rather than centimetres, and use the high-speed camera to resolve the shorter, gentler
deceleration. Better still, sweep `D/V` deliberately and show the crossover, since the whole
mechanism predicts that hardening is controlled by that one group.

**This is the single most useful thing in this document**, because it is a change to the rig
specification that costs nothing now and would be expensive to discover after the sand is in
the lab.

---

## 3. Q4's ceiling is reachable, but the simple model cannot span Strader's range

**The reachable part.** Suction cannot exceed about one atmosphere, so `V^1.6` must saturate
at `V_cav = k(p_atm - p_vap)/(η ψ D)`, which scales as `d²`:

| d | V_cav (ψ = 0.1, D = 20 mm) | drop height |
|---|---|---|
| 0.10 mm | 0.49 m/s | 0.01 m |
| 0.15 mm | 1.10 m/s | 0.06 m |
| 0.20 mm | 1.96 m/s | 0.19 m |
| 0.30 mm | 4.40 m/s | 0.99 m |
| 0.50 mm | 12.2 m/s | 7.6 m |
| 1.00 mm | 48.9 m/s | 122 m |

So the ceiling is inside an ordinary bench drop for `d ≲ 0.3 mm` and out of reach above
`0.5 mm`. **That contrast is the test**: fine sand should show the departure, coarse sand
cannot, and the boundary should move as `d²`. Q4 is cheaper and sharper than the protocol
credits it for, and it needs a shelf rather than a tower.

**The part that does not work.** Strader et al. report an unbroken `a_max ∝ d^-0.8 η^0.4 V^1.6`
up to about 2.5 m/s across beads from 53 to 850 µm. Two requirements follow, and they
contradict each other:

- the fine end must not have cavitated, or the law would have broken there:
  `ψD ≤ k(53 µm)(p_atm - p_vap)/(η V) = 1.10 × 10⁻⁴ m`
- the coarse end must still have been pore-pressure controlled, or the same law could not be
  dilatancy physics there: `ψD ≥ ρ_s V k(850 µm)/η = 4.73 × 10⁻³ m`

**The window is empty by a factor of 43.** Since `k ∝ d²`, the window closes as
`(d_fine/d_coarse)² = 1/257`, so any wide grain-size sweep does this.

Taken alone, the cavitation constraint bounds the protocol's loosest number:
`ψ ≤ 0.0055` for `D = 20 mm`, **at least 18× below the assumed 0.1**.

Read this as a statement about the estimate, not about their measurements. A single constant
`ψ`, with the intruder diameter as the only length scale, cannot span a 16× range in grain
size. Candidate resolutions, all testable: `ψ` depends on the packing state rather than being
a constant; the relevant length is the dilating shear-band thickness rather than `D`; or
partial cavitation still transmits stress, so the ceiling is soft rather than hard.

**This reframes the project, and upward.** `PROTOCOL.md` calls Q4 "speculative, worth one
afternoon, not worth building the project on". On these numbers Q4 is where the existing
model demonstrably fails, and mapping *where* it fails is a better question than confirming
Hazen. The honest expectation in the protocol ("Q1 most likely confirms Hazen, a modest
paper") can be upgraded: a measured `ψ(φ₀, d)` that closes the 43× gap would be a real
contribution, and the beach hands you the `φ₀` spread for free.

---

## What to change

1. **`PROTOCOL.md` §2**: write `n` for porosity, keep `φ` for packing fraction only.
2. **Rig spec**: sweep `D/V`, and include a 0.2 m/s drop so the footstep regime is actually
   reached. Do not buy a 200 mm sphere; drop the 20 mm one from 2 mm.
3. **Promote Q4** from an afternoon to a named aim, framed as measuring `ψ` rather than as
   confirming a saturation.
4. **Grading curves, not percentiles**: the Sauter mean needs the full distribution, and it
   is the discriminating quantity against Hazen's `d₁₀`.
5. Everything above rests on `ψ`, `D` and the Kozeny-Carman constant. Sweep them rather than
   quoting single numbers; every function here takes them as arguments for that reason.
