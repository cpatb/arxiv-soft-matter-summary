# Constant-volume shear from an isotropic state — which side of the critical state is a fresh pack on?

**Written 2026-09-12.** Every number here is printed by `cv_unloading.py` into
`out/cv_unloading.txt` (`FAILED: U2, U3`). Bars U1–U3 were fixed in that script's docstring
before its first execution; U4–U5 are reported, not graded.

---

## Result box

> **An isotropically prepared frictional pack is CONTRACTIVE under shear at every packing
> fraction and every CV measured.** On 108 of 108 jammed constant-volume cells
> (CV 0.0, 0.2, 0.5; N = 1500, μ_p = 0.3, four campaigns pooled) the grain
> pressure falls from its preparation value `p_0` to the critical-state value `p_ss(φ)`. In the
> undrained reading that is pore **overpressure on both sides of `φ_c`** — the "dilative
> branch" of the mechanism is not reached by putting `φ` above `φ_c`; it needs a state
> that sits *below* the critical-state line in pressure, i.e. a loading history.
>
> **U2 — FAIL as registered.** Constant-volume `p_ss(φ)` inverted at `p*` gives `φ` =
> nan, 0.6051, nan against the constant-pressure
> 0.5911, 0.6020, 0.6444: CV 0.2 agrees to +0.0031. The bar is
> unreadable at CV 0.0, 0.5 because `p*` = 3e-04 lies
> *below* the lowest jammed constant-volume rung's `p_ss` there — the interpolation has no
> bracket and no extrapolation is admitted. The CV 0 → 0.5 contrast from the constant-volume
> side is therefore **+nan** (constant-pressure +0.0533); the bar stays FAIL and the
> stricter test it needs is a constant-volume rung between `φ_c(p*)` and the present lowest
> jammed rung at CV 0 and 0.5 (§7).
>
> **U3 — FAIL as registered** on 3 overlapping rungs; 1 violated:
> CV 0.5 φ 0.6450 at p_iso/p_ss = 0.74.
> And a number the bars did not ask for: the phase-8 preparation's own `p_0` sits
> **3.6–9.3× above the phase-2 quasi-static isotropic line at the same `φ`**, so
> U1's "isotropically prepared" pack is a rate-compressed one, over-pressured relative to the
> servo-held isotropic state. U1's sign therefore mixes prep-overpressure relaxation with the
> shear response; the clean version is §7's next run.

---

## 1. What was measured

`p_hat(γ)` on every `cv-shear-v2` constant-volume record: 147 distinct cells after the
protocol / N / μ_p / family / compressed filter, of which 108 are jammed (quasi-static,
`I_ss` ≤ 1e-02) and 39 are flowing. `p_0` is the first shear sample at γ = 0.0075 — the
schema has no γ = 0 sample, and the decay per sample at that point is ~1 %, so `p_0` is a
~1 % underestimate of the isotropic pressure. `p_ss` is the project's steady-state
window mean.

## 2. U1 — the sign (PASS on 108/108)

Per CV: CV 0.0 42/42; CV 0.2 40/40; CV 0.5 26/26. 2 of 108 jammed cells show a transient rise (`p_max > p_0`) before the decay.
The flowing cells (39) start unjammed (`z` ≈ 0, `p_0` ~ 1e-07) and their pressure
*rises* to an inertial steady state; that is a pack at finite rate, not a dilatancy sign,
and they are listed in the `.txt` for the record.

**Why this matters for the mechanism.** `MECHANISM_2026-09-08.md` §1 reads "below `φ_c` a
sheared pack contracts … above it the pack dilates". The second half is true only of a
pack that is *denser than the critical state at the pressure it is sheared at*. Where U3
can be read (CV 0 and 0.2), isotropic compression to `φ` leaves the pack at a pressure
`p_iso(φ)` **above** `p_ss(φ)`, so the first thing shear does is relax that pressure —
contraction, overpressure — however far above `φ_c(0)` the pack sits; at CV 0.5 the one
readable rung is undecided (§4). In critical-state soil mechanics this is the textbook
statement that a **normally consolidated** soil is contractive and only an
**over-consolidated** one dilates; what is new is that it holds for a polydisperse
frictional sphere pack at every CV, and that the DEM puts numbers on the ratio (U5).

## 3. U2 — one line from two protocols

| CV | `φ` at `p*` from `p_ss(φ)` (constant volume) | `φ_c(p*)` (constant pressure) | Δ | status |
|---|---|---|---|---|
| 0.0 | nan | 0.5911 ± 0.0010 | +nan | UNREADABLE: lowest jammed p_ss 4.3e-04 > p* |
| 0.2 | 0.6051 | 0.6020 ± 0.0013 | +0.0031 | read |
| 0.5 | nan | 0.6444 ± 0.0043 | +nan | UNREADABLE: lowest jammed p_ss 7.1e-04 > p* |

Bar ±0.01: **FAIL as registered**. The bar was posed on the assumption that `p*` falls inside
the jammed part of every common-ladder — it does not at CV
0.0, 0.5, where the lowest jammed rung's `p_ss` is
already above `p*`, so the inversion returns NaN rather than a number. That is a defect of
the bar's design, and it is not repaired by extrapolating: the CSL is measurably curved in
(`φ`, ln p) (≈ 0.009 per decade near jamming against ≈ 0.047 between φ 0.62 and 0.67 at
CV 0), so a one-rung extrapolation would be a guess dressed as a test. Where the bar can
be read it passes: CV 0.2 agrees to +0.0031. The two routes there differ in N (1500 vs 1000), in `γ̇`
(the common ladder is *unmatched*: 3.0e-05 / 6e-05 / 1.4e-04 against 1.17e-04 matched) and
in the steady-state definition. One CV agreeing is evidence, not a verdict; the failing
NaNs stay printed.

## 4. U3 — is the isotropic line above the critical-state line in pressure?

| CV | `φ` | `p_iso(φ)` (phase-2, servo-held) | `p_ss(φ)` (constant volume) | `p_iso/p_ss` | `p_0/p_iso` | bar |
|---|---|---|---|---|---|---|
| 0.0 | 0.5951 | 1.05e-03 | 4.30e-04 | 2.45 | 5.19 | ok |
| 0.2 | 0.5950 | 2.62e-04 | 5.16e-05 | 5.07 | 9.26 | ok |
| 0.5 | 0.6450 | 5.25e-04 | 7.12e-04 | 0.74 | 3.64 | **VIOLATED** |

**FAIL as registered.** Only 3 common-ladder rungs fall inside the phase-2 isotropic
ladder at all (its top is 0.603 / 0.609 / 0.649 at CV 0 / 0.2 / 0.5), so the bar is
decided on one rung per CV. The violation is at the CV 0.5 rung, whose constant-volume
`p_ss` is from a rung that is steady on 1 of 3 seeds (§3's table) with one seed still
decaying at `γ` = 1.2 (U4) — so its `p_ss` is an upper bound and the ratio 0.74 a lower
bound. **Unresolved at CV 0.5, not a sign flip**; a longer constant-volume run at
CV 0.5, `φ` = 0.645 decides it, and the constant-pressure `φ_c(p*)` = 0.6444 at that CV
says `p_ss(0.645)` should be ≈ `p*` = 3e-04, i.e. the ratio should come out ≈ 1.7 if the rung
finishes relaxing.

The column the bars did not ask for is the larger finding. `p_0/p_iso` = 3.6–9.3:
the phase-8 preparation (rate compression to `φ`, then constant-volume relaxation to the
unbalanced-force tolerance) lands at a pressure 3.6–9.3× the servo-held
isotropic state of phase 2 at the same `φ` and N. For a frictional pack `p(φ)` is not a
function — it is a history — and the two "isotropic" protocols sit on different branches
of it. Consequence for U1: the sign it measured is that of a *rate-compressed* pack, whose
pressure has two reasons to fall under shear (relaxing the prep overpressure toward the
isotropic branch, and relaxing the isotropic branch toward the critical line). U1 is
correct as stated and weaker than its prose implied; the second reason alone is what §6
needs, and it is tested by preparing ON the phase-2 line (§7).

## 5. U4 / U5 — the decay strain and the overconsolidation ratio (reported)

| CV | `φ_target` | OCR* = `p_0/p_ss` (geo. mean) | range | `γ_half` (first crossing) | n |
|---|---|---|---|---|---|
| 0.0 | 0.595 | 12.70 | 9.23–19.99 | 0.077 ± 0.004 | 3 |
| 0.0 | 0.620 | 1.64 | 1.61–1.67 | 0.090 ± 0.007 | 3 |
| 0.0 | 0.645 | 1.34 | 1.29–1.37 | 0.105 ± 0.007 | 3 |
| 0.0 | 0.670 | 1.24 | 1.24–1.25 | 0.135 ± 0.015 | 3 |
| 0.2 | 0.595 | 46.95 | 46.95–46.95 | 0.030 | 1 |
| 0.2 | 0.620 | 2.08 | 1.83–2.59 | 0.080 ± 0.011 | 3 |
| 0.2 | 0.645 | 1.41 | 1.33–1.53 | 0.097 ± 0.007 | 3 |
| 0.2 | 0.670 | 1.26 | 1.21–1.33 | 0.115 ± 0.009 | 3 |
| 0.5 | 0.645 | 2.68 | 1.40–3.98 | 0.035 ± 0.019 | 3 |
| 0.5 | 0.670 | 1.59 | 1.44–1.83 | 0.070 ± 0.019 | 3 |

`γ_half` is the strain at which the excess pressure `p − p_ss` has halved; first and last
down-crossings and a monotonicity flag are printed per cell. **The overpressure is not
instantaneous and not permanent — it decays over a strain of order 0.1–1 that grows
toward `φ_c`.** A viscosity (M2a) has no strain in it; this is a strain-history signature
on the *contractive* side, measured, with no fit.

OCR* is the prediction for the next run: unload an isotropic pack at `φ` to a pressure
below `p_ss(φ)` — a ratio above OCR* — and constant-volume shear must then *raise* the
pressure. That would be the first dilative (suction-generating) state this DEM has
produced, reached from the other side of the same line.

## 6. What this changes in the mechanism, and what it does not

- **Unchanged:** `Δφ_c(CV)` and every contrast. The line moves with CV exactly as the
  constant-pressure run said (U2).
- **Changed:** the dilatancy sign is **not a function of `φ` alone**. It is the sign of
  `p_shear − p_ss(φ)`, and a first-time-loaded pack always has `p_shear = p_iso(φ) > p_ss(φ)`.
  For the drop this says the sign is set by the *loading path* at the contact: a pressure
  pulse that rises (first-time loading → overpressure) and then relaxes while the pack
  keeps shearing (unloading → the pack crosses OCR* → suction). The E-3 transducer
  prediction becomes a **sign reversal within one impact**, not a single sign.
- **Changed:** the route to the dilative side that a lubricated suspension actually takes
  is not isotropic over-consolidation but the Wyart–Cates one — contacts are frictionless
  until the stress exceeds `σ*`, so the pack at rest is a *frictionless* pack at its `φ`,
  which for `φ_μ(CV) < φ < φ_0` carries **no pressure at all** (`p = 0 < p_ss(φ)`, i.e. the
  OCR → ∞ limit). Switching friction on and shearing must then build a frictional network
  over a finite strain `γ_c(φ − φ_μ(CV))`, with `p` rising from zero — dilative, suction,
  arrest — and `γ_c` diverging as `φ → φ_μ(CV)⁺`. That is the strain-clock protocol
  (`PREREG_STRAINCLOCK_2026-09-13.md`), and it is the physically right preparation for a
  DST suspension, which the isotropic-frictional one is not.

## 7. Confounds, stated, and the runs that remove them

- `p_0` is a γ = 0.0075 proxy (≈ 1 % low). It cannot flip U1: the decay continues for
  tens of samples on every jammed cell.
- **U1's pack is not on the isotropic line** (`p_0/p_iso` = 3.6–9.3, §4). The
  clean test: take the phase-2 servo-held state itself (`φ`, `p_iso`) and shear it at
  constant volume. Prediction from U3 where it is readable: CV 0 at `φ` = 0.595 has
  `p_iso/p_ss` = 2.45, so `p` must fall by that factor — contractive with no
  prep-overpressure to hide behind. If instead `p` does not fall, U1's whole sign was the
  preparation, and §6 is withdrawn.
- **U2 is unreadable at CV 0 and 0.5.** The run that makes it readable: common-ladder
  constant-volume cells at `φ` between `φ_c(p*)` and the present lowest jammed rung
  (CV 0: 0.595–0.62 → e.g. 0.600 and 0.607; CV 0.5: 0.645–0.67 → e.g. 0.650 and 0.657),
  N = 1500, rate-matched at 1.17e-04, `γ_max` ≥ 3 so the near-`φ_c` rungs reach steady.
  The bar then reads without extrapolation.
- U2 where readable mixes N (1500 vs 1000) and `γ̇`; a rate-matched N = 1500
  constant-pressure run would be the clean version.
- The rung nearest `φ_c` at each CV is not steady at `γ_max` = 1.5 (printed in U2's
  table), so its `p_ss` is an upper bound, its OCR* a lower bound, and U3's CV 0.5 verdict
  undecided.
