# The `p̂` ladder runs out of dynamic range, and it does so exactly where today's new claims live

**2026-09-19.** Every number printed by `ladder_dynamic_range_2026-09-19.py` into
`ladder_dynamic_range_2026-09-19.txt` (rc 0, ALL CHECKS PASSED). Read-only; the live
`out/phase2_wcbidi*.jsonl` campaign is excluded from the load.

> ## RESULT
>
> `φ_J` is never measured — it is **extrapolated** from a four-rung `p̂` ladder. An
> extrapolation needs the ladder to move by more than the noise, and nobody had asked how
> far it moves. One number per cell group settles it:
>
> ```
>     SPAN/SD  =  [ φ(p̂ = 1e-3) − φ(p̂ = 3e-5) ]  /  seed scatter of φ(p̂ = 3e-5)
> ```
>
> **12 of 34 cell groups on disk fall below SPAN/SD = 2**, and **2 of them have a
> physically inadmissible NEGATIVE span** — compressing a packing harder cannot lower its
> packing fraction, so those ladders run backwards and their extrapolated `φ_J` is
> meaningless rather than merely imprecise. That second statement needs no threshold.
>
> ### Which committed numbers this touches
>
> | quantity | worst SPAN/SD | |
> |---|---|---|
> | `paper`: `\dPhiJLogn`, the headline shift | 3.8 | ok |
> | `paper`: `\dPhiJUnif` | 8.0 | ok |
> | `paper`: `\dPhiJBidi` | 6.8 | ok |
> | `paper`: `\phiJFrictionless` | 7.6 | ok |
> | `paper`: `\phiJMonoMuThree` | 14.4 | ok |
> | `wc_endpoints`: the frictionless CV shift | **1.2** | 🔴 unreadable |
> | `wc_friction`: the `μ_p = 1.0` window | **−0.2** | 🔴 inadmissible |
> | `wc_friction`: the `μ_p = 3.0` endpoint | **−0.1** | 🔴 inadmissible |
>
> **The manuscript is untouched. All three of the claims I made today rest on groups that
> fail.** Both of today's derivation pages carry a correction pointing here.
>
> **Verdict: SUPPORTED, and it is a protocol result rather than a physics one.** It says
> which cells can be read, not what they say.

---

## 1. The mechanism, and it is the same one twice

SPAN/SD averaged over families:

| `μ_p` \ CV | 0.0 | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 |
|---|---|---|---|---|---|---|
| 0.0 | 7.6 | 8.7 | 15.2 | 7.0 | 8.8 | **4.1** |
| 0.3 | 14.7 | 13.9 | 20.6 | 5.5 | **3.5** | 6.2 |
| 1.0 | 3.6 | — | **0.8** | — | — | **−0.2** |
| 3.0 | **−0.1** | — | — | — | — | — |

Mean SPAN/SD: **9.9** at `μ_p ≤ 0.3` against **1.0** at `μ_p ≥ 1.0`; **11.9** at CV ≤ 0.2
against **5.0** at CV ≥ 0.4. Both axes degrade it, and for the same reason:

> **Friction and polydispersity both RAISE `φ_J`.** At a fixed pressure ladder, a packing
> with a higher `φ_J` sits further above its own jamming point and is stiffer, so it
> compresses less between `p̂ = 1e-3` and `3e-5`. The ladder has less to say precisely where
> the physics being studied is strongest.

This is an instance of the laboratory's standing detector-versus-its-own-quantum problem, in
a new place: the *estimator's input range* collapses, rather than the *observable's*
resolution.

## 2. Two checks that could have failed

**(a) Does SPAN/SD predict an independent symptom?** Ladder non-monotonicity is a separate
observable. Below the bar, **1 of 12** ladders is monotone; above it, **44 of 90** — a
6× difference.

⚠ **A bar I registered here FAILED and is not loosened.** I wrote "Spearman(SPAN/SD,
monotone fraction) > 0.3" without computing it; it is **0.274**. It was also the wrong
statistic: the effect is a **floor**, not a gradient, so a rank correlation across the whole
range dilutes it. The failing number stays printed and the split above is reported beside it.

**(b) Do the convention and the inadmissibility agree?** The two groups with negative spans
(`lognormal, μ_p = 1.0, CV = 0.5`, span −0.00090; `lognormal, μ_p = 3.0, CV = 0`, span
−0.00058) are both below the SPAN/SD bar. They did not have to be.

## 3. What it means for the two pages written today

- **`wc_friction_endpoint_2026-09-19.md`** used `phi_J_lin` at `μ_p = 1.0` and `3.0`. Those
  are the two inadmissible groups. **Its channel must change to `φ` at the matched lowest
  rung `p̂ = 3e-5`, which is fit-free.** On that channel the page's conclusion survives:
  `μ_p = 0.3` captures **79 %** of the friction drop rather than 77 %, and the `μ*` ratio
  between CV 0 and CV 0.5 is **1.327** rather than 1.243. The numbers move; the verdict does
  not.
- **`wc_endpoints_move_together_2026-09-19.md`** used the frictionless CV = 0.5 cells
  (SPAN/SD = 1.2). Its `Δφ_0` is therefore read off a ladder spanning about one seed
  scatter. `paper/main.tex` already warned about this qualitatively — *"the frictionless
  ladders are not monotone in `p̂` at this seed count"* — and this is that warning with a
  number on it.

## 4. What it means for the protocol

The four-rung ladder `(1e-3, 3e-4, 1e-4, 3e-5)` was chosen at `μ_p ∈ {0, 0.3}` and
CV ≤ 0.5, where it works. It does **not** transfer to higher friction or wider
distributions, and the existing fit-quality gate does not catch that: a ladder can have an
acceptable `r²` while spanning nothing.

**The cheapest fix is a higher top rung.** The ladder's top is `p̂ = 1e-3` and its purpose is
to give the extrapolation leverage; `phase2_phiJ.py` records that a `3e-3` rung was dropped
on a timing argument — *"`analyse_phase2.P_HAT_FIT_MAX` excludes it from the extrapolation
anyway"*. For high-friction or high-CV cells that reasoning should be revisited, because
those are exactly the cells where the accessible range is too small.

⚠ **This is a suggestion, not a measurement.** Whether a higher rung restores the span at
`μ_p = 1.0` is testable in three cells and has not been tested.

## 5. What this diagnostic cannot do

1. **SPAN/SD = 2 is a convention.** Every group's value is printed so another bar can be
   applied; the *ordering* does not depend on it. The negative-span statement does not
   depend on it at all.
2. **The `sd` is a 3-seed estimate**, 2 dof, itself noisy by roughly 40 %. A group at 1.5 and
   one at 2.5 are not reliably apart.
3. **It says an extrapolation is unreliable, not that `φ_J` is wrong.** A degenerate ladder
   whose four points happen to agree still *brackets* the jamming point; it cannot locate it
   better than that bracket.
4. **It does not replace the fit-quality gate.** A ladder can have good span and still be
   non-monotone from one bad rung.
5. **It is computed on `φ`, not on `phi_J_lin`**, so it is independent of the fit — which is
   the point, but it also means it cannot tell a bad fit from a bad ladder.

## 6. Referee

⚠ **Unrefereed.** The diagnostic itself is arithmetic on committed data and the two checks in
§2 were written to be able to fail (one did). The *interpretation* in §1 and the protocol
suggestion in §4 have not been attacked.
