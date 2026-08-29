# Overnight 2026-08-06/07: what changed

Two tracks. First, implementation plans for the top of the drop-impact open-problems
slate. Second, a pass through the projects folder making physics and statistics
contributions. **102 tests across four projects, all passing.** Four headline claims in
existing work turned out to be wrong or unsupported, and all four are corrected in place
with the evidence.

---

## Track 1: implementation plans

`Desktop\Papers\STATE_OF_THE_FIELD\IMPLEMENTATION_PLANS.md`

Bench-ready protocols for P1, A1, A4 and A2 from `OPEN_PROBLEMS_MASTER.md`, in the order
they should run. Each carries the same five headings so they can be compared: what gets
decided, the pre-registration, the run sheet, the analysis, the stop rules.

The pre-registrations are the point. Each one fixes the estimator, the prediction and the
failure branch **before** data arrive, because every candidate this programme has killed
died to a confound that was decided after the fact. A1's power calculation is in there
(12 drops per cell gives ±0.15 on the exponent, enough to separate 0 from 1), as is A4's
resolution floor rule, which is the guard that would have prevented the retracted 80×
crater span.

---

## Track 2: project contributions

### `uroflow_splashback` — rewritten, 44 tests where there were none

Five defects, two reversed conclusions. Full derivations in `PHYSICS_NOTES.md`.

1. **The target geometry was self-inconsistent.** v1 integrated the trajectory to a
   vertical plane, then scored splash with `sin(angle from horizontal)`, which for a
   vertical wall is the *tangential* fraction. It scored a head-on impact as the driest
   configuration. The verbal reasoning in the README was right for a horizontal target;
   the code integrated to a vertical one, and the two were never reconciled.
2. **The jet was treated as parallel.** A falling jet accelerates and thins by
   continuity; over 0.75 m the diameter falls to 63% of exit, so v1 over-stated the
   impact Weber number by 1.6×.
3. **Breakup was uncalibrated and 3-4× too long.** v1 printed 63 cm next to its own note
   that real streams break at 15-20 cm, then used 63 cm to draw the coherent region its
   headline depended on. Now an accumulated-disturbance integral along the thinning jet,
   with `Λ = ln(r₀/δ₀)` calibrated to the observation: `Λ = 2.7-3.6`, implying a 3-7% exit
   disturbance. That is a falsifiable side-prediction.
4. **Post-breakup drop size used the jet diameter**, under-stating `We` by 1.89×.
5. **The splash index was invented.** Replaced by Mundo (dry) and Cossali (wetted)
   thresholds on the normal component, so the model now predicts *whether* it splashes,
   not only which case is worse.

**The cleanest new result:** at a vertical wall the splash-driving velocity is `v₀cos θ`,
which contains no launch height; at a horizontal pool it is the fall speed, which contains
no aim angle. The two geometries have orthogonal levers, and conflating them is exactly
what v1 did. Practical number: standing 30 cm from a wetted wall, only aim angles steeper
than ~37° below horizontal stay under the threshold.

### `bloodstain_bpa` — an error rate, 17 new tests

The 2009 NAS and 2016 PCAST reports both criticise BPA for reporting conclusions without
error rates. `uncertainty.py` supplies one.

- **The stain-size route to impact speed is a ~3.5× error amplifier**
  (`d ln v / d ln ξ = 3.1-3.9`). A 5% stain measurement is an 18% speed error; with 10%
  parent-drop uncertainty it is 38%. This was not stated anywhere.
- **That makes the height estimate biased, not merely noisy**, and attribution shows
  parent-drop-size uncertainty carries −17.2 cm of the −22.1 cm total.
- ⚠ **The README claimed the ballistic estimate was "near-unbiased" under noise while the
  project's own output showed −18 cm.** The output was right; the README is corrected.
- **Calibration curve:** unbiased only while `σ_D₀ ≲ 4%`; −13 cm at 10%, −35 cm at 20%.
- **The deliverable is now an interval:** 95% credible interval [0.46, 1.33] m against a
  true 1.15 m. Honest, and 87 cm wide.
- A first-order Jensen correction was tried and **made it worse**; kept in the code, off
  by default, so the negative result stays reproducible.

### `frisbee` — the direction ordering is only half resolved, 13 tests

`P5_Synthesis/q_significance.py`, findings in `SIGNIFICANCE_FINDING_2026-08-07.md`.

The synthesis reports "best = 210°, worst = 270°". Paired over the 12 common launch
points, Holm-corrected, with bootstrap CIs:

| comparison | ΔQ | 95% CI | verdict |
|---|---|---|---|
| 90 vs 210 | −0.0040 | [−0.0095, +0.0018] | **not resolved** (4/12 points) |
| 90 vs 270 | +0.0086 | [+0.0035, +0.0148] | resolved (11/12) |
| 210 vs 270 | +0.0127 | [+0.0056, +0.0205] | resolved (10/12) |

**"270 is the worst of the three" is supported. "210 is the best" is not** — 210 and 90
are indistinguishable. `\QbestHeading` should be withdrawn and the H3 verdict softened to
partially supported. The sheltering headline (24-28% of the 10 m reference) is unaffected.

Separately, the Monte Carlo standard error is 9-14× smaller than the spatial spread, so
the sweep is adequately sampled — but `Q_min-Q_max` is a spatial range and must not be
presented as an error bar.

**A provenance trap, with a guard added.** `Q.npz` was built with `--exclude-throws hammer`
while `runDownstream.sh` passes no such flag, and the shipped artifact predates the
provenance keys. Using all three throw types would have inflated the binomial error by 1.5×
and made the sweep look under-sampled. The analysis now recovers the subset by matching and
says which route it took; it recovers `right_backhand, right_forehand` exactly.

### `beach_sand_dilatancy` — the protocol's own §2, made computable, 22 tests

`dilatancy_design.py`, findings in `DESIGN_FINDINGS_2026-08-07.md`. This project had no
code. Three results, before any sand is collected.

1. **A notation collision.** `PROTOCOL.md` uses `φ` for porosity in the Kozeny-Carman line
   and for packing fraction everywhere else. Reading one for the other changes `k` by 7.6×.
   It cost exactly that error while writing the module. The code never uses the bare name.
2. **The rig as specified is in the wrong regime, and the fix is free.** The hardening
   ratio depends on the intruder only through `D/V`. A footstep is `D/V = 0.1 s`; the
   proposed 20 mm sphere at 2 m/s is `0.01 s`, a factor of 10 into the drained side.
   Matching needs either a 200 mm sphere or **the same sphere dropped from 2 mm instead of
   0.2 m**. This is the most useful line in the document, because it costs nothing now and
   would be expensive to find later.
3. **Q4 is reachable, and the simple model cannot span Strader's range.** The cavitation
   ceiling `V_cav ∝ d²` sits inside a bench drop for `d ≲ 0.3 mm` and out of reach above
   0.5 mm, which is the test. But requiring Strader's clean law to be uncavitated at 53 µm
   *and* pore-pressure controlled at 850 µm gives an **empty window, short by 43×**. Alone,
   the cavitation constraint bounds `ψ ≤ 0.0055`, at least 18× below the assumed 0.1. So Q4
   is where the existing model demonstrably fails, and mapping that is a better question
   than confirming Hazen. Q4 should be promoted from "one afternoon" to a named aim.

Q1 meanwhile survives comfortably: at `C_u = 2` the `d₁₀` and `d₅₀` predictions for `a_max`
differ by 59%, and the discrimination threshold is `C_u ≈ 1.23`, which every real sand
exceeds.

---

## The pattern across all four

Three of the four contributions are the same move: **compare a claimed effect against its
own noise, or a detector against its own quantum.** The BPA README claimed unbiased where
its output said −18 cm. The frisbee synthesis claimed a best direction that a paired test
cannot resolve. The uroflow model drew a coherent region 3-4× too large from a number its
own README admitted was wrong. The sand protocol's design table rested on a parameter that
the literature already bounds 18× lower.

None of these needed new data. They needed the existing numbers to be asked whether they
were bigger than their own uncertainty.

---

## Suggested next steps

- **uroflow:** one phone-camera clip of the first 20 cm would pin `Λ` and check the
  predicted threshold crossing. An afternoon, and it converts the model from
  literature-anchored to measurement-anchored.
- **bloodstain:** the spread law is the smooth-surface correlation. Replacing it with a
  blood-specific law from a Basilisk drop-impact run is the stated next step and now has a
  clear target: reduce `k`, the amplifier, or at least measure it for real blood.
- **frisbee:** regenerate `Q.npz` with the current `throwability.py` so the throw
  provenance is recorded, and make `runDownstream.sh` pass `--exclude-throws hammer`.
- **beach sand:** revise the rig spec to sweep `D/V` before ordering anything.
