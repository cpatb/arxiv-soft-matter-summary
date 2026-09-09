# NEXT STEPS — polydispersity, particle size, and the impact regime

**Written 2026-09-09, end of an autonomous session that worked through
`NEXT_STEPS_2026-09-09.md`.** Feed this file to a new session; it supersedes that one.

**Read first, in this order (about 10 minutes):**

1. This file, §0 and §1.
2. `PREREG_ROOT_PINNING_2026-09-09.md` §6–7 — the pre-registered test that decided the
   session, and my own corrected motivation for it.
3. `SYNTHESIS_2026-09-08.md`, **with §2(b), (d) and (e) read against §0 below.**

`NEXT_STEPS_2026-09-09.md` is now history. Two of its statements were already stale
when it was written: P5 *was* in `MECHANISM` §8 and in the derivation's P-table.

---

## 0. What changed, and it is more than the defect the last session flagged

The handoff asked for one repair (quality gates counted but not applied). Doing it
turned up two further defects, and the second is larger than the first.

### ✅ Done — the gate defect (task 0 of the old file, all four parts)

`analyse_phase8.py` has a `GATES` registry, `apply_gates()`, and
`analyse_group(recs, gate_names)` so **every** estimator honours the same gate. `main`
prints every contrast under **all four gate sets, side by side, never one alone**. The
ladder was extended by +0.10 and +0.12 at **all three** CV levels (18 cells, 0 failures)
so every gate set resolves n = 3 seeds everywhere:

| gate | Δ`φ_c` (CV 0 → 0.5) | slope | sd |
|---|---|---|---|
| ungated | **+0.0607** | 0.1214 | 42.6 |
| quasi-static only | **+0.0446** | 0.0892 | 6.0 |
| steady only | **+0.0471** | 0.0941 | 4.6 |
| both | **+0.0479** | 0.0957 | 6.8 |

**Quote +0.045 to +0.061 with the gate named.** `test_shear.py` fails if a gate is
defined and no estimator honours it, and asserts the gate *moves* the primary estimator
rather than merely dropping a row. 42 tests pass.

### ✅ Controlled — the shear rate was never a controlled variable

`gdot = γ_max/(steps·dt)` and `dt = DT_SAFETY·PWaveTimeStep()` is set by the **smallest
particle**, so widening the PSD raises the rate: **3.000e-05 at CV = 0 against 1.437e-04
at CV = 0.5, a factor 4.78**, against a within-CV seed spread of at most 1.10. Nobody had
recorded this. 30 cells were re-run with `steps` chosen per cell so every arm sits at
`gdot` = 1.430e-04 (match **checked** from the realised rate: 1.000×/1.000×/1.002×).

Removing it moves Δ`φ_c` by **+0.0002 / +0.0038 / +0.0107** across the gate sets — every
one **positive**, so the confound made the published number **conservative**, which is
the direction `φ(I)` rheology predicts. At most 21 % of the contrast.
`grade_rate_matched.py`, `out/rate_matched_verdict.txt`.

### 🔴 THE ONE THAT MATTERS — E2 does not measure `φ_c`

Pre-registered in `PREREG_ROOT_PINNING_2026-09-09.md` **before** the deciding cells ran,
graded in its §7. Outcome **SPLIT**, the third branch the pre-registration named.

- **The root follows the bottom of its own fitting window.** Move the bottom by one rung
  (which the quasi-static gate does): the root moves **1.15** of that. Add two rungs at
  the **top**: it moves **0.06** of that.
- **On a ladder that actually brackets the root, E2 refuses on all three seeds**, with
  its exponent saturated at the bound `b` = 8.00; the roots it would have returned are
  0.4511 / 0.4522 / 0.4594, far *below* the ladder. A power law cannot represent a rung
  below its own root, where `p → 0` and `log p → −∞`, so the fit runs the exponent to the
  ceiling instead of moving the root. **E2 returns a number only when the whole ladder
  sits above the root — and then the number is the bottom rung minus a small offset.**
- **Every ladder is centred on the project's own published isotropic `φ_J`**
  (`phase8_criticalstate.PHI_J_ISO`). So Δ`φ_c` and Δ`φ_J` have to agree, and they do, to
  2.6 %.

**Withdrawn:** that the critical-state slope 0.1214 "matches the independently measured"
isotropic 0.1246 (one number read twice, not corroboration); and that `φ_c` ≈ `φ_J` at
`μ_p` = 0.3.

**Survives:** the *sign and rough size* of the CV effect, as a ladder-relative root; and
**the critical state is genuinely bracketed** — a fit-free readout at CV = 0 gives
`p_ss` = 1.16e-07 / 2.72e-07 / 7.65e-07 / 3.74e-06 / 4.73e-04 and `z_nr` = 2.38 / 2.84 /
3.30 / 3.94 / 4.58 at φ = 0.5145 … 0.5945, with rattlers 0.796 → 0.066. **φ_c(CV = 0) is
between 0.5545 and 0.5945**, and E2's 0.5695 sits inside that bracket.

**Untouched:** **P3**, the matched-`φ` pressure collapse (5.4×–43×; no fit, no root, no
ladder placement) — this is now the project's strongest simulation result. And **P4**,
the frictionless anchor, where E4's crossing is a genuine bracketed interpolation and
`Z` = 2d = 6 is exact.

### 🔴 And the obvious replacement has a mis-set constant

`Z_ISO_FRICTIONAL = 4.0` is the `μ_p → ∞` limit. The production arm runs at `μ_p` = 0.3,
where nothing fixes the count. `measure_z_iso.py` extrapolates `z_nr` to zero pressure on
the project's **own** isotropic ladders (same records, same construction `phase2_phiJ.py`
uses) and gets **4.535 ± 0.072 over 18 cells** — 0.535 contacts out, worth 0.019–0.024 in
`φ_c` at the ladders' own `dz/dφ`. And it **drifts with CV** (lognormal 4.609 → 4.377),
putting a CV-dependent bias of ≈ +0.0045 into E4's own contrast. `E4m` (the crossing
against the measured count) is now reported beside `E4`, never instead of it.

⚠ **The consequence is structural, not a constant.** A corrected E4 needs `z_iso` read
off the isotropic ladder — which is the measurement the shear ladder is already centred
on. **At `μ_p` = 0.3 this protocol has no `φ_c` estimator independent of `φ_J`.**

### ✅ Two questions from the old file, answered

- **§2 item 2, "does P5 break the degeneracy?" — NO.** `p5_degeneracy.py`. The shape
  effect is present in `φ_J`, measured with **no shear at all** on an imposed-pressure
  ladder: 2.735× between families against `φ_c`'s 2.823×, 3.2 % apart. So the shape
  effect is a property of the packing geometry, not of the critical state. To use P5
  against the rheological rival one would have to show that `φ_m` responds to skewness
  differently from every packing-fraction measure this project has. ⚠ A referee notes
  `φ_m = φ_J` is itself an unverified identification (Wyart–Cates-type `φ_m` is
  stress-dependent) [UNVERIFIED], which is why the claim is the weaker "P5 does not
  **establish** a discrimination".
- **§2 item 5's blocker — CLEARED.** `check_regime_flag.py`: the experiment's `arrested`
  flag is **not** stale. It reproduces as `d_arrest < 1.2` on all 160 rows, `d_arrest` is
  already the dimensionless spread ratio, and the cut sits in an **empty band 0.738
  wide** — it could move anywhere in ±0.37 without relabelling a drop. The 27 rows with a
  finite `d_arrest` and `arrested = False` are simply the non-arrested arm.

---

## 1. Where the science stands now

**Established, and not touched by any of the above:**

- **P3, at matched `φ`:** raising CV 0 → 0.5 collapses the constant-volume grain pressure
  **5.4×–43×** (quote the range), with rattlers 0.004 → 0.448 and `Z_nr` 6.1 → 4.4. No
  fit, no estimator, no ladder placement. In both pairs the higher-CV cell is the denser,
  so the ratios are conservative.
- **`φ_J` rises with CV**, +0.0623 over CV 0 → 0.5, measured on an imposed-pressure
  ladder — a different protocol and code path from anything above. **This is the sound
  packing measurement**, and the one the "widening raises the jamming fraction" claim
  should rest on.
- **P5, the shape effect:** 2.735× (isotropic) / 2.823× (shear) between the lognormal and
  uniform families at matched `d₃₂` (~9 %). Real; not a discriminator (see §0).
- **P4, the frictionless anchor:** `φ_c(μ_p = 0)` = 0.6347 ± 0.0010 against an
  independently measured `φ_J(μ_p = 0)` = 0.6332. ⚠ A reviewer notes this is a *post-hoc
  estimator choice*: 1 of 4 estimators passes, and the pre-registration named the
  quantity, not the estimator. The physical argument for E4 there is independent of the
  outcome (Z = 2d = 6 is the definition of frictionless isostaticity, and E2's exponent
  is at its bound on 2/3 seeds), but it should be quoted with that caveat.
- **`μ_s`** = 0.359 (μ_p = 0.3) vs 0.105 (μ_p = 0), factor 3.42. Internal check, not an
  anchor.

**Open, in order of how much it matters:**

1. 🔴 **There is no trusted `φ_c` estimator at `μ_p` = 0.3.** See §2 item 1 — this now
   blocks every quantitative `φ_c` statement, and it is the first thing to fix.
2. **The mechanism is still degenerate with plain suspension rheology**, and P5 does not
   break it (§0). The three separators remain experimental: pore-pressure **sign**,
   **strain history**, and a boundary in **(U₀, D₀)**.
3. **The lab's silica is on the wrong side of the measured `φ_c`** (0.44 against
   ~0.57–0.59) and friction cannot close the gap (floor 0.5574, gap 0.117).
4. **The drainage length is unmeasured and moves `ID` by 10⁵.**
5. **`α`, the dilatancy coefficient, is assumed CV-independent and untested** — needs
   constant-*pressure* shear.

---

## 2. The programme, re-ordered by what today changed

### 1. 🔴 AN ESTIMATOR THAT IS NOT A PROPERTY OF ITS OWN LADDER — do this first

Everything quantitative about `φ_c` waits on it. Three parts, in order:

**(a) Put every CV level on the SAME ABSOLUTE `φ` ladder.** The whole inheritance
problem is that the ladder is centred on `φ_J(CV)`. A common grid — say φ = 0.55, 0.58,
0.61, 0.64, 0.67, 0.70 — cannot transmit Δ`φ_J` into Δ`φ_c`, because the window is
identical for every CV. 3 CV × 3 seeds × 6 rungs = 54 cells. This is the single
highest-value run left. Pre-register that Δ`φ_c` on a common ladder stays within, say,
0.015 of the centred-ladder value; if it collapses, the centred-ladder contrast was
inheritance.

**(b) An estimator that tolerates a rung below the root.** E2 cannot: `log p → −∞`. Fit
`p` in **linear** space with `p = c·max(φ − φ_c, 0)^b`, or locate the crossing of a
declared small pressure, or use the `z_nr` crossing. Whatever is chosen, grade it by the
**window test in `grade_root_pinning.py`** — move the window and require the answer not
to move — which is now the repo's standard for a `φ_c` estimator.

**(c) The E4 repair ladder is ALREADY RUNNING** — see §5. `downlad2`, 30 cells, cv = 0.2
and 0.5 with offsets −0.06 … +0.02, which brackets `z_nr` = 4 at every CV so E4 and E4m
can report the full contrast for the first time. Analyse with
`--pattern "phase8_downlad*.jsonl"` and read the `E4` and `E4m` rows together.

### 2. E-6, the nozzle test — still the cheapest bench experiment, and now relatively more valuable

`drop_size_test.py`. The three drainage readings move the boundary velocity in three
**directions** — `U₀* ∝ D₀`, `∝ 1/D₀`, or flat — 3.6 / 14.4 / 7.2 m/s at 2 mm against a
shared 7.2 at 4 mm. Decided by the **sign**, so it needs no band thickness, no `ID`
prefactor and no boundary location. ⚠ **Cornstarch, not silica** (silica's boundary is
278 m/s); score the **regime**, not the spread factor; Newtonian control at matched
viscosity and density.

It rises in priority because it does not touch `φ_c` at all.

### 3. The salt re-analysis — NOW UNBLOCKED

`check_regime_flag.py` cleared the flag. The salt ladder (0/10/50/100 mM) is the only
designed axis already on disk, and the effective-diameter model has `φ_c^eff` moving
0.23 → 0.56 across it, so the sign of `φ − φ_c` should flip inside data already taken.

⚠ **The naive fraction is a trap and the script says so.** `n_arrested/n` and
`n_arrested/n_scored` run in **opposite directions** across the ladder (0.76→0.58 against
0.79→1.00), because scoring coverage falls with salt; and the `U₀` ranges are not matched
(0 mM reaches 4.95 m/s, 100 mM only 2.63). **Match `U₀` before comparing anything.**

### 4. Constant-PRESSURE shear — closes two open items at once

A servo on the isotropic part of `velGrad` holding `p` constant while shearing gives
(i) `α` directly, testing the derivation's untested assumption, and (ii) `φ_c(p)`, the
critical-state line. ~15 lines on `shear.py`; the risk is servo stability, so shake it
out on one cell. **And note it sidesteps §2 item 1 entirely**: at constant pressure the
critical state is where dilatancy vanishes, which is not a root fitted from below.

### 5. A repulsive contact law — demoted, but not dropped

Was item 1 in the old file. It decides whether a charge-stabilised pack has a critical
state anywhere near the lab's `φ` = 0.44 — but it is a `φ_c` measurement, so it inherits
§2 item 1's estimator problem and should wait for it. When it runs, the cheap route is
inflating interaction radii by the double-layer thickness while computing `φ` from the
**true** radii (Yade `aabbEnlargeFactor` + `interactionDetectionFactor`), and it has a
pre-registered prediction attached: `CV_eff/CV = d₁₀/(d₁₀+2h)` = 0.74 at low salt, so
Δ`φ_c` should come out ~26 % smaller at the same true CV.

### 6. E-2, E-0, and the salt 2×2

Unchanged from the old file. E-2 (raise `φ` at high CV until the solid branch returns)
needs `φ` steps ≲0.03 and drift <0.015, and the existing `φ` ladder is 2.1× too short.
E-0 (raw size distributions of the two bench samples) is still the one thing only the
author can do, and it decides whether "800 nm" is a number mean at all. The salt 2×2
(CV low/high × salt low/high) is the design that separates the two opposite salt effects
that a one-dimensional ladder confounds.

---

## 3. Standing hazards in this codebase

Carried forward, plus what today added.

- 🔴 **A `φ_c` estimator must be graded by moving its window.** `grade_root_pinning.py`
  is the standard. E2 fails it; E4 passes it structurally (it interpolates between the
  two bracketing rungs) rather than empirically.
- 🔴 **`Z_ISO_FRICTIONAL = 4.0` is the `μ_p → ∞` limit** and is wrong at the production
  friction. Use `E4m` / `Z_ISO_MEASURED` and say which.
- 🔴 **`gdot` is an OUTPUT of the packing, not a control.** Anything that changes the
  smallest particle changes the shear rate. `analyse_phase8` flags a contrast whose
  realised `gdot` moves more than 1.20×.
- **The ladders are centred on `PHI_J_ISO`.** Any agreement between a shear-derived
  packing fraction and the isotropic one is suspect until the window test is run.
- **Estimator choice depends on friction.** E2 is unusable at `μ_p` = 0 (exponent at
  bound) and now unusable at 0.3 (window-dependent); E4 finds no crossing at CV = 0.5 on
  the *centred* ladder. Say which estimator, on which ladder.
- **E1 is retracted** — on a fixed-rate ladder `I → 0` is the densest rung. A test pins
  it; do not reinstate it.
- **`sync_back.sh` moves only raw jsonl and logs.** Derived `.csv`/`.txt` are re-derived,
  never copied. And **sync the WSL tree before grading**: a referee caught a verdict
  graded on a partial copy this session.
- **`E:\CPB\Drop Impact Analyses…` ends in an untypeable U+F028.** Glob for it.
- **Never grade a knob at one setting of the knob it interacts with.**
- Raw data on `E:` is **read-only**.
- **Do not run two `scripts/launch_phase8.sh` invocations concurrently** — it re-copies
  `*.py` into `~/pj` on every call, and a concurrent `cp` can hand a starting worker a
  half-written module. `scripts/launch_0909b.sh` shows the standalone pattern.

---

## 4. Verification commands

```
cd C:\Users\compl\Desktop\Projects\polydisperse_jamming
python -m pytest test_shear.py test_poro.py -q          # 42 tests
python derivations/dilatancy_sign_mechanism_check.py    # 21 checks
python analyse_phase8.py                                # CV arm, all four gate sets
python analyse_phase8.py --pattern "phase8_p5uniform*.jsonl" --out p5_verdict.txt
python grade_root_pinning.py                            # the window test
python grade_rate_matched.py                            # the shear-rate control
python measure_z_iso.py                                 # the isostatic count
python p5_degeneracy.py                                 # P5 vs the rheological rival
python check_regime_flag.py                             # the experiment's regime label
python grade_p4_anchor.py ; python phi_axis_offset.py
python drop_size_test.py ; python experiment_coverage.py
python make_mechanism_figures.py --check                # byte-identical regeneration
```

Yade runs go through WSL. Use a standalone launcher against the copy already in `~/pj`
if anything else is running:
`wsl bash -lc 'bash ~/pj/scripts/launch_0909b.sh'`

---

## 5. In flight when this was written

- **`downlad2`, 30 cells, RUNNING.** cv = 0.2 and 0.5 with offsets −0.06 … +0.02, so
  `z_nr` = 4 is bracketed at every CV and E4/E4m can report the full CV contrast. Sync
  with `cp ~/pj/out/phase8_downlad2*.jsonl <repo>/out/` and analyse with
  `python analyse_phase8.py --pattern "phase8_downlad*.jsonl" --out downlad_verdict.txt`.
  **Read the `E4` and `E4m` rows together** — they differ by the isostatic-count
  correction, and neither is "the" answer.
- Nothing else is running.
