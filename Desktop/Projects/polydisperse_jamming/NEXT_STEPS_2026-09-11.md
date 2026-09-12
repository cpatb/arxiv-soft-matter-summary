# NEXT STEPS — after the common ladder landed and a fourth mechanism appeared

**Written 2026-09-09, end of a session that worked through `NEXT_STEPS_2026-09-10.md`.**
Feed this to a new session; it supersedes that file. `NEXT_STEPS_2026-09-09.md` and
`NEXT_STEPS_2026-09-10.md` are now history.

**Read first, in this order (about 10 minutes):**

1. This file, §0 and §1.
2. `PREREG_COMMON_LADDER_2026-09-09.md` **§8**, the outcome of the deciding run.
3. `PREREG_DST_WINDOW_2026-09-09.md` §1–3 and §5 — the new mechanism and the four
   pre-registered hypotheses.

---

## 0. What changed

### ✅ THE BIG ONE — `Δφ_c` is a material effect, not ladder inheritance

`NEXT_STEPS_2026-09-10.md` §2 item 1(a) ran. Pre-registered in
`PREREG_COMMON_LADDER_2026-09-09.md` before any cell launched; graded by
`grade_common_ladder.py` into `out/common_ladder_grade.txt`. **54 cells, 0 failures,
3 seeds at every CV.**

All three CV levels were put on an **identical absolute** `φ` grid — 0.545, 0.570, 0.595,
0.620, 0.645, 0.670 — which **cannot** transmit `Δφ_J` into `Δφ_c`, because the window is
the same six numbers at every CV.

| estimator | gate | `Δφ_c` (CV 0 → 0.5) |
|---|---|---|
| **E4** | ungated | **+0.0554** |
| E4 | steady | +0.0521 |
| **E4m** | ungated | **+0.0497** |
| E4m | quasistatic | +0.0506 |
| E4m | steady | +0.0494 |

**H_real on every estimator and every resolved gate set**, against a centred-ladder
reference of +0.041…+0.051. **Quote `Δφ_c` = +0.049 to +0.055 on the common ladder.**

Per level (E4, ungated): 0.5745 ± 0.0015 / 0.5850 ± 0.0010 / 0.6299 ± 0.0048.

**This restores what 2026-09-09 withdrew.** The doubt was that the ladders were centred
on `φ_J`, so `Δφ_c` could be `Δφ_J` read twice. On a window that cannot do that, the
answer is the same. `Δφ_c` is measured.

### ✅ E4's window-independence is now MEASURED, not argued

- The pre-registered §5a test (drop the bottom rung) passed with a move of **exactly
  0.0000 on all 9 groups** — and that exact zero is the proof it is **vacuous**: E4
  interpolates between the two rungs bracketing `z_nr` = 4, which is never the bottom
  rung here, so any correct implementation passes it and so would a constant. The bar
  stayed and a stricter test was added beside it.
- **§5a′, the placement test.** Same cells, ladders at genuinely different absolute `φ`
  (common 0.545–0.670 vs centred 0.515–0.695 at CV 0 and 0.577–0.757 at CV 0.5):
  **worst |E4(common) − E4(centred)| = 0.0051 over 9 cells**, bar 0.010, range
  0.0003–0.0051 so it had room to fail. Moving a ladder by 0.06 in `φ` moves E4 by
  0.005.

### ✅ P3 sharpens, and the fit-free bracket is now on a common grid

At six *shared* `φ` (realised mismatch ≤ 8.2e-05, against the centred ladders' 0.0023):
**3.2× to 177.5×** at and above `φ` = 0.595, all > 1 as pre-registered. `z_nr` crosses 4
between **0.570 and 0.595** at CV 0 and CV 0.2, and between **0.620 and 0.645** at
CV 0.5 — a full grid interval of separation with no estimator at all.

### 🆕 A FOURTH mechanism: the WINDOW, not the edge

`PREREG_DST_WINDOW_2026-09-09.md`, `dst_window.py` → `out/dst_window.txt`,
`MECHANISM_2026-09-08.md` §7b. **No new simulation** — a re-reading of ladders already
run, and the reason it is new is that every mechanism before it asked where `φ` sits
relative to **one** packing fraction.

Wyart–Cates turns on **two**, and the position *between* them fixes the **type** of
thickening: CST below `φ_μ`, DST / shear jamming between, jammed above `φ_0`.

```
   mu_p = 0.3 : W = 0.0575   dphi_mu/dCV = +0.128 +- 0.013  ->  2.23 windows/CV
   mu_p = 1.0 : W = 0.0733   dphi_mu/dCV = +0.115 +- 0.020  ->  1.56 windows/CV
   => dCV = 0.45 to 0.64 traverses the WHOLE window at fixed phi
```

`mu_p` = 1 is the SATURATED friction, i.e. Wyart-Cates' `phi_mu` as this protocol reaches
it. **REPAIRED numbers** (`capfix`, 180 cells, `out/capfix_grade.txt`); the pre-repair
estimate was 0.50-0.65, so only the third significant figure moved.

**M4:** the window is narrow, and polydispersity slides it bodily upward about twice as
fast as it is wide, so a modest widening carries a material clean across it — a
**threshold in CV**, with no pore fluid. **Its own discriminator is the SHAPE of the flow
curve** (S-shaped ⇄ monotonic), which needs the rheometer the lab has rather than the
pressure transducer M2′ needs.

🔴 **It does not rescue the lab's working point; it prices it.** At `φ` = 0.44 the
material is **1.6–2.3 window-widths BELOW** the window at every CV, and raising CV moves
it *further* below. So at the stated `φ` the bench observation reads as **M2a**, and
impact must supply **31 % more packing fraction (23 % volumetric compaction)** merely to
reach the bottom edge, 44 % to be jammed at all stresses.

### 🔴 Two defects in the isotropic campaign, both on the axis under test

**(a) The frictionless arm is censored, and asymmetrically.** `analyse_phase2` drops
`converged = False` cells; split by arm — which its **pooled** SELECTION CHECK cannot do —
the rate runs **0.92 → 0.33** across CV at `μ_p` = 0 against **1.00 → 0.75** at
`μ_p` = 0.3. So the window differences a heavily censored arm against a light one. All
82 non-converged cells stopped at exactly `MAX_ITER` (slow, not wrong); 18 met the
physical tolerance and were rejected on the iteration count alone. The bias is −0.0047
per unit CV in the frictionless arm, predicted in advance to be ≈28 % of the measured
window narrowing with the same sign — **measured after the repair: 31 %.**

**(b) 🔴 The isotropic protocol is NOT deterministic at CV = 0.** `load_all`'s docstring
asserted duplicates are redundant because "the protocol is deterministic given the seed".
**15 duplicated keys return a different `φ` and a different iteration count** from the
same (cv, μ_p, seed, `p̂`, N, protocol) — worst 5.0e-03, median 2.1e-03, against a seed
scatter of 2.2e-04. **All 15 are at CV = 0.0.** Which one survives is set by the
**filename sort order**, and enumerating the choices puts `W`(CV = 0) anywhere in
**0.05594–0.05911, a span of 0.0032 = 11× the quoted standard error.**

Found because two scripts disagreed (`grade_capfix` gave `dW/dCV` = −0.0071 where
`dst_window` gave −0.0169) and the disagreement was chased rather than papered over.
`analyse_phase2` now counts the conflicts and prints the worst; the CSVs regenerate
byte-identically, so no number moved when the detection was added.

⚠ **Reach, PROPAGATED rather than guessed.** `Δφ_J` and `dφ_μ/dCV` survive (≤3.5 %
move). For `dW/dCV` I first guessed "~40 % of the total change" from the leverage arm;
computed properly (`grade_capfix`, moving `W`(CV = 0) across its 0.0032 spread and
refitting) it is **±0.0023 against a fit error of ±0.0035** — *smaller* than the fit
error, so the slope survives at t ≈ −2.8 combined. The guess was wrong by ~4×, which is
why it is a line of code now. **And none of it touches the shear campaign**: different
protocol, different code path, `grade_common_ladder` reports 0 duplicate-state pairs.

### ✅ `capfix` LANDED — 180 cells, 0 nonzero exits

`out/capfix_grade.txt`, graded by `grade_capfix.py` against bars quoted from
`PREREG_DST_WINDOW_2026-09-09.md` §5; full account in that file's **§7**.

- **Convergence PASSES**: lognormal `mu_p` = 0, CV 0.5 went **0.33 -> 0.83** (bar 0.80);
  every level now resolves 3 seeds. `phi_J`(`mu_p` = 0, CV 0.5) = **0.690401 +- 0.004665
  on 3 seeds**, against 0.683324 on ONE seed before.
- **`dW/dCV` at `mu_p` = 0.3 = -0.0117 +- 0.0035 (fit) +- 0.0023 (systematic)** -- the
  declared THIRD BRANCH. **-0.0169 is retracted as an over-estimate.** Censoring
  accounted for 31 %, against ~28 % predicted in advance.
- 🔑 **At the SATURATED friction the window is FLAT** (+0.0040 +- 0.0030, t = +1.34, not
  resolved), which is exactly the rigid-window premise M4 rests on.
- 🔴 **The control FAILED and was MIS-SPECIFIED**, which is the finding. 2 of 12 cells sat
  AT the old cap and legitimately moved (`converged` False -> True at 205k/212k iters).
  Of the 10 real controls, 5 bit-identical and 5 not -- and **`phi` is identical EXACTLY
  when the iteration count is, on all 10**. So the TRAJECTORY is deterministic and the
  STOPPING is not. Pre-dates the knob; the knob's own effect is at or below the 5.0e-03
  CV = 0 spread.
- ⚠ **A DECISION IS PENDING and is deliberately not taken.** `PHI_J_ISO` and
  `Z_ISO_MEASURED` were measured under the DEFAULT cap and every `phi_c` already reported
  rests on them; adopting the repaired set would move E4m everywhere, including the
  common-ladder verdict. The two tests that pin them now read the default-cap rows and
  say why, and `test_a_REPAIRED_phi_J_set_exists_and_differs` fails if the repaired rows
  are deleted or if the two sets stop differing.

### (superseded) RUNNING: `scripts/launch_capfix.sh`

Repairs both defects at `PJ_MAX_ITER` = 1e6 (a new `packing.max_iter()` env knob,
recorded as `max_iter`/`tol_met` in every record, appended to the cell key when
non-default, 8 tests exercising **both** settings). CV = 0 is the **bit-identity
control**. It also adds the `μ_p` = 1 CV ladder the window's own definition needs.
**Grade with `python grade_capfix.py`** — written before the run, bars quoted from
`PREREG_DST_WINDOW` §5 H-1.

---

## 1. Where the science stands now

**Established:**

- 🔑 **`Δφ_c`(CV 0 → 0.5) = +0.049 to +0.055**, on a fitting window identical at every
  CV. Not inheritance. `out/common_ladder_grade.txt`.
- 🔑 **P3 at matched `φ`: 3.2×–177.5×** collapse of the constant-volume grain pressure,
  no fit, no root, no ladder placement.
- **The critical state is bracketed with no estimator**: `z_nr` = 4 between 0.570–0.595
  (CV 0, 0.2) and 0.620–0.645 (CV 0.5).
- **E4 is placement-independent to 0.0051** over a 0.06 ladder shift.
- **`φ_J` rises +0.0623 with CV**, isotropic, a different protocol.
- **`W` = 0.0586–0.0757 and `dφ_μ/dCV` = +0.117**, giving a window traverse of
  ΔCV = 0.50–0.65.

**Open, in order of how much it matters:**

1. **The mechanism is still degenerate**, but the tree is now 2-bit and one branch is
   cheap: **H-4** (a boundary in `(U₀, D₀)`, one nozzle) separates pore-fluid {M2′} from
   packing {M2a, M4}; **H-2** (flow-curve *shape* on a CV ladder, a rheometer the lab
   has) separates M4 from M2a. `PREREG_DST_WINDOW` §5.
2. **`Z_ISO` is still structural.** E4 and E4m split on H-3 (offset constant in CV at
   0.0069 vs moving at 0.0126) and they differ *only* in the isostatic target, whose
   measured value drifts with CV (4.609 → 4.377). Nothing here is settled until that is.
3. **The lab is on the wrong side**, now priced: 1.6–2.3 window-widths below, needing
   23 % volumetric compaction to reach the bottom edge.
4. **The drainage length is unmeasured and moves `ID` by 10⁵.**
5. **`α`, the dilatancy coefficient, is untested** — needs constant-*pressure* shear.

---

## 2. The programme, re-ordered

### 1. Grade `capfix` — it is running, do this first

`python grade_capfix.py`. Three things to read, in order:

- **The control.** CV = 0, `μ_p` = 0 must come back **bit-identical**. ⚠ Given §0(b),
  this is now a *stronger* test than when it was written: it is also the first direct
  measurement of whether the CV = 0 non-determinism reproduces under a re-run. **If the
  control fails, do not read it as "the knob is broken" without first checking whether
  the difference is the same size as the 15 duplicate conflicts (2.1e-03 median).** Those
  are two different diagnoses and the grader's message says only the first.
- **Convergence** at lognormal `μ_p` = 0, CV 0.5 must exceed 0.80.
- **`dW/dCV`**, graded against H-1a (< 0.006, it was censoring), H-1b (≤ −0.012 with
  t < −2, it is real), or the declared third branch.

### 2. H-2, the flow-curve shape ladder — the cheapest thing that separates M4 from M2a

On a rheometer the lab already has, at fixed `φ` and fixed mean size, over a CV ladder
spanning ΔCV ≈ 0.5. Score the **shape** by a pre-registered statistic (the maximum local
`d log η / d log σ`, bootstrapped over the measured points), **not by eye**, and require
the two candidates' predictions to be ≥1.6× apart on the ladder actually run before
calling the test informative. ⚠ Guy *et al.* report Wyart–Cates failing for a predominant
fraction of large particles, i.e. the broad-lognormal end — so the ladder must be dense
enough to locate the change *before* it reaches that end.

### 3. E-6, the nozzle test (H-4) — unchanged and still the cheapest bench experiment

`drop_size_test.py`. The three drainage readings move `U₀*` in three **directions**
(3.6 / 14.4 / 7.2 m/s at 2 mm against a shared 7.2 at 4 mm), so it is decided by the
**sign** and needs no band thickness or `ID` prefactor. ⚠ Cornstarch, not silica; score
the regime, not the spread factor; Newtonian control at matched viscosity and density.

### 4. The salt re-analysis — still unblocked, still untouched

`check_regime_flag.py` cleared the flag last session. ⚠ **Match `U₀` before comparing
anything** — `n_arrested/n` and `n_arrested/n_scored` run in *opposite* directions across
the salt ladder because scoring coverage falls with salt.

### 5. Constant-PRESSURE shear — closes `α` and `φ_c(p)` at once

~15 lines on `shear.py`; the risk is servo stability, so shake it out on one cell.

### 6. A repulsive contact law — still demoted

It is a `φ_c` measurement, and although §0 repaired the estimator problem it did not
repair `Z_ISO`. Pre-registered prediction attached: `CV_eff/CV = d₁₀/(d₁₀+2h)` = 0.74 at
low salt, so `Δφ_c` should come out ~26 % smaller at the same true CV.

---

## 3. Standing hazards in this codebase

Carried forward, plus what today added.

- 🔴 **The isotropic protocol is not deterministic at CV = 0** (§0(b)). Never quote a
  CV = 0 isotropic number with a seed-scatter error bar alone.
- 🔴 **Two scripts computing the same quantity must be made to agree before either is
  believed.** Today's non-determinism was found exactly that way, and the culprit was a
  missing dedup in the newer script — so the *first* thing to check on a disagreement is
  the record selection, not the physics.
- 🔴 **A test that returns EXACTLY zero has probably not run.** §5a passed at 0.0000 on
  9 of 9 groups because it cannot fail by construction. Ask what would have had to be
  true for a bar to fail before quoting its pass.
- 🔴 **`Z_ISO_FRICTIONAL = 4.0` is the `μ_p → ∞` limit** and is wrong at the production
  friction (measured 4.535 ± 0.072, and it drifts with CV). Use `E4m` beside `E4`, never
  instead.
- 🔴 **`gdot` is an OUTPUT of the packing, not a control.** The common ladder was run
  deliberately UNMATCHED, at the same setting as the centred ladders, so the comparison
  differs in one knob only; the rate correction (+0.0002…+0.0107, always conservative)
  is applied on top, never re-derived.
- **E2 is unusable** — at `μ_p` = 0 (exponent at bound), at `μ_p` = 0.3 on a centred
  ladder (window-dependent), and on the common ladder (unresolved at every gate set,
  because a power law cannot represent rungs below its own root).
- **The quasi-static gate removes nearly everything on a ladder that brackets the root**,
  because the rungs below it run at `I` ~ 0.1. Report all four gate sets; quote ungated
  and steady.
- **Match a claim's own subset before extending it.** The window analysis pairs the two
  frictions **per seed** (`cv_real` agrees to 1.1e-16 on all 15 pairs); the unpaired
  difference-of-slopes gives a different error bar and a different verdict.
- **`sync_back.sh` moves only raw jsonl and logs.** Derived `.csv`/`.txt` are re-derived,
  never copied.
- **Never edit a shell script while bash is executing it**, and never build a WSL script
  through a double-quoted `wsl bash -lc "..."` — the outer shell eats every `$`, and
  `bash -n` still reports SYNTAX OK on the wreckage. (It did today.) Write the file, strip
  CR, copy it, and check the md5 on both sides.
- **Give concurrent invocations DIFFERENT tags.** Three `phase8` runs sharing one tag
  would have appended to one jsonl; a record here is ~30 kB, far above `PIPE_BUF`, so
  those appends are not atomic.
- Raw data on `E:` is read-only; `E:\CPB\Drop Impact Analyses…` ends in an untypeable
  U+F028, so glob for it.

---

## 4. Verification commands

```
cd C:\Users\compl\Desktop\Projects\polydisperse_jamming
python -m pytest test_shear.py test_poro.py test_psd.py test_cap.py test_dst_window.py -q
                                                        # 70 tests
python grade_common_ladder.py                           # the deciding run, all bars
python dst_window.py                                    # the window and the traverse
python grade_capfix.py                                  # H-1, once capfix lands
python analyse_phase8.py --pattern "phase8_commonlad*.jsonl" --out common_ladder_verdict.txt
python analyse_phase2.py                                # now prints convergence BY ARM
                                                        # and the duplicate conflicts
python grade_root_pinning.py ; python grade_rate_matched.py ; python measure_z_iso.py
python make_mechanism_figures.py --check                # byte-identical regeneration
```

Yade runs go through WSL against the copy in `~/pj`:
`wsl bash -lc 'bash ~/pj/scripts/launch_capfix.sh'`
