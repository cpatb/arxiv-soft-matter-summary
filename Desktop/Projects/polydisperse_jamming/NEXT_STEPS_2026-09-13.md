# NEXT STEPS — after the Z_ISO audit and the first constant-pressure runs

**Written 2026-09-12.** Supersedes `NEXT_STEPS_2026-09-11.md`, which is now history
(its §2 items 1, 2 and 5 are all discharged below).

**Read first, in this order (about 10 minutes):**

1. This file, §0 and §1.
2. `PREREG_CONSTP_2026-09-12.md` **§7** — the constant-pressure shakeout's outcome.
3. `PREREG_ZISO_2026-09-12.md` **§8** — why the `Z_ISO` blocker turned out not to be one.
4. `PREREG_CONSTP_PROD_2026-09-12.md` — the production ladder, and §6 for its outcome if
   the run has landed.

---

## 0″. 2026-09-13 — §0′.1 IS DISCHARGED. Read this first; it supersedes §0′.1 entirely.

The strain-clock campaign was pre-registered as **`PREREG_STRAINCLOCK_2026-09-13.md`**
(not `-09-12`; the two documents citing the old name were corrected in the same commit),
ran 45 cells with 0 failures, and is graded in that file's §7 by `grade_strainclock.py`
→ `out/strainclock.txt`.

**The answer to §0′.1 prediction (iii) is COLLAPSE.** `γ_c` at matched distance from the
critical line is the same at every CV: `R` = 1.069 ± 0.092 (CV 0.2) and 0.914 ± 0.103
(CV 0.5) against a pre-registered [0.67, 1.50], with the candidates 2.0–3.0× apart, so
the test was informative and **M4 — "CV changes the network-building kinetics" — is
demoted on its own discriminator**. M2a is the surviving packing channel.

Three more things landed with it:

- **`Δφ_c` now has two protocols.** `p_ss(φ)` from constant-volume shear crosses
  `p*` = 3e-04 at 0.5876 / 0.5987 / 0.6399, giving `Δφ_c`(0→0.5) = **+0.0523** against the
  constant-pressure **+0.0533** — agreement to 0.0011 between methods sharing no machinery.
- **§0′.1 prediction (ii) is REFUTED**: `γ_c` does not grow toward `φ_μ`. It is flat at
  ≈ 0.15 across `φ` and CV, and near the critical line the observable stops existing
  (the pressure fluctuates by its own size). Whether a divergence is hidden under those
  fluctuations is the one gap left and it needs a different observable, not more rungs.
- **Two registered bars were wrong when written** (a quasi-static `I` bar that no rate can
  meet at an unjammed rung; an isostatic bracket at `d+1` = 4.0 where this project has
  measured 4.535 at `μ_p` = 0.3). Neither was moved; both failing numbers stay printed
  with a sharper companion beside them.

Two analyses of data already on disk landed the same night and are in the Status card:
`window_shearedge.py` (the DST window's width trend is a property of which lower edge you
use, and is not established under any of them; the TRAVERSE is robust at 2.20–2.45
windows per unit CV) and `cv_switch.py` (the switchable band is 0.0533 ± 0.0018 wide in
`φ`; which channel of Eq. (11) dominates is set by `Δφ`, with the crossover computed).

**The programme now is the Status card §3 list**, whose item 1 is the ensemble-averaged
`p(γ)` near `φ_μ` that closes the gap above.

---

## 0′. EVENING UPDATE, 2026-09-12 — read this before §1 and §2; it supersedes §1 items 1, 2, 6 and §2.1

Three things landed after this file was written. Every number below is printed by the
script named beside it.

**(a) P1 is RESOLVED: `Δφ_c`(CV 0 → 0.5) = +0.0533 at matched `I` ≈ 6–7e-03**
(`grade_constp_rm.py` → `out/constp_ratematch.txt`; `PREREG_CONSTP_RATEMATCH_2026-09-12.md`
§5). Rate spread 1.001 (was 4.70); `φ_c` = 0.5911 ± 0.0010 / 0.6020 ± 0.0013 /
0.6444 ± 0.0043; CV 0.2 bit-identical to production on `φ_ss`, `p_ss_hat` and `α` (the
constant-pressure protocol is deterministic, unlike the isotropic one). §2.1 below is
discharged. ⚠ The move is not linear in CV: +0.0109 over 0 → 0.2, +0.0425 over 0.2 → 0.5.
**The H-3 "sign-flip candidate" is withdrawn**: `φ_c − φ_iso(p*)` = +0.0078 / +0.0063 /
+0.0028, never negative; the −0.0072 was the rate.

**(b) `α` is RETRACTED** (`derivations/alpha_from_dem.md`, evidence in
`alpha_from_dem_numbers.txt`): the fitted rate read the servo clamp (86–90 % of the dense
transit at `SERVO_EDOT_CLAMP`), the loose branch started unjammed and was 60–100 % servo
compression, and the branch labels were swapped. §1 item 6 is closed the other way —
`α` is not measurable by that protocol at all. `poro_impact.ALPHA_DIL` = 1.0 is an
unmeasured import again; E-3 is a sign test again.

**(c) The undrained SIGN of a fresh pack, from the 147 constant-volume records already on
disk** (`cv_unloading.py` → `out/cv_unloading.txt`, `cv_unloading.md`; bars U1–U3 fixed in
the docstring before the first run; exit `FAILED: U2, U3`):

| bar | result |
|---|---|
| U1 — `p_ss < p_0` on every jammed cell | **PASS 108/108**, every CV, both sides of `φ_c(0)` |
| U2 — CSL from two protocols agree at `p*` to 0.01 | **FAIL as registered**: NaN at CV 0 and 0.5 (`p*` = 3e-04 lies *below* the lowest jammed constant-volume rung, 4.3e-04 / 7.1e-04 — the bar was posed outside the data); CV 0.2, the one readable level, agrees to **+0.0031** |
| U3 — `p_iso(φ) > p_ss(φ)` on overlapping rungs | **FAIL as registered**: 2.45× (CV 0, 0.595) and 5.07× (CV 0.2, 0.595) hold; **0.74× at CV 0.5, 0.645** on a rung steady 1 of 3 seeds — undecided, not a sign flip |
| U4 reported — `γ_half` of the overpressure | grows with `φ`: 0.077 → 0.135 (CV 0), 0.030 → 0.115 (CV 0.2), 0.035 → 0.070 (CV 0.5); at matched `φ` the CV 0.5 pack relaxes 2–3× faster |
| U5 reported — OCR* = `p_0/p_ss` | 1.24–1.64 at `φ` ≥ 0.62 (CV 0), 12.7 at 0.595; CV 0.5 is **2.00×** CV 0 at `φ` = 0.645 and 1.28× at 0.670 |

**And the number no bar asked for, which is the session's most important one:
`p_0/p_iso` = 3.6–9.3.** The phase-8 preparation (rate compression, then constant-volume
relaxation) sits 3.6–9.3× above the phase-2 servo-held isotropic state at the *same* `φ`
and N. `p(φ)` for a frictional pack is a history, not a function, and U1's "isotropically
prepared pack is contractive" is true of a *rate-compressed* pack. The clean test is
§0′.2 below.

**What this does to the mechanism.** The dilatancy sign is not a function of `φ`; it is
`sign(p_shear − p_ss(φ))`, set by the loading path. A DST suspension at rest is a
*frictionless* pack at its `φ`, which for `φ_μ(CV) < φ < φ_0` carries no frictional
pressure at all — the OCR → ∞ limit — so the first frictional shear must *build* `p`
from zero: dilative, suction, arrest, over a finite strain `γ_c(φ − φ_μ(CV))` that
diverges at `φ_μ(CV)⁺`. That is the strain-clock protocol, and it is the physically right
preparation for the drop; the isotropic-frictional one was not. The drop-side prediction
E-3 becomes a **sign reversal within one impact** (loading → overpressure; unloading while
still shearing → suction), not a single sign.

### The programme now, in order

**0′.1 — The strain-clock campaign (frictionless prep → `utils.setContactFriction(0.3)` →
constant-volume shear).** Pre-register as `PREREG_STRAINCLOCK_2026-09-12.md` (referenced
by `alpha_from_dem.md` §4 and `cv_unloading.md` §6 — **it does not exist yet; write it
before any cell runs**). N = 1000, lognormal, CV {0, 0.2, 0.5}, 3 seeds, `φ` ladder from
just below `φ_μ(CV)` [= `φ_c(p*)`: 0.591 / 0.602 / 0.644] to just below `φ_0(CV)`
[frictionless `φ_J`], 5 rungs, rate-matched at 1.17e-04, `γ_max` ≥ 3. Predictions to
write down: (i) `p` rises from ~0 at every rung above `φ_μ(CV)` — the first dilative,
suction-generating state this DEM has produced; (ii) `γ_c` (strain to `p` = half its
steady value) grows as `φ → φ_μ(CV)⁺` and no rung below `φ_μ(CV)` jams; (iii) **the
mechanism bit: `γ_c` plotted against `φ − φ_μ(CV)` COLLAPSES across CV if polydispersity
acts only through the critical line, and does NOT if it changes the network-building
kinetics** (the 0.23 rattler fraction at CV 0.5 is a reason it might not). (iii) is the
test that separates "CV shifts `φ_c`" (M2a-degenerate) from "CV changes the clock" (M4),
which no steady-state measurement can. Launch: copy `scripts/run_constp_ratematch.sh`
into `scripts/run_strainclock.sh`; a `--friction-on-at-shear` flag in `shear.py::_prepare`
(prep at `μ_p` = 0, set friction before stage C); invoke by path from PowerShell in a
background Bash call; cells must fit ~600 s each. Also register there the OCR test: unload
an isotropic pack at `φ` to `p₂ < p_ss(φ)` and require `p` to *rise* — graded against U5's
OCR* per rung.

**0′.2 — The clean U1: shear the phase-2 state itself.** Servo-consolidate to (`φ`, `p_iso`)
with the phase-9 servo (it already exists), then constant-volume shear. 6 cells: CV 0 at
`φ` = 0.595 (prediction from U3: `p` falls 2.45×, contractive, no prep overpressure to hide
behind) and CV 0.5 at 0.645 (undecided; the constant-pressure line puts `p_ss` ≈ `p*`, i.e.
a fall of ≈ 1.7×; a *rise* would make CV 0.5 dilative at its own `φ_c` and U3's violation
real). If CV 0's `p` does not fall, U1's sign was the preparation and `cv_unloading.md` §6
is withdrawn.

**0′.3 — Make U2 readable.** Constant-volume rungs between `φ_c(p*)` and the present lowest
jammed rung: CV 0 at 0.600 and 0.607, CV 0.5 at 0.650 and 0.657, 3 seeds, rate-matched,
`γ_max` ≥ 3 (the near-`φ_c` rungs are the ones that were not steady at 1.5). 12 cells.
⚠ N = 1500 to match the ladder cannot run from an agent session (600 s cap); either the
author runs it from a terminal, or it runs at N = 1000 with the finite-size caveat stated
in the pre-registration, not after.

**0′.4 — Documents that are now stale and must catch up in the same commit as 0′.1's
pre-registration:** `MECHANISM_2026-09-08.md` §1 ("above `φ_c` the pack dilates" —
history-dependent now), §6a/E-3 (α magnitude → sign reversal within one impact), R4 H-3
candidate (withdrawn), `Δφ_c` = +0.053 rate-matched; this file's §1 items 1, 2, 6 and §2.1.

**0′.5 — Lab side, unchanged and still the cheapest discriminator the author owns:** H-2
(flow-curve shape on a CV ladder, §2.3) and E-6 (§2.4). Nothing in (a)–(c) touches them.

**Verification for this update (all ran 2026-09-12 evening):** `python grade_constp_rm.py`
(exit 0), `python derivations/alpha_from_dem_numbers.py` (`ALL CHECKS PASSED`),
`python cv_unloading.py` (exit 1, `FAILURES: ['U2', 'U3']` as reported above),
`ruff check cv_unloading.py` clean.

---

## 0. What changed

### 🔑 `φ_c` now has an estimator that takes NO input from `φ_J`

`shear.shear_at_constant_pressure` shears at a fixed grain pressure with a servo on the
volumetric strain rate, so the steady-state solid fraction **is** `φ_c(p)` — no ladder, no
power-law root, no fitting window, no isostatic contact target. That closes the standing
complaint from `PREREG_ROOT_PINNING_2026-09-09.md` §7 (restated in
`PREREG_ZISO_2026-09-12.md` §8): at `μ_p` = 0.3 every previous route read the isotropic
`φ_J` the shear ladders were centred on.

**Shakeout: 10 cells, 0 failures** (`out/constp_shakeout.txt`, graded by
`grade_constp.py` against bars fixed before launch).

- `φ_ss`(p = 3e-04, CV 0, N = 600) = **0.5935 ± 0.0012** on the seven valid cells.
- It moves **0.0023 across a 5× range of servo gain** and **0.0032 across preparations
  0.09 apart in `φ`** (bars 0.005). Those were the two ways this protocol could have been
  meaningless.
- Against the **isotropic** ladder at the same pressure (0.583314 ± 0.000725, N = 1500,
  n = 9) it is **denser by +0.0102**, the predicted direction: steady shear reaches a
  denser state than isotropic compression.
- Against **`E4m`'s `φ_c`(CV 0) = 0.5902** on the common ladder: **|Δ| = 0.0033.**
  ⚠ N = 600 vs 1500, one seed, one pressure. A corroboration, not a measurement.
- **`α` exists with the right sign structure** (fitted rate negative 10/10, 0.66–4.08
  decades of travel), lands in 1.5–7.3, and the dilating and contracting branches
  disagree by up to **4.41×**. Not quotable on one seed; the production run tests whether
  the asymmetry beats the seed scatter.

🔴 **Two failures worth carrying, both diagnosed rather than excused.**

- **S1 (the servo holds) failed at 4 of 10**, and it is not droop: above `K ≈ 0.25` the
  servo drives the cell **through unjamming** (`p_hat` visits 1e-08 against a 3e-04 set
  point, `z_nr` falls to 0.02), and the mean then reads **+20 % HIGH** because the
  excursions are one-sided. Below the boundary the opposite limit bites — at `K` = 0.05
  the servo is too slow and the offset returns at +8.7 %. Usable window **K ∈ [0.1,
  0.25]**, and the residual offset is priced on the project's own ladder:
  `dφ/d ln p` = 0.00906, so 20 % in pressure is **0.0017 in `φ`**.
- **S2 as written had passed on a different question.** It averaged `φ_ss` over every
  cell, including the two gains that unjam. S2′/S3′ restrict to the cells that passed
  S1's own pre-registered criterion — no new bar — and **on the original six cells S2′
  would have been UNREADABLE**, because only one stable gain existed. The four low-gain
  cells (prescribed in advance by §5's own S1 remedy) are what made it readable. S3′
  still rests on a single shared held gain.

### ✅ `Z_ISO` is measured, and it is not load-bearing

`PREREG_ZISO_2026-09-12.md`, `out/ziso_audit.txt`, **8 of 8 gates**. Two defects were
found in `measure_z_iso.py` *before* anything was re-measured:

- **D1** it read the shared glob `phase2_*.jsonl`, which `capfix` began writing into on
  09-09, with no cap filter and no dedup — re-running it unchanged reads 4.549 instead of
  4.535 and takes the `μ_p` = 1 arm from 1 cell to 2. Fixed with `--cap` and dedup.
- **D2** `N_LOW_RUNGS = 5` was documented as *rungs* and implemented as *records*, so with
  three seeds per rung the extrapolation ran on **two distinct pressures** on ladders that
  offer four. Fixed with an explicit `--rule`; the shipped behaviour is kept as `broken`
  so the historical artifact can be reproduced on purpose.

**And the correction does not matter for the contrast.** Moving the whole lognormal target
table to the rung-correct values (CV 0: 4.609 → 4.557, CV 0.5: 4.377 → 4.343) moves
`E4m`'s `Δφ_c` by **+0.0003**, 0.0497 → 0.0501, and the per-level absolutes by ≤ 0.0012.
A uniform shift of the target cancels out of a contrast exactly; only its CV-dependence
survives. **So the E4/E4m split on H-3 (0.0069 vs 0.0126) cannot be a target effect of
this size** — that has to be something else, and naming it is now the open question.

🔑 **New, from an observable with no `φ_J` in it:** the CV drift in `z_iso` is a
**lognormal** effect (−0.213) and nearly absent on the zero-skew families (uniform
−0.056, bidisperse −0.012). Uniform's headline-looking −0.193 under the shipped rule
rested entirely on a two-pressure ladder. That is the project's CV/skewness separation
arriving from a third direction.

**Decisions taken, both rules pre-registered:** do **not** adopt a repaired-cap `Z_ISO`
table (`capfix` covers `μ_p` = 0.3 only at CV 0.4/0.5 and nothing at CV 0 or 0.2, so a
half-adopted table would put the **iteration cap on the CV axis under test**); and
**`PHI_J_ISO` is not updated under any outcome** — it records where the ladders were
physically placed, not an estimate. The "pending decision" was two decisions.

### ✅ The `_prepare` extraction is verified, not asserted

Stages A–C were lifted verbatim out of `shear_at_constant_volume` so both drivers share
the preparation. `check_prepare_refactor.py` runs the **pre-refactor** module (out of git,
via a sidecar, because `~/pj` is a file copy and not a checkout) and the current one on
the same cell: **bit-identical on all five continuous channels**.

---

## 1. Where the science stands

**Established:**

- 🔑 `Δφ_c`(CV 0 → 0.5) = **+0.049 to +0.055** on a φ grid identical at every CV.
- 🔑 `φ_c` measured at fixed pressure with no `φ_J` input agrees with `E4m` to **0.0033**
  at CV 0 (N = 600).
- P3 at matched `φ`: **3.2×–177×** collapse of the grain pressure.
- The critical state bracketed with no estimator: `z_nr` = 4 between 0.570–0.595 (CV 0,
  0.2) and 0.620–0.645 (CV 0.5).
- `φ_J` rises **+0.0623** with CV, isotropic; `W` = 0.0586–0.0757 and `dφ_μ/dCV` = +0.117,
  giving a window traverse at ΔCV = 0.45–0.64.
- `z_iso`(`μ_p` = 0.3) = 4.535, stable to **0.024 contacts** across three fit rules, and
  its CV drift lives on the tail family only.

**Open, in order of how much it matters:**

1. **The production constant-pressure ladder** (§2.1) — 18 cells, launched 2026-09-12.
   It is the first `Δφ_c`(CV) from a route with no `φ_J` in it.
2. **H-3's split is now unexplained.** `E4` and `E4m` differ only in the target, and the
   target is worth 0.0003 on the contrast. Whatever splits them at 0.0069 vs 0.0126 is
   not that, and nothing here says what it is.
3. **The mechanism is still degenerate**, tree now 2-bit: **H-4** (a boundary in
   `(U₀, D₀)`, one nozzle, cornstarch) separates pore-fluid {M2′} from packing
   {M2a, M4}; **H-2** (flow-curve *shape* on a CV ladder, on a rheometer the lab has)
   separates M4 from M2a. `PREREG_DST_WINDOW` §5.
4. **The lab is on the wrong side**, priced: 1.6–2.3 window-widths below, needing 23 %
   volumetric compaction to reach the bottom edge.
5. **The drainage length is unmeasured and moves `ID` by 10⁵** — E-6 attacks it.
6. **`α` is measurable now but not measured.** The branch asymmetry (up to 4.41×) is
   either physics or noise and the production run is the first chance to tell.

---

## 2. The programme, re-ordered

### 1. Re-run the CV 0.5 arm quasi-statically — that is the whole of what P1 needs

The production ladder **landed** (18 cells, N = 1000, `out/constp_prod.txt`,
`PREREG_CONSTP_PROD_2026-09-12.md` §6): servo held 18/18, steady 18/18, and

- `φ_c` = **0.593254 ± 0.000641** (CV 0) and **0.601989 ± 0.001313** (CV 0.2), from a
  route with no `φ_J` in it;
- **P4 at 4.8e-05 / 3.9e-05** — the critical state forgets a preparation 0.09 apart in
  `φ` to five decimals, on 3 seeds each side;
- the admitted CV 0 → 0.2 contrast is **+0.00874** against the common ladder's +0.0105,
  agreeing to 0.0018;
- 🔴 but **all six CV 0.5 cells fail the quasi-static bar** at `I` = 1.4–1.5e-02, because
  `γ̇` is an output of the packing and runs **4.70×** across the CV ladder. So
  `Δφ_c`(CV 0 → 0.5) is **UNRESOLVED** and the raw +0.0412 is not quotable.

**What to run:** the same grid at CV 0.5 with `--shear-steps 300000` (`I` × 0.68 → 1.0e-02,
just inside the bar) — and note that **matching** `γ̇` to CV 0's 6e-05 needs 940 000 steps,
≈3900 s/cell, which this harness cannot hold (§3). Clearing the bar is what is available
here; rate-matching is the better science and is an author-run job. Say which one the
number came from.

**Two rule repairs to make in the same commit**, both diagnosed in §6:

- The admission `z` floor is **absolute (2.0) and refuses polydisperse packs for being
  polydisperse**: `z_series` is `z_ALL`, and at CV 0.5 the rattler fraction is 0.42–0.52,
  so `z_all` ≈ 2 by construction (the project's own M3 result). Read relatively,
  `z_win_min/z_ss` = 0.77–1.00 on **all 18** cells and nothing unjammed. Use
  `z_nr_series` (now recorded) or a relative floor.
  ⚠ **Repairing it does NOT rescue P1** — the rate refusal applies to all six CV 0.5
  cells independently. Do not report the repair as having moved the verdict.
- `grade_constp_prod.py` needs its explicit `N_REQUIRED`: the abandoned N = 600 / 800
  attempts wrote refusal records under the same tag, and 8 of them reached the analysis
  directory because a cleanup chain beginning with `pkill` aborted on pkill's own exit
  status before its `rm` ran. The filter is in place and prints what it excluded.

🆕 **And carry the flag:** `φ_ss` − `φ_iso` runs +0.0099 → +0.0064 → **−0.0072** across
CV 0 / 0.2 / 0.5 (the last on refused cells). A sign flip in the shear-minus-compression
offset would be a new statement about polydispersity and the critical state, and it is the
first candidate explanation for **H-3**.

### 2. If P1 lands, the constant-pressure ladder replaces the estimator argument

`φ_c(p)` at the four isotropic rungs {3e-05, 1e-04, 3e-04, 1e-03} gives the
shear-minus-compression offset at four pressures **with no estimator on either side** —
which is H-3's question asked in a way that cannot be about a contact-count target. That
is the natural next campaign, and it needs a servo with an integral term (§3).

### 3. H-2, the flow-curve shape ladder — still the cheapest thing the author can do

Fixed `φ`, fixed mean size, ΔCV ≈ 0.5 on the lab's rheometer. Score the shape by a
pre-registered statistic (max local `d log η / d log σ`, bootstrapped), never by eye, and
require the two candidates' predictions to be ≥ 1.6× apart on the ladder actually run
before calling the test informative.

### 4. E-6, the nozzle test (H-4) — unchanged

`out/drop_size_test.txt`. Decided by the **sign** of three drainage readings, so it needs
no band thickness or `ID` prefactor. Cornstarch, not silica (`U₀*` = 278 m/s at 2.5 mm on
the lab's silica); score the regime, not the spread factor; Newtonian control at matched
viscosity and density.

### 5. The salt re-analysis — still unblocked, still untouched

⚠ Match `U₀` before comparing anything: `n_arrested/n` and `n_arrested/n_scored` run in
*opposite* directions across the salt ladder because scoring coverage falls with salt.

---

## 3. Standing hazards in this codebase

Carried forward, plus what today added.

- 🔴 **NOTHING SPAWNED BY A `wsl.exe` CALL SURVIVES THAT CALL EXITING**, and this is the
  binding constraint on every campaign from an agent session. Measured 2026-09-12: not
  under `nohup`, not under `setsid nohup`, and **not while a second WSL session is held
  open** — tested with a bare `sleep 300`, which was gone by the next call. Combined with
  the harness's 600 s cap on one call, **the longest possible cell is ~600 s**, which is
  N = 600 at 200 000 steps and NOT N = 1500 (~900–1000 s). The first `N` = 1500 production
  attempt was killed at 600 s having produced nothing, and the run was re-registered at
  N = 600 (`PREREG_CONSTP_PROD_2026-09-12.md` §2a).
  **What unblocks N = 1500** — an author decision, not an agent's: `vmIdleTimeout = -1`
  in `%USERPROFILE%\.wslconfig` with a service-style launch, or simply running the
  campaign from a terminal the author leaves open. Everything else in the pipeline is
  ready for it; only the session lifetime is missing.
- 🔴 **Invoke WSL scripts from PowerShell, by file path.** Git Bash rewrites a
  `/mnt/c/...` argument before `wsl.exe` sees it (it became
  `C:/Program Files/Git/mnt/c/...`), and `wsl bash -lc "..."` with double quotes eats
  every `$`.
- 🔴 **A gate computed over every cell can pass on a different question.** S2 passed while
  averaging in cells that periodically unjam. Restrict to the cells that passed the
  admission criterion *the pre-registration already contains* — and if that makes the gate
  unreadable, say so and run more cells.
- 🔴 **A servo is a knob and its answer must not depend on it**, and the check is a gain
  sweep, not an argument. Report `z_win_min` per cell: a cell can hold `⟨p⟩` on average
  while visiting an unjammed state.
- 🔴 **`phase2_*.jsonl` is a shared glob two campaigns write into.** Any analysis reading
  it must select a cap explicitly and dedup by cell key.
- 🔴 **The isotropic protocol is not deterministic at CV = 0** — 15 duplicated keys return
  a different `φ` (worst 5.0e-03). Never quote a CV = 0 isotropic number with a
  seed-scatter error bar alone.
- 🔴 **`Z_ISO_FRICTIONAL = 4.0` is the `μ_p → ∞` limit** (measured 4.535 ± 0.072, drifting
  with CV on the lognormal family). Use `E4m` beside `E4`, never instead — and note that
  the target is worth only 0.0003 on the contrast.
- **`gdot` is an OUTPUT of the packing, not a control.**
- **E2 is unusable** at every setting tried.
- **Two scripts computing the same quantity must be made to agree before either is
  believed**, and the first thing to check on a disagreement is the RECORD SELECTION.
- **A test that returns EXACTLY zero has probably not run.**
- **Give concurrent invocations different tags** — a record is ~30 kB, far above
  `PIPE_BUF`, so appends to one shard file are not atomic.
- **`sync_back.sh` moves only raw jsonl and logs.** Derived `.csv`/`.txt` are re-derived.
- Raw data on `E:` is read-only; `E:\CPB\Drop Impact Analyses…` ends in an untypeable
  U+F028, so glob for it.

---

## 4. Verification commands

```
cd C:\Users\compl\Desktop\Projects\polydisperse_jamming
python -m pytest test_shear.py test_poro.py test_psd.py test_cap.py \
                 test_dst_window.py test_ziso.py test_constp.py -q      # 102 tests
python grade_constp_prod.py        # the production ladder, P1-P5
python grade_constp.py             # the shakeout, S1-S7
python grade_ziso.py               # the Z_ISO audit, G1-G7
python grade_common_ladder.py      # the deciding constant-volume run
python dst_window.py               # the window and the traverse
python grade_capfix.py             # H-1
python analyse_phase2.py           # convergence BY ARM + duplicate conflicts
python measure_z_iso.py            # --cap default --rule r3 by default now
python make_mechanism_figures.py --check     # byte-identical regeneration
```

Yade runs go through WSL against the copy in `~/pj`, launched **from PowerShell**:

```
wsl bash /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/run_constp_prod.sh
wsl bash /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/run_prepare_check.sh
```
