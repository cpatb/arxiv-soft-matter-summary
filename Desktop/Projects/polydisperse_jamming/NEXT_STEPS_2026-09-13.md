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
