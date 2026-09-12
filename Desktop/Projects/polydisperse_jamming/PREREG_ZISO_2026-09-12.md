# PRE-REGISTRATION — the `Z_ISO` decision, and an audit of the estimator that sets it

**Written 2026-09-12, before any `z_iso` number was re-measured.** Every bar below is
fixed here. Outcome goes in §8 of this file; artifact `out/ziso_audit.txt`, grader
`grade_ziso.py`.

Supersedes nothing. It discharges the item left open in `NEXT_STEPS_2026-09-11.md` §2
item 2 ("`Z_ISO` is still structural") and the **pending decision** recorded in
`PREREG_DST_WINDOW_2026-09-09.md` §7 ("adopt the repaired `PHI_J_ISO` / `Z_ISO_MEASURED`
set, or not").

---

## 1. Why this is the highest-value analysis item

`E4`, the primary `φ_c` estimator, returns the `φ` at which `z_nr` crosses a target
count. The target is a *constant* in the code (`Z_ISO_FRICTIONAL = 4.0`) and a
*measurement* in `measure_z_iso.py` (`4.535 ± 0.072` at the production friction, with a
CV drift of `−0.232` on the lognormal family). By the artifact's own arithmetic
(`out/z_iso_measured.txt`, `dz/dφ` = 28.1 / 23.0 / 21.9 at CV 0 / 0.2 / 0.5) a target
error of +0.535 contacts is worth **+0.019 to +0.024 in `φ_c`** — between a third and a
half of the CV contrast the project's headline quotes. Nothing about `Δφ_c` is settled
until the target is.

## 2. Three defects found in the estimator BEFORE re-measuring anything

Recorded here so that the bars below cannot be read as having been chosen after the
fact. All three are properties of the code and of the file inventory, not of any result.

**D1 — `measure_z_iso.py` reads a SHARED GLOB that two campaigns now write into.**
`_load_phase2()` globs `out/phase2_*.jsonl` and filters on `converged` alone. Since the
`capfix` campaign landed (2026-09-09) that glob returns **537 default-cap records
(`max_iter` absent) plus 180 at `max_iter = 1e6`**, with no cap filter and **no dedup by
cell key** — so a ladder is now fitted through two protocols at once, and the low-pressure
rungs are double-weighted. This is the 2026-09-06 lesson ("appending to a list a gate
iterates is an edit to every gate that iterates it") relocated from a Python list to a
filesystem glob. `out/z_iso_measured.txt` as committed therefore describes a file
inventory that no longer exists.

**D2 — `N_LOW_RUNGS` counts RECORDS, not RUNGS.** The constant is documented as "how many
lowest-pressure rungs the extrapolation uses" and is set to 5. The code is

```
keep = x <= np.sort(x)[min(N_LOW_RUNGS - 1, len(x) - 1)]
```

which keeps the 5 lowest **records**. Every ladder carries 3 seeds per rung, so 5 records
is **2 distinct pressures** — the intercept that the whole `E4m` correction rests on is a
straight line through two pressures, with zero residual degrees of freedom, on ladders
that have 4 pressures available. (Same family as the 2026-09-03 entry: a 3-point
Richardson fit cannot self-diagnose. This one is a 2-point fit wearing a 5-point label.)

**D3 — the repaired cap does not cover the CV levels the headline uses.** `capfix`
targeted the cells that were failing convergence, so at `μ_p` = 0.3 it ran **lognormal
CV 0.4 and 0.5, and uniform CV 0.4 and 0.5 — and nothing at CV 0 or CV 0.2.** The common
ladder's contrast is CV 0 → 0.5. Therefore **"adopt the repaired `Z_ISO_MEASURED` set" is
not an available option**: any table built from it would be default-cap at CV 0/0.2 and
repaired-cap at CV 0.4/0.5, i.e. **the cap would be aligned with the CV axis under test.**
Feasibility, established before the bars; it is not a gate.

**D4 — the pending decision is two decisions, and one of them has already been answered.**
`PHI_J_ISO` is not an estimate of anything: it records **where the shear ladders were
physically placed**, which was the default-cap measurement, and `test_shear.py`'s comment
already says so. Updating it would make the constant disagree with every ladder on disk.
`Z_ISO_MEASURED` *is* an estimate. So the decision splits:

| constant | kind | what the repair can do to it |
|---|---|---|
| `PHI_J_ISO` | protocol record (where the ladders sat) | **nothing — it must not move** |
| `Z_ISO_MEASURED` | measurement of a material property | could move, but see D3 |
| the H-3 comparison `φ_c − φ_J` | comparison against the best current `φ_J` | should use the repaired `φ_J` |

## 3. What is being run

No new DEM. Three things, all analysis:

1. `measure_z_iso.py` gains an explicit `--cap {default,repaired,all}` and dedup by cell
   key, and reports the counts it selected and any duplicate-state conflicts.
2. The rung-selection rule is made explicit and run three ways:
   - **R-broken** — the shipped rule (5 lowest records ⇒ 2 distinct pressures).
   - **R2** — every record at the **2 lowest distinct rungs** (same pressure range as
     R-broken, balanced seeds).
   - **R3** — every record at the **3 lowest distinct rungs** (wider range, more points,
     more curvature).
3. The selected table is propagated into `E4`/`E4m` on the common ladder and the contrast
   re-graded.

## 4. The rule is selected by an ANCHOR OUTSIDE the code under test

The fit rule must not be chosen by which answer it gives at `μ_p` = 0.3, where nothing is
known. Two arms have external targets:

- `μ_p = 0`: **`z_iso = 2d = 6` exactly** for frictionless spheres (isostaticity).
- `μ_p = 3`: `d + 1 = 4` is the `μ_p → ∞` limit and the project's own anchor measures
  3.98–4.07 there.

**Selection rule, fixed now:** the rule with the smaller `|z_iso(μ_p = 0) − 6|` wins, and
it must also keep `|z_iso(μ_p = 3) − 4| ≤ 0.15`. If the two arms disagree about which rule
wins, no rule is selected and the fit-rule spread is reported as a systematic on
`z_iso(0.3)` instead.

## 5. Pre-registered predictions and bars

**G1 — the committed artifact is stale.** Re-running the shipped script unchanged returns
a `μ_p = 0.3` mean that differs from the committed 4.535, **and** a `μ_p = 1.0` count
above 1 (`capfix` added a `μ_p` = 1 CV ladder). *Predicted: both true.* If the numbers
come back identical, D1 is wrong and the glob does not reach those files.

**G2 — reproduce before extending.** With `--cap default` and dedup on, the per-arm means
must reproduce the committed artifact to **≤ 0.005 contacts** at `μ_p` = 0 and 0.3.
*Predicted: reproduces.* A miss means the dedup changed the sample, and the dedup is then
the finding rather than the repair.

**G3 — the fit rule is a systematic, and it is measured, not argued.** Report
`span = max − min` of `z_iso(μ_p = 0.3)` over {R-broken, R2, R3}. Bar: **span < 0.10
contacts** for the quoted `4.535 ± 0.072` to be a number at the precision it is quoted.
*No directional prediction is offered at `μ_p = 0.3` — nothing is known there.* At
`μ_p = 0`, where the answer IS known, the directional prediction is: **R2 lands closer to
6.0 than R3**, because R3 widens the pressure range and `z = z_c + a√p` is an asymptotic
law. *Predicted: |R2 − 6| < |R3 − 6|.*

**G4 — the CV drift survives the repair.** On the lognormal family at `μ_p` = 0.3 the
default-cap drift is `z_iso(CV 0.5) − z_iso(CV 0)` = **−0.232**. Under the selected rule,
predict the drift stays **negative** with **|drift| ≥ 0.10**. *Predicted: survives.* If it
collapses, the CV-dependent bias in `E4` was a fit artefact and `E4m` loses its reason to
exist.

**G5 — propagation: does the headline move?** Recompute `E4m`'s `Δφ_c`(CV 0 → 0.5) on the
common ladder with the selected `z_iso` table. Published: **+0.0497** (`E4m`, ungated).
Bar: **|Δ(Δφ_c)| < 0.010** for the +0.049…+0.055 headline to stand unqualified. Above
0.010 the headline is requalified in the same commit. *Predicted: < 0.010*, because a
**uniform** shift in the target cancels out of a contrast exactly and only the
CV-dependence of the shift survives — and G4 predicts that dependence is the same drift
already in the table.

**G6 — the refactor control, and it can fail.** `E4` with the theory target 4.0, computed
through the new code path, must reproduce the published **+0.0554** to **1e-9**. This is
the check that the cap filter and dedup did not change the default path. ⚠ It is *not*
vacuous: the new path re-selects records, so a wrong selection shows up here.

**G7 — the adoption decision, with its rule fixed in advance.**
Adopt a repaired `Z_ISO_MEASURED` table **only if** the repaired cap covers *every* CV
level used in a contrast that quotes it. By D3 it does not, so the pre-registered outcome
is: **do not adopt; keep the default-cap table, and record that the repaired cap cannot
re-measure the target at CV 0 or 0.2.** `PHI_J_ISO` is not adopted under any outcome (D4).
What *is* adopted, if G3 selects a rule: the corrected rung-selection.

## 6. What would make this whole item uninformative

If G3's span is large **and** the anchor arms disagree about the rule (G4's selection
fails), then `z_iso(0.3)` is not measurable to better than that span from these ladders,
and the honest end state is that `E4m` carries a systematic of that size — which is a
result about the protocol, not a number. Say so rather than picking the rule that flatters
`E4m`.

## 7. Verification

```
python measure_z_iso.py --cap default            # reproduces the committed artifact
python measure_z_iso.py --cap all                # what the stale glob was doing
python grade_ziso.py                             # G1-G7, all bars
python -m pytest test_ziso.py test_shear.py -q
```

## 8. OUTCOME — 2026-09-12, `out/ziso_audit.txt`, 8 of 8 gates pass

**The headline survives the target correction, and the surviving reason is stronger than
the one it had: `Δφ_c` is not sensitive to `z_iso` at all.**

| gate | bar | measured | verdict |
|---|---|---|---|
| G1 stale artifact | any move, and n(μ 1) > 1 | `μ` 0.3 → 4.549 (was 4.535), n(μ 1) = 2 (was 1) | **PASS** |
| G2 reproduce first | ≤ 0.005 contacts | worst 0.0038 (at `μ_p` = 0) | **PASS** |
| G3 fit-rule systematic | span < 0.10 | **0.0243** (4.535 / 4.533 / 4.510) | **PASS** |
| G4 CV drift survives | negative, ≥ 0.10 | **−0.213** (lognormal, R3) | **PASS** |
| G5 propagation | \|Δ(Δφ_c)\| < 0.010 | **+0.0003** | **PASS** |
| G6 refactor control | 1e-9 | 1.8e-05 — **the bar was unachievable** | replaced |
| G6a (replacement) | 5e-05, the artifact's own precision | 1.8e-05 | **PASS** |
| G6b (replacement) | the patch is reversible, 1e-12 | exact | **PASS** |
| G7 adoption | pre-registered rule | **DO NOT ADOPT** | **PASS** |

### What the numbers say

- **The `Z_ISO` blocker is not load-bearing on the contrast.** Moving the whole lognormal
  target table from the shipped values to the rung-correct measurement (CV 0: 4.609 →
  4.557; CV 0.2: 4.513 → 4.527; CV 0.5: 4.377 → 4.343) moves `E4m`'s `Δφ_c` by
  **+0.0003**, from +0.0497 to +0.0501. Per level the absolutes move by at most 0.0012.
  The reason is arithmetic and was pre-registered: a **uniform** shift of the target
  cancels out of a contrast exactly, so only the CV-*dependence* of the shift survives —
  and that dependence is the drift already in the table.
- **D2 is a real defect of the code with a small numerical reach.** The estimator was
  fitting two pressures while claiming five rungs; correcting it moves `z_iso(μ_p = 0.3)`
  by 0.025 contacts, i.e. **0.0011 in `φ`**. Both halves of that have to be said: the
  code was wrong, and the number it produced was right to within a twentieth of its own
  error bar. (Blast radius is a property of the caller, 2026-09-06.)
- **The fit rule was selected by an anchor outside the code**, `z_iso(μ_p = 0) = 2d = 6`.
  **R3 wins** (|R3 − 6| = 0.132 against R2's 0.148, both admissible at `μ_p` = 3), and
  ⚠ **my directional prediction was wrong**: I predicted the narrower window would win on
  a curvature argument. It lost. The explanation is in the scatter, not the curvature —
  R2's sd at `μ_p` = 0 is 0.060 against R3's 0.027, so with only four rungs the two-rung
  fit is **noise-limited, not curvature-limited**. `test_ziso.py` pins the distinction
  with a synthetic ladder carrying real curvature, where R2 does beat R3 as predicted.
- 🔑 **The CV drift in `z_iso` is a LOGNORMAL effect, and that is new.** Under the
  rung-correct rules the drift is −0.213 on lognormal, **−0.056 on uniform and −0.012 on
  bidisperse**. And uniform's headline-looking −0.193 under the shipped rule rested
  *entirely* on its CV 0.5 cell, a ladder with only two pressures; over CV 0 → 0.4 it is
  −0.059 / −0.058 / −0.056, stable across all three rules. **So the isostatic contact
  count drifts with width on the tail family and not on the zero-skew families** — which
  is the project's CV/skewness separation (`REPORT_2026-08-31.md`) arriving from an
  observable that has nothing to do with `φ_J`.
- **G6 was the wrong test, not a failed one, on two counts.** Its 1e-9 bar was compared
  against a number the artifact prints to four decimals, so no correct implementation
  could have passed; and `E4` uses the *theory* target and never calls `measure_z_iso`,
  so a correct and a broken refactor of that module both pass it. The real refactor
  control is **G2**. Both replacements are stricter in kind and the failing number stays
  printed.

### The decision, taken

1. **`Z_ISO_MEASURED` is NOT replaced by a repaired-cap table.** The repaired cap measures
   `z_iso` at `μ_p` = 0.3 only for lognormal CV 0.4/0.5 and uniform CV 0.4/0.5 — **nothing
   at CV 0 or CV 0.2**, which is where the contrast is read. A half-adopted table would
   put the **iteration cap on the CV axis under test**. This was G7's pre-registered
   outcome and it is what the data says.
2. **`PHI_J_ISO` is not updated, under any outcome.** It records where the shear ladders
   were physically placed, not an estimate of anything; `test_shear.py` already pins it to
   the default-cap row and says why. The pending decision was two decisions and this half
   was already answered.
3. **The rung-correct rule (R3) is adopted in `measure_z_iso.py`** as the default, with
   `--rule broken` retained so the historical artifact can be reproduced on purpose.
   `out/z_iso_measured.txt` is left as the committed historical artifact; G2 reproduces
   against it and nothing overwrites it.
4. **`Z_ISO_MEASURED` in `analyse_phase8.py` is left at the shipped (broken-rule) values**,
   with this outcome recorded beside it, because changing it moves `Δφ_c` by 0.0003 while
   invalidating the reproduction of every artifact already committed. The correction is
   *measured and bounded*, which is worth more than the third decimal place.

### What is still open on `Z_ISO`

- **H-3 is unresolved and this does not settle it.** `E4` and `E4m` split on whether the
  shear-minus-compression offset is constant in CV (0.0069 vs 0.0126); they differ only in
  the target, and the target correction moves `Δφ_c` by 0.0003, so **H-3's split cannot be
  a target effect of this size.** It has to be something else, and naming it is the next
  `Z_ISO` question rather than this one.
- **`E4m` still is not independent of `φ_J`.** Reading a per-cell target off the isotropic
  ladder is reading the measurement the shear ladders are centred on. Nothing here changes
  that; at `μ_p` = 0.3 this protocol still has no `φ_c` estimator independent of `φ_J`,
  and the fit-free bracket (`z_nr` crossing 4 between 0.570–0.595 and 0.620–0.645) remains
  the reading that needs no target at all.
