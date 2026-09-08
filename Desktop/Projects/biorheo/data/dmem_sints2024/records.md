# DMEM + FBS: the registry's first commercially licensable records

> **CORRECTED 2026-09-08 (item 23). The reference values graded against on this page are the paper's FIGURE CAPTION numbers -- the flow indices it fed to COMSOL, rounded and shared across conditions -- not its Table 2, "Fit results for `K` and `n` values in (4) for the experimental rheometry data shown in Fig. 1", which resolves these four conditions one row each and disagrees with ALL FOUR of them (two are the right number against the wrong row). Re-graded against Table 2 the anchor FAILS the same 0.05 bar at 0.0599. The page below is left as it was written, with its numbers, and `data/dmem_sints2024/anchor_recheck.md` carries the correction and what survives it.**

Written by `scripts/records_from_dmem.py` from `data/dmem_record_prereg.md`. Source, licence and traps: `data/dmem_sints2024/SOURCE.md`. Every number here is printed by that script.

**Eight `power_law` records, CC BY 4.0, `commercial_use: true`.** They are the first records in this registry that may be used in a public demo, a webinar figure, a short course or a paid deliverable. Every other record is `commercial_use: false`.

**They cannot carry an injectability claim.** The measured rates stop at 1 1/s (10 1/s on one interval) against 1e4-1e5 in a fine needle: four to five decades of extrapolation, against 1.42 for the blood records.

## G22.1 — the fit reproduces the PAPER'S OWN published flow indices

The authors fit the same Ostwald-de Waele power law and publish their indices in the figure captions. Those numbers were published before this repo existed and are not computed by any code here, which is what makes this an anchor rather than a self-check.

**FAILED**, 7 intervals graded, worst |n - paper n| = 0.06511 against a 0.05 bar. 1 interval is reported without a bar: it is the only one running to 10 1/s, and the paper's indices are quoted for the 0.01-1 1/s range.

**G22.1's bar was MIS-SPECIFIED and has not been moved.** It compares each INTERVAL's index against a value the paper publishes once per CONDITION. Those are different objects, and this run's own numbers prove no fit could have satisfied it: the per-interval spread within a FULLY GRADED condition reaches 0.0753, larger than the 0.05 bar itself. Its number stays on this page. (Over every condition including the one with a single graded interval it is 0.1044, which is the number this page printed before 2026-09-08; a condition with one graded interval cannot illustrate a claim about BOTH of its intervals.)

## G22.1' — the same anchor, with the interval-choice excuse removed

One published index per condition means a fit to that condition, so both intervals of a condition are fitted TOGETHER and that is what is graded. This is stricter in kind rather than looser in size: it removes an excuse instead of widening a tolerance, and it is still graded against numbers this repo did not compute.

**PASS**, 4/4 conditions, worst |n - paper n| = 0.02986 against a 0.05 bar.

| condition | pts | pooled `n` | paper `n` | &#124;dev&#124; |
|---|---|---|---|---|
| 10 % FBS, 20 C | 29 | **0.3195** | 0.32 | 0.0005 |
| 5 % FBS, 20 C | 30 | **0.3310** | 0.31 | 0.0210 |
| 1 % FBS, 20 C | 30 | **0.5699** | 0.54 | 0.0299 |
| 10 % FBS, 37 C | 25 | **0.5120** | 0.5 | 0.0120 |

**biorheo's fitting layer reproduces an independent paper's published flow indices to within 0.03, from the raw deposited files, having never seen those indices.** On the 10 % FBS 20 C condition the agreement is 0.0005. That is the strongest external anchor this repo has for its identification layer: every other check compares biorheo against biorheo, a closed form, or a solver.

One defect in this grading code is recorded rather than quietly fixed: the first run keyed the pooling on the published index as well as the condition, which split the 37 C condition (whose second interval has no published index) and pooled only its 5-point half. That was a bug in the grader, not a finding about the data.

## Per-record fits

| record | FBS | T (C) | pts | dropped | rate range (1/s) | `k` | `n` | paper `n` | &#124;dev&#124; | verdicts |
|---|---|---|---|---|---|---|---|---|---|---|
| `dmem_sints2024_10pct_20C_a` | 10 % | 20 | 14 | 1 | 0.0714-1 | 0.02796 | **0.2945** | 0.32 | 0.0255 | k=measurable, n=measurable |
| `dmem_sints2024_10pct_20C_b` | 10 % | 20 | 15 | 0 | 0.01-1 | 0.02718 | **0.3134** | 0.32 | 0.00655 | k=measurable, n=measurable |
| `dmem_sints2024_5pct_20C_a` | 5 % | 20 | 15 | 0 | 0.01-1 | 0.04165 | **0.3687** | 0.31 | 0.0587 | k=measurable, n=measurable |
| `dmem_sints2024_5pct_20C_b` | 5 % | 20 | 15 | 0 | 0.01-1 | 0.03237 | **0.2934** | 0.31 | 0.0166 | k=measurable, n=measurable |
| `dmem_sints2024_1pct_20C_a` | 1 % | 20 | 15 | 0 | 0.01-1 | 0.01575 | **0.5394** | 0.54 | 0.000569 | k=measurable, n=measurable |
| `dmem_sints2024_1pct_20C_b` | 1 % | 20 | 15 | 0 | 0.01-1 | 0.0185 | **0.6003** | 0.54 | 0.0603 | k=measurable, n=measurable |
| `dmem_sints2024_10pct_37C_a` | 10 % | 37 | 5 | 0 | 0.01-1 | 0.005005 | **0.4349** | 0.5 | 0.0651 | k=measurable, n=measurable |
| `dmem_sints2024_10pct_37C_b` | 10 % | 37 | 20 | 0 | 0.01-10 | 0.004572 | **0.5393** | - | - | k=measurable, n=measurable |

## The protocol these records inherit

Each point is a **160 s** measurement taken after **2000 s of rotation at 1 1/s**. The authors found that up to 2000 s was needed for the viscosity to become stationary, and state that they take the plateau rather than the transient startup value to be the one relevant to flow in their device. **These are post-pre-shear plateau viscosities, not as-loaded ones**, and every record says so in its `notes`. A reader who forgets this is reading a different measurement from the one that was made.

## The one exclusion, and it is stated rather than silent

Points below **0.001 1/s** are dropped. Exactly one point in the whole dataset meets that: `DMEM_and_10pct_FBS.txt` interval 1 opens at 2.24e-7 1/s with an apparent viscosity of 8710 Pa s and a torque of 0.0641 uN m, which is the instrument's zero-rate point rather than a measurement. The dropped count is a field on every record.
