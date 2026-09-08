# Item 22 pre-registration — the registry's first commercially licensable records

**Committed before any generator code.** Bars are not moved. A bar that turns out to be the
wrong bar is recorded with its number and a *stricter* test is added beside it. A
pre-registration is never rewritten; corrections go in an addendum at the foot.

## 0. What this item is, and why it is worth the disruption

Every one of the 115 records in this registry is `commercial_use: false` (item 18). That
does not block a paid study on a client's own curve, but it blocks every public use of this
work: a demo page, a webinar figure, a short course, a conference slide. All four worked
deliverable pages are stamped NOT LICENSABLE.

`data/dmem_sints2024/` is a **CC BY 4.0** rheometer dataset of DMEM supplemented with FBS
(Zenodo 13828992, `10.5281/zenodo.13828992`; paper arXiv:2409.13650). Item 20 downloaded
it, scanned it, and read the paper. `SOURCE.md` records what the paper settles. This item
turns it into records.

**It cannot carry an injectability claim and no page built on it may make one.** The
measured rates stop at 1 1/s (10 1/s on one interval) against 1e4–1e5 in a fine needle:
four to five decades of extrapolation, against 1.42 for the blood records.

## 1. What is being written

**Family: `power_law` only.** That is the form the authors themselves fit (Ostwald–de Waele,
their Eq. 4), and the scan already shows it beats a Bingham alternative on 7 of 8
intervals. Fitting the registry's other families here would be inventing a description the
source does not make and would need its own justification.

**Eight records, one per interval**, not four per file. Each `.txt` holds more than one
`Meas. Pts.` table and an importer reading the first one discards half the data — the same
shape as item 7's concatenated ramps. The intervals are separate measurements at different
rate ranges and are not replicates: the two 10 % FBS 20 °C intervals run 0.0714–1 and
0.01–1 1/s, and the two 37 °C ones run 0.01–1 and 0.01–10.

Each record carries:

- `temperature_C`: **20 or 37**, from the paper, not guessed. Results §*Medium & its
  viscosity* and Methods §*Viscosity measurements* both state `T = 20 °C` for "room
  temperature".
- `geometry`: cone-plate, Anton Paar MCR 502, CP50-2, `d = 0.208 mm`, Peltier temperature
  control — from the deposited header and the paper's Methods.
- **the pre-shear protocol in `notes`**, because a record inherits it: *2000 s of rotation
  at 1 1/s before measurement, 160 s per point*, and the authors' stated position that the
  plateau rather than the startup value is the relevant one. A later reader must not take
  these for as-loaded viscosities.
- `licence`: CC BY 4.0, `commercial_use: true`. **These will be the first such records.**

## 2. The exclusion rule, stated before the fit

One point is dropped and only one: `DMEM_and_10pct_FBS.txt` interval 1 opens at
`gd = 2.24e-7` 1/s with an apparent viscosity of 8710 Pa s and a torque of 0.0641 µN·m.
That is the instrument's zero-rate point, not a measurement. **The rule is `gd > 1e-3` 1/s**,
it is the same rule `dmem_candidate_scan.py` already uses, and the dropped count is printed
by the generator and asserted by a test. No other point is excluded for any reason.

## 3. Gates

| gate | statement | bar | prediction |
|---|---|---|---|
| **G22.1** | the fitted flow index reproduces the paper's own published values | `n` within **0.05** of the paper's for the six intervals it names (0.32, 0.31, 0.31, 0.5/0.54); the two it does not name are reported, not barred | PASS on 6/6 |
| **G22.2** | the generator's own round trip: write, read the record BACK, and require its envelope to contain both extremes of its source curve | 8/8, exit non-zero otherwise | PASS (the existing `round_out` gate, applied to a new source) |
| **G22.3** | the DOI resolves | `10.5281/zenodo.13828992` in `doi_cache.json`, resolved with `allow_network=True` | PASS |
| **G22.4** | the records are the FIRST licensable ones | `commercial_use` true on 8/8 new, false on all 115 existing | PASS |
| **G22.5** | a deliverable page built on one is NOT stamped NOT LICENSABLE, and IS stamped extrapolated at four-plus decades | `non_commercial` absent, `extrapolated` present, `decades_past` > 4 | PASS |
| **G22.6** | **the two temperatures do NOT collapse under a single shift factor** | `collapse_shift` rms log residual on the 10 % FBS 20 °C vs 37 °C pair is **above 0.05**, i.e. an order of magnitude above the 1.4e-05 a real shift gives | PASS, and see §3.1 |

### 3.1 G22.6 is the one worth writing down in advance

A one-parameter time–temperature shift rescales viscosity and rate together and **cannot
change the flow index**: under `eta -> a eta(a gd)` a power law maps `k -> a^n k` with `n`
fixed (item 21's derivation). The paper's own indices are **0.32 at 20 °C and 0.5 at
37 °C** for the same 10 % FBS medium. So the two curves differ in a parameter the shift
cannot touch, and **the collapse must fail**.

This is the first time item 21's residual machinery meets real data, and the prediction is
that it reports a *large* residual. If it instead reports a small one, either the shift
implementation or the fitted indices are wrong, and finding out which would be the item.
**A small residual here would be a failure, not a success** — which is the opposite of how
a "does it work" test usually reads, and is why it is pre-registered.

## 4. The consumers that must move, enumerated BEFORE the change

This repo has been bitten by exactly this: item 11 appended two records to a shared list and
silently re-decided two closed gates, because the suite asserted verdicts that stayed true
*on a different question*. Every place that depends on the registry's contents is listed
here, found by grep, and each must be handled deliberately.

| file | what it asserts | how it must be handled |
|---|---|---|
| `tests/test_offer_claims.py` | `len(recs) == 115`; all `commercial_use: false`; the OFFER's "All 115 records" sentence | **The claim CHANGES.** The OFFER's §6 item 1 must be rewritten: the blocker is lifted for a demo. Both must move together. |
| `tests/test_design_space.py` | `g["n"] == 115`; `n_nc == len(rows) == 115` | Item 18's artifact reads its OWN frozen record list, not `load_all()`, so its numbers do not move. The `commercial_use` assertion must become "all 115 of THAT list", not "all records". |
| `tests/test_openfoam_export_nonyield.py` | `len(recs) == 115` | A CLOSED item-6 gate. It must read a frozen list or be restated as `>= 115`, and which was chosen must be recorded. |
| `tests/test_openfoam_export_null.py` | `n == 115`; the `.md` string "115 records, of which 20 `newtonian` and 20 `power_law`" | A CLOSED item-11 gate. **Its G11.5 separations are computed over the registry**, so new `power_law` records would enter a closed gate's sample. It must be frozen. |
| `tests/test_report.py` | every record triggers `non_commercial` | **This is the test that should FAIL**, and its failing is the result. Restate it as: every record from the pre-2026-09-08 registry is non-commercial, and the new ones are not. |
| `tests/test_profiling_wrapper.py` | `>= 115` | already written as an inequality; nothing to do. |
| `case_studies/*/run.py` | select by model name or an explicit donor list | verify each selects by an explicit list or a family filter that excludes `dmem`, and record which. |

**Rule for this item: no existing gate may change its VALUE. Any gate whose sample would
grow gets a frozen list instead**, and the frozen list is committed in the same change as
the records. A gate that silently re-decides is the failure this table exists to prevent.

## 5. What this item does NOT do

- No injectability number, no force, no needle. §0.
- No temperature MODEL. G22.6 tests that a shift does *not* describe this pair; it does not
  fit one.
- No other family. `power_law` only, because that is what the source fits.
- No change to `dst_fit`, and no change to any existing record.
- **No claim that this fluid is representative of anything.** It is cell-culture medium,
  chosen because it is licensable and biological, not because a buyer asked about it.

---

## Addendum A, 2026-09-08 — this pre-registration named the wrong published numbers

**Not edited above, per the repo's standing rule.** Three reviewer passes on the shipped
item found that §3's reference indices `0.32, 0.31, 0.54, 0.50` are **figure-caption
values used as COMSOL inputs**, rounded and shared across conditions, and that the paper
also carries **Table 2, "Fit results for `K` and `n` ... for the experimental rheometry
data shown in Fig. 1"**, which resolves the four deposited conditions one row each and
disagrees with three of the four values graded against here. Table 2 also publishes `K`,
<!-- CORRECTED: it is FOUR of four; see dmem_anchor_correction_prereg.md addendum A.4. -->
which this pre-registration never asked anyone to compare.

Three further defects in the grading built on §3, all confirmed by rerun:

- G22.1's pass expression tests `len(graded) == 6` while 7 intervals are graded, so it
  could not have passed whatever the fits returned;
- G22.1' pools the 37 C interval that §3 itself excludes as out of the paper's rate range,
  and that interval supplies 20 of the 25 pooled points;
- the "per-interval spread reaches 0.076" argument is printed as 0.1044, a max taken over a
  condition with only one graded interval.

The correction, its gates and the restated headline are in
**`data/dmem_anchor_correction_prereg.md`** (item 23). `G221_TOL` is not moved; the
corrected anchor FAILS it on one of four conditions and the failure is reported.
