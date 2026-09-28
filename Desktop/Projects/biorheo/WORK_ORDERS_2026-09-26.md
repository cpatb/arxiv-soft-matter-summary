# Work orders 2026-09-26 — biorheo, dst_fit, polydisperse_jamming (in silico only)

**Dispatch sheet.** Same contract as `WORK_ORDERS_2026-09-08.md`: one package per agent, `owns` sets never intersect, pre-register before running, never loosen a bar, every prose number script-printed, tests pin every result, a failed gate is a recorded result. `NEXT_STEPS.md` "Open work" items 29–40 point here in order; `SESSION_PROMPT.md` is unchanged and still works.

**Scope decision that this sheet inherits (founder, 2026-09-26, `biofluids_venture\09_INSILICO_PLAN_2026-09-26.md`):** no bench. Nothing below asks for a measurement. Anything that needs one is an *input with a sensitivity*, not a default.

**Still open from 2026-09-08 and dispatched first:** B1 (device budget), B3 (anchor audit), B4 (omitted terms), B5 (referee 24–26), B2 (licensable protein curve), B6 (re-grade OFFER). Their briefs are in the 09-08 sheet and are not restated; only their order changes (see waves).

**Model choice**: Opus for judgement / anything client-facing or solver-facing; Sonnet for sharp, rule-checked execution.

---

# Part 1 — biorheo (`Desktop\Projects\biorheo`)

## N1 — Packaging, CI, and the one-command deliverable
**Sonnet.** ~2–3 h. No physics.

The page a client receives is produced by `case_studies/deliverable_example/run.py`, which a client-shaped input cannot reach. Make the chain one command.

1. `src/biorheo/cli.py` + `[project.scripts] biorheo = "biorheo.cli:main"`. Subcommands:
   - `biorheo fit CURVE.csv --family herschel_bulkley --out record.yaml` — columns `shear_rate_1_s, viscosity_Pa_s` (or `shear_stress_Pa`), optional `replicate`; writes a record YAML **and** its identifiability sidecar via `identify.write_identifiability` (the only writer).
   - `biorheo report record.yaml --needle-id 0.21e-3 --needle-len 12e-3 --barrel-id 2.3e-3 --flow 1.67e-9 --friction 0 --out page.md` — calls `apps.report.build` + `render`; refuses a gauge number (existing rule).
   - `biorheo screen --eta 0.0123 --eta-sd 0.0009 --licence ... <same device args>` — `build_from_screen`.
   - `biorheo export record.yaml --solver openfoam --rho 1050 --out transportProperties` (Fluent/COMSOL added by N3).
2. `.github/workflows/ci.yml`: `pip install -e ../dst_fit -e .[dev]`, `pytest -q -m "not slow"`, plus every `derivations/*.py` (each exits non-zero on a failed identity). dst_fit's own suite runs in a second job.
3. `python -m build` produces a wheel; a test installs it into a temp venv and runs `biorheo report` on `data/dmem_sints2024` to reproduce `deliverable_example/01_clean.md` byte-for-byte.

**Pre-register**: the CSV loader must refuse (a) a column with no unit in its name, (b) a rate below the record's envelope lower bound after fitting, (c) fewer than `n_params + 2` distinct rates. Gate N1.1: the three refusals fire with a named error. N1.2: byte-identical page. N1.3: CI green on the first push of a private GitHub repo (`cpatb/biorheo`, private; licence file deferred to founder D9 — do not add one).

**owns**: `src/biorheo/cli.py`, `pyproject.toml`, `.github/`, `tests/test_cli.py`, `tests/test_wheel.py`.

---

## N2 — L5: the credibility manifest (ASME V&V 40) as code
**Opus.** ~6–8 h over two nights. **The in-silico product that the 09-26 plan puts first (P2).**

Nothing exists. Build `src/biorheo/credibility/`:

- `factors.py`: the V&V 40 objects — `QuestionOfInterest`, `ContextOfUse`, `ModelRisk(influence, consequence)`, and the credibility factors: verification (code: SQA, numerical code verification; calculation: discretisation error, numerical solver error, use error), validation (model form, model inputs, test conditions, assessment — with UQ on both sides), applicability (relevance of QOI, relevance of validation activities to COU). Each factor has a **gradation** (V&V 40's 0–4-style levels, named not numbered) and an `Evidence` list.
- `manifest.py`: `Manifest.from_study(...)` collects evidence from a study directory: the identifiability sidecar → *model inputs* evidence; `solver/openfoam.export_report` round trip → *numerical code verification*; a case study's `results.json` with intervals → *validation assessment*; the extrapolation flag → *applicability*. **Every factor defaults to `MISSING`**; nothing is satisfied by default; the default risk is the highest level.
- `render.py`: a report section that prints the factor table, the achieved vs required gradation, and a *missing evidence* list.

Wrap two studies retroactively, each as a script + committed `manifest.json` + `manifest.md`:
1. `case_studies/mab_injectability/` — expected result: model-inputs evidence PRESENT (measurement repeatability, not a parameter interval — the manifest must carry the screen/fitted distinction from item 26), validation MISSING (no bench, by decision), applicability FLAGGED (0.74 decades).
2. The bidet grid study (`Desktop\Projects\bidet\`, read-only) — design-set unit 1.5: the manifest reproduces the observed orders of accuracy from the stored fields and flags the out-of-band design-model evaluation as an applicability gap.

**Pre-register**: G-N2.1 a study with no evidence renders every factor MISSING and risk HIGH. G-N2.2 the mAb manifest never uses the words "fitted" or "validated". G-N2.3 the bidet manifest's observed orders match `bidet/cfd_fit_summary.json` to printed precision. G-N2.4 (predicted FAIL, record it): the bidet study reaches the "required" gradation for a medium-risk COU — it will not, and the gap list is the deliverable.

Cite V&V 40-2018 and FDA's Nov-2023 guidance by title in docstrings; **no gradation wording copied verbatim from the standard** (it is copyrighted) — paraphrase and say so.

**owns**: `src/biorheo/credibility/`, `case_studies/*/manifest.*`, `tests/test_credibility*.py`. Read-only on `bidet/`.

---

## N3 — Exporters for the solvers clients own: Fluent and COMSOL, with round trips
**Opus.** ~5–6 h. Same contract as `solver/openfoam.py`.

- `solver/fluent.py`: for `newtonian, power_law, carreau_yasuda (Fluent: Carreau with a=2 only — refuse a≠2 with a named error), cross, herschel_bulkley (Fluent's HB has its own low-rate branch: name it, compute its induced error like `regularise.py` does), casson (Fluent lacks it: emit a `DEFINE_PROPERTY` UDF in C from the SymPy expression)`. Outputs: a `.jou`/scheme snippet setting the viscosity panel, or the UDF source. `parse_*` back.
- `solver/comsol.py`: a piecewise-cubic **interpolation table** `eta(gammadot)` over the measured envelope only (a table cannot extrapolate silently — write the envelope bounds as the table bounds and `extrap = none`), plus the analytic coefficient block for Carreau/Cross/HB in COMSOL's naming (`mu_inf, mu_0, lambda, n`, `mu_p, tau_y` with its own `gamma_dot_ref` regularisation — name it and compute the induced error).
- Round trip: 1-D plane Couette and Poiseuille reference solutions already exist in `derivations/`; the harness compares the **parsed-back** law against `record.eta_at` over the envelope (no solver run needed for Fluent/COMSOL — the point is that the dictionary is exact and its regularisation error is quantified; the OpenFOAM item already proved the solve is where it breaks).

**Pre-register**: N3.1 parse(export(record)) reproduces `eta_at` to 1e-12 rel for every exportable family. N3.2 the HB regularisation error curve for Fluent and COMSOL matches `regularise.induced_error` to 1e-6 at 50 rates. N3.3 a record whose sidecar says `degenerate` is refused without `allow_unidentified=True`. N3.4 the COMSOL table's bounds equal the record envelope exactly.

**owns**: `src/biorheo/solver/fluent.py`, `solver/comsol.py`, `tests/test_fluent*.py`, `tests/test_comsol*.py`, one `data/solver_export_matrix.md` (which family → which solver → exact / regularised / UDF / refused).

---

## N4 — The terms that default to zero become inputs with sensitivities
**Opus.** ~3–4 h. Depends on B1 and B4 having reported.

`apps/injectability.injection_force` carries `friction=0.0`; barrel drop, entrance loss and the kinetic term are absent (item 28, B4). Add `DeviceTerms(friction_N, friction_sd_N, barrel_len, barrel_id, entrance_coeff, include_kinetic)` to `apps/report.DeviceSpec`; **none defaulted** except by an explicit `DeviceTerms.hydrodynamic_only()` constructor whose name appears on the page. The page then prints: each term's contribution, the log-sensitivity of the total to each, the separability loss from `friction_sd` (item 28's `s*` recomputed per cohort — this is B1's budget function, reused not rewritten), and a one-line verdict "which term dominates at this operating point".

Temperature: `constitutive/shift.py` already carries `a`; add an **Arrhenius input** `E_a_J_mol` (user-supplied, never looked up) → `a(T)`, and print the superposition residual beside it as item 21 does.

**Pre-register**: N4.1 with all terms zero every existing page reproduces byte-for-byte (a test diff). N4.2 with the Funke 2016 glide bracket (doi:10.1002/jps.23945, 1.3–1.5 N) the mAb-vs-blood ratio prints 2.07–1.95× and the difference stays 1.70 N. **[SOURCE MIS-ATTRIBUTED — see Addendum A1 at the foot of this file, 2026-09-28; the number stands, the citation does not.]** N4.3 sensitivities sum to 1 in log space at every operating point tested (a homogeneity check that can fail).

**owns**: `src/biorheo/apps/injectability.py` (additive), `apps/report.py` (after B1 lands — coordinate: run B1 first), `constitutive/shift.py`, `tests/test_device_terms.py`.

---

## N5 — Registry to 25 licensable records via DataCite
**Sonnet.** Repeated 1–2 h packages, **one at a time** (exclusive on records, as B2).

Source list: `data/candidate_search_2026-09-08.md` (107 licensable DOIs). Priority order: any curve reaching ≥1e4 1/s; mucus and mucus simulants; IDDSI thickeners; synovial/HA; then anything CC BY with >8 rates. Per record: the B2 rules (no digitising, fit table not caption, licence from deposit metadata, grep `RECORD_IDS` consumers before appending), and a founder-review budget line: **≤1 founder-hour per record** (R3) — print the agent's own wall time in the record's `.md`.

Gate per package: N5.x every new DOI resolves offline-cached; envelope equals the deposit's rate range; `write_identifiability` sidecar present; fast suite green with the count updated in `NEXT_STEPS.md`.

**owns**: `data/<deposit>/`, `registry/records/`, `scripts/records_from_<deposit>.py`, `tests/test_<deposit>_records.py`.

---

## N6 — Viscoelastic family (L1) — **DO NOT DISPATCH** until a mucus/synovial/saliva record exists from N5
Placeholder. The 09-26 plan (C12) builds this only on demand.

---

# Part 2 — dst_fit (`Desktop\Projects\dst_fit`)

## D1 — Retire gate G4c with the reason recorded (founder decision D5, taken 2026-09-26)
**Sonnet.** ~1 h.

G4c (cross-curve consistency of `phi_m`) fails 1 of 6 on PMMA and the package's own result is that `phi_m` is not cross-curve consistent in any material (N4, N6, G4, the gauge-covariance finding). Mark it `RETIRED` in `scripts/gate_g4.py` (skip with a reason string, exit 0), keep the number in `README.md` "Status" verbatim, add a line to `NEXT.md`. `biorheo/tests/test_gate_0_1_dstfit_unchanged.py` pins dst_fit's suite — update its expected output in the same commit and say so in the message.

**owns**: `dst_fit/scripts/gate_g4.py`, `dst_fit/README.md`, `dst_fit/NEXT.md`, `biorheo/tests/test_gate_0_1_dstfit_unchanged.py`.

---

## D2 — `dst-fit` 0.1.x to TestPyPI, then a `biorheo` extra
**Sonnet.** ~1 h. **After D1.** Publish to **TestPyPI only** until founder D9 (licence) is decided; the repo has a `LICENSE` file already — check what it says and report, do not change it.

**owns**: `dst_fit/pyproject.toml`, a `dst_fit/PUBLISH.md`.

---

# Part 3 — P3: dense-suspension jamming risk (`Desktop\Projects\polydisperse_jamming` → new `biorheo.apps.jamming`)

## J1 — Freeze the φ_J relation as an importable function with a validity domain
**Opus.** ~3–4 h. Read-only on the DEM outputs.

`ANSWER_SHAPE_SIZE_2026-09-25.md` gives φ_J = Z/(Z + 2√3·w) with a population filling factor `w` (1.00 mono frictionless, 0.99 at μ 0.3; lognormal CV 0.5 → 0.76; doublet → 1.37) and a POST HOC accuracy of ±0.014, bias +0.006, at one pressure rung (p̂ 3e-5). Turn it into `polydisperse_jamming/phij_model.py`:

- `filling_factor(cv, shape, mu)` interpolated **only inside the trained set** (list every trained (cv, shape, mu, p̂) triple in a table in the module; outside → `OutsideTrainedDomain`).
- `phi_j(z, w)`, `chi(z_m, z_0, w)`, `x(phi, phi_0)`, and `verdict(phi, cv, shape, mu) -> "thickens_only" | "shear_jams" | "jammed"` with the 0.007-margin caveat at μ 0.3 carried as a flag, not swallowed.
- `predict_from_psd(d, weights)`: CV and skew from a discrete PSD (`check_discrete_psd.py` exists — reuse).

**Pre-register**: J1.1 held-out DEM packs (choose them by seed BEFORE fitting; list the seeds) recover φ_J within ±0.014 on ≥ 90 %; J1.2 the verdict reproduces the DNS `cvmap` result (CV 0 jams, CV 0.5 cannot, at fixed weighed φ); J1.3 (predicted FAIL) the relation holds across pressure rungs — it will not (w drifts with p̂); record the drift and pin the rung.

**owns**: `polydisperse_jamming/phij_model.py`, `polydisperse_jamming/test_phij_model.py`, `polydisperse_jamming/PHIJ_MODEL_2026-09-26.md`.

---

## J2 — `biorheo.apps.jamming`: the formulation-risk page
**Opus.** ~4 h. **After J1 and N1.**

Same shape as the injectability page. Inputs: PSD (or CV + shape class), weighed φ, particle density, a flow curve if any. Outputs on one page: φ₀, φ_J, χ, x, the verdict, σ* with its profile interval from `dst_fit` (gate G10's route) when a curve is supplied, and a **needle criterion**: at the needle wall stress from `apps.injectability`, is `x > χ`? plus `d_needle / d_90` against a stated bridging threshold that is **tagged UNVERIFIED with the literature it needs** (do not invent the number).

**Pre-register**: J2.1 a Newtonian record produces a page that says "no thickening physics applies" and refuses σ*; J2.2 the DMEM/Carbopol records (no particles) are refused; J2.3 a synthetic monodisperse hard-sphere input at φ = 0.55 returns `shear_jams` and at φ = 0.66 `jammed`.

**owns**: `src/biorheo/apps/jamming.py`, `tests/test_jamming.py`, `case_studies/jamming_example/`.

---

## J3 — Preprint draft from `ANSWER_SHAPE_SIZE_2026-09-25.md` + J1
**Opus.** ~3 h. Draft only, `polydisperse_jamming/preprint/`. Every number via `*_generated.tex` macros with the `\MISSINGNUM` fallback; `manuscript_ready.py` must pass. **Not posted** — founder D6/D7 first.

**owns**: `polydisperse_jamming/preprint/`.

---


---

# Part 4 — Phase 3 lab (post-PhD form of the company): the high-shear rig

Decided 2026-09-26: post-PhD the company is a specialized rheology computation + measurement lab. Pre-PhD nothing physical is bought; these two packages make the rig **costed** and its **analysis pipeline built and tested on open data** so that day one of Phase 3 is a purchase, not a design.

## L1 — Capillary / needle rheometer: bill of materials from vendor list prices
**Sonnet.** ~3 h. Research only; no purchase.

Replace the card's `~$3–8k UNVERIFIED` with a BOM whose every line carries a vendor URL and a list price (or `QUOTE REQUIRED` with the vendor named). Three configurations: (a) syringe-pump + pressure transducer + interchangeable capillaries/needles (the wedge rig); (b) (a) plus a load-cell force stand; (c) a commercial benchtop capillary rheometer for comparison (name the models). For each: rate range reachable, wall-stress range, sample volume per point, and the one measurement each cannot make. Deliverable `biofluids_venture/RIG_BOM_2026-XX.md`; numbers in a `.json` beside it. **Do not invent a price**; a price with no URL is `UNVERIFIED`.

**owns**: `biofluids_venture/RIG_BOM_*`.

## L2 — The rig's analysis pipeline, built now on open data
**Opus.** ~5 h. `src/biorheo/apps/capillary.py`.

From (pressure drop, flow rate, radius, length) series: apparent wall shear rate and stress; Weissenberg–Rabinowitsch correction; Bagley end correction from ≥2 lengths; Mooney slip from ≥2 radii **with the monotonicity-in-1/R check that item 27 (G27.6) showed can fail**; an uncertainty budget propagated from transducer, pump and geometry tolerances (user-supplied, no defaults); and the hand-off into `identify.fit_rate_curve` so the output is a record with a sidecar. Positive control: the Koponen deposit (`data/foam_koponen2023/`, 3 radii, imaged slip) — the pipeline must reproduce the anchor's reduced curve and its compressibility factor.

**Pre-register**: L2.1 Newtonian synthetic data returns n = 1 and zero slip to 1e-10; L2.2 Koponen reduction reproduces the anchor's numbers (median compressibility factor 1.574; 0 rows with slip above mean velocity after correction); L2.3 (predicted FAIL, record it) the Mooney construction on Koponen returns non-negative slip on every point — it will not, and the pipeline must refuse rather than report; L2.4 a series with one radius and one length returns `NOT_CORRECTABLE` flags, never a corrected curve.

**owns**: `src/biorheo/apps/capillary.py`, `tests/test_capillary.py`, `data/capillary_koponen_control.md`.

# Dispatch waves

| wave | packages | notes |
|---|---|---|
| **1 (now)** | B3 (Sonnet), B5 (Opus), N1 (Sonnet), D1 (Sonnet), J1 (Opus) | no file collisions; N1 and D1 both touch `test_gate_0_1_dstfit_unchanged.py` expectations → run D1 first, then N1 |
| **2** | B1 (Opus), B4 (Sonnet), N2 (Opus), B2 (Opus, exclusive on records) | N2 reads case studies only; B1 owns `report.py` |
| **3** | N4 (Opus, after B1+B4), N3 (Opus), D2 (Sonnet), J2 (Opus, after J1+N1) | N3 and N4 are disjoint |
| **4** | B6 (Sonnet, last), N5 (Sonnet, repeated, one at a time), J3 (Opus) | N5 never runs concurrently with B2 |
| **hold** | N6 | until N5 lands a viscoelastic biofluid record |

**Session sizing**: one package per night. A package that will not fit a night is split at its numbered gates, never stretched.

---

# Addenda (corrections to pre-registrations above; the originals are not rewritten)

## A1 — N4.2 cites the wrong source for its 1.3–1.5 N glide bracket (found by B3, 2026-09-28)

N4.2 as written attributes "1.3–1.5 N glide" to Funke et al. 2016 (doi:10.1002/jps.23945).
That is not where the number comes from. `case_studies/device_friction/results.json`
(`glide_note`) and `data/injection_force_terms_lit_2026-09-08.md` §1 both record the glide
bracket as **1.28–1.45 N from a *BioProcess International* trade article with no DOI,
carried as UNVERIFIED**; Funke 2016 publishes **breakloose** (~3 N silicone-oil-free,
~9 N siliconised), which `results.json` carries separately as `breakloose 3.0-9.0 N` with
its own ratio 1.18–1.52×. Full row: `data/external_anchor_audit_2026-09-08.md` row 7.

**What N4's owner must do before running N4.2**, choosing one and saying which:

- (a) keep the bar's number (1.3–1.5 N → 2.07–1.95×, difference 1.70 N) and re-cite it as
  "trade-article glide bracket, no DOI, UNVERIFIED (device_friction `glide_note`)"; or
- (b) re-register N4.2 on Funke's breakloose 3–9 N, whose predicted ratio is the 1.18–1.52×
  already printed on the device-friction page, and say the bar changed BEFORE any run.

Neither option moves a number that has been computed; the defect is the label on the
source. A bar graded later against "Funke 2016" would otherwise be graded against 3–9 N
and fail for a reason with no physics in it.
