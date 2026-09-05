# 04 — Asset map: what already exists, and where it lands in the stack

Every entry was opened, not read from its README alone. Readiness scores are the July 2026
sweep's where one exists; "state" is what the code does today. Layers refer to
`03_ARCHITECTURE.md`.

## 1. Direct contributions

| project | state today | lands in | what carries over | what must change |
|---|---|---|---|---|
| `dst_fit` | Packaged (`pyproject`, `src/`, `tests/`), profile-likelihood intervals, identifiability verdicts on real deposited data (Edinburgh DataShare), a strain-clock measurement with a named definition dependence, and a permanently failing gate G4c kept on purpose | **L2 identify** (the seed of the whole layer) | the verdict vocabulary (measurable / bounded / degenerate), the profile-likelihood machinery, the steady-subset lesson, the "a fail is data" gate discipline | generalise the model interface from Wyart–Cates to all of L1; the `si_fit` withdrawn preset must not migrate; the G4c decision (retire with reason, or ship failing) is the founder's |
| `tier0_wyart_cates_closure` | Tested closure with reversibility gate; PyTorch autodiff calibrator recovers synthetic truth; `calibration/` documents what a flow curve can and cannot hand a solver | **L1 constitutive** (suspension family), **L3** (the calibration handoff argument) | model, state kinetics, the calibration write-up as the first L3 design note | move autodiff from PyTorch to JAX/diffrax for consistency with L2; parameters are illustrative until the registry supplies real ones |
| `bio_fluids_bridge` | Runs, five gates pass; Moore–Coussot structural closure for blood, thixotropic loop | **L1** (structural family), **L0** (first blood record once calibrated) | the lambda-kinetics model and the silica-to-blood transferability argument | calibrate to published viscometry with provenance; the readiness note is right that alone it is textbook, and its value is as one leg of the closure story |
| `thrombosis_jamming` | Runs, signature validation passes; occlusion threshold and divergence; the discriminating test against Ku-group data is named but not run | **L6 thrombosis**, flagship paper | the closure, the falsifiable critical-point prediction, the prior-art positioning | run the named test on digitised published occlusion data; resolve the arXiv 2602.20886 precedent flagged as "verify"; fix gamma_c and the modulus scale from independent literature before calling anything a prediction |
| `uroflow_splashback` | v2, 44 tests, Weber long-wave breakup, Mundo/Cossali splash thresholds with an explicit normal, wall-vs-pool result | **L6 urology**, **L4** jet rig analysis, consumer tier case study | the whole physics chain; the wall-versus-pool result is the first public case study | the calibration measurement (breakup length from a phone clip) has not been done; splash index magnitudes are ranking-only until it is |
| `bloodstain_bpa` | 6 + 17 tests; ballistic origin reconstruction; the 3.5x error-amplifier finding; bias attribution | **L1/L6** drop-impact spread model for a bodily fluid; forensic vertical if ever wanted | spread-factor inversion, uncertainty propagation pattern, the amplifier analysis | prior art (Attinger group) unchecked and possibly fatal to the novelty claim; as a company asset it is a demonstrator of the UQ style, not a product |
| `bidet` | Journal manuscript; four-grid convergence with iterative convergence fixed first; OpenFOAM impinging-jet cases; `gen_case.py` as the single source of truth for case names | **L3** (OpenFOAM tooling), **L5** (the worked example of calculation verification done right), consumer-tier flagship | the grid-study discipline including the lesson that the fine grid was the unconverged one; the case-generation pattern; the manuscript as a sales asset | the design model evaluates its fitted closure outside its own fit band (r_max/d* = 23–26 against [2, 20]); that must be fixed or bounded before it is shown to a client |
| `carbopol` | HB fits to published flow curves; bead-drop yield-number bracket with the three-convention trap documented | **L0** (surrogate records), **L1** (yield-stress regularisations), **L4** (mucus simulant) | the fits, the convention warning (Beris 0.143 is the middle of three definitions), the literature file | the L0 record must carry which yield-number convention was used |
| `polydisperse_jamming` | DEM packing and impact study; rattler vs backbone result; the confinement-force and reduced-units bugs found and fixed; finite-size flag split by family | **L1** (polydispersity correction to phi_J), particle rung of the fidelity ladder | the packing-to-phi_J pipeline; the energy-admissibility gates | the dynamic verdicts still flip between E = 1 and E = 10 at CV = 0.5; only the static observables are scale-free, and the registry must say so |
| `paper1`, `axidrop`, `Drop Impact — Simulation` | Basilisk solvers (`wc_drop_impact.c`, `hb_drop_on_bath.c`), the regularisation sweeps, dropout-flag lesson, the `EPS_MIN` and `MUMAX` exposure work | **L3 basilisk exporter**, **L4** drop rig, free-surface rung | the solvers, the `getenv`-parameter convention, the regularisation-sweep protocol as a product feature | `axidrop` is still scaffold; the drop-impact solvers are research code and need the pluggable constitutive interface `axidrop` was designed for |
| `impact_stress` | Wagner and shear-stress toolkit, pore-pressure relaxation, `conventions.py`, `materials.py` | **L1/L6** (impact-loading of suspensions and gels), the stress-gauge convention lesson | the toolkit and the explicit convention file | `materials.py` ships the withdrawn `si_fit`; strip before reuse |
| `rehab_smartfluid` | Layers 2–3 built: activation number, actuator bandwidth closed form, six tests; field survey identifying the physics-to-device gap | **L1** (transient thickening), a device vertical if pursued (protective and orthotic devices are medical devices) | the activation-number result and the literature map of transient constitutive models | pending decision on whether this is a vertical or a paper |
| `datacenter_cooling` | ROM plus 2-D CFD calibration; the Lewis-number sweep that exposed a diffusion artefact as a headline | **fidelity ladder** pattern: ROM calibrated to CFD, with a transport-coefficient sweep as the mechanism test | the ROM-calibration pattern and the sweep-a-transport-coefficient lesson | domain is unrelated; only the pattern transfers |
| `hydrogels` | Validated measurement layer (36/36) with no data; calibration from fiducial; indentation modulus; crack metrics | **L4** imaging pipelines, **L6 absorbents** (swelling kinetics), wound-care and depot-gel verticals later | calibration module, swelling kinetics (D, Biot), the review discipline that caught five silent bugs | needs footage; unrelated to fluids until a gel vertical opens |
| `glymphatic_transport` | Péclet regime map reproducing the field synthesis | **L6** scoping template; a CSF record in L0 | the dimensional-analysis scoping format, which is exactly the shape of a paid feasibility study | none; it is a template, not a product |
| `slime mold`, `catfish school`, `frisbee` | PTV, collective motion, urban CFD with a full verification pipeline | **L4** PTV/PIV code; the frisbee verification and retraction discipline | image-analysis pipelines and the test-from-repo-root lesson | none of the domains transfer |

## 2. Infrastructure that transfers wholesale

- `~/.claude/rules/`, the reviewer agents (`physics-referee`, `units-auditor`,
  `science-code-reviewer`, `numerics-verifier`, `citation-checker`, `calibration-auditor`) and
  the skills `derive-and-check`, `convergence-study`, `lit-verify`, `dimensional-audit`,
  `paper-figure`. These are the company's quality system in prototype form. A V&V 40 package
  is what they produce when pointed at a client study, and that should be said out loud in
  the sales material once it has been demonstrated once.
- The generated-macro pattern (`*_generated.tex` with `\MISSINGNUM`) for reports.
- The Status-card and ledger discipline for project state.
- `lessons.md`: a large fraction of it is directly about the failure modes the credibility
  layer must catch mechanically (regularisation readouts, quantised observables, unconverged
  grids, collinear predictors, pooled-versus-within-cell correlations).

## 3. What does not exist yet, in order of how much it blocks

1. **L0 registry.** Nothing. Schema, loader, provenance checks, first records.
2. **L3 exporters and round trip.** Nothing productised; the Basilisk `getenv` convention and
   the bidet OpenFOAM cases are the raw material.
3. **L5 credibility engine.** Nothing as code; the discipline exists in skills and reviewers.
4. **L4 rigs as build sheets.** The imaging pipelines exist; the rigs are lab equipment
   used ad hoc and undocumented as products.
5. **L1 viscoelastic family.** Nothing; needed for mucus, synovial, saliva, thickeners.
6. **`libbiorheo` for OpenFOAM.** Nothing; the bidet cases use built-in Newtonian.
7. **Poroelastic rung.** Nothing; needed for W1 depot back-pressure and absorbents.

## 4. What must be retired before it is shown to anyone

- `si_fit` in any path (`impact_stress/materials.py`, any preset in the Basilisk solvers).
- The bidet design model's out-of-band evaluation, until bounded.
- Uroflow splash-index magnitudes; only rankings until the calibration clip is analysed.
- `bloodstain_bpa` novelty language, until the Attinger prior art is checked.
- Any registry-style number in the existing READMEs marked "illustrative" or "verify".
