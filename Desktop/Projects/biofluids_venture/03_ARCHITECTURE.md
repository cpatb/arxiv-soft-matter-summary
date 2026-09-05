# 03 — Technical architecture

The company's product is a chain (characterise, identify, carry, predict, validate, certify).
The software is that chain as a set of packages with hard boundaries, so that any single link
can be sold, tested, or replaced on its own. This document fixes the layers, their contracts,
and the simulation capability each one needs. Package names are placeholders until D1 in
`07_DECISIONS.md` is settled; `biorheo` is used below.

## 0. Design rules, inherited from the last two years of mistakes

1. **A number without provenance does not enter the system.** Every parameter carries its
   source (DOI, figure, digitisation method, or a measurement manifest), the conditions it was
   measured under, and the range of shear rate, stress, frequency and time the data actually
   covered. Evaluating a model outside that envelope raises, never warns.
2. **Every gate is code.** Acceptance criteria are written before the run and enforced by a
   test that exits non-zero. A failed gate is stored as a result with the same status as a pass.
3. **Round trip or it did not happen.** A solver export is accepted only when the target solver,
   given the export, reproduces the rheometer measurement within a stated tolerance.
4. **Regularisation is a declared parameter.** Any clamp, floor, Papanastasiou exponent,
   bi-viscosity ratio or minimum gap is a named input, swept in every study that quotes a
   magnitude, and reported beside the result. An observable whose threshold is itself the
   regularisation parameter is rejected at design time (`yfrac-is-a-clamp-readout`).
5. **Fidelity ladder, lowest rung first.** Closed form, then 1-D, then reduced-order, then
   full CFD or particles. A higher rung is used only when the lower one cannot answer the
   question of interest, and the lower rung stays in the report as the sanity anchor.
6. **Units at the boundary.** `pint` on every public function argument that carries a
   dimension. SI internally. Camera and rheometer conversions happen once, in a manifest.
7. **Anchors that do not pass through the code under test.** Each physical result has at
   least one hand-computed point, analytic limit, published number, or discriminating negative.
8. **Reproducible by one command.** Every figure, table and report regenerates from a script;
   byte-identical regeneration is the refactor test.

## 1. The layers

```
L6  apps/          injection · airway · thrombosis · urology · dysphagia · absorbents
L5  credibility/   V&V 40 workflow engine, report generator, evidence manifest
L4  bench/         rig designs, imaging pipelines, calibration audit, benchtop UQ
L3  solvers/       export to OpenFOAM/Basilisk/COMSOL/Fluent/SimVascular; round-trip
L2  identify/      fitting, identifiability, model selection, protocol design
L1  constitutive/  model forms, state evolution, regularisations, limits
L0  registry/      cited, versioned fluid property records + schema
     ─────────── cross-cutting: units · provenance · gates · reporting · seeds ───────────
```

Dependencies point downward only. `apps` may import everything; `registry` imports nothing
above it. Each layer is a separately installable package with its own tests and its own
acceptance gate.

### L0 `biorheo.registry` — the fluid property registry

The moat asset. A git-versioned collection of records, one per (fluid, condition, source):

```yaml
fluid: human_whole_blood
condition: {T_C: 37, hematocrit: 0.45, anticoagulant: EDTA, donor: pooled, pathology: none}
protocol: {geometry: cone_plate, gap_um: 50, mode: steady_shear, gammadot_range_1_s: [0.1, 1000]}
model: carreau_yasuda
params: {eta0_Pa_s: ..., eta_inf_Pa_s: ..., lambda_s: ..., n: ..., a: ...}
uncertainty: {type: profile_likelihood_95, ...}
identifiability: {eta0: bounded_only, eta_inf: measurable, lambda: measurable, n: measurable, a: degenerate_with_n}
envelope: {gammadot_1_s: [0.1, 1000], time_s: null, frequency_Hz: null}
provenance: {doi: 10.xxxx/..., figure: "3b", extraction: "WebPlotDigitizer 4.6, 2026-09-xx", verifier: "lit-verify"}
notes: "..."
```

Rules: no record without a resolved DOI or a measurement manifest; every record has an
envelope; the identifiability field is filled by L2, not typed by hand; records are immutable
and superseded, never edited, so a report can cite a commit hash.

Contents at launch (target): whole blood and plasma across hematocrit and temperature; a
mucus simulant family and human sputum from the 2025 multimode BMP paper; synovial fluid;
saliva; urine (Newtonian, but with the surface-tension and contamination data the splash work
needs); CSF; vitreous; representative high-concentration mAb formulations from published
viscosity-concentration curves; IDDSI thickener families; and the lab's own surrogates
(Carbopol, xanthan, glycerol, PEG, silica suspensions) with in-house measurements.

Business shape: a core subset published open (CC-BY) with the open benchmark suite; the
curated, identifiability-graded, envelope-checked tier is licensed.

### L1 `biorheo.constitutive` — model forms

One interface for every model:

```
eta(gammadot, state, T)        apparent viscosity, generalised Newtonian
stress(D, state, T)            full tensor where the model has one
dstate_dt(state, D, sigma, T)  structure kinetics; empty for time-independent models
limits()                       envelope and regularisation the model requires
regularise(scheme, **knobs)    returns a regularised variant with the knobs recorded
```

Families to implement, in order of use:

- *Generalised Newtonian:* power law, Carreau, Carreau–Yasuda, Cross, Herschel–Bulkley,
  Casson, Quemada, Bingham. With every regularisation in circulation for the yield-stress
  members (Papanastasiou, bi-viscosity, `mu_max` clamp, Bercovier–Engelman) implemented
  side by side, because the export layer must be able to emit whichever the target uses.
- *Structural / thixotropic:* first-order structure kinetics (Moore, Houska), the
  lambda-closure already in `bio_fluids_bridge`, multimode BMP for mucus, and the
  Wyart–Cates frictional-contact closure with the strain clock from `tier0_wyart_cates_closure`
  and `dst_fit`.
- *Suspension:* Maron–Pierce and Krieger–Dougherty with a migrating jamming point;
  the polydispersity correction from `polydisperse_jamming` where the particle-size
  distribution is known.
- *Viscoelastic:* Oldroyd-B, Giesekus, PTT, FENE-P, multimode Maxwell, for mucus, synovial
  fluid, saliva, vitreous and the extensional behaviour of dysphagia thickeners. These are the
  hardest to carry into a solver and the ones where the round-trip rule earns its keep.
- *Damage and activation kernels:* the power-law hemolysis kernel with its calibration
  provenance made explicit, and the VWF activation function from `thrombosis_jamming`.
  Both are shipped as "relative, not absolute" and the report says so.

Every model ships with its SymPy definition, analytic limits, and a numeric spot-check at
three parameter points, per `derive-and-check`.

### L2 `biorheo.identify` — fitting and identifiability

Generalisation of `dst_fit` from dense suspensions to all of L1. Inputs are rheometer records
in a common schema (steady shear, amplitude and frequency sweeps, step shear and shear
reversal, creep, capillary and slit, CaBER/DoS extensional). Outputs:

- parameters with confidence intervals from profile likelihood, not Wald;
- the parameter correlation matrix and its condition number;
- a verdict per parameter: *measurable*, *bounded only* (hit a bound, quoted as an interval),
  *degenerate* (enters only as a combination with another), *not identifiable from this
  protocol*;
- model selection across families with AIC/BIC and a residual-structure check, because a
  better AIC with structured residuals is the wrong model winning;
- **protocol design:** given the degeneracies found, the additional measurement that would
  break them (a reversal transient for the strain clock; a second temperature for
  Arrhenius parameters; an extensional test for the viscoelastic branch). This is a sellable
  output on its own.

Built-in corrections, each with a test: cone truncation and edge effects, wall slip
(gap-dependence test), Taylor and inertial instability limits, evaporation during the
measurement (the drop-impact project's own control-parameter drift), instrument compliance,
and low-torque limits.

Implementation: `scipy.optimize` for the base fits; JAX with `diffrax` for models with state
evolution so that gradients of a transient fit are exact rather than finite-differenced.
This is the inverse-problem tooling the landscape memo identified as "published, validated,
never shipped" (RhINNs), done with adjoints rather than a PINN because adjoints are easier
to validate.

### L3 `biorheo.solvers` — the measurement-to-solver bridge

For each target, an exporter and a round-trip harness.

| target | what is emitted | round-trip geometry |
|---|---|---|
| OpenFOAM (current release) | `transportProperties` / `momentumTransport` entries for built-in GNF models; a maintained `libbiorheo` for the models OpenFOAM lacks (structural, multimode, migrating-jamming) | cone-plate or Couette case, steady sweep |
| Basilisk | header with the viscosity closure and its regularisation knobs as `getenv`-read parameters, following the existing drop-impact solvers | the same Couette case in axisymmetric form |
| COMSOL | material-model expressions and a parameter table; later a user-defined model | Couette |
| Ansys Fluent | UDF source for `DEFINE_PROPERTY` viscosity with regularisation | Couette |
| SimVascular | Carreau–Yasuda parameter block, with a warning when the fit was not Carreau–Yasuda | Poiseuille |

The `libbiorheo` OpenFOAM library is the one piece of solver-side code the company owns.
It is the modern replacement for the gap rheoTool left when it fell behind, scoped to the
fluids the company sells, and it is what makes the OpenFOAM path fully open and scriptable.

The round-trip harness runs the target solver on the rheometer geometry, extracts the flow
curve, and compares it with the measurement and with the L1 evaluation. Three outcomes are
reported: agreement within tolerance; disagreement traced to the target's regularisation
(the induced error is then quoted with the export); disagreement not traced (the export is
rejected). This is where "every solver regularises Herschel–Bulkley differently" becomes a
number in a report instead of a war story.

### L4 `biorheo.bench` — the simplified experimental model layer

The company's answer to "how do you know the simulation is right." Reference rigs, each with a
build sheet, a calibration manifest, an analysis pipeline and a stated uncertainty budget.

| rig | measures | serves | already have |
|---|---|---|---|
| syringe–needle force rig (load cell, syringe pump or drop-weight, high-speed camera on the needle exit) | injection force vs rate, jet formation, drip | W1 | camera, calibration audit, drop-impact tracking |
| capillary and slit viscometer | high-shear viscosity beyond rheometer range, wall slip | W1, W2 | Buckingham–Reiner in closed form |
| jet and spray rig (nozzle, pressure vessel, backlit imaging, PDPA optional) | breakup length, drop size distribution, plume geometry, splash on wet and dry targets | W2, W4 | `uroflow_splashback`, `bidet`, drop-impact rigs |
| drop and splash rig | spread, splash threshold, retention on textured or wetted surfaces | W4, consumer tier | `paper1` pipeline, `bloodstain_bpa` spread model |
| flow loop with transparent test section (PIV/PTV) | velocity fields in a device surrogate, pressure drop | W3, catheters, dialysers | PTV code from hydrogels and slime-mold work; glymphatic PTV plan |
| absorbent and porous test cell | imbibition rate, front shape, saturation | consumer tier | Washburn/Richards closed forms to be written |
| sedimentation and settling column | suspension stability, phi profile | W1 depots, W5 | `polydisperse_jamming` |

Pipelines reuse the calibration-auditor and units-auditor discipline; every image-derived
quantity carries px/mm and fps provenance, and camera-native observables (ratios of pixel
counts and frame counts) are preferred where a calibration constant could be the confound.

### L5 `biorheo.credibility` — the V&V 40 engine

Objects mirror the standard: question of interest; context of use; model risk as model
influence times decision consequence; credibility goals per factor (code verification,
calculation verification, validation for computational model and comparator, applicability).
Evidence is attached from the pipeline, not typed:

- calculation verification: grid and timestep ladders with fitted observed order and Richardson
  extrapolation from the `convergence-study` skill; the ladder is refused if any rung failed its
  own iterative criterion, because the bidet study showed the fine grid was the unconverged one;
- validation: benchtop comparison with uncertainty propagated on both sides, using paired
  differences rather than differenced variances (frisbee lesson);
- applicability: the registry envelope against the device's operating range; extrapolation is
  listed as a gap.

Output: a report (LaTeX via the lab's existing generated-macro pattern, with `\MISSINGNUM`
rendering visibly) and a machine-readable manifest listing every credibility factor with its
evidence or its absence. The manifest is what a reviewer or a client's regulatory consultant
reads first.

### L6 `biorheo.apps` — vertical packages

Thin. Each composes L0–L5 for one class of question and ships with one open case study.

- `injection`: needle flow for GNF fluids (closed form → axisymmetric CFD), injection-force
  prediction with rheological UQ, shear exposure and residence-time histograms for
  aggregation risk, depot back-pressure coupling to a poroelastic tissue model
  (FEniCSx or an OpenFOAM poro-solver; the published porohyperelastic work is the benchmark).
- `airway`: nozzle spray and plume, deposition on a nasal cast, mucus-layer clearance under
  cough-like transients with a multimode BMP mucus.
- `thrombosis`: the jamming closure in a stenosis and in a benchmark device, occlusion time vs
  shear rate, the relative-only hemolysis kernel, and the FDA nozzle as the Newtonian anchor.
- `urology`: uroflow splashback (done), catheter and urethral flow, bidet and urinal splash.
- `dysphagia`: the IDDSI syringe flow test as a forward model for GNF and viscoelastic
  thickeners, with extensional viscosity as the parameter the test cannot see.
- `absorbents`: imbibition into layered porous media for menstrual and incontinence products.

## 2. Simulation capability required, by fidelity rung

| rung | tool | what it answers | validation anchor |
|---|---|---|---|
| closed form | SymPy-derived, `derive-and-check` | Poiseuille for power-law and HB (Buckingham–Reiner), jet breakup (Weber long-wave), Taylor–Aris dispersion, Washburn imbibition, lubrication, Windkessel | textbook limits, hand-computed points |
| 1-D / 0-D transient | `scipy.integrate`, `diffrax` | structure kinetics under a prescribed history, injection with a compliant device, syringe flow test, occlusion in a throat | analytic step responses at limiting parameters |
| reduced-order | ROM + calibration to CFD (the `datacenter_cooling` pattern) | design sweeps in minutes for client optimisation | the CFD it was calibrated to, out of sample |
| continuum CFD, single phase | OpenFOAM + `libbiorheo` | device flows, deposition, shear-exposure fields, pressure drop | grid ladder with observed order; benchtop flow loop |
| free-surface CFD | Basilisk (existing solvers), OpenFOAM VoF as the second opinion | drops, jets, sprays, splash, filling, drip | the drop-impact benchmarks already validated in `axidrop` |
| particle / discrete | YADE-style DEM (`polydisperse_jamming`), LF-DEM for reference | suspension microstructure, jamming point vs polydispersity, depot settling | the packing benchmarks already run |
| poroelastic tissue | FEniCSx (or OpenFOAM poro-solver) | subcutaneous depot, back-pressure, absorbents | confined-compression data from the literature |
| differentiable inverse | JAX + `diffrax` | constitutive calibration from transients and from partial flow measurements | synthetic recovery with known truth, then real data |
| uncertainty | polynomial chaos or GP surrogate over the L2 posterior, Monte Carlo through the rung in use | the rheological uncertainty on the quantity of interest | convergence of the surrogate error against direct samples |

Solvers are not written from scratch. OpenFOAM and Basilisk are the workhorses; commercial
targets are served by export. The only owned solver code is `libbiorheo` and the Basilisk
headers, and both are material models, not solvers.

## 3. Repository layout

```
biofluids_venture/            this design; becomes the company handbook
biorheo/                      the monorepo (created at Phase 0, see 05_ROADMAP.md)
  packages/
    registry/                 L0   records/, schema/, loader, provenance checks
    constitutive/             L1
    identify/                 L2
    solvers/                  L3   openfoam/ (libbiorheo C++), basilisk/, comsol/, fluent/, simvascular/
    bench/                    L4   rigs/, pipelines/, calibration/
    credibility/              L5
    apps/                     L6
  benchmarks/                 the open suite: geometry, fluid record, reference solution, tolerance
  case_studies/               public write-ups, one per app, regenerated from scripts
  studies/                    client work, one directory per engagement, never in the public repo
  docs/
```

Python 3.12, `pint`, `numpy`, `scipy`, `sympy`, `jax`/`diffrax`, `pandas`, `matplotlib` with
the lab's house style. C++ only inside `libbiorheo`. Tests with `pytest`; gates are tests. CI
runs the round-trip harness against OpenFOAM in a container so the bridge cannot drift from
the solver the way the WSL launch scripts drifted from the repo.

## 4. Data contracts between layers

- `registry → constitutive`: a `FluidRecord` (the YAML above, validated by schema).
- `constitutive → identify`: a `Model` object exposing `eta`, `stress`, `dstate_dt`, `limits`.
- `identify → registry`: a `FitResult` with parameters, intervals, correlation, verdicts,
  residual diagnostics; writing it back creates a new record, never edits one.
- `constitutive → solvers`: `Model` plus a chosen `Regularisation`; the exporter returns an
  `Export` with the emitted files, the target version, and the round-trip result.
- `bench → credibility`: a `Measurement` with value, uncertainty budget, calibration manifest.
- `solvers/apps → credibility`: a `SimulationResult` with QoI, convergence ladder, and the
  regularisation sweep.
- `credibility → report`: an `EvidenceManifest`; the report is a view of the manifest.

Each contract is a dataclass with a schema, a version, and a test that round-trips it
through JSON.

## 5. What is deliberately not built

- No GUI. Scripts, notebooks-as-scripts, and generated reports. A web front end is a
  Phase 3 question at the earliest and only if a licensing product materialises.
- No PINN forward solver. Adjoint calibration covers the inverse problem; classical solvers
  win the forward one.
- No in-house Navier–Stokes solver. The company's code is material models and bridges.
- No patient data. Registry records are literature or in-house surrogate measurements;
  client data stays in `studies/` under the client's terms.
