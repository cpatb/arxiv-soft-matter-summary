# biorheo

Cited fluid records, identifiable constitutive models, and their handoff to solvers,
for fluids that are not water. Phase 0 of the biofluids venture
(`../biofluids_venture/`, design spec and plan under `docs/superpowers/`).

The package exists to make three questions answerable from a flow curve before any
simulation is run:

1. Which constitutive parameters does this data actually determine? (`identify`)
2. Where did the numbers come from, and over what range were they measured? (`registry`)
3. What does a device calculation rest on when it leaves that range? (`apps`)

## Layout

```
src/biorheo/
  units.py            one pint registry; boundaries convert once, everything inside is SI
  constitutive/       RateModel (a dst_fit Model in shear-rate space) and the families:
                      gnf.py (newtonian, power_law, herschel_bulkley, casson,
                      carreau_yasuda, cross), structural.py (Moore-Coussot kinetics),
                      regularise.py (papanastasiou / bercovier_engelman / clamp, with the
                      knob named in the model and its induced error reported)
  identify/           fit_rate_curve -> IdentifiedModel: profile-likelihood intervals and a
                      per-parameter verdict (measurable / bounded_only / not_identifiable /
                      degenerate, the last from |corr| > 0.99 or a Fisher null direction);
                      synthetic.py recovery studies; record.py writes the sidecar
  registry/           YAML fluid records (parameters with units, measured envelope, DOI
                      sources, licence); identifiability lives only in a sidecar the
                      identify layer writes
  apps/               injectability.py: Herschel-Bulkley needle flow (closed form for any
                      n, quadrature for any monotone model), plunger force with profile
                      and covariance intervals, log-sensitivities, and an extrapolation
                      flag against the fitted envelope
derivations/          SymPy scripts that regenerate their .md; each exits non-zero on a
                      failed identity or hand spot
scripts/              build the Carbopol records from the ETD2020 deposit and identify them
case_studies/         one script per study, writing results.md and a figure beside it
tests/                fast suite; `-m "not slow"` skips re-running dst_fit's own suite
```

## Install and test

```
pip install -e ../dst_fit -e .[dev]
python -m pytest -q tests/ -m "not slow"
python derivations/gnf_limits.py && python derivations/structural_kinetics.py && python derivations/buckingham_reiner.py
```

`dst_fit` is consumed unchanged; `tests/test_gate_0_1_dstfit_unchanged.py` runs its
suite to prove it.

## What the data currently supports

Fifteen Carbopol records (doi:10.25835/hjn2qcvy, ETD2020, smooth plates, 20 C) with
measured envelopes of [1e-3, 200] 1/s and every Herschel-Bulkley parameter graded
`measurable`; 13 of the 15 are flagged `residual_structure.structured`, so the HB form
is a systematic misfit at rmse 0.2-17 % even where every parameter is determined. The
one case study so far, needle injectability of two of those gels, puts the needle wall
2.9 decades above the top of that envelope and finds the force ordering of the two gels
reverses between rheometer and needle; see
`case_studies/injectability_surrogates/results.md`.

Twenty human whole-blood records (Horner 2020 thesis workbook via Farrington,
Jariwala, Armstrong 2023, doi:10.1007/s00397-023-01402-2; 37 C, EDTA, 16 rates in
[0.1, 700] 1/s, one donor each) as Casson fits with `commercial_use: false` (author
posted, no licence stated). The authors' Casson table is reproduced to 1e-8 before
anything is built on it, and it is a convention: an unweighted line in sqrt space,
whose yield stress is 0.4-71 % (median 33 %) above biorheo's log-viscosity fit of the
same 16 points, because Casson misfits every donor but one (3 runs of residual sign
against 9 expected, under either objective). Cross and Carreau-Yasuda fit the same curves 5x better; the table
is `data/horner2020/model_comparison.md`. The workbook's step-shear and LAOS sheets are
unused transient data for the structural family.

## Conventions

- Models return viscosity `eta(gd)` in Pa s; `fit_rate_curve` takes `eta = sigma/gd`.
- Parameters are SI floats internally. A record stores `value` and `unit`, never `si`.
- A flow-curve fit reports its envelope; anything computed outside it is flagged
  `extrapolated` and the flag is carried, not dropped.
- Every derivation has a hand-computed anchor that does not pass through the code
  under test.
