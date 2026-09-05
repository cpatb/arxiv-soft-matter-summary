# Design spec — biofluids venture foundation (2026-09-05)

Status: **awaiting founder review.** No product code has been written. This spec fixes what
the first implementation plan will build (Phase 0 of `../../../05_ROADMAP.md`) and points to
the longer documents for the reasoning.

## Problem

Simulation sold to medical-device, drug-product and personal-care companies treats bodily
fluids and viscous formulations as Newtonian or as textbook Carreau. The error is unquantified
and no vendor sells the material layer (measurement → identifiable model → solver → benchtop
validation → credibility package) as one chain. Landscape and evidence: `01_LANDSCAPE.md`.

## Product

A seven-layer stack (`03_ARCHITECTURE.md`): registry (L0), constitutive (L1), identify (L2),
solvers bridge (L3), bench (L4), credibility (L5), vertical apps (L6). Sold as consulting
studies first, licensed software and registry second. Entry wedge: injectability of viscous
drug products (W1); second: airway and mucus (W2); flagship paper: thrombosis as jamming (W3);
public case studies: urology and splash (W4). Reasoning: `02_NICHE_AND_POSITIONING.md`.

## Phase 0 scope (what the first plan builds)

1. Monorepo `biorheo/` with `packages/{registry,constitutive,identify}`; `dst_fit` migrated
   into `identify` behind a generalised `Model` interface with its tests and figures unchanged.
2. Registry schema (`FluidRecord`), loader, provenance validator refusing records without a
   resolved DOI or measurement manifest; envelope required.
3. Records: Carbopol surrogates (own data), Guy–Hermes–Poon PMMA, whole blood from two
   independent published sources. Identifiability field written by L2 output only.
4. L1 generalised-Newtonian family with Papanastasiou, bi-viscosity and clamp regularisations
   side by side; SymPy definitions, limits, three-point spot checks via `derive-and-check`.
5. L1 structural family (Moore–Coussot lambda closure; Wyart–Cates with strain clock); L2
   recovers synthetic truth on 20 seeds and returns "degenerate" on a deliberately degenerate
   protocol.
6. Case study 1: needle injectability of a yield-stress and a shear-thinning surrogate,
   closed form (Buckingham–Reiner) with UQ from the L2 posterior; hand-computed anchor at three
   parameter sets.
7. Ten discovery conversations (founder). Kill test R1.

Out of scope for Phase 0: any solver exporter, `libbiorheo`, the credibility engine, rigs,
viscoelastic models, a GUI, a company name.

## Data contracts

`FluidRecord`, `Model`, `FitResult` as dataclasses with schema, version, and JSON round-trip
tests (`03_ARCHITECTURE.md` §4).

## Testing and gates

Every unit has a gate written before the run and enforced by a `pytest` test that exits
non-zero. Gates are listed per unit in `05_ROADMAP.md`. Reviewers (`science-code-reviewer`,
`units-auditor`, `physics-referee`, `numerics-verifier`) run before any unit is called done.

## Error handling

Out-of-envelope evaluation raises. Missing provenance raises at load. A fit that hits a bound
returns an interval and the verdict `bounded_only`, never a point estimate. Regularisation
knobs are mandatory named arguments with no defaults in L3-facing code.

## Open decisions blocking implementation

D1 (package name, placeholder `biorheo` is unregistered on PyPI as of today), D2 (wedge),
D5 (G4c). All others can be deferred past Phase 0. See `07_DECISIONS.md`.

## Risks

`06_RISKS_AND_KILL_CRITERIA.md`. The one that decides Phase 1 is R1.
