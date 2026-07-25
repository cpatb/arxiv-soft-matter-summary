# FINISH BIDET — drive to a CFD-validated journal paper

**Paste into Claude Code with `C:\Users\compl\Desktop\Projects\bidet` as the working directory.**
Environment: Windows + WSL2 with **OpenFOAM v2606** installed (`source /usr/lib/openfoam/openfoam2606/etc/bashrc`), Python 3 (numpy/scipy/matplotlib), MiKTeX (`pdflatex`). Use your engineering/code-review skills and git as you go.

## Where this stands (read first)
The paper was reframed after a validated recompute (see `model/`, `bidet_optimization_1_CORRECTION_DOSSIER_2026-07-23.md`). Current state: `bidet_optimization_1.CORRECTED.tex` is a correct, honest **scaling/framework** paper. Its central result — the laminar Watson r⁻⁵ closure is inapplicable; a turbulent closure reverses the design conclusions — is sound, but the **turbulent closure prefactor C_P and exponent n_t are calibrated to a single commercial datum, not computed.** The reviewers' one hard blocker for a journal is exactly this: the "reversal" is asserted from a one-point fit, not demonstrated. **This brief removes that blocker by computing the turbulent wall-shear law with CFD, then turns every estimated number into a computed one.**

Target venue: **European Journal of Physics** (pedagogical/methodological) is the honest home; with the CFD in hand, *Physics of Fluids* becomes defensible. Decide at the end based on what the numbers support.

## The validated model you're extending
`model/bidet_model.py` (8/8 tests in `model/test_bidet_model.py`) already implements the Watson laminar closure, a general power-law closure `tau_w = C_tau·rho·V0²·Re_d^p·(d/r)^n`, the cleaning-criterion inversion, and the lever calculations. **You will replace the calibrated turbulent (C_P, n_t) with CFD-computed values, then re-run everything.** Do not rewrite the model; extend it.

---

## STAGE 1 — CFD: compute the turbulent wall-shear law (the blocker)

Goal: an axisymmetric OpenFOAM simulation of a round jet impinging normally on a flat plate, giving τ_w(r) at the bidet operating point, from which you fit the closure `(τ_w/ρV0²)·Re_d^{1/2} = C_P·(r/d)^{-n_t}` (the Phares 2000 form) — i.e. pin **C_P and n_t** directly instead of calibrating to commercial flow.

1. **Build the case** in `cfd_openfoam/` (create it). Axisymmetric wedge (`wedge` patches, ~5° wedge, 1 cell thick) is the right geometry for a round jet on a plate — far cheaper than 3D and standard for this problem.
   - Solver: **`simpleFoam`** (steady, incompressible — water; use `nu = 8.9e-7`) with **k-ω SST**. (A single liquid phase is fine; the free surface is not needed to get the wall shear under the impinging core — state this assumption.)
   - Geometry: nozzle diameter d = 1 mm; plate at standoff H; radial extent to r ≥ 20 mm; resolve the wall with y⁺ ≈ 1 (low-Re SST) or wall-functions y⁺ 30–100 — pick one and be consistent. Inlet: uniform V0 from Re_d.
   - Cases: **Re_d ∈ {3e3, 1e4, 3e4}** × **H/d ∈ {1, 6, 20, 60}** (start with H/d = 6 and Re_d = 1e4 to shake the case out, then sweep).
2. **Version note:** the case dicts must be v2606-correct (ESI). `simpleFoam`, `kOmegaSST`, `snappyHexMesh`/`blockMesh` are all present. If a dict keyword errors, fix it (usual suspects at a new version) and record the fix.
3. **Extract τ_w(r)** on the plate patch (`wallShearStress` function object → integrate/sample the radial profile). Write `tau_w_r.csv` per case.
4. **Fit the closure:** for each case fit `n_t` (the radial decay exponent) and `C_P` (the prefactor in the Phares form). Report how `n_t` and `C_P` vary with Re_d and H/d. Expect `n_t ≈ 1–1.3` in the developed turbulent wall jet (vs the laminar 5). **Validate** the case first: (a) reproduce the known free-jet centerline decay Vc/V0 = B·d/H (B≈6) for H/d≳6; (b) check grid independence (refine once, τ_w within ~10%); (c) sanity-check τ_w magnitude against the Phares 2000 correlation the paper already cites.
5. Save a figure `fig_cfd_tau_w.png`: CFD τ_w(r) vs the fitted power law vs the laminar r⁻⁵, at the nominal case.

**Compute reality:** an axisymmetric wedge case is small (~10⁴–10⁵ cells) and each run is minutes, not hours — the whole sweep is feasible on the workstation in an afternoon. If any case won't converge, drop to the single nominal case (Re_d=1e4, H/d=6); one computed (C_P, n_t) still removes the blocker.

## STAGE 2 — feed CFD back into the model and recompute

1. In `model/`, replace the commercial-point calibration of the turbulent closure with the **CFD-fitted (C_P, n_t)** (and their Re_d/H dependence if you captured it). Keep `recirc`-style backwards-compat: leave the old calibration available behind a flag; make CFD the default.
2. Re-run `run_bidet.py`: recompute d*, Q*, and **every lever** (pressure exponent, standoff sign+magnitude, oblique, oscillation) with the computed closure. The lever *directions* should hold (that's the robust result); the *magnitudes* are now computed, not estimated.
3. Extend `test_bidet_model.py` with a test that the CFD-fitted closure gives the physical standoff sign and a pressure exponent in the CFD-derived range. All tests must pass.
4. Regenerate `fig_Qstar_vs_n.png` and the closure figure with the CFD point marked.

## STAGE 3 — upgrade the manuscript from "framework" to "CFD-validated"

Edit a **new copy** `bidet_optimization_1.JOURNAL.tex` (don't overwrite CORRECTED). Changes:
1. Add a **CFD methodology + results subsection** (the axisymmetric jet case, validation, the fitted C_P/n_t, `fig_cfd_tau_w.png`). This is the new evidence.
2. Replace every "≈ΔP⁻² (exponent-sensitive)" / "order 10× (estimate)" hedge with the **computed** value + its CFD-based uncertainty. The reversal is now *demonstrated*, not asserted.
3. Keep the honest scoping where it's still true (τ_c is still an assumed detachment shear — see Stage 4), but drop the "calibration circularity" caveat: C_P now comes from CFD, independent of the commercial datum.
4. Update the abstract, discussion, and conclusions to reflect a validated result. Re-evaluate the venue claim (EJP vs PoF) in the conclusions.
5. **Verify the two flagged citations** (`Bhagat & Wilson 2016`, `Wilson et al. 2014`) — confirm author list / volume / pages (search the web / the DOIs) and remove the `% FLAG FOR AUTHOR: verify` once done.

## STAGE 4 — τ_c, the remaining honest gap
The cleaning shear τ_c is an *assumed* detachment stress (the cited Woolley data are bulk fecal rheology, a different quantity). You cannot compute this from CFD. Do the honest thing: (a) search for an actual wall-detachment/adhesion-shear measurement for soft biological soils and cite it if found; (b) otherwise keep τ_c explicitly labeled an assumption, and lean on the fact that the *design levers are ratios* largely insensitive to τ_c (Q* ∝ τ_c^{1/3}). State this clearly. Leave a `% FLAG FOR AUTHOR` if a real source can't be found.

## STAGE 5 — review and finalize
1. Run the multi-tier review harness on the journal draft: `Read ../CLAUDE_CODE_PREPRINT_REVIEW_PIPELINE.md and follow it, TARGET = bidet/bidet_optimization_1.JOURNAL.tex`. It runs validation → postdoc → professor subagents and writes a dossier. Resolve what it raises.
2. **Compile:** `pdflatex` twice (inline bibliography, no BibTeX). Confirm it builds, figures resolve, no undefined refs, no `[PLACEHOLDER]`/`[CITE`/`a^{?}` survive. Sane page count.
3. Commit at each stage (`git add -A && git commit`) so the CFD, the model update, and the manuscript changes are separable.

## Definition of done
- `cfd_openfoam/` runs and produces a fitted (C_P, n_t) with a validation check; `fig_cfd_tau_w.png` exists.
- `model/` uses the CFD closure by default; all tests pass; levers recomputed.
- `bidet_optimization_1.JOURNAL.tex` compiles, presents the CFD as the evidence, and survived the review pipeline; citations verified; τ_c honestly handled.
- A short `BIDET_JOURNAL_STATUS.md` at the top level: what the CFD found, which numbers changed vs the estimates, the recommended venue, and any remaining author decision.

## Guardrails
- Don't fabricate CFD results — if a case didn't converge, say so and use what did.
- Preserve `CORRECTED.tex` and `model/` originals; work on copies/extensions.
- No claim that isn't backed by a passing test or a converged, validated CFD case.
- If the CFD *contradicts* the paper's qualitative story (e.g. a wrong-sign lever), report it prominently — that's a finding, not a failure.
