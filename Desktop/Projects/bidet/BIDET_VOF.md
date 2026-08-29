# BIDET VOF — does jet breakup restore standoff as a strong lever?

**Paste into Claude Code with `C:\Users\compl\Desktop\Projects\bidet` as the working directory.**
Environment: Windows + WSL2 with **OpenFOAM v2606** (`source /usr/lib/openfoam/openfoam2606/etc/bashrc`),
Python 3 (numpy/scipy/matplotlib), MiKTeX (`pdflatex`). Use your engineering/code-review skills and git.

## Why this exists (read first)
The completed journal drive (`BIDET_JOURNAL_STATUS.md`) computed the turbulent wall-shear closure with a
**single-phase, steady, axisymmetric** `simpleFoam` case and found — correcting the paper's own analytic
claim — that **standoff is a weak lever** (τ_w(r_max) varied only ~1.14× over H/d 1–20). That result is firm
**only for a coherent jet**, because single-phase RANS structurally cannot see jet breakup. The reviews named
this the single highest-value follow-up: the weak-lever finding is extrapolated at the toilet-scale standoff
(H/d ≈ 60), exactly where a real water jet may have **broken into ligaments and droplets**, collapsing the
momentum flux delivered to the plate. If breakup collapses the wall shear at large H, **standoff becomes a
strong lever again in the real regime** — the opposite design guidance — and the crossover standoff is a
genuine, PoF-defensible result.

**This brief tests that hypothesis with two-phase VOF CFD.** It does NOT overwrite the single-phase result;
it asks where that result stops being valid.

## The operating point (already established)
d = 1 mm nozzle; at **Re_d = 10⁴**, V₀ = Re·ν/d ≈ **8.9 m/s**. Water: ρ = 1000 kg/m³, ν = 8.9e-7 m²/s,
σ = 0.072 N/m. Air: ρ = 1.2, ν = 1.5e-5. Dimensionless: **We = ρV₀²d/σ ≈ 1100**, **Oh ≈ 0.003** —
turbulent, inertia-vs-surface-tension breakup regime. Breakup length is uncertain (order 10–100 d); whether
it is shorter or longer than the H/d = 60 standoff is precisely the open question.

## The physics change vs the last run
- Solver: **`interIsoFoam`** (VOF water+air, geometric isoAdvector interface — sharper than `interFoam` for
  ligament/drop pinch-off; present in v2606). **Transient** (PIMPLE), NOT steady.
- Surface tension explicit (σ = 0.072). Courant-limited: `maxCo ≈ 0.2`, `maxAlphaCo ≈ 0.2`,
  `adjustTimeStep yes`. Integrate long enough (several flow-through times, t ≳ 5·H/V₀) for the wall-shear
  statistic to reach a quasi-steady time-average.
- Turbulence: **URANS k-ω SST** first (robust, cheap). Note in the writeup that **LES (WALE)** is the fidelity
  upgrade for breakup — state URANS-smears-breakup as a limitation, do not hide it.

---

## STAGE A — axisymmetric screen (cheap, decisive-enough), do this first
Goal: find out whether even an axisymmetric jet loses its coherent core (and its wall-shear delivery) before
H/d = 60. Axisymmetry captures the varicose (Rayleigh–Plateau) mode and gives a **conservative lower bound**
on breakup — real 3D breakup is at least this strong.

`cd cfd_openfoam` (extend the existing single-phase case tree; make a `vof/` subdir — do NOT overwrite the
`simpleFoam` cases).

1. **Geometry:** axisymmetric **wedge** (~5°, 1 cell thick), nozzle d = 1 mm on the axis, radial extent
   R ≥ 25 mm, axial length = the standoff H. Build cases at **H/d = 6, 20, 60** at Re_d = 10⁴.
   - **Mesh: static, graded — NOT AMR here.** 2D is cheap, and `dynamicRefineFvMesh` is finicky on wedges;
     avoid that trap. Put a fine band around the jet column: **≥ 20 cells across the 1 mm diameter**, graded
     coarser toward R. (Reserve AMR for the 3D Stage B.)
2. **Initial/boundary conditions:**
   - Nozzle inlet: `alpha.water = 1`, `U = 8.9 m/s` uniform (or a 1/7-power profile — state which).
   - Plate (the impingement wall): no-slip wall; this is where you sample τ_w.
   - Far/atmosphere patches: `alpha.water` `inletOutlet`, `U` `pressureInletOutletVelocity`, `p_rgh`
     `totalPressure` (p₀ = 0) — the standard open-boundary VOF trio so entrained air can come and go.
   - Domain initially air (`setFields`: `alpha.water = 0` everywhere, jet introduced through the inlet).
   - **v2606 gotcha carried over:** any `surfaceFieldValue` function object needs an explicit
     `writeFields false;` entry (this bit the data-center case). `pRefCell` is NOT needed here — the
     atmosphere patch fixes the pressure datum.
3. **Run** transient to quasi-steady wall shear. Watch `alpha.water` stays bounded in [0,1] (MULES) and the
   interface isn't parasitically shredding from spurious currents (a sign of under-resolution / too-large Co).

## STAGE A — validation (gate before trusting anything)
Two independent checks, both required:
1. **Coherent-limit cross-check:** at **H/d = 6** the core is intact, so the VOF time-averaged **τ_w(r) must
   agree with the single-phase Stage-1 result** (`sweep_results.jsonl` / `tau_w_r.csv` at H/d = 6,
   Re_d = 10⁴). If it doesn't, the VOF case is wrong — fix before proceeding. This is the hard gate.
2. **Breakup length:** measure the intact/breakup length L_b (axial distance where the connected liquid core
   first pinches) and compare to a turbulent-jet breakup correlation (e.g. Grant & Middleman, or a
   Weber-number turbulent-breakup fit). Order-of-magnitude agreement; record the correlation used.

## STAGE A — the deliverable figure
For each H/d, from the **time-averaged** field after the jet establishes:
- **Axial momentum flux delivered to the plate** (∫ ρ u_z² dA over the plate-adjacent region, water phase)
  **vs H/d.** Single-phase conserves this (hence the weak lever); breakup + air entrainment + droplet drag
  should make it decay. This is the headline.
- **τ_w(r) vs H/d**, and the resulting **cleaning radius r(τ_c) vs H/d** at τ_c = 1, 3, 5, 10 Pa (the band
  the journal paper already uses).
- Save `fig_vof_momentumflux_vs_H.png` and `fig_vof_tauw_vs_H.png`.

**Decision point:** if the delivered momentum flux (and hence τ_w at the cleaning radius) **collapses by
H/d = 60**, the weak-standoff-lever finding breaks in the real regime → standoff is a **strong** lever there,
and you have the crossover. If it does NOT collapse even with breakup, the weak-lever result is *strengthened*
(now robust to breakup in axisymmetry). **Either outcome is a publishable result — report what you find.**

## STAGE B — 3D confirmation (only if Stage A shows the crossover, and only at the crossover case)
Axisymmetry suppresses the azimuthal modes that dominate real breakup, so a positive Stage-A result must be
confirmed in 3D before it goes in a fluids journal.
1. Build **one** 3D case at the crossover standoff (likely H/d = 60), Re_d = 10⁴: a periodic **sector**
   (e.g. 30–90°) if the mean is near-axisymmetric, else full 3D.
2. **Here use AMR:** `dynamicRefineFvMesh` refining on the `alpha.water` interface (`field alpha.water`,
   `lowerRefineLevel 0.01`, `upperRefineLevel 0.99`, 2–3 refinement levels) to keep the cell count sane —
   otherwise ~20 cells across a 1 mm jet over a 60 mm domain in 3D is millions of static cells.
3. Confirm the Stage-A momentum-flux collapse survives 3D (expect it to be **at least as strong**, since 3D
   breaks the jet faster). Report the 3D delivered-momentum-flux number against the axisymmetric one.
4. If Stage B is out of compute reach, **stop after Stage A and say so** — the axisymmetric result stands as
   a conservative lower bound, honestly labeled. Do not fake the 3D case.

## STAGE C — feed back into the model and manuscript
1. In `model/`, add a **breakup-aware standoff correction**: a delivered-momentum-flux factor f(H/d) from the
   VOF (behind a flag, single-phase remains the default coherent-core closure). Recompute the standoff lever
   with it. Add a test that the VOF-corrected lever has the sign/magnitude the CFD showed. All tests pass.
2. Edit a **new copy** `bidet_optimization_1.VOF.tex` (do NOT overwrite `.JOURNAL`). Add a two-phase CFD
   subsection: the VOF method, the two validations, the momentum-flux-vs-H figure, and the crossover finding.
   **Rewrite the standoff-lever conclusion** to state the coherent-jet result AND its breakdown at large
   standoff — the honest two-regime picture. Update abstract/discussion/conclusions and re-evaluate the venue
   claim (if a clean crossover is demonstrated in 3D, PoF/PRF becomes arguable; if Stage A only, keep EJP/AJP
   and scope 3D as future work).

## STAGE D — review and finalize
1. Run the multi-tier review harness on the new draft:
   `Read ../CLAUDE_CODE_PREPRINT_REVIEW_PIPELINE.md and follow it, TARGET = bidet/bidet_optimization_1.VOF.tex`.
   Resolve what validation → postdoc → professor raise.
2. `pdflatex` ×2 (inline bib). Confirm it builds, figures resolve, no undefined refs, no placeholders.
3. Write `BIDET_VOF_STATUS.md`: what the VOF found (breakup length, whether the momentum flux collapsed and
   by how much, the crossover standoff if any), which numbers/conclusions changed vs the single-phase paper,
   Stage-A-only vs Stage-B-confirmed, the revised venue call, and remaining author decisions.
4. Commit per stage (case → validation → 3D → model → manuscript → review), scoped to `bidet/` only.

## Compute reality
- Stage A (2D axisymmetric, static mesh, 3 standoffs): transient VOF, hours per case on a workstation — an
  afternoon total. This is the screen; do it fully.
- Stage B (one 3D case + AMR): a day or more. Only if Stage A earns it. Check `nproc`; run in parallel
  (`decomposePar` → `mpirun -np N interIsoFoam -parallel`).
- If a case won't stay bounded, the usual fixes: lower `maxAlphaCo`, refine the interface band, check the
  `cAlpha`/isoAdvector settings, reduce inlet turbulence intensity. Record the fix.

## Honest limitations (state all four plainly in the paper)
1. VOF resolves **primary** breakup (column → ligaments → large drops), **not** fine atomization — droplets
   below cell size are unresolved; claim only the primary-breakup momentum budget.
2. **URANS smears turbulent breakup**; LES/WALE is the fidelity upgrade (scope as future work).
3. **Axisymmetry (Stage A) under-predicts breakup** — it is a conservative lower bound; the real jet breaks
   sooner. Stage B is what removes this caveat.
4. The real toilet **air environment** (confinement, cross-drafts, the actual receiving surface) is not
   modeled — this is a canonical jet in still air.

## Guardrails
- Don't fabricate CFD results — if breakup doesn't collapse the flux, that is the finding; report it.
- Preserve `bidet_optimization_1.JOURNAL.tex`, the single-phase `cfd_openfoam/` cases, and `model/` originals;
  work on copies/extensions (`vof/`, `.VOF.tex`).
- The H/d = 6 coherent-limit cross-check against the single-phase τ_w is the hard validation gate — no
  downstream number is trustworthy until it passes.
- Every headline number traces to a bounded, time-averaged, validated run. State Stage-A-only vs 3D-confirmed
  wherever the crossover is claimed.
