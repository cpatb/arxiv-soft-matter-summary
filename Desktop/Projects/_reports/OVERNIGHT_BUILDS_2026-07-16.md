# Overnight builds — 2026-07-16

Four small, self-contained, validated computational projects built overnight (with Claude),
each mapping to a research/portfolio thread. See the Obsidian vault's *Frontiers 2026-07 MOC*
and *Simulations MOC* for the surrounding plan.

## 1. `tier0_wyart_cates_closure/`  ⭐ portfolio piece #1
Reduced-order **Wyart–Cates suspension closure** (the Tier-0 rung of the drop-impact sim plan).
- `wyart_cates.py` — closure physics (steady flow curve + two-state reversibility latch).
- `run_sanity_gate.py` — reproduces DST/shear-jamming (silica) vs CST+yield (corn starch) and the
  retention-vs-recovery reversibility switch. All gates PASS.
- `test_wyart_cates.py` — 6 invariant tests (pass, also runs under pytest).
- `diff_calibrate.py` — **PyTorch-autodiff calibration**: recovers ground-truth closure parameters
  from noisy synthetic viscometry by gradient descent (the "make the closure differentiable"
  frontier move). Validation PASS.

## 2. `bio_fluids_bridge/`  — suspension physics → biological fluids
Reuses the same state-variable closure to model **blood** (shear-thinning via RBC rouleaux
aggregation) beside silica (shear-thickening) — one framework, opposite sign of the
structure→viscosity coupling — plus a blood thixotropic hysteresis loop.
- `blood_rheology.py`, `run_blood_bridge.py` (validation PASS).
- Entry rung for the biological-fluids research bridges (thrombosis, sickle cell, tissue).

## 3. `thrombosis_jamming/`  — the closure as arterial thrombosis (third fluid)
Shear-induced platelet aggregation (SIPA) as a stress-activated jamming transition (VWF unfolding
switches on adhesive contacts → bonded network jams). Produces the **falsifiable signature** — a
critical occlusion threshold + divergence of occlusion time — that distinguishes a jamming
mechanism from the incumbent deposition-kinetics thrombosis models.
- `thrombosis_jamming.py`, `run_thrombosis.py` (signature validation PASS).
- The seed of a wet-lab-free first paper (see vault: *Frontiers/Thrombosis as Jamming*).

## 4. `glymphatic_transport/`  — brain clearance regime map (a different flavour)
A Péclet-number dimensional-analysis map (NOT a sim) operationalizing the brain-clearance
advection-vs-diffusion controversy: perivascular = advective (Pe≈7.5), deep parenchyma =
diffusive (Pe≈0.05), matching the field synthesis. The entry rung for a "fluid dynamics and
the brain" track (see vault: *Frontiers/Fluid Dynamics & the Brain*).

## Common thread
The first three are *one state-variable jamming/rheology closure, many fluids* — silica (drop impact)
→ blood (rheology) → platelets (thrombosis), all from the same skeleton; the fourth is a transport
(Péclet) map for the brain. The unifying research vision is in the vault:
**Frontiers/Research Vision — Jamming Closures Across Fluids**.

## Reproduce
Each folder has a `requirements.txt`. All run on a CPU workstation:
```
cd tier0_wyart_cates_closure && python run_sanity_gate.py && python test_wyart_cates.py && python diff_calibrate.py
cd ../bio_fluids_bridge && python run_blood_bridge.py
```

## Portfolio next steps (from the frontier notes)
- Package `tier0-closure` as an installable, tested library; publish to GitHub with the figures.
- Extend differentiable calibration through the transient ODE (JAX/diffrax) + uncertainty quantification.
- Calibrate against real silica/corn-starch flow curves; then couple into the Basilisk drop sim.
