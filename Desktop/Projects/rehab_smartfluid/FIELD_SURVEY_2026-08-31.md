# Field survey: transient shear thickening, and the device literature that ignores it

**Date:** 2026-08-31 · arXiv (cond-mat.soft, physics.flu-dyn) plus web search for the
engineering-journal side. First pass, not exhaustive.

---

## 1. The headline

There are two literatures here and **they do not cite each other**.

The **physics** side has spent six years building transient constitutive models for
shear-thickening suspensions, largely out of the Cates, Mari and Forterre groups. The
**device** side designs shear-thickening dampers, orthoses and protective pads using
steady flow curves and Bouc-Wen hysteresis fits.

Nobody has carried the first into the second. That gap is the company, and it is
narrower and better defined than the design document assumed.

## 2. What the physics side already has

**Take this seriously: the constitutive model is not the novel part.**

- **Gillissen, Ness, Peterson, Wilson & Cates**, *Constitutive model for time-dependent
  flows of shear-thickening suspensions*, arXiv:1911.00280. A tensorial model coupling a
  microstructure evolution equation to stress-activated friction. This is the grown-up
  version of the 0-D strain clock in `layer2/activation.py`, and it is the natural
  layer 2 for phase 1. **Read before writing any more constitutive code.**
- **Athani, Metzger, Forterre & Mari**, *Transients in shear thickening suspensions:
  when hydrodynamics matters*, arXiv:2408.15130. Particle-based, pressure-imposed. Finds
  that approaching shear jamming the transient stress does **not** diverge the way the
  steady-state flow curve predicts. This both supports the shortfall result in
  `SCALING_2026-08-31.md` §4 and warns that a 0-D model will get its magnitude wrong.
- **Acharya & Mari**, *Modelling the evolution of flow-induced anisotropy of concentrated
  suspensions*, arXiv:2506.05222. Steady and transient, 2025.
- **Han, Wyart, Peters & Jaeger**, *Constitutive relations for shear fronts in
  shear-thickening suspensions*, arXiv:1711.02196. Sudden driving at a boundary, which is
  the impact-loading case a fall-protection device sees.
- **Nakanishi, Nagahiro & Mitarai**, arXiv:1010.5063 and arXiv:1112.4593. Phenomenological
  state-variable fluid dynamics for dilatant fluids. The precedent for treating the
  microstructure as one scalar state, which is what `activation.py` does.
- **Mobaseri, Gordillo, Burton, Yoon, Lee, Kumar, Driscoll & Cheng**, *Inertia-Dilatancy
  Interplay Governs Shear-Thickening Drop Impact*, arXiv:2601.12642. Already known to this
  vault as the inertial-dilatant number. Note the author's own collaborators are on it.

**Consequence for the novelty claim.** The strain-clocked element is the simplest
possible version of published work, not a new model. The defensible claim is the
*design map* and the *application*: nobody has asked what these transient results imply
for a device operating at gait rates. Do not write a paper claiming the constitutive
model.

## 3. What the device side has, and what it assumes

Searched arXiv for shear-thickening fluid combined with wearable, damper, protective or
impact: **7 results in total**, and none of them are wearables. The device literature is
in engineering journals, not preprints. From the web search side:

- STF variable-stiffness prosthetic foot, *Smart Mater. Struct.* 2020, with critical
  shear rates measured at **10-80 1/s**. This is phase 0's reproduction target.
- STF-embedded foam protective wearables, strain reduction **30% (foam) to 50% (TPU)**.
- STF dampers characterised by force-velocity hysteresis loops; squeeze mode giving
  38.48 kPa normal stress at 71 wt% silica against a 6.38 kPa peak in shear.
- MR dampers characterised across 0-2 A, 0.5-4 Hz, 1-2.5 mm amplitude.
- **Montenegro & Galindo-Rosales**, arXiv:2305.05306, entry flow of continuous
  shear-thickening fluids in microtubes, explicitly motivated by energy-dissipating
  devices (Rheinforce). The nearest thing to a bridge found so far, and it is a flow
  problem rather than a device design.

Every one of these is quasi-steady or phenomenological. None carries a packing fraction.

## 4. What this changes about the plan

1. **Do not build a constitutive model from scratch in phase 1.** Port Gillissen et al.
   The 0-D strain clock stays as the cheap screening tool and as the thing the tensorial
   model must reduce to.
2. **Athani et al. is a validation target and a threat.** If transient stresses near
   jamming are lower than steady state predicts for reasons a 0-D model cannot capture,
   then `SCALING` §4's amplification is directionally right and quantitatively wrong.
   Reproducing their result is a sharper phase-1 gate than anything currently in the
   design document.
3. **The corpus for the extrapolation gate is thinner than hoped on one axis.** The
   transient literature is largely simulation, not the multi-formulation experimental
   ladders the gate needs. The weight-fraction ladders found on the material side
   (9/12/15/17/20 wt%, and the two-particle-size study) remain the corpus of record.
4. **The two-literature gap is worth a preprint on its own**, independent of any device.
   That is a cheap, publishable phase-0 output that also functions as recruiting material.

## 5. Not yet searched

- Non-English and patent literature on STF orthoses.
- The MR side's transient literature, which may be further ahead than the STF side given
  it is a commercial technology.
- Whether anyone has measured `gamma_0` for a silica/PEG STF specifically. `SCALING` §6
  names this as the single measurement that would most change the conclusion.
