# 02. Niche and positioning

## 1. The niche, in one sentence

**Credible simulation and benchtop testing for products whose working fluid is a bodily fluid
or a drug formulation that does not behave like water.**

The important word is *credible*. The niche is not "CFD for medical devices"; that
market is served. It is the material-model layer underneath it, which nobody productises, plus
the evidence chain that turns a simulation into something a reviewer can accept.

## 2. What the company sells

One chain, sold whole or in parts:

1. **Characterise.** Rheometry on the client's fluid or a validated surrogate, with the
   protocol chosen to make the parameters identifiable rather than chosen from habit.
2. **Identify.** A constitutive model, its parameters with confidence intervals, the
   correlation structure, and a verdict per parameter: measurable, bounded only, or not
   identifiable from this data. Plus the protocol that would resolve what is not.
3. **Carry.** The model exported into the client's solver (OpenFOAM, Basilisk, COMSOL, Fluent,
   SimVascular) with the solver's own regularisation named, its induced error computed, and a
   round-trip check that the solver reproduces the rheometer.
4. **Predict.** The device or process simulation itself, at the fidelity the question needs:
   closed-form and 1-D where that answers it, full CFD or particle simulation where it does not.
5. **Validate.** A simplified physical model in-house: syringe-and-needle force rig, jet and
   spray rig, drop and splash rig, capillary viscometer, flow loop; high-speed imaging with
   audited calibration; uncertainty on both sides of the comparison.
6. **Certify.** A V&V 40 credibility package assembled from the pipeline's own manifest:
   question of interest, context of use, model risk, verification with observed orders,
   validation with uncertainty, applicability. Missing evidence is listed, not hidden.

## 3. What makes it unique, and why each item is hard to copy

| feature | why competitors do not have it |
|---|---|
| Identifiability verdicts ("you cannot have phi_m from this curve, and here is why") | Requires having done the degeneracy analysis on real data. `dst_fit` has done it for dense suspensions; nobody has shipped it for any fluid class. |
| Regularisation-aware solver export with round-trip proof | Requires having been burned by it. The `yfrac-is-a-clamp-readout` episode (a month of a yield-stress observable measuring the regularisation ceiling) is where this requirement came from. |
| Benchtop surrogates plus imaging plus simulation from one shop | CFD shops do not own rigs; testing labs do not own solvers. The founder has run both for five years. |
| Pre-registered gates enforced in code, failures reported as results | A process, not a feature. Most consultancies cannot afford to tell a client the model cannot answer the question. It is hard to copy because it costs revenue in the short run. |
| Structural state-variable closures spanning silica, blood, thrombosis and smart fluids | The jamming and thickening physics is the founder's research; the transfer to blood and clotting is already prototyped. |
| A cited, versioned fluid registry | Nobody has curated it because it is slow, unglamorous, and needs both rheology judgement and literature discipline (`lit-verify`, DOI resolution) that the existing tooling here already enforces. |

## 4. Three regulatory tiers, one stack

The same physics and code serve three markets with very different burdens. This is the
deliberate structure of the business: cash from the light tier funds the heavy one.

| tier | examples | burden | what is sold | why it fits now |
|---|---|---|---|---|
| **Consumer and personal care** | bidets, toilets and urinals (splash), water flossers, nasal rinse, menstrual and incontinence absorbents, diapers, feeding bottles, breast pumps, thickened drinks, skincare with yield stress | none beyond the client's own QA | feasibility studies, design optimisation, splash and aerosol assessments, absorbency modelling | the bidet paper and uroflow model are finished marketing assets; sales cycles are short; no QMS needed |
| **Drug development, internal decision support** | injectability of viscous biologics, fill-finish drip and splash, needle shear and protein aggregation, nasal and inhaled deposition, dissolution hydrodynamics (USP II), TFF of high-concentration mAb, depot injection back-pressure | GMP-adjacent; simulations inform decisions but are not themselves submitted, or are submitted as supporting evidence | characterisation plus simulation studies; solver exports; UQ on formulation variability | largest budgets; the fluid is genuinely non-Newtonian and the industry knows it; surrogates and client-supplied rheometry make it workable without biologics on site |
| **Medical devices, regulatory** | catheters, ECMO and oxygenators, blood pumps, dialysers, autoinjectors as combination products, ophthalmic injection, viscosupplements | V&V 40, FDA credibility guidance, ISO 13485 environment | credibility packages, validation studies, expert reports; initially as subcontractor to established consultancies | highest value per study and the hardest to copy, but requires a track record and insurance the company does not yet have |

## 5. Entry wedges, ranked

Ranked on: transfer from existing work, whether non-Newtonian behaviour actually changes the
answer, budget, whether it can be validated on a benchtop without access to biologics or
blood, and regulatory pull.

**W1. Injectability and drug-product fluid handling (primary).**
Viscous biologics through needles, autoinjectors, needle-free, depot gels, fill-finish.
Pressure-driven flow through a capillary is analytically tractable for power-law and
Herschel-Bulkley fluids (Buckingham-Reiner), so the first product can be closed-form with a
CFD check, which is exactly the fidelity ladder in `03_ARCHITECTURE.md`. The benchtop rig is a
syringe, a load cell, a needle and a camera. Surrogates (glycerol, PEG, Carbopol, xanthan,
silica) are in hand; real formulation rheometry comes from the client under NDA. Existing
budget line: ISO 11608 injection-force testing is already routine, so the company is selling
into a cost centre that exists.

**W2. Airway, nasal and mucus (second).**
Sprays, nebulisers, nasal rinse, mucus clearance, cough. Shares the jet and spray rig with W1
and the founder's high-speed imaging skill directly. Mucus simulants (Carbopol, HEC, mucin
solutions) are standard, and the Carbopol work is the yield-stress anchor. The FDA has
already signalled that CFD deposition can carry bioequivalence weight, which makes the
credibility package sellable.

**W3. Device-induced thrombosis and hemolysis (flagship research asset, not first revenue).**
The `thrombosis_jamming` closure is novel as far as the prior-art check has gone (an arXiv
precedent is still flagged for verification in `04_ASSET_MAP.md`): occlusion as a stress-activated jamming
transition with a critical shear rate, where incumbent models produce arrest by deposition
kinetics. The field's own benchmark work concedes that absolute hemolysis predictions fail by
three decades. This is the paper and the open benchmark that give the company a scientific
identity; it is not the first invoice, because validation needs blood work the company cannot
do alone and the buyers are large medtech with long cycles.

**W4. Urology, continence, surgical splash and aerosol (public case studies).**
Uroflow splashback, bidet jets, drop impact on wet surfaces, catheter flow. Finished or
nearly finished. Small budgets, short cycles, publishable, and a natural consumer-tier entry.
These become the free case studies and the first two or three paid feasibility studies.

**W5. Ophthalmic, joint, GI, dysphagia (later, same stack).**
Vitreous and intravitreal injection; synovial fluid and viscosupplements; gastric fluid and
dissolution hydrodynamics; IDDSI thickened liquids. All are constitutive-identification
problems the registry and identification layers will handle once they exist. Dysphagia in
particular is small, open, and has both a medical and a food-industry buyer.

Recommendation: **W1 primary, W2 second, W3 as the flagship paper, W4 as the case-study
line.** They share the benchtop rig, the registry, and the identification code.

## 6. What the company is not

- It does not compete as a general CFD shop. Newtonian device flows are already served, and
  competing there is a price war.
- It does not do preclinical image analysis for hire. That was declined on record
  (`monetization-direction`).
- It is not branded as "rigour for AI-generated science". The rigour is a feature of the named
  work rather than the identity.
- It does not sell a solver. Solvers are open (OpenFOAM, Basilisk) or the client's; the product
  is the material layer, the bridge, the validation, and the credibility trail.

## 7. Working positioning line

"When the fluid isn't water." Alternatives and the naming decision are in `07_DECISIONS.md`.
