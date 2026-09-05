# 01. Landscape: who does what, and where the gaps are

Surveyed 2026-09-05. Web-checked claims carry a link. The August 2026 software survey in
memory (`computational-rheology-software-landscape`) is folded in and not repeated in full.

## 1. The regulatory frame that makes this a business

- **ASME V&V 40-2018** is the FDA-recognised, risk-based framework for deciding how much
  credibility a computational model needs before it can support a device decision. It was
  written jointly by device manufacturers, simulation vendors and the FDA.
  ([ASME](https://www.asme.org/codes-standards/find-codes-standards/assessing-credibility-of-computational-modeling-through-verification-and-validation-application-to-medical-devices),
  [overview](https://dmd.umn.edu/2023/vv-40))
- The FDA's 2023 guidance *Assessing the Credibility of Computational Modeling and Simulation
  in Medical Device Submissions* is written to be consistent with V&V 40 and points sponsors to
  the Q-Submission programme for feedback on credibility plans.
  ([FDA PDF](https://www.fda.gov/media/154985/download))
- For nasal and inhaled drug products, the FDA has funded and published CFD-based deposition
  work, and CFD-linked deposition is positioned as a route to bioequivalence with reduced
  clinical burden. ([FDA CDER spotlight](https://www.fda.gov/drugs/spotlight-cder-science/computational-simulations-shed-light-factors-affecting-nasal-spray-distribution),
  [GDUFA report](https://www.fda.gov/industry/generic-drug-user-fee-amendments/fy2016-regulatory-science-report-locally-acting-orally-inhaled-and-nasal-drug-products))
- In-silico evidence for devices was formally recognised in the FDA's 2021 reporting guidance;
  vendor market reports put the in-silico-trials market around $2B in 2026. Treat the dollar
  figure as marketing, not data. ([systematic review](https://link.springer.com/article/10.1007/s43441-025-00871-2))

What this means commercially: a simulation with a credibility trail is worth more than the
same simulation without one, and the trail is a specifiable, auditable deliverable. That is
the thing this company should be built to produce by default.

## 2. Who sells simulation to these industries today

**Full-service engineering consultancies.** Veryst Engineering is the template for the shape
of this company: in-house materials testing, simulation, medical-device focus, and a
productised material-model library (PolyUMod) sold into COMSOL/Abaqus alongside consulting.
Their published capability lists "non-Newtonian fluids" under CFD, but their material-model
product is for polymers and solids. There is no fluid PolyUMod.
([Veryst simulation](https://www.veryst.com/services/simulation-analysis),
[Veryst medical](https://www.veryst.com/services/specialized-expertise/medical-devices))
Others in the same shape: Exponent, Thornton Tomasetti's applied-science group, SimuTech
(Ansys channel), CFD Direct (OpenFOAM). All have solvers; none have rheology as a product.

**In-silico platform companies.** Dassault (Living Heart, BIOVIA), Certara (PK/PBPK), ELEM
Biotech, HeartFlow, Materialise. These operate at organ or population scale with a fixed
material model, usually Newtonian or Carreau blood. They are customers or partners for a
material layer, not competitors to it.

**Contract rheology and testing labs.** Instrument vendors' application labs, contract
analytical labs, university core facilities. They return a flow curve and a fitted
Herschel-Bulkley triple. They do not say whether the parameters are identifiable, and they do
not simulate.

**Academic codes.** rheoTool (OpenFOAM, behind current releases), LF-DEM (dense suspensions,
unpackaged), SimVascular (Newtonian/Carreau blood, vascular focus), Basilisk (free surface, no
ecosystem). Each proves feasibility of one layer; none is a product.

**Fitting software.** pyRheo, hermes-rheo, RepTate own viscoelastic model fitting for
polymers. None handle dense suspensions, thixotropy with identifiability, or the handoff to
a solver. (Detail in the August memory note.)

## 3. The state of the science where the money is

**Blood damage.** Absolute CFD predictions of hemolysis in the FDA benchmark nozzle are wrong
by roughly three orders of magnitude; relative predictions between geometries are usable.
The power-law coefficients come from uniform-shear Couette devices and do not transfer to
device flows. ([PMC 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10101913/),
[FDA nozzle non-Newtonian study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3965442/),
[UQ approach 2024](https://arxiv.org/pdf/2407.18757)). This is a field that already knows it
needs identifiability verdicts and does not have a vendor for them.

**High-concentration antibody products.** Every biologic is moving to subcutaneous delivery,
and the formulations at 100 to 300 mg/mL are viscous, often shear-thinning, and the viscosity
limits device design (injection force, needle gauge, dose time). Tissue back-pressure and
depot mechanics are being modelled with porohyperelastic FE; the fluid is still usually taken
as Newtonian. ([Biophysical Reviews 2025](https://link.springer.com/article/10.1007/s12551-025-01346-2),
[mAbs 2025](https://www.tandfonline.com/doi/full/10.1080/19420862.2025.2577159),
[porohyperelastic SC injection](https://pubmed.ncbi.nlm.nih.gov/36529050/),
[autoinjector tissue study](https://link.springer.com/article/10.1007/s10237-023-01746-x))

**Mucus.** A 2025 *Scientific Reports* paper characterises human mucus as
thixo-viscoelastoplastic with a multimode BMP model. This is a constitutive form nobody
outside a few rheology groups can fit, and airway, nasal, GI and reproductive-health products
all sit on top of it.
([Sci Rep 2025](https://www.nature.com/articles/s41598-025-16178-z))

**Dysphagia thickened liquids.** IDDSI levels are defined by a gravity flow test through a
syringe. The literature now agrees that shear viscosity at 50 1/s is not sufficient and
extensional viscosity governs bolus cohesion, yet product classification still runs on the
syringe test. The modelling problem is open, nobody prominent is working on it, and it has a
medical and a consumer market attached.
([Sci Rep 2019](https://www.nature.com/articles/s41598-019-52269-4),
[Food Hydrocolloids 2023](https://www.sciencedirect.com/science/article/abs/pii/S0268005X23001492))

**Biofluid material data.** There is no consolidated, cited, versioned database of
constitutive parameters for human fluids. Parameters are scattered across papers, with
protocols, temperatures and donor conditions that rarely match. A search for one returns
review chapters. ([Springer chapter](https://link.springer.com/chapter/10.1007/978-1-4419-8119-6_4))

## 4. The empty slots

1. A **fluid material-model library and registry** sold the way PolyUMod is sold for solids.
2. **Identifiability verdicts** on rheological parameters as a deliverable.
3. The **measurement-to-solver handoff** with each solver's regularisation made explicit and
   checked by round trip.
4. **Credibility packages** (V&V 40) generated from the analysis pipeline rather than written
   after it.
5. **Simplified physical models** (benchtop surrogates with high-speed imaging) as the
   validation arm, offered by the same people who run the simulation.
6. An **open benchmark suite** for non-Newtonian biofluid device flows, the way the FDA
   nozzle and blood pump became benchmarks for Newtonian blood.

Slots 1 to 4 are software. Slot 5 is the founder's laboratory skill. Slot 6 is how a company
with no track record acquires one.
