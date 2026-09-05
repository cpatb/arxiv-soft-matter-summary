# biofluids_venture: foundation for a bodily-fluid simulation and testing company

**Created 2026-09-05.** Working name only; the company has no name yet (see `07_DECISIONS.md`, D1).

This directory is the design of a company, not yet its code. It answers four questions the
founder asked: what niche the company fills, how it is architected, what simulation capability
it needs, and how to build the software so that it is robust rather than merely running. The
raw material is the existing portfolio in `Desktop/Projects`, most of which was built with
this company in mind whether or not it said so at the time.

## The thesis in one paragraph

Nearly every simulation sold to medical-device, drug-product and personal-care companies
treats the working fluid as water, or as water with a viscosity multiplier. Blood, mucus,
synovial fluid, high-concentration antibody formulations, thickened liquids for dysphagia,
urine, saliva and semen are none of those things. They shear-thin, yield, build and lose
structure in time, and jam. The error this introduces is not small, and it is almost never quantified. The company sells the missing link: measure the fluid; identify a
constitutive model with its uncertainty and its limits stated; carry that model into the solver
without the silent distortions each solver adds; validate against a simplified physical model
built in-house; and deliver all of it inside the credibility framework regulators already
recognise (ASME V&V 40 and the FDA's 2023 guidance on computational modelling in device
submissions). Nobody currently sells that chain as one service, and the closest analogue in
solid mechanics (Veryst's PolyUMod material library sold alongside consulting) has no fluid
counterpart.

## What is in here

| file | what it settles |
|---|---|
| `01_LANDSCAPE.md` | who already does what, with sources; which slots are empty and why |
| `02_NICHE_AND_POSITIONING.md` | the niche statement, what makes the service unique, the three regulatory tiers, and the ranked entry wedges |
| `03_ARCHITECTURE.md` | the technical stack: seven layers, their contracts, and the simulation implementations each needs |
| `04_ASSET_MAP.md` | every existing project mapped onto the stack, with its readiness and what must be redone |
| `05_ROADMAP.md` | phases, gates and budgets under the 5 to 10 hr/week constraint |
| `06_RISKS_AND_KILL_CRITERIA.md` | what would make this plan wrong, and the test that decides it |
| `07_DECISIONS.md` | the decisions only the founder can make, each with a recommendation |
| `docs/superpowers/specs/` | the design spec in the repository's standard form |

State card: `Obsidian Vault/Status/Biofluids Venture.md`.

## Provenance rules for this directory

Market and regulatory claims carry a source link or the tag `[UNVERIFIED]`. Numbers from
vendor market reports are marked as such; they are advertising, not measurements. Nothing in
here about the existing projects was taken from their READMEs alone; each was opened.
