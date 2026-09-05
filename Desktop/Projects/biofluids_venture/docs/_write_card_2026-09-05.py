"""One-shot: write the Obsidian Status card, router row, dormant pointer, README entry.
Kept for provenance; safe to delete after the vault validates."""
from pathlib import Path

STATUS = Path(r"C:\Users\compl\Documents\Obsidian Vault\Status")

card = r'''---
type: status
project: biofluids-venture
tier: 1
status: active
health: yellow
updated: 2026-09-05
tags: [status, biofluids, rheology, company, monetization]
---

# Biofluids Venture — Bodily-Fluid Simulation and Bench Consulting

## 1. One-liner

The post-PhD company: simulation and simplified benchtop work for products whose working
fluid is a bodily fluid or a drug formulation that does not behave like water — medical
devices, drug development, personal-care and human-fluid products — sold as a
measurement → identifiable model → solver → bench validation → credibility-package chain.
Design set written 2026-09-05 in `Desktop\Projects\biofluids_venture\`. **No product code
exists yet.**

## 2. Where it stands

- **Eight design documents written, zero code.** `README.md`, `01_LANDSCAPE`, `02_NICHE`,
  `03_ARCHITECTURE`, `04_ASSET_MAP`, `05_ROADMAP`, `06_RISKS`, `07_DECISIONS`, plus the
  Phase-0 spec at `docs\superpowers\specs\2026-09-05-biofluids-venture-design.md`
  (status: awaiting founder review). Health is yellow because the whole thing is a
  hypothesis until the founder reads it.
- **Niche sentence:** credible simulation and benchtop testing for products whose working
  fluid is not water. Working line: *"When the fluid isn't water."* The sharpest verified
  gap: Veryst's PolyUMod is a cited, validated material library for solids; **no fluid
  counterpart exists** and no consolidated biofluid constitutive database exists either.
- **Seven-layer stack** (placeholder name `biorheo`): L0 registry of cited, versioned
  fluid records with envelope + identifiability; L1 constitutive (GNF, structural,
  suspension, viscoelastic, damage kernels, every regularisation side by side); L2
  identify (dst_fit generalised); L3 solver bridge (OpenFOAM, Basilisk, COMSOL/Fluent
  exports, round trip mandatory); L4 bench rigs; L5 V&V 40 credibility engine; L6
  vertical apps. Three regulatory tiers (consumer / drug-dev internal / device) on one
  stack.
- **Wedges ranked:** W1 injectability (primary, closed-form first product, syringe rig);
  W2 airway/mucus; W3 thrombosis-as-jamming (flagship paper + open benchmark, never an
  absolute predictor); W4 urology/splash (public case studies); W5 later.
- **Existing assets mapped to layers** (`04_ASSET_MAP.md`): `dst_fit` seeds L2; `tier0`,
  `bio_fluids_bridge`, `carbopol`, `polydisperse_jamming` seed L1; the Basilisk solvers
  and bidet OpenFOAM cases seed L3; reviewers + skills are the quality system in
  prototype. Nothing exists for L0, L3 exporters, L5, rig build sheets, viscoelastic
  family, poroelastic rung.
- **Kill test is R1:** ten discovery conversations; fails if fewer than 3/10 describe a
  situation where the fluid model was the issue. Then licensing is dropped and the
  company narrows to consumer tier + bench.

## 3. Next action

Founder reads the design set (start with `README.md`, then `02`, then `07`) and answers
the three decisions that block Phase 0: **D1** package name (`biorheo`, `biofluid`,
`biorheology` unregistered on PyPI 2026-09-05; GitHub and trademark unchecked), **D2**
first wedge (recommended W1 injectability), **D5** what to do with `dst_fit`'s failing
gate G4c (recommended: retire with reason recorded). On approval, the next agent session
invokes `writing-plans` on the spec to produce the Phase 0 implementation plan; unit 0.1
is the monorepo with `dst_fit` migrated behind a generalised `Model` interface, gate =
its tests green from the new location and byte-identical figures.

## 4. Locations

| role | path |
|---|---|
| design set | `C:\Users\compl\Desktop\Projects\biofluids_venture\` |
| Phase-0 spec | `C:\Users\compl\Desktop\Projects\biofluids_venture\docs\superpowers\specs\2026-09-05-biofluids-venture-design.md` |
| predecessor product spec (rheology line) | `C:\Users\compl\Desktop\Projects\rheo_products\SPEC_2026-08-27.md` |
| L2 seed | `C:\Users\compl\Desktop\Projects\dst_fit\` |
| L1 seeds | `C:\Users\compl\Desktop\Projects\tier0_wyart_cates_closure\`, `C:\Users\compl\Desktop\Projects\bio_fluids_bridge\`, `C:\Users\compl\Desktop\Projects\carbopol\` |
| L6 seeds | `C:\Users\compl\Desktop\Projects\thrombosis_jamming\`, `C:\Users\compl\Desktop\Projects\uroflow_splashback\`, `C:\Users\compl\Desktop\Projects\glymphatic_transport\` |
| intent record | `Logs (Raw)\2026-08-04.md` (human-authored; quoted in [[Dormant & Prototypes]]) |

## 5. Run / verify

Nothing runs yet. When unit 0.1 exists: `pytest` from the monorepo root must reproduce
`dst_fit`'s current 205-green count and its figures byte-identically.

## 6. Load-bearing facts

- **Verified 2026-09-05 (web):** ASME V&V 40-2018 and the FDA 2023 credibility guidance
  are the frame for the device tier; FDA's nasal-spray CFD-linked bioequivalence pathway
  exists; hemolysis CFD absolute predictions are off by ~3 decades across benchmarks
  (relative rankings usable) — every damage kernel ships **relative-only**.
- **Vendor market figures in `01_LANDSCAPE.md` are marketing, not measurements.**
- `si_fit` is withdrawn ([[Drop Impact — Simulation]]); it must not migrate into any
  layer. `impact_stress/materials.py` still ships it.
- The Beris sphere yield number 0.143 belongs to the **middle** of three conventions
  (factor-of-three spread); every L0 record carrying a yield number must name its
  convention.
- The August spec's headline (corn-starch `phi_m` drift) was retired within a week of its
  runs because it rested on disowned curves — `06_RISKS` R7: every headline here is a
  hypothesis until its gate has run.

## 7. Open questions

D1 name, D2 wedge, D5 G4c (blocking). D3 bench work yes/no, D4 consulting-first, D6 PhD
end date and advisor sign-off before anything public, D7 thrombosis paper ownership,
D8 `rehab_smartfluid` vertical or paper, D9 licences — all in `07_DECISIONS.md` with
recommendations.

## 8. Landmines

- **The design set is agent-written in one night** and has not been read by the founder.
  Treat every market claim as unverified beyond the cited links; treat the roadmap's
  month numbers as guesses until D6 anchors them.
- `bloodstain_bpa` has published prior art (Attinger group) — no novelty language.
- Uroflow splash-index magnitudes are ranking-only until the calibration clip is analysed.
- The bidet design model evaluates its closure outside its fit band (`r_max/d*` = 23–26
  against [2, 20]) — bounded before it is shown to a client.
- The brainstorming skill's gate: **no implementation before founder approval of the
  spec.** Do not start unit 0.1 on the strength of this card.

## 9. History

- 2026-09-05 — Card created. Eight design documents + Phase-0 spec written; promoted from
  the [[Dormant & Prototypes]] "bodily-fluid dynamics" intent entry.
- 2026-08-04 — Author states intent to make bodily-fluid dynamics a research area
  (`Logs (Raw)\2026-08-04.md`).

## 10. Links

[[_Router]] · [[dst-fit — Suspension Constitutive Fitting]] · [[Drop Impact — Substrate Grant]]
(HB solver) · [[Tier-0 Closure — Calibration and Onset]] · [[Dormant & Prototypes]]
'''
(STATUS / "Biofluids Venture.md").write_text(card, encoding="utf-8", newline="\n")

# Router row: insert before Dormant & Prototypes; update dormant row + updated date
r = STATUS / "_Router.md"
t = r.read_text(encoding="utf-8")
row = ("| [[Biofluids Venture]] | 🟡 | **The post-PhD company.** Bodily-fluid simulation + bench consulting for devices, drug products, human-fluid products. Design set only, no code. "
       "| 🆕 **05 SEP: eight design docs + Phase-0 spec written overnight** in `Desktop\\Projects\\biofluids_venture\\`. Founder reads `README` → `02` → `07` and answers D1 (name), D2 (wedge; rec. W1 injectability), D5 (dst_fit G4c). **Nothing is built until the spec is approved** |\n")
marker = "| [[Dormant & Prototypes]] |"
assert t.count(marker) == 1
t = t.replace(marker, row + marker)
old = "**Bodily-fluid dynamics named as a research direction 04 AUG** — intent recorded, nothing done. Needs a defined question before it earns a card"
assert t.count(old) == 1, t.count(old)
t = t.replace(old, "Bodily-fluid dynamics promoted to its own card 05 SEP → [[Biofluids Venture]]")
t = t.replace("updated: 2026-09-03", "updated: 2026-09-05", 1)
r.write_text(t, encoding="utf-8", newline="\n")

# Dormant note pointer
d = STATUS / "Dormant & Prototypes.md"
t = d.read_text(encoding="utf-8")
old = "> Promoting this to a Tier-1 card needs a defined question, not a subject area. Until then\n> it stays here. Per [[_Design]], a Tier-2 entry is promoted by *writing it a card*."
assert t.count(old) == 1, t.count(old)
t = t.replace(old, old + "\n>\n> 🆕 **2026-09-05 — promoted.** The card is [[Biofluids Venture]]; the design set is in\n> `Desktop\\Projects\\biofluids_venture\\`. This entry stays as the record of the intent.")
d.write_text(t, encoding="utf-8", newline="\n")

# Projects README
p = Path(r"C:\Users\compl\Desktop\Projects\README.md")
t = p.read_text(encoding="utf-8")
old = "| `rheo_products/` | Product spec for the computational-rheology line. |"
assert t.count(old) == 1, t.count(old)
t = t.replace(old, old + "\n| `biofluids_venture/` | Design set for the post-PhD company: bodily-fluid simulation + bench consulting. Supersedes `rheo_products/` as the business plan; no code yet. |")
p.write_text(t, encoding="utf-8", newline="\n")
print("ok")
