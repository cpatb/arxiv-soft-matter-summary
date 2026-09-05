# 07 — Decisions only the founder can make

Each has a recommendation. None has been acted on. Where a decision changes what gets built
first, the roadmap unit it gates is named.

## D1 — Name and package name

Gates: 0.1 (repository name), everything public.

The company name is a founder decision and a trademark search. The package name matters
sooner because it is in every import. `biorheo` is used throughout these documents as a
placeholder. On 2026-09-05 `biorheo`, `biofluid` and `biorheology` returned 404 on PyPI (unregistered); GitHub and trademark were not checked. Do not adopt from these documents alone.
Working positioning line: "When the fluid isn't water." Alternatives worth keeping:
"Fluids the body makes"; "Constitutive models for biological fluids, with the uncertainty
attached."

## D2 — Which wedge first

Gates: 0.6 (first case study), 2.1–2.2 (which rig is built first).

**Recommendation: W1 injectability**, for the reasons in `02_NICHE_AND_POSITIONING.md` §5: the
first product is closed-form, the rig is a syringe and a load cell, surrogates are in hand,
and the industry already spends on injection-force testing. **Alternative: W2 airway**, which
uses the imaging skill more and has clearer regulatory pull, but needs the spray rig and
free-surface CFD before it produces anything sellable. Choosing W2 first swaps units 0.6 and
1.6 and moves the spray rig ahead of the needle rig.

## D3 — Is this a company that does bench work

Gates: the whole of Phase 2; whether L4 exists as a product.

The differentiation argument rests on owning rigs. That means space, equipment and the
founder's hands, not only overnight compute. **Recommendation: yes**, at the scale of a bench
in a rented lab or a shared space, with the high-speed camera as the only expensive item. If
the answer is no, the company is a software and analysis firm, L4 becomes "we analyse your
bench data," and the positioning line has to change, because "validated against our own
measurements" is the sentence that separates it from the CFD shops.

## D4 — Consulting-first or product-first

Gates: Phase 2.5 versus a licensing tier in Phase 2.

**Recommendation: consulting first, product second**, the Veryst order. Studies fund the
registry and the exporters; the licensing tier sells what studies have already used. The
reverse order requires a sales channel the founder does not have and a product that has
never been used in anger. The August spec reached the same conclusion for the narrower
rheology line.

## D5 — What to do with `dst_fit`'s gate G4c

Gates: 0.1, because the package migrates with its tests.

Still open from the August spec. The gate requires cross-curve consistency of phi_m and fails
at 1 of 6 on PMMA. **Recommendation: retire it with the reason recorded** — the work has shown
phi_m is not cross-curve consistent in any material tested, so the gate measures the physics,
not the fitter, and a permanently failing test that everyone learns to ignore is worse than a
documented retirement. Keep the number in the README.

## D6 — The PhD end date, and whether any of this is public before it

Gates: phase boundaries; whether case studies carry the founder's name and affiliation.

The roadmap is in elapsed months and needs an anchor. Separately, publishing case studies and
an open benchmark under the founder's name while still a student may need advisor and
institutional sign-off, particularly where lab equipment or lab-collected data is involved.
**Recommendation:** raise it with the advisor before Phase 0.6 is public; use only own-time
surrogate data in anything public until then.

## D7 — The thrombosis paper: company asset or thesis chapter

Gates: 2.6.

If it becomes part of the thesis, its authorship and timing are the university's. If it is
the company's flagship, it is written on own time and the company's name goes on the
benchmark. **Recommendation:** decide by whether the advisor wants it; the science is the same
either way, and the open benchmark can carry the company's name regardless.

## D8 — `rehab_smartfluid`: a vertical or a paper

Gates: none in Phases 0–2.

Protective and orthotic devices with shear-thickening elements are medical or sports
devices, and the activation-number result is a device-design rule. It could be a W6. **Recommendation:** paper first; revisit as a vertical only if a device company appears in the discovery calls.

## D9 — Open subset of the registry and benchmark: licence

Gates: 2.3, 2.4.

**Recommendation:** CC-BY 4.0 for records and benchmark definitions; a permissive software
licence (BSD-3 or MIT) for the open layers; the identifiability-graded tier and the exporters
proprietary. Decide before anything is published, because it cannot be narrowed afterwards.
