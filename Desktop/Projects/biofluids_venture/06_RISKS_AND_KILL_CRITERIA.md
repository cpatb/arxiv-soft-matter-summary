# 06 — Risks and kill criteria

Each risk names the test that decides it and what happens if the test fails. Ordered by how
early it can be decided, because the cheap tests come first.

## R1 — Nobody pays for constitutive rigour

The whole thesis is that "we tell you which parameters you can have" is worth money. The
counter-hypothesis: buyers passed review last time with Carreau and textbook constants, and
will again.

**Test:** Phase 0.7, ten discovery conversations. Ask specifically: has a simulation ever been
questioned on the fluid model; has anyone been asked for uncertainty on a rheological
parameter; would a report that says "not identifiable from your data" be welcome or a
problem. **Fails if** fewer than three of ten describe a situation where the fluid model was
the issue. **Then:** the stack survives as a research tool and an open benchmark, the
company narrows to the consumer tier plus benchtop testing, and the licensing tier is dropped.

## R2 — Device-tier work needs a track record and a QMS a solo founder cannot supply

Regulatory submissions are reviewed for the sponsor's process as well as the model. A first
company with no ISO 13485 environment and no prior submissions will not be a prime contractor.

**Test:** whether an established consultancy will subcontract device-tier work on the strength
of a credibility manifest (Phase 3). **Fails if** two attempts return "we need to see prior
submissions." **Then:** the device tier stays a subcontract line indefinitely and the company's
revenue is the drug-development and consumer tiers, which do not require it. This is not a
kill; it is a ceiling.

## R3 — The registry cannot be curated at useful density

Literature parameters are scattered, protocols mismatch, and digitising figures is slow.
Facts are not copyrightable, but wholesale extraction of a paper's data may still breach
publisher terms, and the company must not build on data it cannot show the provenance of.

**Test:** Phase 2.4, 25 records with every DOI resolved and every envelope checked. Time it.
**Fails if** a record averages more than a founder-hour of review, or if a meaningful fraction
of needed sources are figure-only with no deposited data and restrictive terms. **Then:** the
registry is limited to deposited datasets (Edinburgh DataShare pattern), own measurements,
and client data under licence, and the open subset is smaller than hoped. The identification
layer is unaffected.

## R4 — Hemolysis and thrombosis prediction is a graveyard

The field's own benchmarks show absolute predictions off by three decades. Selling into it
with a new closure invites the same fate.

**Mitigation, not a test:** W3 is positioned as a paper and a benchmark, never as an absolute
predictor. All damage kernels ship labelled relative-only. The jamming closure's claim is a
critical point that discriminates mechanisms, which is a scientific claim with a named test
(Phase 2.6), not a product claim. **If the discriminating test finds the data too sparse,** the
paper says so and the closure stays a hypothesis in the registry with that status.

## R5 — The bridge is a maintenance treadmill

rheoTool fell behind OpenFOAM releases and stayed there. `libbiorheo` could do the same.

**Test:** whether the round-trip harness in CI catches a solver release change within a week
of it. **Mitigation:** scope `libbiorheo` to the models the company sells; support one
OpenFOAM release at a time; the harness is the definition of "supported."

## R6 — Time

5–10 hours a week is the binding constraint until the PhD ends, and the PhD comes first.

**Test:** Phase 0 exit at three months. **Fails if** fewer than five of seven units are
gated by then. **Then:** drop Phase 1.6 and 2.3 (the second case study and the benchmark
suite) rather than compressing everything; the chain matters more than its breadth.

## R7 — Reproducing the last plan's failure modes

The August product spec's own headline (the corn-starch phi_m drift) was retired within a
week of the runs that tested it, because it rested on curves the source labs had disowned.
Every headline in these documents is a hypothesis about the market or the physics until its
gate has run.

**Mitigation:** no claim in a case study or sales document that has not passed through the
same reviewers the research does; the `04_ASSET_MAP.md` §4 retirement list is enforced
before anything is shown.

## R8 — Liability

A simulation used in a device decision carries liability for the modeller. Insurance is a
Phase 3 prerequisite, listed there, and until it exists no deliverable may be described as
suitable for a regulatory decision. Reports before that point say so in the applicability
section.

## What would make the whole plan wrong

- A vendor ships a biofluid material library with identifiability grading into COMSOL or
  Ansys. Watch Veryst, Ansys Granta, and the instrument vendors. If it happens, the company is
  the benchtop validation arm and the open benchmark, which is still a business.
- The founder does not want to run a company that does bench work. The rigs are the
  differentiator; without them this is a software company competing with pyRheo and the
  solver vendors on their ground. That is a founder decision (`07_DECISIONS.md`, D3), and it
  should be made early.
