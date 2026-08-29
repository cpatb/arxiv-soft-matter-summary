# Skill Observation Log

Observations captured during task-oriented work. Each entry identifies a
potential skill improvement or new skill opportunity.

**Status key:** OPEN = not yet actioned | ACTIONED = skill updated/created |
DECLINED = user decided not to pursue

---

### Observation 1: Driver docs may reference files that don't exist — verify repo state before executing

**Date:** 2026-07-25
**Session context:** Executing FINISH_BIDET_JOURNAL.md end-to-end (CFD-validate a bidet paper).
**Skill:** New skill candidate: execute-driver-doc (or a checklist addition to any "follow this brief" workflow)
**Type:** open-source
**Phase/Area:** Intake / pre-flight

**Issue:** The driver doc named `bidet_optimization_1.CORRECTED.tex` and a
`CORRECTION_DOSSIER` as the canonical base; neither existed. The newest real
manuscript was `bidet_optimization_1.REVIEWED.tex` (blockers still flagged inline)
and the dossier was dated differently. Proceeding on the doc's assumptions verbatim
would have failed or fabricated a base. Verifying actual repo state first, adopting
the real newest file, and surfacing the discrepancy in the final report was the
correct move.

**Suggested improvement:** Any "follow this brief/driver doc" workflow should begin
with an intake step that verifies every referenced input file/path exists, and on
mismatch adopt the closest real artifact + record the deviation, rather than trust
the doc. Fabricating a missing input is never acceptable.

**Principle:** A brief describes intent, not ground truth. Reconcile the brief against
the actual repository before acting; surface deviations instead of papering over them.

---

### Observation 2: Adversarial review subagents catch verifiable bugs a single editing pass introduces

**Date:** 2026-07-25
**Session context:** Same session; a subagent edited a 1300-line manuscript, then three
independent review subagents (validation/postdoc/professor) reviewed it.
**Skill:** review (preprint/code review workflows)
**Type:** open-source
**Phase/Area:** Verify stage

**Issue:** The editing subagent produced a clean-compiling manuscript, but the postdoc
review tier caught three verifiable errors it had carried/introduced: a sensitivity
exponent taken from the superseded (laminar) closure instead of the operative (CFD)
one — inverting a robustness claim (tau_c^{1/3} vs tau_c^{1.43}); a load-bearing number
quoted at the wrong operating point (5.9 Pa at Re=1e4 ≈ 40 kPa, not the stated 70 kPa);
and an algebra factor-of-2 in an exponent. A "compiles + reads fine" pass would have
shipped all three.

**Suggested improvement:** For any generated technical artifact carrying numbers/algebra,
add an independent adversarial verify pass whose explicit job is to recompute exponents/
units/load-bearing numbers from first principles — not just proofread. Give the verifier
the source-of-truth data files so it can cross-check, not just the prose.

**Principle:** Generation and verification should be separate passes with separate context.
The verifier's mandate must be "recompute and refute," not "read for sense" — clean
compilation and fluent prose are not evidence of numerical correctness.

---

### Observation 3: Compute a nominal case to shake out a pipeline before parallelizing a sweep

**Date:** 2026-07-25
**Session context:** Building an OpenFOAM axisymmetric impinging-jet sweep from scratch.
**Skill:** New skill candidate: cfd-parameter-sweep (reproducible OpenFOAM case generation)
**Type:** open-source
**Phase/Area:** Execution

**Issue:** The wedge-apex mesh needed two debugging iterations (unterminated header comment;
degenerate axis faces needing a shared collapsed vertex). Doing this on a single nominal
case first — then parallelizing the 13-case sweep only once the generator, solver, post-fit,
and figure all worked end-to-end — avoided 13× the debugging churn. A Python case-generator
(one function per dict) made the sweep reproducible and the fixes single-point.

**Suggested improvement:** A CFD-sweep skill should mandate: (1) a parametric case generator
(no hand-copied case dirs), (2) shake out one nominal case fully (mesh→solve→post→figure)
before any fan-out, (3) report per-case convergence honestly and flag any case that hit an
iteration cap rather than the residual target.

**Principle:** De-risk the whole pipeline on one instance before scaling; scale only a
working pipeline. Parametric generation makes fixes single-point and the sweep reproducible.
