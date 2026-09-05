# 05 — Roadmap: phases, gates, budgets

Two resources, spent differently. Human hours (5–10 per week until the PhD ends) buy
decisions, bench measurements, and buyer conversations. Overnight agentic sessions buy
execution against gates written in advance. Each unit below names which resource it spends
and what test decides it. A failed gate is a recorded result, not a blocked phase.

The PhD end date sets the phase boundaries and is not known here; phases are written in
elapsed months from a start and should be re-anchored once it is (`07_DECISIONS.md`, D6).

## Phase 0 — Foundation (months 0–3). Software only, no company yet.

Goal: the stack exists as a repository with three layers running end to end on one fluid, and
one public case study demonstrates the chain.

| unit | resource | gate |
|---|---|---|
| 0.1 Create the `biorheo` monorepo; move `dst_fit` into `packages/identify` behind the generalised `Model` interface without changing its test results | overnight | `dst_fit` tests green from the new location; byte-identical figures |
| 0.2 L0 schema, loader, provenance validator; refuse any record without a resolved DOI or measurement manifest | overnight | a synthetic record with a fake DOI is rejected; a real one loads |
| 0.3 First registry records: Carbopol surrogates (own data), Guy–Hermes–Poon PMMA (already fitted in `dst_fit`), whole blood from two independent published sources | overnight, then founder review | each record's envelope matches the source figure's axis range; identifiability field written by L2, not by hand |
| 0.4 L1 generalised-Newtonian family with every yield-stress regularisation side by side, SymPy definitions, limits, spot-checks | overnight | `derive-and-check` outputs in `derivations/`; Papanastasiou, bi-viscosity and clamp agree on the flow curve where the regularisation is inactive and disagree by the predicted amount where it is not |
| 0.5 L1 structural family: port `bio_fluids_bridge` and `tier0` closures; L2 fits them from synthetic transients with known truth | overnight | recovery within the interval on 20 seeded synthetic sets; a deliberately degenerate protocol returns "degenerate" rather than a number |
| 0.6 First public case study: injectability of a yield-stress and a shear-thinning surrogate through a needle, closed form (Buckingham–Reiner) with full UQ from the L2 posterior | overnight for the code, founder for the write-up | matches a hand-computed point at three parameter sets; the report states which parameters the injection force is insensitive to |
| 0.7 Ten discovery conversations with people who buy or commission this work (formulation scientists, device engineers, regulatory consultants) | founder, 1 hr each | written notes; the question is whether anyone has paid for constitutive rigour or would; this is the kill test in `06_RISKS` R1 |

Exit: the chain runs from a registry record through identification to a closed-form prediction
with uncertainty, on one fluid, with a public write-up. Discovery notes exist.

## Phase 1 — Bridge and credibility (months 3–7)

Goal: a model goes into a real solver and comes back; one existing study is wrapped in a V&V 40
package retroactively.

| unit | resource | gate |
|---|---|---|
| 1.1 OpenFOAM exporter for built-in GNF models plus the Couette round-trip harness in a container | overnight | round trip within 2 % on a Carreau fluid; the induced error of the Papanastasiou export of an HB fluid is reported and matches the L1 prediction |
| 1.2 `libbiorheo` v0: structural closure with state transport as an OpenFOAM viscosity model | overnight, multiple | reproduces the L1 0-D transient in a homogeneous-shear case to 1 % |
| 1.3 Basilisk exporter using the existing `getenv` convention; round trip on the axisymmetric Couette | overnight | same bar as 1.1 |
| 1.4 L5 skeleton: the objects, the evidence manifest, the report generator | overnight | the manifest for a study with no evidence lists every factor as missing; nothing renders as satisfied by default |
| 1.5 Wrap the bidet grid study in L5 retroactively | overnight, founder review | the manifest reproduces the study's observed orders from the stored fields and flags the out-of-band design-model evaluation as an applicability gap |
| 1.6 Second case study: nasal or airway spray plume with a mucus simulant, Basilisk free surface, benchtop imaging from the existing rig | founder bench day plus overnight | drop-size distribution from imaging with audited calibration; simulation versus bench compared with paired uncertainty |

Exit: two solvers round-trip; one V&V 40 manifest exists for a real study; two case studies public.

## Phase 2 — Bench and first revenue (months 7–12)

Goal: the rigs are documented products; first paid work.

| unit | resource | gate |
|---|---|---|
| 2.1 Syringe–needle force rig build sheet, calibration manifest, pipeline; commissioning on glycerol against the Newtonian closed form | founder bench, overnight analysis | measured force within the stated uncertainty budget of Hagen–Poiseuille over a decade of rate |
| 2.2 Same rig on Carbopol and a xanthan solution; L2 identifies from the rheometer, L6 predicts the force, rig measures it | founder bench, overnight | prediction and measurement agree within combined uncertainty, or the disagreement is traced to a named regularisation or slip effect and reported |
| 2.3 Open benchmark suite v0: three geometries (needle, FDA nozzle, impinging jet), three fluids (Newtonian, Carreau, HB), reference solutions with tolerances | overnight | an independent OpenFOAM run by someone other than the author reproduces the reference within tolerance |
| 2.4 Registry to 25 records across blood, mucus and simulants, mAb viscosity curves, IDDSI thickeners, synovial fluid | overnight, founder review of every record | `lit-verify` resolves every DOI; every envelope checked against its figure |
| 2.5 Two paid feasibility studies at the consumer tier (splash, absorbency, or dispensing), fixed price, using W4 assets | founder sales, overnight execution | delivered inside the fixed scope; a client-readable 90-second summary; the client is asked whether the identifiability verdict was worth paying for |
| 2.6 Thrombosis paper: run the named discriminating test on digitised published occlusion data | overnight, founder | the exponent is either resolved or the data is shown too sparse to resolve it, and either result is written up |

Exit: a rig with a stated uncertainty budget; a benchmark suite others can run; two invoices;
one flagship manuscript submitted or the reason it cannot be written down.

## Phase 3 — Company (from PhD end)

- Incorporate; professional liability and errors-and-omissions insurance before any device-tier
  work; a minimal quality manual written against the V&V 40 manifest so the QMS question has an
  answer when a client asks.
- Subcontract into two established consultancies for device-tier work; the credibility manifest
  is the pitch.
- First drug-development-tier engagement in W1 with client-supplied rheometry.
- Licensing tier of the registry and `biorheo` pro (identifiability, exporters, credibility);
  the PolyUMod model.
- Viscoelastic family and the poroelastic rung, driven by the first W2 or depot engagement.

## Budget split, per week, during Phases 0–2

| resource | spend |
|---|---|
| founder hours | 2 on reading and approving overnight output; 2 on bench or calls; 1–3 on writing that only the founder can do (case-study prose, discovery notes, gate verdicts) |
| overnight sessions | one unit per night, two or three subagents at most, launch compute first; a unit that will not fit a night is split, never stretched |

## Rules that carry from `rheo_products/SPEC_2026-08-27.md` §0

Gates written before runs and enforced in code; subagent budget is the limiter, not compute;
prefer agents that return a conclusion to agents that edit; long compute launched first; done
means tests green, gate verdict recorded (pass or fail), and one artefact a human can read in
90 seconds.
