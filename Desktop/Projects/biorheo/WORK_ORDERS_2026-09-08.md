# Parallel agent work orders — biorheo and Drop Impact — Simulation

Written 2026-09-08. **This is a dispatch sheet, not a plan of record.** Each package below
is sized to be handed to ONE agent in a fresh context with no other package running on the
same files. The `owns` line is the contract: two agents must never be given packages whose
`owns` sets intersect.

Copy of this file lives in `C:\Users\compl\Github\Drop Impact\WORK_ORDERS_2026-09-08.md`.

## How to dispatch one

Paste the package's **brief** into a new agent. Every brief inherits, without restating:

> Read `~/.claude/CLAUDE.md` and `~/.claude/lessons.md` first and follow them. Pre-register
> before running. Never move a bar that fails — add a stricter test beside it and keep the
> original number. Never tune to a picture. Say "unverified" rather than guessing.
> Reproduce the claim you are extending BEFORE extending it. Report what failed.

**Model choice.** Opus where the work is a judgement about what the right question is, or
where a wrong answer would be written into a client-facing page or a solver. Sonnet where
the question is already sharp and the work is execution, search, or audit against a rule
that is written down.

**The one collision that has bitten this repo repeatedly**: appending to a shared list
(`RECORD_IDS`) silently re-decides closed gates that iterate it. Any package that adds a
record must grep every consumer of the name before appending. Only ONE record-adding
package may run at a time.

---

# Part 1 — biorheo (the consulting company)

Repo `C:\Users\compl\Desktop\Projects\biorheo`, branch `portfolio-audit-2026-08-29`.
Ledger `NEXT_STEPS.md`; items 1–28 closed except 17's founder decision (D5).

**Shared-file rule.** `NEXT_STEPS.md`, `OFFER.md` and `README.md` are written by the
INTEGRATOR (you), not by package agents. A package agent writes its own item file under
its own directory and reports the ledger paragraph it wants added; you paste it. This is
the cheapest way to run four agents on one repo without a merge conflict.

---

## B1 — Turn item 28 into product: the device budget the client actually needs
**Opus.** ~3–4 h. **Highest commercial value on this list.**

Item 28 (`case_studies/device_friction/`) measured that a formulation screen stops being
decidable once syringe-to-syringe friction repeatability exceeds **s\* = 0.1397 N**, and
that **nobody publishes a friction SD for any device** — the literature reports means.
So the quantity that decides whether a client's screen can distinguish two formulations is
unmeasured, and measuring it is the cheapest bench task in the whole programme: N
empty-barrel injections on one syringe lot.

Build two things.

1. **A computed device budget on the deliverable page.** `apps/report.py` should be able to
   answer, for the client's OWN cohort and their own decision, *"your syringe friction must
   repeat to better than X newtons for this comparison to hold"*. That is `s*` recomputed
   per study, not item 28's number reused — item 28's `s*` belongs to 58 mAb formulations
   at one flow rate, and a bug's blast radius is a property of the caller.
2. **A bench protocol page** for the measurement that fills it: empty-barrel injection at
   the study's own plunger speed, how many syringes, what to record, and what the result
   licenses you to say.

**Pre-register**, before writing either: what would make the budget section WRONG to ship.
At least these — (a) a cohort whose forces are already indistinguishable at `s = 0`, where
the budget is undefined rather than zero and must say so; (b) a single-formulation study,
where there is no pair and the section must refuse rather than print a number; (c) the
screen path (`build_from_screen`), where the interval is a measurement repeatability and
the device term is the SAME KIND of quantity, so adding it in quadrature may be
double-counting — decide that explicitly and say which.

**Do not** default a friction value on the caller's behalf, and do not let the page say
"the force" about a number that omits the device. The five existing pages each carry the
sentence *"Syringe friction is a property of the device, not the fluid"* — that sentence is
correct and must stay; what changes is that the page now says what the friction has to be
repeatable TO.

**owns**: `src/biorheo/apps/report.py`, `case_studies/deliverable_example/`,
`case_studies/device_friction/` (extending), `tests/test_report*.py`,
`tests/test_device_friction.py`.

---

## B2 — Close the wedge-fluid data gap: a licensable high-concentration protein FLOW CURVE
**Opus.** ~4–6 h. The registry's standing open item.

The registry is 125 records, 10 licensable, and **not one of them is the wedge fluid at
device rates**. Item 25 demonstrated the chain on 58 real mAb formulations but they are
single viscosities — one point cannot determine even a one-parameter model, so no record
exists. What is needed is a licensable protein flow curve with more than one point that
reaches rheometer-to-needle rates.

Start from `data/candidate_search_2026-09-08.md` (**107 licensable dataset DOIs**, from a
route that works — Zenodo's REST API returns 504; use DataCite). Three named candidates
are unopened:

- `10.5281/zenodo.15469989` — PNAS, *Electrostatics and viscosity are strongly linked in
  concentrated antibody solutions*. Ships one `.tar.gz`, never opened.
- `10.17877/resolv-2025-mht4btcr` — *Predicting the Dynamic Viscosity of
  High-Concentration Protein Solutions*.
- `10.6084/m9.figshare.32084472` — *Disentangling mucus rheology and transport efficiency
  in human airways*. **A second biofluid family**, and worth taking on those grounds even
  if it does not close the protein gap.

**Rules that are not negotiable**: do not digitise a figure; do not take a table of fitted
parameters without the curve behind it; resolve the licence from the deposit's own metadata
and record it. Before appending any record, **grep every consumer of `RECORD_IDS`** — this
has silently re-decided closed gates twice, and the suite stayed green both times because
every verdict was still true, on a different question.

**Grade against the source's own FIT TABLE, never a figure caption or an abstract.** That
mistake has now been made twice here (item 22, and the Singh polydispersity comparison).
And compare parameter sets **as viscosity curves over the measured range**, not parameter
by parameter — `k` and `n` trade off, so a deviation in one alone is uninterpretable.

**owns**: `data/`, `src/biorheo/registry/records/`, `scripts/datacite_candidate_search.py`,
`tests/test_registry*.py`. **Exclusive**: no other package may add a record while this runs.

---

## B3 — Audit every external anchor against its source's fit table
**Sonnet.** ~2 h. Read-only until the findings file is written.

Item 23 found the repo's then-strongest external anchor had been graded against the
paper's **figure-caption** values — what the authors fed a simulation, rounded and shared
across conditions — while the same paper carried a Table 2 of fit results, and all four
values disagreed. Two of the four were the right number against the WRONG ROW, which is
how a total mismatch reads as a partial one.

That check has been applied to one anchor. Apply it to **every external comparison in the
repo**: for each, name the source, name exactly which table/figure/caption the compared
number came from, and say whether the source publishes a fit table that was not used.
Candidates include the DMEM/FBS power-law records, the Carbopol-in-propylene-glycol shift
anchor, the Koponen foam pipe anchor, the Horner blood Casson parameters and the mAb
concentration law.

Deliverable: `data/external_anchor_audit_2026-09-08.md`, one row per anchor, with a
verdict column: GRADED AGAINST A FIT TABLE / GRADED AGAINST A CAPTION / NO FIT TABLE EXISTS
/ UNVERIFIED. Do not fix anything; rank what needs fixing.

**owns**: `data/external_anchor_audit_2026-09-08.md` only. Read-only elsewhere.

---

## B4 — Size the omitted terms for the mAb cohort, and reconcile two routes that disagree
**Sonnet.** ~2–3 h.

`blood_injectability/results.md` sizes the entrance (Hagenbach) loss at **0.43–0.99 % of
`dp`** for the 27G geometry at 1 mL/min. The literature survey run 2026-09-08
(`data/injection_force_terms_lit_2026-09-08.md`) reaches **1.4–4.6 %** for the same
`L/D ≈ 60`, by a different route (capillary-rheometer end-correction coefficients, Boger,
Gupta & Tanner 1978, doi:10.1016/0377-0257(78)80006-8). **Those two numbers disagree by
about 3× and nobody has reconciled them.** Do that, and then size the term for the mAb
cohort, which has never had it done.

Both routes must be stated as what they are: the first is a kinetic-energy (Hagenbach)
estimate at the mean velocity; the second is a shear-thinning end correction whose
coefficient is a function of the flow index. Say which one the mAb formulations' own
Newtonian classification licenses.

Also open and separate: the **hub/taper** contraction between barrel and cannula, which no
paper found in the survey isolates. If it cannot be sized from published data, that is the
answer — say so and say what measurement would settle it.

**owns**: a new `case_studies/omitted_terms/`. Do NOT edit `blood_injectability/` or
`mab_injectability/`; their artifacts are pinned by tests and by item 28's C1 reproduction.

---

## B5 — Referee items 24, 25 and 26
**Opus, `physics-referee` subagent type.** ~2 h. Read-only.

Items 21–23 and 27 have been refereed. Items 24 (the temperature-shift external anchor),
25 (the mAb wedge fluid) and 26 (the screen deliverable page) have not. Attack, in
particular: item 24's `a_fit = 0.13984` against `a_pred = 0.11765` — is a ratio of 1.1887
against a bar of 1.30 a pass or a bar too wide to fail? Item 25's concentration-normalised
7.33× — does the exponential concentration law hold over 159–212 mg/mL for THESE
formulations, or is the normalisation extrapolating? Item 26's screen page — what does it
claim that a single viscosity cannot support?

**owns**: one report file. Read-only elsewhere.

---

## B6 — Re-grade `OFFER.md` against everything found since it was written
**Sonnet.** ~1–2 h. **Run this LAST**, after B1–B5 report.

Every number in `OFFER.md` is recomputed from an artifact by `tests/test_offer_claims.py`.
Item 28 changes what the page may claim about a force: any RATIO of two forces is a
property of the friction being zero, and the invariant statement is the DIFFERENCE. Check
the page for ratio claims, and check the capability-gap list against B1–B5's findings.

**owns**: `OFFER.md` (in `Desktop\Projects\biofluids_venture\`),
`tests/test_offer_claims.py`.

---

# Part 2 — Drop Impact — Simulation (Wyart–Cates f-field DNS)

Repo `C:\Users\compl\Github\Drop Impact`, branch `overnight/2026-06-29-autorun`.
Solvers in `basilisk-drop-impact\`; Basilisk in WSL, **serial only** (OpenMP hangs on this
machine), 22 cores, parallelism via concurrent background subshells with
`OMP_NUM_THREADS=1`.

**Compute budget is the scarce resource, not agents.** Two campaigns already run as of
2026-09-08 (`armb_hbath_2026-09-08`, `armc_ladder_2026-09-08` + `armd_ladder_2026-09-08`),
holding 18 of 22 cores. **Do not dispatch S1 or S5 until those report.**

**Standing constraints for every compute package.** New campaign directory only; never
`rm -rf`; never edit a shell script while bash is running it; LF line endings only (a CRLF
kills a WSL script on a bare `do`); never chain a command onto a backgrounded launch line;
every poll loop echoes the value it actually read; `timeout` on every cell; a source-md5
guard that refuses to start if the solver moved; per-cell `RUN_STATUS.json` with rc,
`wall_s`, `stop_reason` and `series_md5`. `rc = 0` proves *a* run finished, not *your* run
— check the banner encodes the parameters you asked for.

Wall-clock, measured: `wc` L9 ≈ 0.03–0.6 h, L10 ≈ 0.6–5 h, L11 ≈ 1.5 h; `hb` bath cells at
L9 ≈ 0.05–0.15 h, an order of magnitude cheaper than arm A.

---

## S0 — Refresh the Status card (do this first; it costs nothing and unblocks everyone)
**Sonnet.** ~1 h. No compute.

`Status\Drop Impact — Simulation.md` §3 is dated 2026-08-31 and **two of its five priority
items are already done**: arm B's clean J→0 control was built, run and graded that evening
(`ARMB_CONTROL_2026-08-31.md`), and the "face hoop term wrong on radial faces" defect was
fixed 2026-08-09 in all five solvers that carry the term (`HOOPFACE = 1` default-on). The
card's own §2 banner says it was three weeks stale at 08-31, and four September sessions
have run since.

Rewrite §2 and §3 from the `NIGHT_REPORT_*` series and `ARMB_CONTROL_2026-08-31.md`,
newest first. Keep §6's load-bearing numbers with their provenance. Mark what is genuinely
still open — and note that the hoop item's SECOND half is open even though its first half
is fixed (see S1).

**owns**: `Status\Drop Impact — Simulation.md`.

---

## S1 — The hoop default, and the half of the fix nobody flipped
**Opus.** ~1 h for the gate + re-run cost. **Best value-per-hour on this list.**

The radial-face defect is fixed: `0.5*(u.y[]+u.y[-1,0])/y` is the radial velocity on
x-faces and the AXIAL velocity on y-faces, inflating γ̇ by up to √3 for pure axial shear
and growing under refinement — and the cell-centred `d33c` fix is present with
`HOOPFACE = 1` default-on in `hb_drop_on_bath.c`, `wc_drop_on_bath.c`, `rheo_substrate.c`,
`sphere_impact_hb.c`, `sphere_settle_hb.c`.

**What is still open is that `HOOP` itself defaults to 0** (`hb_drop_on_bath.c:335`,
`wc_drop_on_bath.c:128`), so the hoop term is omitted from the face invariant on every
default bath run. The direction of the bias is known and it points at the measured
quantity: omitting `D_tt` under-estimates γ̇ → over-estimates `τ_y/γ̇` → biases the crater
**shallower, on the axis, exactly where it is read**. Four further solvers
(`power_law_drop.c`, `stf_drop_on_bath.c`, `rheo_matrix.c`, `hb_drop_on_bath_papa.c`) have
no hoop term and no knob at all.

Three-tier gate, pre-registered before any run:
- **A** — byte-identity at `HOOP=0`, to prove the change is inert where it should be. Read
  it on a CONTINUOUS channel, not on a grid-quantized one: this project has twice called
  two different runs identical because `beta`, `h/D0` and `jam_frac` cannot express a 1 %
  difference that `umax` shows plainly.
- **B** — a manufactured velocity field, checked for direction symmetry on BOTH face
  families. **This is the gate that can only pass if the fix is real**; the old code
  returns different answers on x- and y-faces for the same physical field.
- **C** — a refinement ladder in which the OLD error grows as `1/Δ`. A defect that grows
  under refinement is the one a convergence study will otherwise attribute to the mesh.

Then decide, and state the consequence plainly: **flipping the default invalidates the back
catalogue.** The honest options are (i) flip and re-run only the cells whose numbers are
still quoted, (ii) leave the default and mark every bath crater number as carrying a known
one-sided bias of measured size. Do not flip silently.

**owns**: the five solvers' `HOOP` default, a new gate campaign, `SOLVER_EXPOSURE_REGISTER`.
**Blocked by**: the two campaigns now running (core budget).

---

## S2 — The augmented-Lagrangian port (AL3)
**Opus, staged.** Weeks. **The item that unblocks the most.**

One change closes three separate dead ends: arm D's HB-drop cannot start (a free-flying
drop has γ̇ ≡ 0 everywhere, so a regularised H-B drop is a ceiling-viscosity plug bounded by
air across one cell — `rc = 136`, Poisson residual 8e48 on the FIRST step, and a
pre-registered 4-cell rescue showed neither `MUMAX`, nor τ_y/K, nor CFL fixes it); arm B's
ceiling never converges; and the L10/L11 `dt`-collapse is a pressure-solve failure that
halving the CFL does not fix (`dt` reaches 7e-09; the card's own verdict is that refinement
above L10 needs a solver change, not a longer queue).

Already done and not to be redone: the algebra, elliptic and saddle spikes all gate
(7/7, 6/6, 7/7); `AL_PORT_DESIGN_2026-09-01.md` is the spec; `r* ≈ max(10, 5·Bn)` is
priced, arrested end ≈ 3200 iterations; `AL2f` showed `r` can stay uniform and `Λ` needs no
advecting; τ_y MAX at corners beats MEAN by 3.4 % at 4.8× cost.

**Do NOT port dual-FISTA.** It was written up as "not optional above Bn ≈ 2" on a
measurement taken at one value of `r`; swept over `r`, FISTA caps at every `r > 1` and never
converges at Bn = 20, while plain ALG2 at its own `r*` wins 4.4×/1.0×/3.1×. The dual step
is `2r` and Nesterov is stable only below a step-size ceiling — momentum and large `r` do
not compose.

Stage it with go/no-go gates: (1) the moving VOF interface, which is the untested part;
(2) a single arm-D cell that must start from rest; (3) one L10 cell that must clear the
`dt`-collapse. Stop at the first stage that fails and report, rather than pushing through.

**owns**: `AL*` sources and specs, a new solver variant. Does not touch production solvers
until stage 3 passes.

---

## S3 — Audit the stranded worktree before rescuing anything from it
**Sonnet.** ~2 h. Read-only.

`worktree-sim-dev-session` carries 11 commits including an air-cushion refutation, an L12
rung and an F1 figure, plus a `wc_drop_impact.c` diff that conflicts with `SNAP_MODE`.

**The first move is `git log -- <path>` on every file it claims to add.** A directly
analogous rescue in a sibling project copied a 400-line module out of an orphan worktree
into master, where a 511-line version had existed for two months — and had the copy landed
it would have reversed the sign of the result by duplicating the data. **An orphan can be
the WORSE copy.** Produce a per-file verdict (ORPHAN IS NEWER / MASTER IS NEWER / GENUINELY
NEW / CONFLICTS) and a merge plan. Merge nothing.

**owns**: an audit report. Read-only.

---

## S4 — `wc_drop_on_bath.c` has never been built in WSL
**Sonnet.** ~1 h.

The repo copy is **4 863 bytes newer** than the newest WSL build (`thetafill_2026-08-18`)
and has never been compiled there. Every other solver is md5-identical between repo and its
newest build. Compile it, run one cell that reproduces a recorded `thetafill` result, and
report what the 4 863 bytes changed and whether any number on record moves.

**owns**: a build/repro campaign for that solver only.

---

## S5 — The dilatancy-ceiling unit defect, and then a finer ladder
**Opus.** ~4 h + compute.

`DILCEIL_DECISION_MEMO_2026-08-19.md` records a blocking defect in the 10 AUG launcher:
**ρU₀² = 7657 Pa against the solver's own 4800 Pa, a 1.595× error**, which must be fixed
before exposure E6 (`DIL_CEIL` resolved only to 1.773×) can be scored on a finer ladder.
Fix it, then re-score.

Read the existing verdict before re-opening the physics: the dilatancy sweep (`a1dil`) is
already scored, and **D2 is a NO** — the ratio rises 0.6370 → 0.6628 but never approaches
its declared 0.70, and a 10× amplitude sweep moves it 4 %. My own D1 prediction was refuted
(the term bites hardest at LOW U₀). So this package is about the CEILING's resolution, not
about re-litigating whether dilatancy answers referee objection #1. It does not, and there
is no Darcy velocity field or pore-pressure equation anywhere in the suite — that is a
solver project, not a campaign, and should stay an acknowledged exposure.

**owns**: the `a1dil` launcher family and its analysis. **Blocked by**: core budget.

---

# Suggested dispatch waves

| wave | biorheo | drop impact | why |
|---|---|---|---|
| **now** | B1 (Opus), B3 (Sonnet), B5 (Opus) | S0 (Sonnet), S3 (Sonnet), S4 (Sonnet) | no file collisions, no compute contention; S0 unblocks every later sim package |
| **next** | B2 (Opus), B4 (Sonnet) | S1 (Opus) | B2 is exclusive on records; S1 needs cores freed by the two running campaigns |
| **then** | B6 (Sonnet) | S5 (Opus) | B6 needs B1–B5's findings |
| **long** | — | S2 (Opus, staged) | weeks; run it alongside everything else, it touches no production solver until stage 3 |

**Founder-only, and no agent substitutes for it**: read
`biofluids_venture\docs\market_research_2026-09-06\00_SYNTHESIS.md` §1/§5/§6/§7 and
`OFFER.md`, then book the first five discovery calls. The research is unambiguous that
R1 — ten calls — is the binding gate on this venture, and none of the packages above moves
it.
