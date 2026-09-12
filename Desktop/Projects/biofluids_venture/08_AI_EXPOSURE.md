# 08 — AI exposure: which layers of this venture survive a frontier model

Date: 2026-09-12. Written in response to the founder's own challenge ("it feels like not a
lot of progress is being made and I'm worried the premises aren't AI-proof as a product or
service"). This file grades the venture's premise against that challenge and says which
parts of `05_ROADMAP.md` / `00_SYNTHESIS.md` change and which do not.

It is adversarial by intent. It does not try to rescue the plan.

---

## 0. Why the existing research does not already answer this

`00_SYNTHESIS.md` §2 F2–F3 addresses AI, and §3's ranking table carries a **durability**
column. But every entry in that column is an *assertion*, not a test:

- rank 1: "high: the verdict and CI are what AI cannot fake"
- rank 5: "highest — provenance + data compounds"

Nothing was run to check either. F3's claim that physics-AI surrogates are "not a threat
for 5+ years" is about *surrogate solvers* — a different product from the one being sold.
The question the founder is asking is narrower and harder: **can a frontier model plus a
competent engineer produce this deliverable?** That question has an answer on disk, and
nobody looked at it.

## 1. The replication test already ran, and the venture ran it on itself

`biorheo` — L0 registry, L1 constitutive, L2 identification with profile-likelihood
intervals and an identifiability verdict, L3 OpenFOAM export, L6 injectability app —
was built by agents.

Measured, not argued (`git log` in `Desktop\Projects\biorheo`, 2026-09-12):

| quantity | value |
|---|---|
| commits touching the package | 91 |
| first commit → last commit | 2026-09-05 → 2026-09-08 (**4 calendar days**) |
| Python in `src/` + `tests/` | 12,765 lines |
| fast tests green | 564 (+11 slow, 1 xfail) |
| founder lines of code written | **0** (`repo-is-agent-written`) |

So the honest reading of the AI-exposure question is not a forecast. **The differentiator
in `00_SYNTHESIS.md` §1 — "the artefact nobody sells, in any industry, at any price" —
was built from nothing in four days by the very technology it is supposed to be defended
against.** Anyone with a Claude/GPT subscription, a rheology background and four days has
the same artefact. The absence of sellers (F1) is therefore weak evidence of a moat; it is
equally consistent with the artefact being *newly cheap* rather than *hard*.

## 2. Layer-by-layer grading

The replication test for each layer is named so it can be disputed.

| layer | effort spent | replication test | survives? |
|---|---|---|---|
| L1 constitutive laws (Casson, Cross, Carreau–Yasuda, HB, power law, structural) | large | Write them from the textbook definitions with `pint` units | **No.** Public formulae. |
| L2 fitting + profile-likelihood CIs + identifiability | largest | Ask a model for `scipy.optimize` + profile likelihood on a flow curve | **No.** `dst_fit`'s own defect (item 13) was found and fixed by an agent in one session. |
| L3 OpenFOAM `transportProperties` export | medium | Write a dictionary writer against the v2606 source | **No.** It is a file format. The *name collisions* (`CrossPowerLaw`'s `m`/`n` reversal) are a gotcha, not a moat — and they are now written down publicly in this repo. |
| L6 closed-form needle force + UQ | medium | Hagen–Poiseuille / Buckingham–Reiner with a sensitivity sweep | **No.** And worse: F4's own note says a 2018 Roche model already gets <12 % error, and buyers own a Hagen–Poiseuille spreadsheet (01 §Summary 8). |
| L0 registry, 125 records / 10 licensable | several items | Re-run `scripts/datacite_candidate_search.py` | **Weakly.** That one script found **107 licensable dataset DOIs** in an afternoon. The catalogue is a head start measured in days, not a barrier. |
| The verdict layer — 17 failed gates never loosened, item 23's caption-vs-fit-table catch, item 27's compressibility factor | large | see §3 | **Contested. Read §3 — the repo's own record refutes the optimistic answer.** |
| A named person who signs a V&V 40 / EDDO credibility argument | **zero** | — | **Yes.** Not compressible. |
| Rheometry on the client's own proprietary formulation | **zero** | — | **Yes.** No model puts a sample on a rheometer. |
| Deciding *what* to measure, and the client relationship that reveals it | **zero** | — | **Yes.** This is F9's whole finding: buyers do not know they need this. |

**The pattern is the finding.** Approximately all effort to date has gone into the rows
graded *No*. The three rows graded *Yes* have had none.

## 3. "AI cannot fake the verdict" is refuted by this repo's own history

This is the claim the plan's durability rests on, so it gets its own section.

The `biorheo` record is genuinely unusual: seventeen pre-registered gates failed or were
mis-specified and **not one was loosened**; a published anchor was caught being graded
against a figure caption instead of the paper's Table 2; a missing compressibility factor
(median 1.574) was caught by an admissibility check rather than an argument; a "415 %
parameter move" was retracted when it turned out to be 9.8e-19 over 2.4e-19.

That is exactly the judgement the plan sells. **And agents produced every one of it.**
Four reviewer passes found item 23's defect; two found it independently. The numerics pass
re-derived the corrected number from a from-scratch fit. No human was in that loop.

So the defensible statement is *not* "AI cannot fake a verdict." It is narrower:

> An adversarial, pre-registered, gate-keeping process produces verdicts a casual user of
> the same model does not get. The process is the asset. The process is also written down
> in `~/.claude/CLAUDE.md`, `lessons.md` and every `prereg.md` in the repo — i.e. it is
> copyable prose, and its value decays as the default behaviour of frontier models moves
> toward it.

That is a real advantage with a shelf life, not a moat. Treating it as permanent is the
premise the founder is right to distrust.

## 4. The larger threat is not AI

Two non-AI risks are better evidenced than the AI one and are already in the research,
unweighted:

1. **The buyer's existing model is good enough.** 01 §Summary 8: "the physics is 80 %
   solved by a 2018 Roche model with <12 % error" and the buyer "already owns Instron/Zwick
   rigs and a Hagen-Poiseuille spreadsheet." An error bar on a number they already trust is
   a hard sell. This is the risk R1 question 4 exists to measure.
2. **No buyer was ever found.** §8 flag 2: *no buyer quote anywhere saying "we wanted an
   uncertainty on the constitutive law."* Four days of code do not change that, and 28
   closed items have not tested it once.

## 5. What this changes

**Unchanged.** The wedge (W1 injectability), the refusal to sell hours, simulants-only, no
bench before an invoice, and the ranking itself. None of that rested on the code being
defensible.

**Changed.**

1. **Stop counting code as progress.** The ledger to watch is calls held, proposals out,
   invoices paid. Current values: 0, 0, 0.
2. **Re-rank on the durable rows, not the built rows.** Rank 2 (V&V 40 credibility
   package) is now arguably *above* rank 1, because its value is a signature and a
   regulator-facing argument — the one row graded *Yes* that a solo founder can supply
   immediately. Rank 1's value is a calculation. This is a recommendation, not a decision;
   it is D-A and it is the founder's.
3. **Teaching moves from channel to product.** F9 already says teaching is the only channel
   with repeated evidence. §2 adds a reason: a short course sells the *process* — how to
   know your extrapolation is wrong — which is the asset with the shortest shelf life, so
   it should be monetised early rather than late. Month 7–8 is too late; it belongs in
   month 2–4.
4. **Every agent work order should be judged by which row of §2 it feeds.** The current
   `WORK_ORDERS_2026-09-08.md` first wave (B1/B3/B5) feeds rows graded *No*. That is not
   worthless — B3 in particular protects against shipping a wrong anchor to a client — but
   it must not be mistaken for progress toward revenue.

## 6. What would falsify this memo

Stated so it can be checked rather than believed:

- **F-1.** Three or more of the ten R1 calls say they have paid, or would pay, for the
  high-shear/identifiability layer specifically. That is the existing R1 bar; passing it
  refutes §4.2 and substantially weakens §1, because it would establish that the artefact
  being cheap to build does not make it cheap to *trust*.
- **F-2.** A named competitor ships a fit with a CI and a solver round trip within 12
  months. That confirms §1 and kills rank 5.
- **F-3.** A client asks for the OpenFOAM export and the record to keep (rank 5's own kill
  test). That would show the catalogue has value independent of who computed it.
- **F-4.** One paid short course sells ≥10 seats. That confirms §5.3.

## 7. The one-line verdict

The code is the least defensible part of this venture and has absorbed essentially all of
the effort; the defensible parts are a signature, a rheometer and a relationship, and have
absorbed none. That is a correctable allocation problem, not a dead premise — but it is
not correctable by writing more code, and the gate that would settle it (R1) has been
sitting unrun since 2026-09-07.
