# Before China: priorities, 06 AUG to 31 AUG 2026

Written 2026-08-06. You fly the week of Aug 31, so this covers about three
weeks. The visa is in hand, so nothing here is waiting on the consulate.

## The one idea that sorts everything

You have 18 working days and roughly 11 of them already have a meeting in
them. That is not enough time to move twelve projects, and trying to is most
of why the list feels crushing.

But you are not short of time in general. You are short of one specific kind
of time. Three weeks in China costs you nothing on a desk task. It costs you
everything on anything that needs the bench, the Zetasizer, the E: drive, a
purchase order, or a person on campus.

So the sort is not by importance. It is: does this need Evanston? If yes, it
happens before Aug 28. If no, it goes on the plane.

Judged that way, most of the list turns out not to be urgent, and the part
that is urgent is short.

---

## Tier 1: needs the lab, and dies if it slips

### 1. Settle the AIBN number, then run TPM batch 1

Top of the list, because the 03 AUG meeting made synthesis the centre of the
program and because particles keep in a bottle. Synthesize before you fly and
you come home to material that is ready for characterization. Skip it and the
synthesis program loses the whole trip.

One number blocks the oven. The sphere-only sheet still says AIBN at
1 mg/100 mL. `TPM_BENCH_CARD_04AUG26.md` argues that is 20 to 67 times below
both van der Wel (67 mg) and Crothers (20 mg), and that its failure mode looks
identical to dead initiator. At 0.5 vol% TPM, 1 mg works out to about
0.2 wt%.

Running the sheet verbatim is a defensible test of whether the recipe is
viable on its own. But if the batch fails on 1 mg you cannot tell whether the
recipe was bad or the initiator was, which is the one thing a first pass is
supposed to settle. Ask Jairo or Michelle and close it. This is a
conversation, not a decision you have to make alone.

Ammonia and TPM are already settled at the 0.5 vol% pair. Size the batch by
SEM before believing either 700 nm or 800 nm, since the sheet's own note says
nominally identical recipes differ by up to a factor of two between labs.

### 2. Start the Zetasizer clock this week

Three separate things, and each carries lead time you cannot compress from
Shanghai: Keck-II access, a training refresher, and confirming the lab has the
cells and glassware or ordering them. Order consumables in week one and they
arrive while you are away. Order them in week three and they arrive after you
get back.

The card treats this as a step. It has been a blocker since 30 JUL.

### 3. Carbopol rig, one commissioning day

Lowest mass fraction plus the Newtonian control, in that order, treated as
commissioning rather than data. Your own stated purpose was dialing in the
joint low-speed and high-speed camera system.

The argument for doing this before you leave is not about the data. A rig you
left dialed in is worth roughly a week of September. A rig you left cold is
not.

If any of it starts counting as data, retain the batches or measure the flow
curve, because a wt% without a τ_y cannot be placed on the J axis afterwards.
And seed the control drop with the same 1 % w/w ground silica. If the control
is clear water, the tracer is confounded with the substrate.

### 4. Mount E:, and stage the laptop

Two jobs on one drive. Process the already-shot 50 mM LiCl cine, which is pure
re-analysis and fills the DLVO midpoint with no new experiment. And check
`E:\rheology\` while the drive is up.

Then the part that makes the trip productive: copy onto your laptop everything
the deferred desk work in Tier 4 will need. E: is removable and frequently
unmounted, and you will not reach it from another continent.

### 5. Launch the long compute on your last day

The loose end is L10 at θ = 20° and 45°. Machines run while you are on a
plane. Anything you queue on Aug 31 is three weeks of free progress, and it is
the cheapest item on this entire list.

---

## Tier 2: owed to people, and it gets worse with time

### 6. The replication README for Alex Davis

Open since 30 JUL. A collaborator is holding code he cannot run, on a paper
with an external deadline. Small job, and the one most likely to turn into an
awkward email while you are abroad.

Start from `SI_image_processing_concise.tex` rather than a blank file. It
needs the entry points, the expected input tree, which figure each script
produces, and the calibration constants that are not in the code: gap_row
29/33, px/mm 6.9034/9.2045, `min_gap_px` 22.091.

While you are in that repo, commit the two untracked files,
`build_velocity_field.py` and `mode_transition_analysis.py`, 747 lines between
them. The mode-transition work currently sits on an orphan branch whose only
commit is empty. That is a real chance of losing it, and ten minutes removes
the risk.

### 7. Jamming mini-lecture curriculum: give Michelle a date

Overdue since June. The honest move before a three-week absence is not to
draft it badly at 11 pm. It is to tell her at the Aug 10 one-on-one when it
lands. Pick a date in October and say it out loud.

### 8. The Aug 17 Drop Impact meeting

You created this on 05 AUG with Xiang Cheng and Michelle, and it repeats. A
new external collaborator's first meeting sets the tone. Decide this week what
you owe it, because it is 11 days out and no card names it.

---

## Tier 3: nearly finished, so ship it

### 9. Post the slime mold preprint

The best value per hour on the board. A complete 9-page manuscript, five
independent verification passes, headline reproduces exactly, zero undefined
references, held back by a line reading `\author{Prepared for Charlie}`.

Three steps: fill the author block, add the one-sentence Table 1 footnote on
the networkx Steiner/MST coincidence plus a pinned `requirements.txt`, then
post to `physics.soc-ph`. Maybe two hours.

Do it before you fly rather than after, so any arXiv hold resolves while you
are still reachable.

### 10. Two ten-minute items

Waive or confirm NU-SHIP. The window runs until 10/1, so it technically
survives the trip, but you will be abroad for a good chunk of it.

Ask Surbhi whether the gels you imaged on 30 JUL are the coumarin cure-depth
thread. One question, and if the answer is yes, that project has data nobody
has looked at.

---

## Tier 4: deliberately deferred to the trip

These are not dropped, they are scheduled. All desk work, all portable, and
moving them off the pre-trip list is what makes the pre-trip list achievable.

Paper 1's scalar-set decision, and the 18 band-C deposit annotations (0, 10
and 100 mM at `u0_3`, velocity spread 1.24 %, zero rebounds). One caveat: the
clicks are the irreplaceable part, but the downstream route is broken, since
`manual_overrides.csv` and `build_times.py` went with the deleted folder. Bank
the clicks and rebuild the consumer later.

Michelle's porous-diffusion estimate, meaning how long an 800 nm jammed
network retains structure against capillary solvent flow. It is a
back-of-envelope calculation that speaks directly to the late-tail result.

Frisbee, resuming the hardening brief at WP1. Still do not post.

Bidet, the three author decisions: venue, τ_c, Tier-B CFD.

Career items, including the portfolio repos and confirming work
authorization.

---

## What the three weeks look like

Week of Aug 10 is decisions and other people's blockers. Settle AIBN. Book SEM
and the Zetasizer refresher and order consumables. Write Alex's README and
commit the stray files. Post the slime mold paper. Ask Surbhi the one
question, and waive NU-SHIP. Bring the TPM plan and a real curriculum date to
Monday's one-on-one.

Week of Aug 17 is the bench. TPM batch 1 and SEM sizing, Carbopol commissioning
day, and the Xiang Cheng meeting on the 17th.

Week of Aug 24 is characterization and staging. Zeta and rheometry on the TPM,
Carbopol at the lowest wt% plus the Newtonian control, then mount E:, process
the LiCl cine, and stage the laptop.

Mon Aug 31: launch every long run you have, then go.

---

## Two things worth saying plainly

This list is shorter than it feels. Tier 1 is four bench jobs and Tier 2 is
three conversations. It reads as overwhelming because eleven projects each
carry a live next action, but nine of those projects need nothing from you
this month, and several are further along than they feel. The slime mold paper
is finished. The bidet paper is journal-ready and waiting on three decisions.
The frisbee preprint is drafted and explicitly on hold.

If the three weeks compress, protect TPM. It is the one item where the advisor
meeting, the §8 rheology exposure, and the timing of the trip all point the
same direction. Carbopol commissioning is second. Everything else can move.
