# Item 14 pre-registration — the case study ships point forces and Phase 0's gate needs UQ

**Committed before any of it was computed.** Bars below are not moved; a bar that turns
out to be the wrong bar is recorded with its number and a *stricter* test is added beside
it, as items 3, 6, 11, 12 and 13 did. A pre-registration is never rewritten; corrections
go in an addendum at the foot.

## 0. Why, and what is already known

The venture's Phase-0 unit 0.6 (`../../../biofluids_venture/05_ROADMAP.md`) requires the
injectability case study to be *"closed form (Buckingham–Reiner) with full UQ from the L2
posterior"*, and its acceptance line requires the report to *state which parameters the
injection force is insensitive to*. `results.md` ships POINT forces.

Until item 13 that was not fixable: `force_with_uq` read profile upper ends that were
bisection-grid artefacts. Item 13 closed that and measured the size of the thing being
added — **the Casson push-through hull is a median 29.1 % of the force at this needle,
against a 12.5 % family-choice effect** (`../../data/profile_interval_g139.md`). The
uncertainty is the larger of the two, so shipping point values understates the honest
error by more than the effect this study's own cross-check is about.

### 0.1 Two things measured before this file was written

Both are diagnosis, not results, and both licensed the approach rather than tested it.

- `results.json` regenerates with only `git_head` moving (verified 2026-09-06).
- The forces in `data/profile_interval_g139.json`, which come from a **refit** of each
  donor curve, agree with the shipped `results.json` forces to a worst **1.2e-06**
  relative over 20 donors. That residual is explained: records store parameters at
  `%.6g` (`records_from_horner.py`), and six significant figures is a part in `1e-6`.

## 1. The approach, fixed here

The per-donor intervals are **read from `data/profile_interval_g139.json`**, not
recomputed. Computing them a second way would create a second source of truth for a
number this repo already ships, which is the failure mode recorded in `NEXT_STEPS.md`
("two pieces of code that construct the same thing by different rules will eventually
disagree"). G14.2 and G14.3 exist to check that reading them across is legitimate.

## 2. The gates

| id | claim | bar | what a failure means |
|---|---|---|---|
| **G14.1** | the study reproduces before it is extended | `results.json` regenerates with only `git_head` changed | do not extend a claim you cannot reproduce |
| **G14.2** | the borrowed intervals belong to these forces | refit vs shipped force agree to **<= 1e-5** relative on **20/20** donors | the two are not the same calculation and the intervals may not be borrowed |
| **G14.3** | both modules describe the SAME device | `run.py` and `crosscheck.py`'s needle ID, needle length, barrel radius and `Q` are **exactly equal** | **a defect to fix with one source of truth, never a tolerance to widen** |
| **G14.4** | an open interval is never printed as a number | **0** donors with a numeric interval printed while `profile_open` is non-empty | the report would quote a limit the data does not set |
| **G14.5** | donors can still be ranked at 95 % | **>= 50 %** of the 190 donor pairs have disjoint hulls | see §3 — the ordinal claim needs restating, not deleting |
| **G14.6** | the insensitivity list is reproduced by an independent route | `sensitivities()` names `tau_y` insensitive on **>= 19/20** donors | either the log-sensitivity route or the earlier `x = tau_y/tau_w` route is wrong; say which |

### Why G14.5 is the gate that matters, and why it is not the obvious one

`results.md`'s surviving ordinal claim is that the donor-to-donor spread (128 %) is
larger than the family-choice effect (30 %). Adding a 29 % uncertainty does not threaten
that arithmetic, so re-checking it would be a gate that cannot fail.

**The comparison that CAN fail is a different one, and item 14's own text in
`NEXT_STEPS.md` flags it: the hull is PER DONOR and the spread is ACROSS donors.** A
customer does not use the spread; they ask whether *this* patient needs more force than
*that* one. That is a pairwise question, and a 29 % per-donor interval against a 128 %
total range does not obviously answer it for the middle of the distribution.

So G14.5 asks it directly: over the 190 unordered donor pairs, how many have disjoint
95 % hulls?

### G14.6 is a discriminating check, not a restatement

`results.md` already says the yield stress is invisible here, on two grounds:
`x = tau_y/tau_w` reaches only 1.7e-4, and zeroing `tau_y` moves the force under 2.9 %.
`sensitivities()` is a **third and independent** route — a central-difference
`d ln F / d ln p` — and it can disagree. If it does, one of the three is wrong and the
report says which rather than quoting the two that agree.

## 3. What happens if G14.5 fails — decided now

**The ordinal claim is restated, not deleted, and not rescued by widening anything.**
If fewer than half the pairs separate, the honest form of the claim is
*"the extremes are separated at 95 %; the middle of the donor distribution is not"*,
with the separating fraction printed and the rank at which separation fails named. The
128 %-versus-30 % arithmetic stays on the page because it is true; what changes is the
sentence that reads it as "which patient matters more", which is a pairwise statement
the arithmetic does not license.

The bars `INSENSITIVE`, and the item-13 interval machinery, are not touched.

## 4. Predictions, so that being wrong is visible

1. G14.1, G14.2, G14.3 pass. (G14.3 is the one I would least like to be wrong about:
   the two modules carry the constants as *separate literals* with the same values.)
2. G14.4 passes trivially for `casson` — the item-13 run reported no open ends on the 20
   blood Casson records. **If it does not, that is the more interesting outcome.**
3. **G14.5 FAILS.** A median hull of 29 % against a total range of 128 % means adjacent
   donors overlap heavily, and 190 pairs is dominated by nearby ones. I expect roughly
   35–50 % of pairs to separate. Written down as a failure prediction so that a pass
   cannot be read afterwards as what was expected.
4. G14.6 passes at 20/20, because `x <= 1.7e-4` makes `d ln F / d ln tau_y` tiny by the
   same algebra.

Prediction 3 is the one this item exists to test. If it passes, §3 does not fire and the
ordinal claim stands as written.
