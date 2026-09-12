# R1 — the call log, and the rule that reads it

Date opened: 2026-09-12. The gate was set on 2026-09-07 and had no log, no target list and
no draft on disk until now.

## Files

| file | what it is |
|---|---|
| `R1_TARGETS.md` | fourteen named, sourced targets across the four §7 segments |
| `R1_OUTREACH.md` | four send-ready drafts, under 150 words each, with the fill-in slots marked |
| `r1_log.csv` | **the data.** One row per target, filled in by hand after each contact |
| `score_r1.py` | reads the CSV and prints the pre-registered verdict; does not change it |

## The rule, pre-registered and not to be adjusted after data exists

From `docs/market_research_2026-09-06/00_SYNTHESIS.md` §7:

> R1 = number with paid-before OR would-pay ≥ 3 (of ten calls).

`score_r1.py` implements it with two guards, both chosen before any call:

1. **A verdict needs ten HELD CALLS.** A no-reply is not a call. An email thread that never
   reached the questions is not a call. Below ten the verdict is **UNDETERMINED**, which is
   a real answer and is not the same as a pass or a fail — three yeses out of the first
   three is still undetermined.
2. **A would-pay counts only if the CONTACT named a deliverable.** "That sounds useful" is
   not a would-pay. Rows where `would_pay=yes` but `deliverable_named` is empty are printed
   separately as *soft yes* and are excluded from R1.

Guard 2 exists because the whole risk in `00_SYNTHESIS.md` §8 is that *no buyer quote
anywhere says "we wanted an uncertainty on the constitutive law"*. A gate that counts
politeness as demand cannot detect that, and a gate that cannot fail is not a gate.

**Consequence, written down now so it cannot be renegotiated later:** R1 < 3 fires D-A's
pre-registered fallback — drop W1, move to rank 4 (battery slurry / coating window). That is
a switch of wedge, not a reason to re-run the calls with a friendlier script.

## Verified behaviour of the scorer (run 2026-09-12)

Five control cases on throwaway copies of the log; the real `r1_log.csv` was not touched:

| case | expected | got |
|---|---|---|
| 10 held, 3 paid-before | PASS, exit 0 | `R1 PASS (3 >= 3)` ✓ |
| 10 held, 2 paid + 1 soft yes | FAIL, exit 1 | `R1 FAIL (2 < 3)` ✓ |
| 4 held, all four yes | UNDETERMINED, exit 2 | `an early 4/4 is not a pass` ✓ |
| `call_held = maybe` | malformed, exit 3 | `'maybe' is not yes/no/blank` ✓ |
| 10 held, 3 would-pay **with** a named deliverable | PASS, exit 0 | `R1 PASS (3 >= 3)` ✓ |

Both branches of the yes-rule and both guards have now fired at least once, so none of them
is untested code.

## Filling in the CSV

| column | how to fill it |
|---|---|
| `target_id` | the row id from `R1_TARGETS.md` (`A1`, `C2`, …) |
| `segment` | exactly one of `formulation`, `device`, `regulatory`, `platform` — anything else is not counted and the scorer says so |
| `date_contacted` | the day the note went out, not the day it was drafted |
| `channel` | `email`, `linkedin`, `intro`, `conference` |
| `replied` / `call_held` | `yes`/`no`. A written answer that covers the questions counts as a held call; say so in `notes` |
| `paid_before` | §7 Q2: has this group paid an outside party for rheology-to-force prediction? |
| `would_pay` | §7 Q7: would the one-page verdict be worth paying for? |
| `price_named_usd` | only if **they** name a number. This is the single most valuable cell in the file — `00_SYNTHESIS.md` §8 records that no public price for this service exists anywhere |
| `deliverable_named` | what *they* said they would want to receive, in their words |
| `referral_name` | §7 Q8. Referrals are the only channel F9 finds evidence for |
| `verbatim_quote` | one sentence in their words. Do not paraphrase into agreement |

## Run it

```
python docs/R1/score_r1.py
```

Current output: `UNDETERMINED -- the log is empty. The gate has not been run.`

## The other number to watch

`08_AI_EXPOSURE.md` §5.1: the ledger that matters is **calls held, proposals out, invoices
paid**. As of 2026-09-12 that is **0, 0, 0**, against 91 commits and 12,765 lines of code.
Nothing in `biorheo` moves any of those three.
