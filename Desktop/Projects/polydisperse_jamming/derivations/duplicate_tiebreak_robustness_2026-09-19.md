# The duplicate tie-break does not move the published numbers — measured, not argued

**2026-09-19.** Every number printed by `duplicate_tiebreak_robustness_2026-09-19.py` into
`duplicate_tiebreak_robustness_2026-09-19.txt` (rc 0, ALL CHECKS PASSED). Read-only; the live
`out/phase2_wcbidi*.jsonl` campaign is excluded.

> ## RESULT
>
> `paper/main.tex` discloses in its Methods that some shard keys carry records that
> **disagree on `φ`** — the same physical cell re-run, returning a different answer — and
> resolves them with a tie-break ("highest `max_iter`, then first in file order"). It quotes
> the worst *input* disagreement, `\conflictWorstDphi` = 0.0040. **It never bounds the
> *output*.** This enumerates all **32** duplicate choices among the 5 affected selected
> cells and recomputes the three published shifts under each:
>
> | family | min | max | range | committed | range as % of the shift |
> |---|---|---|---|---|---|
> | lognormal | 0.06128 | 0.06252 | **0.00123** | 0.0621 | **2.0 %** |
> | uniform | 0.02420 | 0.02420 | **0.00000** | 0.0242 | **0.0 %** |
> | bidisperse | 0.01389 | 0.01389 | **0.00000** | 0.0139 | **0.0 %** |
>
> and the headline ratio spans **4.41× to 4.50×** against the committed **4.47×**, on the
> correlated pairing (same choice for both families, which is the correct one because they
> share their CV = 0 cells).
>
> **Verdict: SUPPORTED. The tie-break is not a threat to any published number.** Two of the
> three families do not move at all, because every disagreeing duplicate in the selected set
> sits under a **lognormal** key at CV = 0.
>
> This is the robustness check a referee asks for and the manuscript did not have. It cost
> nothing and it can now be quoted.

---

## 1. Why the disclosure was not enough

The Methods bound the *input*: the worst two records for one cell differ by 0.0040. That
says nothing about the *output*, because the published quantity is a **paired difference of
extrapolated `φ_J`** over three seeds — pairing cancels shared error and the four-rung fit
averages. The propagation could have been larger than the input (an extrapolation can
amplify) or smaller. **Measured, it is smaller: at most 0.00123.**

Which is also why the protocol's own per-rung re-run scatter — median **0.00234** over the
14 disagreeing keys across the whole tree — cannot be read directly as an error bar on the
shift. It is a different quantity. The enumeration is the test; the scatter is context.

## 2. 🔴 A correction to this page's own first run

The first run of this module reported that the disagreeing duplicates are **not** all at
CV = 0, contradicting `main.tex`'s Methods and `extract_numbers.py`'s own assertion, and
found a spread of **0.125** in a uniform CV = 0.5 cell (`φ` = 0.603 against 0.478).

**My filter was wrong, not their claim.** I had omitted the `converged` cut that
`sweep_descriptive.load_shards` applies. Those `φ ≈ 0.48` records are non-converged cells the
paper correctly excludes. With the cut restored the selection reproduces the paper's exactly
— **717 records → 607 at N = 1500 → 270 at `μ_p` = 0.3 → 257 converged → 217 cells** — and
every disagreeing duplicate does sit at CV = 0, 5 of them, all under a lognormal key.

The standing rule this violated: *reproduce the selection before contradicting what it
supports.* The failing first-run numbers are printed in the script's comment rather than
deleted.

⚠ **A second self-correction in the same module.** The check "the committed value lies
inside the enumerated range" initially FAILED for the uniform family — `\dPhiJUnif` = 0.0242
is a macro printed to four decimals against a computed 0.0242013, so the comparison failed on
the fifth decimal and meant nothing. It now compares at the precision the manuscript prints,
with the tolerance stated (5e-5, half the last digit).

## 3. What this cannot decide

1. **It bounds the tie-break, not the protocol.** The 14 disagreeing keys are the cells that
   *happen* to have been run twice. Cells run once could disagree with a re-run by as much
   and there is no way to know from this data.
2. **It says nothing about why the re-runs disagree.** The affected records differ in
   `iters` (e.g. 27831 against 59566 for the same key), so the likeliest cause is the
   iteration count at which the relaxation was stopped — but that is a hypothesis this page
   does not test.
3. **Only `μ_p = 0.3` is enumerated**, because that is the selection the published shifts
   use. The frictionless block has its own duplicates and they are not covered here.
4. **The three families' CV = 0 cells are shared**, so the "correlated pairing" in the ratio
   is the right one and the worst-case pairing quoted beside it is an over-estimate by
   construction. Both are printed.

## 4. Referee

⚠ **Unrefereed.** The enumeration is exhaustive and mechanical; the two self-corrections in
§2 were caught by checks written to be able to fail, and both failing states are recorded.
