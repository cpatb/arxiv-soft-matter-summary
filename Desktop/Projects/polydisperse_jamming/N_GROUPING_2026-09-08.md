# System size was being pooled into the ladder (2026-09-08)

Every number below is printed by `diag_n_pooling.py` (`out/n_pooling.txt`,
`out/n_pooling.json`) or by `analyse_phase2.py` (`out/phase2_verdict.txt`). Pinned by
`test_n_grouping.py` (12 tests).

## What was wrong

`analyse_phase2.py` de-duplicates on a record key that encodes N — its own docstring says
*"the cell key encodes cv, mu, seed, p_hat and N"* — but the three keys it then **grouped**
by did not:

| purpose | key, before | contained N? |
|---|---|---|
| ladder fit → `phi_J`, `b_free` | `(family, cv, mu_p, seed)` | no |
| seed scatter → the gate's denominator | `(family, cv, mu_p, p_hat_target)` | no |
| seed aggregate → the published `phi_J` table | `(family, cv_nom, mu_p)` | no |

`out/` holds cells at **N = 1500 and N = 3000**, so a ladder could carry the same `p_hat`
twice at two system sizes and a "seed scatter" cell could be measuring a system-size
difference. **18 of 95 ladders and 20 of 114 scatter cells were mixing.**

The seed scatter is the denominator of the gate that decides which ladders are informative,
and the ladder fit produces `b_free`, the project's pre-registered exponent anchor. So the
defect sat under both halves of that anchor.

## What it was worth — measured, not argued

- median seed scatter in `phi`: **0.001811 pooled → 0.001776 split**, a 1.9 % move in the
  denominator.
- per-ladder `phi_J`: **median |Δ| = 0.00104, max 0.00521 = 2.88 seed scatters.** Above
  one seed scatter is this project's own threshold for calling something N-sensitive, so
  the correction is a correction and not housekeeping.
- **The pooling did NOT rescue the exponent anchor.** Median `b_free` is **0.518 pooled**
  and **0.452 split**, against a band of `[0.6, 1.6]` pre-registered before today: **FAIL
  either way.** That is the useful outcome — a defect that *could* have explained the
  anchor's failure, ruled out by measurement.

**The honest cost of the fix, stated because it cuts the other way.** Splitting N takes
ladders from 8 points to 4, or from 6 to 3, and `fit_phi_J` returns `b` as NaN below 4
points. So separating N does not merely change `b` — it **disqualifies ladders**
(31 → 48 with no free fit). Pooling had been quietly buying the free fit its degrees of
freedom.

## The fix, and the gate it immediately re-decided

N is now in all three keys. **Adding it broke anchor A4 on the first run** — and not
because any physics changed. A4 asks whether `phi_J` falls monotonically in `mu_p`; the
summary table now carried two rows at `mu_p` = 0.3, one per N, so the sequence read
`['0.0:0.6332', '0.3:0.5745', '0.3:0.5764', '1.0:0.5590', '3.0:0.5574']` → **FAIL**, on a
different question than the one A4 asks. This is the same failure mode as the defect being
fixed, one level down. A4 is now graded **per N**: **PASS at N = 1500**, "not gradeable" at
N = 3000 (one `mu_p` value there).

A5 had the same exposure and worse: it sorts by `cv_real` and takes the **first and last**
row, so without N in the loop key it could weld the low-CV row of one system size to the
high-CV row of another. It is now keyed on N, and every block names its N.

**A guard was added rather than a comparison assumed.** Two blocks at different N are only
a finite-size comparison if they span the same CV range, and they need not: the uniform
family has N = 3000 cells only at `cv` 0.4–0.5 against 0.0–0.5 at N = 1500. The verdict
file now prints

```
!! CV SPAN DIFFERS BETWEEN SYSTEM SIZES: N=1500: cv 0.00-0.50, N=3000: cv 0.40-0.50
!! These blocks are NOT a matched finite-size comparison.
```

so the uniform N = 3000 block's 0.99× cannot be read as contradicting the 0.73× at
N = 1500. They measure different intervals.

## Four siblings had the same blindness, and their numbers had already moved

Fixing `analyse_phase2.py` is not the same as fixing the defect. Grepping for the property
rather than for the shape of the bug (`n_requested` references per file) found **four
consumers with zero**:

| file | what it carries | had `n_requested`? |
|---|---|---|
| `analyse_skewness.py` | claim **C4**, lognormal vs zero-skew at fixed CV | no |
| `analyse_families.py` | claim **C3**, families agree at fixed CV | no |
| `diag_convergence.py` | is the non-converged set a biased sample? | no |
| `diag_ladder.py` | the raw `phi(p_hat)` ladders | no |
| `analyse_exponent_constrained.py` | the constrained exponent | yes |
| `analyse_finite_size.py` | the N ladder itself | yes |

Their printed numbers **moved the moment the new cells were ingested** — `analyse_skewness`
alone shifted every row of its comparison table — without anyone choosing that.

C3 and C4 are between-**family** comparisons at fixed CV, and the N = 3000 coverage does not
span the CV range they need (bidisperse: one CV; uniform: two). So the honest form is not a
per-N split but a **stated restriction to the size every claim is quoted at**, through a
shared `restrict_to_primary_N(recs, why)` in `analyse_phase2.py` that **prints the count it
drops and the sizes present** on every run. A later N sweep therefore cannot quietly stop
being included.

**The restriction restored the published numbers exactly.** `families_verdict.txt` and
`skewness_verdict.txt` are byte-identical to their committed state again, and
`skewness_comparison.csv` differs only in the 15th significant figure (a summation-order
change), with all six agree/differ verdicts unchanged. That is the check that says the
restriction is the right one rather than merely a choice.

## What came out of it for free: a finite-size check on the headline verdict

The lognormal family at `mu_p` = 0.3 is the one pair that **does** span the same CV range
at both sizes, so the M2 verdict is now graded independently at each:

| | `phi_J`(cv≈0) | `phi_J`(cv≈0.5) | Δ | Δ/scatter | effect at φ = 0.44 | verdict |
|---|---|---|---|---|---|---|
| N = 1500 | 0.5745 | 0.6369 | +0.0623 | 13.11 | **0.467×** | INSUFFICIENT |
| N = 3000 | 0.5764 | 0.6367 | +0.0603 | 11.14 | **0.481×** | INSUFFICIENT |

**The two sizes agree to 3.0 % on the effect and 3.2 % on Δ`phi_J`, and reach the same
verdict.** That is a finite-size corroboration of the project's headline M2 result, and it
came out of fixing the grouping rather than from a new run.

## The data that prompted this

Three `wsl_n3k_*.jsonl` files had been sitting **untracked in the project root** since
2026-09-04, where no loader could see them (every consumer globs `out/phase2*.jsonl`).
They are now `out/phase2_n3k_bidisperse.jsonl` and `out/phase2_n3k_uniform.jsonl`.

⚠ **`wsl_n3k_lognormal.jsonl` was a byte-for-byte duplicate** of cells already in
`out/phase2_n3000_shard*.jsonl` — all 12 keys already present. It is kept as
`out/phase2_n3k_lognormal_DUPLICATE.jsonl.txt`, renamed out of the glob so it cannot be
loaded twice, rather than deleted, because the fact that one sync landed twice under two
names is worth knowing.

**Ingesting them before the grouping fix would have made things worse, not better**: the
largest per-ladder `phi_J` move (0.00521, 2.88 seed scatters) is a *uniform* ladder, i.e.
it comes from exactly the cells that were about to be added.

## What this does not say

- Not that N = 3000 is wrong, or that `phi_J` is wrong. The N ladder itself was measured
  separately (`analyse_finite_size.py`) and the extrapolated `phi_J` moved 0.0033 when N
  doubled.
- Not that the exponent anchor is rescued. It fails at `b` = 0.452 against `[0.6, 1.6]`,
  and the bar was not touched.
- Not a finite-size check on bidisperse or uniform. Bidisperse has N = 3000 at one CV only,
  so it produces no A5 block at all; uniform's CV span does not match.
