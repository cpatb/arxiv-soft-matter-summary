# PREREG — the lognormal large-N replicate that would close the finite-size flag

**Registered 2026-09-17. Nothing below has been run.** This document is written *before* any
cell is launched and before any number it predicts exists. Launching it is the author's
call (§6).

> **Filename note.** This was commissioned as `PREREG_LOGN_N3000_<date>.md`. It is not
> called that, because §1 shows **N = 3000 does not clear the flag** — it is already on
> disk and still flags on 18 of 24 draws. A file named for an N that its own contents
> refute is the kind of mislabelling this project has been bitten by before. The
> registered size is **N = 6000**.

Evidence base: `out/diag_finitesize_n.txt`, printed by `python diag_finitesize_n.py`.
**Every number in this document is quoted from that file or from `paper/numbers.json`**;
none is typed from memory. Re-run the diagnostic before reading the prediction table — if
its §0 reproduction fails, this registration is void.

---

## 0. The claim this exists to close

`paper/main.tex` §III B discloses, and does not resolve:

> `\fabFlaggedTopLognOf` = **12 of 12** lognormal fabric cells on the primary rung at
> CV = 0.5 carry `d_max/L_dense > 0.25`, against **none** in the two zero-skewness
> families; and dropping the wholly flagged rung changes the lognormal fabric slope from
> `+0.197 ± 0.225` to `-0.070 ± 0.265` — the verdict survives but the point estimate
> changes sign.

and in the Discussion:

> "The lognormal family at CV ≥ 0.3 carries the finite-size flag at N = 1500, which
> touches the fabric slope and, **through the same tail, may touch its φ_J**; a larger-N
> repeat of that family is the obvious next cell to spend."

The φ_J half of that sentence is what this registration tests. The fabric half is **not**
covered here: the fabric campaign has its own 12-seed grid and closing it needs its own
registration.

---

## 1. Premise check — done first, and it moved the design

Two premises were checked before writing any gate. Both came back differently from the way
the Status card states them.

### 1a. The N ≥ 3000 replicate has partly been run, and it is not a ladder

`out/phase2*.jsonl` already holds a lognormal block at **N = 3000, μ_p = 0.3, seeds 1–3,
all four p̂ rungs**, at CV = 0.0 and CV = 0.5. It cannot be read:

| ladder | rungs | converged | usable? |
|---|---|---|---|
| CV = 0.0, seed 1 / 2 / 3 | 4 / 4 / 4 | 4 / 4 / 4 | USABLE |
| CV = 0.5, seed 1 | 4 | 3 | **NOT A LADDER** |
| CV = 0.5, seed 2 | 4 | 2 | **NOT A LADDER** |
| CV = 0.5, seed 3 | 4 | 3 | **NOT A LADDER** |

4 of 24 cells stopped at `MAX_ITER = 200 000` without converging, and the bottom rung
(p̂ = 3e-05) is among them **on all three seeds** — the rung the zero-pressure
extrapolation leans on hardest. `(N, seed)` pairs with both ends of the paired CV = 0 → 0.5
difference usable: **0**. So no paired Δφ_J at N = 3000 exists and none can be recovered
from disk.

This is item 4 of §3 of the Status card — *`MAX_ITER` is an N-dependent criterion
project-wide* — showing up in the cells that would answer this question. The paper's own
N = 1500 lognormal CV = 0.5 cells were run at the **raised** cap of 1 000 000; the N = 3000
block was not.

### 1b. N = 3000 would not clear the flag anyway

`d_max/L_dense` is **not a simulation output**. Reading `packing.py:525-547` it is

```
d_max/L_dense = d_max / ( Σ_i (π/6) d_i³ / PHI_DENSE_BOUND )^(1/3),   PHI_DENSE_BOUND = 0.72
```

a deterministic function of (family, CV, N, seed, truncate) evaluated on the drawn
diameters before a single time step runs. `diag_finitesize_n.py` §0 regenerates it from
the project's own PSD generator and matches all 9 on-disk draws to **4.4e-13** — so the
column below is a calculation, not an extrapolation.

| N | seed 1 | seed 2 | seed 3 | worst of 3 | median of 24 | max of 24 | flagged / 24 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1500 | 0.3214 | 0.3067 | 0.3364 | 0.3364 | 0.3084 | 0.3857 | **24/24** |
| 3000 | 0.2511 | 0.2975 | 0.2651 | 0.2975 | 0.2715 | 0.3020 | **18/24** |
| 4000 | 0.2468 | 0.2687 | 0.2423 | 0.2687 | 0.2565 | 0.2778 | 15/24 |
| 5000 | 0.2299 | 0.2489 | 0.2522 | 0.2522 | 0.2412 | 0.2582 | 6/24 |
| **6000** | **0.2164** | **0.2345** | **0.2385** | **0.2385** | 0.2319 | 0.2471 | **0/24** |
| 8000 | 0.2013 | 0.2118 | 0.2164 | 0.2164 | 0.2117 | 0.2257 | 0/24 |

The flag falls *slower* than N^(−1/3) because the largest grain is itself a sample maximum
and grows with N until it saturates at the truncation ceiling 5.0 d₁₀: measured medians go
0.3084 → 0.2715 over the 1500 → 3000 doubling, ratio 0.8804 against N^(−1/3) = 0.7937, i.e.
the drawn maximum grew by 1.1092.

**N = 6000 is registered** as the smallest size on this grid that clears the flag on all
three registered seeds *and* on all 24 sampled draws. Margins, both printed above:
4.6 % below threshold on the worst **registered** seed (0.2385), **1.2 %** on the worst of
24 (0.2471). The 1.2 % figure is the honest one to quote about the family; the 4.6 % is
what applies to the cells actually registered here, because the seeds are fixed in §2 in
advance and their values are already on the table. **No seed may be substituted after a
cell fails.** If a seed has to be dropped it is recorded as an amendment with the reason,
and the family runs at the reduced n with its own t.

---

## 2. The cells — fixed here, in advance

| | |
|---|---|
| runner | **`phase2_phiJ.py`, unchanged** — only its arguments move. It is resumable (a cell already in the shard file is skipped), so an interrupted campaign restarts without re-running anything |
| invocation | `--family lognormal --n 6000 --mus 0.3 --cvs 0.0,0.5 --seeds 1,2,3 --tag <new> --shard i --nshards K`, then the same with `--family uniform` |
| families | **lognormal** (the flagged one) and **uniform** (the control, §3 G5) |
| CV | 0.0 and 0.5 — both ends, because the claim is a *paired* difference |
| N | **6000** |
| seeds | **1, 2, 3** — the manuscript's own seeds |
| μ_p | 0.3 |
| p̂ rungs | 1e-03, 3e-04, 1e-04, 3e-05 (4 rungs, the campaign ladder) |
| `MAX_ITER` | **1 000 000**, set through the env var `PJ_MAX_ITER` (`packing.max_iter()`), matching the paper's own N = 1500 lognormal CV = 0.5 cells. `phase2_phiJ.py:130` appends `_it{cap}` to the record key at any non-default cap, so these records **cannot pool** with a 200 000-cap ladder — the separation is structural, not a convention |
| truncate | (0.2, 5.0) |
| tag | a **new** tag matching none of `phase2*`, `shp*`, `mixr*`, `extr*`, `tern*`, `kl*`, `ug*` — so these records cannot pool into any existing ladder |

**Cell count: 2 families × 2 CV × 3 seeds × 4 rungs = 48.**

The uniform arm is not optional. Without it a shift in Δφ_J at N = 6000 cannot be
attributed to the flag rather than to N, and the uniform family is the clean control
because it carries **no** flag at any size on this grid (`d_max/L` = 0.1496 / 0.1485 /
0.1502 at N = 1500 and 0.1200 / 0.1188 / 0.1199 at N = 3000).

---

## 3. Gates — bars quoted verbatim, predictions blind

The registered comparison values, all from `paper/numbers.json` as committed at 5573617:

```
dPhiJLogn   = +0.0621     dPhiJLognSe   = 0.0026
dPhiJUnif   = +0.0242     dPhiJUnifSe   = 0.0008
dPhiJBidi   = +0.0139     dPhiJBidiSe   = 0.0025
phiScale    = 2.1e-03     (run-to-run scatter on repeated monodisperse cells)
ratioPhiJ   = 4.47        (logn / bidi; quoted in the abstract as "three to five")
tCritTwoDof = 4.303
```

### G1 — WIRING, not a gate

The runner's own `d_max/L_dense` at N = 6000 must equal **0.2164 / 0.2345 / 0.2385** for
seeds 1 / 2 / 3 at CV = 0.5, to the precision the shards store.

⚠ **This cannot fail on the physics** — §1b already computed it from the same generator.
It is recorded as a **wiring check**: what it tests is that the runner's insertion RNG
order and truncation at N = 6000 are the ones `diag_finitesize_n.py` modelled. A
disagreement means the draw changed, and everything below is void.

### G2 — LADDER HEALTH (can fail; a failure at N = 3000 is on record)

All 48 cells converge under `MAX_ITER = 1 000 000`, and each of the 12 ladders has
r² ≥ 0.85 with a non-negative slope (`r2Min = 0.85`, the campaign's own bar).

**Prediction: 48/48 converge.** Basis: the paper's N = 1500 lognormal CV = 0.5 cells all
converge at this cap, and the 4 failures at N = 3000 were all at the 5× lower cap.

If fewer than 48 converge, the non-converged cells are **listed**, and φ_J is not quoted
for any ladder missing a rung. A ladder is not repaired by dropping its bottom rung.

### G3 — THE PRIMARY: does Δφ_J move? (can fail)

Let `D6 = Δφ_J^logn(N=6000)` be the seed-paired CV = 0 → 0.5 change, computed by exactly
the estimator `sweep_descriptive.py` uses at N = 1500.

**G3a — statistical.** `D6` and `+0.0621` agree under a Welch t at **the design's own
degrees of freedom**, with the critical value taken from that dof — not 2, and not 4.303
unless the dof works out at 2.

**G3b — material.** `|D6 − 0.0621| ≤ 3 × phiScale = 0.0063`.

**Predictions (blind, made now): G3a PASS, G3b PASS.**

What the G3b bar is worth, stated so a reader does not have to guess:

* it is **10.1 %** of the effect (0.0063 / 0.0621), i.e. far stricter than this project's
  existing finite-size convention of "safe if the shift is under a third of the effect"
  (which would be 0.0207);
* on the upside it is the bar that **keeps the abstract's "three to five"**: with the
  bidisperse arm held at 0.0139, a ratio of 5 needs Δφ_J^logn ≤ 0.0695, i.e. a rise of
  0.0074. So 0.0063 is slightly stricter than the claim needs on that side, and much
  stricter on the other;
* it is **1.25 ×** the worst stopping-rule irreproducibility this protocol has measured
  (5.03e-03 in φ at CV = 0, `out/capfix_grade.txt` §(ii)/(iii)). **A failure only just
  above the bar therefore cannot be distinguished from stopping-rule noise on one seed
  triple**, and must be read against the uniform control (G5) before it means anything.

**Power guard — so the bar is not unachievable, and the guard itself is not vacuous.**
The standard error of the difference is `se_diff = sqrt(0.0026² + se6000²)`, with `se6000`
unknown until the run. **Guard: if `se_diff > 0.0063`** — the design cannot resolve its own
bar even at one standard error — the verdict is **UNDETERMINED**, not FAIL. `se_diff` is
printed beside the verdict either way.

⚠ The first version of this guard read `se_diff > 0.0063 / 2` and is **retracted here**,
before the run, because it could only ever have returned UNDETERMINED: substituting the
best available prior for `se6000` (the measured `se1500 = 0.0026`) gives
`se_diff ≈ 0.0037`, already above `0.0063 / 2 = 0.00315`. A guard that fires on every
admissible input is not a guard. The number is left in the record rather than deleted.

**Expected behaviour of G3b under a true zero shift**, registered so a marginal result is
not over-read: at `se_diff ≈ 0.0037` the bar sits at **1.7 se_diff**, so a genuinely
unmoved Δφ_J clears it roughly 5 times in 6 (≈ 0.84 on a t with ~4 dof). **A single FAIL
landing just above 0.0063, with the uniform control steady, is therefore a marginal
result and not a refutation** — it is reported as such, with `se_diff` and the z value,
and it does not by itself re-quote the headline.

### G4 — the flag rate (can fail only through G1)

Flagged cells among the 24 lognormal N = 6000 records: **0 of 24 predicted**. Reported as
a count, not a rate.

### G5 — THE CONTROL: is any shift specific to the flagged family? (can fail)

`D6_unif = Δφ_J^unif(N=6000)` against `+0.0242 ± 0.0008`, same two arms as G3.

**Prediction: G3 and G5 both pass, i.e. neither family moves.** The discriminating
outcome, and the reason the control is here: **if the lognormal arm moves and the uniform
arm does not**, the shift is attributable to the flag; **if both move**, it is a generic
N effect and says nothing about the tail. A lognormal-only run could not tell these apart.

⚠ The candidate-separation caveat applies: if neither family moves, the control has
**discriminated nothing** and is reported as UNINFORMATIVE rather than as corroboration.

### G6 — what already-known N dependence predicts (context, not a gate)

`out/finite_size_verdict.txt` grades the φ_J channel **SAFE** on 30 matched N = 1500 vs
3000 cells (max N shift 0.0091 against an effect of 0.0623, ratio 0.15) and finds the CV
trend holds (+0.0609 → +0.0600, ratio 0.98). **That file does not settle this question**,
for two reasons that must be restated with any use of it: it compares φ at *fixed p̂*,
not the extrapolated φ_J; and it is computed on the same partly non-converged N = 3000
block that §1a disqualifies.

---

## 4. What each outcome does to the manuscript

Registered in advance so the write-up cannot be chosen after the fact.

| G3a | G3b | G5 | consequence for `paper/main.tex` |
|---|---|---|---|
| pass | pass | pass | The Discussion sentence "through the same tail, may touch its φ_J" is **replaced** by the measurement, with `D6` and its s.e. quoted and the N = 6000 flag count printed. The §III B fabric caveat is untouched — this run does not address it. |
| pass | fail | either | **Not resolved.** The caveat stays, and the measured shift is printed beside it. Both numbers ship. |
| fail | — | uniform steady | Δφ_J^logn is N-dependent. The headline is re-quoted at N = 6000 or withdrawn; the ratio and the abstract's "three to five" are recomputed. |
| fail | — | uniform moves too | A generic N effect, not a tail effect. The whole φ_J ladder needs an N study before any family number is quoted. |
| UNDETERMINED | — | — | The design is under-powered. Report the bound, keep the caveat, and say what seed count would resolve it. |

**No bar in §3 may be loosened after a failure.** If a bar turns out to be the wrong test
rather than a failed one, the replacement must be stricter *in kind*, and the original
number stays printed — the standing rule this project has applied five times.

---

## 5. Cost, and why it cannot run from an agent session

Measured, from `out/diag_finitesize_n.txt` §3:

| cv | N | n | iters med | iters max | wall med (s) | wall max (s) | non-converged |
|---|---|---|---|---|---|---|---|
| 0.0 | 3000 | 12 | 21 946 | 70 176 | 151 | 494 | 0 |
| 0.5 | 1500 | 12 | 118 598 | 367 411 | 280 | 1 492 | 0 |
| 0.5 | 3000 | 12 | 156 771 | 200 000 | 500 | 657 | **4** |

Wall time at CV = 0.5 scales as N^0.84 across the two measured sizes, giving **≈ 893 s**
for a median N = 6000 cell. ⚠ **That is a floor, not an estimate**: the N = 3000 row is
**censored at the 200 000 cap** on 4 of 12 cells, so both its median iteration count and
its median wall time are under-reported, and the true exponent is larger. At the raised
cap a CV = 0.5 cell that needs the full 10⁶ iterations would run **≈ 1.5–2 h**.

Budget: **48 cells, ≈ 12–30 h serial**, ≈ 1.5–3 h wall at 12 workers, dominated by the 12
lognormal CV = 0.5 cells.

### 6. Launcher requirement — **AUTHOR'S CALL, and it blocks the run**

**Nothing spawned by a `wsl.exe` call survives that call exiting** — not `nohup`, not
`setsid nohup`, and not while a second WSL session is held open. With a 600 s cap on one
agent call, the longest cell an agent session can complete is ~600 s, and §5 says the
binding cells are 900–7 200 s. **This campaign cannot be launched from an agent session at
all.**

Running it needs both of:

1. `vmIdleTimeout = -1` in `.wslconfig`, and
2. a **service-style launch** (a Windows scheduled task or service that owns the WSL
   process), not an inline `wsl bash -lc` from a session.

Both are changes to the author's machine configuration. **This is the author's decision and
this registration does not presume it.** Launch discipline that has worked before, once
that decision is made: write the script with an editor, strip CR, copy it into `~/pj`,
compare md5 **on both sides**, and invoke it **by file path from PowerShell** — never
inline, and never through a double-quoted `wsl bash -lc "cat > f <<'EOF' ..."`, which lets
Git Bash expand `$` first and silently empties every shell variable.

---

## 7. What this registration does *not* cover

* The **fabric** flag. `T_nr`'s lognormal slope sign-flips when the flagged rung is
  dropped; that is a 12-seed campaign on its own grid and needs its own prereg.
* The **exponent** of the approach to jamming. The ladder is fitted linearly by choice and
  the paper says so; nothing here tests the free exponent.
* Whether **N = 1500 results should be withdrawn**. Even a G3 failure re-quotes the
  number at N = 6000; it does not retract the N = 1500 measurement, which is correctly
  reported as a measurement at its own N.

---

*Registered 2026-09-17, before any cell. Evidence: `out/diag_finitesize_n.txt`
(`python diag_finitesize_n.py`). Predictions in §3 are blind and must not be edited after
the run; a changed prediction is an amendment with a date and a reason.*
