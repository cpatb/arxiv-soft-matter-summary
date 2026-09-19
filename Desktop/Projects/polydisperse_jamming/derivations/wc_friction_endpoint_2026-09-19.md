# μ_p = 0.3 is not the frictional endpoint, and the Wyart–Cates window's flatness was tested with the wrong yardstick

**2026-09-19.** Every number is printed by `wc_friction_endpoint_2026-09-19.py` into
`wc_friction_endpoint_2026-09-19.txt` (rc 0, ALL CHECKS PASSED). No new cells: this runs on
`out/phase2_fits.csv`, which has held the `μ_p = 1.0` and `μ_p = 3.0` blocks since before
this session. **No analysis in this repository had ever read them.**

> ## RESULT
>
> **Two independent errors in `wc_endpoints_move_together_2026-09-19.md` were cancelling.**
> Correcting either one alone changes its conclusion; correcting both brings the conclusion
> back, for a better reason than it was first given.
>
> ### 1. 🔴 `φ_J(μ_p = 0.3)` is not `φ_m`
>
> | `μ_p` | `φ_J` (lognormal, CV = 0) |
> |---|---|
> | 0.0 | 0.63237 ± 0.00042 |
> | 0.3 | 0.57453 ± 0.00038 |
> | 1.0 | 0.55904 ± 0.00067 |
> | 3.0 | 0.55742 ± 0.00097 |
>
> `μ_p = 0.3` captures **77.2 %** of the `μ_p = 0 → 3.0` drop. The remaining **22.8 %**
> (`+0.01711 ± 0.00120`) is resolved at **\|t\| = 14.29**. `μ_p = 1.0` is not resolved from
> 3.0 (`+0.00161 ± 0.00058`, \|t\| = 2.78), so it is the best `φ_m` on disk.
> ⚠ The campaign chose 0.3 for timing, and said so: `phase2_phiJ.py` — *"`mu_p = 1.0`
> dropped from the production grid for the same timing reason."*
>
> ### 2. 🔴 The window is **not** flat at `μ_p = 0.3`, on the right statistic
>
> The published page tested flatness as a **range** over six rungs against twice a
> **two-point** propagated s.e. That yardstick does not cancel the shared seed draw. The
> per-rung **paired** statistic does, and its s.e. is up to 25× smaller.
>
> | statistic | `φ_m` at 0.3 | `φ_m` at 1.0 (converged) |
> |---|---|---|
> | range over rungs | 0.00646 | 0.00240 |
> | worst per-rung \|t\| | **27.70 — RESOLVED** | 2.78 |
> | **trend `dW/dCV`, fitted per seed** | **−0.01217 ± 0.00155, \|t\| = 7.87 — RESOLVED** | +0.00398 ± 0.00178, \|t\| = 2.24 |
> | as a fraction of `W(0)` over the ladder | **10.5 %** | 2.7 % |
>
> ### 3. 🔑 And the mechanism is identifiable: **the residual friction gap grows with CV**
>
> `W(0.3-based) ≡ W(1.0-based) − [φ_J(0.3) − φ_J(1.0)]`, so the drift must live in that gap.
> Measured, paired on seed:
>
> | CV | `φ_J(0.3) − φ_J(1.0)` |
> |---|---|
> | 0.0 | +0.01550 ± 0.00073 |
> | 0.2 | +0.01542 ± 0.00068 |
> | 0.5 | +0.02136 ± 0.00146 |
>
> Growth over CV 0 → 0.5, paired both ways: **+0.00586 ± 0.00077, \|t\| = 7.65.**
>
> **Reading: a wider size distribution needs more friction to reach its frictional limit.**
> That is a statement about the *approach* to the limit, not about the limit itself, and it
> is exactly why a campaign that stopped at `μ_p = 0.3` sees a CV-dependent window where the
> converged one does not.
>
> ### 4. What survives, and it is the load-bearing claim
>
> **The CV shift is friction-independent.** Paired CV 0 → 0.5, lognormal:
> **+0.05803 ± 0.00268** (`μ_p = 0`), **+0.06207 ± 0.00255** (0.3), **+0.05621 ± 0.00332**
> (1.0). Worst pairwise \|t\| = **1.40** against 4.303, while `φ_J` itself moves 0.07333 over
> the same friction range — **13× larger than the spread in the shift.**
>
> **Verdict: SUPPORTED as a bound.** The resolving power on the worst pair is 0.01835, i.e.
> **30 % of the shift**, so "friction-independent" means "to within 30 %", not "identically".

---

## 1. Why nobody had looked

The friction ladder was cut for cost. `phase2_phiJ.py`'s own comment:

> *"`mu_p = 1.0` dropped from the production grid for the same timing reason. The A4 anchor
> (`phi_J` falls with friction) is still testable on 0.0 vs 0.3, and a supplementary
> `mu = 1.0` run can be added with `--tag` if time allows."*

The supplementary run *was* added — 12 cells at `μ_p = 1.0` (CV 0, 0.2, 0.5) and 3 at
`μ_p = 3.0` (CV 0). They then sat unread, because every downstream analysis filters to the
production grid `MU_LADDER = (0.0, 0.3)`. The A4 anchor they were run for only needs the
*sign*, which 0.0 vs 0.3 already gives.

So the situation was: the data that shows `μ_p = 0.3` is not the frictional limit was
generated, committed, and never loaded — while a page built on the assumption that it *is*
the limit was written this morning.

## 2. Additive or multiplicative separability?

If polydispersity and friction act on different things, `φ_J(CV, μ)` should factor. Two
one-parameter forms:

```
   additive        φ_J(CV, μ) = F(CV) + G(μ)   ⇒  the ABSOLUTE shift is μ-free
   multiplicative  φ_J(CV, μ) = F(CV) · G(μ)   ⇒  the RELATIVE shift is μ-free
```

| `μ_p` | absolute shift | relative shift |
|---|---|---|
| 0.0 | +0.05803 | 0.09177 |
| 0.3 | +0.06207 | 0.10804 |
| 1.0 | +0.05621 | 0.10055 |

Fractional spread: **10.0 % absolute, 16.3 % relative.** The additive form fits better.

⚠ **This is a ranking of two forms on three points, not a test that either is right.**
Neither is rejected. What it does say is that if a separable description is sought, the
additive one is the better starting point — and an additive separation is the natural
statement of "polydispersity changes the geometry of space filling while friction changes
the contact network", because the two contributions add in `φ` rather than scaling it.

## 3. What this cannot decide

1. **`μ_p = 1.0` has only three CV rungs** (0, 0.2, 0.5), so the converged window's flatness
   rests on two paired comparisons, not five.
2. **Only the lognormal family has any friction above 0.3** — checked, 12 cells for
   lognormal and **0** for the other two. The uniform and bidisperse families cannot be
   tested for this at all.
3. **`μ_p = 3.0` exists at CV = 0 only**, so "`μ_p = 1.0` is converged" is established at
   zero width and *assumed* at finite width. Given §3 of the result box — the residual gap
   grows with CV — that assumption is the one most likely to be wrong, and it points the
   same way: at CV = 0.5 the true `φ_m` may be further below `φ_J(1.0)` than at CV = 0.
4. **"Not resolved" at n = 3 is a bound.** `μ_p = 1.0` is converged to within 0.00250, which
   is 3.3 % of the total friction drop.
5. `paper/main.tex`: the frictionless ladders are not monotone in `p̂` at this seed count.
   That caveat attaches to every `φ_0` here.

## 4. What would settle it

- **A `μ_p` ladder at fixed CV = 0.5** — the cheapest thing that would turn §3's mechanism
  from identified into measured. Three or four `μ_p` values at one CV, one family.
- **`μ_p = 3.0` at CV = 0.5**, three cells, which converts limit 3 above from an assumption
  into a measurement.
- A `μ_p = 1.0` block for the uniform family, which would say whether the
  friction-independence of the CV shift is a property of the lognormal family or of the
  protocol.

## 5. Referee

⚠ **Unrefereed.** Written in the same session as the finding. Its two retractions of the
companion page are arithmetic on committed data and were verified before being written down;
the *interpretation* in §2 and §3 of the result box has not been attacked.
