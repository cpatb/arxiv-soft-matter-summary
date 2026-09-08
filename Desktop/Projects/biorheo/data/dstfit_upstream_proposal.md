# Item 17 — a PROPOSAL: should the profile-interval defect be fixed in `dst_fit`?

Written by `scripts/dstfit_upstream_proposal.py`. **This document changes nothing.** `dst_fit` is the founder's own package and this repo's standing rule is that it is consumed unchanged. The decision is D5's and this is the case, not a patch.

## The defect, in one line

`dst_fit.profiling.profile_likelihood` sizes its outward step floor from the parameter's declared **box bound** rather than from its scale:

```python
step = max(wald / 4.0 if finite else span / 200.0,
           span / 5000.0)          # span = |bound - best|
```

so the returned endpoint can be a point on the grid `best + m*span/5.12e6` rather than a likelihood crossing. Measured in `biorheo` (item 13): **201 of 296 closed upper ends** on that grid, worst **1639x the best-fit value** above it, and the shipped endpoints sat at a median **710x `f_crit`** over a 70-end pre-registered sample. A confidence limit has `F = f_crit` there, by definition.

## What is NOT in dispute

- The endpoints are not confidence limits (item 13, G13.2, n = 70).
- `biorheo`'s corrected walk matches a **closed form that touches none of its code** to 2.4e-9 on 20 records (G13.4).
- The blast radius is computable per parameter without running anything, and is now a field in every `biorheo` sidecar (item 16, G16.4).
- `dst_fit`'s own suite is green and its one figure is byte-identical today.

## The number nobody had computed: what it does to `dst_fit`'s OWN results

The defect acts through

```
exposure = (span / 5000) / (wald / 4)          span = |bound - best|
```

so **the box floor wins exactly when `wald < span / 1250`**. That depends only on a parameter's declared box and the precision the fit reports.

### Step 1 — which shipped code paths reach the profile walk at all

- `scripts/` calling `profile_likelihood` directly: **2** line(s).
- `scripts/` passing `profile=` to `assess()`: **0**.
- `tests/` passing `profile=`: **2**.

      scripts\gate_n7.py:106: pr = profile_likelihood(WYART_CATES, s, eta, "phi_m", fit=fit, level=LEVEL)
      scripts\gate_n7.py:147: pr = profile_likelihood(model, s.stress, s.viscosity, "phi_m",

**No shipped `dst_fit` script passes `profile=` to `assess()`.** The only shipped consumer of the walk is `scripts/gate_n7.py`, and it profiles **`phi_m` only**. The one place `sigma_star` is profiled is a test, and a test asserting behaviour is not a published number.

### Step 2 — `phi_m`'s exposure, which is the whole question

`phi_m`'s declared box is **[0.3, 0.639]** — a packing fraction, so it is O(1) wide. `gate_n7.py`'s own docstring quotes **phi_m = 0.5876 +/- 0.0013** on the real series. With that precision:

| end | `span` | `span/5000` | `wald/4` | exposure | box floor wins? |
|---|---|---|---|---|---|
| lower | 0.2876 | 5.752e-05 | 3.250e-04 | **0.177** | no |
| upper | 0.0514 | 1.028e-05 | 3.250e-04 | **0.0316** | no |

**Both ends are below 1**, by factors of 6 and 32. The box floor never sets the first step for `phi_m` at the precision `dst_fit` reports, so **the defect is inert for the only number `dst_fit` currently publishes through the profile walk.**

That is the same mechanism items 15 and 16 found downstream: the defect's size scales with the DECLARED BOX, and a parameter whose box is O(1) is barely touched while a viscosity bounded at 1e6 is destroyed.

### Step 3 — what it WOULD do to the wide-boxed parameter

`sigma_star`'s box is **[0.0001, 100000]**. At an order-of-magnitude fitted value of 44.2 Pa with a standard error of 4 Pa — **illustrative numbers, not a fit** — the upper end's exposure is

    exposure_hi ~ 20

i.e. the first outward step would be about **20x** the Wald-scaled one. In absolute terms that step is `span/5000` = **20 Pa** on a parameter whose value is 44.2 Pa — roughly **45 % of the parameter itself, in a single step**. Nothing shipped profiles it today. **That is the risk the fix buys against: not a wrong number now, but a silently wrong one the first time anyone profiles a parameter with a wide box.**

## The case FOR fixing it upstream

1. The failure is **silent**. It returns a plausible-looking interval, wider than the truth, with no warning. A too-wide interval looks like an honest conservative one, which is why it survived this long.
2. The fix is **one line** (size the floor from the parameter, not its box) plus the endpoint verification item 13 already wrote and tested.
3. The correction has an **independent anchor**: a closed form that touches none of the walk, agreeing to 2.4e-9 on 20 records.
4. It is **cheap to verify per parameter** — the exposure ratio above says which of `dst_fit`'s numbers could move, before any code runs. This document is that check and the answer is *none of them*, which makes the change low-risk as well as low-value.

## The case AGAINST, which is not weak

1. **`dst_fit` is consumed unchanged on purpose**, so that `biorheo` cannot quietly reshape the founder's package to suit itself. That rule has more value than one correct interval.
2. **`biorheo`'s wrapper already gives `biorheo` correct intervals.** Patching upstream buys correctness for OTHER callers of `dst_fit` — of which there is currently one, `dst_fit` itself, and Step 2 shows it is not affected.
3. `tests/test_gate_0_1_dstfit_unchanged.py` pins a figure hash that **will move**, and a gate whose whole purpose is to detect change cannot be left to fail quietly. It must be updated in the same commit, which means the change is not one line in practice.
4. **`dst_fit`'s published numbers were computed with the defective walk.** Step 2 says they do not move, but 'do not move' is a claim that has to be re-established after the change, not assumed from before it.

## Recommendation, for the founder to accept or reject

**Fix it, but not because anything is currently wrong.** Step 2 establishes that no shipped `dst_fit` number moves; the argument is entirely about the *next* caller, and the failure mode is a wrong answer with no error. If it proceeds:

1. **Pre-register what 'unchanged results' means before touching code**, per parameter, via the exposure ratio — the table in Step 2 is that pre-registration for `phi_m` and predicts NO movement at all. A gate that predicts no movement and then sees some has found something.
2. `test_gate_0_1_dstfit_unchanged.py` moves **in the same commit**.
3. **Do not delete `biorheo`'s wrapper.** It carries the endpoint verification, the second nuisance start and the `step_rule` / `exposure` reporting, none of which is part of the one-line fix and all of which caught something real — including two records where the profile is **not a function** of the parameter and a root finder returns a discontinuity as if it were a limit.

**If it is rejected, that is a defensible position** and the thing to do is record it here, so item 17 stops being reopened. What is not defensible is leaving the question undecided for a fifth item.
