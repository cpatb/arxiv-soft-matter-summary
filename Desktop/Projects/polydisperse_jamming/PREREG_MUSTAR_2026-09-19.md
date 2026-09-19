# PREREG — does polydispersity move the FRICTION SATURATION SCALE μ\*, or just the amplitude?

**Written 2026-09-19, BEFORE any cell of this campaign is launched.** Every prediction below
is computed from data committed before this session and is written down here so the
comparison cannot be made afterwards. Companion: `derivations/wc_friction_endpoint_2026-09-19.md`.

---

## 1. The question, and why it is now well posed

`wc_friction_endpoint_2026-09-19.md` measured, on the lognormal family, that the residual
friction gap above `μ_p = 0.3` **grows with CV**:

| CV | `φ_J(0.3) − φ_J(1.0)` |
|---|---|
| 0.0 | +0.01550 ± 0.00073 |
| 0.2 | +0.01542 ± 0.00068 |
| 0.5 | +0.02136 ± 0.00146 |

growth over CV 0 → 0.5, paired both ways: **+0.00586 ± 0.00077, \|t\| = 7.65**.

That page reported this as "a wider size distribution needs more friction to reach its
frictional limit". **There are two mechanisms that produce it and the page could not
separate them.** Write the saturation as

```
   φ_J(μ) = φ_∞ + A exp(−μ/μ*)
```

- **(M1) SCALE.** `μ*` grows with CV — the approach to the frictional limit is genuinely
  slower for a wider distribution.
- **(M2) AMPLITUDE.** `μ*` is CV-independent and `A` grows with CV — the whole friction drop
  is simply bigger, so a fixed fraction of it is a bigger absolute gap.

These are physically different. (M1) says polydispersity changes how friction mobilises;
(M2) says it only changes how much there is to mobilise.

---

## 1b. 🔴 AMENDMENT 1, written BEFORE launch: the channel changes from `phi_J_lin` to a matched rung

A referee pass on the companion page found that **the high-friction ladders fail the
project's own ladder-health criterion**, and I verified it before writing this:

| `μ_p` | non-monotone ladders | worst `r²` |
|---|---|---|
| 0.0 | 7 / 18 | — |
| 0.3 | 7 / 18 | — |
| **1.0** | **9 / 9** | **0.003** |
| **3.0** | **3 / 3** | 0.005 |

At `μ_p = 1.0`, CV = 0.5, seed 1 the ladder reads 0.61740 / 0.61388 / 0.62075 / 0.61264 —
that is noise, not a ladder. **`phi_J_lin` at `μ_p ≥ 1.0` is a fit to noise**, and §2's
original predictions were built on it.

**The mechanism, measured:** the ladder's dynamic range collapses. `φ` span across the four
rungs, against the seed scatter at the bottom rung:

| `μ_p` \ CV | 0.0 | 0.2 | 0.5 |
|---|---|---|---|
| 0.0 | 16.8 | 7.8 | **1.3** |
| 0.3 | 17.7 | 48.7 | 4.0 |
| 1.0 | 6.3 | **1.7** | **1.5** |
| 3.0 | **2.1** | — | — |

A high-friction packing is nearly incompressible over the accessible pressure range, so
there is nothing left to extrapolate. **Wherever span/sd ≲ 2 the extrapolation is not
measuring a jamming point.** Note this also catches the *frictionless* CV = 0.5 cells.

**What changes in this prereg, all of it before any cell runs:**

1. 🔴 **The primary channel is now `φ` at the matched lowest rung `p̂ = 3e-5`, fit-free.**
   `phi_J_lin` is still recorded and reported, but it is not graded at `μ_p ≥ 1.0`.
2. **The §2 fits and the §4 predictions are recomputed on that channel** and replaced below.
   The old `phi_J_lin` values are kept struck through so the two can be compared.
3. **A NEW GATE, F0, is added**: every `(μ_p, CV)` group must report its span/sd, and any
   group below 2.0 is reported as **not readable on the extrapolated channel** regardless of
   what its `phi_J_lin` says.

**Recomputed on the matched rung** (`φ` at `p̂ = 3e-5`, mean over 3 seeds):

| CV | `μ_p` = 0.0 | 0.3 | 1.0 | 3.0 |
|---|---|---|---|---|
| 0.0 | 0.63182 | 0.57535 | 0.56148 | 0.56054 |
| 0.2 | 0.64202 | 0.58755 | 0.57316 | — |
| 0.5 | 0.69084 | 0.63667 | 0.61486 | — |

| CV | `φ_∞` | `A` | `μ*` |
|---|---|---|---|
| 0.0 | 0.56105 | 0.07077 | **0.1877** |
| 0.5 | 0.61346 | 0.07737 | **0.2491** |

`μ*` ratio **1.327**, `A` ratio **1.093** — against 1.243 and 1.031 on the discredited
`phi_J_lin` channel. **The conclusion is robust to the channel; only the numbers move.**

**Revised blind predictions, on the matched rung:**

| prediction | value |
|---|---|
| `φ(μ_p = 0.6, CV = 0.0)` | **0.56394** |
| `φ(μ_p = 0.6, CV = 0.5)` | **0.62042** |
| `φ(μ_p = 3.0, CV = 0.5)` | **0.61346** |
| `φ(μ_p = 3.0, CV = 0.0)` | 0.56105 predicted vs **0.56054** measured — *not blind*; the form already misses by 0.00051, which sets the scale |

> **REGISTERED: each blind prediction lands within ±0.0030 of the measured mean**, four times
> the 0.00051 the form misses by where it can be checked. Unchanged from §4; only the values
> move.
>
> **REGISTERED, unchanged: `μ*(0.5)/μ*(0)` in [1.05, 1.45]** — the matched-rung estimate 1.327
> sits inside the band that was written down from the `phi_J_lin` estimate 1.243, so the band
> is not being moved to fit the new channel.

---

## 2. What the existing points already say, and why that is not enough

⚠ **Superseded by §1b** — kept because §1b's band was set from it and that has to remain
visible. The numbers here are on the `phi_J_lin` channel, which §1b discredits at `μ_p ≥ 1.0`.

Fitting the three-parameter form to what is on disk:

| CV | points | `φ_∞` | `A` | `μ*` |
|---|---|---|---|---|
| 0.0 | 4 (`μ_p` = 0, 0.3, 1.0, 3.0) | 0.55816 | 0.07421 | **0.1986** |
| 0.5 | 3 (`μ_p` = 0, 0.3, 1.0) | 0.61392 | 0.07649 | **0.2469** |

`μ*` ratio **1.243**; `A` ratio **1.031**. That points at (M1) and away from (M2).

🔴 **But the CV = 0.5 fit has THREE points and THREE parameters — zero residual degrees of
freedom. It cannot self-diagnose, and any three points fit it exactly.** This is the same
defect the bidet grid study was caught on (a 3-point Richardson fit has no residual dof).
The whole purpose of this campaign is to add points so the CV = 0.5 fit becomes
over-determined and its residual becomes readable.

## 3. What is being run

| | |
|---|---|
| runner | `phase2_phiJ.py`, **unchanged** |
| family | `lognormal` (the only family with any friction level above 0.3) |
| friction | `μ_p ∈ {0.6, 3.0}` — the two values missing from the ladder |
| CV | 0.0, 0.2, 0.5 |
| seeds | 1, 2, 3 |
| `p_hat` | the production ladder 1e-3, 3e-4, 1e-4, 3e-5 |
| N | 1500 |
| tag | **`mustar`** — matches no existing tag |
| records | **72** = 2 `μ_p` × 3 CV × 3 seeds × 4 rungs |

`(μ_p = 3.0, CV = 0)` already exists and is **deliberately re-run** as the wiring gate.

**Cost.** High friction is cheap on this protocol: the existing `μ_p = 3.0, CV = 0` cells
have a median `wall_s` of **67 s** against 487–724 s at `μ_p = 0`. Budgeting the new cells at
the `μ_p = 1.0` rates (107 s at CV 0, 350 s at CV 0.2, 561 s at CV 0.5) gives ≈ 4 CPU-hours,
or **≈ 25–40 min** on 12 workers. **This campaign is cheap because friction makes packings
jam sooner, which is the very effect being measured.**

**Launch discipline:** `scripts/launch_mustar.sh`, by file path from PowerShell, held open in
a background Bash call for the whole run (`wsl.exe -e bash -lc '... &'` returns immediately
and WSL kills the group). **Not launched until the `wcbidi` campaign finishes** — running
both would put 24 workers on 16 physical cores and the wall times above would stop meaning
anything.

## 4. Gates

### F1 — WIRING, and it can fail

The 12 new `(μ_p = 3.0, CV = 0)` records must be **bit-identical** on `phi`, `z_all`,
`z_nr`, `rattler_frac` and `n_contacts` to the existing ones at the same (seed, rung).
Existing extrapolated value to reproduce: **`φ_J(μ_p = 3.0, CV = 0) = 0.55742`**.

### F2 — ladder health

**Registered: 72 of 72 converge.** High friction jams sooner, so if anything caps here it is
a surprise and is reported as an amendment, not re-run at a raised cap.

### F3 — **PRIMARY: (M1) or (M2)?**

Refit both CV levels with the new points. The CV = 0.5 fit then has **5 points, 3
parameters, 2 residual dof** and its residual is readable for the first time.

> **REGISTERED: `μ*(CV = 0.5) / μ*(CV = 0)` lands in [1.05, 1.45].**
> (M1) predicts > 1; (M2) predicts = 1. The band is centred on the 1.243 the existing
> points give, widened to cover the extra scatter three new points will bring.
>
> **REGISTERED: `A(CV = 0.5) / A(CV = 0)` lands in [0.90, 1.20].**
> This is the discriminating half. If `A` moves as much as `μ*` does, the two mechanisms are
> not separated by this design and the verdict is UNDETERMINED regardless of what `μ*` does.

**Resolving power, stated in advance.** The `μ*` ratio's uncertainty will come from the fit
covariance with seed-level resampling; if its 95 % interval contains 1.0 the verdict is
**UNDETERMINED**, not (M2). A bootstrap over the three seeds is registered as the error
estimate, with the per-seed fits reported individually so a single bad seed is visible.

### F4 — blind point predictions

From the fits in §2, before the cells exist:

| prediction | value |
|---|---|
| `φ_J(μ_p = 0.6, CV = 0.0)` | **0.56178** |
| `φ_J(μ_p = 0.6, CV = 0.5)` | **0.62065** |
| `φ_J(μ_p = 3.0, CV = 0.5)` | **0.61392** |
| `φ_J(μ_p = 3.0, CV = 0.0)` | 0.55816 — *not blind*, this point is in the CV = 0 fit; its measured value is 0.55742, so the form is already 0.00074 off there, which sets the scale for what "agreement" means |

> **REGISTERED: each of the three blind predictions lands within ±0.0030 of the measured
> mean**, i.e. four times the 0.00074 the form already misses by at the one point where it
> can be checked. A miss is reported as a miss; the band is not widened afterwards.

⚠ The exponential form is a **choice**, not a derivation. Nothing in the literature resolved
so far pins the functional form of `φ_J(μ_p)`. F4 therefore tests the form and the
parameters together, and a failure does not distinguish them. F3's ratio is the more robust
statement because a common form error largely cancels between the two CV levels.

### F5 — the form itself, now testable

With 5 points and 3 parameters the CV = 0.5 fit gains 2 residual dof.

> **REGISTERED: the residual RMS of the 5-point fit is reported at both CV levels, and if it
> exceeds the mean seed s.e. of the points (≈ 0.0015) the exponential form is recorded as
> INADEQUATE and F3's ratio is quoted as form-dependent.**

A two-parameter alternative (`φ_∞ + A/(1 + μ/μ*)`, algebraic rather than exponential) is
named **now** as the comparison, and whichever fits better is reported — not chosen.

### F6 — the CV = 0.2 rung, registered as an interpolation check

CV = 0.2 gets the same two new `μ_p` values, so all three CV levels end with 5-point ladders.

> **REGISTERED: `μ*(CV)` is monotone in CV across the three levels.** The gap measurement
> that motivated this campaign was *non*-monotone in its raw form (+0.01550 / +0.01542 /
> +0.02136, i.e. flat then rising), so a monotone `μ*` would be a real strengthening and a
> non-monotone one is a genuine outcome worth reporting rather than a failure.

## 5. What would make this campaign worthless

- F1 fails → the runner does not reproduce its own committed cells and nothing here is read.
- `A` moves as much as `μ*` → the two mechanisms are not separated and the verdict is
  UNDETERMINED. This is registered as a likely outcome, not an unlikely one.
- The exponential form is inadequate on residuals at both CV levels → F3's ratio is a
  statement about a form nobody has justified, and the honest verdict is that the campaign
  measured `φ_J(μ_p)` on a 5-point ladder and did not settle the mechanism.

## 6. What this cannot decide, whatever it returns

- **One family.** Only lognormal has friction above 0.3, and this campaign does not change
  that. Whether `μ*(CV)` is a property of the lognormal family or of polydispersity as such
  is not addressed.
- **`μ_p → ∞` is still not reached**, only approached; `μ_p = 3.0` is the largest value on
  this protocol and F1's own reproduction is the only evidence it is converged.
- **N = 1500 throughout**, so nothing here bears on the finite-size flag.
