# PRE-REGISTRATION — the strain clock: frictionless prep → friction on → constant-volume shear

**Written 2026-09-13, before any cell of this campaign ran.** Cited by
`derivations/alpha_from_dem.md` §4 and `cv_unloading.md` §6 as
`PREREG_STRAINCLOCK_2026-09-12.md`; those two references are corrected to this filename
in the same commit. The friction-switch code it depends on was written 2026-09-12 and its
smoke test is §3 below.

Supersedes nothing. It is the first campaign in this project whose preparation is the
physically right one for a dense suspension at rest.

---

## 1. Why this protocol and not the one we have been running

Every constant-volume campaign on disk prepares the pack **frictionally** — grains are
inserted, compressed to `phi_target` at `mu_p` = 0.3, relaxed, then sheared. That pack
already carries a frictional force network before the shear starts, and
`cv_unloading.py` measured how much: `p_0/p_iso` = **3.6–9.3**, i.e. the preparation sits
three to nine times above the servo-held isotropic state at the *same* `phi`. `p(phi)`
for a frictional pack is a history, not a function, and U1's "the pack is contractive"
is a statement about a **rate-compressed** pack.

A dense suspension at rest is not that. Below `sigma*` the contacts cannot transmit
tangential force, so the pack at rest is **frictionless** at its own `phi`. For
`phi_mu(CV) < phi < phi_0(CV)` a frictionless pack is *below its own jamming point* and
carries **no** grain pressure at all — the `OCR -> infinity` limit. The first frictional
shear must therefore **build** `p` from zero. That building takes a finite strain, and
that strain is the quantity this campaign measures.

The measured premise, not an assumption: the smoke cell (§3) prepared frictionless at
`phi` = 0.605, CV 0 — above `phi_mu` = 0.591, below `phi_0` = 0.632 — carries
`p_hat_prep` = **6.66e-11** with **5** real contacts and a rattler fraction of **1.0**,
and frictional shear takes it to `p_ss` = 3.76e-03. **7.75 decades**, on the channel
every prediction below is read from.

### What it decides that no steady-state measurement can

`Delta phi_c`(CV) is resolved (+0.0533, `PREREG_CONSTP_RATEMATCH_2026-09-12.md` §5), and
it is **degenerate**: "polydispersity shifts the critical line" (M2a) and "polydispersity
changes how fast a force network assembles" (M4) predict the *same* steady state. They
differ in the **transient**. If CV acts only through `phi_c`, then the clock measured at
matched distance from the critical line must be the same at every CV. If CV also changes
the network-building kinetics, it will not be. G4 below is that test, and it is the
first two-way discriminator this project has had between M2a and M4.

There is a prior pointing at "does not collapse": `cv_unloading.py`'s U4 found that at
matched `phi` the CV 0.5 pack relaxes its overpressure **2–3x faster** than CV 0. That is
a different protocol (relaxation from above, not building from zero) so it is a prior and
not a prediction, but it is the reason the G4 bar is set where it is.

---

## 2. The grid, fixed here

45 cells: lognormal, `mu_p` = 0.3 during shear, `mu_prep` = 0.0, **N = 1000**, 3 seeds,
3 CV levels, 5 rungs.

`phi_mu(CV)` is the **constant-pressure** `phi_c(p* = 3e-04)`, rate-matched
(`out/constp_ratematch.txt`); `phi_0(CV)` is the **isotropic frictionless** `phi_J`
(`out/dst_window.txt`, lognormal, N = 1500). The rungs sit at fixed
`x = (phi - phi_mu)/(phi_0 - phi_mu)` = **-0.15, 0.10, 0.35, 0.60, 0.85**:

| CV | `phi_mu` | `phi_0` | `W` | rungs |
|---|---|---|---|---|
| 0.0 | 0.591127 | 0.632371 | 0.0412 | 0.5849, 0.5953, 0.6056, 0.6159, 0.6262 |
| 0.2 | 0.601989 | 0.645392 | 0.0434 | 0.5955, 0.6063, 0.6172, 0.6280, 0.6389 |
| 0.5 | 0.644447 | 0.690401 | 0.0460 | 0.6376, 0.6490, 0.6605, 0.6720, 0.6835 |

`W` moves only 11 % over the ladder, so the **absolute** offsets `phi - phi_mu` are
matched across CV to within 11 % as well (-0.0062/0.0042/0.0145/0.0248/0.0351 at CV 0
against -0.0068/0.0046/0.0161/0.0276/0.0391 at CV 0.5). Both candidate collapse
coordinates are therefore testable on this one grid and **both will be reported**.

**Rate.** `gamma_max` = 3.0 at `gdot` = 1.17e-04, the rate the constant-pressure ladder
was matched at, obtained by setting `shear_steps` per CV from the measured `dt`:
CV 0 -> 102 000, CV 0.2 -> 200 000, CV 0.5 -> 480 000. `gdot` is an OUTPUT of the packing
in this code and `dt` varies by seed at CV 0.5 (0.0522–0.0546), so the realised spread
will be reported and G0 requires it below 1.2x.

**Sampling.** `n_samples` = 1500, i.e. `Delta gamma` = 0.002. The clock is expected at
`gamma_c` ~ 0.01–0.2 (U4's `gamma_half` ran 0.030–0.135 on the *relaxation* problem), so
the smallest resolvable clock carries ~5 samples. Sampling is passive in this driver — it
reads the state, it does not feed the timestep — so unlike `MEASURE_NSAMP` in the
Basilisk solvers it cannot drive the trajectory; the only cost is wall time.

### Three gauge caveats, stated before the run rather than after

1. **`phi_mu` and `phi_0` are measured in different pressure gauges.** `phi_mu` is
   `phi_c` at a finite `p*` = 3e-04; `phi_0` is a `p -> 0` extrapolation. `W` is
   therefore a mixed-gauge width and `x` is a convenience coordinate, not a measured
   fraction of a physical window. Everything graded below is graded on `phi - phi_mu`
   as well, which uses one gauge only. (Whether the `dphi/dln p` slope is itself
   CV-dependent is a separate, cheap question answered from the isotropic ladders in the
   same commit; if it is, `W(CV)` inherits it.)
2. **`phi_0(CV = 0.5)` carries `sd` = 0.0047**, so the top rung (0.6835, nominally
   0.0069 below `phi_0`) may in fact sit **above** the frictionless jamming point, where
   the prep is jammed and `p_hat_prep` is not ~0. That is detectable per cell rather than
   arguable: G1 reports `p_hat_prep` for every cell and any cell above the threshold is
   **excluded from the collapse and reported as excluded**, not counted as a failure.
3. **N = 1000, not 1500.** It matches the constant-pressure ladder `phi_mu` came from,
   and it is what a ~600 s cell allows. No finite-size claim is made.

---

## 3. The blocker, cleared before this was written

`check_prepare_refactor.py --sidecar shear_head_src.txt --old-has-prepare --fprep-smoke`,
run in WSL 2026-09-13 (`out/fprep_check.log`, exit 0):

- **the default path is BIT-IDENTICAL** to HEAD on all five continuous channels
  (`p_ss_hat`, `tau_ss_hat`, `mu_macro`, `steady_split`, `i_ss`), so the `mu_prep` and
  `n_samples` additions changed no trajectory that is already on disk;
- the **non-default arm is live and verified**: `friction_switch_ok = True`
  (materials and all real interactions read back `tan(atan(mu_p))`),
  `protocol = "cv-shear-v3-fprep"`, 607 samples recorded at `n_samples` = 600, and
  `p_ss_hat` differs from the default arm's (ratio 1.246) — a knob that returned the
  default's number would be inert rather than verified.

---

## 4. Definitions, fixed here

- **`p_prep`** — `p_hat` after stage C, before the friction switch and before any shear.
- **`p_ss`** — the driver's existing steady-state mean over the last `STEADY_FRAC` of the
  strain. Unchanged from every other campaign.
- **`gamma_c`** — the strain at which `p_hat` **first** reaches `p_prep + 0.5*(p_ss - p_prep)`,
  by linear interpolation between the two bracketing samples. **First** crossing, and the
  series is required to be monotone up to it within a tolerance (below) — the 2026-09-03
  `t_half` lesson: a crossing time whose monotonicity is assumed is a measurement of the
  rule, not of the material. Reported beside it: `gamma_c_last` (last up-crossing) and
  the count of down-crossings before the argmax. If `gamma_c` and `gamma_c_last` differ
  by more than 25 %, the cell's clock is **not readable** and is excluded with its
  numbers printed.
- **`gamma_90`** — the same at 0.9 of the rise, reported as a shape control: a clock that
  moves under the choice of level is a readout of the level.
- **the collapse statistic** — `R(CV) = geometric mean over the four rungs above
  `phi_mu` of `gamma_c(CV)/gamma_c(CV = 0)` at matched `x`, with the seed-to-seed
  scatter propagated.

---

## 5. Predictions, with bars. Written before the first cell.

**G0 — the campaign is what it says it is.** 45/45 records with
`protocol = "cv-shear-v3-fprep"`, `friction_switch_ok = True`, `mu_prep` = 0.0,
`n_samples` = 1500; realised `gdot` spread across CV below **1.2x**; max `I` below
**1e-02**. Any failure here voids the campaign rather than being reported as physics.

> **Amendment, written 2026-09-13 03:25, after the FIRST cell of 45 landed and before
> any other.** That cell is CV 0, `x` = -0.15, seed 1 (`phi` = 0.5849), i.e. a rung
> *below* the frictional critical line, and it returns `I` = **1.26e-02**, over G0's
> 1e-02 bar. This is not a rate problem and no rate fixes it: `I = gdot*d/sqrt(p/rho)`
> and `p` at an unjammed rung is ~0 by construction (this cell's `p_hat` = 8.7e-05
> against the 3e-03 of a jammed one), so `I` diverges as the rung unjams. **G0's `I` bar
> is therefore unmeetable at `x` = -0.15 at any rate, and it stays exactly as written.**
> The grader will print **both** numbers: the all-cells max `I` as registered (expected
> FAIL), and the max `I` over the cells that carry the clock (`z_nr_ss` >= 4, the
> criterion G3 already contains, no new bar). The second is the stricter statement,
> because a below-`phi_c` rung failing a quasi-static bar is the campaign's own G3
> prediction restated, not a defect in it. No claim below rests on an `x` = -0.15 cell
> except G3 and G5's lower bracket, and G5's bracket uses `p_ss`, not `I`.

**G1 — the premise: the prep carries no pressure.** `p_hat_prep` < **1e-06** (three
decades below the `p_ss` scale of ~3e-03) on every rung with `x` <= 0.60, on all three
CVs. PASS requires 36/36. The `x` = 0.85 rungs are **reported, not graded** (caveat 2).
A failure here is the interesting kind: it would mean the frictionless pack jams below
the isotropic frictionless `phi_J`, which contradicts `phi_0`'s own construction.

**G2 — the clock exists.** `p_ss/max(p_prep, 1e-12)` > **100** on every jammed cell, and
`p_hat` rises monotonically (allowing a 10 % non-monotone tolerance on the sample-to-sample
difference) up to `gamma_c`. PASS requires the rise on >= 90 % of admitted cells.

**G3 — below `phi_mu` nothing jams.** At `x` = -0.15, `z_nr_ss` < **4.0** (the project's
own isostatic bracket) on all three CVs, and at `x` >= 0.35, `z_nr_ss` >= **4.0**. This is
the internal check that the ladder is placed where it claims to be; if it fails, `phi_mu`
measured at constant pressure is not the constant-volume critical fraction and that is a
finding about the two protocols, not about the clock.

> **Amendment 2, written 2026-09-13 03:55, with 31 of 45 cells landed and the G4
> statistic not yet computed.** G3's isostatic bracket is written above as **4.0**, which
> is `d + 1` — the **`mu_p -> infinity`** limit. This campaign shears at `mu_p` = 0.3,
> where this project has *measured* the isotropic contact number at **4.535 ± 0.072**
> (`out/ziso_audit.txt`). Using the frictionless-limit constant at a finite friction is a
> lesson already in this lab's notes and I repeated it. **The bar is not moved.** 4.0
> stays, its verdict stays printed (it FAILS: the `x` = -0.15 rungs return `z_nr` up to
> 4.46 at CV 0 and 4.32 at CV 0.2, i.e. above 4.0 but below 4.535), and the companion at
> the measured value is printed beside it, where both halves hold on every CV with cells.
> The companion is the stricter test of the half that matters — it makes "`x` >= 0.35 is
> jammed" harder to pass, and that is the half the clock is read on.
>
> Recorded at the same moment, from the same partial data and *not* a change to anything:
> the clock is **not readable** at `x` = -0.15 and `x` = +0.10 (first and last
> half-crossings differ by 5-15x, tens to hundreds of down-crossings before the argmax),
> and it is cleanly readable at `x` >= 0.35 (first crossing == last, 0 span). That is the
> readability guard of §4 doing exactly what it was written for, and it means **G4 will
> run on three rungs, not four**. Whether the near-critical unreadability is physics (the
> pressure fluctuates by its own size near `phi_c`) or a sampling limit is a question for
> G6, not a reason to admit those cells.

**G4 — THE MECHANISM BIT. Does the clock collapse on the distance from the critical
line?** Computed at matched `x` and, separately, at matched `phi - phi_mu`.

- **COLLAPSE (M2a — CV acts only through `phi_c`)** if `R(0.5)` and `R(0.2)` both lie in
  **[0.67, 1.50]** and neither trends monotonically with rung.
- **DOES NOT COLLAPSE (M4 — CV changes the network-building kinetics)** if `R(0.5)` lies
  outside that band. U4's prior says the direction would be `R(0.5)` < 1 (faster).
- **UNDETERMINED** if the propagated seed scatter puts the 1.50 bar within **2 SE** of
  `R(0.5)`. This is a real verdict and will be reported as one rather than rounded to
  either side. The discriminator's own separation will be printed beside the verdict:
  if the two candidates' predictions are not at least **1.6x** apart on the data actually
  collected, the test is declared NOT INFORMATIVE regardless of which side `R` falls.

**G5 — the cross-protocol anchor, and it is not circular.** `p_ss(phi)` interpolated
across the five rungs crosses `p*` = 3e-04 at some `phi`; that `phi` must agree with
`phi_c(p*)` measured by the **constant-pressure** driver to **0.010** at each CV. The
ladder placement used `phi_mu`, but the crossing is read from measured `p_ss` values and
is **bracketed** by rungs 1 and 2 by construction (`x` = -0.15 and 0.10), so it is an
interpolation and not the window-bounded extrapolation the 2026-09-09 root-pinning
pre-registration was about. If the crossing falls outside the ladder it is REFUSED and
reported as refused.

**G6 — reported, not graded.** `gamma_90/gamma_c` per cell (shape); `p_ss` against
`cv_unloading.py`'s frictionally-prepared `p_ss(phi)` at the same `phi` — the clean
version of U2, since these two protocols reach the same `phi` from opposite preparations
and the critical state is *defined* by forgetting which; `rattler_frac_prep` and
`z_nr_prep` per rung (both expected to be 1.0 and NaN below `phi_0`, which is not a bug).

---

## 6. What each outcome buys, written now so it cannot be chosen later

| outcome | what it means | what it changes |
|---|---|---|
| G4 COLLAPSE | polydispersity enters the transient only through `phi_c` | M4 is not needed; `ID`'s `phi - phi_c` pole carries the whole CV dependence and `MECHANISM_2026-09-08.md` §6a can be written on one axis |
| G4 NOT COLLAPSE, `R(0.5)` < 0.67 | a broad pack builds its network in less strain at matched distance | M4 is live and separable from M2a; the drop-side consequence is that the CV 0.5 drop reaches its jammed state EARLIER in the impact, which is a time the high-speed camera can see |
| G4 NOT COLLAPSE, `R(0.5)` > 1.50 | a broad pack is SLOWER | same structurally, opposite sign; and it would contradict U4's relaxation prior, which would itself be worth reporting |
| G4 UNDETERMINED | the seed scatter is the size of the effect | more seeds at the two rungs with the most leverage, decided by the measured scatter, not more rungs |
| G1 FAIL | the frictionless prep is not the OCR -> infinity limit | the whole mechanism paragraph in `cv_unloading.md` §6 is withdrawn |
| G3 FAIL | constant-pressure `phi_c` is not the constant-volume critical fraction | a protocol result, and it would put `Delta phi_c` = +0.0533 back in question |

---

## 7. OUTCOME

**Landed 2026-09-13 04:20. 45 cells, 0 failures, ~2100 s per 3-cell shard on 15 workers.**
Every number below is printed by `grade_strainclock.py` -> `out/strainclock.txt`; the
verdict line is `FAILURES: ['G0', 'G2', 'G3']`, and all three failures are the
pre-registered bars catching something real rather than the campaign misfiring.

### 🔑 G4 — the mechanism bit: COLLAPSE. CV acts through `phi_c`, not through the clock.

| | `R` = `gamma_c(CV)/gamma_c(0)` at matched `x` | per-rung | verdict |
|---|---|---|---|
| CV 0.2 | **1.069 ± 0.092** | 1.214, 0.987, 1.019 | COLLAPSE |
| CV 0.5 | **0.914 ± 0.103** | 1.013, 0.795, 0.948 | COLLAPSE |

and the test is **informative**, which had to be established separately: the two
candidates predict `R` = 1.00 (M2a) against `R` = 0.33–0.50 (M4, using `cv_unloading`'s
U4 prior that a CV 0.5 pack relaxes 2–3x faster at matched `phi`), i.e. **2.0–3.0x apart**
against the 1.6x rule. The measurement sits **0.8 se** from M2a and **4.0 se** from the
M4 band. **A kinetic effect of the size U4's relaxation result suggested is excluded.**
One smaller than the seed scatter is not: the 2 se bound on `R`(0.5) is **[0.71, 1.12]**,
so a CV effect on the clock larger than about 30 % in either direction is ruled out and
anything smaller is not resolved by three seeds.

⚠ **A code defect found while reading this, and it is reported rather than quietly
fixed.** §5's separation rule says "the two candidates' PREDICTIONS must be at least
1.6x apart". The first coding computed the *span of the per-rung ratios* (1.27x) — a
measure of rung-to-rung scatter, not of candidate separation — and on that number
declared the test NOT INFORMATIVE. Both numbers are now printed; the registered
*sentence* is graded on the candidates' predictions, which is what it says.

⚠ **Where the collapse is NOT established.** The `x` = +0.10 rung — nearest the critical
line, and the rung where a kinetic difference was predicted to be *largest* — was lost at
all three CVs to unreadability (below). So the collapse rests on `x` = 0.35, 0.60, 0.85
and says nothing about the approach to `phi_mu`.

### 🔑 G5 — the same `phi_c`, from a protocol with nothing in common but the material

`p_ss(phi)` interpolated across the ladder crosses `p*` = 3e-04 at

| CV | crossing | constant-pressure `phi_c(p*)` | \|d\| |
|---|---|---|---|
| 0.0 | 0.5876 | 0.591127 | 0.0035 |
| 0.2 | 0.5987 | 0.601989 | 0.0033 |
| 0.5 | 0.6399 | 0.644447 | 0.0046 |

**PASS at a bar of 0.010**, and the contrast is better than any absolute:
**`Delta phi_c`(CV 0 -> 0.5) = +0.0523 from these crossings against +0.0533 from the
constant-pressure ladder, agreeing to 0.0011.** The two protocols share no machinery — one
is constant PRESSURE with a volumetric servo, the other constant VOLUME from a
FRICTIONLESS preparation — and the residual offsets are one-sided (all three crossings sit
*below* their constant-pressure value), which is why they cancel in the contrast. The
interpolation is bracketed by rungs 1 and 2 by construction, so this is not the
window-bounded extrapolation `PREREG_ROOT_PINNING` was about.

### G1 PASS — the premise holds, by 7 to 8 decades

`p_hat_prep` < 1e-06 on **36/36** graded rungs, running down to **2.0e-11**, with
`p_ss/p_prep` up to **3.1e+08**. A frictionless pack between `phi_mu(CV)` and `phi_0(CV)`
carries no grain pressure, and frictional shear builds the whole of it. `z_nr_prep` is NaN
and `rattler_frac_prep` = 1.0 on 34 of 45 cells, which is what "unjammed" looks like and is
not a defect. Of the 9 ungraded `x` = 0.85 rungs, **2 were prepared JAMMED**
(`p_prep` up to 1.96e-04) — exactly the `phi_0(CV 0.5)` uncertainty §2 caveat 2 flagged in
advance, detected per cell rather than argued about.

### G0 FAIL — on the bar that could not be met, exactly where it was predicted to fail

45/45 records carry `protocol = "cv-shear-v3-fprep"`, `friction_switch_ok = True`,
`mu_prep` = 0, `n_samples` = 1500, N = 1000; `gdot` spread across CV is **1.001** against
a 1.2 bar. The failure is the `I` bar: **1.46e-02** against 1e-02, and every offending
cell is an `x` = -0.15 rung whose `p_ss` is ~1e-04, i.e. the amendment's point made
by the data. Over the cells above the MEASURED `z_iso` = 4.535 the max `I` is reported
beside it.

### G2 FAIL — and the failure is a physical result, not a defect

Readable on **27 of 43** jammed cells, monotone rise on 31 of 43 against a 90 % bar. The
split is the whole story:

| rung | readable | monotone |
|---|---|---|
| `x` = -0.15 | 0/7 | 0/7 |
| `x` = +0.10 | 0/9 | 5/9 |
| `x` = +0.35 | 9/9 | 8/9 |
| `x` = +0.60 | 9/9 | 9/9 |
| `x` = +0.85 | 9/9 | 9/9 |

(the `x` = -0.15 row is out of 7, not 9: two of its nine cells fall below the `z_nr` = 4.0
admission the G2 table is taken over — they are counted in G3, which is where they belong)

Near the critical line the grain pressure fluctuates by its own size: the first and last
half-crossings differ by 5–15x with tens to hundreds of down-crossings before the argmax,
so **a half-rise strain is not a well-posed observable there**. Away from it the series is
monotone to its crossing and first == last exactly. The readability guard of §4 is what
separated the two, and it was written before any cell ran.

### G3 FAIL as registered, PASS on the measured constant

At the registered bracket of 4.0 the `x` = -0.15 rungs return `z_nr` up to 4.46 (CV 0),
4.32 (CV 0.2) and 4.26 (CV 0.5) — above 4.0. At the **measured** `z_iso` = 4.535 for
`mu_p` = 0.3 both halves hold on all three CVs: every `x` < 0 cell below it, every
`x` >= 0.35 cell above it. Amendment 2 records why the bar was wrong and why it was not
moved.

### The prediction that did NOT survive, and it was one of the campaign's own

`NEXT_STEPS` §0'.1 prediction (ii) said `gamma_c` would **grow** as `phi -> phi_mu+`.
On the readable rungs it does not. The nine cell means run 0.1747/0.2121/0.1769 at
`x` = 0.35, 0.1750/0.1726/0.1390 at 0.60 and 0.1668/0.1699/0.1582 at 0.85 (CV 0/0.2/0.5) —
a total span of **1.53x with no systematic trend in either `x` or CV**, against a
within-cell seed sd of up to **37 %** (0.065 on 0.175). The spread is consistent with seed
scatter and there is no growth toward the critical line. What happens on the approach to
`phi_mu` instead is that the observable stops existing. The
divergence may still be there and be hidden by the fluctuations, and distinguishing those
two needs a different observable (an ensemble-averaged `p(gamma)` over many seeds, or a
first-passage time to a fixed `p`), not more rungs.

### What this changes

1. **M4 is demoted, on its own discriminator.** The clock does not depend on CV at matched
   distance from critical, so "polydispersity changes the network-building kinetics" has
   no support at the size U4 suggested. `MECHANISM_2026-09-08.md`'s decision tree loses a
   branch and **M2a is the surviving packing channel**.
2. **`Delta phi_c` is now a two-protocol number**, +0.0523 / +0.0533, which is the
   strongest form the project's central measurement has taken.
3. **The strain to build a force network from nothing is ~0.15 and universal here** across
   `phi` and CV. That is a number the drop side can use directly: at an impact strain rate
   of `U_0/D_0`, 0.15 of strain is 0.15 `D_0/U_0` — tens of microseconds at 1 m/s and
   2.5 mm, i.e. fast compared with the impact, so the network is built essentially
   instantly and the *steady* critical state is the right thing to compare `phi` against.
