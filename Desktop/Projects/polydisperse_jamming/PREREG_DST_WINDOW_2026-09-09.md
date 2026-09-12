# The shear-thickening WINDOW — what sets the *type*, and four testable hypotheses

**Written 2026-09-09.** Sections 1–4 are analysis of artifacts already on disk
(`dst_window.py` → `out/dst_window.txt`, reading `out/phase2_fits.csv`). Section 5 is a
**pre-registration** written before its cells were launched. Nothing here is a new
simulation result yet.

Companion: `PREREG_COMMON_LADDER_2026-09-09.md`, the run that decides whether `φ_c` is a
material root or a property of its ladder.

---

## 0. The question, sharpened

The project's question has been *"why does raising CV move an impacting drop from
solid-like to spreading?"*, and every mechanism weighed so far answers it with **one**
packing fraction and the sign or size of `φ − φ_something`:

| | turns on | what it predicts |
|---|---|---|
| M1 permeability | `d₃₂(CV)` | a smooth factor in `ID` |
| M2a plain rheology | `φ_m − φ` | viscosity falls, smoothly |
| M2′ dilatancy sign | sign of `φ − φ_c` | pore pressure changes sign |

The Wyart–Cates account of shear thickening turns on **two**, and it is the position of
`φ` *between* them that fixes the **type** of thickening rather than its size:

```
    φ  <  φ_μ            CST      thickens, viscosity finite at every stress
    φ_μ <  φ  <  φ_0     DST / shear jamming: viscosity diverges at a FINITE stress
    φ  >  φ_0            jammed at every stress
```

`φ_0` is the frictionless (low-stress) jamming fraction, `φ_μ` the frictional
(high-stress) one. So there is a **window** `W = φ_0 − φ_μ`, and a reduced coordinate

```
    x = (φ − φ_μ) / W
```

This is the first time this project has asked what polydispersity does to `W` and to `x`
rather than to a single jamming point. It is the natural home for the phrase *"shear
thickening **type**"*, and it is answerable from ladders already run.

⚠ **The identification is an assumption and is flagged in the code.** Wyart–Cates' `φ_0`
and `φ_μ` are jamming points of the *sheared* suspension; what is measured here is the
**isotropic** `φ_J` at `μ_p` = 0 and `μ_p` = 0.3, because that is the pair with CV
coverage at both frictions. The frictionless leg is the safer of the two (frictionless
isotropic compression goes to random close packing, which is what `φ_0` means). The
frictional leg is **not** safe: this project's own shear campaign puts the frictional
critical state *above* the isotropic frictional `φ_J`. Every number below is therefore
*a window between two isotropic ladders*, and the Wyart–Cates reading of it is the
hypothesis, not the measurement.

---

## 1. The measured window, and the one number that matters

`out/dst_window.txt`, lognormal, N = 1500, three seeds, formed **per seed and then
averaged** — the two frictions were run on the *same* seeds and `cv_real` matches to
1.1e-16 (machine epsilon) on all 15 pairs, so the arms are paired and the paired
difference is the correct statistic. (The unpaired difference-of-slopes is printed beside it and is
conservative by the shared-realisation variance; both are in the artifact.)

| CV | `φ_0` | `φ_μ` | `W` | ±se |
|---|---|---|---|---|
| 0.0 | 0.633156 | 0.574533 | **0.058623** | 0.000284 |
| 0.1 | 0.637495 | 0.578113 | 0.059382 | 0.000795 |
| 0.2 | 0.645392 | 0.588050 | 0.057342 | 0.000653 |
| 0.3 | 0.658335 | 0.603985 | 0.054350 | 0.001512 |
| 0.4 | 0.672854 | 0.620182 | 0.052671 | 0.000229 |
| 0.5 | 0.683324 | 0.632565 | 0.050759 | single seed |

🔑 **The window is narrow and polydispersity slides it about twice as fast as it is
wide.**

```
    W(CV = 0)      = 0.0586
    dφ_μ/dCV       = +0.117 +- 0.016      = 2.00 windows per unit CV
    ΔCV to traverse the whole window at fixed φ = 0.50
```

That is the mechanism statement for *type*: a material sitting anywhere inside the DST
window is carried clean out of it — from shear-jamming-capable to CST-only — by a
**width change of ΔCV ≈ 0.5**, at fixed `φ`, fixed chemistry and fixed mean size. It is a
**threshold**, which is what the mechanism document says a smooth permeability factor
cannot produce, and it is reached without invoking pore pressure at all.

The window itself also narrows, and only for the skewed family:

| family | `dW/dCV` | t | verdict |
|---|---|---|---|
| lognormal | **−0.0169 ± 0.0037** | −4.6 | resolved: −12 % over CV 0 → 0.4 — ⚠ **but see §3(c)**: the CV = 0 anchor alone moves `W` by 0.0032 between runs |
| uniform (zero skew) | +0.0050 ± 0.0058 | +0.9 | **not resolved**; \|dW/dCV\| < 0.017 at 2 se |

⚠ **§5 pre-registers a run that may remove this row**, for reasons in §3. The traverse
number above does **not** depend on it: `dφ_μ/dCV` is an absolute slope with no window in
its denominator, and it is quoted alongside `x` for exactly that reason.

### ⚠ The window is a DEFINITION, and it costs a factor 1.29

Wyart-Cates' `φ_μ` is the **fully frictional** jamming point. Everything above reads it
at `μ_p` = 0.3, the production friction, which is not the saturated one. The project's own
friction ladder at CV = 0 prices the choice:

| `μ_p` | `φ_μ` | `W` | lab `x` | ΔCV to traverse |
|---|---|---|---|---|
| 0.3 | 0.574533 | 0.0586 | −2.29 | **0.50** |
| 1.0 | 0.559036 | 0.0741 | −1.61 | 0.63 |
| 3.0 | 0.557423 | 0.0757 | −1.55 | 0.65 |

`φ_J` saturates in friction (the `μ_p` 1 → 3 step is inside seed scatter), so the window is
**bounded between 0.0586 and 0.0757** and every conclusion here survives that spread *in
kind*: the traverse is 0.50–0.65 in CV either way, the lab sits 1.6–2.3 windows below.
**No single number should be quoted without naming the `μ_p` it was read at.**

⚠ **Not measured:** `dφ_μ/dCV` at the saturated friction. The `μ_p` = 1 and 3 ladders exist
only at CV = 0, so the traverse above pairs a possibly-saturated width with a
`μ_p` = 0.3 slope. An isotropic ladder at `μ_p` = 1 across CV is the cheapest remaining
gap in this section and should be run with H-1.

---

## 2. What this says about the lab's own material, and it is not comfortable

At the bench `φ` = 0.44:

| CV | `φ − φ_μ` | `x` | class | densification needed to ENTER / to JAM |
|---|---|---|---|---|
| 0.0 | −0.1345 | −2.29 | CST only | 1.306× / 1.439× |
| 0.2 | −0.1481 | −2.58 | CST only | 1.336× / 1.467× |
| 0.5 | −0.1926 | −3.79 | CST only | 1.438× / 1.553× |

Two things follow, and they pull in opposite directions.

**(a) It sharpens `MECHANISM` §6b from a gap into a requirement.** That section says the
lab is 0.13 below its critical state and needs 1.29× densification. The window version
says the same thing in units that carry more information: the lab is **2.3 window-widths
below the window**, and impact must supply **31 % of extra packing fraction (23 %
volumetric compaction) merely to reach the bottom edge, 44 % to be jammed at all
stresses**. The window is the natural yardstick because it is the only length in the
problem that is set by the *material* rather than by an arbitrary offset.

**(b) It is a point AGAINST reading the bench observation as a change of type.** At
`φ` = 0.44 the material is CST-only at every CV measured, and raising CV moves it
*further* below the window (−2.29 → −3.79), not across it. On these numbers the CV effect
at the lab's stated `φ` is **M2a — the viscosity falling as the material retreats from
jamming — and nothing more exotic**. A type change requires resolution (1), (2) or (3) of
`MECHANISM` §6b to be true; it is not licensed by the window analysis on its own.

Stated plainly: **the window analysis does not rescue the lab's working point. It prices
it.**

---

## 3. What is wrong with the frictionless arm, and it was never looked at

The frictionless ladders are the weak leg of every number above, in two ways that were
found by checking rather than assumed away.

**(a) Every frictionless cell is a THIN extrapolation.** `φ_J` is a linear fit in `p̂`
over the ladder; at CV ≥ 0.3 the frictionless fits use **2 points**, which has zero
residual degrees of freedom, reports `r²` = 1.000 by construction and cannot
self-diagnose. Those are exactly the rungs that carry the `dW/dCV` slope.

**(b) The censoring is on the axis under test and is ASYMMETRIC between the arms.**
`analyse_phase2` drops cells with `converged = False`. Split by arm — which the script's
own pooled SELECTION CHECK cannot do — the drop rate is:

| family | `μ_p` | CV 0 | CV 0.3 | CV 0.4 | CV 0.5 |
|---|---|---|---|---|---|
| lognormal | 0.0 | 0.92 | 0.58 | 0.42 | **0.33** |
| lognormal | 0.3 | 1.00 | 1.00 | 0.75 | 0.75 |
| uniform | 0.0 | 0.83 | 0.67 | 0.58 | **0.17** |
| uniform | 0.3 | 1.00 | 1.00 | 0.92 | 0.50 |

So `W = φ_0 − φ_μ` at CV = 0.5 differences an arm censored 67 % against one censored
25 %. Three further facts, all measured:

- **All 82 non-converged cells stopped at exactly `MAX_ITER` = 200 000.** They are slow,
  not wrong. Median `unbalanced` = 0.046 against a tolerance of 0.02, i.e. a factor 2.3.
- **18 of the 82 met the physical bar** (`unbalanced` < 0.02) and were rejected *only* on
  `O.iter < MAX_ITER`; **12 of those 18 are in the starved frictionless arm.** A
  fail-closed criterion rejecting valid states, which is a shape this repo has met before.
- **The censoring biases `φ` in the direction that inflates the finding.** Per (cv, p̂)
  cell, keeping only converged members shifts mean `φ` by a median 0.0007, worst 0.0042 —
  and the shift's slope in the frictionless lognormal arm is **−0.0047 per unit CV**.
  Over CV 0 → 0.4 that manufactures **0.0019 of the measured 0.0068 narrowing, about
  28 %, with the same sign**.

**(c) 🔴 And the protocol is NOT deterministic at CV = 0 — found 2026-09-09 by a
disagreement between two scripts.** `analyse_phase2.load_all` deduplicates shard files
by cell key, with a docstring saying duplicates are "genuinely redundant" because "the
protocol is deterministic given the seed". It is not. **15 duplicated keys return a
different `φ` and a different iteration count from the same (cv, μ_p, seed, `p̂`, N) and
the same `protocol` string** — worst `|Δφ|` = **5.0e-03**, median 2.1e-03, against a seed
scatter of 2.2e-04 at that cell. **All of them are at CV = 0.0**, the monodisperse and
most degenerate packing, and none of them disagrees about `converged`.

Which record survives is set by the **filename sort order**. Enumerating every choice:

```
   phi_J(CV 0, mu_p 0  ) spans 1.7e-03 to 2.2e-03 per seed
   phi_J(CV 0, mu_p 0.3) spans 3.0e-04 to 1.8e-03 per seed
   => W(CV = 0) ranges 0.05594 to 0.05911     span 0.0032
      against the quoted standard error of   0.00028   -- a factor 11
```

**So `W`(CV = 0) = 0.0586 ± 0.0003 understates its own reproducibility by an order of
magnitude**, and CV = 0 is the leverage end of every slope in §1. `analyse_phase2` now
counts the disagreements and prints the worst instead of silently keeping the last; the
CSVs regenerate byte-identically, so no number moved when the detection was added.

⚠ **What this does and does not reach.** `Δφ_J`(CV 0 → 0.5) = +0.0623 moves by at most
0.0022, ~3.5 %, and survives. `dφ_μ/dCV` = +0.117 likewise. **`dW/dCV` carries the leverage at CV = 0**, and the
size of that was GUESSED here at "~40 % of the total change" and then computed in §7:
moving `W`(CV = 0) across the 0.0032 spread and refitting gives **±0.0023 against a fit
error of ±0.0035**, i.e. smaller than the fit error and about 4× less than the guess. The
slope survives. Left in place as written because the guess was the wrong way round and
the correction is the point. The `capfix` re-run re-measures CV = 0 at both frictions and will supply
a run-to-run replicate rather than a seed replicate. **And it does not touch the shear
campaign at all** — different protocol, different code path, and
`grade_common_ladder`'s duplicate-state check reports 0 offending pairs there.

**None of this touches `dφ_μ/dCV` or the traverse**, which live entirely in the
frictional arm (censored 0–25 %) and are the numbers §1 quotes as the mechanism.

---

## 4. Where the four mechanisms now stand, and what separates them

| | mechanism | moves with CV through | signature nothing else has |
|---|---|---|---|
| M1 | permeability | `d₃₂²`, a smooth factor | none — excluded as the weaker channel |
| M2a | plain rheology | `φ_m − φ`, smooth | **no threshold, no boundary in (U₀, D₀), no strain history** |
| M2′ | dilatancy sign | sign of `φ − φ_c` | **the SIGN of the pore pressure** |
| **M4** | **window traverse** | **`x = (φ − φ_μ)/W`, a threshold at `x` = 0 and 1** | **the SHAPE of the flow curve: S-shaped ⇄ monotonic** |

M4 is not a rival to M2′ so much as its dry twin: both need the material to sit near a
packing threshold, and both predict a *sharp* change rather than a factor. What separates
them is where you look — M2′ needs a pressure transducer in the substrate, **M4 needs a
rheometer the lab already has.**

---

## 5. 🔒 PRE-REGISTRATION — four hypotheses, in cost order

Written before the cells in H-1 were launched. `scripts/launch_capfix.sh` and
`grade_capfix.py` are named here and did not exist before this file.

### H-1 (simulation, running next) — is the window narrowing real, or is it the censoring?

**Run.** `scripts/launch_capfix.sh`: both families, `μ_p` = 0 at CV 0.0/0.3/0.4/0.5 and
`μ_p` = 0.3 at CV 0.4/0.5, seeds 1–3, all four `p̂` rungs, at `PJ_MAX_ITER` = 1 000 000.
Tag `capfix`; the cell key carries `_it1000000` so these records cannot merge with the
capped ones they replace, and `max_iter` / `tol_met` are written into every record.

**The control, and it is not padding.** CV = 0.0 converged at ~11 000 iterations, three
orders below either cap, so raising the cap **cannot** change it. Pre-registered: `φ`
returns **bit-identical** on all 3 seeds × 4 rungs. If it does not, the protocol is not
deterministic under this knob and nothing else in the run is interpretable.

**Predictions.**
- Convergence at lognormal `μ_p` = 0, CV = 0.5 rises from **0.33 to > 0.80**. If it does
  not, the cells are not merely slow and the diagnosis in §3 is wrong.
- **H-1a, the narrowing is censoring:** `|dW/dCV| < 0.006`, i.e. no longer resolved at
  2 se.
- **H-1b, the narrowing is real:** `dW/dCV ≤ −0.012` and still `t < −2`.
- **Third branch, declared now:** a value in (−0.012, −0.006) means it is partly real;
  the repaired number is quoted and the −0.0169 is retracted as an over-estimate. This is
  the outcome the §3(b) arithmetic (28 % of the effect) actually predicts, and it will be
  reported as a *reduced* effect, not as either hypothesis passing.

**Separation.** 0.0169 against 0, on a slope whose se is 0.0037 → **4.6**, above the 1.6
informativeness floor.

**What it cannot decide.** Nothing in §1's traverse. `dφ_μ/dCV` and `ΔCV` = 0.50 live in
the frictional arm and are not at stake.

### H-2 (bench, a rheometer the lab has) — does the flow curve change SHAPE with CV?

The discriminator M4 owns and M2a does not. At fixed `φ` and fixed mean size, measure
steady flow curves `η(σ)` on a CV ladder, and score the **shape**, not the magnitude:

- **M4:** the curve is S-shaped / has a region of `d log η / d log σ ≥ 1` at low CV, and
  becomes monotonic and continuous at high CV. The change is a **threshold in CV**.
- **M2a:** every curve keeps the same shape; only the prefactor falls, smoothly and with
  no special CV.

⚠ **Score the shape by a pre-registered statistic, not by eye** — the local logarithmic
slope's maximum, with its own bootstrap over the measured points — and require the two
candidates' predictions to be at least 1.6× apart on the ladder actually run before
calling the test informative. ⚠ And note the standing literature caveat: Guy *et al.*
report Wyart–Cates *failing* for suspensions with a predominant fraction of large
particles, i.e. the broad-lognormal end of this ladder
(`LITERATURE_2026-08-31.md` §1). A shape change there is ambiguous between "left the
window" and "Wyart–Cates does not apply", so the ladder must be dense enough to locate
the change *before* CV reaches that end.

**This experiment is worth more than E-2** (raise `φ` until the solid branch returns),
because it is on a bench instrument rather than in a drop tower and it separates M4 from
M2a directly, where E-2 is consumed by both.

### H-3 (from data already on disk) — is `φ_μ` the same number under shear and compression?

The window's lower edge is measured here by *isotropic compression*. The concurrently
running common-ladder campaign (`PREREG_COMMON_LADDER_2026-09-09.md`) measures the
frictional critical state by *constant-volume shear* on a `φ` grid that is identical at
every CV. Those are different protocols and different code paths.

- **Prediction:** `φ_c(common ladder, CV)` − `φ_J(iso, μ_p = 0.3, CV)` is a **constant**
  offset in CV, to within 0.010. A constant offset means the two protocols measure the
  same CV dependence and the window's lower edge is protocol-independent up to a shift.
- **If the offset moves with CV**, the window has no single lower edge and `x` is not
  well defined; M4 would then need the *shear* pair, which does not exist at `μ_p` = 0
  above CV = 0.

⚠ **This comparison may NOT be made on the centred production ladders.** Those are
centred on the very `φ_J(iso, μ_p = 0.3)` values it would be comparing against, so on
them the agreement is forced. It is exactly why the common ladder was run.

### H-4 (bench, one nozzle) — is there a boundary in `(U₀, D₀)` at all?

Unchanged from `MECHANISM` §7a / E-6 and restated here because M4 shares M2a's silence on
it: **a viscosity has no `(U₀, D₀)` boundary and neither does a window position.** Only
M2′ does, through the drainage number. So H-4 is the test that separates the *pore-fluid*
family {M2′} from the *packing* family {M2a, M4}, while H-2 separates M4 from M2a inside
the packing family. Run on cornstarch, score the regime, Newtonian control at matched
viscosity and density.

**The three experiments are not alternatives — they are a 2-bit decision tree:**

```
                    H-4: boundary in (U0, D0)?
                     yes /                \ no
                 M2' pore-fluid       packing family
                                            |
                             H-2: does the flow-curve SHAPE change with CV?
                              yes /                        \ no
                              M4 window traverse            M2a plain rheology
```

---

## 6. What would kill M4

| observation | verdict |
|---|---|
| `W` measured under shear at both frictions is not ≈ 0.06 | the window identification is wrong; M4 has no stage |
| the flow-curve shape does not change over a CV ladder that spans `ΔCV` = 0.5 | dead — the traverse predicts a type change over exactly that |
| the shape change happens but at a CV unrelated to `(φ − φ_μ)/W` | mechanism survives, the quantitative version does not |
| `φ_μ` under shear moves with CV differently from `φ_μ` under compression (H-3 fails) | `x` is not well defined; M4 must be restated on shear quantities only |
| the lab's `φ` is confirmed at 0.44 with no densification at impact | M4 is irrelevant *to the lab's own experiment*, though it may still be right about the material class |

---

## 7. OUTCOME of H-1 — graded 2026-09-09, `out/capfix_grade.txt`

**180 cells, 0 nonzero exits**, `PJ_MAX_ITER` = 1e6, `grade_capfix.py`.

### Convergence — **PASS**, and every level now resolves three seeds

| family | `μ_p` | CV | old cap | new cap |
|---|---|---|---|---|
| lognormal | 0.0 | 0.3 | 0.58 | **1.00** |
| lognormal | 0.0 | 0.4 | 0.42 | **0.92** |
| lognormal | 0.0 | 0.5 | **0.33** | **0.83** ← the bar (0.80) |
| uniform | 0.0 | 0.5 | 0.17 | **0.75** |
| lognormal | 0.3 | 0.4 / 0.5 | 0.75 | **1.00** |
| uniform | 0.3 | 0.5 | 0.50 | **1.00** |

The censoring was slowness, as diagnosed. `φ_J`(μ_p = 0, CV 0.5) is now
**0.690401 ± 0.004665 on 3 seeds**, against 0.683324 on **one** seed before — so
`Δφ_J`(μ_p = 0, CV 0 → 0.5) = **+0.0580**.

### The control — **FAILED, and it was MIS-SPECIFIED, which is the finding**

It assumed every CV = 0 cell stopped "at ~11k iters, three orders below either cap".
That is the median, not the set. Split three ways:

- **(i) 2 of 12 sat at EXACTLY the old cap.** Raising it is *supposed* to move those, and
  it did: both went `converged` **False → True** at 205 276 and 211 806 iterations — 3–6 %
  more than the old cap. Lumping them into a bit-identity control makes a working knob
  look broken.
- **(ii) The real control, 10 cells that stopped early: 5 bit-identical, 5 not**, worst
  `|Δφ|` = **5.028e-03**.
- **(iii) `φ` is identical EXACTLY WHEN the iteration count is identical, on all 10.**

🔑 **So the trajectory is deterministic and the STOPPING is not** — the convergence check
fires at a different iteration between runs. That is a stopping-rule artefact, it
**pre-dates this knob**, and the worst value here (5.028e-03) is the *same cell and the
same size* as the worst pre-existing duplicate conflict from §3(c). This run is a third
independent draw confirming it.

⚠ **The knob is NOT exonerated**, because the control cannot separate the two effects at
that cell. What can be said: **the knob's effect, if any, is at or below the pre-existing
run-to-run spread of 5.0e-03 at CV = 0.** "Void" would be the wrong word and the grader
says so.

### `dW/dCV` at the production friction — **THIRD BRANCH, as the arithmetic predicted**

```
   dW/dCV = -0.01168 +- 0.00351 (fit)  +- 0.00229 (the CV = 0 systematic, propagated)
   t = -3.33 on 6 levels, all now 3 seeds     pre-repair value: -0.01690
```

Between H-1a (< 0.006) and H-1b (≤ −0.012), so the **declared third branch**: partly
real. **`−0.0169` is retracted as an over-estimate; quote `−0.0117`.** §3(b) predicted
censoring would account for "≈28 %" of the narrowing; measured, it is
(0.0169 − 0.0117)/0.0169 = **31 %**.

The systematic was *computed*, not asserted: moving `W`(CV = 0) across its 0.0032
run-to-run spread takes the slope to [−0.01397, −0.00939], i.e. **±0.0023 against a fit
error of ±0.0035** — smaller than the fit error, so the slope survives at t ≈ −2.8
combined. (An earlier back-of-envelope guess of ±0.009 was wrong; this is why it is a
line of code.)

### 🔑 At the SATURATED friction — the window is FLAT, which is what M4 needed

`μ_p` = 1 is Wyart–Cates' `φ_μ` as this protocol can reach it (`φ_J` stops moving between
`μ_p` 1 and 3). Measured across CV for the **first time**:

| CV | `φ_0` | `φ_μ` | `W` |
|---|---|---|---|
| 0.0 | 0.632371 | 0.559036 | 0.073334 |
| 0.2 | 0.645392 | 0.572635 | 0.072757 |
| 0.5 | 0.690401 | 0.615247 | 0.075154 |

```
   dW/dCV = +0.00398 +- 0.00298 (fit) +- 0.00295 (systematic), t = +1.34  -- NOT RESOLVED
```

**So the narrowing is a `μ_p` = 0.3 phenomenon, and at the friction the Wyart–Cates
`φ_μ` actually means, the window is a rigid object.** That is M4's premise, and it is
now measured rather than assumed.

### The repaired traverse — M4's headline barely moves

| friction | `dφ_μ/dCV` | `W`(CV 0) | windows per unit CV | traverse ΔCV |
|---|---|---|---|---|
| `μ_p` = 0.3 | +0.1284 ± 0.0127 | 0.0575 | 2.23 | **0.448** |
| `μ_p` = 1.0 (saturated) | +0.1148 ± 0.0203 | 0.0733 | 1.56 | **0.639** |

**Quote ΔCV = 0.45 to 0.64** (was 0.50–0.65 pre-repair). §1's headline survives the
repair; only its third significant figure moved.

### The zero-skew control also sharpens

`uniform`, `μ_p` = 0.3: `dW/dCV` = **+0.00036 ± 0.00176**, t = +0.21 — bounded to **3 %
of `W` over CV 0 → 0.5**, against 12 % before the repair. The skewed/zero-skew contrast
in §1 therefore stands, and stands on a tighter bound.

### One decision left to the author, and it is deliberately not taken here

`Z_ISO_MEASURED` and `PHI_J_ISO` were measured under the **default** cap, and every
`φ_c` already reported rests on them. The repaired isotropic set is better, and adopting
it would move E4m everywhere — **including the common-ladder verdict graded the same
day**. The two tests that pin those constants now read the default-cap rows explicitly
and say why; `test_a_REPAIRED_phi_J_set_exists_and_differs` fails the day the repaired
rows are deleted *or* the day the two sets stop differing, so the decision cannot be
lost.
