# Overnight 2026-08-08/09 — what changed

Read §1 and §7. Everything else is detail.

---

## 1. The short version

**The solver's weakest parameter is gone.** `tau_relax` and `tau_up` — two guessed
numbers that gated whether the drop was jammed when it mattered — are replaced by a
measured law with two independently-measured constants. Swapping them moves peak
`jam_frac` from 0.12 to ~0.98. The old law was holding the drop's state variable
**frozen through the entire impact**.

**Three solvers were running at the wrong CFL, and the repo's list of which ones was
wrong in both directions.** All five files recorded as still broken had been fixed;
three that were never on the list had not. One of them is the instrument that certifies
the penalization kernel used elsewhere.

**Two long-standing nulls turned out to be detector artifacts**, both now measured:
the theta = 20/45 jamming dip, and Stage 5b's "the deposit never respreads".

**One thing you should know before writing anything on bloodstains:** the core idea is
already published, twice, by the people named in your own TODO as the prior-art check
that was never run. Details in §6.

**I hit the session limit at 12:40am** and five subagents were killed mid-task. §7 lists
exactly what is unfinished.

---

## 2. Solver: what is genuinely better now

### 2.1 The un-jamming clock is measured, not guessed  (commit `9fbad8a`)

The old law relaxed `fhat` on `tau_up` while rising and `tau_relax` while falling, both
PLACEHOLDER in every preset — necessarily so, since steady flow curves cannot fix a
transient. The asymmetry between loading and unloading existed **because two different
numbers were typed in**.

`WC_CLOCK=1` uses the measured law (Richards et al., PRL 123, 038004):

```
dfhat/dt = -( gammadot/gamma0 + 1/t_r ) [ fhat - fss(sigma) ]
```

One rate, both directions. The asymmetry now comes *out*: while flowing, the strain term
dominates; at rest it vanishes. `gamma0 ≈ 0.2` is a strain — Han et al. (PRF 3, 073301)
measure the same quantity independently at **0.197 ± 0.002** — and `t_r = 0.24 ± 0.05 s`,
with four flow-cessation studies bracketing it at 0.1–2 s.

Four gates pass, including an **exact-limit gate**: send `gamma0 → ∞` and the new law
provably reduces to the old one at a matched single time constant, and reproduces it
byte for byte. That gate tests the new branch rather than merely observing that it
differs — it would catch a units error, a sign error, or a mis-set rate.

**Why it matters quantitatively.** Impact shears the drop at `gammadot ~ U0/D0 ~ 670 s⁻¹`,
so the effective clock is ~0.3 ms against the 2 s the silica preset assumed — a factor of
~7000. At L8, peak `jam_frac` goes 0.116 → 1.000.

**The resolution check on that number** (`wcclock_ladder.sh`, 8 cells) came back:

| | L8 | L9 | L10 |
|---|---|---|---|
| strain clock, `jam_frac` at t/τ = 2 | 0.913 | 0.997 | 0.975 |
| legacy | 0.096 | 0.325 | *(did not reach t/τ = 2)* |

and the `gamma0` row at L9: **0.912 / 0.997 / 0.671** at `gamma0` = 0.05 / 0.2 / 0.8.

Read together, these say the result is **not** cap-limited: the pre-registered failure
mode was "flat in `gamma0`", and it is not flat — 0.671 against 0.997 is a large,
correctly-signed response. What *is* true is that at the measured `gamma0` the observable
is saturated near its ceiling, so **quote the ordering, not the value**. And note the
legacy baseline is itself unconverged (0.096 → 0.325 between two levels), so the *ratio*
between the two laws is not a quotable number either.

⚠ **Practical cost, and it is large.** The strain clock drives the viscosity to its
ceiling, which stiffens the implicit solve badly: the L10 cell took 3.8 h to reach
t/τ = 2.55 and then hit dt-collapse, against 23 min for legacy. Three of four
strain-clock cells ended in dt-collapse. This is now the binding numerical problem, and
it points straight at the augmented-Lagrangian item in §8.

### 2.2 A deposit observable that is not quantized  (commit `6f8fd22`)

`beta`, `beta_wet` and `h/D0` are all extremal (max y, max x over `f > 0.5`), so they only
move when the interface crosses a whole cell. Stage 5b measured deposit motion of 2–3
cells on those columns and concluded nothing moves — which is the smallest thing those
columns can express, not a measurement.

`SHAPE_MODE=1` adds five f-weighted integrals: `vol, zc, Rg_r, Rg_z, AR`. The useful one
is **`AR = Rg_r/Rg_z`, which equals √2 exactly for a sphere of any radius** — a scale-free
analytic anchor that does not depend on drop size or code units. Measured at t = 0:
1.41336 against 1.41421, a relative error of 6.0e-4. `AR > √2` reads as an oblate
(spread) deposit, `< √2` as a beaded one.

### 2.3 Three solvers were running at the wrong CFL  (commit `ad01d1b`)

`CFL = CFLnum;` in `main()` is a no-op — `centered.h` overwrites it inside `run()`. The
repo found this in July and recorded five files as still carrying it. **All five had been
fixed. Three others had not:**

| solver | asked | actually ran at |
|---|---|---|
| `ve_drop_on_bath.c` | 0.10 | **0.5** |
| `sphere_impact_hb.c` | 0.2 | **0.5** |
| `sphere_settle_hb.c` | 0.4 | **0.8** — no `vof.h`, so nothing clamped it |

Fixed and gated per-solver against *its own* pre-fix value. Two consequences worth your
attention:

- **The viscoelastic arm ran log-conformation at 5× its requested CFL** — the whole Γ
  ladder, the De rows, the τ_p peak-and-decay. Orderings are safe (every cell had the
  same wrong number); **magnitudes are not yet at the timestep their launchers claim.**
- **`sphere_settle_hb.c` is the instrument that licenses the Brinkman kernel** via the
  Stokes drag gate (1.029 / 1.013 / 0.992), and that gate was therefore taken at twice
  its intended timestep on a quantity already known to be dt-sensitive. **I checked it
  rather than leaving it as a worry: all runs take an identical 31 steps at both 0.4 and
  0.8, so their dt is set by DTMAX and not by CFL, and the drag numbers survive.**

---

## 3. Results from compute

### 3.1 The theta = 20/45 loose end is closed: it was an L9 artifact

Eight cells, four angles, two levels, with theta = 5 and 90 as reproduction anchors.

| theta | L9 peak `jam_frac` (at t/τ) | L10 peak `jam_frac` (at t/τ) |
|---|---|---|
| 5 | 0.3591 @ 1.21 | 0.4092 @ 1.18 |
| 20 | 0.0204 @ **9.46** | 0.3867 @ 1.33 |
| 45 | 0.0184 @ **9.25** | 0.3725 @ 1.13 |
| 90 | 0.2070 @ 1.33 | 0.4043 @ 1.18 |

At L10 all four agree to 9.3% of the mean, every peak lands at t/τ ≈ 1.2, and `beta_max`
is identical at three angles and one quantum away at the fourth.

**The give-away is the *time* of the L9 peak.** At L9 the intermediate angles peak at
t/τ ≈ 9.3–9.5 with a value of 0.02 — they never develop an impact-phase jamming peak at
all, and `argmax` is landing on late-time noise across a flat, near-zero series. The
"non-monotonic peak jamming" was two cells that failed to jam. **The theta-insensitivity
claim now covers 5–90° continuously**, which removes the caveat sitting under Stages 3–5a.

⚠ **One genuine discrepancy, not swept up.** `beta_max` reproduces the 29 JUL L10
reference to one quantum (2.19531 vs 2.1797) but peak `jam_frac` does not: 0.409/0.404
here against a recorded 0.9135–0.9331, a factor of 0.44. **The flow reproduces and the
state variable does not**, which is the more worrying way round. The conclusion above is
unaffected — it rests only on comparisons within this campaign — but do not quote these
absolute `jam_frac` values beside the 29 JUL ones until the configuration difference is
found. The analyzer records this as `A2_jamming_anchor: MISMATCH` rather than letting a
beta-only anchor pass the campaign.

### 3.2 The un-jamming clock DOES have a deposit signature — and both fixes mattered

14 cells: `tau_relax` over four decades × theta ∈ {5°, 120°}, read on `AR` at a common
t/τ = 6.0. Two cells excluded for not reaching the read time (never read at their own
last frame — that rule is what the retracted pool `cavity_max` was about).

**Control passes exactly.** The Newtonian arm (`sigma* = 1e9`, so `fhat ≡ 0`) gives
`AR = 5.28087` at both ends of the ladder at theta = 5, and `2.10723` at both ends at
theta = 120 — identical to nine decimals. There is no state variable to relax, and
nothing responds. The campaign is valid.

| `tau_relax` (s) | AR at theta = 5 | AR at theta = 120 |
|---|---|---|
| 2e-4 | *(excluded)* | 2.78665 |
| 2e-3 | 4.48045 | 4.09335 |
| 2e-2 | 4.66911 | 4.56281 |
| 2e-1 | 4.90796 | *(excluded)* |
| 2e0 | 4.95713 | 4.80078 |
| **span** | **10.0%** | **49.6%** |

**Both ladders are monotone**, and the response at the experimental angle is **4.95×**
larger.

**This corrects my own pre-registration, which is the useful part.** I predicted theta = 5
would reproduce Stage 5b's null (span < 2%). It does not — it shows a clean monotone 10%
response. So the honest reading is that **Stage 5b's null was mostly the detector, not the
contact angle**: the continuous observable finds a monotone signal at *both* angles where
the extremal ones found none, and the correct contact angle then multiplies it by five.
Both fixes were needed; the observable was the bigger one.

The physics is the predicted mechanism: `AR` *rises* with `tau_relax`. A deposit that
un-jams quickly is free to retract under its own capillarity (AR falls toward the sphere
value √2); one that never un-jams stays a pancake.

---

## 4. The pore-size question you raised

Two independent investigations were run without sight of each other; they disagreed, and
the disagreement is the result. Full adjudication:
`Desktop\Projects\tier0_wyart_cates_closure\PORE_UNJAMMING_ADJUDICATION_2026-08-08.md`.

**The derivation.** To un-jam, the network must dilate; the drop is incompressible and its
free surface capillary-sealed, so the new pore volume must be fed by Darcy flow through
the packing. That gives `tau ∝ eta_f L² / (k sigma*)` with `k ∝ a²`, and it predicts
silica 1.5 s / cornstarch 0.013 s, **ratio 111 against the solver's assumed 100** — with no
free parameter in the ratio, robust across three permeability models, and a one-parameter
cross-prediction good to 10% across a factor of 100. It also shows that *every*
contact-scale mechanism gets the sign wrong, because they all reduce to `eta_f/sigma*`
times something without `a` in it.

**The literature says it is wrong, and identifies why it looked right.** Release times have
been measured repeatedly (Maharjan & Brown 2017; Cho 2022; Richards 2019; Barik &
Majumdar 2022) and cluster at **0.1–2 s across a 5× range of particle diameter**.
Cornstarch measures at 0.2–2 s — the solver's 0.02 s placeholder is ~10× too fast. **The
derivation agreed with a guess.** Once the target moves to the measured value, silica and
cornstarch are not two decades apart, and a mechanism whose job was to produce two decades
is solving a problem that does not exist.

> **The resolution: pore-scale drainage sets the STRESS, not the CLOCK.** Mobaseri,
> Driscoll & Cheng build exactly this Kozeny–Carman permeability into a drop-impact model
> — so the picture is current practice in this geometry — but their dilatancy time is
> **~9 µs**, five decades below the measured release times. Athani et al. reach the same
> structure: Darcy bounds the stress while the state variable relaxes on a strain clock.

Your steer was pointing at real physics on the wrong observable. Keep the mechanism, move
it to the stress.

**Honest caveat:** the measured times *do* span ~20×, comparable to the 25× that `a⁻²`
predicts over that size range. What is missing is any *correlation* with size, and **no
single study has varied particle size systematically**. So this is strong enough to stop
us building `a⁻²` into the solver and not strong enough to call it refuted. §8 has the
experiment that would close it.

**On electrostatics — half right.** For the *solver*: right, and for a better reason than
"small" — whatever electrostatics does, it does through `sigma*`, a single calibrated
scalar, so there is no separate double-layer physics to model. For the *material*: not
small, and specifically not for your own system. Ishibashi et al. (2026) fit silica
110–2580 nm with a double-layer term carrying `F* = 0.3 nN` at 0.1 mM and **`F* = 0` at
100 mM**; it is the whole story at low salt and gone above the critical coagulation
concentration. Treating it as minor for **low-salt colloidal silica would discard Paper
1's own independent variable.** Against that, Bourrianne et al. change *only* surface
chemistry on identical 300 nm silica and **extinguish DST entirely** — so chemistry and
roughness are larger where they act.

---

## 5. Complex substrates: a new solver, and an honest partial result

`wc_drop_on_texture.c` — a Wyart-Cates drop onto a **concentric-groove** substrate. Grooves
are the one texture that is axisymmetric by construction, so this needs no 3-D, and they
are a real bench geometry (a lathe cuts them from ~50 µm up). The question it exists to
ask: **texture sets a threshold by geometry and a shear-thickening drop sets one by
rheology — which wins?** Texture is the only substrate threshold you can sweep on a bench
without a rheometer.

Reduction gates **pass**: `TEXTURE=0` reproduces the parent solver on all nine columns, and
a zero-depth groove is inert.

The go/no-go did not pass, and the failure is informative:

- **Viscosity penalization (the spec's Route C as written) is dead.** Every cell dies with
  a floating-point exception at t = 0, identically at `MU_SOLID` = 1e2, 1e4 and 1e6.
  Insensitivity to the penalization magnitude is the signature of a structural failure,
  and this repo has met it before: the 2026-08-04 Stage 1 diagnosis found that a
  high-viscosity body suspended in **air** is a disconnected island of enormous viscosity
  contrast on which the first implicit `u.x` solve fails at t = 0. **Empty grooves are
  exactly that geometry** — contrast 6.7e8 even at the weakest setting, against the 6.7e7
  that already killed Stage 1.
- **Brinkman drag instead — the kernel `sphere_settle_hb.c` already validates — runs.**
  Cells reach t/τ ≈ 1.0–1.14 instead of dying at t = 0.
- **It leaked at first — `u_solid_max` ≈ 4× the impact speed — and the leak was FLAT
  across `ETA_PEN` = 1e-6 to 1e-10.** Ten thousand times the penalization strength with
  no change cannot be a weak penalization; at 1e-10 the implicit update divides the
  velocity by ~1e6 per step, so the velocity *was* being zeroed and something was putting
  it back.
- **The something is the projection**, which enforces only that the field is
  divergence-free — it knows nothing about the solid and is free to drive flow straight
  through it. Penalizing before the projection is the right place to shape what the
  pressure solve sees and the wrong place to have the last word on the velocity.
  **Doing both passes fixes it:**

| `ETA_PEN` | `u_solid_max` | rc |
|---|---|---|
| 1e-6 | 0.01099 | 0 |
| 1e-8 | 0.0001103 | 0 |

The leak falls by **99.6× for a 100× change in `ETA_PEN`** — linear, which is the
signature of a volume penalization working as intended rather than of a number that
merely happens to be small. At 1e-8 the substrate holds to 1e-4 of the impact speed, and
`impaled_frac` agrees to four figures between the two settings (0.3735 / 0.3736), so the
result no longer depends on the knob. **GO.**

**No impalement number is trustworthy yet** and none is quoted. Still owed, per the spec:
the Wagner check on the flat wall (the programme's one external analytic anchor, C → √6)
and a groove-resolution pair — 250 µm is only 5.3 cells at L9 and should be expected to
move. Also recorded in the solver header rather than discovered later: `contact.h` imposes
the contact angle on the domain boundary only, so the groove walls carry **no wetting
condition** — this model tests inertial impalement, not the capillary threshold
`P_c = 2σ|cosθ|/w`.

---

## 6. Bloodstain project: stop before writing

The `[CRITICAL — check before writing anything]` prior-art item in that project had never
been done. It is done now, and the answer is decisive.

**The core idea is published, twice, independently — by the authors named in the flag.**

- **Laan, de Bruin, Slenter, Wilhelm, Jermy, Bonn (2015)**, *Sci. Rep.* 5:11461. Verbatim
  from the abstract: investigators use a straight-line approximation "ignoring effects of
  gravity and drag and thus overestimating the height of the source… by including gravity
  and drag… **roughly four times more accurately** than with the straight-line
  approximation." Also inverts a drop-impact spread law to get velocity from stain
  geometry, and quantifies the stringing bias (up to 45 cm).
- **Attinger, Comiskey, Yarin, De Brabanter (2019)**, *FSI* 298:323–331. Backward ballistic
  integration with drag and gravity, automated stain segmentation, propagated measurement
  uncertainty, and **nested 90/99/99.9% likelihood region-of-origin volumes** — error
  < 10 cm at P = 99%.

Both were verified from source (abstracts plus the full NIJ final report), and a check of
the OSAC BPA bibliography through 2025 found no newer drag-corrected origin method.

**What may survive** is narrow and worth knowing: neither paper states the *conditioning*
of the stain-size→velocity inversion, and both sidestep it — Attinger by carrying a range
of drop-size/velocity pairs rather than one velocity, Laan by *measuring* D₀ with OCT
instead of assuming it. So a quantitative sensitivity result is a legitimate
**methodological note on an existing method**, not a new method. Two things would need
fixing first: the work is synthetic-only when public benchmark data exists (Attinger et
al., *Data in Brief* 18:648–654 and 22:269–278), and the README's novelty claim needs
those two papers read in full.

Separately, the audit found the README's "overestimates by 1.3–2.2×" is wrong — the
figure plots +127% to +217%, i.e. 2.3–3.2×, and the run output gives 2.44×. Percentages
were written down as factors. **Not yet fixed** (§7).

---

## 7. What is NOT done — I hit the session limit at 12:40am

Five subagents were killed mid-task. None of these is blocked on anything except time:

| item | state |
|---|---|
| **bidet** headline fix | ✅ **DONE and verified** — I checked this one by hand after the agent was killed, and it had finished before the kill. The τ_c = 1 Pa contour genuinely does lie outside the R = 25 mm domain (`case_Re10000_H6` ends at r = 24.77 mm with τ_w = **1.957 Pa**), and the manuscript now says so: it gives both extrapolations (39 mm from the full closure, 35 mm from an outer-half refit, r² = 0.9998), **labels them as extrapolations rather than computed values**, rescopes the "within a factor of ~2" claim to τ_c ≳ 3 Pa where it actually holds (ratios 1.4 / 1.1 / 0.8), states the τ_c = 1 Pa case as an extrapolated 2.3–2.6× over-provisioning bounded only by r > 25 mm, and says what would settle it (re-run at R ≈ 50 mm). Abstract, results, figure caption and design implication all updated. **Recompiles clean: 12 pp, 0 errors, 0 undefined references.** |
| **beach_sand** protocol | **Not started.** Five verified corrections in `DESIGN_FINDINGS_2026-08-07.md` never reached `PROTOCOL.md`, which still specifies 0.5–4 m/s (D/V ≤ 0.04 s against a ~0.1 s footstep) and still collides the symbol φ. Has a hardware cost if it slips. |
| **frisbee** correction | **Partially applied — check state before continuing.** The agent was mid-edit. `significance_generated.tex` is orphaned (never `\input`), `\QbestHeading` still asserted, and it had just found that `fig4_throwability` **prints "best"/"worst" into the pixels**, where no macro guard can see it. |
| **catfish** fringe | **Partially applied.** Every number in `fringe_analysis.json` is 1.57× too small (hardcoded 21.70 px fallback instead of the measured 13.86), and it is the only script skipping the camera-pan correction. |
| **datacenter** figures | **Not started.** Both figures predate the Jul 25 switch to the CFD recirculation model; the containment sweep is the one most affected. |
| **solver landscape / industrial pitch** | **Not started.** This was your "next-gen / industrial" ask and it is the one substantive item with nothing to show. |

**Recommendation for the next session:** **beach_sand first** — under an hour, it is a wrong
claim rather than a missing one, and it has a hardware cost if it slips. Then frisbee and
catfish, both of which are in a *partially edited* state and should be checked before being
continued.

⚠ **Treat `LIMBO_AUDIT_2026-08-08.md` as leads, not findings.** Every item I verified by
hand was directionally right, but the bidet entry described the manuscript's claim in terms
that do not match what is in the file, and the audit was itself written while agents were
editing the same folders. Verify before acting on any entry — which is what the audit
itself says to do about this repo's other summary documents.

---

## 8. Decisions and open questions for you

1. **Should `WC_CLOCK=1` become the default?** It is measured where the current default is
   guessed, but it changes every impact-phase number and makes runs 3–10× more expensive
   with frequent dt-collapse. My recommendation: keep the default at 0, re-run the Stage
   3–5 headline cells at 1, and report both — the difference *is* a result.
2. **The `jam_frac` mismatch against the 29 JUL campaign (factor 0.44) needs a cause.**
   Flow reproduces, state variable does not. Until it is found, no `jam_frac` number from
   this programme should be quoted across campaigns.
3. **The stiff solve is now the binding constraint**, not the physics. The capability
   audit ranks augmented Lagrangian as the highest-value single item (40–60 h) — it
   removes the regularization parameter entirely rather than tuning it.
4. **The one experiment that would settle the size scaling**: one material, three particle
   sizes, same chemistry, same φ/φ_M, stress-drop release time for each. The whole
   literature is missing it. **A null result is publishable on its own.**
5. Cheaper and immediate: **sweep `t_r` over the measured 0.1–2 s bracket** and report the
   sensitivity instead of one value. If impact-phase results are flat across it — likely,
   since the strain term dominates during impact — the clock stops being a threat to any
   conclusion.

---

## 9. Commits

`6f8fd22` shape observable · `36a8d61` two campaigns · `ad01d1b` CFL fix in three solvers
+ texture solver · `9fbad8a` strain-clock law · `1fd3039` theta loose end closed

Also written: `SOLVER_CAPABILITY_AUDIT_2026-08-08.md` (all 15 solvers, verified against
source — it is what found the CFL and hoop defects), `NEW_PAPERS_SYNTHESIS_2026-08-08.md`,
`UNJAMMING_TIMESCALE_LITERATURE_2026-08-08.md`, `LIMBO_AUDIT_2026-08-08.md`,
`PORE_UNJAMMING_ADJUDICATION_2026-08-08.md`.

⚠ **One more defect found and NOT yet fixed**, from the capability audit: the axisymmetric
hoop term in the face strain-rate invariant is written inside a generic `foreach_face()`,
which rotates field components but not coordinates — so on radial faces it computes
`u_z/r` instead of `u_r/r`. The spurious term scales as **1/Δ**, i.e. it *grows* under
refinement (~64 at L9 against a true γ̇ of 1–10). Live and unconditional in
`rheo_substrate.c`, `sphere_impact_hb.c` and `sphere_settle_hb.c`. `hoop_gate.sh` cannot
catch it — its sign test passes for the wrong reason.
