# The measured clock is now the default — and it costs two headline results

**2026-08-09, second session.** Read §1 and §5.

---

## 1. The short version

`WC_CLOCK=1` is the default. Then the Stage 3–5 headline cells were re-run under both
laws, in single campaigns so the comparison is internal, and **two of the project's three
headline claims do not survive.**

| claim | legacy (guessed clock) | measured clock |
|---|---|---|
| **Stage 4** two-clocks retention: silica retains more | **8/8**, median 2.43× | **0/8**, median 0.171× — **inverted** |
| **Stage 5a** drive discrimination: b-ratio > 2 | 2.277, PASS | 1.542, **5.1σ below threshold, FAIL** |
| **Stage 3** convergence structure | — | unchanged in kind; cost is 3–10× |

Neither was wrong arithmetic. Both were correct consequences of a parameter nobody had
measured. The control arm reproduces the original result in each case, which is what makes
the comparison mean anything.

Separately, the landscape research came back and **three things this project believes about
its own novelty are false.** §4.

---

## 2. Stage 4: the two-clocks result inverts

32 cells, matched (φ, We), φ ∈ {0.44, 0.47} × U₀ ∈ {1.5, 2, 3, 4}, both clocks, read at a
common t/τ = 5.5. Zero cells excluded. The solver's own `CLOCK` banner was checked against
what the launcher requested in all 32.

| | si retains more | median si/cs | range |
|---|---|---|---|
| WC_CLOCK=0 (legacy) | **8/8** | 2.43 | 1.92 – 3.32 |
| WC_CLOCK=1 (measured) | **0/8** | 0.171 | 0.041 – 0.409 |

**The ranges do not overlap.** The ordering does not weaken, it reverses, in every pair, at
both volume fractions and all four speeds.

**The control passes**, which is the precondition for reading any of this: the legacy arm
returns Stage 4's own 8/8 at median 2.43 against the recorded 3.0. Count and direction
reproduce exactly; the magnitude sits ~20% low, which is the cross-campaign magnitude
instability this programme has documented repeatedly.

**Why it inverts, and it was pre-registered before the run.** The result is *named* the
two-clocks asymmetry because it was built on cs and si being handed `tau_relax` of 0.02 s
and 2.0 s — a factor of 100 typed in by hand. Under the measured law both materials share
γ₀ and t_r, so retention is set by `fss = exp(-(sigma*/sigma)^delta)` and **the material
with the lower onset stress recruits contacts more readily.** The calibrated sets put
cornstarch at 10.52 Pa against silica's 44.23 — a factor 4.2 the other way, which the 100×
clock difference had buried.

> The surviving physics is a **σ\*/φ_M asymmetry, not a clock asymmetry**, and it points the
> opposite way. The claim has to be rewritten, not renamed.

⚠ **Caveat that must travel with it.** Under the measured clock the absolute retentions are
tiny — 1e-4 to 3e-2, against 0.11–0.37 under the legacy law — because the strain term lets
`fhat` decay as soon as shearing stops. A ratio between two small numbers is less robust
than one between two large ones, and the smallest cell (si = 9.4e-5) may be a handful of
cells. **The 8/8 ordering is the claim; the median 0.171 needs a resolution check before it
is quoted as a magnitude.**

---

## 3. Stage 5a: the drive discrimination collapses

20 cells, σ\* ∈ {11.06, 22.12, 44.23, 88.47, 176.93} Pa × two drives × two clocks, at
U₀ = 2.0 — the only column Stage 3 found converged to ~3%. Zero excluded.

| | b (stress) | b (inertial) | ratio | verdict |
|---|---|---|---|---|
| legacy | 0.718 ± 0.090 | 0.316 ± 0.019 | **2.277** | PASS (0.87σ above 2) |
| measured | 3.553 ± 0.155 | 2.303 ± 0.089 | **1.542** | **FAIL (5.09σ below 2)** |

Control reproduces: 0.718 / 0.316 / 2.277 here against Stage 5a's recorded 0.623 / 0.278 /
2.459, within the window-dependence that recompute already documented (the ratio was
recorded as moving 1.99–2.79 depending which rungs are fitted).

Both exponents rise 5–7× under the measured clock, which is expected — `<fhat>_final` is a
*retained* quantity, so it reads the clock directly, and a state variable that can actually
decay responds far more steeply to the onset stress than one held frozen. What does not
survive is the **separation**, and the separation was the gate.

**This is not a noisier fit washing out a real effect.** The power-law assumption is
rejected in all four arms, but the local log-log slope spans **1.52× and 1.62×** under the
measured clock against **3.32×** for the legacy stress arm. The exponents are more
power-law-like than they have ever been, and the ratio still lands 5σ the wrong side.

---

## 4. The landscape research: three novelty claims are false

Full document: `SOLVER_LANDSCAPE_2026-08-09.md`.

**What cannot be claimed any more:**

- ❌ *"First to treat the frictional contact fraction as a dynamical variable."* Baumgarten
  & Kamrin (PNAS 2019) already do — their *f* obeys `ḟ/(K₀γ̇) = H(f_m−f) − Sf` and sets both
  φ_m(f) and the dilation parameter.
- ❌ *"First advected Wyart–Cates field."* Angerman, Sandnes, Seto & Ellero (Phys. Fluids
  2024; JFM 1016 A33, 2025) solve `Df/Dt = k_f γ̇(f̂−f) + α∇²f` in SPH — advected, WC steady
  state, **plus a non-locality we do not have.**
- ❌ *"Nobody has modelled shear jamming under impact."* Athani, Metzger, Forterre & Mari
  (PRF 10, 043301, 2025) built the WC + Reynolds-dilatancy closure and state its purpose is
  inhomogeneous flows "where shear jamming is local, e.g. under impact".

**What is genuinely unoccupied is the intersection.** Each competitor holds at most two of
{WC field with its own dynamics, sharp capillary interface, free-surface impact}. Mobaseri,
Kumar & Cheng (PNAS 122 e2500163122, 2025) run cornstarch drop impact on Basilisk VOF and
state in print that the DST/shear-jamming transition "is not captured by the current study".

**And the sharpest technical differentiator is the thing built last night.** Both prior
WC implementations relax *f* on a **pure strain clock** — so neither can un-jam material
that has stopped moving. The `1/t_r` term is the differentiator, and it only exists as of
yesterday.

**Two referee exposures, and the first one connects to the pore-size question:**

1. **No Reynolds dilatancy / Darcy backflow.** This is Athani's central result and should be
   expected as referee objection #1. It is also exactly where last night's pore-size
   adjudication landed — *drainage sets the stress, not the clock* — so the user's original
   instinct and the competitive analysis point at the same missing term, independently.
2. **No fracture.** Slutzky, Pelosse, van der Naald & Jaeger (PRE 2026) show the jammed
   state cracks, with likelihood rising **as surface tension is reduced** — which makes it
   both an exposure and the best available validation target for a solver whose whole point
   is that it has capillarity.

**Commercial CFD, honestly:** nothing off the shelf does this — Fluent, STAR-CCM+ and
OpenFOAM offer only algebraic η(γ̇), and COMSOL's Houska thixotropy has the wrong sign and no
φ_J. **But do not claim commercial CFD has never touched it:** Veryst Engineering publishes
a COMSOL case study tracking "the fraction of frictional contacts" for nozzle blockage and
injection clogging. That is a channel partner and a list of paying customers, not a
competitor to dismiss.

**Best-fit applications** (fit × pain): direct ink writing of dense pastes (needs wall
slip); high-solids coating and battery slurries (needs 3-D — ribbing is not axisymmetric);
STF impact protection (best physics match, but fielded armour is confined, so the free
surface is not the governing feature — say so first); small-bore clogging. Poor fits worth
naming as poor: proppant screen-out and STF dampers, both confined.

**Funding:** NSF **CBET Particulate and Multiphase Processes PD 23-1415** is the closest
match and names suspensions explicitly. Cheapest shots on goal are rolling DoD BAAs —
DEVCOM Soldier Center W911QY-25-R-0023 and ARO W911NF-23-S-0001.

---

## 5. What I recommend, and what needs you

1. **Rewrite Stage 4's claim rather than renaming it.** The retention ordering is real, it
   is a σ\*/φ_M effect, and it points from cornstarch to silica — the opposite of what is
   currently written. This is the biggest single writing job.
2. **Stage 5a's gate should be retired, not re-tuned.** It was already weak three ways
   before today (1.09σ, window-dependent, both fits rejected as power laws). It now fails by
   5σ. Reporting the two exponents and their local slopes is more honest than a threshold.
3. **Implement Reynolds dilatancy + Darcy.** It is simultaneously referee objection #1, the
   correct home for the pore-size physics, and the thing that would let the solver claim the
   Athani closure in a free-surface flow — which nobody has done.
4. **The `1/t_r` term is the paper's technical hook.** Nobody else's WC field can un-jam at
   rest. That should be stated in the first paragraph, and it is a stronger and more
   defensible claim than any of the three novelty claims that just died.
5. ⚠ **Still open from last night:** the `jam_frac` factor-0.44 mismatch against the 29 JUL
   campaign, where flow reproduces and the state variable does not. Nothing here resolves it,
   and it is the reason no `jam_frac` number should be quoted across campaigns.

---

## 6. Other projects

| project | state |
|---|---|
| **frisbee** | ✅ **Done.** Best-direction claim withdrawn in text, figure and JSON. 26 pp, 0 undefined refs, 0 `\MISSINGNUM`. The paired test resolves only that 270° is worst; 90° vs 210° is not separable (p = 0.206, CI spans zero, favoured at 4/12 launch points). Needs your author line. |
| **catfish** | ✅ **Done.** Fringe numbers were 1.566× too small — a silent fallback to `working_body_px` (21.70) instead of the measured 13.86. Now raises instead of substituting. **One conclusion reversed:** fringe NND is 1.33–1.39 BL, *over* one body length, where the docs said 0.85–0.95 and concluded "under one". The section's actual conclusion survives untouched, because it is a ratio of two speeds sharing the length unit. |
| **beach_sand** | ✅ **Done.** All five corrections verified from code and applied, 22/22 tests pass. Speed range corrected to **0.20–2.5 m/s** (derived from the 0.1 s footstep, not guessed). Found an **extra** defect: §2's footstep ratio is 95, not "a factor of a few" — which independently supports promoting Q4. |
| **datacenter** | ✅ **Done and independently verified.** I recomputed the sweep myself: cfd 1.2276 → 1.2049, linear ΔPUE = 0.0396, matching the README's "1.227→1.205" and "0.039" exactly. |
| **bidet** | ✅ Done (verified last night). |
| **bloodstain** | 🔴 **Do not write.** Core idea published twice — Laan et al. Sci. Rep. 5:11461 (2015) and Attinger et al. FSI 298:323 (2019). |

⚠ **Two `LIMBO_AUDIT_2026-08-08.md` entries did not hold up on inspection** — the bidet
claim description and the catfish fig6 mismatch. Treat that audit as leads, not findings.
