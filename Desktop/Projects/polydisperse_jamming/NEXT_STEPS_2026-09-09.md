# NEXT STEPS — polydispersity, particle size, and the impact regime

**Written 2026-09-09, end of an autonomous session.** Feed this file to a new session.
It is the single instruction sheet; everything else it names is a supporting document.

**Read first, in this order (about 10 minutes):**

1. `SYNTHESIS_2026-09-08.md` — what is established, what is not, one page.
2. This file, §0 and §1.
3. `Documents\Obsidian Vault\Status\Drop Impact — Polydispersity Jamming Model.md` §1–3.

Do **not** read `REPORT_2026-08-31.md` to get oriented; it predates all of this.

---

## 0. 🔴 START HERE — one defect is open and it moves the headline by 30 %

A numerics reviewer found this at the very end of the session and it is **not fixed**,
only measured and documented.

**`analyse_phase8.py` counts its quality gates but never applies them.** `e1_phi_of_I`
filters on the quasi-static bar; `_pressure_root` (E2 — **the primary estimator**) and
`e4_z_crossing` filter on nothing but a pressure floor. The module docstring claimed
cells were excluded. They were not.

Measured on the production run:

| | Δ`φ_c` (CV 0 → 0.5) |
|---|---|
| as shipped, no gate applied | **+0.0623** |
| quasi-static cells only | **+0.0444** |
| steady cells only | +0.0485 |
| both gates | +0.0438 *(but n = 1 seed at CV = 0)* |

**The headline falls ~30 %.** Sign and qualitative conclusion survive; the number does
not. Worse, the exclusions are **CV-asymmetric** — 3/3/1 non-quasi-static cells at
CV = 0/0.2/0.5 — which is the direction that inflates the contrast.

**Why it cannot simply be switched on:** applying both gates leaves **one seed** at
CV = 0. The low-CV end of the ladder is where cells fail the gates, because those rungs
sit nearest `φ_c` where the pressure is small and `I` is large.

### Task 0 (do this before anything else)

1. Make `analyse_phase8.py` report **both** the gated and ungated contrast, always, side
   by side. Never one alone.
2. Re-run the low-CV end with **more rungs on the jammed side** so that gated cells
   survive: at CV = 0 and 0.2, add offsets `+0.10, +0.12` (the ladder is currently
   `0, 0.02, …, 0.08` from the published `φ_J`). ~12 cells, under an hour.
3. Re-grade P1/P2/P3 on the gated numbers and **update every document that quotes
   +0.0623** — `SYNTHESIS`, `MECHANISM` §8, the derivation's §4 and §5 table, the Status
   card, and `salt_narrows_cv.py`'s `DPHI_C_DCV = 0.1247`.
4. Add a test that fails if a gate is defined and not applied.

⚠ **Everything below is written against +0.0623 and will need the same rescaling.**

---

## 1. Where the science stands

### Established (measured this session, 90 DEM cells, 0 failures)

- **`φ_c` rises with CV.** Δ`φ_c`(0→0.5) = +0.0623 ungated / +0.044 gated. Slope
  0.125 ungated, matching the independently-measured isotropic `φ_J` slope 0.1246.
- **🔑 P5, run last and the strongest discriminator: the effect is a distribution
  SHAPE effect, not a moment effect.** Uniform (zero-skew) family gives Δ`φ_c` =
  **+0.0215 ± 0.003** against lognormal's +0.0623 — a ratio of **2.90×** against a
  pre-registered bar of 1.5. The two families have `d₃₂` within ~9 % at matched CV, so
  **a permeability (Sauter-mean) account cannot produce this.** The isotropic analogue
  was 2.73×; the shear measurement gives 2.90×. *(This result is not yet written into
  `MECHANISM_2026-09-08.md` §8 — do that.)*
- **The headline, needing no fit:** at *matched* `φ`, raising CV 0 → 0.5 collapses the
  constant-volume grain pressure **5.4×–43×** (two overlapping rungs; quote the range),
  with rattlers 0.004 → 0.448. Robust to a late-strain window (4–5 %).
- **External anchor:** `φ_c(μ_p=0)` = 0.6347 ± 0.0010 vs an independently measured
  `φ_J(μ_p=0)` = 0.6332 — off by 0.0015, 7× inside the bar.
- **`φ_c` ≈ `φ_J`** to within 0.005 at both frictions.
- **`μ_s`** = 0.359 (μ_p = 0.3) vs 0.105 (μ_p = 0).

### The mechanism

`ID = ρκ(d₃₂(CV))U₀ / (A μ_w (φ − φ_c(CV)) D₀)`. The size distribution enters **twice**;
the permeability channel is `O(CV)` and vanishes at CV → 0, the critical-state channel is
`O(1)` and has a **pole** where `φ_c(CV) = φ`. Past it the pack is contractive, the pore
pressure changes **sign**, and the drop liquefies.

### Open, in order of how much it matters

1. **Degenerate with a simpler rival.** Plain suspension rheology (CV raises `φ_m`, so
   `η ∝ (φ_m−φ)⁻²` falls) consumes the same number and predicts the same direction. **No
   simulation separates them.** ⚠ But see §2 item 2 — P5 may weaken it.
2. **The lab's silica is on the wrong side of the measured `φ_c`** (0.44 vs 0.5695) and
   friction cannot close the gap (floor 0.5574, gap **0.117**).
3. **The drainage length is unmeasured and moves `ID` by 10⁵.**
4. **`α` (the dilatancy coefficient) is assumed CV-independent and is untested** — needs
   constant-*pressure* shear, which the protocol does not do.

---

## 2. The ambitious programme, in priority order

### 1. 🔬 A REPULSIVE CONTACT LAW — the highest-value simulation left

**Why.** It decides resolution (2) of the `φ_c` gap: whether a charge-stabilised pack has
a critical state anywhere near 0.44. Right now the mechanism cannot say which side the
lab's own material is on, which blocks every quantitative statement about the experiment.

**What.** Two routes, in increasing cost:

- **(a) Cheap, ~1 day.** Inflate interaction radii by the double-layer thickness while
  computing `φ` from the **true** radii. Yade: `Bo1_Sphere_Aabb(aabbEnlargeFactor)` +
  `Ig2_Sphere_Sphere_ScGeom(interactionDetectionFactor)`. This *is* the effective-diameter
  model, run properly.
  🔑 **And it tests something the algebra already predicts**: a constant shell narrows the
  effective PSD, `CV_eff/CV = d₁₀/(d₁₀+2h)` = 0.74 at low salt (`salt_narrows_cv.py`), so
  Δ`φ_c` should come out **~26 % smaller** than the hard-sphere value at the same true CV.
  If the DEM reproduces that factor, the effective-diameter model is validated; if not, it
  is not, and that is worth knowing before anything is built on it.
- **(b) Proper, ~1 week.** A genuine screened-Coulomb law with a finite-range force
  (custom `Law2_`), giving a **stress-dependent** `φ_c` — the critical-state *line*, which
  is also the closure `poro_impact.py` says its transient model is missing. Same gap
  surfaces in two independent places, which is evidence it is real.

**Pre-register before running:** `φ_c` falls below the hard-sphere value; the shift scales
as `(d/(d+2λ_D))³`; and Δ`φ_c`/CV falls by the narrowing factor.

### 2. Re-examine whether P5 breaks the degeneracy

P5 showed the CV effect is a **shape** effect (2.90× between families at matched `d₃₂`).
**Does the rheological rival predict that too?** `φ_m` is a packing-geometry quantity, so
it plausibly does — but *check it*, do not assume. If `φ_m` and `φ_c` respond differently
to skewness, P5 becomes a simulation-only discriminator and §1's "no simulation separates
them" is wrong. This is a literature question plus one extra DEM arm (bidisperse, which
already matches `d₃₂` to 4 decimals — the family exists in `psd.py`).

### 3. Constant-PRESSURE shear — closes two open items at once

A servo on the isotropic part of `velGrad` holding `p` constant while shearing. Gives
(i) the dilatancy coefficient `α` directly, testing the derivation's untested assumption,
and (ii) `φ_c(p)`, the critical-state line. ~15 lines on top of `shear.py`; the risk is
servo stability, so shake it out on one cell first.

### 4. E-6, the nozzle test — the cheapest bench experiment, and it attacks the biggest uncertainty

`drop_size_test.py`. The three drainage readings move the boundary velocity in three
**directions**: `U₀* ∝ D₀`, `∝ 1/D₀`, or flat — **3.6 / 14.4 / 7.2 m/s at 2 mm** against
a shared 7.2 at 4 mm. Decided by the **sign**, so it needs neither the band thickness nor
the `ID` prefactor nor where the boundary sits.

⚠ **Run it on cornstarch**, not silica (silica's boundary is 278 m/s). Score the
**regime**, not the spread factor — `We` and `Re` move with `D₀` on their own — and use a
Newtonian liquid of matched viscosity and density as the control.

### 5. Re-analyse the EXISTING salt series

`experiment_coverage.py` established that the salt ladder (0/10/50/100 mM, 4 levels) is
**the only designed axis already on disk**; `φ` spans 0.030 (2.1× too short for E-2) and
`d₀` is nozzle scatter. The effective-diameter model has `φ_c^eff` moving 0.23 → 0.56
across that ladder, so **the sign of `φ − φ_c` should flip inside data already taken**.

⚠ **First establish which regime criterion is current.** The project records that the
`arrested` flag "stayed True for 22 drops that never jammed after the classification
moved". In `phase_drops.csv`, 27 drops have a finite `d_arrest` but `arrested = False` —
find out whether that is the stale flag before building anything on it.

### 6. The 2×2 the salt axis actually needs

Two salt effects point **opposite ways** in the same model: screening moves `φ_c^eff`
itself (large, toward contractive) and *narrows* the effective PSD (small, weakening the
CV lever). On a one-dimensional ladder they are confounded. **Design: (CV low/high) ×
(salt low/high).** The interaction term is the narrowing; the salt main effect is the
`φ_c^eff` shift.

### 7. E-2, which is NOT blocked

Raise `φ` at high CV until the solid branch returns; the increment should equal Δ`φ_c`.
Needs `φ` steps ≲0.03 and drift <0.015 — an order of magnitude looser than the ±0.0025
this was wrongly written up as needing. ⚠ **The existing `φ` ladder is 2.1× too short**;
extend to ≈0.42–0.56 with ≥6 levels.

### 8. E-0, still the one thing only the author can do

Raw size distributions of the two bench samples. Decides (i) whether they are matched on
the **number** mean, (ii) whether a pure size threshold can work (it needs CV ≈ 1.0–1.5),
and (iii) whether "800 nm" is a number mean at all — `psd.py` notes laser diffraction
reports `d₄₃`, and if 800 nm is volume-weighted then `d₃₂` **falls** with CV and the sign
of the permeability channel reverses.

---

## 3. Standing hazards in this codebase

- **Estimator choice depends on friction.** E2 (pressure root) is primary at `μ_p` = 0.3;
  **E4 (isostatic crossing) at `μ_p` = 0**, because E2's exponent hits its bound there and
  E4 finds no crossing at 0.3. One rule is silently wrong in one arm.
- **E1 is retracted** — on a fixed-rate ladder `I → 0` is the *densest* rung, not the
  critical state. A test pins this; do not reinstate it.
- **A fit on a bound is not a measurement.** Guards exist for both the root and the
  exponent; both have fired on real data.
- **`sync_back.sh` moves only raw jsonl and logs.** Derived `.csv`/`.txt` are re-derived,
  never copied — it once clobbered fifteen committed outputs.
- **`E:\CPB\Drop Impact Analyses…` ends in an untypeable U+F028.** Glob for it; a typed
  path reads as missing data.
- **Never grade a knob at one setting of the knob it interacts with.** This session
  committed that error twice (the size-threshold operating point; the isostatic count).
- Raw data on `E:` is **read-only**.

## 4. Verification commands

```
cd C:\Users\compl\Desktop\Projects\polydisperse_jamming
python -m pytest test_shear.py test_poro.py -q          # 32 tests
python derivations/dilatancy_sign_mechanism_check.py    # 21 checks
python analyse_phase8.py                                # CV production arm
python analyse_phase8.py --pattern "phase8_p5uniform*.jsonl" --out p5_verdict.txt
python grade_p4_anchor.py                               # frictionless anchor
python phi_axis_offset.py ; python drop_size_test.py ; python experiment_coverage.py
python make_mechanism_figures.py --check                # byte-identical regeneration
```

Yade runs go through WSL:
`wsl bash -lc 'bash /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/launch_phase8.sh 14 <args>'`

## 5. Unfinished from this session

- The numerics reviewer **did not complete** (session rate limit). Its one delivered
  finding is §0. **Re-run it** — the remaining claims it was asked to check are
  independent re-fits of Δ`φ_c`, the matched-`φ` pairs, the P4 estimator choice, and
  whether the gate exclusions bias the headline (partly answered in §0).
- P5's result is **not yet written into** `MECHANISM_2026-09-08.md` §8 or the derivation's
  P-table. It currently says "NOT RUN".
