# A geometric account of the jamming window: size fills, shape and friction count

2026-09-25. Numbers: `derivations/geometric_window_2026-09-25_check.py` → `_check.txt` (algebra, lever arms,
polyhedra, drainage), `derivations/filling_factor_posthoc_2026-09-25.py` → `.txt` (existing records, POST HOC),
`grade_rollchi.py` → `out/grade_rollchi.txt` (the registered campaign `PREREG_ROLLCHI_2026-09-25.md`).
Every number below is printed by one of those three.

> **Result.** For every frictional or frictionless DEM pack in this project, the jamming packing fraction
> is well described by one geometric relation,
>
> **(G)  φ_J = Z / (Z + 2√3 · w)**
>
> (Song, Wang & Makse, Nature 453, 629 (2008), doi:10.1038/nature06981, where w = 1 for monodisperse
> spheres). **Z** is the mechanical coordination the contacts must supply for rigidity; it is set by
> **constraint counting** — friction, rolling resistance, and the rotational freedom that only a non-sphere
> has. **w** is a Voronoi **filling factor**; it is set by the **grain population** — the size distribution
> lowers it, elongation raises it — and is nearly friction-independent. So:
> - **size (CV) is a filling knob**: it lowers w (1.00 → 0.76, lognormal CV 0 → 0.5), which raises φ₀ and,
>   by the algebra of (G), also *narrows* the window slightly. At fixed weighed φ a polydisperse sphere
>   suspension drops out of the window from both ends.
>   (Referee pass 1 corrected the shape bullet below; see §5.)
> - **shape worsens filling on both branches and counting compensates on only one.** An elongated grain
>   raises w (1.00 → 1.37, the same at μ_p 0 and 0.3). Frictionless, a non-sphere's rotations are real degrees
>   of freedom that must be constrained (counting: Z_iso 6 → 10–12; measured Z₀ 5.97 → 7.95, hypostatic), and
>   the extra contacts offset the worse filling, so φ₀ barely moves (−0.9 %). Frictional counting is shape-blind
>   (Z_iso = 4 for any shape), so no compensating contacts are required; the measured Z_m rise (+1.0) is
>   geometric crowding that counting does not demand (Z counts DISTINCT neighbours — `packing.contact_graph`
>   deduplicates by grain pair — so it is not contact multiplicity either; see §5), and φ_m falls (−5.8 % for doublets
>   at μ_p 0.3). Net Z_m/Z₀ 0.775 → 0.705, window ×1.49.
> - **Every shape family tested here (doublets, α 2 clumps, their mixtures) has a wider window than spheres.**
>   Sphericity is not what lets silica jam into a ball; in these packs non-sphericity makes shear jamming
>   *easier*. Mild asphericity and real facets are open (§6, §9). External, weaker corroboration: same-paper
>   fits give χ = 0.833 for 1.5 µm silica spheres (Royer, Blair & Hudson, PRL 116, 188301 (2016),
>   doi:10.1103/PhysRevLett.116.188301) and 0.762 for cornstarch (Han, Wyart, Peters & Jaeger, PRF 3, 073301
>   (2018), doi:10.1103/PhysRevFluids.3.073301; φ convention not verified) — different labs and conventions.
>
> ⚠ (G) is an empirical-geometric closure: SWM derive it for monodisperse spheres with Z ∈ [4, 6]; its use
> with w ≠ 1 and Z outside [4, 6] is this project's extension and is tested, not assumed (§4, §5).
> **Registered campaign (§6, partly complete):** a new aspect-ratio ladder puts the window at 1 − χ = 0.072 for
> spheres and 0.10–0.11 for every clump from α 1.1 to 2.0 (p̂ 1e-3; the registered 3e-5 grade is INCOMPLETE).
> Mild asphericity (α ≤ 1.5) RAISES the frictional onset φ_m (more neighbours, same w); a closing neck lowers it,
> continuously (N1 GRADUAL). Spheres with a saturated contact moment limit reproduce the doublet's φ_m at μ_p 0.3
> (R3) but through the opposite geometry (Z ↓, w unchanged) and without the size interaction (R7 not clump-like);
> rolling resistance saturates near η ≈ 0.3, not 0.1 (R2 FAIL).

---

## 1. Frame

- Quantities: φ₀ (frictionless jamming, μ_p = 0), φ_m (frictional jamming; this project's clean realistic
  value is μ_p 0.3), χ = φ_m/φ₀, x = φ/φ₀. Convention (project CLAUDE.md / Status card): a suspension can
  shear-jam only if χ ≤ x < 1 (Wyart–Cates window); χ at a stated DEM pressure rung p̂.
- DEM: Yade, periodic, isotropic compression to reduced pressure p̂ = p/E, linear contacts, N = 1500–5400.
  "φ at p̂ 3e-5" stands in for the jamming point (`PREREG_KMER` convention); both rungs 3e-5 and 1e-3 are
  quoted wherever a result could depend on the rung.
- Z = z_nr, the mean contact number over non-rattlers (grain level for clumps).
- Assumptions: (A1) jamming is rigidity: Z reaches the isostatic count at φ_J; (A2) the relation between
  mechanical Z and φ is a local Voronoi-volume balance (SWM); (A3) the grain population enters that balance
  through one number w; (A4) friction and rolling resistance change Z, not w. A3–A4 are the testable content.

## 2. Counting (C1, SymPy)

For bodies with n_f degrees of freedom and n_c independent constraints per contact, rigidity needs
N n_f = (N Z/2) n_c, so

(1) Z_iso = 2 n_f / n_c.

| body / contact | n_f | n_c | Z_iso |
|---|---|---|---|
| frictionless sphere (rotations are null modes) | 3 | 1 | 6 |
| frictionless non-sphere, axisymmetric | 5 | 1 | 10 |
| frictionless non-sphere, general | 6 | 1 | 12 |
| frictional, sliding blocked (any shape) | 6 | 3 | 4 |
| + rolling blocked | 6 | 5 | 12/5 |
| + rolling and twisting blocked | 6 | 6 | 2 |

With a fraction x_s of contacts at the Coulomb limit and x_r, x_t at the rolling/twisting limits,
(2) n_c = 6 − 2x_s − 2x_r − x_t (checked: x_r = x_t = 1, x_s = 0 → 3; all zero → 6).

The row pair that matters: **frictionless, a sphere needs 6 and a non-sphere 10–12; frictional, both need 4.**
(The frictional shape-blindness is `geom_shape_2026-09-19` E1; the frictionless jump is Donev et al.,
Science 303, 990 (2004), doi:10.1126/science.1093010, who measure Z ≈ 10 for spheroids.)

## 3. From counting to packing fraction (C2, SymPy)

Take (G) with a = 2√3. The window is

(3) χ = φ(Z_m)/φ(Z₀) = Z_m (Z₀ + a w) / [Z₀ (Z_m + a w)],

(4) ∂χ/∂w = −a Z_m (Z₀ − Z_m) / [Z₀ (Z_m + a w)²] < 0 whenever Z_m < Z₀.

Limits (SymPy): w → 0 ⇒ χ → 1 (perfect filling closes the window); w → ∞ ⇒ χ → Z_m/Z₀; Z_m → Z₀ ⇒ χ → 1.
At w = 1 with ideal counts: spheres (6, 4) χ = 0.845; the measured sphere pack is hyperstatic at finite p̂
(Z₀ 5.97, Z_m 4.63 at 3e-5) and gives 0.9035 from (3) vs 0.9083 measured.

Two consequences follow from (3)–(4) alone:
- **Better filling narrows the window** (smaller w ⇒ larger χ), in addition to raising φ₀.
- **The window is set by the counting ratio Z_m/Z₀, tempered by w.** A knob that lowers Z_m/Z₀ opens it.

## 4. Size is a filling factor (POST HOC, existing records)

Invert (G) per cell: w = Z(1 − φ)/(a φ). Over the 70 existing (family, CV, rung) groups with both μ_p 0 and
μ_p 0.3 at N 1500 (`filling_factor_posthoc_2026-09-25.txt`):
- monodisperse spheres satisfy SWM with **no fitted number**: w = 1.001 (μ 0) and 0.989–0.992 (μ 0.3) at p̂ 3e-5.
- w falls monotonically with CV: lognormal 1.00 / 0.95 / 0.90 / 0.84 / 0.76 at CV 0 / .2 / .3 / .4 / .5
  (3e-5); uniform falls less (0.908 at CV 0.5, p̂ 1e-4 — no 3e-5 row exists) — the family dependence is the project's known "filling
  functional" result (`jamming-shift-is-a-filling-functional`), now expressed as one number per population.
- w is nearly friction-independent: median spread across μ_p 0 / 0.3 (/ 3) 0.026.
- **As a prediction it is approximate, not exact**: w from the μ_p 0 cell plus the measured Z at μ_p 0.3
  predicts φ_m with residual +0.0057 ± 0.0039 (max 0.014); only 31 % of groups within 2·PHI_SCALE. The sign
  is systematic (frictional packs fill ~1–3 % better per contact than (G) with the frictionless w says).
- μ_p 3 is where (G) is worst even for monodisperse spheres (w 0.910 at 3e-5): SWM's RLP end (0.536) sits
  below this protocol's 0.560.

**Why polydisperse spheres fail to jam (the geometric statement).** At fixed weighed φ, x = φ/φ₀ and
φ₀ = Z₀/(Z₀ + a w). CV 0 → 0.5 takes w 1.00 → 0.76, so φ₀ 0.632 → 0.694 and x falls by 9 %; by (4) the
window's lower edge rises (χ 0.9109 → 0.9170, both N 1500, lognormal, p̂ 3e-5; referee pass 1 caught a splice with the N 1700 value 0.908). A suspension sitting at x 0.95 as monodisperse spheres
(inside [0.908, 1)) sits at x ≈ 0.865 as a CV-0.5 population — below its χ — and can then only thicken.
Physically: small grains sit in the Voronoi cells of large ones and add solid volume without adding the
contacts rigidity is counted in (the rattler result, rattling ∝ a², `rattler_closure_2026-09-24`); the
same weighed φ is a looser network. The DNS already shows the drop consequence (`cvmap`: CV 0 jams, CV 0.5
spreads 35–67 % further, 99.5–101.8 % through the window).

## 5. Shape: counting, and why the sphere is special (POST HOC, existing records)

Same inversion on the doublet ladder and α = 2 clumps (N 1700, CV 0):

| grain | p̂ | Z₀ | w (μ 0) | Z_m | w (μ .3) | Z_m/Z₀ | χ meas | χ from (3) |
|---|---|---|---|---|---|---|---|---|
| sphere | 3e-5 | 5.97 | 1.001 | 4.63 | 0.989 | 0.775 | 0.9083 | 0.9035 |
| doublet w₂ 0.5 | 3e-5 | 6.66 | 1.128 | 4.96 | 1.145 | 0.745 | 0.8814 | 0.8874 |
| doublet w₂ 1 | 3e-5 | 7.95 | 1.364 | 5.61 | 1.371 | 0.705 | 0.8631 | 0.8652 |
| clump α 2 | 3e-5 | 7.92 | 1.359 | 5.63 | 1.381 | 0.711 | 0.8620 | 0.8685 |

- **w is set by the grain, not the friction**, for shape too (1.36 vs 1.37–1.38).
- **Z₀ rises by 2.0** (toward the frictionless DOF count 10, not reaching it: hypostatic, as for ellipsoids),
  **Z_m by only 1.0** — so Z_m/Z₀ falls 0.775 → 0.705 and the window opens ×1.49 (the `kmer` result).
- The two changes on the FRICTIONLESS branch cancel (w up: worse filling; Z₀ up: more contacts, as counting
  demands for a body with real rotational DOF), which is why φ₀ barely moves (−0.9 %, shapechi).
- On the FRICTIONAL branch counting is shape-blind (Z_iso = 12/(3 − x_m) for any shape, E1), so nothing
  demands extra contacts. The worse filling is therefore not compensated and φ_m falls.
- ⚠ Referee pass 1 (FATAL, acted on) rejected an earlier draft that explained the frictional +1.0 as partially
  released rotational DOF. Its proposed replacement ("the +1.0 is double-contact multiplicity, E1") is ALSO wrong
  on this code: Z = z_nr counts DISTINCT non-rattler neighbours (`packing.contact_graph` deduplicates every grain
  pair; `dup_contact_frac` is reported beside it and never enters Z). A doubly-touching pair supplies MORE
  constraints per neighbour, which, if anything, lowers the count rigidity needs. So the frictional rise is
  geometric crowding — an elongated grain touches more distinct neighbours — that counting neither demands nor
  forbids. What counting explains is the ASYMMETRY (the frictionless branch MUST add contacts; the frictional
  one need not); the SIZE of each rise (Z₀ +2.0, short of the +4 DOF count; Z_m +1.0) is measured, not predicted.

**Sphericity is a singular point (hypothesis).** Counting says Z₀ jumps from 6 to ≥ 10 at *infinitesimal* asphericity
(Donev's singularity), while frictional Z_iso = 4 does not move. If the DEM follows the counting, the window
should open most steeply right at α → 1, i.e. a slightly non-spherical grain (cornstarch-like) already has most
of the doublet's window. If it follows the lever arm (the eccentricity of contacts grows ∝ (α − 1)), the
opening is gradual. That is `PREREG_ROLLCHI` A2 — [§6].

## 6. The registered campaign (`PREREG_ROLLCHI_2026-09-25.md` → `out/grade_rollchi.txt`)

Run on a VM shared with ~26 solver processes of another session, so the grid did not complete in the session;
every test whose cells had not landed is INCOMPLETE, never graded on a subset. Execution order was changed
twice and the grid extended once (Amendments 1–3), all before grading and all committed.

**Arm R — spheres given a contact moment limit η (rolling + twisting), N 1700.**
- R0 inert half (μ_p 0: η 1 vs η 0 bit-identical): holds on every pair that has landed; live half: η 1 moves φ
  by −0.028 / −0.035 / −0.030 (CV 0 at 1e-3, 3e-5; CV 0.5 at 1e-3) against a 0.0042 bar. Registered verdict
  INCOMPLETE on pair count only.
- R1 (CohFrict η 0 = FrictMat): |Δφ| ≤ 0.0016 on every row with n = 3.
- **R2 FAIL (as registered).** Rolling resistance does NOT saturate by η = 0.1 at μ_p 0.3: φ(0.1) − φ(1) =
  +0.0037 ± 0.0004 (1e-3) and +0.0062 ± 0.0034 (3e-5), both above PHI_SCALE; η 0.3 is within 0.0025 of η 1. The
  2026-09-13 scoping probe (N 600, one seed, p̂ 3e-4) had said "saturates by 0.1". The half-effect scale is
  η ≈ 0.01–0.03: at p̂ 1e-3, φ(0.01) sits 0.0244 (n 3) and φ(0.03) 0.0157 (n 2) above φ(1) of a 0.0282 total
  (the 3e-5 η 0.01 / 0.03 cells are single seeds; referee pass 2).
- **R3 PASS (confirmatory — seen at seed 99 before registration).** Saturated rolling spheres land on the pure
  doublet's φ_m: 0.5674 vs 0.5693 (1e-3), 0.5399 vs 0.5414 (3e-5); the rolling channel carries 1.05–1.08 of the
  doublet's shift.
- **But not by the same geometry (POST HOC, (G) decomposition).** Rolling spheres keep w ≈ 1 (0.94–1.00) and
  reach that φ by LOWERING Z below the sliding-only count (Z 3.82 at η 1, p̂ 3e-5; rattlers 29 %). Doublets reach it
  with w 1.37 and Z 5.61. Same φ, opposite changes in Z. R3 therefore decides nothing about mechanism: φ alone
  cannot tell the two apart, and (G) can.
- **R7 at p̂ 1e-3: NOT CLUMP-LIKE** (the registered discriminator). ρ = D(η 1)/D(η 0) = 0.966, 95 % CI
  [0.809, 1.123] (n 3; informative, half-width 0.157 < 0.165), against 0.62 for α 2 clumps and a separable
  prediction of 0.953. Rolling resistance leaves the size effect intact; shape does not. 3e-5: INCOMPLETE.
- R5 (SWM with w = 1 on CV-0 spheres, all μ and η): 92.8 % of the 69 landed cells within 0.015 (mean residual
  +0.0004, sd 0.0085), including Z down to 3.15 — i.e. below the [4, 6] range SWM state it for. Registered verdict
  INCOMPLETE (69/156 cells). R4, R6 INCOMPLETE; R6's post-hoc per-rung variant repeats §4's systematic
  (+0.009 to +0.016: frictional packs sit denser than (G) with the frictionless w predicts).
  **Reading:** "shape is blocked rolling" is TRUE AS A NUMBER at realistic friction and FALSE AS A MECHANISM.

**Arm A — the aspect-ratio ladder, CV 0, N 1700 (0 of 34 checked cells K4-excluded).**

φ_m at μ_p 0.3 (seed means, n in brackets where < 3):

| α | 1.0 | 1.1 | 1.2 | 1.5 | 1.75 | 1.9 | 1.95 | 2.0 |
|---|---|---|---|---|---|---|---|---|
| p̂ 1e-3 | 0.5940 | 0.6013 (1) | 0.6151 (1) | 0.6289 (1) | 0.6086 | 0.5861 | 0.5769 | 0.5665 |
| p̂ 3e-5 | 0.5747 | 0.5792 | 0.5852 (1) | 0.5993 (2) | 0.5834 | 0.5606 | 0.5517 | 0.5407 |

and at p̂ 1e-3 the window 1 − χ is **0.072 (sphere) → 0.106 / 0.112 / 0.101 / 0.101 / 0.108** at α 1.1 / 1.2 /
1.5 / 1.75 / 2.0 (frictionless φ₀ 0.672 / 0.693 (1) / 0.700 (1) / 0.677 / 0.635; Baule et al.'s mean-field maximum
for dimers is 0.707, Nat. Commun. 4, 2194 (2013), doi:10.1038/ncomms3194 — the abstract value; the peak location
was not checked in the primary text).
- **N1 (Amendment 3): GRADUAL** at both rungs. φ_m falls continuously as the neck closes (r_neck/r_s 0.66 → 0.44
  → 0.31 → 0 at α 1.75 → 2.0), every step > PHI_SCALE. The collapse toward α = 2 is a smooth function of how
  concave the grain is, not an event at tangency.
- **The (G) decomposition splits shape into two regimes** (POST HOC): for α ≤ 1.5, w stays at the sphere value
  (0.995–1.008 at μ_p 0.3) and Z rises (4.63 → 5.22 at 3e-5) — denser packing through more neighbours at unchanged
  filling. For α ≥ 1.75, w climbs with neck closure (1.12 → 1.24 → 1.38) while Z saturates (5.4–5.6) — worse
  filling, looser packing. The window opens in both regimes (Z_m/Z₀ falls on the frictionless-branch rise alone).
- ⚠ **Neck confound.** At α 1.1, 32 % of frictionless neighbour pairs touch both lobes (10 % frictional): grains
  nestle into the shallow neck. A convex near-sphere (a faceted grain) has no neck. The 1.1 point therefore shows
  that a *clump* at 10 % asphericity already has most of the window-opening; it does not show that a convex
  faceted grain does.
- A1, A2, A3 (graded at p̂ 3e-5) are INCOMPLETE: no frictionless clump cell at 3e-5 finished in the session
  (they run 2–4 h each under the shared load). At 1e-3, A3's eccentricity map sits within ±0.0134 of every α
  (report only; not the registered rung).

**What the campaign changes in the argument.** (i) The "sphericity is singular" hypothesis of §5 is supported
at p̂ 1e-3 for clumps: every α from 1.1 to 2.0 gives 1 − χ clustered at 0.10–0.11, ~1.5× the sphere's. The
point-to-point differences (≤ 0.011) are the size of the propagated PHI_SCALE noise and α 1.2 / 1.5 rest on single
seeds, so the ladder is NOT resolved point-to-point: "opens by α 1.1" is supported, "then flat" is not tested (the
registered 3e-5 grade A2 is INCOMPLETE; referee pass 2). (ii) What lowers the ONSET
φ_m is concavity (a closing neck), continuously; mild asphericity RAISES φ_m. (iii) A contact moment limit
(rolling resistance) reproduces the doublet's φ_m at μ_p 0.3 but not its geometry and not its interaction with
size — so a faceted grain's lever arms (§7) are a separate mechanism from a clump's shape, not a proxy for it.

## 7. Faceted grains (cornstarch proxy; C4)

A faceted grain's contact normals miss its centroid. For vertex-on-face contacts over isotropic
orientations, the lever arm ⟨e⟩/r_eq is 0.35 (dodecahedron), 0.44 (icosahedron), 0.65 (cube), 0.76
(octahedron); rounding the edges with a sphere of radius equal to 50 % / 75 % / 90 % of the grain radius
leaves 0.22 / 0.11 / 0.04 (icosahedron). A two-sphere clump reaches 0.075 at α 1.1 and 0.62 at α 2 (C3,
Monte Carlo agrees with (π/4)(α − 1)(d_s/d_eq) to < 5e-4).

Against the measured rolling scales (§6, μ_p 0.3): half the rolling effect by η ≈ 0.01–0.03, saturation by η ≈ 0.3.
Sharp polyhedra (0.35–0.76) and an icosahedron rounded to half its radius (0.22) sit at or near saturation; one
rounded to 75 % (0.11) is part-way. **If** a faceted grain's lever arms act like a contact moment limit, cornstarch's
frictional φ_m would drop by up to the saturated rolling shift (−0.028 to −0.035 at μ_p 0.3). R7 says that channel
leaves the size effect alone, and the (G) decomposition says it works by removing contacts, not by changing filling
— a different mechanism from the clumps'. This is an inference from spheres with a moment law, not a measurement on
polyhedra; the project has no polyhedral DEM.

## 8. Drainage: why size, not shape, carries ball vs hybrid (C5)

Kozeny–Carman (form not verified here; used for RATIOS only) κ ∝ ψ² d_v² (1 − φ)³/φ²: cornstarch
(d_v 10–14 µm, ψ 0.8–1) against 0.8 µm Stöber silica is κ ×100–220, of which sphericity is at most ×1/1.56.
`ball_hybrid_scaling_2026-09-25` puts silica at De ≈ 7×10⁵ (undrained) and cornstarch at De ≈ 7, with the
contact modulus K_p (EDL vs Hertzian, ×1870) the larger lever. Shape moves neither by more than O(1).

## 9. Open
- (G)'s w for mild asphericity (arm A) and for real facets (no polyhedral DEM in this project).
- Why frictional packs sit systematically ~1 % denser than (G) with the frictionless w.
- The DEM φ axis is not the bench's; only ratios (x, χ) transfer (ANSWER §6).
- **w drifts with the DEM pressure rung** (lognormal CV 0 μ 0: 1.001 / 1.010 / 1.022 / 1.040 / 1.064 at
  p̂ 3e-5 … 3e-3; referee pass 1): (G) is a jamming-point relation, so every comparison here is at a FIXED
  rung, with 3e-5 the jamming stand-in. The drift (0.063 over two decades) is larger than the CV 0 → 0.2 step,
  so a cross-rung comparison of w would be meaningless.
