# Singlets, doublets, multiplets of ONE primary — and every knob of a dense-suspension drop impact

Written 2026-09-22 (fourth block). Checks: `multiplet_knobs_2026-09-22_check.py` →
`..._check.txt` (**9/9 PASS**; K6 is an arithmetic cross-check against `predict_kmer.py`,
not a physics test). Predictions: `predict_kmer.py` → `out/predict_kmer.txt`. Test:
`PREREG_KMER_2026-09-22.md` (126 DEM cells, graded; §11).
Status: **derived, checked; referee pass 1 returned 1 FATAL + 4 MATERIAL + 2 MINOR, all acted on below (marked "(pass 1)"); pass 2 confirmed them and returned 1 MATERIAL + 2 MINOR, acted on ("pass 2"). Pass 3 (after the data) returned 1 FATAL + 5 MATERIAL + 6 MINOR: the result box was rewritten from §11, §10 scoped to the μ_m 3 convention, stale pre-data sections marked, P-Z for clusters set UNDETERMINED ("pass 3"). Pass 4 returned 3 MATERIAL (net bound 0.58 → 0.59 over both rungs; Θ 7.88 → 7.78 from a mis-copied χ; triangles fail K4 at μ_p 0.3 too), acted on ("pass 4"). Pass 5 recomputed every result-box, §10 and §11 number from the artifacts and dumps: **nothing new**; the referee loop is closed.**

> ## RESULT BOX (rewritten from §11 after referee pass 3; the pre-data version is in git history)
>
> For a suspension of ONE primary sphere (diameter `d = 2a`) aggregated into rigid, tangent
> k-clusters, at the realistic friction μ_p = 0.3 and p̂ 3e-5 (the clean arm; `grade_kmer.txt`):
>
> 1. **Aggregation is a WINDOW knob, not a size knob.** Tangent clusters keep `S/V = 6/d`
>    exactly (K1) and compact k ≤ 4 clusters enclose no topologically closed pore (K2), so the
>    effective-sphere levers (Sauter diameter, trapped fluid) are absent. What moves is χ:
>    **measured 0.9083 → 0.8631 from singlets to pure dumbbells (window ×1.49)**, on the
>    ideal-mixing curve (12/12 clean cells within 2·PHI_SCALE; the largest gap is −0.0041 in φ at
>    w₂ 0.75, −0.0061 in χ; the μ 0.3 p̂ 3e-5 deviations run +0.0012 / −0.0022 / −0.0041 across
>    w₂ 0.25 / 0.5 / 0.75, so a small deficit growing with w₂ is not excluded). φ₀ moves 0.9 %, so the doublet fraction barely moves `x`.
> 2. **The composition that opens the window** (K4, measured clean endpoints):
>    `w*(x)` = 0.84 / 0.61 / 0.39 / 0.18 at x = 0.87 / 0.88 / 0.89 / 0.90. At fixed weighed φ in
>    `x ∈ [0.863, 0.908)` singlets cannot shear-jam and a sufficient doublet fraction can.
> 3. **Onset.** Every contact is primary–primary, so `F_rep = κa/2` is uniform and there is no
>    diameter-convention ambiguity. Measured `σ*(dumbbells)/σ*(singlets)` = **0.783** (isotropic)
>    and **0.786** (sheared, `grade_shearkmer.txt`).
> 4. **Sharpness.** Isotropic δ is unchanged at p̂ 3e-5 (+0.004, t 0.05) but falls −0.119
>    (t −3.79, below t_crit 4.303) at p̂ 1e-3; SHEARED it falls **1.254 → 1.063** (−15 %).
>    Under stress-controlled activation (§3 scope) the net stress to jam at FIXED composition is
>    ×0.67 at f_c 0.1 … ×1.05 at f_c 0.9; at FIXED weighed φ (the window opening lowers f_c too)
>    it is **at most ×0.59 anywhere in the window** (`net_fixed_x.py` → `out/net_fixed_x.txt`,
>    worst 0.579 at p̂ 3e-5, 0.593 at 1e-3).
> 5. **Isostatic counting.** Frictionless Z is linear in doublet number fraction (3/3) and
>    dumbbells are hypostatic (7.92 against 10). The rough Z = 4 prediction for CLUSTERS is
>    **UNDETERMINED**: every μ_p 3 cluster cell fails the K4 overlap gate (singlets: 4.04).
> 6. **In the drop** (§10): the only DNS test so far (`kmerdrop`) used μ_m 3 χ values, which are
>    the predicted, K4-excluded convention; at the clean μ 0.3 χ its x = 0.85 lies OUTSIDE the
>    window at every w₂. The clean test, `kmermap`, **passes R1**: doublets jam at x 0.88 where singlets cannot (J up to 0.196); see §10.
>
> Scope: everything in 4 and 6 assumes **stress-controlled activation** (M6), which rests on
> one informative DNS pair (`jamwin` P2) and censored bounds (`jamwin2`).

---

## 0. Frame and assumptions

**Sought:** how the k-mer distribution `{w_k}` (solid-volume fractions; `w₁ + w₂ + w₃ + … = 1`)
of one primary moves the reduced Wyart–Cates knobs `(x, χ, σ*, δ)` and the drop-impact
groups, stated as closed forms where geometry allows and as registered DEM predictions where
it does not.

- **(M1) Rigid, tangent clusters.** Members touch at a point; no necks. Sintered necks are
  treated in K1b and change one conclusion (permeability) in a known direction.
- **(M2) One primary**, diameter `d = 2a`, CV = 0. A polydisperse primary is the CV axis,
  already measured (`three_knobs.py`), and composes with this one through `x`.
- **(M3) Compact small clusters**: dumbbell (k = 2), equilateral triangle (k = 3). Open
  fractal aggregates (large k) are the `multiplet_psd` model's domain and are NOT this file's.
- **(M4) No breakage** in §1–§3; breakage is §4's separate knob.
- **(M5) The reduced framework** (`three_knobs.py`, `jam_stress_2026-09-22.md`): jammed iff
  `f̂ ≥ f_c = (1 − x)/(1 − χ)`, onset scale `σ*`, activation shape `δ`.
- **(M6) Stress-controlled activation at impact.** `jamwin` (graded 2026-09-22) found the
  solver's drop jamming more with sharper activation at `f_c` = 0.75 — the stress-controlled
  sign, not the rate-driven one (one informative pair). §3 uses the stress-controlled form
  and flags where the rate-driven one would differ.

## 1. Exact geometry of one-primary clusters

**K1 — specific surface is invariant.** A k-cluster of tangent spheres has `V = k·πd³/6` and
`S = k·πd²` (point contacts remove zero area), so `S/V = 6/d` for every k: the Sauter
diameter is `d` whatever the composition. With Carman–Kozeny, `κ = d_s²(1−φ)³/(180φ²)`, the
permeability at fixed φ is independent of the singlet/doublet ratio **only if the Kozeny
constant (180 here) is shape-independent** — it encodes tortuosity, and a dumbbell's concave
waist ("as corrugated as a body can be", `shape.py`) can raise it. S/V is exact; κ is not
(pass 1). (For a SINTERED
dumbbell at centre separation `s < d`, the removed caps lower `S` faster than `V`, so
`d_s` = 1.08 d at `s = 0.8d` and 1.12 d at `0.6d` — K1b, checked against Monte-Carlo volume.
Necks raise κ.) *This is the opposite of polydisperse spheres, where κ ∝ d₃₂² moves with CV.*

**K2 — no trapped fluid for compact small clusters.** The throat between three mutually
tangent primaries has inscribed radius `(2/√3 − 1)a = 0.155a > 0`, so the pore inside a
triangle or tetrahedron of tangent spheres connects to the outside: there is no
TOPOLOGICALLY closed pore, so **Λ ≈ 1** for k ≤ 4 compact clusters. ⚠ The hydrodynamic
criterion (below) is not computed: whether fluid in a tetrahedral interstice drains through
a 0.155a throat on the impact timescale is an open estimate (pass 1). The
trapped-fluid lever needs OPEN aggregates whose interior is screened hydrodynamically (the
Brinkman screening length of the aggregate interior smaller than its size), i.e. large k.

**K3 — isostatic counting.** A rigid body has 6 dofs; a frictionless body feels only normal
forces, and a tangent dumbbell's spin about its own axis is unconstrained, so its effective
dofs are 5 (a sphere 3, a triangle 6). Isostatic frictionless packs have `Z = 2⟨dof⟩`:
`Z = 2(3n₁ + 5n₂ + 6n₃)` in number fractions — 6, 10, 12. Perfectly rough contacts carry 3
constraints each, so `Z = 2·6/3 = 4` for **every** rigid body. Measured (`shapechi`, p̂ 3e-5):
frictionless singlets 5.97, dumbbells 7.92 (**hypostatic** by 21 %), μ_p = 3: 4.04 and 3.88.

**K3′ — the uniform threshold.** Every contact between clusters is a contact between two
primaries, so its reduced radius is `R* = a/2` and its Derjaguin threshold `F_rep = κa/2` is
the same everywhere. The size-convention trap that dogged the CV onset (σ\* shift 0.736 at
`d₁₀`, 1.145 at `d₃₂`) **does not exist here**.

## 2. The reduced knobs as functions of composition

Write `A = 1/φ` (specific volume). **H_mix** (ideal mixing): `A_J(w) = (1−w)A_J,s + wA_J,d` for
each branch (φ₀ at μ 0; φ_m at μ_m). Justification, and its limit (pass 1): all grains share one primary and the body-size ratio
is only 2^{1/3} = 1.26. This project's own SPHERE binaries at size ratio 1.25
(`phase12_mixr125`, μ_p 0.3, p̂ 1e-3) show a filling excess of at most **+0.0036** over the
monodisperse value (0.5976 at fines volume fraction X_s = 0.4 vs 0.5940) — but that single
point is **not resolved**: seed-paired t = 2.06 against t_crit 4.303 (n = 3). Its SIGN is
corroborated (20 of 24 seed-paired differences over the eight X_s rungs are positive; not
independent, the same three seeds). So size-ratio filling is small, of either sign within
noise at a single composition, and the non-spherical excess is untested (pass 2). (An
earlier version cited a "~2.4 filling threshold" that is not traceable in this repository;
withdrawn.) The rival — a positive excess from mixing shapes — is registered beside it.

- **Position.** `x(w) = φ/φ₀(w)`. φ₀ moves only 0.9 % from singlets to dumbbells, so **the
  doublet fraction barely moves x.** (Contrast CV, which moves x by ~a window.)
- **Contrast.** `χ(w) = A₀(w)/A_m(w)` — endpoints 0.886 → 0.788 (μ_m = 3). Window
  width `1 − χ` ×1.85. ⚠ (pass 3) The μ_m 3 dumbbell endpoint is a K4-excluded clump pack
  (§11); the clean measurement is μ_p 0.3: 0.9083 → 0.8631, window ×1.49.
- **Window opening (K4).** Setting `χ(w) = x` is LINEAR in `w`:
  `w* = (xA_m,s − A₀,s)/[(A₀,d − A₀,s) − x(A_m,d − A_m,s)]` (SymPy, K4a);
  0.87 / 0.54 / 0.34 / 0.15 at x = 0.80 / 0.83 / 0.85 / 0.87 (μ_m = 3, excluded convention).
  At the realistic μ_p = 0.3 (measured, clean) the window goes 0.092 → 0.137 and `w*` =
  0.84 / 0.61 / 0.39 / 0.18 at x = 0.87 / 0.88 / 0.89 / 0.90; it exists only for x ≥ 0.863.
- **Required activation** inside the window: `f_c(x, w) = (1 − x)/(1 − χ(w))` FALLS as `w`
  rises — a doublet-rich suspension needs a smaller frictional fraction to jam.
- **Onset (K6).** At fixed p and φ the mean contact force is `⟨f⟩ ∝ p/(n_c ℓ)` (Love–Weber:
  n_c contacts per volume, ℓ the mean contact lever arm about the body centroid). External
  contacts per primary: 4.63 (singlets) vs 2.82 (dumbbells) at μ_p 0.3 → `n_c` ratio 0.61;
  the lever arm of a lobe contact is between 1 and √2 times a sphere's (a HEURISTIC bracket:
it assumes contacts uniform over the lobe surface and ignores the exclusion zone next to the
sibling lobe — pass 1 MINOR; the `kmer` dumps measure the ratio directly). So
  `σ*(w=1)/σ*(0) = (n_c ℓ)_d/(n_c ℓ)_s ∈ [0.61, 0.86]` — the onset stress FALLS with doublets.
  (pass 3) **Measured 0.783 isotropic, 0.786 sheared** — inside the bracket.
- **Sharpness.** δ(w). ⚠ **The first version of this bullet predicted δ(dumbbells) <
  δ(singlets) from the frictional-dumbbell force CV of 1.3–8.8. That CV is an OUTLIER
  statistic, not a width**: the `kmer` per-contact dumps show 2–7 contacts per cell out of ~3500–4700
  carrying 50–143× the mean force in frictional clump packs (2–7 per cell across the dumps, pass 2) (μ_p 3 doublets, μ_p 0.3
  triplets), while the 99th percentile is 2.2–3.2× the mean — the same as spheres.
  Whether those few contacts are real wedging in a dumbbell's neck or an overlap artefact is
  open. The quantile estimators (δ₁₀/₉₀, σ\* at the 1/e quantile) are insensitive to five
  outliers in 4000, so the registered H_δ is graded as registered, but its REASON is
  withdrawn: there is no longer a mechanism on record predicting its sign.
  (pass 3) **Measured:** isotropic δ unchanged at p̂ 3e-5 (+0.004, t 0.05) and −0.119 (t −3.79,
  unresolved at n 3) at p̂ 1e-3; sheared, 1.254 → 1.063 (§11).

## 3. What a drop sees — every channel, one composition

With (M6), a material point jams where the impact stress exceeds
`σ_jam(w) = σ*(w)·[ln 1/f_c(x,w)]^(−1/δ(w))`. Doublets lower `σ*` and lower `f_c`; both lower `σ_jam` unconditionally. A change in `δ`
has an `f_c`-dependent sign: by `jam_stress_2026-09-22.md` eq. (9), a LOWER δ lowers `σ_jam`
for `f_c < 1/e` and RAISES it for `f_c > 1/e` (an earlier version had this backwards — pass 1
FATAL). The doublet direction of δ was unmeasured when this was written; sheared, it FALLS (§11),
so near the window edge the δ channel offsets the σ\* channel at FIXED composition (×1.05 at
f_c 0.9). At FIXED weighed φ the window opening also lowers `f_c`, and the net is ≤ ×0.59
everywhere in the window (`out/net_fixed_x.txt`; pass 3). All of this is under (M6). **Outside** the window the drop cannot jam but
still thickens by `Θ = [(1−x)/(1−x/χ)]^n`, which RISES as χ falls (K5). And the
inertial-dilatant number `ID ∝ κ/(φ − φ_c)` FALLS: κ is unchanged (K1) while the critical
packing φ_c tracks the frictional branch, which dumbbells lower by 12 %.

So the window, onset, thickening and dilatancy channels move the same way — toward jamming,
earlier activation, stronger thickening, stronger dilatancy — at fixed primary size,
chemistry and weighed φ; the activation-shape channel is the exception whose sign is open.

## 4. The complete knob atlas for a dense-suspension drop impact

Reduced variables: position `x = φ/φ₀`, contrast `χ = φ_m/φ₀`, onset `σ*`, sharpness `δ`,
permeability `κ` (sets dilatancy/Darcy times), and the drop groups
`We = ρU₀²D₀/γ`, `Re = ρU₀D₀/η`, Stokes `St = ρ_p γ̇ b²/η_s`, and the stress ratio
`Π = σ_imp/σ*` (gauge named per use).

| particle / process knob | x | χ | σ\* | δ | κ | status | the decisive test |
|---|---|---|---|---|---|---|---|
| packing fraction φ | ✔ | — | — | — | ✔ | definition | — |
| **size distribution CV** | ✔ (~a window per ΔCV 0.5) | ✗ (paired \|t\| 0.56) | ✔ (convention-dependent sign) | ✔ (δ −19 % at CV 0.5, μ 0.3) | ✔ d₃₂² | MEASURED | bench at fixed φ, two CVs, same d₅₀ |
| **singlet/doublet ratio (this file)** | ~0 (0.9 %) | ✔✔ (window ×1.49 at μ_p 0.3, clean) | ✔ −22 % measured (0.783 iso., 0.786 sheared) | ✔ sheared −15 %; isotropic unresolved | ✗ if the Kozeny constant is shape-free (K1) | MEASURED (μ_p 0.3); μ_m 3 excluded (K4) | bench: fractionate or grow doublets, same primary, same φ |
| **open aggregates (large k)** | ✔ via Λ ≥ 1 | ✔ (shape, unmeasured) | ✔ | ? | ✔ ↑ | modelled (`multiplet_psd`) | salt/pH-aggregated suspension vs dispersed |
| friction μ_p | — | ✔✔ (window 0.113 at 0→3) | — | ✔ (δ −0.22 at 0→0.3) | — | MEASURED | coated vs bare particles |
| shape (elongation) | ~0 | ✔✔ | ? | ? | ✔ (via S/V) | χ MEASURED | rods vs spheres, same volume |
| surface potential κ | — | — | ✔ (∝ κ) | — | — | algebra | salt at fixed particle |
| primary size a | — | — | ✔ (∝ a^(m−2), m ∈ {0,1}) | — (A4b untested) | ✔ ∝ a² | algebra + one size series | size ladder at fixed chemistry — degenerate with κ unless κ is varied |
| stiffness E | — | — | — | — | — | MEASURED scale-free | — |
| particle density ρ_p | — | — | — | — | — | via St only | density-matched vs not |
| **breakage of doublets** (§5) | — | ✔ reverses during impact | ✔ | ? | — | criterion only | sintered vs salt-bridged doublets |
| drop U₀, D₀, fluid η_s, γ | — | — | — | — | — | set Π, We, Re, the film thickness | the impact ladder |
| **Brownian motion (Péclet)** | ✔ (thermal φ₀ shift) | ? | ✔ (entropic repulsion adds to κ) | ✔ | — | breaks A4b (a size-free ĝ); not in the DEM | Pe ladder via solvent viscosity at fixed a |
| **lubrication film vs particle size** (h_film/a) | — | — | ✔ (contact vs lubricated) | ? | — | not modelled; the DNS has no particle scale at DIL_COEFF = 0 | drops of one D₀ at two a |
| **particle inertia** (St = ρ_p γ̇ b²/η_s; doublets ×2^{2/3}) | — | — | — | — | — | not modelled | density-matched vs not |

**What the atlas says about the question.** Only three knobs move χ — friction, shape, and
(through shape) the singlet/doublet ratio — and χ is what decides whether jamming is
possible at all. Of those, friction is a surface property the bench cannot easily vary, shape
needs a different synthesis, and **the singlet/doublet ratio is set by the very process that
makes the suspension** (drying, salt, sintering, sonication). That is the argument for
aggregation being where the answer lies, and it is testable.

## 5. Breakage — the knob that can reverse mid-impact

If doublets are held by a finite bond force `F_b` (van der Waals or a salt bridge, not a
sintered neck), the impact stress can split them. The dimensionless group is
`Π_b = σ_imp a²/F_b` (K7). For `Π_b ≪ 1` doublets survive and §3 applies; for `Π_b ≫ 1`
they split during loading, `w` falls, the window NARROWS while the drop is being loaded, and
a drop that began jamming can un-jam — a history the reduced framework's static `χ` cannot
represent. **No number is put on `F_b` here**: a van der Waals estimate needs a Hamaker
constant and a contact separation, neither measured for the lab's particles
**[UNVERIFIED]**. The test is binary and cheap: the same composition prepared with sintered
vs salt-bridged doublets must behave identically if `Π_b ≪ 1` and differently if not.

## 6. Predictions, and what would falsify each

| # | prediction | falsified by | where |
|---|---|---|---|
| P-mix | φ(w) follows H_mix within 2·PHI_SCALE | a systematic excess > 0.0042 | `kmer` |
| P-χ | χ(w) 0.886 → 0.788 near-linearly; `w*(0.85)` ≈ 0.34 — (pass 3) the μ_m 3 arm is UNDETERMINED (K4); the clean μ_p 0.3 arm falls 0.9083 → 0.8631, monotone | χ flat in w, or non-monotone | `kmer` |
| P-Z | rough Z = 4 ± 0.25 at every composition | Z drifting with w at μ_p 3 | `kmer` |
| P-σ\* | σ\*(dumbbell)/σ\*(singlet) ∈ [0.55, 0.95] at μ_p 0.3 | ratio ≥ 1 | `kmer` dumps |
| P-δ | registered as δ(dumbbell) < δ(singlet), but its mechanism is WITHDRAWN (§2): graded as registered, read as a measurement of an open sign, not as a test of a mechanism | — | `kmer` dumps |
| P-bench | at fixed φ and primary, raising doublet volume fraction past `w*(x)` switches the drop from spreading to arrest | no switch across the whole w range at an x inside [χ_d, χ_s) | bench |
| P-κ | tangent aggregation does NOT change the Darcy time at fixed φ | a measured permeability change at fixed φ with tangent (non-sintered) doublets | Darcy cell |

## 7. Open

- H_mix passed 12/12 clean cells; a small deficit growing with w₂ at μ_p 0.3 p̂ 3e-5 is not excluded (pass 3).
- Isotropic packs are not sheared packs — (pass 3) now measured sheared (`shearkmer`): σ\* reproduces, δ does not.
- Triplets: the chain trimer widens the μ 0.3 window no more than a dumbbell (`kchain`); the compact triangle is K4-excluded at μ 3 AND at μ_p 0.3 (pass 4).
- The lab's working φ (0.44) sits below every DEM window (09-09 record); the transferable
  statements are in `x` and in RATIOS, not in raw φ.

## 8. The bench protocol, in quantities a rheometer measures (scale-free by construction)

The DEM φ axis does not transfer to the bench (the lab's φ = 0.44 sits below every DEM window;
09-09 record). What transfers are **ratios measured on one material**, because the reduced
framework is exactly covariant in `(φ, φ₀, φ_m)`:

1. **χ is measurable without knowing φ's absolute scale.** From steady flow curves at several
   φ: `φ₀` = divergence of the LOW-stress viscosity branch, `φ_m` = divergence of the
   HIGH-stress branch (both fits of `η ∝ (1 − φ/φ_J)^(−n)`, same n). `χ = φ_m/φ₀`.
   **Prediction: `χ(w)/χ(0)` falls from 1 to 0.95 (μ_p = 0.3, measured clean: 0.8631/0.9083 =
   0.950; the μ_m = 3 convention's 0.89 is K4-excluded) as the doublet volume fraction goes
   0 → 1**, near-linearly (H_mix), independent of
   whatever sets the absolute φ scale of the real particles.
2. **The onset ratio.** `σ*(w)/σ*(0)` read at the thickening onset of the same flow curves:
   0.786 at w = 1 in sheared DEM packs (pass 3; predicted [0.61, 0.86]).
3. **The drop.** With `x = φ/φ₀(w)` measured on the rheometer, a drop at fixed weighed φ
   switches from spreading to arrest when the doublet fraction crosses `w*(x)` (§2), computed
   from the SAME material's measured endpoint χ's — no DEM number enters.
4. **Controls that make it a test and not a correlation.** Vary w at FIXED chemistry: fractionate
   one batch (centrifugation separates singlets from doublets) and remix, rather than
   aggregating with salt — salt moves `κ` (and so σ\*) at the same time, which is the
   size–potential degeneracy of `onset_size_degeneracy_2026-09-19.md` in another form. Count
   w by microscopy on the remixed batches. Run a sintered-doublet batch beside it to read the
   breakage number `Π_b` (§5) off the difference.

## 9. What the lab's own salt series already says (existing data, no new experiment)

Salt screens charge-stabilised silica — which lowers the repulsive threshold κ — **and** lets it
aggregate. The two predict different things, and the bench has measured both channels:

- **Onset.** `impact_stress/results.py` LAB records the fitted onset at 0 / 10 / 50 / 100 mM:
  γ̇_c = 631.6 / 478.8 / 470.6 / 195.5 s⁻¹ and u₀\* = 1.473 / 1.120 / 1.138 / 0.455 m/s. At
  100 mM the onset has fallen to **×0.310** (γ̇_c) or ×0.309 (u₀\*; ×0.095 in u₀\*²).
- **Arrested state.** lessons.md (2026-08-30): among arrested drops at matched velocity, salt
  moves the footprint by only **×1.28 against a ×6.35 total range** — "salt shifts the ONSET,
  not the arrested state".

Read through this file:

1. **Doublets alone cannot produce the onset shift.** Their predicted onset factor is ≥ 0.61
   even at w₂ = 1 (§2), against a measured 0.31 on the most generous (viscous) gauge. So if
   aggregation drives the salt onset, it is not by doublets; it needs larger aggregates.
2. **Aggregation of any size must also move the window** (shape → χ, K3/§2), and hence the
   arrested state. The arrested state barely moves.
3. **So the lab's salt effect reads as SCREENING (κ ↓, σ\* ↓ at fixed χ), not aggregation.**
   Breakage does not rescue the aggregation reading: aggregates that broke under impact
   (`Π_b ≫ 1`, §5) would stop moving the arrested state, but they would ALSO stop lowering
   the onset, which the data show falling ×3. ⚠ This is a reading of two recorded numbers
   against a model; no data on disk measures the aggregation state itself.
4. **The decisive bench measurement is therefore the doublet fraction itself** (microscopy or
   analytical centrifugation) at 0 and 100 mM. If w₂ is ≈ equal, the salt series is a κ
   ladder, and it is the cleanest σ\* lever the lab already owns (it moves the onset by ×3
   without moving the window). If w₂ rises strongly with salt, the arrested-state constancy
   is a puzzle for this theory and must be reported as such.

## 10. The DNS consequence (`PREREG_KMERDROP`, graded 2026-09-22)

Feeding the DEM-measured χ(w₂) into the Wyart–Cates drop solver at a FIXED weighed φ (x = 0.85,
σ\* 300 Pa, δ 1.3169 — only the window channel moves): the jammed volume goes **0 / 0 / 0.049 /
0.220 / 0.385** at w₂ = 0 / 0.25 / 0.5 / 0.75 / 1 (D2 PASS; the zeros at w₂ 0 and 0.25 are by construction, since x = 0.85 < χ there puts the drop outside the window), and the maximum spread falls from
β 1.39 to 1.08 (continuous radial Rg 0.460 → 0.348, strictly monotone; the registered β test
FAILED on a one-grid-quantum tie between w₂ 0.75 and 1). So, within the closure, **raising the
doublet fraction of one primary at fixed φ switches the drop from spreading to arrest** — the
author's hypothesis in its continuum form. ⚠ It inherits every DEM → DNS caveat (χ from the
μ_m = 3 convention; the DEM φ axis is not the bench's; stress-controlled activation assumed
from jamwin's one informative pair).

🔴 **Scope (pass 3 FATAL).** The χ(w₂) fed to `kmerdrop` were the μ_m 3 PREDICTIONS
(0.886/0.859/0.834/0.811/0.788), and the only measured μ_m 3 cluster endpoint is K4-excluded.
At the clean μ_p 0.3 χ (0.9083 → 0.8631), x = 0.85 lies **outside the window at every w₂**, so
this run cannot be quoted as "doublets switch a drop to arrest" for a realistic-friction
material. It shows what the solver does IF χ moves by that much. The clean-χ test is
`PREREG_KMERMAP` (x 0.88/0.92, U₀ 1/2/4). **Graded 2026-09-22 (after the pass-5 close):** with
χ, σ\* and δ all from the clean μ_p 0.3 DEM, pure doublets jam at x 0.88 (J 0.058 / 0.196 at
U₀ 2 / 4) where singlets cannot (R1 PASS), and at x 0.92 they jam more (0.636 vs 0.138 at U₀ 4;
other pairs censored by L9 deaths). The largest spread cut, however, is at x 0.85, where nothing
jams: the footprint does not discriminate jamming from thickening. Stress-controlled scope (M6) as
before.

**Which channel (`PREREG_KMERDECOMP`, graded 2026-09-23).** Holding the singlets' σ\* and δ and
swapping in only the doublets' χ reproduces the x 0.88 jamming (J 0.0591 vs 0.0579 all-channel, U₀ 2)
and the whole x 0.85 footprint cut at U₀ 4 (ln Rg −0.372 vs −0.371). Swapping in only the doublets' σ\* and δ
moves the footprint by 0.08 % at U₀ 4, growing to 2.9 % at U₀ 1 (Amendment 1). So **the doublets
act on the drop almost entirely through the window**: the onset and sharpness channels measured in
the DEM are real, but minor in this DNS (≤ 3 % in Rg over U₀ 1–4). Grid-checked: at L10 the activation cuts are −0.0009 (U₀ 4) and −0.0280 (U₀ 1) against L9 −0.0008 / −0.0286, and the window cut is −0.390 (`PREREG_KMERDECOMP` Amendment 2, F1 PASS).

**Why spreading falls before any jamming (K5 made quantitative).** Outside the window the
activated drop can only thicken, by `Θ = [(1−x)/(1−x/χ)]^MP_EXPO` (MP_EXPO = 2 in the solver).
At x = 0.85 with the μ_m 3 PREDICTED χ: Θ = 13.6 at w₂ = 0 (χ 0.886) and 205 at w₂ = 0.25
(χ 0.859), ×15. ⚠ (pass 3) That step is ill-conditioned — 0.859 sits 0.009 from x, so Θ is
set by a predicted, K4-excluded χ. With the CLEAN μ_p 0.3 χ at the same x: Θ = **5.46 → 7.78
(×1.43)** for w₂ 0 → 0.25 and 5.46 → 97.7 for w₂ 0 → 1 (χ 0.8631, still outside the window).
That this thickening is what lowered the spread in the K5 panel between w₂ 0 and 0.25 is a
HYPOTHESIS, not measured. ⚠ Θ is quoted WITH its exponent (lessons.md); at MP_EXPO = 5 the
predicted-χ step is ×877.

## 11. What the DEM then measured (`PREREG_KMER` RESULT, 126 + 36 cells, graded 2026-09-22)

| prediction (§6) | outcome |
|---|---|
| P-χ (window opens with doublets) | ✅ at realistic friction μ_p 0.3, clean: χ 0.908 → 0.863 (window ×1.49), on the ideal-mixing curve within ±0.005. At μ_m 3 UNDETERMINED: those clump packs interpenetrate (below). |
| P-mix (ideal specific-volume mixing) | ✅ 12/12 clean cells within 2·PHI_SCALE (as registered: INCOMPLETE, the μ 3 half excluded). No excess; at μ 0.3 p̂ 3e-5 the deviation runs +0.0012 / −0.0022 / −0.0041, so a small deficit is not excluded. |
| P-Z (isostatic counting) | frictionless Z linear in doublet number fraction (3/3 PASS); frictionless dumbbells hypostatic. Rough Z = 4: singlets 4.04; **clusters UNDETERMINED** — every μ_p 3 cluster cell is K4-excluded (the grader earlier auto-passed the p̂ 1e-3 rows; fixed, pass 3). |
| P-σ\* (onset falls) | ✅ σ\*(pure doublets)/σ\*(singlets) = 0.783 at μ_p 0.3 (isotropic packs). |
| P-δ (sharpness) | isotropic: unchanged (+0.004, t 0.05). ⚠ **SHEARED (`PREREG_SHEARKMER`): δ falls 1.254 → 1.063** at μ_p 0.3 (−0.191 ± 0.008) — so the §3 cancellation IS real under shear, but only at the window edge: σ_jam ratio ×0.67 at f_c 0.10, ×0.94 at 0.75, ×1.05 at 0.90 (`out/shearkmer_sjam.txt`). σ\* ratio reproduced under shear (0.786). |
| chains vs dumbbells | a straight trimer widens the μ 0.3 window no more than a dumbbell at the quoted rung (0.8638 vs 0.8631). |

🔴 **Frictional (μ_p 3) clump packs lock in genuine interpenetration** (penetration depth up to
0.74 radius, not a time-step artefact, worse at lower pressure). ⚠ (pass 4) Measured from the dumps as u/r = fn/0.25, the worst is 0.92 in μ 3 dumbbell mixtures and **1.31** in triangle-bearing cells; and **triangle compositions fail K4 at μ_p 0.3 too** (worst 1.09). Dumbbells are the only shape clean at μ_p 0.3, and even their K4-clean cells hold a few contacts (under the 0.1 % gate) up to u/r 0.62. Every μ_m 3 clump number —
including `shapechi`'s window ×1.85 — carries that caveat; the overlap-volume sensitivity says the
χ values survive (≤ 0.0017), the forces do not.

**Net, for the question:** at realistic friction (μ_p 0.3, clean), turning a suspension of
singlets into doublets of the same primary widens the shear-jamming window by half (×1.49) and
lowers the onset stress by a fifth (×0.78). Under SHEAR it also softens the activation (δ −15 %),
which at FIXED composition lowers the stress to jam deep in the window (×0.67 at f_c 0.1) and
slightly raises it at the edge (×1.05 at f_c 0.9). At FIXED weighed φ the window opening also
lowers f_c, and the net stress to jam is ≤ ×0.59 everywhere in the window (`out/net_fixed_x.txt`,
recomputed pass 3). ⚠ These σ_jam statements assume **stress-controlled activation** (M6),
supported by one informative DNS pair; under rate-driven activation they would differ. The DNS
has shown the window channel switching a fixed-φ drop to jamming only with the μ_m 3 predicted χ
(§10 scope); the clean-χ DNS test is `kmermap`.
