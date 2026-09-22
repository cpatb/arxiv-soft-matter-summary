# Singlets, doublets, multiplets of ONE primary — and every knob of a dense-suspension drop impact

Written 2026-09-22 (fourth block). Checks: `multiplet_knobs_2026-09-22_check.py` →
`..._check.txt` (**9/9 PASS**; K6 is an arithmetic cross-check against `predict_kmer.py`,
not a physics test). Predictions: `predict_kmer.py` → `out/predict_kmer.txt`. Test:
`PREREG_KMER_2026-09-22.md` (126 DEM cells, running at time of writing).
Status: **derived, checked; unrefereed.**

> ## RESULT BOX
>
> For a suspension of ONE primary sphere (diameter `d = 2a`) aggregated into rigid, tangent
> k-clusters (singlets, dumbbells, triangles, …):
>
> 1. **Aggregation is not a size knob. It is a WINDOW knob.** Tangent clusters leave the
>    solid's specific surface exactly unchanged (`S/V = 6/d`, K1), so the Sauter diameter
>    and the Carman–Kozeny permeability at fixed φ do not move; and a compact cluster of
>    k ≤ 4 encloses no fluid (throat radius 0.155 a > 0, K2), so the trapped-fluid factor is
>    **Λ = 1 exactly** — the lever the effective-sphere model (`multiplet_psd_2026-09-22.md`)
>    is built on is absent for singlet/doublet/triplet populations. What survives is SHAPE:
>    measured pure doublets have `χ` = 0.788 against 0.886 for singlets (μ_m = 3; 0.862 vs
>    0.908 at μ_m = 0.3), i.e. **a doublet fraction widens the shear-jamming window.**
> 2. **The composition that opens the window** at fixed position `x` is, under ideal
>    specific-volume mixing, a closed linear-fractional form (K4):
>    `w*(x) = (x/φ_m,s − 1/φ_0,s) / [(1/φ_0,d − 1/φ_0,s) − x(1/φ_m,d − 1/φ_m,s)]`
>    — 0.87 / 0.54 / 0.34 / 0.15 at x = 0.80 / 0.83 / 0.85 / 0.87 (μ_m = 3). **A suspension
>    that cannot shear-jam as singlets CAN at the same weighed φ once ~a third of its volume
>    is doublets** (x = 0.85). This is the tractable, bench-testable form of the hypothesis.
> 3. **Every contact is primary–primary**, so the repulsive threshold `F_rep = κR*`
>    (R* = a/2) is UNIFORM. The onset then moves only through the force per contact, which
>    rises because a dumbbell exposes fewer contacts per primary: **σ\*(doublets)/σ\*(singlets)
>    ∈ [0.61, 0.86]** (Love–Weber bracket, μ_p = 0.3). And, unlike polydisperse spheres,
>    **there is no diameter-convention ambiguity in this onset** — the primary IS the size.
> 4. **Isostatic counting fixes the contact number**: rough (μ → ∞) clusters of ANY shape
>    jam at Z = 4 (measured 4.04 singlets, 3.88 doublets); frictionless dumbbells are
>    hypostatic (7.92 against an isostatic 10).
> 5. For a drop, doublets push the same way on **every** channel: wider window (jam), lower
>    onset stress (earlier activation), larger thickening ratio outside the window (K5), and a
>    smaller inertial-dilatant number through a lower critical packing at unchanged
>    permeability. **The singlet/doublet ratio is the one particle property that moves the
>    window, the onset and the dilatancy together at fixed primary size and chemistry.**

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
permeability at fixed φ is **independent of the singlet/doublet ratio**. (For a SINTERED
dumbbell at centre separation `s < d`, the removed caps lower `S` faster than `V`, so
`d_s` = 1.08 d at `s = 0.8d` and 1.12 d at `0.6d` — K1b, checked against Monte-Carlo volume.
Necks raise κ.) *This is the opposite of polydisperse spheres, where κ ∝ d₃₂² moves with CV.*

**K2 — no trapped fluid for compact small clusters.** The throat between three mutually
tangent primaries has inscribed radius `(2/√3 − 1)a = 0.155a > 0`, so the pore inside a
triangle or tetrahedron of tangent spheres connects to the outside. **Λ = 1 exactly** for
k ≤ 4 compact clusters: the effective volume fraction equals the solid fraction. The
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
each branch (φ₀ at μ 0; φ_m at μ_m). Justification for "no excess": all grains share one
primary, so no size-ratio filling (body-size ratio 2^{1/3} = 1.26, far below the ~2.4 at which
filling appears in this project's binaries). The rival — a positive excess from mixing
shapes — is registered beside it (`PREREG_KMER` §3).

- **Position.** `x(w) = φ/φ₀(w)`. φ₀ moves only 0.9 % from singlets to dumbbells, so **the
  doublet fraction barely moves x.** (Contrast CV, which moves x by ~a window.)
- **Contrast.** `χ(w) = A₀(w)/A_m(w)` — measured endpoints 0.886 → 0.788 (μ_m = 3). Window
  width `1 − χ` ×1.85.
- **Window opening (K4).** Setting `χ(w) = x` is LINEAR in `w`:
  `w* = (xA_m,s − A₀,s)/[(A₀,d − A₀,s) − x(A_m,d − A_m,s)]` (SymPy, K4a);
  0.87 / 0.54 / 0.34 / 0.15 at x = 0.80 / 0.83 / 0.85 / 0.87 (μ_m = 3). At the realistic
  μ_m = 0.3 the window only goes 0.092 → 0.138 and `w*` exists only for x ≥ 0.862.
- **Required activation** inside the window: `f_c(x, w) = (1 − x)/(1 − χ(w))` FALLS as `w`
  rises — a doublet-rich suspension needs a smaller frictional fraction to jam.
- **Onset (K6).** At fixed p and φ the mean contact force is `⟨f⟩ ∝ p/(n_c ℓ)` (Love–Weber:
  n_c contacts per volume, ℓ the mean contact lever arm about the body centroid). External
  contacts per primary: 4.63 (singlets) vs 2.82 (dumbbells) at μ_p 0.3 → `n_c` ratio 0.61;
  the lever arm of a lobe contact is between 1 and √2 times a sphere's. So
  `σ*(w=1)/σ*(0) = (n_c ℓ)_d/(n_c ℓ)_s ∈ [0.61, 0.86]` — the onset stress FALLS with doublets.
- **Sharpness.** δ(w): the frictional dumbbell records already show a force CV of 1.3–8.8
  against 0.7 for spheres, so **δ(dumbbells) < δ(singlets)** is predicted; the CV value is
  large enough to need checking on per-contact data before it is believed (a few huge forces
  would do it) — which the `kmer` dumps now record.

## 3. What a drop sees — every channel, one composition

With (M6), a material point jams where the impact stress exceeds
`σ_jam(w) = σ*(w)·[ln 1/f_c(x,w)]^(−1/δ(w))`. Doublets lower `σ*`, lower `f_c` and (predicted)
lower `δ`. The first two lower `σ_jam` unconditionally; the third lowers it only for
`f_c > 1/e` (`jam_stress_2026-09-22.md` §3). **Outside** the window the drop cannot jam but
still thickens by `Θ = [(1−x)/(1−x/χ)]^n`, which RISES as χ falls (K5). And the
inertial-dilatant number `ID ∝ κ/(φ − φ_c)` FALLS: κ is unchanged (K1) while the critical
packing φ_c tracks the frictional branch, which dumbbells lower by 12 %.

So **every channel moves the same way**: toward jamming, earlier activation, stronger
thickening, stronger dilatancy — at fixed primary size, chemistry and weighed φ.

## 4. The complete knob atlas for a dense-suspension drop impact

Reduced variables: position `x = φ/φ₀`, contrast `χ = φ_m/φ₀`, onset `σ*`, sharpness `δ`,
permeability `κ` (sets dilatancy/Darcy times), and the drop groups
`We = ρU₀²D₀/γ`, `Re = ρU₀D₀/η`, Stokes `St = ρ_p γ̇ b²/η_s`, and the stress ratio
`Π = σ_imp/σ*` (gauge named per use).

| particle / process knob | x | χ | σ\* | δ | κ | status | the decisive test |
|---|---|---|---|---|---|---|---|
| packing fraction φ | ✔ | — | — | — | ✔ | definition | — |
| **size distribution CV** | ✔ (~a window per ΔCV 0.5) | ✗ (paired \|t\| 0.56) | ✔ (convention-dependent sign) | ✔ (δ −19 % at CV 0.5, μ 0.3) | ✔ d₃₂² | MEASURED | bench at fixed φ, two CVs, same d₅₀ |
| **singlet/doublet ratio (this file)** | ~0 (0.9 %) | ✔✔ (window ×1.85) | ✔ predicted −14…−39 % | ✔ predicted ↓ | ✗ (K1, exact) | endpoints MEASURED, mixtures RUNNING | bench: fractionate or grow doublets, same primary, same φ |
| **open aggregates (large k)** | ✔ via Λ ≥ 1 | ✔ (shape, unmeasured) | ✔ | ? | ✔ ↑ | modelled (`multiplet_psd`) | salt/pH-aggregated suspension vs dispersed |
| friction μ_p | — | ✔✔ (window 0.113 at 0→3) | — | ✔ (δ −0.22 at 0→0.3) | — | MEASURED | coated vs bare particles |
| shape (elongation) | ~0 | ✔✔ | ? | ? | ✔ (via S/V) | χ MEASURED | rods vs spheres, same volume |
| surface potential κ | — | — | ✔ (∝ κ) | — | — | algebra | salt at fixed particle |
| primary size a | — | — | ✔ (∝ a^(m−2), m ∈ {0,1}) | — (A4b untested) | ✔ ∝ a² | algebra + one size series | size ladder at fixed chemistry — degenerate with κ unless κ is varied |
| stiffness E | — | — | — | — | — | MEASURED scale-free | — |
| particle density ρ_p | — | — | — | — | — | via St only | density-matched vs not |
| **breakage of doublets** (§5) | — | ✔ reverses during impact | ✔ | ? | — | criterion only | sintered vs salt-bridged doublets |
| drop U₀, D₀, fluid η_s, γ | — | — | — | — | — | set Π, We, Re, the film thickness | the impact ladder |

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
| P-χ | χ(w) 0.886 → 0.788 near-linearly; `w*(0.85)` ≈ 0.34 | χ flat in w, or non-monotone | `kmer` |
| P-Z | rough Z = 4 ± 0.25 at every composition | Z drifting with w at μ_p 3 | `kmer` |
| P-σ\* | σ\*(dumbbell)/σ\*(singlet) ∈ [0.55, 0.95] at μ_p 0.3 | ratio ≥ 1 | `kmer` dumps |
| P-δ | δ(dumbbell) < δ(singlet) | δ(dumbbell) ≥ δ(singlet) | `kmer` dumps |
| P-bench | at fixed φ and primary, raising doublet volume fraction past `w*(x)` switches the drop from spreading to arrest | no switch across the whole w range at an x inside [χ_d, χ_s) | bench |
| P-κ | tangent aggregation does NOT change the Darcy time at fixed φ | a measured permeability change at fixed φ with tangent (non-sintered) doublets | Darcy cell |

## 7. Open

- H_mix is a hypothesis; the rival (excess) is live.
- Isotropic packs are not sheared packs; σ\* and δ from `kmer` carry the caveat the CV lead did.
- Triplets carry no prediction (unmeasured endpoint).
- The lab's working φ (0.44) sits below every DEM window (09-09 record); the transferable
  statements are in `x` and in RATIOS, not in raw φ.

## 8. The bench protocol, in quantities a rheometer measures (scale-free by construction)

The DEM φ axis does not transfer to the bench (the lab's φ = 0.44 sits below every DEM window;
09-09 record). What transfers are **ratios measured on one material**, because the reduced
framework is exactly covariant in `(φ, φ₀, φ_m)`:

1. **χ is measurable without knowing φ's absolute scale.** From steady flow curves at several
   φ: `φ₀` = divergence of the LOW-stress viscosity branch, `φ_m` = divergence of the
   HIGH-stress branch (both fits of `η ∝ (1 − φ/φ_J)^(−n)`, same n). `χ = φ_m/φ₀`.
   **Prediction: `χ(w)/χ(0)` falls from 1 to 0.89 (μ_m = 3 convention) or 0.95 (μ_p = 0.3)
   as the doublet volume fraction goes 0 → 1**, near-linearly (H_mix), independent of
   whatever sets the absolute φ scale of the real particles.
2. **The onset ratio.** `σ*(w)/σ*(0)` read at the thickening onset of the same flow curves:
   predicted in [0.61, 0.86] at w = 1 (⚠ from isotropic packs; `kmer` dumps test it).
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
3. **So the lab's salt effect reads as SCREENING (κ ↓, σ\* ↓ at fixed χ), not aggregation**,
   unless the aggregates break under impact (`Π_b ≫ 1`, §5), which would hide them from the
   arrested state while still having lowered the onset… no — broken aggregates would lower the
   onset less, not more. ⚠ This is a reading of two recorded numbers against a model; it is
   not a measurement of the aggregation state, which no data on disk contains.
4. **The decisive bench measurement is therefore the doublet fraction itself** (microscopy or
   analytical centrifugation) at 0 and 100 mM. If w₂ is ≈ equal, the salt series is a κ
   ladder, and it is the cleanest σ\* lever the lab already owns (it moves the onset by ×3
   without moving the window). If w₂ rises strongly with salt, the arrested-state constancy
   is a puzzle for this theory and must be reported as such.
