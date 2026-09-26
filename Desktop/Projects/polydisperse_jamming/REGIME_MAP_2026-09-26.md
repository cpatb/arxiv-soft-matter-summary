# What decides the TYPE of shear thickening an impacting drop shows: size, shape, size distribution, potential

2026-09-26 (overnight). The author's question: *how do size, shape, size distribution/variance and surface potential
each affect the type of shear thickening seen in suspensions in drop impact?* This page is the map; the derivations
and campaigns behind every row are linked. DEM = this project's Yade packs (dry, frictional, isotropic or sheared);
DNS = the Wyart–Cates drop solver `wc_drop_impact.c`; scaling = order-of-magnitude estimates on the project's
poroelastic model. **No new bench data.** Status per row: measured / derived / DNS / scaling / open.

## 1. Three gates, in the order an impact meets them

**Gate A — which TYPE the material can show at its packing fraction** (`derivations/thickening_type_2026-09-25.md`).
With x = φ/φ₀, χ = φ_m/φ₀ and the activation sharpness δ, divergence exponent n:

| x range | steady-shear type | drop impact, if the impact stress reaches the thickened branch |
|---|---|---|
| x < x_DST(χ, nδ) | CST — continuous thickening, bounded (Θ = [(1−x)/(1−x/χ)]^n) | stiffens, spreads, retracts |
| x_DST ≤ x < χ | DST — S-shaped flow curve, discontinuous jump | stiffens abruptly; still cannot jam |
| χ ≤ x < 1 | SJ — shear jamming | a jammed solid region forms (arrests within t/τ ≈ 0.1 in the DNS) |
| x ≥ 1 | jammed at rest | a paste |

A DST band exists only if nδ > 1. x_DST depends on n and δ only through nδ.

**Gate B — whether the impact actually gets there.** A material point in the SJ band jams where the local stress
exceeds σ_jam = σ*[ln(1/f_c)]^(−1/δ), f_c = (1 − x)/(1 − χ) (`derivations/jam_stress_2026-09-22.md`,
stress-controlled; the DNS drop behaves stress-controlled, `jamwin` P2). σ* ∝ F_rep/a² (onset stress) carries size
and surface potential. For the lab silica the wall-shear activation gauge is G_act = τ_wall/σ* ≈ 4–70 across the
literature σ* bracket (`ball_hybrid_scaling`), i.e. the impact is far past onset.

**Gate C — what KIND of solid the jammed drop is: ball or hybrid** (`derivations/ball_hybrid_scaling_2026-09-25.md`,
scaling). Drainage Deborah number De = τ_pore/τ_impact and dilatant-strength ratio S = K_p Δφ/(φ ρU₀²):
De ≫ 1 with S ≪ 1 is undrained with a weak contact network (the 800 nm silica: De ≈ 7×10⁵, S ≈ 0.004 — the ball);
De ~ 1–30 with S ~ 1–20 has drainage, dilatancy and inertia all O(1) (10 µm cornstarch: De ≈ 7, S ≈ 4 — the hybrid).

## 2. The knobs

Arrows: ↑/↓ = raises/lowers the quantity; "→ CST" = pushes the material toward continuous thickening (away from
jamming), "→ SJ" = toward shear jamming, "→ ball"/"→ hybrid" = Gate C.

| knob | Gate A: x, χ, nδ (type) | Gate B: σ* (onset stress) | Gate C: De, S (ball vs hybrid) | status |
|---|---|---|---|---|
| **mean size a** | none: the dry DEM is scale-free in size (only ratios transfer) | σ* ∝ a^(m−2): a⁻² at fixed repulsive force, a⁻¹ Derjaguin → bigger grains thicken at lower stress | κ ∝ a² → De ↓ (∝ a⁻² with Hertz contacts) → hybrid; S: size-free (Hertz) or ↓ (EDL) | algebra + scaling |
| **size variance (CV)** | better filling: w 1.00 → 0.76 (lognormal CV 0.5) ⇒ φ₀ +9.9 %, x ×0.91; χ 0.911 → 0.917; δ ×0.81 ⇒ x_DST ↑. **All three → CST** | dry DEM: ×0.73–0.79 at fixed d₁₀ via Zφ; lubricated LF_DEM: 0.96–1.17 — **contested** | κ via d₃₂: ×(1+CV²)⁴ ⇒ De ×0.41 at CV 0.5 → toward drained | measured (A); DNS: at fixed weighed φ CV 0 jams, CV 0.5 cannot and spreads 35–67 % further (window carries 99.5–101.8 %) |
| **distribution SHAPE** (family at fixed CV) | w at CV 0.5: lognormal 0.76, uniform 0.91 — the family moves x as much as CV does; two binaries matched in CV pack 0.0635 apart (filling functional) | via Zφ | via d₃₂ | measured |
| **mild asphericity** (convex-ish, α ≤ 1.5) | window opens: 1 − χ 0.072 → ~0.10 (clumps, p̂ 1e-3); φ₀ AND φ_m rise (φ₀ ≈ 0.70 near α 1.2–1.5) ⇒ at fixed weighed φ, x falls. Net: wider SJ band, starting higher. [convex faceted grains: level-set campaign, §4] | not measured | κ × ψ² (≥ 0.64) — minor | measured (clumps; neck-confounded at α 1.1) |
| **concavity / interlocking** (necks, doublets, aggregates) | φ_m falls continuously as a neck closes (w ↑ 1.12 → 1.38); pure doublets window ×1.49 at the jamming rung; δ ×0.85 sheared → **→ SJ at lower φ** | ×0.78 (doublets, every contact primary–primary) | — | measured; DNS: doublets jam where singlets cannot (`kmermap` R1), entirely through χ |
| **facets / contact lever arms** | (spheres with a contact moment limit as the stand-in) φ_m −0.028 … −0.035 at μ_p 0.3, half by η ≈ 0.04–0.05, saturated by ≈ 0.3; by removing contacts, not by filling; leaves the CV effect intact (R7 ρ 0.97) | — | — | measured on the stand-in only; polyhedra [§4] |
| **friction μ_p** (roughness, surface chemistry) | window 0.09 (μ 0.3) → 0.11 (μ 3) at CV 0; δ −0.22 from μ 0 → 0.3 | — | — | measured |
| **surface potential** (charge, Debye length, salt) | none directly; indirectly if it aggregates (→ concavity row) or lubricates/roughens contacts (→ friction row) | σ* ∝ F_rep: the lab salt series (onset ×0.31 at 100 mM, arrested state ×1.28 of ×6.35) reads as screening | **K_p**: a double-layer network is soft (silica K_p ≈ 590 Pa) → S ≪ 1 → ball; bare Hertzian contacts (cornstarch K_p ≈ 1.1×10⁶ Pa) → S ~ 4 → hybrid. The largest single lever between the two lab materials (×1870) | scaling; bench reading |
| **fluid viscosity** μ_f | η_s only | — | De ∝ μ_f (×10–1000 with glycerol–water); S unchanged | scaling |
| **impact** U₀, D₀ | — | stress scale ρU₀², τ_wall | τ_impact = D₀/U₀; S ∝ U₀⁻² | — |

## 3. Reading the map, per question
- **Size (mean).** Does not change what the packing CAN do (Gate A is scale-free); it shifts WHEN (σ*) and HOW the
  jammed state drains (De ∝ a⁻²). Bigger grains: easier onset, more drained → toward the hybrid.
- **Size variance.** Acts on filling. Every channel it touches (x, χ, x_DST) pushes toward CST; it is the reason a
  polydisperse suspension at the same weighed φ thickens and spreads where a monodisperse one jams.
- **Shape.** Acts on counting and filling. Spheres are the one shape with free frictionless rotation, so they have
  the narrowest window; any asphericity widens it. Concavity (necks, aggregates) additionally lowers the jamming
  onset; mild convex asphericity raises it. Shape barely touches drainage.
- **Surface potential.** Does not decide the type; decides the onset stress σ* and, through the double layer's
  softness, whether the jammed drop is a weak undrained ball (silica) or a stiff dilatant hybrid (cornstarch).
- **Across the two lab materials:** Gate A favours cornstarch (non-spherical, wider window; published χ 0.76 vs 0.83);
  Gate C separates them by contact mechanics first (K_p ×1870), size second (κ ×190).

## 4. Convex faceted grains (level-set superellipsoids) — [in progress]

