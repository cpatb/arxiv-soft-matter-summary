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
| x < x_DST(χ, nδ) | CST — continuous thickening, bounded (Θ = [(1−x)/(1−x/χ)]^n) | stiffens by Θ, spreads, retracts |
| x_DST ≤ x < χ | DST — S-shaped flow curve, discontinuous jump | **the same as CST in the drop**: stiffens by Θ, cannot jam |
| χ ≤ x < 1 | SJ — shear jamming | a jammed solid region forms (arrests within t/τ ≈ 0.1 in the DNS) |
| x ≥ 1 | jammed at rest | a paste |

A DST band exists only if nδ > 1; x_DST depends on n and δ only through nδ. **But the drop does not see the CST/DST
line** (`PREREG_DSTTYPE_2026-09-26` RESULT: 20 drops, peak spread strictly monotone in Θ, no jump at x_DST — D4, the cleanest evidence; D1 THETA; the confound-free within-set D1′ is UNDETERMINED, the DST drops sitting slightly ABOVE the CST trend, i.e. saturation, never a type drop).
So for drop impact the classification is two-way: **jammable (x ≥ χ) or stiffening by Θ(x, χ)**, and with the cvmap
finding that, ON THE CV AXIS, activation (σ*, δ) carries < 0.05 % of the effect, the particles act on a drop that is
past onset through x and χ (at the fixed n = 2 of every DNS cell here — n was never varied; one closure, L9). Whether
friction, shape and potential also collapse onto (x, χ) in the drop was not tested directly.
Spheres at x 0.90 (χ 0.911): Θ = 74; the same weighed φ at CV 0.5 (x 0.819, χ 0.917): Θ = 2.9.

**Gate B — whether the impact actually gets there.** A material point in the SJ band jams where the local stress
exceeds σ_jam = σ*[ln(1/f_c)]^(−1/δ), f_c = (1 − x)/(1 − χ) (`derivations/jam_stress_2026-09-22.md`,
stress-controlled; the DNS drop behaves stress-controlled, `jamwin` P2). σ* ∝ F_rep/a² (onset stress) carries size
and surface potential. For the lab silica the wall-shear activation gauge is G_act = τ_wall/σ* ≈ 4–70 across the
literature σ* bracket (`ball_hybrid_scaling`), i.e. the impact is far past onset.

**Gate C — what KIND of solid the jammed drop is: ball or hybrid** (`derivations/ball_hybrid_scaling_2026-09-25.md`,
scaling). Drainage Deborah number De = τ_pore/τ_impact and dilatant-strength ratio S = K_p Δφ/(φ ρU₀²):
De ≫ 1 with S ≪ 1 is undrained with a weak contact network (the 800 nm silica: De ≈ 7×10⁵, S ≈ 0.004 — the ball);
De ~ 1–30 with S ~ 1–20 has drainage, dilatancy and inertia all O(1) (10 µm cornstarch: De ≈ 7, S ≈ 4 — the hybrid).

⚠ **Scope, from the final registered grade (`PREREG_ROLLCHI` R5 FAIL, 2026-09-26 19:45):** with w = 1 the relation
holds on 83 % of 156 CV-0 sphere cells (registered bar 90 %); it fails where friction and rolling resistance push Z
below 4 — residual +0.018 at μ_p 1, +0.014 at μ_p 3, +0.013–0.014 at rolling η ≥ 0.3, up to +0.08 — i.e. outside the
[4, 6] range SWM state it for. The per-(CV, rung) separability test predicts within 0.015 on only 62–92 % of cells.
So (G) with one filling factor is a good description at REALISTIC friction (μ_p ≤ 0.3, no rolling resistance),
not a general law.
## 2. The knobs

Arrows: ↑/↓ = raises/lowers the quantity; "→ CST" = pushes the material toward continuous thickening (away from
jamming), "→ SJ" = toward shear jamming, "→ ball"/"→ hybrid" = Gate C.

| knob | Gate A: x, χ, nδ (type) | Gate B: σ* (onset stress) | Gate C: De, S (ball vs hybrid) | status |
|---|---|---|---|---|
| **mean size a** | none: the dry DEM is scale-free in size (only ratios transfer) | σ* ∝ a^(m−2): a⁻² at fixed repulsive force, a⁻¹ Derjaguin → bigger grains thicken at lower stress | κ ∝ a² → De ↓ (∝ a⁻² with Hertz contacts) → hybrid; S: size-free (Hertz) or ↓ (EDL) | algebra + scaling |
| **size variance (CV)** | better filling: w 1.00 → 0.76 (lognormal CV 0.5) ⇒ φ₀ +9.9 %, x ×0.91; χ 0.911 → 0.917; δ ×0.81 ⇒ x_DST ↑. **All three → CST** | dry DEM: ×0.73–0.79 at fixed d₁₀ via Zφ; lubricated LF_DEM: 0.96–1.17 — **contested** | κ via d₃₂: ×(1+CV²)⁴ ⇒ De ×0.41 at CV 0.5 → toward drained | measured (A); DNS: at fixed weighed φ CV 0 jams, CV 0.5 cannot and spreads 35–67 % further (window carries 99.5–101.8 %) |
| **distribution SHAPE** (family at fixed CV) | w at CV 0.5: lognormal 0.76, uniform 0.91 — the family moves x as much as CV does; two binaries matched in CV pack 0.0635 apart (filling functional) | via Zφ | via d₃₂ | measured |
| **mild asphericity** (convex-ish, α ≤ 1.5) | window opens: 1 − χ 0.092 → 0.131 at α 1.1 (clumps, jamming rung p̂ 3e-5, n 3); φ₀ AND φ_m rise (φ₀ ≈ 0.70 near α 1.2–1.5 at p̂ 1e-3, where α 1.5's frictionless value and α 1.2–1.3's frictional values are n = 1–2) ⇒ at fixed weighed φ, x falls. Net: wider SJ band, starting higher. [convex faceted grains: level-set campaign, §4] | not measured | κ × ψ² (≥ 0.64) — minor | measured (clumps; neck-confounded at α 1.1) |
| **concavity / interlocking** (necks, doublets, aggregates) | φ_m falls continuously as a neck closes (w ↑ 1.12 → 1.38); pure doublets window ×1.49 at the jamming rung; δ ×0.85 sheared → **→ SJ at lower φ** | ×0.78 (doublets, every contact primary–primary) | — | measured; DNS: doublets jam where singlets cannot (`kmermap` R1), entirely through χ |
| **convex facets** (rounded cubes, level-set DEM) | **φ₀ falls 0.026 (resolved, t 17)** — flat faces resist tipping without friction, so the frictionless branch jams at lower density; φ_m and hence the window's change are NOT resolved at N 100 (W narrower on 2/4 seeds). ⚠ The earlier stand-in (spheres with a contact moment limit, which lower φ_m by up to 0.035 at μ_p 0.3) is NOT what convex facets do | — | κ × ψ² — minor | pilot (N 100, 2–3 seeds, one rung; §4) |
| **friction μ_p** (roughness, surface chemistry) | window 0.09 (μ 0.3) → 0.11 (μ 3) at CV 0; δ −0.22 from μ 0 → 0.3 | — | — | measured |
| **surface potential** (charge, Debye length, salt) | none directly; indirectly if it aggregates (→ concavity row) or lubricates/roughens contacts (→ friction row) | σ* ∝ F_rep: the lab salt series (onset ×0.31 at 100 mM, arrested state ×1.28 of ×6.35) reads as screening | **K_p**: a double-layer network is soft (silica K_p ≈ 590 Pa) → S ≪ 1 → ball; bare Hertzian contacts (cornstarch K_p ≈ 1.1×10⁶ Pa) → S ~ 4 → hybrid. The largest single lever between the two lab materials (×1870). Across the lab salt range (λ_D 30 → 1 nm) silica's K_p rises ×10–31 but S stays 0.001–0.04 and De ≥ 7×10⁴ (`derivations/potential_ladder_2026-09-26.txt`): salt moves the onset, not the ball | scaling; bench reading |
| **fluid viscosity** μ_f | η_s only | — | De ∝ μ_f (×10–1000 with glycerol–water); S unchanged | scaling |
| **impact** U₀, D₀ | — | stress scale ρU₀², τ_wall | τ_impact = D₀/U₀; S ∝ U₀⁻² | — |

## 3. Reading the map, per question
- **Size (mean).** Does not change what the packing CAN do (Gate A is scale-free); it shifts WHEN (σ*) and HOW the
  jammed state drains (De ∝ a⁻²). Bigger grains: easier onset, more drained → toward the hybrid.
- **Size variance.** Acts on filling. Every channel it touches (x, χ, x_DST) pushes toward CST; it is the reason a
  polydisperse suspension at the same weighed φ thickens and spreads where a monodisperse one jams.
- **Shape.** Acts on counting and filling, and its SIGN depends on which geometric feature it adds. Features that
  add rotational freedom a frictionless contact must pin (elongation, a concave neck: clumps, doublets) WIDEN the
  window (×1.4–1.7 at the jamming rung, already at 10 % elongation). Features that let a frictionless contact pin
  rotation by itself (flat faces) lower the frictionless edge φ₀ (convex rounded cubes, resolved); whether the window
  narrows is not resolved (4 seeds, N 100). Concavity additionally lowers the
  jamming onset φ_m; mild convex elongation raises it. Spheres sit between the two families, not at an extreme.
  Shape barely touches drainage.
- **Surface potential.** Does not decide the type; decides the onset stress σ* and, through the double layer's
  softness, whether the jammed drop is a weak undrained ball (silica) or a stiff dilatant hybrid (cornstarch).
- **Across the two lab materials:** Gate A does NOT clearly favour either on shape alone — cornstarch is both
  faceted (lowers φ₀; window change unresolved) and mildly elongated/irregular (widens); the published χ 0.76 vs 0.83 says wider but mixes lab,
  size, friction and φ convention;
  Gate C separates them by contact mechanics first (K_p ×1870), size second (κ ×190).

## 4. Convex faceted grains (level-set superellipsoids) — PILOT (`PREREG_LSPILOT_2026-09-26`, `out/grade_lspilot.txt`)

The clump ladder cannot speak for cornstarch, because a two-sphere clump has a concave neck. A superellipsoid with
equal axes goes continuously from a sphere (ε = 1) to a rounded cube (ε → 0) while staying convex. Level-set DEM, N 100,
p̂ 1e-3, 1000 surface nodes. The method control reproduces the FrictMat sphere to 0.0003 (μ_p 0.3) and 0.008 (μ_p 0),
after a 102-node version failed it by +0.02 (node under-resolution of penetration, diagnosed and fixed).

Seed means (sphere 5 seeds, rounded cube ε 0.5 4–5 seeds; Welch t):
| branch | sphere | ε 0.5 | difference |
|---|---|---|---|
| frictionless φ₀ | 0.6370 ± 0.0031 | 0.6108 ± 0.0011 | **0.026, t 17.5 — resolved** |
| frictional φ_m (μ_p 0.3) | 0.5914 ± 0.0147 | 0.5814 ± 0.0062 | 0.010, t 1.4 — NOT resolved |
| window W = 1 − φ_m/φ₀ (seed-paired, 4 pairs) | | | +0.022 ± 0.018; narrower on 2/4 seeds — NOT resolved |

- **What is established: convex faceting lowers the FRICTIONLESS jamming point**, by 0.026 (node-limited: at 2000
  surface nodes one ε-0.5 seed gives 0.6197 against 0.6105 at 1000, so ±0.009 on this number; Amendment 4 NODE-LIMITED).
  Counting gives a geometric reading: a flat face resists tipping even without friction, because the contact point jumps
  to the face's edge, as for a box on a table. Face contacts then supply rotational constraints that only friction
  supplies to spheres, so a faceted pack jams at a lower density already at μ_p 0 (Z₀ ≈ 6.4, far below the Z_iso = 12
  that point contacts on a general body would need). A reading, not a measured contact census.
- **What is NOT established: whether the shear-jamming WINDOW narrows.** That needs φ_m, and the sphere's frictional
  branch scatters by 0.015 between seeds at N 100 (0.576–0.613). The first two seeds said "narrows by half"; seeds 4–5
  did not. The registered 3-seed + N-200 reversal test is INCOMPLETE (one cell never landed); the N-200 pair points the
  narrowing way (0.084 vs 0.044). ε 0.7 is within scatter of the sphere on every branch.

**What this does to the cornstarch expectation.** Elongated and concave grains (clumps, doublets) widen the window
×1.4–1.7 at the jamming rung (replicated). Convex faceting moves the frictionless edge down by a resolved 0.026. It
therefore pushes the other way if φ_m holds, but that is unresolved. Cornstarch granules are faceted polyhedra with
modest irregularity, so on geometry alone this project cannot say whether their window is wider or narrower than
Stöber silica's. The published same-paper fits (cornstarch χ 0.76 vs 1.5 µm silica spheres 0.83) say wider, but they
differ in lab, size, friction and φ convention. What the pilot does predict robustly is a LOWER frictionless jamming
fraction φ₀ for faceted grains — measurable on a rheometer as the high-shear-rate (frictionless-branch) divergence.
