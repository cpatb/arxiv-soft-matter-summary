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
| **convex facets** (rounded cubes, level-set DEM) | **the window NARROWS**: W 0.07–0.10 (sphere) → 0.04–0.05 (ε 0.5), because φ₀ falls (−0.027) while φ_m holds — flat faces resist tipping without friction, so the frictionless branch jams 'frictional-like'. ⚠ The earlier stand-in (spheres with a contact moment limit, which lower φ_m by up to 0.035 at μ_p 0.3) is NOT what convex facets do | — | κ × ψ² — minor | pilot (N 100, 2–3 seeds, one rung; §4) |
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
  window (×1.4–1.5 at the jamming rung, already at 10 % elongation). Features that let a frictionless contact pin
  rotation by itself (flat faces) narrowed it in the pilot's 2/2 available seeds (convex rounded cubes; the registered
  replication is incomplete). Concavity additionally lowers the
  jamming onset φ_m; mild convex elongation raises it. Spheres sit between the two families, not at an extreme.
  Shape barely touches drainage.
- **Surface potential.** Does not decide the type; decides the onset stress σ* and, through the double layer's
  softness, whether the jammed drop is a weak undrained ball (silica) or a stiff dilatant hybrid (cornstarch).
- **Across the two lab materials:** Gate A does NOT clearly favour either on shape alone — cornstarch is both
  faceted (narrows) and mildly elongated/irregular (widens); the published χ 0.76 vs 0.83 says wider but mixes lab,
  size, friction and φ convention;
  Gate C separates them by contact mechanics first (K_p ×1870), size second (κ ×190).

## 4. Convex faceted grains (level-set superellipsoids) — PILOT (`PREREG_LSPILOT_2026-09-26`, `out/grade_lspilot.txt`)

The clump ladder cannot speak for cornstarch, because a two-sphere clump has a concave neck. A superellipsoid with
equal axes goes continuously from a sphere (ε = 1) to a rounded cube (ε → 0) while staying convex. Level-set DEM, N 100,
p̂ 1e-3, 1000 surface nodes. The method control reproduces the FrictMat sphere to 0.0003 (μ_p 0.3) and 0.008 (μ_p 0),
after a 102-node version failed it by +0.02 (node under-resolution of penetration, diagnosed and fixed).

| ε | seeds | φ₀ | φ_m | window W = 1 − φ_m/φ₀ |
|---|---|---|---|---|
| 1.0 (sphere) | 3 | 0.6416 / 0.6374 / 0.6370 | 0.5797 / 0.5756 / 0.5927 | 0.097 / 0.097 / 0.070 |
| 0.7 | 3 | 0.6297 / 0.6399 / 0.6293 | 0.5781 / 0.5808 / 0.5775 | 0.082 / 0.092 / 0.082 |
| 0.5 (rounded cube) | 2 | 0.6105 / 0.6124 | 0.5871 / 0.5826 | 0.038 / 0.049 |

**In the pilot's two available seeds, convex faceting narrows the window** — the opposite of elongation or concavity;
the registered replication (3 seeds + N 200, Amendment 3) and the ε-0.5 node-convergence check (Amendment 4) decide
whether it stands (status in `out/grade_lspilot.txt`). The node check came back NODE-LIMITED (φ₀(ε 0.5) = 0.6105 at
1000 nodes, 0.6197 at 2000; bar 0.005): the sign survives — every ε-0.5 value at ≥ 400 nodes is ≥ 0.017 below the
lowest sphere seed — but the size of the narrowing is uncertain by about a third. The frictionless branch falls
(φ₀ −0.027 at ε 0.5) while the frictional branch holds. Counting gives a geometric reason: a flat face resists tipping
even without friction, because the contact point jumps to the face's edge, as for a box on a table. Face contacts
therefore supply the rotational constraints that only friction supplies to spheres, and a faceted pack jams
"frictional-like" already at μ_p 0 (Z₀ ≈ 6.4, far below the Z_iso = 12 point contacts on a general body would need).
ε 0.7 (mild faceting) is within the seed scatter of the sphere. **Status: the registered replication (3 seeds + N 200)
is [see grade_lspilot.txt: REVERSAL line]; one rung, N 100–200; the face-contact mechanism is a counting reading,
not a measured contact census.**

**What this does to the cornstarch expectation.** Two DEM shape families now point opposite ways. Elongated and
concave grains (clumps, doublets) widen the window ×1.4–1.5; convex faceted grains narrow it. Cornstarch granules
are faceted polyhedra with modest elongation, so geometry alone does not decide whether its window is wider or
narrower than Stöber silica's. The published same-paper fits (cornstarch χ 0.76 vs silica spheres 0.83) say wider,
but they differ in lab, size, friction and φ convention, and cannot isolate shape.

