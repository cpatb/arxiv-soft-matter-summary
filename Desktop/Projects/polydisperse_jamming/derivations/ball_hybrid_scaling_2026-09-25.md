# Ball → inertia-dilatant hybrid: which dimensionless group can carry the size and CV trend?

2026-09-25. Numbers: `derivations/ball_hybrid_scaling.py` → `out/ball_hybrid_scaling.txt` (every number below is
printed there; the physics modules are `impact_stress/{pore_relaxation, inertial_dilatant, stress}.py`, reused,
not retyped). Companion registration: `PREREG_DILSIZE_2026-09-25.md`.

> **Result (a scaling estimate, not a measurement).**
> 1. Four groups could separate a solid ball from an inertia-dilatant hybrid: shear activation G_act = τ_wall/σ\*,
>    drainage Deborah number De = τ_p/τ_imp, dilatant-strength ratio S = p_dil,max/ρU₀², and Mobaseri et al.'s
>    ID. The two lab materials sit in **different corners of (De, S)**: 800 nm silica is deeply undrained
>    (De ≈ 7×10⁵) with a negligible dilatant strength (S ≈ 0.004), while 10 µm cornstarch is near the
>    drained/undrained crossover (De ≈ 7, bracket 1.6–34) with S ≈ 4 (0.9–19), i.e. **drainage, dilatancy and
>    inertia are all O(1) only for cornstarch** — which is what "inertia-dilatant hybrid" means.
> 2. **Size alone at fixed silica chemistry does not get there.** From 0.4 to 20 µm, De stays ≥ 10³ and S ≤ 0.007
>    under both σ\* conventions. The silica → cornstarch jump in (De, S) is dominated by the network modulus K_p
>    (double-layer 590 Pa vs Hertzian 1.1×10⁶ Pa, ×1870), with permeability (κ ×190) second. So on this model
>    the bench trend "larger size → hybrid" is really "larger size AND stiffer contact mechanics".
> 3. **The WC solver cannot represent the strength route at all.** `wc_drop_impact.c` adds the dilatant pressure
>    only to the stress that RECRUITS frictional contacts (`sigma += p_dil`); nothing makes the jammed state
>    stronger or weaker with drainage. Where shear activation is already saturated (silica: G_act = 4–70), the
>    dilatancy term is expected to be inert, and every solver-representable size handle (σ\* ∝ a⁻¹ or a⁻²) pushes
>    toward MORE activation — the opposite of the bench trend. That is a registered negative prediction
>    (PREREG_DILSIZE P1–P2), and if it holds, the transition lives outside the current closure.
> 4. CV enters three ways, only one of them tested: the window (x falls, `cvmap`, shown), permeability
>    (κ ∝ d₃₂² → ×2.44 at CV 0.5, De ×0.41, ID ×2.44), and σ\* (×0.735 at fixed d₁₀, which pushes toward jamming).

## 1. The groups and where each comes from
- **G_act = ⟨τ_wall⟩(t\*=1)/σ\*** — the shear-stress activation gauge (`dilatancy.lambda_shear`; wall shear from
  `stress.wall_shear_mean`, a Philippi/Wagner field). The project's memory records that ρU₀² is the wrong gauge by
  ~2 decades; this one is used instead. Needs η₀ of the unthickened suspension; the lab-silica value is from a
  retracted fit, so η₀ is a declared bracket 0.01–0.1 Pa s.
- **De = τ_p(L)/τ_imp**, τ_p = L²/D_p, D_p = κ/(μ S_stor), S_stor = (1−φ)χ_f + 1/K_p (`pore_relaxation`). L = the
  contact-patch radius 0.49 D₀, the length `PORE_PHYSICS.md` §3 chose for exactly this question. De ≫ 1: undrained,
  dilation frustrated; De ≪ 1: drained, no dilatant stress.
- **S = K_pΔφ/(φ ρU₀²)** — the largest dilatant pressure the network can carry (`dilatant_pressure_ceiling`,
  validated against Jerome et al.'s measured suctions in `test_pore_relaxation.py`) over the impact stress.
- **ID** — [M26] Eq. 6, printed for reference only: the project already recorded that Eq. 6 does not evaluate to 1
  at the boundary it is drawn at (memory: inertial-dilatant-number), and it contains neither K_p nor σ\*.

## 2. Numbers (lab operating points: silica D₀ 2.5 mm, U₀ 2.3 m/s; cornstarch [M26] D₀ 4 mm, U₀ 5.1 m/s)
| | silica φ 0.44 | silica φ 0.50 | cornstarch φ 0.42 |
|---|---|---|---|
| K_p (bracket) | 590 Pa (197–5.9e3) | 590 Pa | 1.1e6 Pa (2.4e5–5.1e6) |
| De | 7.3e5 (7.3e4–2.2e6) | 1.3e6 | **7.2 (1.6–34)** |
| S | 0.0035 (0.0012–0.035) | 0.012 | **4.1 (0.9–19)** |
| ID | 5.5e-5 | 7.1e-6 | 4.7e-3 |
| ⟨τ_wall⟩, η₀ 0.01–0.1 Pa s | 425–1340 Pa | 434–1370 Pa | 1030–3240 Pa |
G_act for silica (σ\* 19 / 44 / 100 Pa, literature bracket): 22–70 / 10–30 / 4–13. No measured cornstarch σ\* is on
disk, so its G_act is not given. With the network's confining pressure set to each σ\* of that bracket instead of
the retracted 44 Pa fit, silica stays deeply undrained: De 1.7×10⁶ / 7.3×10⁵ / 3.2×10⁵ and S 0.0015 / 0.0035 /
0.008. ⚠ `silica_poro` sets the confining pressure equal to σ\*, which equates a SHEAR onset stress with a particle
PRESSURE; they differ by an O(1) friction coefficient, which moves K_p, De and S by that factor, not by decades.

## 3. Size ladder at fixed silica chemistry (anchor 0.8 µm, φ 0.44), both σ\* conventions
Exponents (d ln/d ln a): κ +2. [F] fixed repulsive force: σ\* −2, K_p (= σ\* a/6λ_D, `silica_poro`'s rule) −1,
De −1, S −1. [D] Derjaguin (force ∝ a): σ\* −1, K_p 0, De −2, S 0. Over 0.4–20 µm: De 1.5e6 → 2.9e4 [F] or
2.9e6 → 1.2e3 [D]; S 0.007 → 0.0001 [F] or constant 0.0035 [D]. **Neither convention reaches De ~ 1, S ~ 1.**
⚠ The K_p rule ties the network modulus to σ\* (P_conf = σ\*); that is `silica_poro`'s declared choice, not a
measurement, and it is the most model-laden line here.

## 4. What this means for "why size and CV move ball → hybrid"
- The candidate that fits the two materials is the (De, S) corner: silica is undrained with a weak network, cornstarch
  sits where drainage, dilatant strength and inertia all compete. But the move between them is mostly contact
  mechanics (K_p), so a size series at fixed chemistry is predicted to stay ball-like. **If the bench shows a
  size-driven transition at fixed chemistry, this model is missing it.**
- Inside the WC solver, the only size handles are σ\* and the dilatant ACTIVATION term, both of which push toward
  jamming as size grows. The solver therefore has no mechanism for the size half of the observation; it does
  reproduce the CV half through the window.

## 5. Open, and what would decide it
- K_p for the lab silica (unmeasured; EDL modulus depends on λ_D — the salt series is a K_p series).
- Whether the bench size trend is at fixed chemistry (it compares silica sizes? or silica vs cornstarch?). The
  answer decides whether (De, S) or something else carries it.
- A closure term for drainage-limited STRENGTH of the jammed state (not just activation): the minimal form is a
  cap on the jammed-state stress at min(p_dil capacity-limited, Darcy) × friction coefficient. Not implemented.
