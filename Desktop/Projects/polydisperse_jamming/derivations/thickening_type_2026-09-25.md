# Which TYPE of shear thickening: CST, DST or shear jamming, in closed form

2026-09-25 (night). Checks: `derivations/thickening_type_2026-09-25_check.py` → `_check.txt` (SymPy identity,
brute-force flow curves at 4 points, limits). Companion: `derivations/geometric_window_2026-09-25.md`.

> **Result.** In the Wyart–Cates closure the drop solver integrates, a suspension at x = φ/φ₀ in steady shear is
> - **continuously shear thickening (CST)** for x < x_DST,
> - **discontinuously shear thickening (DST, S-shaped flow curve)** for x_DST ≤ x < χ,
> - **shear-jamming (SJ)** for χ ≤ x < 1,
>
> where χ = φ_m/φ₀ and x_DST(χ, nδ) is the smallest x at which
>
> **(T)  H(u) = nδ · x(1−χ) · u e^{−u} / [J (J − x)] > 1** for some u > 0, J = 1 − (1−χ)e^{−u}.
>
> Two exact consequences: **x_DST depends on the activation sharpness δ and the divergence exponent n only through
> the product nδ**, and **a DST band exists at all iff nδ > 1** (as x → χ, sup_u H → nδ). So every particle property
> acts on the TYPE of thickening through three numbers: x (position), χ (the SJ band's lower edge) and nδ (the DST
> band's lower edge).
> On the project's measured inputs at n = 2 (μ_p 0.3): spheres CV 0 have DST band 0.029 and SJ band 0.092;
> **CV 0.5 shrinks BOTH (0.016, 0.083) while lowering x by 9 %**; pure doublets keep DST (0.027) and widen SJ (0.137).
> At n = 5 (the cornstarch fit's exponent) the DST band is ~4× wider (0.12 for spheres).

**Literature status (citation check 2026-09-26).** The three-band CST / DST / SJ structure in φ is stated by Singh,
Mari, Denn & Morris, J. Rheol. 62, 457 (2018), doi:10.1122/1.4999237, Fig. 1 (their φ_C ≈ φ₀x_DST, φ_J^μ = φ_m).
The closed form (T), the nδ-only dependence and the "DST iff nδ > 1" criterion were **derived here** and were NOT
found in Wyart & Cates, PRL 112, 098302 (2014), doi:10.1103/PhysRevLett.112.098302 (which fixes n = 2 and δ = 1 and
gives a related large-stress bound on 1 − f) nor in Singh et al. 2018; Guy/Hermes/Poon, Royer et al., Mari et al.,
Morris 2020 and Ness et al. 2022 were NOT checked. Treat novelty as unestablished.

## 1. Frame
Closure (`wc_drop_impact.c` header, steady limit fhat = fss): f(σ) = exp[−(σ*/σ)^δ]; φ_J = f φ_m + (1 − f) φ₀;
η = η_s (1 − φ/φ_J)^−n (n = MP_EXPO; 2 in the shipped presets, ~5 fitted to cornstarch, solver header). Steady simple
shear is stress-controlled here: γ̇(σ) = σ/η(σ). DST ⇔ γ̇(σ) non-monotone (a range of σ with dγ̇/dσ < 0).
Assumptions: steady, homogeneous shear; the closure's φ_J(f) is linear in f (as in the solver); no inertia.
This is a statement about the CLOSURE, i.e. about what the drop solver can produce; whether real suspensions follow
it is the closure's own validation question (the shipped cs_fit has rms(log η) = 0.675).

## 2. Derivation (C1: SymPy)
With x = φ/φ₀, φ_J/φ₀ = J(σ) = 1 − (1−χ)f, so γ̇ = (σ/η_s)·g^n, g = 1 − x/J.
(1) ln γ̇ = ln σ + n ln(1 − x/J) + const.
(2) d ln γ̇/d ln σ = 1 + n (x/J²)(dJ/d ln σ)/g.
(3) dJ/d ln σ = −(1−χ) df/d ln σ, and df/d ln σ = f·δ(σ*/σ)^δ = δ u e^{−u} with u = (σ*/σ)^δ.
(4) ⇒ d ln γ̇/d ln σ = 1 − n δ x (1−χ) u e^{−u} / [J² g] = 1 − H(u), since J² g = J(J − x).
SymPy confirms (2)–(4) reduce to 1 − H exactly (C1). DST ⇔ max_u H > 1. H increases with x at fixed u (∂/∂x of
x/(J − x) > 0), so the DST set in x is an interval [x_DST, χ).

## 3. Checks
- **Brute force (C2).** At (χ, δ, n) = (0.908, 1.254, 2), (0.863, 1.063, 2), (0.908, 1.254, 5), (0.80, 2.0, 2) the
  directly computed flow curve is monotone at x_DST − 0.002 and non-monotone at x_DST + 0.002, all four.
- **Limits (C4).** δ → 0 or n → 0: no DST (H ∝ nδ). χ → 1: H → 0, no DST (no frictional drop to discontinue).
  x → χ: sup_u H → nδ (numerically 0.799 / 1.198 / 2.504 / 1.498 for nδ 0.8 / 1.2 / 2.508 / 1.5), so for nδ < 1 the
  suspension passes CST → SJ with no DST band. nδ scaling: x_DST(χ 0.9; δ 2, n 2) = x_DST(δ 1, n 4) = 0.829005.
- **Dimensions.** H is dimensionless (u, x, χ, n, δ all are); σ* sets only WHERE on the stress axis the S appears.

## 4. Where the particle properties enter (inputs and their provenance)
| input | spheres μ_p 0.3 CV 0 | CV 0.5 | pure doublets | source |
|---|---|---|---|---|
| χ (p̂ 3e-5) | 0.908 | 0.917 | 0.863 | `grade_kmer.txt`; CV: `cvmap` inputs (ANSWER_2026-09-23 §3b) |
| δ (sheared, Derjaguin δ₁₀/₉₀, p̂ 3e-4) | 1.254 | 1.017 (×0.81) | 1.063 | `grade_actshape.txt` H4; `grade_shearkmer.txt` |
| x shift at fixed weighed φ | — | ×0.91 (φ₀ +9.9 %) | ≈ ×1.01 (φ₀ −0.9 %) | `filling_factor_posthoc`, shapechi |
⚠ δ is the DEM contact-force distribution's Weibull exponent, identified with the closure's activation exponent —
the project's standing convention (`jam_stress_2026-09-22`), not a derivation. The clumps' δ is not measured.

Bands (C3), in x:
| | n = 2: x_DST | DST band | SJ band | n = 5: x_DST | DST band |
|---|---|---|---|---|---|
| spheres CV 0 | 0.879 | 0.029 | 0.092 | 0.788 | 0.120 |
| spheres CV 0.5 | 0.901 | 0.016 | 0.083 | 0.832 | 0.086 |
| pure doublets | 0.836 | 0.027 | 0.137 | 0.733 | 0.130 |
| clumps α 1.1–2 (χ at 1e-3, δ assumed sphere) | 0.861 | 0.033 | 0.106 | 0.761 | 0.133 |

## 5. Reading
- **Size distribution** pushes toward CST three ways at once: x falls (better filling), the SJ band narrows (χ ↑), and
  the DST band narrows (δ ↓ raises x_DST). A monodisperse suspension that is DST or SJ at a given weighed φ can be CST
  at CV 0.5. This is the closure-level statement of the founding bench observation.
- **Shape (doublets)** widens SJ and leaves the DST band about the same (its δ drop is offset by χ falling).
- **The exponent n matters as much as either** (it multiplies δ). The solver's shipped n = 2 gives a DST band of ~0.03
  in x; the cornstarch fit's n ≈ 5 gives ~0.12. Quoting any DST boundary without n is meaningless (the 2026-09-22
  lesson on Θ, again).
- **Surface potential and mean size do not appear**: they set σ* (where on the stress axis), not the type. They decide
  whether an impact's stress reaches the S or the jam at all (`jam_stress_2026-09-22` §(8)).

## 6. What the DROP sees (registered DNS test, `PREREG_DSTTYPE_2026-09-26` RESULT)
The steady-shear type above is a rheometer classification. In the drop solver (same closure, n = 2, L9) it has **no
signature**: 20 drops straddling x_DST on two input sets show peak spread strictly monotone in the thickening ratio
Θ = [(1−x)/(1−x/χ)]^n, with no jump at x_DST (D2 FAIL 0/4) and matched-Θ drops of different type equal within the
same-type control (D1 THETA). With the cvmap result that activation (σ*, δ) carries < 0.05 % of the CV effect, the
drop's outcome in this closure depends on the particles ONLY through (x, χ) and n: **shear-jams if x ≥ χ; otherwise
stiffens by Θ(x, χ)**. Example: spheres at a weighed φ giving x = 0.90 (χ 0.9106) have Θ = 74; the same weighed φ as
a CV-0.5 population (x 0.819, χ 0.917) has Θ = 2.9 — a ×26 weaker thickening, which is the founding observation's
extra spreading. Scope: one solver, one resolution, n = 2; the rheometer CST/DST distinction is real in the closure,
it just does not reach the drop's spread.
