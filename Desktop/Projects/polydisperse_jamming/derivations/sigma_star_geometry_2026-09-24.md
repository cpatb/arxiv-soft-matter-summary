# σ\* from contact geometry: why polydispersity lowers the onset stress

2026-09-24. Checks: `derivations/sigma_star_geometry_check.py` → `out/sigma_star_geometry.txt`
(every number below is printed there).

> **Result.** For a Derjaguin-scaled threshold (`F_rep = κR*`), the Wyart–Cates onset stress of a
> sheared polydisperse packing is, EXACTLY,
>
> **σ\*/κ = (N_c/V)·⟨aᵢaⱼ⟩_c · C / (3 q Λ)** (eq. 7)
>
> and, with a pair-statistics closure that holds to ≤ 2 % (μ_p 0) and ≤ 8 % (μ_p 0.3) on 52 sheared DEM
> cells, **σ\*/κ ≈ Z φ / (8π⟨a⟩) · C/q**. At fixed NUMBER-mean radius, polydispersity therefore
> acts on σ\* through the contact density Zφ alone (plus two O(5 %) corrections). CV 0 → 0.5 cuts σ\*
> by ×0.73–0.79; the geometric factor Zφ alone is ×0.60–0.71, and the force-size correlation C and the
> distribution shape q each push BACK by 5–10 %. **The reduction is geometric (fewer load-bearing
> contacts per volume), not a flattening of the force tail** — the tail-shape factor moves σ\* the
> other way.

## 1. Frame

- Quantity: σ\* in the solver's convention, `σ*/κ = 1/Q_{1−1/e}(g)`, `g = f_n/(R* p)` (`actshape.py`
  eqs. 2, 4), i.e. the stress at which a fraction 1/e of contacts carries a load above threshold.
  `R* = aᵢaⱼ/(aᵢ+aⱼ)`, κ a force per length (Derjaguin, `PREREG_ACTSHAPE` eq. 1).
- Data: `actshape` dumps (sheared, constant-pressure Yade packs, N 1000/4000, per-contact f_n, r₁, r₂,
  10 steady snapshots per cell), healthy cells only (steady ∧ servo_held ∧ quasistatic): 52.
- Assumptions: (A1) contacts are between touching spheres, branch vector ℓ = (aᵢ+aⱼ)n up to overlap;
  (A2) tangential forces are ⟂ n (Cundall–Strack), so they carry no trace. Nothing else is assumed
  until eq. 12.
- Convention: σ\* ratios are quoted at FIXED number-mean diameter d₁₀ (the `predict_cvmap` "d₁₀"
  convention); since σ\*/κ ∝ L⁻¹, a length rescale by λ multiplies it by λ.

## 2. Derivation

1. Love–Weber (Christoffersen–Mehrabadi–Nemat-Nasser) mean stress, exact for a static assembly:
   `p = (1/3V) Σ_c f_c · ℓ_c`. Under A1–A2, `f_c·ℓ_c = f_n,c (aᵢ+aⱼ)`. Define the closure residual
   `Λ ≡ Σ_c f_n,c(aᵢ+aⱼ)/(3Vp)` (= 1 up to overlap and averaging). (eq. 1–2)
2. Identity: `R*(aᵢ+aⱼ) = aᵢaⱼ` (SymPy: 0). (eq. 3)
3. Substituting `f_n = g R* p`: `Σ_c g_c aᵢaⱼ = 3VΛ`, so `⟨g aa⟩_c = 3VΛ/N_c`. (eq. 4)
4. Define the force–size correlation `C = ⟨g aa⟩/(⟨g⟩⟨aa⟩)` (= 1 if g is independent of pair size)
   and the shape factor `q = Q_{1−1/e}(g)/⟨g⟩`. Then `⟨g⟩ = 3VΛ/(N_c C ⟨aa⟩)`. (eq. 5–6)
5. `σ*/κ = 1/(q⟨g⟩) = (N_c/V)⟨aa⟩_c C/(3qΛ)`. **Exact** given the definitions. (eq. 7)
6. `N_c/V = (Z/2)·φ/⟨v⟩ = 3Zφ/(8π⟨a³⟩)` with Z = 2N_c/N over ALL grains (rattlers included). So
   `σ*/κ = Zφ ⟨aa⟩_c /(8π⟨a³⟩) · C/(qΛ)`. (eq. 8)
7. Pair closures for ⟨aa⟩_c, with aᵢ, aⱼ drawn independently from the number PSD and pair weight w:
   - A: w ∝ (aᵢ+aⱼ)² (geometric cross-section) → `⟨aa⟩ = (m₃m₁+m₂²)/(m₂+m₁²)` (eq. 9); lognormal at
     ⟨a⟩ = 1: `⟨aa⟩/⟨a³⟩ = 1/(1+CV²)` (eq. 10).
   - B: w ∝ aᵢ²aⱼ² (surface × surface) → `(m₃/m₂)²`; lognormal: `1+CV²` (eq. 11).
   - **C: w ∝ aᵢaⱼ(aᵢ+aⱼ) = R*(aᵢ+aⱼ)²** → `⟨aa⟩ = 2m₃m₂/(2m₂m₁) = m₃/m₁`, so
     **`⟨aa⟩_c/⟨a³⟩ = 1/⟨a⟩` for ANY PSD** (eq. 12). A, B, C are the three terms of the
     Boublík–Mansoori–Carnahan–Starling–Leland hard-sphere-mixture contact value
     (g_ij ∝ const + R*·… + R*²·…) `[UNVERIFIED — structure recalled, not looked up; used here only as
     motivation, the closure is tested on data, not on that paper]`.
8. With C: **`σ*/κ = Zφ/(8π⟨a⟩) · C/(qΛ)`** (eq. 13). Monodisperse: `Zφ/(8πa)·C/q`, with C = 1.
9. Z decomposes as `Z = (1−f_r) Z_nr` (f_r rattler NUMBER fraction, Z_nr the backbone count) (eq. 14).
   `geom_isostatic_backbone_2026-09-19.md` shows z_iso contains no length, so Z_nr is forced nearly
   flat in CV; eq. 13 then predicts σ\*(CV)/σ\*(0) ≈ (1−f_r) ratio × φ ratio × C/q ratio.

## 3. Checks (outputs in `out/sigma_star_geometry.txt`)

- **SymPy**: eq. 3 → 0; eq. 9 → 0; monodisperse ⟨aa⟩ = a²; lognormal eq. 10 = 1/s. PASS.
- **Monte Carlo** (2×10⁶ weighted pairs): eq. 9 to ≤ 1.6e-3 at CV 0.2/0.35/0.5. PASS.
- **External anchor (mechanics identity, not this code):** Λ = 0.9994–1.0042 on all 52 cells
  (snapshot range 0.9935–1.0103). The residual rises with p̂ (1.0003 → 1.0029 at μ 0, 1e-4 → 1e-3), the
  sign and scaling expected from contact overlap shortening ℓ below aᵢ+aⱼ. This also validates V,
  the units of p̂·E, and the dump bookkeeping.
- **Reconstruction:** eq. 7 from snapshot means reproduces the measured σ\*/κ to ≤ 0.3 % except two
  low-p̂ μ 0.3 cells (0.8 %, 1.9 %) where the mean of a ratio differs from the ratio of means.
- **Limiting case:** CV 0 → C = 1 identically, ⟨aa⟩ = a², eq. 13 = Zφ/(8πa)/q: 0.1932 predicted vs
  0.1932 measured (μ 0.3, p̂ 3e-4, seed 1).
- **Dimensions:** N_c/V [L⁻³] · ⟨aa⟩ [L²] = L⁻¹ = σ\*/κ ([F L⁻²]/[F L⁻¹]). PASS.

## 4. What the data say (seed-paired, (seed, N) matched, d₁₀-rescaled)

| μ_p | p̂ | CV | R_meas | geometric Zφ·⟨aa⟩/⟨a³⟩ | C ratio | 1/q ratio | ⟨aa⟩_meas/⟨aa⟩_C |
|---|---|---|---|---|---|---|---|
| 0.3 | 1e-4 | 0.5 | 0.794 | 0.701 | 1.063 | 1.070 | 1.077 |
| 0.3 | 3e-4 | 0.5 | 0.733 | 0.640 | 1.070 | 1.070 | 1.060 |
| 0.3 | 1e-3 | 0.5 | 0.788 | 0.677 | 1.093 | 1.066 | 1.048 |
| 0 | 3e-4 | 0.5 | 0.785 | 0.719 | 1.053 | 1.036 | 1.009 |
| 0.3 | 3e-4 | 0.2 | 0.955 | 0.935 | 1.016 | 1.006 | 0.991 |

(all 11 rows in the output file). The product reproduces R_meas to ≤ 1.1 % (exact identity).
Closures A and B fail (⟨aa⟩/⟨a³⟩ ratio 0.80 and 1.25 at CV 0.5 against a measured 1.01–1.08); C is
within 2 % for μ 0 and within 8 % for μ 0.3.

**The Z split (qualitative only).** From the records' own `rattler_frac_ss`, `z_nr_ss`: at CV 0.5 the
(1−f_r) ratio is 0.50–0.68 while Z_nr moves 0.85–1.01 — the rattler fraction carries the drop
(sheared μ 0.3 CV 0.5 packs have f_r = 0.40–0.68 by number). ⚠ The records' (1−f_r)Z_nr·φ does NOT
reproduce the dump-exact Zφ (0.461 vs 0.650 in the worst row): the record fields are sampled or
thresholded differently from the dumps (not traced). Only the dump Zφ enters the exact statement.

## 5. Consequences

1. **σ\* is naturally a d₁₀ quantity**, and eq. 13 says why: with closure C the PSD enters only through
   ⟨a⟩ and Zφ. The "d₃₂ convention" bracket (`predict_cvmap` H2, ×1.145) has no geometric basis here.
2. **The earlier reading "polydispersity flattens the force tail, hence lowers σ\*"
   (`three_knobs` §6, card 2026-09-22) is the wrong mechanism.** The tail factor q moves σ\* UP by
   4–7 % at CV 0.5; the whole reduction is contact density.
3. Transfer to the bench: σ\*(CV)/σ\*(0) ≈ [Zφ]_CV/[Zφ]_0 at fixed number-mean size, and Z is set by
   the rattler population under shear. A measurement of Z (or rattler fraction) in a sheared
   polydisperse pack is a σ\* measurement.

## 6. Registered for `lfsig` (P4, before any lfsig cell is read)

In LF_DEM the activated fraction at stress σ is the fraction of interacting pairs with gap h < 0, and
the threshold is exactly Derjaguin (`RepulsiveForce.cpp`: contact force = F_R·R*). Prediction: the
ratio of rheological onsets R_LF = σ\*_η(CV 0.5)/σ\*_η(CV 0.2) (PREREG_LFSIG P1) equals
[Z_c φ]_{0.5}/[Z_c φ]_{0.2} × 1.05 within ×1.2, with Z_c the LF_DEM contact number per particle
(`data_*` col 19) averaged over the σ = 30 and 100 cells at the same x_m, and 1.05 the Yade C/q ratio
for CV 0.5 vs 0.2 at μ 0.3 (≈ (1.070·1.070)/(1.016·1.006)·… rounded; recomputed exactly in the grader
from `out/sigma_star_geometry.txt`).

## 7. Open

- A derivation of the rattler fraction f_r(CV) under shear (the only non-geometric input left).
- Whether closure C's 5–8 % residual at μ 0.3 is the rattlers being small (the backbone is a
  size-biased subset); a per-particle dump would decide it.
- The BMCSL attribution in step 7 is unverified.
