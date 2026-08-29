# MOC physics flags surfaced while building drying_analysis (2026-07-06)

Per the lab CLAUDE.md ("Do not fix physics silently; flag and let the author decide"),
these are discrepancies noticed while implementing the roadmap observables. Code was
written to the *correct* physics; the MOC prose needs the author's eye.

## Flag 1 — `ℓ ≈ l_fc for brittle low-stretch` is inverted (hydrogel_MOC.md §5)

**MOC text (§5, "Key relation"):**
> `ℓ = Γ/E = l_fc · f(λ_c)`, with neo-Hookean `f = (λ_c²+2/λ_c−3)/6`. For brittle
> low-stretch PEGDA (λ_c ≈ 1.1–1.5), `ℓ ≈ l_fc`; for stretchy gels `ℓ` is 4–16× larger.

**Problem:** `f = W_c/E` is the ratio of the network length `ℓ=Γ/E` to the
fractocohesive length `l_fc=Γ/W_c` (since `ℓ/l_fc = (Γ/E)/(Γ/W_c) = W_c/E = f`).
`f → 0` as `λ_c → 1`, so for a **brittle low-stretch** gel `ℓ` is **much smaller**
than `l_fc` (e.g. f(λ_c=1.15) ≈ 0.013 pure-shear, ≈ 0.010 uniaxial), **not** `≈ l_fc`.
`ℓ` grows *toward and past* `l_fc` only for stretchy gels. The "brittle → ℓ≈l_fc"
and the direction of the inequality read backwards.

**Also (reviewer catch):** the MOC writes the uniaxial neo-Hookean `f=(λ_c²+2/λ_c−3)/6`,
but the toughness pipeline and the notch series are **pure shear**, whose
`f=(λ_c²+1/λ_c²−2)/6`. `ell_from_lfc()` uses the pure-shear form for consistency.

**Resolution in code:** `toughness.ell_from_lfc` computes the correct pure-shear
`f` (validated exactly). The direction (`ℓ << l_fc` for brittle, growing with stretch)
is asserted in `validate.py`. **Author action:** correct the §5 prose: likely intended
"`ℓ` is a small fraction of `l_fc` for brittle gels, approaching `l_fc` for stretchy ones."

## Flag 2 — Biot-number length convention (hydrogel_MOC.md §5A vs §6.3)

§5A writes `Bi = k_m ℓ / D` with `ℓ` the **full** disc thickness; §6.3 writes
`Bi(t) = k_m H(t)/D_p` with `H` the side-view **half**-thickness. These differ by 2×.
`swelling_kinetics.fit_D_km_robin` reports `Bi = k_m·(L/2)/D` (half-thickness = the
physical diffusion length of a symmetric two-face sheet) and documents it. **Author
action:** pick one convention across the docs so `Bi_c` comparisons are consistent.
