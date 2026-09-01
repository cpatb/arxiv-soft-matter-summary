# Polydispersity, shear thickening and jamming: what is already done

Written 2026-08-31, during the build. **This corrects the novelty claim in
`SPEC_2026-08-31.md` §1.3.**

That claim was checked only against `OPEN_PROBLEMS_MASTER.md`, which is a **drop-impact**
slate. It is silent on the suspension-rheology literature, and the rheology half of this
project's question turns out to be actively worked. The honest position is below.

---

## 1. What is already done, and by whom

**Singh, Tyagi & Santra, arXiv:2510.16464 (Oct 2025), "Rigidity transition in polydisperse
shear-thickening suspensions".** `[web-verified 2026-08-31, abstract read]`
2D LF-DEM at high φ and large dimensionless stress. Rigidity of stress-induced clusters
analysed with the **pebble game**; order parameter `f_rig` (fraction of particles in rigid
clusters) scales as `(φ − φ_c)^β` with 2D-percolation exponents, susceptibility
`χ_rig ~ |φ − φ_c|^(−γ)`, finite-size `ν ≈ 1.33`.

**This is the closest published work to this project's Phase 2/3, and it carries two
findings that bear directly on my anchors:**

1. **`φ_c` and `φ_J^μ` vary NON-MONOTONICALLY with the polydispersity index.** My anchor
   **A5** assumed φ_RCP rises monotonically with polydispersity. That assumption is
   contradicted for the *frictional, shear-jamming* φ_J in 2D. A5 must be demoted to a
   direction-only anchor for the **frictionless** case (where Desmond & Weeks applies) and
   must not be asserted for `φ_J^μ`.
2. **`φ_c` and `φ_J^μ` depend on particle stiffness.** This project's control parameter is
   the reduced pressure `p_hat = p/E`, which is exactly that stiffness axis, and φ_J is
   taken as the `p_hat → 0` extrapolation. Their result predicts that limit may not be
   cleanly attained, **and the non-monotonic frictionless φ(p̂) ladders observed here
   (0.6434, 0.6319, 0.6352, 0.6313 over 1.5 decades) are consistent with exactly that.**

3. Polydisperse systems behave identically to **statistically equivalent bidisperse** ones,
   for both `f_rig` and `χ_rig`. If that holds in 3D it is a strong simplification, and it
   is a testable prediction against this project's uniform-vs-lognormal comparison.

**Guy, Ness, Hermes, Sawiak, Sun & Poon, arXiv:1901.02066 (2019), "Testing the
Wyart–Cates model for non-Brownian shear thickening using bidisperse suspensions".**
`[web-verified 2026-08-31]` Uses bidispersity as a deliberate probe of the WC model's
contact-fraction `f`. Same group whose Edinburgh DataShare deposits this lab already uses
(see the `edinburgh-datashare-rheology` memory). Must be cited by anything here that
touches WC and size distribution together.

**Malbranche (CUNY thesis) and the bidisperse-DST line more generally.**
`[web-verified 2026-08-31, secondary source]` The reported consensus is that **increasing
the particle size ratio weakens thickening, DST gives way to CST**, and a 2024 study
reports a CST→DST transition accompanied by an abrupt rise in the **onset stress**.

⚠ **That direction is the author's bench observation.** More polydisperse → weaker
thickening → more spreading and inertia. So the observation may already be explained by
published bulk rheology, and the novel content of this project is *not* "polydispersity
weakens thickening", that appears to be known.

**Desmond & Weeks, *PRE* **90**, 022204 (2014), arXiv:1303.4627.** `[web-verified]`
φ_RCP as a function of polydispersity **and skewness**, with an empirical formula.
Skewness can outweigh polydispersity; kurtosis barely matters. This is what motivated
adding the zero-skewness `UniformPSD` family here.

**Silbert, *Soft Matter* **6**, 2918 (2010), arXiv:1108.0012.** `[web-verified]`
φ_J falls monotonically 0.64 → 0.55 and `Z` 6 → 4 as friction rises over decades. This
project's measurements (φ_J 0.6332 ± 0.0002 at μ = 0, 0.5744 ± 0.0004 at μ = 0.3;
`Z` 6.37–6.47 → 4.94) sit on that published curve, a genuine external anchor, passed.

---

## 2. What this project can still claim

**Not novel:** that polydispersity weakens thickening; that φ_J depends on the size
distribution; rigidity/percolation analysis of polydisperse thickening suspensions.

**Still open, and where the value is:**

1. **The free-surface / impact geometry.** Every item above is bulk rheology or a
   rheometer geometry. `OPEN_PROBLEMS_MASTER.md` §3 still lists "intrinsic bulk transition
   against boundary-mediated dilatancy, unresolved **for a free surface**, where no rigid
   wall supplies confinement" as contested. Phase 4, free vs capillary vs caged on the
   *same* carved aggregate, is the piece nobody has done.
2. **3D, and the monotonicity question.** Singh et al. is **2D**. Whether `φ_J^μ` is
   non-monotonic in polydispersity in **3D** is not settled here, and this project's
   Phase-2 sweep measures it directly.
3. **Separating polydispersity from skewness.** Desmond & Weeks showed skewness can
   dominate for frictionless RCP; nobody appears to have separated the two axes for
   *frictional* `φ_J`. The lognormal-vs-uniform comparison at matched CV does exactly
   that, and it is cheap.
4. **The moment-matching point itself.** The `d₃₂/d₁₀ = (1+CV²)²` trap, that "same
   particle size" names one moment while the permeability reads another, is a
   methodological result that applies to every experimental comparison of a monodisperse
   and a polydisperse suspension, and I have not found it stated in this literature.

---

## 3. Actions this forces

- **`SPEC_2026-08-31.md` §1.3's novelty claim is narrowed** to items 1–4 above. Done.
- **Anchor A5 is demoted**: monotonic φ_RCP rise is asserted only for the frictionless
  case; for `φ_J^μ` the sweep must *measure* monotonicity rather than assume it, and the
  analysis reports the ordering rather than grading against an assumed direction.
- **Read arXiv:2510.16464 in full before writing anything up.** Only the abstract has been
  read. Its `f_rig` order parameter is close to this project's largest-backbone-fraction
  observable and the two should be compared explicitly.
- Both Poon-group references belong in any write-up that touches Wyart–Cates.
