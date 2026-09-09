# Figures for the dilatancy-sign mechanism

Built by `make_mechanism_figures.py`. Each figure has a sibling `.csv` in `figures/`
carrying every plotted point, so no number here is read off a picture. Regeneration is
byte-identical (`python make_mechanism_figures.py --check`).

House style: minimal text on the axes; the explanation is here, not on the plot.

---

## M1 — the sign, and the length nobody has measured

`figures/M1_regime_map.{pdf,png,csv}`

**(a) The pore pressure changes sign at the critical state.** `Δp_f` against `φ − φ_c`
for 800 nm silica at `φ` = 0.44, `U₀` = 3 m/s, `D₀` = 2.5 mm. Positive is
**overpressure** (the pack contracts and drives fluid out, effective stress falls,
the drop liquefies); negative is **suction** (the pack dilates and draws fluid in,
effective stress rises, the drop arrests). The vertical dotted line is `φ = φ_c` and
crossing it is the whole mechanism.

Three curves, because the honest answer is the inner envelope of two:

- **solid** — the network-modulus ceiling, `p_dil^max = K_p |Δφ|/φ`
  (`PORE_PHYSICS.md` Eq. 4.5). Linear in `Δφ`, which is why capping does not disturb
  the sign structure.
- **dashed** — the uncapped Darcy expression at the drop-scale drainage length.
- **dash-dot** — the same at a grain-scale drainage length, `10 d₃₂` = 8 µm.

🔑 **The ceiling lies below BOTH Darcy curves everywhere on the plot.** So although the
drainage length moves `Π` by 9.8 × 10⁴ (panel b), **it does not move the pore-pressure
magnitude at all** — the packing runs out of dilation strain first either way. The
section-6 ambiguity is an ambiguity about `ID`, the force *ratio*, not about the
pressure.

⚠ The Pa scale carries `K_p` from `impact_stress`'s double-layer route, whose `P_conf`
is the **withdrawn** `si_fit` `σ*` = 44.2 Pa. The *sign*, the *shape* and the *ordering*
of the three curves do not depend on it; the absolute Pa values do. The ratio form
`p_dil^max/P_conf = (a/6λ_D)(Δφ/φ)` avoids it entirely.

The `y` axis is **symlog** (linear within ±1 Pa, logarithmic outside): the quantity
spans six decades *and* changes sign, and neither a linear nor a log axis can show both.

**(b) `Π` against the drainage length, at `Δφ` = 0.02.** `Π ∝ ℓ²`, dashed line at
`Π = 1`. The two markers are the grain-scale and drop-scale readings. They straddle the
line, and the crossing is at

> **`ℓ*` = 21.2 µm = 26.5 grain diameters**, i.e. `ℓ*/D₀` = 8.5 × 10⁻³.

That is squarely inside the range a granular shear band occupies (5–15 grains is the
usual quoted width, so `ℓ*` is a factor ~2 above it). **This system sits on the
drained/undrained boundary, and which side it is on is decided by a length that has
never been measured for an impacting drop.** Both axes are logarithmic and each spans
more than a decade.

---

## M2 — the transient, and the strain clock

`figures/M2_transient.{pdf,png,csv}`

From `poro_impact.py` at `Λ` = 2, `S` = 1, 201 nodes. Both panels are linear and both
include zero, because zero is the quantity of interest in each.

**(left) Pore pressure at the impermeable substrate against strain.** The three curves
are `φ > φ_c`, `φ = φ_c` and `φ < φ_c`. The dilatant and contractive branches are exact
mirror images at this level of the model — the mechanism is a sign flip, not an
asymmetry — and the `φ = φ_c` case is identically zero, which is the null the sign
claim needs.

The shape matters as much as the sign: the pressure **peaks near `γ ≈ 0.55` and then
decays**. It is not a steady state.

**(right) Why.** `φ − φ_c` relaxes toward zero as the pack dilates (or contracts), so
the source shuts itself off. Two consequences:

1. **The closed form is an upper bound, not a prediction.** The transient always
   undershoots `Π/2`; `test_poro.py` pins that it must.
2. **The effect has a lifetime of a strain of order `1/α`.** This is the model's own
   account of the author's bench phrase *"inertial-dominated for some time"* — a strain
   clock rather than a permanent state, and it is a shape a rival mechanism built on the
   suspension viscosity does not predict, because a viscosity is a state function.

---

## M3 — `φ_c(CV)`, and where it crosses `φ`

`figures/M3_phi_c_vs_CV.{pdf,png,csv}` — **built only when `out/phase8_summary.csv`
exists.** The script prints a reason and skips otherwise rather than drawing an empty
panel.

Four estimators (`analyse_phase8.py`: the `φ(I → 0)` extrapolation, two pressure roots
and the coordination crossing), mean ± seed sd.

**No `φ` = 0.44 reference line is drawn**, deliberately.

⚠ `y` does **not** start at zero here, against the house rule, because `φ_c` is a
packing fraction near 0.6 while the effect under test is ~0.06. The suppressed range is
stated instead: the axis shows only the interval actually occupied, and the differences
between CV values — not their ratio to zero — are what the figure is for.

⚠ **The absolute level of `φ_c` from a frictional-sphere DEM (~0.58) is not
commensurable with the `φ_c` = 0.37–0.42 that `impact_stress` uses for cornstarch and
silica.** Only the *shift* `dφ_c/dCV` transfers, which is why the figure carries no
`φ` = 0.44 line: an intersection with it would be a `CV*` computed from two numbers that
are not on the same axis.
