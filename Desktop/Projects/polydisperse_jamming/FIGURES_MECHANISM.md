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

## M2 — the sign, and undrained locking

`figures/M2_transient.{pdf,png,csv}`

⚠ **REBUILT 2026-09-08.** The first version of this figure showed the pore pressure
peaking near `γ ≈ 0.55` and decaying, and a caption calling that a *strain clock*. That
decay was an **artefact of a missing term** in the `φ`-evolution law (`poro_impact.py`
eq. iv), found by a referee. The claim is withdrawn and the corrected behaviour is the
opposite. Both panels are linear and include zero.

**(left) Pore pressure at the impermeable substrate, in the undrained limit
(`Λ` = 10⁻⁶, which is where the lab's silica sits: `Λ` ≈ 10⁻⁷).** The three curves are
`φ > φ_c`, `φ = φ_c` and `φ < φ_c`. Three things to read off it:

- the **sign flips** across `φ_c` — suction above, overpressure below;
- the dilatant and contractive branches are **exact mirror images** at this level of the
  model, so the mechanism is a sign flip and not an asymmetry;
- the `φ = φ_c` curve is **identically zero**, which is the null the sign claim needs;
- and the pressure **grows linearly and without bound**. It does not peak and it does
  not decay.

**(right) Why — and why only when undrained.** `φ − φ_c` against strain, at the same
two `Λ`.

- **Undrained (solid):** the pack **locks**. As it dilates plastically the pore pressure
  falls, so the effective stress rises, so the skeleton compresses elastically by
  exactly what it just dilated. Net volume change zero — which is what "constant volume"
  *means*. With `α` = 1 the cancellation is exact and `φ − φ_c` is pinned at its initial
  value to machine precision. This is **dilatancy hardening**, the standard reason a
  dense saturated sand is strong under fast loading.
- **Drained (dashed, `Λ` = 2):** the fluid carries pressure away, the cancellation is
  incomplete, and the pack really does relax toward `φ_c`. **The strain clock exists —
  but only on the drained side, which is not where this system is.**

What actually bounds the undrained growth is *not* in this model: the critical state is
a **line**, not a point, so `φ_c` falls as the effective stress rises and undrained shear
ends when `φ_c(σ')` comes down to meet the fixed `φ`. With `φ_c` constant there is no
shut-off at all, and the network-modulus ceiling of M1(a) is the stand-in. Closing that
needs a compression index nothing here measures.

---

## M3 — `φ_c(CV)`

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

---

## M4 — the whole phase-8 result in one panel

`figures/M4_pss_vs_phi.{pdf,png,csv}`

Constant-volume steady-state grain pressure against packing fraction, one ladder per CV,
`μ_p` = 0.3, lognormal, 3 seeds, protocol `cv-shear-v2`.

**The three ladders are the same curve displaced to the right.** That displacement is
Δ`φ_c` = +0.062 over CV 0 → 0.5, and it is the measurement the mechanism turns on.

The headline is the **vertical gap where two ladders overlap in `φ`**: near φ = 0.635 the
CV = 0 ladder sits at `p̂` ≈ 1.1e-2 and the CV = 0.5 ladder at 2.5e-4 — **43×**, at the
same packing fraction, with the CV = 0.5 point marginally *denser*. Read the range 5.4×
to 43× rather than one number: the two overlapping rungs tie on `φ` mismatch and the
ratio grows as `φ_c` is approached.

`y` is **logarithmic** — `p̂_ss` spans four decades. `x` is linear and does **not** start
at zero, because `φ` is a fraction near 0.6 while the ladders span 0.09; the range shown
is the range occupied, and it is stated here rather than implied by the axis.
