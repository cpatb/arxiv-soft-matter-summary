# The aggregation figures — what A1 and A2 show

Built by `make_multiplet_figures.py` (`--check` rebuilds and compares md5; both are
byte-identical). Plotted data ships beside each figure as `.csv`. Derivation:
`derivations/multiplet_psd_2026-09-22.md`. Numbers: `out/multiplet_levers.txt`.

Axes carry almost no text by house rule; this file is the caption.

> 🔴 **THIS FILE WAS STALE AND SHIPPED TWO RETRACTED CLAIMS.** It was written before
> referee pass 2 and not revisited, so the figure it describes used an admissibility
> ceiling pass 2 had already rejected, and this caption restated a precision claim the
> derivation's own result box retracts. Both are corrected below, and both are recorded
> rather than quietly replaced — a downstream artifact going stale when an upstream
> correction lands is a defect this repository already knows about, and it happened here
> inside one session.

---

## A1 — the aggregation plane

`log₂(Λ/G)` over `(d_f, ⟨k⟩)`, where `d_f` is the multiplets' mass-fractal dimension
(3 = compact, lower = more open) and `⟨k⟩` the mean aggregation number (1 = unaggregated).

- **`Λ` is the trapped-fluid factor**, `⟨k^γ⟩/⟨k⟩` with `γ = 3/d_f`: the factor by which
  aggregation inflates the volume fraction the *packing* sees, at fixed solid loading.
  Closed form, no fitted parameter.
- **`G` is the polydispersity gain**, `φ_J(multiplet PSD)/φ_J(mono)`, from this project's
  calibrated filling functional.
- **Red = aggregation promotes jamming** (`Λ > G`); **blue = it suppresses** jamming.

### The two contours are the whole figure, and they are contours of the SAME field

| line | condition | meaning |
|---|---|---|
| **solid** | `Λ/G = 1` | aggregation starts **promoting** jamming |
| **dashed + hatched** | `Λ/G = φ₀ᵐᵒⁿᵒ/φ_solid = 1.2162` | **assumption A3 fails** beyond it |

The second one is new. Rigid, non-interpenetrating aggregates cannot be packed denser
than where that effective-sphere pack **jams**, so A3 fails when
`Λ·φ_solid > G·φ₀ᵐᵒⁿᵒ` — which rearranges to a contour of `Λ/G`, the very quantity
plotted. At `φ_solid = 0.52` (the solver's default) and `φ₀ᵐᵒⁿᵒ = 0.6324` (the paper's
own `\phiJFrictionless`, independently reproduced to 1.2 × 10⁻¹⁴ this session), that is
`Λ/G = 1.2162`.

> 🔑 **THE ADMISSIBLE WINDOW IS NARROW, and that is the honest reading of this figure.**
> Between the solid line and the dashed one there is only a thin band, around
> `d_f ≈ 2.6 – 2.9`, in which aggregation *promotes* jamming **and** the rigid-aggregate
> model still applies. Everything further left — including the whole `d_f ≈ 1.8 – 2.1`
> region where reported colloidal aggregates sit, and where the effect looks largest — is
> **hatched**: there the aggregates must deform or interpenetrate, which is physics this
> model does not contain.
>
> ⚠ The previous version of this figure hatched `φ_eff > 1` instead, a trivial
> space-filling bound (`Λ ≤ 1.923` rather than `Λ/G ≤ 1.216`, **1.58× looser**), which
> drew a large wedge of inadmissible parameter space as valid coloured data.

### What the crossing does and does not establish

The `Λ = G` contour is **nearly vertical**: `d_f*` moves only from 2.893 to 2.946 across
a sevenfold range in `⟨k⟩`.

> 🔴 **THAT FLATNESS IS NOT ESTABLISHED AT THAT PRECISION**, and an earlier version of
> this caption stated it as the result. The plateau is 0.053 wide, and one `PHI_SCALE` of
> protocol noise moves the *root* by 152 % / 80 % / 36 % of it at `⟨k⟩` = 1.1 / 1.2 / 1.5
> (`out/multiplet_levers.txt` §2). The low-`⟨k⟩` roots are less well pinned than the
> spread they are being compared across.
>
> **What survives is ORDINAL, and it is enough:** every `d_f*` sits near **2.9** — far
> below the compact limit 3 and far above the reported colloidal 1.8–2.1 — on all three
> closure calibrations. The conclusion is about which *side* of the crossing a real
> aggregate falls on, and the worst sensitivity (0.081) is an order of magnitude smaller
> than the 0.7–1.1 gap to those values.

⚠ The blue (suppressing) region is a narrow sliver at the right edge: over the whole
plane `log₂(Λ/G)` runs **−0.062 to +4.850**, so the suppressing side never reaches a
factor of 1.05. The colour is therefore **clipped at ±2** with the colorbar marked
`extend`; the `.csv` carries unclipped values.

## A2 — the two levers, and the closure-free form of the argument

`Λ` (solid) and `G` (dotted) against `d_f`, at three aggregation numbers. Log `y`, because
`Λ` spans more than a decade; the house rule is linear-from-zero *unless* the data spans a
decade, and it does. Markers are the `Λ = G` crossings.

**The dashed horizontal line is measured, not fitted**: `G_max = 0.6489/0.5947 = 1.0911`,
the widest pack this project has ever made (lognormal `CV = 0.5`) over the measured
monodisperse `β`. Every dotted curve stays under it.

That is the argument without the closure:

> `G` is **penned into a band measurement already bounds — about +9 % at most.** `Λ` is
> not. So the question is only *where `Λ` passes 1.0911*, which is a closed form.

It passes above `⟨k⟩ =` **1.08 / 1.14 / 1.37** at `d_f =` 1.8 / 2.1 / 2.5. ⚠ **Below
those thresholds the closure-free test is uninformative** — the fitted closure still gives
the same *sign*, but by a per cent rather than a factor.

---

## What these figures do NOT show

1. **Λ is an UPPER BOUND, not a value.** Both figures assume aggregates are impenetrable
   at their effective diameter. Real fractal aggregates interpenetrate, which lowers `Λ`.
   A1's hatching is now the *quantitative* statement of where that assumption has already
   broken, rather than a remark.
2. **No shape.** Each multiplet is replaced by one sphere. A real doublet is a dumbbell.
   The clump control is registered in `PREREG_MULTIPLET_2026-09-22.md` §6 and is **blocked
   on volume bookkeeping**, not on cells.
3. **No `σ*`.** Neither figure says anything about the onset stress; which area it reads
   is a convention owned by two other files (`multiplet_wc.py` §4).
4. **`G` is the fitted half**, and its form is on record as too rigid. `d_f*` moves ≤ 0.01
   between the two refit calibrations and ~0.1 to the textbook one — which is why A2 draws
   the *measured* ceiling rather than leaning on the fit.
5. **A1's A3 boundary uses the FRICTIONAL `G`** as a stand-in for the frictionless one,
   because that is the gauge the closure is calibrated in. `multiplet_wc.py` §2 finds no
   resolved friction dependence of `G` (0 of 4 CV levels) — which is a **bound, not a
   proof of equality**, and the `mustarcap` campaign is what widens it.
