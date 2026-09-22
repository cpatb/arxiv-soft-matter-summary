# The aggregation figures — what A1 and A2 show

Built by `make_multiplet_figures.py` (`--check` rebuilds and compares md5; both are
byte-identical). Plotted data ships beside each figure as `.csv`. Derivation:
`derivations/multiplet_psd_2026-09-22.md`. Numbers: `out/multiplet_levers.txt`.

Axes carry almost no text by house rule; this file is the caption.

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

**The solid contour is the `Λ = G` crossing, and it is almost VERTICAL.** That is the
result: `d_f*` sits at **2.89–2.95 across a seven-fold range in `⟨k⟩`**, so

> **the SIGN of the aggregation effect is set by how OPEN the aggregates are, not by how
> much the suspension has aggregated.**

⚠ **The blue region is a narrow sliver against the right-hand edge.** Suppression happens
only for nearly-perfectly-compact aggregates. Over the whole plotted plane
`log₂(Λ/G)` runs **−0.062 to +4.850** — the suppressing side never reaches even a factor
of 1.05 — which is why the colour is **clipped at ±2** with the colorbar marked `extend`.
The `.csv` carries the unclipped values.

**The dashed line and hatching are the model refuting itself.** `φ_eff = Λ φ_solid` is a
volume fraction and cannot exceed 1, so at the Wyart–Cates solver's own `φ = 0.52` the
bound is `Λ ≤ 1.923`. Inside the hatching that bound is **breached**: open aggregates at
moderate loading *require* interpenetration, so assumption A3 (impenetrable at the
effective diameter) is false there and **no cell in the hatched region may be run or
quoted**. It is the sharpest thing in the figure — a model bounding its own domain from
the inside — and it eats exactly the corner where the effect looks largest.

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

It passes above `⟨k⟩ =` **1.08 / 1.14 / 1.37** at `d_f =` 1.8 / 2.1 / 2.5. **Above those
thresholds no width effect this project has ever measured can offset the trapped-fluid
lever.** ⚠ **Below them the closure-free test is uninformative** — the fitted closure still
gives the same *sign* there, but by a margin of a per cent, not a factor of several. An
earlier draft quoted `Λ = 1.22–4.11` as if it bounded "any realistic aggregate"; that is
the range over five tabulated `⟨k⟩`, and a lightly aggregated suspension is equally real.

---

## What these figures do NOT show

1. **Λ is an UPPER BOUND, not a value.** Both figures assume aggregates are impenetrable
   at their effective diameter. Real fractal aggregates interpenetrate, which lowers `Λ`.
   That matters most where `Λ` is largest — the left of A1, which is hatched out anyway —
   and least at the crossing, where aggregates are nearly compact. So the *location* of
   the boundary sits in the most defensible part of the model; the *magnitude* deep in the
   red does not.
2. **No shape.** Each multiplet is replaced by one sphere. A real doublet is a dumbbell.
   The clump control that would test this is registered in `PREREG_MULTIPLET_2026-09-22.md`
   §6 and is **blocked on volume bookkeeping**, not on cells.
3. **No `σ*`.** Neither figure says anything about the onset stress; which area it reads is
   a convention owned by two other files (`multiplet_wc.py` §4).
4. **`G` is the fitted half**, and its form is on record as too rigid. `d_f*` moves by
   ≤ 0.01 between the two refit calibrations and ~0.1 to the textbook one — which is why
   A2 draws the *measured* ceiling rather than leaning on the fit.
