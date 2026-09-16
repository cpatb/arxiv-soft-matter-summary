# PREREG — the shape axis, 2026-09-18 (design registered; no cell run)

The author's call on the design fork (2026-09-18): **clumps now, LevelSet as the
control** — build the cheap axis, then measure the interlocking confound instead of
caveating it. This file registers the parameterisation, what it holds fixed, the two
defects a capability probe found before anything was built on them, and the gates.

**Nothing has been run.** `shape.py` and `test_shape.py` (19 anchors) are the geometry;
`scripts/probe_clump.sh` is a capability probe that writes no shard.

---

## 1. The parameterisation, written out — because a family is a PATH

The 2026-09-13 lesson: a one-parameter family is a path through a multi-parameter plane,
and that is exactly how this project manufactured its "bidisperse is non-monotone in CV"
negative — the family was 50:50 **by number**, so its fines carried 0.500 → 0.036 of the
**volume** as CV rose and it never visited the composition where the effect peaks. So the
shape axis is written out before it is run.

A grain is a rigid clump of two identical spheres of diameter `ds`, centres `s` apart,
with `u = s/ds` and **aspect ratio α = 1 + u**, α ∈ [1, 2]. The union volume is
`V = ds³·g(u)` with

    g(u) = π/3 − (π/12)(2+u)(1−u)²   (u ≤ 1),      g(u) = π/3   (u ≥ 1)

so the sub-sphere diameter that makes the grain's **volume-equivalent diameter** exactly
`d_eq` is `ds = d_eq·((π/6)/g(u))^(1/3)`.

**HELD FIXED along the α axis, by construction and not by hope:** `d_eq` is drawn from the
*same* `psd.make_psd`, so CV, `d32`, `d43`, skewness, span and every other PSD descriptor
are **identical at every α**. `mu_p`, `N`, the pressure ladder and the protocol are
unchanged. **MOVED:** α, and nothing else about the size distribution.

That makes **(CV, α) a grid, not a path**. The registered grid is
`α ∈ {1.0, 1.25, 1.5, 1.75, 2.0}` × the existing CV ladder, and **α = 1 is the first rung
deliberately**: at α = 1 the construction returns `ds = d_eq` exactly — the grain *is* the
sphere the whole project has run — so the shape campaign must reproduce the existing
sphere results at its own α = 1 before anything at α > 1 is believed.

**MOVED AS A SIDE EFFECT, and this is the confound:** a clump of overlapping spheres has a
bumpy surface. The waist depth relative to the equivalent radius is
`eps(u) = (ds/d_eq)·[1 − √(1−u²)]`, which runs **0 at α = 1 to 2^(−1/3) = 0.7937 at α = 2**,
where the neck radius is exactly zero. This is not a caveat, it is a second knob: at the
top of the ladder the grain's surface is as corrugated as a two-sphere body can be.
`eps` is recorded on every cell.

**THE CONTROL THAT TURNS IT INTO A MEASUREMENT.** A LevelSet body — smooth, matched α,
matched volume, `eps = 0` by construction — run at the same (CV, α) points. The
clump-minus-LevelSet difference *is* the interlocking term. Until that control exists,
**no result from the clump axis may be attributed to elongation rather than to
interlocking**, and the artifacts must say so in those words.

---

## 2. What the capability probe found, before anything was built on it

`scripts/probe_clump.sh`, α = 1.5, 60 clumps in a periodic cell, 50 steps. Four questions,
two of them answered badly, which is why the probe existed.

1. **Clumps work in the protocol.** 60/60 built in a periodic cell; the collider accepts
   them with `allowBiggerThanPeriod = False`; `PeriTriaxController`-style integration ran.
2. 🔴 **`utils.porosity()` DOUBLE-COUNTS the member overlap.** Measured: φ from the
   analytic union volume 0.06135923, φ from `utils.porosity()` **0.07272205**, ratio
   **1.185185**, which is exactly `2(ds/d_eq)³`. **This ratio is α-dependent and runs
   2.000 at α = 1 down to 1.000 at α = 2** — so an uncorrected clump campaign would have
   reported φ **falling by a factor of two along the aspect ladder from the bookkeeping
   alone**, in a project whose central observable is φ. `packing.measure` must compute φ
   from the analytic union volume for clump runs; `utils.porosity()` must not be called.
3. **No internal contacts.** 61 real interactions, 0 of them between members of the same
   clump, so `n_contacts`, `z` and the fabric need no filter.
4. 🔴 **Mass carries the same double count.** Clump mass 0.62056151 against a union mass
   of 0.52359878 — the same 1.185185. The material density must be scaled by
   `(π/6)/(2(ds/d_eq)³·(π/6))` per α, which is a constant independent of grain size, so
   one density per α fixes it.

Both defects are of the same kind, and neither would have raised an error.

---

## 3. Gates, registered

**S0 — α = 1 REPRODUCTION, runs first, exits non-zero.** At α = 1 the clump path must
return `phi`, `z_all`, `z_nr` and `rattler_frac` **bit-identical** to the existing sphere
records at the same (family, CV, seed, `p_hat`, N). `phi` is the gate (continuous
channel). This is the strongest available check on the whole construction: it tests the
volume bookkeeping, the density scaling, the insertion RNG order and the collider path in
one number, and an off-by-anything in §2 breaks it.

**S1 — the PSD must not move.** The realised `cv_realised`, `d32_realised`, `d43_realised`
recorded at α > 1 must equal those at α = 1 to 1e-12, at every (family, CV, seed). If the
construction is volume-preserving this is exact; it is checked rather than assumed.

**S2 — the cell must still be legal.** The longest extent is `span_factor(α)·d_eq`, which
is `α·((π/6)/g(α−1))^(1/3)` and reaches **2^(2/3) = 1.5874** at α = 2. N = 1500 clears the
periodic collider at α = 1 and **does not automatically clear it at α = 2**; the
clearance check is re-run per α and the required N is reported, never worked around by
tightening the truncation (which would move `d32`, the one moment this project turns on).

**S3 — the primary claim.** Per family, the paired CV 0 → 0.5 shift in φ_J at each α. The
registered question: **does the shape-of-the-distribution effect (the existing 4.47×
across families) survive, grow or vanish as grains stop being spheres?** No direction is
predicted. The grid is 2-D by construction, so the answer is read as a surface and not as
a path.

**S4 — the confound, quantified not caveated.** Every α > 1 result is reported with its
`eps`, and every S3 number is labelled **"elongation + interlocking"** until the LevelSet
control exists. A single-number claim about "aspect ratio" from clumps alone is
pre-registered as **not quotable**.

---

## 4. What this design cannot decide

Nothing about **absolute size** — the static DEM is proved scale-free, and a size axis
means a different model. Nothing about **concave or faceted** grains: a two-sphere clump
is a fixed shape family with one parameter, and α > 2 is a different family (the spheres
separate) needing more sub-spheres. And nothing about the **direction** of any fabric
anisotropy: orientations are drawn uniformly on the sphere (`shape.random_axes`, checked
against a θ–φ sampler that scores `n·aniso²` above 50 where uniform scores under 6), so
the pack has no built-in preferred axis to find.
