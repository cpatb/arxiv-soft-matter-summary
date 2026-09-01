# F3 — structure at matched mechanical state

**What is plotted.** For each structural observable, the change between the lowest and
highest CV at *matched* reduced pressure `p̂` and friction `μ_p` — i.e. after the packings
have been put in the same mechanical state. Bars are grouped by `p̂`; error bars are the
seed standard deviation.

**Why matching on `p̂` is enough.** A DEM packing of hard spheres is scale-invariant, so at
fixed distribution *shape* it cannot distinguish `d₁₀` from `d₃₂` — matching on either is
the same experiment. Mechanism M1 (permeability, which reads `d₃₂`) lives in the
interstitial fluid, which this model does not have. **Every quantity in this figure is free
of M1 by construction rather than by matching.**

**What it shows.** The residual CV dependence survives matching, and it has a specific
shape:

- the **rattler fraction rises steeply** — grains drop out of the load-bearing network
- **Z over all grains falls** in step, because the rattlers dilute it
- **Z over the backbone barely moves**, and the backbone still spans at every CV

So polydispersity does not *restructure* the force network. It *removes grains from it*.
The network that remains has almost the same coordination and still percolates; it simply
contains a smaller fraction of the material.

**The threshold is a definition, and it is swept.** Identifying rattlers by iterated removal
of grains with fewer than `k` contacts is k-core decomposition of the contact graph. At
frictional mean `Z ≈ 4.9` the `k = 4` core is *empty* — the whole pack reads as rattlers —
while a single non-iterated pass removes only 15%. All of `k` = 2, 3, 4 are recorded;
`k = 2` (true rattlers, which make no friction-dependent stability assumption) is primary.
