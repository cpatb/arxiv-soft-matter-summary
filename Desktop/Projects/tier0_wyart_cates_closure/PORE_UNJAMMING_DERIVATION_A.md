# Deriving `tau_relax` for a Wyart-Cates frictional-fraction field

**Pore-scale drainage as the un-jamming clock**

Derivation A. Written 2026-08-08. Derived from first principles without consulting the
project's own notes; the hand-set placeholder values were used only in §8, after the
prediction was fixed.

---

## 0. Summary of the result

The un-jamming clock is set by **how fast interstitial fluid can be dragged through the
packing's pores into the space the frictional network must open in order to release**. It is
a poroelastic (Darcy) drainage time, not a contact-scale time:

```
tau_relax = eps_dil * eta_f * L_d^2 / ( k(phi, a) * sigma* )
k(phi, a) = (1 - phi)^3 a^2 / (45 phi^2)        [Carman-Kozeny]
          = (1 - phi) r_h^2 / 5,   r_h = (1 - phi) a / (3 phi)
```

Predicted: **1.5 s (silica), 0.013 s (cornstarch), ratio 111.** The empirical placeholders
are 2.0 s, 0.02 s, ratio 100.

The ratio is the load-bearing prediction and it contains no adjustable parameter: it follows
from `a`, `phi`, `eta_f` and `sigma*` alone. Every other mechanism I examined either gives
the ratio the **wrong sign** (lubrication, strain-control: silica un-jams *faster*) or is
absurd by orders of magnitude (Brownian, inertial).

The physical reason the smaller particle un-jams *more slowly* is a single line:
**permeability scales as `a^2`.** Silica's pores are 148 nm across, cornstarch's are 2.9 µm.
Same dilation, same driving pressure, 360× less permeable packing.

---

## 1. What "un-jamming" has to mean microscopically

In the model as written, `phi` is fixed (or advected) and `phiJ` moves because `fhat` moves.
That is a bookkeeping convenience; it hides the actual mechanics. A frictional contact is
frictional because two asperities are pressed together and are transmitting a tangential
load. To make it frictionless you must **physically separate the two surfaces** by at least
the asperity height `h_r`, so the contacts disengage and only lubrication films remain.

Separating every contact in a region by `h_r` is a **dilation**. The volumetric strain is

```
linear strain per contact = h_r / (2a)          (centre-to-centre distance is 2a)
eps_dil = 3 * h_r / (2a) = 1.5 * (h_r/a)
```

with `h_r/a ~ 0.01` for both systems (colloidal silica: 1-5 nm asperities on a 350 nm
sphere; cornstarch granules: 50-200 nm on a 7.5 µm granule), giving **`eps_dil ~ 0.015`**.

The suspension is incompressible and the drop's free surface is pinned by capillarity (§7),
so that new pore volume has to be **filled by liquid flowing in from elsewhere in the
sample**. That flow is resisted by the packing's permeability. *That* is the clock.

This is exactly the steer: the pore size controls both how fast a jammed region releases and
how long it survives.

---

## 2. Geometry: the hydraulic pore size, explicitly

Per unit total volume, void volume is `(1 - phi)` and wetted particle area is
`phi * (3/a)` (a sphere has surface/volume = 3/a). The hydraulic radius is their ratio:

```
r_h = (1 - phi) a / (3 phi)                                            [m]
```

Carman-Kozeny with Kozeny constant `k0 = 5` and specific surface `S_v = 3/a`:

```
k = eps^3 / (k0 S_v^2 (1-eps)^2)  with eps = 1 - phi
  = (1 - phi)^3 a^2 / (45 phi^2)
  = (1 - phi) r_h^2 / 5                     <-- pore size enters as r_h^2
```

(The second identity was checked numerically to machine precision.)

| | silica (a = 0.35 µm, phi = 0.440) | cornstarch (a = 7.5 µm, phi = 0.465) |
|---|---|---|
| hydraulic pore radius `r_h` | **148 nm** (= 0.424 a) | **2.88 µm** (= 0.384 a) |
| permeability `k` | 2.47e-15 m² (= 0.0202 a²) | 8.89e-13 m² (= 0.0158 a²) |
| Brinkman screening `sqrt(k)` | 49.7 nm | 943 nm |

Cornstarch's `phi` was not supplied; I set it by matching the distance to jamming,
`phi/phiM = 0.8494`, giving `phi = 0.4646`. Sensitivity is in §8.3.

### Why Carman-Kozeny, and is it good enough near jamming?

CK is derived for a bundle-of-tubes with hydraulic radius `r_h` and is known to be accurate
to roughly ±30% for random sphere packings over `0.3 < phi < 0.6`; it fails in the *dilute*
limit, not the dense one. I cross-checked against two independent closures:

| phi | Carman-Kozeny `k/a²` | Happel cell model | Beetstra-van der Hoef LB |
|---|---|---|---|
| 0.44 | 0.02016 | 0.02043 | 0.01081 |
| 0.518 | 0.00927 | 0.00988 | 0.00438 |
| 0.64 (RCP) | 0.00253 | 0.00282 | 0.00091 |

CK and Happel agree to ~7%; the lattice-Boltzmann correlation is ~2× lower but has the same
`a²` scaling and the same mild `phi` dependence. **The predicted silica/cornstarch ratio is
111, 113, 107 under the three laws respectively.** The result is not an artefact of CK.

One genuinely important fact: **permeability does not vanish at jamming.** At random close
packing `phi = 0.64` the porosity is still 0.36 and `k = 2.5e-3 a²`. There is no permeability
divergence at `phiJ`. This matters for §9.

---

## 3. (M1) Poroelastic / Darcy drainage: the derivation

### 3.1 Direct volume-balance derivation

A region of linear size `L_d` must dilate by volumetric strain `eps_dil`. Required liquid
volume:

```
dV = eps_dil * L_d^3
```

The liquid enters through the region's surface, area `~L_d^2`, driven by Darcy flow. The
driving pressure difference is the particle-phase (osmotic) pressure the network releases
when it lets go. In the Wyart-Cates framework the only intrinsic particle-phase stress scale
is the onset stress `sigma* = F*/a^2`, the repulsive force scale that has to be exceeded to
make a contact frictional in the first place. So `dP ~ sigma*`, and

```
Darcy flux         q     = (k / eta_f) * (dP / L_d)              [m/s]
volumetric rate    Q     = q * L_d^2 = k L_d dP / eta_f          [m^3/s]
time               tau   = dV / Q = eps_dil L_d^2 eta_f / (k dP)
```

```
================================================================
   tau_poro = eps_dil * eta_f * L_d^2 / ( k(phi,a) * sigma* )
================================================================
```

### 3.2 Cross-check as a Terzaghi/Biot consolidation time

Poroelastic consolidation gives `tau = L^2 / c_v` with diffusivity `c_v = k K / eta_f`,
`K` the drained constrained modulus. My expression is identical with

```
K = sigma* / eps_dil      (secant modulus: releasing stress sigma* over strain eps_dil)
```

which for `eps_dil = 0.015` gives `K = 2.9 kPa` (silica), `0.75 kPa` (cornstarch). Sanity
check by an independent route, an affine spring-network estimate. Elastic energy density under
isotropic strain `epsilon`, contacts compressing by `2a*epsilon`, `nZ/2` contacts per volume:

```
U = (1/2) k_c (2a eps)^2 (nZ/2) = k_c a^2 eps^2 n Z,  n = 3 phi/(4 pi a^3)
K = 2U/eps^2 = (3 phi Z / 2 pi) * k_c / a
```

with contact stiffness `k_c = F*/lambda = sigma* a^2 / lambda`, `lambda` the range of the
repulsion. Taking `Z = 4` (frictional isostatic), `phi = 0.44`, `lambda = h_r = 0.01a`:
`K = 0.84 * sigma* * (a/lambda) = 0.84 * 44 * 100 = 3.7 kPa` for silica. Same order as the
secant estimate. The two routes agree to within a factor 1.3, and I use the secant
form because it is the one tied to a measurable (`sigma*`) with no extra parameter.

**The L-free, geometry-free material prediction** is therefore the poroelastic diffusivity:

```
c_v = k sigma* / (eta_f eps_dil)
    = 1.54e-9  m^2/s   (silica)
    = 1.72e-7  m^2/s   (cornstarch)
```

This is what an experiment should actually measure (§12). `tau_relax = L_d^2 / c_v` then
follows for whatever drainage length the flow imposes.

### 3.3 The drainage length `L_d`

This is the one input the derivation cannot supply from material properties, and it is the
weakest link in the whole argument. It is *the size of the region that un-jams coherently*.

- **Lower bound, cluster branch.** If clusters release one at a time, `L_d = xi`, the
  force-network correlation length, `xi = a (1 - phi/phiJ)^-nu`. At the stated states this
  is 0.9-2.3 µm (silica) and 19-50 µm (cornstarch) for `nu = 0.5-1`. Because `xi ∝ a`,
  **this branch makes `tau` a-independent** and yields 0.5-3.5 ms (silica), 2-14 ms
  (cornstarch), a ratio of 0.24, wrong sign.
- **Upper bound, coherent branch.** If the release is driven by a *field*, the impact
  pressure decaying below `sigma*` across a whole region at once, then the entire region
  dilates simultaneously and `L_d` is the region's smallest dimension. In a spreading
  lamella that is the film thickness.

The physical situation the solver models is the coherent one: `sigma(x,t)` falls below
`sigma*` over a macroscopic patch, not cluster by cluster. So `L_d` = jammed-layer
thickness. Independent estimate for `D0 = 3 mm`, `U0 = 1-5 m/s`, water:

```
Re = 3000-15000
h_lamella ~ D0 Re^(-2/5)   = 122, 79, 64 um   (U0 = 1, 3, 5 m/s)
h_lamella ~ D0/(2 sqrt(Re)) = 27, 16, 12 um
```

so **`L_d ~ 15-120 µm`, geometric mean ≈ 45 µm.** I take `L_d = 50 µm` as the default. This
was fixed *before* looking at the placeholders; §8 shows what happens next.

### 3.4 Scaling exponents (M1)

```
tau_poro  ~  a^-2                          (through k ∝ a^2)
          ~  L_d^2
          ~  eta_f / sigma*
          ~  phi^2 (1 - phi)^-3            (through k)
          ~  h_r / a                       (through eps_dil)
```

Note `eps_dil ∝ h_r/a` is roughly `a`-independent if roughness scales with particle size, so
the net particle-size scaling stays `a^-2`.

---

## 4. (M2) Lubrication pull-apart at a single frictional contact

Two equal spheres, radius `a`, gap `h << a`, reduced radius `a/2`. Lubrication resistance to
separation along the line of centres:

```
F_lub = 6 pi eta_f (a/2)^2 (dh/dt)/h = (3/2) pi eta_f a^2 (dh/dt)/h
```

Once the confining load is removed, the restoring force is the stabilising repulsion, and
in Wyart-Cates this force is *definitionally* `F* = sigma* a^2`. (Brownian force `~kT/a` is
0.012 pN for silica vs `F* = 5.4` pN; gravity/buoyancy is ~1e-16 N; capillarity is absent in
bulk. The repulsion wins.) Balancing:

```
(3/2) pi eta_f a^2 (dh/dt)/h = sigma* a^2
dh/dt = [2 sigma* / (3 pi eta_f)] h            -> exponential opening
tau_e = 3 pi eta_f / (2 sigma*)                 <-- a CANCELS EXACTLY
tau_lub = tau_e * ln(h_end / h_0)
```

| | `tau_e` | `tau_lub` (ln ≈ 2.3) |
|---|---|---|
| silica | 5.03e-4 s | **1.16e-3 s** |
| cornstarch | 1.63e-3 s | **3.75e-3 s** |

**Ratio 0.31: silica un-jams three times *faster*. Wrong sign.**

This is not a quirk of my prefactors. It is a theorem about *any* purely local mechanism:
the only dimensionally-admissible time built from local quantities is `eta_f/sigma*` times a
dimensionless function of `phi` and `h_r/a`, and `a` drops out because the lubrication
resistance and the repulsive force both scale as `a^2`. **A contact-scale mechanism cannot
produce a size-dependent `tau_relax` at fixed `sigma*`.**

### 4.1 Van der Waals caution

Is the released silica contact actually repulsive? Sphere-sphere vdW with `A_H = 1e-20 J`,
`F_vdW = A_H a / (12 h^2)`:

| | `h` | `F_vdW` | `F* = sigma* a²` | `F*/F_vdW` |
|---|---|---|---|---|
| silica | 1 nm | 2.9e-10 N | 5.4e-12 N | **0.018** |
| silica | 5 nm | 1.2e-11 N | 5.4e-12 N | **0.46** |
| cornstarch | 75 nm | 1.1e-12 N | 6.4e-10 N | 572 |

For silica at contact, vdW attraction is comparable to or larger than the repulsion. If the
suspension sat in the primary minimum, contacts would be adhesive (`U_adh ~ A_H a/(12 h) ~
70 kT`) and would essentially never release. That *would* give the right sign, but it would
also make the suspension irreversibly aggregated and its thickening irreversible. Reversible
DST flow curves say it is not. So `sigma*` should be read as the *net* repulsion at the
stabiliser-set separation, vdW already folded in. **But if the silica is only marginally
stable, an adhesive contribution is a live alternative explanation for its slowness, and it
is not captured by any formula here.** See §13, assumption A7.

---

## 5. (M3) Brownian, and (M4) inertial/acoustic: both ruled out

### M3 Brownian

`D_0 = kT/(6 pi eta_f a)`; time to diffuse a distance `l` is `l^2/2D_0`.

| | `D_0` | diffuse `h_r` | diffuse `a` | Pe at 1e3 /s |
|---|---|---|---|---|
| silica | 1.33e-13 m²/s | 4.6e-5 s | 0.46 s | 9.2e2 |
| cornstarch | 7.5e-15 m²/s | 0.38 s | **3.8e3 s** | 7.5e6 |

Structural Brownian time `~ 6 pi eta_f a^3 / kT` scales as **`a^3`**, so the *bigger*
particle relaxes more slowly (wrong sign), and the cornstarch number (an hour) is
physically absurd. `Pe ~ 1e3-1e9` confirms neither system is Brownian-controlled during
impact. Note the near-coincidence that silica's structural time (0.46 s) is within 4× of the
2 s placeholder; treat that as coincidence, since the same formula fails cornstarch by 5
orders. Also, Brownian motion is a *parallel* release channel, so the fastest channel wins:
`min(0.46 s, 1.2 ms) = 1.2 ms`. **Ruled out.**

### M4 Inertia and acoustics

```
Stokes particle relaxation  tau_p = 2 rho_p a^2/(9 eta_f)
  silica     1.3e-8 s      cornstarch  4.8e-6 s
inertio-viscous  rho a^2/eta_f (silica)  2.6e-8 s
acoustic transit D0/c = 3 mm / 1500 m/s  2.0e-6 s
```

Five to eight orders of magnitude faster than anything else on the table, and `tau_p ∝ a^2`
also has the wrong sign. **Ruled out for completeness.** (Acoustics matters for the *impact
pressure pulse*, not for the un-jamming clock.)

### M5 Strain-controlled release (shear reversal), worth noting

Shear-reversal experiments show frictional contact networks are destroyed after a reverse
strain `gamma_c ~ 0.05-0.2`, independent of rate. That is a genuinely competitive closure:
`tau_relax = gamma_c/|gammadot|`. At impact rates `1e3-1e5 /s` it gives 1e-6 to 1e-3 s, and
it is **`a`-independent**, so again the wrong sign. It matters as an *upper bound while flow
continues*: contacts cannot survive longer than the time to accumulate `gamma_c` of reversed
strain. It says nothing about release at rest, which is the case `tau_relax` is meant to
cover.

---

## 6. Ranked verdict

| Rank | Mechanism | `tau` silica | `tau` corn | ratio | sign correct? |
|---|---|---|---|---|---|
| **1** | **M1 poroelastic drainage, coherent `L_d`** | **1.5 s** | **0.013 s** | **111** | **yes** |
| 2 | M2 lubrication pull-apart | 1.2e-3 s | 3.7e-3 s | 0.31 | no |
| 2b | M1 with cluster `L_d = xi` | 0.5-3.5e-3 s | 2-14e-3 s | 0.24 | no |
| 3 | M5 strain-controlled | 1e-6 - 1e-3 s | same | 1 | no (neutral) |
| 4 | M3 Brownian | 0.46 s | 3800 s | 1.2e-4 | no (inverted) |
| 5 | M4 inertial/acoustic | 1.3e-8 s | 4.8e-6 s | 2.7e-3 | no |

**M1 wins because it is the only mechanism with a length in it that does not scale with the
particle.** Everything else reduces to `eta_f/sigma*` times an `a`-independent number, and
`eta_f/sigma*` happens to be 3× *larger* for cornstarch, so every local mechanism predicts
the opposite of what is observed, by a factor 3, i.e. off by 300× on the ratio.

The others are subdominant for concrete reasons:
- **M2** is real, it just is not rate-limiting: it is the *last* step. A contact cannot open
  until the fluid to fill the gap has arrived. M2 (1 ms) and M1 (1.5 s) are in **series**,
  and M1 is 1000× slower for silica. For cornstarch they are within 4× (3.7 ms vs 13 ms), so
  cornstarch sits near the crossover, which is itself a testable statement.
- **M3** is in **parallel** and slower than M2 for silica, so it never controls.
- **M4** is nine orders too fast.

### What would falsify M1

1. **`tau_relax` independent of rheometer gap.** M1 predicts `tau ∝ L_d^2`; going from a
   100 µm to a 2 mm gap must lengthen the relaxation by ~400×. If the measured decay time is
   gap-independent, M1 is dead and M2 is right. **This is the decisive experiment and it is
   cheap.**
2. **Matching `eta_f/sigma*` between two suspensions of different `a` and finding equal
   `tau`.** M1 says `tau` should still differ as `a^-2`.
3. **A measured poroelastic diffusivity far from `k sigma*/(eta_f eps_dil)`** (§12, E2).
4. Silica turning out to be weakly aggregated (adhesive contacts), which would explain the
   slowness without any drainage physics.

---

## 7. Why capillarity does not rescue the free surface (and why it supports M1)

For a drop, one might hope the dilating network can simply pull air in, or push the surface
out, avoiding long-range drainage. It cannot. The capillary pressure at the pore scale is

```
P_cap = 2 gamma / r_h,   gamma = 0.072 N/m
      = 9.7e5 Pa (silica, r_h = 148 nm)
      = 5.0e4 Pa (cornstarch, r_h = 2.9 um)
```

versus a particle-phase driving pressure of only `sigma* =` 44 and 11.3 Pa,
**`sigma*/P_cap = 4.5e-5` and `2.3e-4`**. Air cannot invade the pores by four orders of
magnitude. The drop is therefore a **closed** poroelastic domain: liquid dilating one region
must be squeezed out of another region inside the same drop. That is precisely what forces
`L_d` to be a flow-scale length rather than a pore-scale one, and it is the second half of
the physical steer. The pore size is also what keeps the region jammed, by making the free
surface an impenetrable wall for the interstitial fluid.

(This is the same capillary pressure that governs granulation and the "cornstarch
monster"; it scales as `1/a`, so it is 20× more restrictive for silica.)

---

## 8. Numbers

### 8.1 Central predictions (`L_d = 50 µm`, `eps_dil = 0.015`, `h_r/a = 0.01`)

| | silica | cornstarch |
|---|---|---|
| `a` assumed | 0.35 µm (mid of 0.2-0.5) | 7.5 µm (mid of 5-10) |
| `phi` | 0.440 | 0.4646 (assumed, `phi/phiM` matched) |
| `k` | 2.47e-15 m² | 8.89e-13 m² |
| `c_v = k sigma*/(eta_f eps_dil)` | 1.54e-9 m²/s | 1.72e-7 m²/s |
| **`tau_relax`** | **1.5 s** | **0.0134 s** |
| **ratio** | **111** | |

### 8.2 Honest uncertainty

Two separable pieces.

**The ratio** depends only on `(eta_f/sigma*)` and `k`, both measured or geometric:

```
tau_s/tau_c = (eta_s/eta_c) (sigma_c*/sigma_s*) (k_c/k_s) = 1.205 * 0.257 * 360 = 111
```

Across permeability laws: 107-113. Across the stated particle-size brackets and `sigma*`
brackets taken to their opposing extremes: **7.5 to 1180**, central 111, with a realistic
band of roughly **40-300**. The `(a_c/a_s)^2 = 459` size factor is doing the work; the
`sigma*` ratio partially cancels it.

**The absolute magnitude** is dominated by `L_d^2`, uncertain by ~1e3:

| `L_d` | silica | cornstarch |
|---|---|---|
| 1 µm | 6.5e-4 s | 5.8e-6 s |
| 16 µm | 0.17 s | 1.5e-3 s |
| **50 µm** | **1.5 s** | **0.013 s** |
| 120 µm | 9.3 s | 0.084 s |
| 1.5 mm (drop) | 1.5e3 s | 13 s |

Material-parameter sensitivity at fixed `L_d = 50 µm`:
- silica, `a ∈ [0.2, 0.5] µm` and `sigma* ∈ [40, 100] Pa` → **0.32 - 5.0 s**
- cornstarch, `a ∈ [5, 10] µm` and `sigma* ∈ [8, 20] Pa` → **0.0043 - 0.043 s**

So: **silica `tau_relax` = 1.5 s, credible 0.3-5 s. Cornstarch = 0.013 s, credible
0.004-0.04 s.** Both assume `L_d ≈ 50 µm`; both scale as `L_d^2`.

### 8.3 The `phiM = 0.438` convention for cornstarch

If cornstarch's `phiM` is 0.438 rather than 0.547 and `phi` is rescaled to keep
`phi/phiM = 0.8494`, then `phi = 0.372`, `k = 2.24e-12 m²` (2.5× higher), and
`tau_corn = 5.3e-3 s`, pushing the ratio to **280**. This is the single largest discrete
sensitivity in the calculation and should be pinned down before the closure is trusted
quantitatively.

### 8.4 Comparison with the placeholders (done last)

| | predicted | placeholder | pred/placeholder |
|---|---|---|---|
| silica | 1.50 s | 2.0 s | 0.75 |
| cornstarch | 0.0134 s | 0.02 s | 0.67 |
| **ratio** | **111** | **100** | **1.11** |

Both are reproduced to within 33% with a *single* drainage length that was chosen from
lamella hydrodynamics, not fitted. Inverting instead, asking what `L_d` reproduces silica's
2.0 s exactly, gives `L_d = 55.5 µm`, which then predicts cornstarch at **0.0180 s** against
a placeholder of 0.020 s. **One parameter, fixed by one datum, predicts the second datum to
10% across a factor of 100.** That is the strongest evidence in this document.

### 8.5 Bonus: the same closure derives the loading/unloading asymmetry

The model already assumes `tau_up < tau_relax` and gets no justification for it. Compaction
(loading) is the same Darcy problem run backwards, driven by the *full applied* stress rather
than by the released repulsive pressure:

```
tau_up / tau_relax = sigma* / sigma_applied
```

With `sigma_applied ~ rho U0^2` (1-25 kPa for 1-5 m/s), this gives `tau_up ≈ 2.6-66 ms`
(silica) and `6 µs-0.15 ms` (cornstarch). **The asymmetry is therefore not an extra assumption. It is `sigma*/sigma`,
and it is automatically large under impact.** Worth adopting.

---

## 9. Behaviour as `phi -> phiJ`

**`tau_relax` does not diverge through the material terms.** This is the most practically
important finding for the solver.

```
tau ∝ phi^2 (1 - phi)^-3
```

| `phi` | `r_h/a` | `k/a²` | `tau/tau(0.44)` |
|---|---|---|---|
| 0.30 | 0.778 | 0.0847 | 0.24 |
| 0.44 | 0.424 | 0.0202 | 1.00 |
| 0.518 (`phiM`) | 0.310 | 0.00927 | 2.17 |
| 0.60 | 0.222 | 0.00395 | 5.10 |
| 0.64 (RCP) | 0.188 | 0.00253 | 7.96 |

Going from the operating point all the way to random close packing lengthens `tau_relax` by
only **8×**, smoothly, with no singularity. The reason is physical: at RCP the porosity is
still 0.36 and the pores are still open. Jamming stops the *particles* moving; it does not
close the *pores*. A closure of this form is numerically benign, with no stiffness blow-up as the
solver approaches or crosses `phiM`.

**The divergence, if any, lives in `L_d`, and it is geometric, not material.** When the
frictional network percolates (`phi > phiM` with `fhat -> 1`), the coherently-releasing
region jumps from cluster-scale to system-scale. Writing `L_d = min(xi, L_flow)` with
`xi = a (1 - phi/phiJ(fhat))^-nu` gives

```
tau ∝ (1 - phi/phiJ)^(-2 nu)     for xi < L_flow      (nu ~ 0.5-1, so exponent -1 to -2)
tau ∝ L_flow^2                    saturated, for xi > L_flow
```

For silica, `xi` reaches `L_flow = 50 µm` only at `1 - phi/phiJ ≈ 0.007`, i.e. within 1% of
jamming; for cornstarch it is already there at the operating point (`xi = 19-50 µm`). This is
an internal tension in the derivation and I flag it as such: at `phi = 0.44 < phiM = 0.518`
the material is thickened (`eta/eta_s = 44`) but **not jammed**, so the cluster branch is
arguably the correct one there, and the cluster branch gives milliseconds with the wrong-sign
ratio. The coherent branch is justified instead by the *mode of release*, i.e. the impact
pressure decaying below `sigma*` across a whole layer at once, not by percolation. See §13,
assumption A3, which is the assumption most likely to be wrong.

---

## 10. Dimensionless group against the impact

```
De = tau_relax / (D0/U0),   D0 = 3 mm
```

| `U0` | `t_imp` | silica `De` (M1) | cornstarch `De` (M1) | silica `De` (M2) | corn `De` (M2) |
|---|---|---|---|---|---|
| 1 m/s | 3.0 ms | 498 | 4.5 | 0.39 | 1.25 |
| 3 m/s | 1.0 ms | 1495 | 13.4 | 1.16 | 3.7 |
| 5 m/s | 0.6 ms | 2492 | 22.4 | 1.93 | 6.2 |

**Silica: `De ~ 5e2 - 2.5e3`.** The frictional network cannot relax at all during an impact.
Practically, the solver can treat silica un-jamming as **frozen**: set `tau_relax = inf`, or
equivalently latch `fhat` once it rises, and check that results are insensitive to any
`tau_relax > ~0.1 s`. This is a real simplification and removes the parameter's influence
entirely.

**Cornstarch: `De ~ 4-22`.** Slow compared with the impact itself, but *not* compared with
retraction and oscillation (5-50 ms). So cornstarch will partially release during retraction
and `tau_relax` genuinely matters for the late-time drop shape. This is where the closure
earns its keep, and where an error of 3× is visible in the output.

A useful corollary: because `De >> 1` for both, **the observed behaviour is controlled by
`tau_up`, not `tau_relax`**, which strengthens the case for adopting §8.5's derived
`tau_up = tau_relax * sigma*/sigma`.

---

## 11. Recommended closure, ready to code

```python
# ---------------------------------------------------------------- units -----
# a       [m]     particle radius
# phi     [-]     local particle volume fraction
# eta_f   [Pa s]  suspending-fluid viscosity
# sigma_s [Pa]    Wyart-Cates onset stress sigma*
# h_r     [m]     asperity height / range of the stabilising repulsion
# L_d     [m]     drainage length = smallest dimension of the coherently
#                 un-jamming region (lamella/jammed-layer thickness)
# returns [s]

def tau_relax(a, phi, eta_f, sigma_s, h_r=None, L_d=50e-6, phi_max=0.639):
    if h_r is None:
        h_r = 0.01 * a                      # roughness ~ 1% of radius
    phi = min(phi, phi_max)                 # clamp: CK is meaningless above RCP
    eps_dil = 1.5 * h_r / a                 # dilation needed to open all contacts
    k = (1.0 - phi)**3 * a**2 / (45.0 * phi**2)      # Carman-Kozeny [m^2]
    return eps_dil * eta_f * L_d**2 / (k * sigma_s)

def tau_up(a, phi, eta_f, sigma_s, sigma_applied, **kw):
    """Loading is the same Darcy problem driven by the full applied stress."""
    return tau_relax(a, phi, eta_f, sigma_s, **kw) * sigma_s / max(sigma_applied, sigma_s)
```

Equivalently, and preferably, store the **material** number and let geometry supply the rest:

```
c_v = k(phi,a) * sigma* / (eta_f * eps_dil)      [m^2/s]     tau_relax = L_d^2 / c_v
```

Three warnings for the implementation:

1. **`L_d` must not be the grid spacing.** If someone wires `L_d = dx`, every result becomes
   mesh-dependent and will "converge" to zero relaxation time under refinement. Set it from
   the physical jammed-layer thickness and hold it fixed under refinement, or make it a
   resolved field via a genuine two-phase model (below).
2. **Clamp `phi`** below RCP; Carman-Kozeny is meaningless at zero porosity.
3. **The right long-run fix is not a scalar `tau_relax` at all.** This derivation says the
   un-jamming clock is a *drainage* time, which means the honest model is a two-phase
   (suspension-balance / Darcy) formulation carrying a relative fluid flux, with `phi`
   advected and dilatancy resolved. A scalar `tau_relax` is a lumped stand-in for that, and
   its geometry-dependence is the symptom. If the drop-impact results turn out sensitive to
   `tau_relax`, that is the signal to upgrade rather than to re-tune.

---

## 12. Two experiments that would measure `tau_relax` directly

### E1. Gap-resolved stress-drop rheometry: decisive and cheap

Load the suspension above `sigma*` to a steady thickened state, then step the stress down to
`sigma_low << sigma*` and record `eta(t)` (or `gammadot(t)`) as it decays. The decay time is
`tau_relax` by construction of the model.

The essential step, and the reason to do this rather than a plain stress-jump: **repeat at
several gaps, 100 µm to 2 mm.** M1 predicts `tau ∝ gap^2`, a 400× change; M2 predicts no
change at all. Nothing else distinguishes them so cleanly.

- Effort: **2-4 days** on a stress-controlled rheometer with roughened parallel plates, plus
  a day of analysis. Cornstarch first: its predicted 4-40 ms is comfortably resolvable.
- Pitfall: instrument inertia. The rheometer's own transient is 1-10 ms and will swamp the
  cornstarch signal at small gaps. Use a low-inertia head, subtract the measured instrument
  response on a Newtonian oil of matched viscosity, and confirm the fitted time exceeds the
  instrument time by 5×. Silica's predicted 0.3-5 s is far outside instrument time and is the
  safer first measurement despite being the less interesting fluid.
- Second pitfall: wall slip and shear banding will mimic slow relaxation. Roughen, and verify
  with two gaps at fixed ratio of expected `tau`.

### E2. Consolidation (oedometer) cell: measures `c_v` directly, geometry-free

Confine the suspension in a filtration cell between two porous frits with a fast pressure
transducer at the base. Apply a step load, let it consolidate, then **release** the load and
record the pore-pressure transient. Standard Terzaghi analysis returns the poroelastic
diffusivity `c_v` directly, the quantity in §3.2, with no reference to any drainage length,
and `tau_relax = L_d^2/c_v` follows for whatever geometry the solver needs.

Predicted values to test against: `c_v = 1.5e-9 m²/s` (silica), `1.7e-7 m²/s` (cornstarch).

- Effort: **~2 weeks**, mostly cell fabrication. The measurement itself is textbook soil
  mechanics. A 5 mm sample column gives silica time constants of ~4 hours and cornstarch of
  ~2 minutes; shorten the silica column to 1 mm (~10 min) or accept the long run.
- Bonus: the same cell measures `k(phi)` independently by steady permeation, which removes
  the Carman-Kozeny assumption entirely and turns the closure into a fully measured object.

*(A third, if instrumentation allows: ultrasound or high-speed X-ray radiography of a rod
withdrawn from a jammed bath, timing the re-fluidisation front. Directly visualises `L_d`.
Effort: months.)*

---

## 13. Every assumption that could be wrong

**Structural, ranked by how much damage they do:**

- **A1. `dP = sigma*`.** I take the pressure released on un-jamming to be the onset stress.
  If instead the relevant pressure is the *residual confining stress* left over from the
  impact (which can be `rho U0^2 ~ 1e4` Pa, 200× larger), `tau_relax` drops by that factor
  and the drop's history enters the closure. **Linear in the answer.**
- **A2. `eps_dil = 1.5 h_r/a ≈ 0.015`, and `h_r/a ≈ 0.01` for both systems.** Roughness is
  the least well known input. If silica is smooth (`h_r/a = 0.003`) and cornstarch rough
  (`h_r/a = 0.03`), the ratio falls from 111 to ~11. **Linear in the ratio.**
- **A3. `L_d` = coherent flow-scale (~50 µm), not the cluster correlation length (~1-20 µm).**
  This is the assumption I trust least; §9 lays out the tension. If the cluster branch is
  right, `tau_relax` is milliseconds for both and the placeholder ratio is *backwards*.
  **Quadratic in the answer, and it flips the sign of the conclusion.** E1 tests it directly.
- **A4. Cornstarch `phi` = 0.4646**, inferred by matching `phi/phiM`. Not given. Together with
  the `phiM = 0.547` vs `0.438` convention this is a 2.5× swing (§8.3).
- **A5. Carman-Kozeny.** Cross-checked against Happel and Beetstra (§2), ratio robust to
  ±6%, absolute to ~2×. Assumes monodisperse spheres; cornstarch granules are polydisperse
  and non-spherical, and polydispersity *reduces* permeability (fines plug throats), so
  cornstarch's `k` is likely overestimated and `tau_corn` underestimated.
- **A6. `sigma* = F*/a^2` with prefactor 1.** Different papers carry factors of `6 pi` or
  `1/pi`. This shifts M2 (§4) but cancels out of M1's *ratio*, since `sigma*` is used as a
  measured stress there.

**Physical:**

- **A7. Silica contacts are net-repulsive on release.** §4.1 shows vdW is comparable to `F*`
  at 1-5 nm. If silica is weakly aggregated, adhesion rather than drainage explains its slowness,
  and the whole `a^-2` story is a coincidence. Test: does the silica flow curve show
  reversible thickening with no yield stress?
- **A8. The dilation is isotropic and affine.** Real dilatancy is anisotropic, concentrated
  in shear bands; a band drains over its *thickness*, which is smaller than my `L_d`.
- **A9. Fluid is incompressible and single-phase; no dissolved air.** Entrained microbubbles
  would provide compliance and short-circuit the drainage entirely.
- **A10. `sigma*` is stress-history independent** and the same on loading and unloading. If
  contacts have any adhesive hysteresis, unloading `sigma*` is lower and `tau` longer.
- **A11. Isothermal, no evaporation, no sedimentation** over the relaxation time. For silica
  at `tau = 2 s` sedimentation is negligible (Stokes velocity ~0.1 µm/s); for cornstarch at
  0.013 s also negligible. Safe.
- **A12. `nu` for the correlation length** is 0.5-1 by analogy with isostatic length scaling;
  I have not verified the exponent for frictional shear jamming.

**Numerical / methodological:**

- **A13.** Both particle radii are stated as ranges and I used midpoints. `tau ∝ a^-2` makes
  this a 6.25× lever on silica and 4× on cornstarch (§8.2).
- **A14.** The empirical targets are hand-guessed placeholders, so agreement to 33% is *not*
  evidence of accuracy. The defensible claim is narrower: **a mechanism with `tau ∝ a^-2`
  reproduces a ratio of ~100 with the smaller particle slower, and no other mechanism
  examined can produce that sign at all.**

---

## 14. Bottom line for the solver

1. Adopt `tau_relax = eps_dil eta_f L_d^2 / (k(phi,a) sigma*)` with `L_d` the jammed-layer
   thickness, and record `L_d` explicitly as a model input rather than burying it.
2. Adopt `tau_up = tau_relax * sigma*/sigma_applied`. The asymmetry then costs no new
   parameter.
3. Because `De >> 1` for silica, **verify the results are insensitive to `tau_relax`** over
   0.1-10 s. If they are, the parameter never mattered and this whole question is moot for
   silica; if they are not, that is surprising and worth understanding.
4. Run E1 (§12) before trusting the absolute magnitude. The gap-scaling test costs days and
   settles the mechanism.
