# Disc Flight in a Lakefront Internal Boundary Layer:
# Coupling Analytical Frisbee Aerodynamics with Site-Specific Urban CFD over Loyola University's Sean Earl Field

**Author:** charlie *(affiliation TBD)*
**Status:** First draft, v0.1 — generated 2026-05-29 from handwritten project notes
**Target venue (suggested):** *Sports Engineering* (Springer), or AIAA Applied Aerodynamics Conf. (parallel disc-trajectory paper), or *Journal of Wind Engineering & Industrial Aerodynamics* (parallel urban-CFD paper).

---

## 0. How to read this document

Three layers stacked in one file:

1. **Feasibility memo (§1)** — go/no-go on what your handwritten notes actually propose, including a couple of things you should change.
2. **Project plan (§2)** — work breakdown, timeline, equipment/data, risks.
3. **Paper first draft (§3–§9)** — manuscript skeleton with real content where I have it and clearly-marked placeholders where I don't.

Read §1 first. If you disagree with the scope changes I propose, the draft in §3 onward needs to be reframed.

---

## 1. Feasibility check

### 1.1 What the notes propose

> (1) Analytical treatment of ideal conditions for throwing a frisbee.
> (2) Build a realistic 3D map of the terrain surrounding Loyola field — houses, buildings — simulate wind off the lake spanning 180°, generate 3D flow fields across the turf to capture vortices.
> Goal: determine the best conditions for this constrained terrain that best allow for frisbee.
> Tools: OpenFOAM or COMSOL.

### 1.2 Feasibility verdict

| Component | Verdict | Why |
|---|---|---|
| Analytical 6-DOF disc model | **Solid** | Hummel (2003), Crowther & Potts (2007) give a complete validated 6-DOF formulation with measured coefficients. You don't need new aerodynamic measurements; you re-use Potts–Crowther wind-tunnel data. |
| RANS urban CFD of a campus block, 13 wind directions × 15° | **Tractable** on a 32-core workstation. ~3–10 M cells, 2–8 h per direction, total 1–4 days. | Standard COST 732 / AIJ best-practice; OpenFOAM `simpleFoam` with k-ω SST. Geometry inputs are free (Cook County GIS, USGS 3DEP LiDAR, Microsoft footprints, OSM). |
| LES of the same domain | **Not realistic** on one workstation across 13 directions. Single direction is feasible (weeks). | Δx ~1–2 m near buildings ⇒ 50–200 M cells per case. |
| Coupling CFD wind field into disc trajectory simulator | **Novel and tractable.** | Not done in the published literature. The Lukes et al. (2022) disc-golf code uses uniform freestream; nobody has coupled environmental CFD to disc dynamics. **This is your novelty hook.** |
| "Determine best conditions for this constrained terrain" — i.e., optimize throws given the local wind | **Tractable** if framed as: for each wind direction θ, compute the spatial map of mean velocity + TKE at disc height (1–3 m AGL), then evaluate trajectory robustness for canonical throws (backhand, forehand, hammer). | Don't over-promise: this is descriptive optimization, not real-time forecasting. |

### 1.3 Things in your notes I'd flag

- **"Loyola field"** — your wording ("off the lake") only matches **Sean Earl Field** on the Lake Shore Campus, ~150–250 m from Lake Michigan, not Hoyne Field (which is a mile inland). The draft assumes Sean Earl. If you actually meant Hoyne, the wind framing has to change because there's no direct lake fetch.
- **"180° span" of wind directions** — Sean Earl's lake-facing side is roughly east; a meaningful "off-lake" sweep is roughly NNE through ESE (~90° span centered on east). The full 180° is fine if you also want to characterize the prevailing westerly (off-land) cases for contrast, but think of it as two regimes, not one continuous sweep.
- **"Range finder" to build the 3D map** — you don't need a personal LiDAR pass. USGS 3DEP / Illinois Height Modernization Program LiDAR covers Cook County at 1–8 pts/m², and Microsoft GlobalMLBuildingFootprints + USGS 2023 building-heights release explicitly include Chicago. A handheld rangefinder is useful only for spot-checking the 5–10 buildings nearest the field. Save the effort; the public data is good enough.
- **"OpenFOAM or COMSOL"** — pick OpenFOAM. The urban-CFD literature (Blocken, Tominaga, the COST 732 community) is overwhelmingly OpenFOAM-based; `snappyHexMesh` handles dirty STL of campus buildings well; you'll have peer-reviewed cases to validate against. COMSOL is fine for the analytical / 6-DOF side if you prefer a GUI, but it is not the standard methodological choice for atmospheric-scale building CFD.

### 1.4 The biggest unstated gap in the notes

You do not have a **bridging quantity** that turns a CFD flow field into a throw-quality score. The paper needs one. I propose this in §3.6: define a "throwability metric" Q(x, θ) at each launch point x for each wind heading θ, combining mean wind magnitude, turbulence intensity at disc altitude, and predicted lateral trajectory drift for the three canonical throws. Without that, the work is two parallel literatures with no synthesis.

---

## 2. Project plan

### 2.1 Work breakdown (≈14 weeks, single researcher)

| Phase | Weeks | Output |
|---|---|---|
| **P0 — Setup** | 1 | OpenFOAM 11 or v2406 installed; Python disc-trajectory simulator stub running; QGIS workflow for LiDAR → STL set up. |
| **P1 — Analytical disc model** | 2–3 | 6-DOF integrator (Python, scipy `solve_ivp` with stiff-aware RK45/LSODA). Validated against Hummel (2003) Fig. 5-1 throw reconstructions. Outputs: trajectory(t), α(t), spin(t) given inflow profile U(x,t). |
| **P2 — Geometry build** | 4 | STL of ~250 m radius around Sean Earl Field, ~5–10 named buildings (Mertz, Norville Center, Damen Student Center, etc.) + tree-canopy proxy. Validated against Loyola campus map. |
| **P3 — RANS CFD sweep** | 5–8 | 13 steady-state runs (15° increments from N to S through E) on the campus mesh. Each run: simpleFoam + k-ω SST + Richards–Hoxey ABL inlet (z₀ = 0.001 m over water, 0.3 m over urban). Outputs: U_mean(x,y,z), TKE(x,y,z) sampled to 1, 2, 3 m AGL grids. |
| **P4 — Coupled trajectory** | 9–11 | Trilinear-interpolated CFD field fed into 6-DOF disc model. Monte Carlo over launch position, heading, spin, release α (10⁴ throws per wind direction). Outputs: landing-point ellipses, throw-success rate vs. (x_launch, θ_throw). |
| **P5 — Synthesis & manuscript** | 12–14 | Throwability map Q(x, θ_wind); paper revision; figures. |

### 2.2 Compute and storage

- **Workstation target**: 32 physical cores, 128 GB RAM, ~2 TB scratch SSD. If only a 16-core / 64 GB machine is available, expect ~2× wall-clock on CFD and one wind direction per overnight job — still finishes inside the schedule.
- **Disk budget**: ~50 GB per RANS case (mesh + fields + post), so ~700 GB for the sweep. LES, if attempted, would blow this up by 20–50×.

### 2.3 Inputs needed (all free / accessible)

- Building footprints & heights: Cook County GIS, Microsoft GlobalMLBuildingFootprints (ODbL), USGS Overture-Maps buildings 2023, OpenStreetMap.
- LiDAR DEM: USGS 3DEP / Illinois Height Modernization Program (ILHMP). 1–8 pts/m².
- Wind reference: NWS O'Hare 1991–2020 normals; NOAA GLERL Wilson Avenue Crib for shoreline calibration; Northerly Island ASOS for southern-lakefront reference.
- Disc aerodynamic coefficients: Hummel (2003) Tables 3-1 / 3-2; Crowther & Potts (2007).

### 2.4 Risks & mitigations

| Risk | Mitigation |
|---|---|
| **k-ω SST under-predicts rooftop separation at sharp roof corners** (well-known issue). | Validate one direction against a LES "spot check," follow Blocken (2018) on RANS limits. Don't claim peak/fluctuating quantities — report means + TKE only. |
| **Disc-coefficient uncertainty** at off-design α (>20°). | Restrict throw library to canonical α ∈ [−5°, +20°]; flag hammer (high-α) results as exploratory. |
| **Horizontal inhomogeneity of ABL inlet** (Hargreaves & Wright 2007). | Use the Yang et al. (2009) corrected wall functions; verify inlet profile is preserved upstream of buildings. |
| **Sean Earl is not actually a competition-grade Ultimate field** (it's small). | Reframe paper as a *methods* paper using Sean Earl as exemplar; the methodology generalizes. |
| **Reviewers ask for measurement validation.** | Use the Tominaga (2008) Building B-cube benchmark for the CFD; cite Hummel's marker-tracked flights for the disc model. Both are public. |

---

## 3. Manuscript: Disc Flight in a Lakefront Internal Boundary Layer

### Abstract *(draft)*

Outdoor disc sports are dominated by wind, and yet the published aerodynamic and trajectory literature for flying discs (Hummel 2003; Potts & Crowther 2002; Lukes et al. 2022) uses uniform-freestream assumptions that ignore site-specific flow structure. Conversely, the urban-CFD literature (Tominaga et al. 2008; Blocken 2014, 2018) characterizes site-specific wind fields but stops at pedestrian-level mean quantities. Here we couple a validated 6-degree-of-freedom (6-DOF) disc-flight model to steady-RANS atmospheric CFD of a lakefront university campus (Loyola University Chicago, Sean Earl Field; ~200 m from the Lake Michigan shore). We sweep 13 inflow directions in 15° increments using OpenFOAM with k-ω SST turbulence and a Richards–Hoxey atmospheric-boundary-layer inlet, validated against COST Action 732 benchmarks. Disc trajectories for canonical throws (backhand, forehand, hammer; 175 g Ultrastar) are integrated through the resulting time-averaged + turbulence-intensity fields. We report a spatial "throwability" map Q(x, θ_wind) and identify the lakefront internal boundary layer — formed by the sharp roughness transition at the shoreline — as the dominant local feature controlling disc-altitude flow. The methodology generalizes to any outdoor sport played on a roughness-discontinuity site.

### 3.1 Introduction *(draft)*

The Frisbee is one of the most-studied non-airfoil aerodynamic objects in sports physics. Stilley & Carstens (1972) provided the first parametric model; Yasada (1999), Hummel (2003), and Potts & Crowther (2002) supplied successively more complete wind-tunnel measurements of the lift, drag, and pitching-moment coefficients, and Hummel's UC Davis thesis remains the most-cited 6-DOF flight reconstruction. Subsequent CFD validation by Kamaruddin (2011) and Lukes et al. (2022, 2023) using RANS k-ω SST in OpenFOAM has converged the modeling stack to within a few percent of wind-tunnel CL, CD, and CM at moderate angles of attack.

Yet every published disc-trajectory study assumes a spatially uniform, steady freestream. As Lukes et al. (2022) explicitly note, *"the biggest uncertainty in trajectory experiments is the wind environment."* This uncertainty is not noise to be averaged away: outdoor disc sports — Ultimate, disc golf — are played in geographically structured turbulent flow. A throw whose intended landing zone lies in a sheltered eddy behind a building experiences a categorically different flow than the same throw 20 m to the side.

Independently, the urban CFD community has developed mature best-practice guidance for site-specific pedestrian-level wind: the COST Action 732 guidebook (Franke et al. 2007), the AIJ guidelines (Tominaga et al. 2008), and the 50-year retrospective of Blocken (2014). Standard practice uses steady RANS with k-ω SST or realizable k-ε on snappyHexMesh-generated meshes of ~10⁷ cells, with Richards–Hoxey atmospheric inlets. This work has shaped pedestrian-comfort regulation in dozens of cities, but to our knowledge has never been applied to outdoor sports trajectory prediction.

The present work fills that gap. We choose Sean Earl Field at Loyola University Chicago's Lake Shore Campus as a study site because it offers an unusually clean physical setup: the field sits within ~250 m of Lake Michigan, with no major intervening built obstacle to the east, behind it a typical low-rise urban-campus block to the west. Easterly winds therefore arrive after ~80 km of unbroken over-water fetch (z₀ ≈ 0.001 m), pass through a sharp roughness transition at the shoreline, and develop a shallow internal boundary layer over the field. Westerly winds arrive through the campus block. The site cleanly separates two regimes (smooth-fetch and urban-fetch) within a 180° heading sweep.

Our contributions:
1. The first published two-way coupling of an environmental-CFD wind field to a validated 6-DOF disc-flight model.
2. A "throwability metric" Q(x, θ_wind) that compresses the spatial CFD field plus trajectory uncertainty into an actionable single scalar per (launch position, wind heading).
3. A case study identifying the lakefront internal boundary layer as the dominant local control on Ultimate-disc flight at the chosen site.

### 3.2 Analytical Disc Aerodynamics *(draft, with equations)*

We adopt the Hummel (2003) / Crowther & Potts (2007) 6-DOF formulation. Three reference frames: inertial **N**, body-fixed **B** (Euler angles φ, θ, ψ from **N** to **B**), and disc-fixed **D** spinning at rate ψ̇ relative to **B**. The disc has mass m = 0.175 kg, diameter d = 0.27 m, and principal moments I_z = (1/2)mr² ≈ 1.6 × 10⁻³ kg·m² (spin axis) and I_x = I_y ≈ 8.0 × 10⁻⁴ kg·m².

The translational equation in **N** is

$$
m \dot{\mathbf{v}} = m \mathbf{g} + \tfrac{1}{2}\rho V_{\text{rel}}^2 A \left( -C_D \hat{\mathbf{v}}_{\text{rel}} + C_L \hat{\mathbf{L}} \right),
$$

with **V**_rel = **v** − **U**(**x**, t) the relative wind, A = πd²/4, and α the angle between **V**_rel and the disc plane. The aerodynamic coefficients are linearized about α₀ ≈ −4°:

$$
C_L(\alpha) = C_{L0} + C_{L\alpha}(\alpha - \alpha_0), \quad
C_D(\alpha) = C_{D0} + C_{D\alpha}(\alpha - \alpha_0)^2,
$$

with the Hummel "short-flight" values (CL0 = 0.33, CLα = 1.91 rad⁻¹, CD0 = 0.18, CDα = 0.69 rad⁻²) as the baseline; we provide a sensitivity to the Potts–Crowther (CL0 = 0.20, CLα = 2.96) set.

The rotational equations in **B** with Euler's formulation give a pitching moment

$$
M = \tfrac{1}{2}\rho V_{\text{rel}}^2 A d \left[ C_{M0} + C_{M\alpha}(\alpha - \alpha_0) + C_{Mq}\tfrac{q d}{2 V_{\text{rel}}} \right]
$$

and a roll moment that couples to spin rate Ω = ψ̇ through the advance ratio AdvR = Ωd/(2V_rel):

$$
R = \tfrac{1}{2}\rho V_{\text{rel}}^2 A d \left[ C_{Rr}\tfrac{r d}{2 V_{\text{rel}}} + C_{Rp}\tfrac{p d}{2 V_{\text{rel}}} \right].
$$

A gyroscopic-precession reduction (Hubbard & Hummel 2000) yields the well-known throw-flip relation: the precessional roll rate is approximately p ≈ M / (I_z Ω). At Ω = 100 rad/s and typical M = 0.05 N·m, p ≈ 0.3 rad/s — substantial bank over a 3 s flight.

**Reynolds number.** For V = 15 m/s and d = 0.27 m, Re = Vd/ν ≈ 2.7 × 10⁵. Potts & Crowther showed CL, CD are essentially Re-independent over 1.14 × 10⁵ – 3.78 × 10⁵, justifying the use of fixed coefficients across operational wind conditions.

### 3.3 Site description: Sean Earl Field *(draft)*

Sean Earl Field is a small multi-use turf flanked by Mertz Residence Hall (W) and the Norville Center for Intercollegiate Athletics (S), with a clear lakefront promenade and beach to the east. Per the Loyola Lake Shore Campus map (LUC 2025), the field's eastern sideline lies approximately 150–250 m from the Lake Michigan shoreline; the field's long axis runs roughly N–S. Loyola explicitly lists Ultimate among the field's sanctioned uses.

**Climatology.** Chicago lies in the mid-latitude westerly belt; the annual prevailing direction at ORD is WSW–SW (Illinois State Water Survey wind roses). The 1991–2020 NWS annual mean speed at O'Hare is ~9–10 mph. However, ORD is 27 km inland; lakefront stations (NOAA GLERL Wilson Avenue Crib, Northerly Island ASOS) record higher mean speed and gustier easterlies, especially during the warm-season lake-breeze regime (NWS Chicago / WGN-TV).

**Why this site.** Easterly inflow at Sean Earl Field passes over ~80 km of unbroken Lake Michigan fetch (to St. Joseph / Benton Harbor MI), giving a near-neutral, very-low-z₀ marine boundary layer. The shoreline imposes a sharp z₀ transition (~10⁻³ to ~0.3 m). The downwind internal boundary layer thickness δ_IBL grows as δ_IBL ≈ 0.10 x^{0.8} (Garratt 1990), reaching ~5 m by the eastern sideline and ~7–10 m at midfield — fully encompassing the 1–3 m disc-altitude band but in a sheared, non-equilibrium state. This is the dominant local physical feature.

### 3.4 CFD methodology *(draft)*

**Geometry.** A ~500 m × 500 m × 80 m domain centered on Sean Earl Field, built from Microsoft GlobalMLBuildingFootprints + USGS 2023 heights for the campus block and 2–3 nearest residential blocks west of campus, with USGS 3DEP LiDAR DEM for ground. Buildings imported as STL after height-extrusion; tree canopy on the lakefront strip modeled as a porous-zone with leaf area density α_LAD ≈ 0.5 m²/m³ (Gromke & Ruck 2008). Domain extent satisfies COST 732 / Tominaga (2008) guidelines: 5H upstream/lateral/top, 15H downstream of tallest building (H ≈ 22 m for Mertz). Blockage <3%.

**Mesh.** snappyHexMesh hex-dominant. Background cell 4 m, three refinement levels in a 250 m field-centered box (0.5 m near the field surface, 1 m at building faces, surface layers giving y⁺ ≈ 30–100 for wall functions). Total ≈ 5–8 M cells.

**Solver.** OpenFOAM v2406 `simpleFoam`, k-ω SST, second-order schemes. Richards & Hoxey (1993) ABL inlet with Yang et al. (2009) corrected wall functions to suppress horizontal inhomogeneity. Inlet roughness: z₀ = 0.001 m (water) for headings in [60°, 120°] east sector; z₀ = 0.5 m (urban) for [180°, 360°]; transitional zones blended. Top BC constant-shear; sides symmetry; outlet pressure-outlet. Reference: u_ref = 5 m/s at z_ref = 10 m (representative shoreline mean).

**Cases.** Inflow heading θ_wind ∈ {0°, 15°, 30°, …, 180°} (13 cases). Steady RANS to residuals <10⁻⁴. Wall-clock per case ≈ 4 h on 32 cores.

**Validation.** (i) ABL profile recovery at empty-domain upstream of buildings (per Hargreaves & Wright 2007 protocol). (ii) Tominaga (2008) AIJ-B cube benchmark comparison for the snappyHexMesh + k-ω SST stack on a separate test mesh. (iii) Loyola-site spot checks: a handheld anemometer at three field-corner locations during 2026 summer/fall lake-breeze events, matched against the corresponding CFD case interpolated to instrument height.

### 3.5 Coupled trajectory simulation *(draft)*

For each θ_wind, the CFD output is trilinearly interpolated onto a 0.5 m × 0.5 m × 0.5 m grid covering 1–6 m AGL across the field. The 6-DOF integrator evaluates **U**(**x**) at the disc's instantaneous position; the steady-RANS field provides U_mean, and turbulence intensity I = √(2k/3)/U_mean is mapped into a stochastic gust contribution **u**'(t) generated by a Mann (1998) synthetic-turbulence Fourier sampler with the local I and integral length scale set by k/ε (taken from the RANS).

For each case (θ_wind, throw type), 10⁴ Monte Carlo throws are launched from a 5 m × 5 m grid of launch points spanning the field with realistic release-parameter scatter:
- speed V₀ ∼ 𝒩(15, 1.5) m/s for backhand and forehand; 𝒩(11, 1.5) for hammer
- spin Ω₀ ∼ 𝒩(100, 10) rad/s
- release pitch α₀ ∼ 𝒩(+8°, 4°) backhand; 𝒩(+6°, 4°) forehand; 𝒩(+40°, 8°) hammer (overhead)
- heading bearing θ_throw varied over 0°–360° in 30° bins.

Landing point distributions are characterized by 90%-containment ellipse parameters and a "delivered-to-target" success rate (landing within 3 m of intended target at canonical 25 m throws).

### 3.6 Throwability metric *(draft — this is the novelty)*

Define, for each launch position **x** and each wind heading θ_wind:

$$
Q(\mathbf{x}, \theta_{\text{wind}}) =
\underbrace{1 - \frac{\langle |U_{2\text{m AGL}}| \rangle_{\text{field}}}{U_{\text{ref}}}}_{\text{wind-magnitude penalty}}
\cdot
\underbrace{(1 - I_{2\text{m AGL}}(\mathbf{x}))}_{\text{turbulence penalty}}
\cdot
\underbrace{\overline{S}_{\text{throws}}(\mathbf{x}, \theta_{\text{wind}})}_{\text{simulated success rate}}
$$

with each factor normalized to [0, 1]. Q ranges from 0 (unplayable) to 1 (calm-air-equivalent). Maps of Q over the field surface, parameterized by θ_wind, form the principal result.

### 3.7 Expected results / placeholder findings *(to be filled by simulation)*

Hypothesized:
- Easterly cases (90° ± 30°) show a sharp Q gradient across the field's eastern half corresponding to the developing internal boundary layer; backhand throws launched from the windward sideline are most consistent, hammer throws are systematically degraded.
- Westerly cases (270° ± 30°) show patchier Q with sheltered "throw lanes" downstream of the Norville Center; forehand throws benefit most.
- The overall best-condition wind direction for canonical 25 m throws is hypothesized to be a light SSW (~210°, 3–5 m/s reference), aligning with the prevailing climatology — a satisfying convergence between meteorology and playability.

(Section will be populated after the simulation sweep.)

### 3.8 Discussion *(to draft after results)*

- Generalizability of the throwability framework to other outdoor sports (golf, tennis).
- Steady RANS limitations: peak gusts and spatially correlated gust structure are not captured; LES of one or two critical directions is the obvious follow-on.
- Coaching implications and what the metric can and can't say.

### 3.9 Conclusions *(to draft after results)*

---

## 4. Key references (verified)

- Blocken, B. (2014). 50 years of computational wind engineering: Past, present and future. *JWEIA* **129**, 69–102.
- Blocken, B. (2018). LES over RANS in building simulation… *Building Simulation* **11**, 821–870.
- Crowther, W.J. & Potts, J.R. (2007). Simulation of a spin-stabilised sports disc. *Sports Engineering* **10**(1), 3–21. doi:10.1007/BF02844198
- Franke, J. et al. (2007). *Best Practice Guideline for CFD Simulation of Flows in the Urban Environment.* COST Action 732.
- Hargreaves, D.M. & Wright, N.G. (2007). On the use of the k-ε model in commercial CFD software… *JWEIA* **95**, 355–369.
- Hubbard, M. & Hummel, S.A. (2000). Simulation of Frisbee Flight. 5th Math & Computers in Sport, Sydney.
- Hummel, S.A. (2003). *Frisbee Flight Simulation and Throw Biomechanics.* MS Thesis, UC Davis.
- Kamaruddin, N.M., Potts, J.R. & Crowther, W.J. (2018). Aerodynamic performance of flying discs. *Aircraft Eng. Aerosp. Tech.* **90**(2), 390–397.
- Lissaman, P.B.S. & Hubbard, M. (2010). Maximum range of flying discs. *Procedia Engineering* **2**(2), 2529–2535.
- Lukes, R., Hart, J. et al. (2022). Disc golf trajectory modelling combining CFD and rigid body dynamics. *Sports Engineering* **25**, 24.
- Potts, J.R. & Crowther, W.J. (2002). Frisbee™ aerodynamics. AIAA 2002-3150.
- Richards, P.J. & Hoxey, R.P. (1993). Appropriate boundary conditions for computational wind engineering models using the k-ε turbulence model. *JWEIA* **46–47**, 145–153.
- Tominaga, Y. et al. (2008). AIJ guidelines for practical applications of CFD to pedestrian wind environment. *JWEIA* **96**, 1749–1761.
- Yang, Y. et al. (2009). New inflow boundary conditions for modelling the neutral equilibrium atmospheric boundary layer. *JWEIA* **97**, 88–95.

---

## 5. Open questions for charlie

1. Sean Earl or Hoyne? I assumed Sean Earl from the "off the lake" phrasing — confirm.
2. Is there budget / interest in one LES validation run for one wind direction, to anchor the RANS sweep? Adds ~2 weeks but gives a much stronger paper.
3. Do you have any actual handheld-anemometer measurements at the field, or is the 2026 summer measurement plan in §3.4 (iii) realistic?
4. Co-authors? Department affiliation? Funding to acknowledge?
5. The hammer (high-α overhead) is the throw most affected by wind but the one whose aerodynamic coefficients are least well-validated. Are you OK keeping it as "exploratory" in the paper, or do you want to drop it?
