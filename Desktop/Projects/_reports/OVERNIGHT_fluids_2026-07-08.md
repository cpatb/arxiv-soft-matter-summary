# Overnight fluids/soft-matter exploration — 2026-07-08

You asked me to develop low-hanging-fruit analytical + simulation projects in fluids /
soft matter, in the spirit of the bidet, frisbee, and bodily-fluids work, either
tangential to current projects or bodily-fluid dynamics. I built **two complete,
validated, self-contained projects** and scoped a ranked backlog of more.

Everything is runnable Python (numpy/scipy/matplotlib), with tests, validation against
known analytic results, and figures. Written to `Desktop/Projects/`.

---

## Project 1 — Bloodstain Pattern Analysis: ballistic region-of-origin
`Desktop/Projects/bloodstain_bpa/`  ·  bodily fluid **and** a drop-impact bridge

**Idea.** A bloodstain is the inverse of your drop-impact problem: reconstruct the impact
(and its 3-D origin) from the splat. Three ellipse measurements decode it: impact angle
`α=arcsin(W/L)`, direction from the major axis, and **impact speed from stain size via a
drop-impact spread factor `ξ=D_stain/D₀=f(We_n)`** (Pasandideh-Fard 1996). That last one
is the bridge to your Basilisk drop-impact sims. With the full impact velocity vector I
integrate the flight **backward** (gravity + Schiller–Naumann drag) and find where all
paths converge.

**What I did.** Forward model (emit from a known origin → synthetic stains) + two
reconstructions: the classical straight-line "stringing" method and a ballistic
back-projection. Full round-trip validation, a bias sweep, the drop-impact bridge figure,
and a measurement-noise robustness study. 6/6 correctness tests pass (free-fall vs
analytic, forward/backward reversibility, spread-model invertibility, origin recovery).

**What I found.**
- With matched physics the ballistic method recovers the known origin to ~0 % (a
  consistency check). The method is self-consistent and correctly coded.
- **The standard forensic straight-line method overestimates the height of origin by
  ~1.3–2.2×** for impact-spatter-onto-floor geometry, and the bias *grows with impact
  speed* (flatter trajectories → more gravitational curvature). The ballistic method
  removes it (mean bias +0.7 % vs +157 %; ~230× reduction).
- Under realistic noise (5 % ellipse, 3° azimuth, 10 % drop size) the ballistic estimate
  keeps a small, honest **−18 ± 5.5 cm** bias driven almost entirely by drop-size
  uncertainty, vs the straight-line's **+212 cm**.

**Where to go next.** Wall+floor combined patterns; a Monte-Carlo *posterior* over the
region of origin rather than a point; **replace the analytic spread law with one fit to
your actual drop-impact simulation data** (blood-specific `ξ(We_n, surface)`); validate
against a public BPA benchmark set. This is genuinely paper-shaped: "gravity/drag bias in
region-of-origin estimation, and a drop-impact-calibrated correction."

---

## Project 2 — Uroflow splashback: Rayleigh breakup + stream tracking
`Desktop/Projects/uroflow_splashback/`  ·  bodily fluid **and** a bidet sibling

**Idea.** The bidet paper aims a jet at a surface to *maximise* cleansing shear; this aims
a jet at a surface to *minimise* splashback: same machinery, opposite objective. The
controlling physics is **Plateau–Rayleigh breakup**: a coherent stream splashes far less
than a drop train (Splash Lab, 2013), so you want to impact *before* breakup length `L_b`.

**What I did.** Computed the fastest-growing varicose mode from Weber's viscous dispersion
relation (validated to 2.5 % against Rayleigh's inviscid `ω_max`), got `L_b`, tracked the
stream as a projectile to the target, and built a **splash-index map over (standing
distance, aim angle)** plus a sit-vs-stand comparison and a ranked recommendation search.

**What I found (all physics-derived, matching the Splash Lab's experimental advice).**
1. **Stand close** — impact before Rayleigh breakup (coherent ≪ drops).
2. **Aim shallow** — grazing impact kills the normal momentum that ejects secondaries.
3. **Sit / lower the source** — drier at every distance.
Driest vs. wettest scanned config span ~30× in splash index.

**Where to go next.** Add a standing-water target (pool-depth vs splash, the real
urinal-design lever); calibrate `ln(r/δ₀)` and the splash index against a phone high-speed
clip; reuse the bidet paper's Watson wall-jet solution for the post-impact spreading film.

---

## Idea backlog (ranked, not yet built)

Low-hanging, self-contained, analytical+sim, bodily-fluid or tangential:

1. **Mucus/saliva filament thinning (CaBER-in-silico).** Capillary-thinning of a
   viscoelastic biopolymer thread (Oldroyd-B / FENE-P 1-D thinning ODE). Extract relaxation
   time from the exponential-thinning regime. Bodily fluid; classic, well-validated against
   the Entov–Hinch analytic law. ~1 evening. **Highest value/effort.**
2. **Cough/sneeze droplet dispersal (Wells–Bourouiba).** Ballistic + evaporation-settling
   of a respiratory droplet cloud; the evaporation-vs-fall competition and the turbulent
   puff range. Reuses Project-1's drag/flight code almost verbatim. Bodily fluid, high
   relevance.
3. **Tear-film breakup time.** Thin-film dewetting / thinning ODE (evaporation + Marangoni)
   predicting non-invasive break-up time. Bodily fluid, soft-matter thin films.
4. **Frisbee flight-envelope solver (extends your frisbee work).** 6-DOF disc with the
   standard `C_L(α), C_M(α)` aero + gyroscopic precession; produce the roll-vs-throw-angle
   stability map. Directly extends existing work.
5. **Sneeze-vs-bidet aerosolization of a toilet plume.** Combines Project-1 flight + the
   bidet jet; the "toilet plume" bioaerosol problem: droplet size → suspension time.
6. **Ketchup/blood yield-stress start-up flow (Bingham/HB).** 1-D transient
   pipe/slot flow of a Herschel–Bulkley fluid; ties to your Carbopol yield-stress bath.

---

## Notes on rigor
Per your CLAUDE.md: library-first (scipy `solve_ivp`, `brentq`, `minimize_scalar`),
seeded RNG (`default_rng`), validation against analytic references (free-fall, inviscid
Rayleigh), boundary/edge tests, and every claim backed by a figure or a printed number.
Both projects report their own limitations honestly (geometry-dependence of the BPA bias
magnitude; ideal-theory over-prediction of `L_b`).

*Built autonomously overnight. Nothing here has been run against real data yet. The
validation is internal (round-trip / analytic-limit) and the next step for both is a
single real measurement to calibrate.*
