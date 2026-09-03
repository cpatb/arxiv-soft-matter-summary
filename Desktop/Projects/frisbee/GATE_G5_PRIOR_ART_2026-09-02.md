# GATE G5 — the prior-art screen (WP4a)

**Run 2026-09-02.** Brief: `CLAUDE_CODE_FRISBEE_HARDENING.md` §4. This is one of the three
gates that can end the programme in its present form. Zero compute; every reference below
was resolved to a DOI or a patent number in this session.

**Verdict in one line: the programme survives, but two of its three stated novelty claims
do not, and the manuscript's own "novelty hook" is prior art from 2012 and 2016.**

The work does not change. The framing does, and it has to change before anything is
posted, because both dead claims are the kind a referee finds in ten minutes.

---

## 0. What was screened

Three claims are made in the tree. They are screened separately because they fail
separately.

| tag | claim | where it is written | verdict |
|---|---|---|---|
| **N1** | "Coupling CFD wind field into disc trajectory … not done in published literature. *This is the novelty hook.*" | `Project1_Frisbee_LoyolaWind_PAPER.REVIEWED.tex:158` | 🔴 **DEAD.** Done in published literature in 2012 and patented for a stadium in 2016 |
| **N2** | "Nobody has written down the equivalent [pedestrian-comfort] criterion for a *projectile*." | `FRISBEE_GENERALIZATION.md` §1 | 🔴 **TOO STRONG.** The urban-air-mobility literature writes CFD-based turbulence-threshold criteria for flying bodies. Survives only in a narrowed form |
| **N3** | The closed-form landing-dispersion law `σ_landing = κ g σ_u T_f²/(2V_d) · √Φ(N)`, `N = L_f/L` | `P6_Theory/dispersion_law.py`, brief §C5 | 🟡 **Novelty survives, as an application not a derivation** — `Φ(N)` is a standard random-vibration response integral, proved below rather than asserted. ⚠ **But the law itself is already rejected by the project's own 6-DOF model** (`THEORY_WIND_SPEED_2026-08-30.md`), so this cannot be the paper's positive claim. See §3a |

---

## 1. N1 is dead, twice

### 1a. Kakimpa, Hargreaves & Owen (2012) — coupled CFD–RBD, in a wind-engineering journal

- *An investigation of plate-type windborne debris flight using coupled CFD–RBD models.
  Part I: Model development and validation.* **J. Wind Eng. Ind. Aerodyn. 111, 95–103
  (2012).** `10.1016/j.jweia.2012.07.008`
- *…Part II: Free and constrained flight.* **JWEIA 111, 104–116 (2012).**
  `10.1016/j.jweia.2012.07.011`

A CFD wind field coupled to six-degree-of-freedom rigid-body dynamics for a body flying
through it. That is the coupling this paper calls its hook, published fourteen years ago,
for a non-lifting plate rather than a lifting disc.

The same construction is standard in exterior ballistics under the identical name — e.g.
Sahu, *Advanced coupled CFD/RBD calculations of free-flight projectile aerodynamics*,
DOD-UGC'05, 86–91, `10.1109/dodugc.2005.7`; Costello, Gatto & Sahu, AIAA 2007-6582,
`10.2514/6.2007-6582`. "CFD/RBD" is a term of art there.

### 1b. US 11,052,285 B2 — the same thing, in a stadium, for a ball

*System and method for using weather applied metrics for predicting the flight of a ball.*
Weather Applied Metrics Inc. (formerly Martin Meter LLC); Martin, Farley, Farley, Jackson,
Sinton & Ozen. **Priority 2016-10-29.** Continuations: US 11,361,132; US 12,485,311.

Claim language, verbatim from the specification:

> "wind measurements by the wind sensors are used as inputs into a CFD model that produces
> 3D wind vectors at grid-points above the field … These 3D wind vectors are used as
> inputs into the trajectory model"

with the vectors linearly interpolated to the ball's position at each timestep, and ANSYS
/ OpenFOAM / Star-CCM+ / COMSOL named as the solvers. This is the pipeline of `P3_CFD` →
`P4_Coupling`, in a sports venue, filed a decade ago.

**What it does not do, and this is the whole of what survives N1:** it is *deterministic*.
There is no gust model, no turbulence variance, no Monte Carlo, no landing-point
distribution — the wind field is a fixed input and the output is one predicted trajectory.
The same is true of the closest sports-engineering analogue, Giljarhus, Gooding &
Njærheim, *Disc golf trajectory modelling combining computational fluid dynamics and rigid
body dynamics*, **Sports Engineering 25, 26 (2022)**, `10.1007/s12283-022-00390-5`, which
uses CFD only to obtain aerodynamic coefficients and then flies the disc in a **uniform**
freestream.

**Required edit.** `PAPER.REVIEWED.tex:158` must stop saying the coupling is novel. The
defensible sentence is that the coupling is established (cite Kakimpa Part I/II and the
patent family), and that what is new here is carrying a **spatially resolved turbulent
field** through it to a **landing-point distribution** rather than a single trajectory.
That is a narrower claim and it is one the tree can actually defend.

## 2. N2 is too strong

The generalization memo's line — that the COST 732 / Lawson / NEN 8100 pedestrian-comfort
tradition has no counterpart for a projectile — was true of the *pedestrian* literature and
is not true of the flying-vehicle literature, which has moved fast since 2020:

- *Urban morphology and wind effects on UAV flight risks: a CFD-based risk mapping
  framework.* Sustainable Cities and Society (2026), ScienceDirect S2210670726002830.
- *Urban wind field effects on the flight dynamics of fixed-wing drones.* **Drones 9(5),
  362 (2025)**, `10.3390/drones9050362`.
- *Ten questions concerning urban wind environments for the safe utilization of urban air
  mobility.* Building and Environment (2026), ScienceDirect S0360132325016026.
- Watkins, Thompson, Loxton & Abdulrahim, *On low-altitude flight through the atmospheric
  boundary layer.* **Int. J. Micro Air Vehicles 2(2), 55–68 (2010)**,
  `10.1260/1756-8293.2.2.55`.

These build CFD-derived, site-resolved **threshold** criteria (wind speed, turbulence
intensity, shear) for flying bodies. So "nobody has written the criterion for a flying
object" is refutable from a single search.

**What survives, and it is a real distinction.** Every one of those criteria is a
*threshold on a controlled vehicle* — a go/no-go on whether a drone can hold station or a
pilot can hold a path. This project's object is **unguided after release**: there is no
control loop, so the quantity of interest is not an exceedance threshold but the
**variance of the terminal position**, and the correct group is `σ_u/V_d` (the disc's own
airspeed) rather than the turbulence intensity `σ_u/U` that a station-keeping criterion
uses. The memo already argues exactly this in §2 and in `dispersion_law.py`'s docstring.

**Required edit.** Rewrite `FRISBEE_GENERALIZATION.md` §1's last sentence to claim the
*unguided-projectile* case against the *controlled-vehicle* literature, and cite the four
references above as the thing being distinguished from rather than the gap being filled.

## 3. N3 survives — and here is the articulation §4 demanded

The brief sets the bar: *"If you cannot articulate what `Φ(N)` adds beyond a von Kármán
response integral, there is no theory paper."*

`Φ(N)` **is** a response integral. That was established here rather than assumed, by
reconstructing it from the physics and never from the closed form:

A high-Stokes body (the brief measures `St ≈ 8.5`, ballistic length 51 m against `L` ≈ 3–6 m)
does not follow eddies; it is *forced* by them, so the lateral displacement is a double
time integral of the gust acceleration,

```
D = ∫₀^T ∫₀^t a(s) ds dt = ∫₀^T (T − s) a(s) ds
```

With `a` proportional to a stationary field of exponential (Dryden / Ornstein–Uhlenbeck)
correlation `e^{−|Δs|/L}`, substituting `s = uT` and `N = T/T_L = L_f/L` gives the
variance-reduction factor

```
Var(D) / Var(D)|_{N→0}  =  4 ∫₀¹∫₀¹ (1−u)(1−v) e^{−N|u−v|} du dv
```

Evaluated numerically against `regrade_gate_g1.phi` at nine values of `N` over eight
decades:

```
N=1e-06   Phi_code=0.9999997333   4*I(N)=0.9999997333   ratio=1.00000000
N=0.05    Phi_code=0.9868043737   4*I(N)=0.9868043744   ratio=1.00000000
N=0.25    Phi_code=0.9366620039   4*I(N)=0.9366620076   ratio=1.00000000
N=1       Phi_code=0.7805956079   4*I(N)=0.7805956222   ratio=0.99999998
N=2       Phi_code=0.6303304085   4*I(N)=0.6303304372   ratio=0.99999995
N=4       Phi_code=0.4450548481   4*I(N)=0.4450549055   ratio=0.99999987
N=4.4     Phi_code=0.4193781383   4*I(N)=0.4193782014   ratio=0.99999985
N=10      Phi_code=0.2274662671   4*I(N)=0.2274664102   ratio=0.99999937
N=100     Phi_code=0.02626674667  4*I(N)=0.02626674668  ratio=1.00000000
```

(`N` = 4.4 is the realised value at Sean Earl. Agreement to 8 significant figures; the
`1e-7` drift at `N` = 10 is `dblquad`'s roundoff on the kink at `u = v`, not a discrepancy.)

**So the identification is exact.** `Φ(N)` is the continuous-turbulence response integral
`σ² = ∫ |H(iω)|² Φ_gg(ω) dω` of Etkin (*Turbulent wind and its effect on flight*, **J.
Aircraft 18(5), 327–345 (1981)**, `10.2514/3.57498`) and Hoblit (*Gust Loads on Aircraft*,
AIAA 1988), evaluated for the one transfer function aerospace never needs — a pure double
integrator, i.e. a body with **no** restoring dynamics — and for the Dryden rather than the
von Kármán spectrum, which is what makes it close in elementary functions at all.

**Three things it adds, and they should be the paper's claim:**

1. **It is closed.** The aerospace tradition evaluates the integral numerically per
   aircraft, because `H` is an aeroelastic transfer function. For an unguided projectile
   `H` is `1/(iω)²` and the integral collapses to
   `Φ(N) = 8[1/(3N) − 1/(2N²) − e^{−N}/N³ + (1 − e^{−N})/N⁴]`, a design formula needing
   two site numbers (`σ_u`, `L`) and three flight numbers (`T_f`, `V_d`, `L_f`).
2. **The controlling group is spatial, not spectral.** Aerospace parameterises by reduced
   frequency `ωc/V`, because the aircraft has an internal timescale. This body has none, so
   the only group left is `N = L_f/L`, *flight length in integral scales*. That is the
   right variable precisely because `V_d/U ≈ 8–10` — the projectile, not the wind, supplies
   the advection, so Taylor's hypothesis applies with the *body's* speed. This is the part
   that generalises across objects and it is the memo's transferable content.
3. **`κ` is not in the aerospace integral.** `κ = √(4 + (C_Lα/C_L)²(σ_w/σ_u)²)` folds the
   *lifting* response to vertical gusts into the same scalar. Etkin's sweeping correction
   handles the kinematics of a body traversing a frozen field; it does not supply this.

**Required edit.** §567 ("Throwability metric (the novelty)") must cite Etkin (1981) and
Hoblit, say in the introduction that the response-integral method is standard in aerospace
and absent from the sports and urban-projectile literature, and confine the novelty to
items 1–3 above. Written that way the claim is small, true, and unattackable.

## 3a. …but N3 cannot be the paper's positive claim, because the project already rejected it

This is the part the brief's §4 did not anticipate, and it changes what the paper is more
than either dead novelty claim does. `THEORY_WIND_SPEED_2026-08-30.md` §"Varying the flight
rejects the ballistic law" tested the law's own shape against the project's 6-DOF model:

| law | fitted `b` | 95 % CI | verdict |
|---|---|---|---|
| ballistic `√Φ(N)` (this law) | **0.671** | [0.522, 0.818] | **shape rejected** (predicts ≈ 1.70) |
| tracer `1/√N` | 0.950 | [0.640, 1.218] | not rejected (predicts ≈ 0.50) |

replicating in 6 of 6 blocks, with the ballistic slope below 1.0 in **12 of 12** sub-blocks
(sign test `p` = 0.0005). The same document is careful that this does not crown the tracer
law either — the tracer slope flips sides between grouping levels (0.955 pooled, 1.306
within heading), and in **11 of 12** fine blocks the measurement is *bracketed*,
ballistic < 1 < tracer.

So the honest reading is: **neither a single time integral nor a double one describes this
disc, and the closed form whose novelty §3 just defended is the rejected endpoint of that
bracket.** A finite-Stokes body with aerodynamic restoring should indeed sit between them,
which is a result — but it is a *negative* one, and it is the paper's theory contribution.

There is no conflict with GATE G1. G1 tests `Φ(N)` as a **ratio at a fixed flight**
(predicted 1.90 against measured 1.868/1.735), where `T_f`, `V_d` and `L_f` do not move and
the partial-averaging factor is the only thing being asked about. The rejection comes from
letting the *flight* vary, which is what exposes the `T_f²/V_d` prefactor. Both statements
are true of different tests and the manuscript must not quote G1's pass as support for the
law as a whole.

**Consequence for the framing.** The paper's positive content is the empirical urban result
— buildings halve the mean wind and raise `σ_u` by half again in two of three directions,
turbulence intensity 12 % → 25–38 % — plus the bracketing negative above. The closed form
belongs in the theory section as *the endpoint that fails*, with its provenance (an Etkin /
Hoblit response integral for a pure double integrator) stated so a referee can see it was
not invented for the occasion.

## 4. What this screen found that the brief did not ask for

**A free, measured validation target that the model has never been run against.** Sasakawa
& Sakurai, *Biomechanical analysis of the sidearm throwing motion for distance of a flying
disc: a comparison of skilled and unskilled Ultimate players*, **Sports Biomechanics 7(3),
311–321 (2008)**, PMID 18972880. Seventeen players, three-dimensional high-speed video:

| | skilled (n = 10) | unskilled (n = 7) |
|---|---|---|
| release speed | 21.7 ± 1.7 m/s | 20.7 m/s (**not significantly different**) |
| spin | 12.9 ± 1.3 rps = **81.1 ± 8.2 rad/s** | 9.4 ± 1.3 rps = **59.1 ± 8.2 rad/s** |
| distance | **51.4 ± 6.6 m** | **29.5 ± 7.6 m** |

Hold `V ≈ 21 m/s` and sweep `ω` from 59 to 81 rad/s: the 6-DOF model must move the range
from ~30 m to ~51 m. Nothing else in the project tests the spin sensitivity the
gyroscopic-stability argument rests on, and this costs no compute.

⚠ **One correction to the brief.** §4 item 7 describes this as a human-throw test and
implies a backhand. It is a **sidearm** (forehand) study. The project's 6-DOF release
conditions are backhand; Winograd & Engsberg (2012) report a precession difference between
the two. So the test is still worth running, but the arm must be matched or the mismatch
stated — do not quote it as a backhand anchor.

## 5. References resolved in this screen

| # | reference | identifier | status |
|---|---|---|---|
| 1 | Kakimpa, Hargreaves & Owen 2012, CFD–RBD Part I, JWEIA 111:95–103 | `10.1016/j.jweia.2012.07.008` | ✅ resolved |
| 2 | Kakimpa, Hargreaves & Owen 2012, CFD–RBD Part II, JWEIA 111:104–116 | `10.1016/j.jweia.2012.07.011` | ✅ resolved |
| 3 | Richards 2012, *Dispersion of windborne debris*, JWEIA 104–106:594–602 | `10.1016/j.jweia.2012.02.026` | ✅ resolved (18 citations); **abstract not retrieved — see §6** |
| 4 | US 11,052,285 B2 (+ 11,361,132; 12,485,311), Weather Applied Metrics Inc, priority 2016-10-29 | patent | ✅ specification read |
| 5 | Giljarhus, Gooding & Njærheim 2022, Sports Engineering 25:26 | `10.1007/s12283-022-00390-5` | ✅ resolved |
| 6 | Crowther & Potts 2007, *Simulation of a spin-stabilised sports disc*, Sports Engineering 10:3–21 | `10.1007/BF02844198` | ✅ resolved |
| 7 | Sasakawa & Sakurai 2008, Sports Biomechanics 7(3):311–321 | PMID 18972880 | ✅ resolved, numbers extracted |
| 8 | Etkin 1981, *Turbulent wind and its effect on flight*, J. Aircraft 18(5):327–345 | `10.2514/3.57498` | ✅ resolved |
| 9 | Watkins et al. 2010, Int. J. Micro Air Vehicles 2(2):55–68 | `10.1260/1756-8293.2.2.55` | ✅ resolved; full text paywalled |
| 10 | *Urban wind field effects on the flight dynamics of fixed-wing drones*, Drones 9(5):362 (2025) | `10.3390/drones9050362` | ✅ resolved |
| 11 | Sahu 2005, coupled CFD/RBD free-flight projectile | `10.1109/dodugc.2005.7` | ✅ resolved |
| 12 | Costello, Gatto & Sahu 2007, AIAA 2007-6582 | `10.2514/6.2007-6582` | ✅ resolved |

Curcin & Harris (1986) — the brief's item 1, and the one it said to close first — **was not
found.** Four searches returned Monte-Carlo synthetic-turbulence flight-simulation work of
that era (Reeves et al., NASA TP-2469, 1985; the AIAA rational-approximation von Kármán
papers) but nothing matching that author pair and year. **Treat the citation as
`[UNVERIFIED]` and do not put it in the `.bib`** until someone produces it; it may be a
mis-transcription. It does not change the verdict: N1 is already dead on items 1, 2 and 4.

## 6. What this screen did NOT settle

- **Richards (2012), *Dispersion of windborne debris*.** Title-level, this is the closest
  published thing to N3 — a *dispersion* law for a body flying through turbulent wind.
  Crossref and Semantic Scholar both return metadata with the abstract elided, and
  ScienceDirect returns 403 from this machine. **This is the single remaining item that
  could still cost N3 its "closed form" claim,** so it must be read before submission.
  Get it through the university proxy; it is 9 pages.
- Potts (2005) PhD thesis and Potts & Crowther (2002) AIAA 2002-3150 remain unacquired, so
  the model still runs on Hummel's transcription of the coefficients rather than the source.
  That is a reproducibility exposure, not a novelty one, and it belongs to WP4b.

## 7. Gate ruling

**GATE G5: PASS, conditional on the four edits in §1, §2, §3 and §3a being made before
posting.** The programme is not ended. Its novelty is smaller than claimed and now has a
defensible perimeter:

> A spatially resolved urban CFD field is carried through a 6-DOF lifting-disc model to a
> **landing-point distribution** rather than a single trajectory; the site is shown to
> remove mean wind while adding fluctuation, so that shelter is not protection; and the
> dispersion is shown to be **bracketed** between the tracer and ballistic limits, so that
> neither a single nor a double time integral describes a finite-Stokes lifting body.

Neither the coupling nor the response-integral method may be claimed as new, and the closed
form must be presented as the rejected endpoint of the bracket rather than as the result.
