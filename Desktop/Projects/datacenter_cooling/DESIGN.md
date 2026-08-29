# Data-Center Cooling Optimization — Project Design & Roadmap

**Author:** Charlie Burton · **Started:** 2026-07-23 · **Status:** Tier-0 model (11/11) + **Tier-1 CFD calibration done** (2D buoyant solver validated vs de Vahl Davis; recirculation φ(η) measured and fed back into the model).
**Intent:** industry-facing engineering demonstrator, open-source simulation throughout, generic layout first → real company plans later.

> **Tier-1 result (2026-07-23):** a de Vahl Davis-validated 2D buoyant-convection solver (`cfd/`) *measured*
> the recirculation fraction vs containment: **φ(η) = 0.007 + 0.228·e^(−η/0.79)**. Two findings the Tier-0
> placeholder missed — a **~7% recirculation floor at full containment** (the placeholder wrongly drove it to
> zero) and **concave, diminishing returns** (containment ROI flattens past η≈0.5). This is now the model's
> default (`recirc_model="cfd"`). See `cfd/README.md`.

---

## 1. The problem, and why it's a good portfolio piece

Cooling is 30–40% of a data center's energy. The whole air-cooling design collapses to **one tension**:

> Raise the supply-air temperature and the chiller gets dramatically more efficient (higher COP) — but hot-air **recirculation** then pushes the hottest rack's inlet past its ASHRAE limit and you cook a server.

Every real operator lives on that knife-edge, plus a second one (blow more cold air to fight recirculation, but fan power scales as flow³). That makes it a clean, honest **constrained optimization** problem with a real objective (PUE / cooling energy), real decision variables (supply temperature, airflow provisioning, containment, later rack placement), and a hard physical constraint (ASHRAE inlet limits). It's the kind of problem a company can immediately recognize as *their* problem — which is exactly what makes it a strong demonstrator and a credible reason to reach out to them for real layouts.

It also rhymes with your existing toolkit: this is a reduced-order closure with an influence matrix and a scipy optimizer, the same shape as your `tier0_wyart_cates_closure` philosophy — cheap sanity model first, CFD only once the layout is chosen.

## 2. Tiered scope (the ladder — we lock each rung before climbing)

**Tier 0 — reduced-order, generic layout.  ← BUILT, this folder.**
Lumped per-rack heat balance + a geometry-driven recirculation influence matrix `B` + a chiller/fan/pump power model, optimized over (supply temp, provisioning) subject to the ASHRAE inlet constraint. Runs in milliseconds, so we can sweep containment, density, and the thermal envelope. This is the object you optimize and the artifact you show first.

**Tier 1 — calibrate the recirculation matrix with CFD.  ← DONE (2D), see `cfd/`.**
The one physically-soft piece in Tier 0 was `B` (how much of each rack's exhaust reaches each inlet). Rather than a placeholder, a **de Vahl Davis-validated 2D buoyant-convection solver** (`cfd/buoyant2d.py`, Ra=1e3→0.3%, 1e4→1.8%) now *measures* the recirculation fraction vs containment on a hot-aisle/cold-aisle cross-section (`cfd/aisle.py`), giving φ(η)=0.007+0.228·e^(−η/0.79). That curve is fed back into the reduced-order model, which stays the optimization engine (the classic capture-index / influence-coefficient reduction; VanGilder & Shrivastava). **The 3D rack-to-rack coupling (off-diagonal `B_ij`) and a hall-Reynolds turbulence model are the remaining upgrade — the `openfoam/` scaffold (`buoyantSimpleFoam`, kOmegaSST) is set up for it.**

**Tier 2 — real layout from a company.**
Swap the generic hall for an actual room: rack coordinates, per-rack loads, CRAH/CRAC positions and capacities, containment type. The Tier-0 model ingests it directly (it's just a load vector + geometry); Tier-1 CFD gives that specific room its `B`. Now the optimization produces recommendations someone can act on — and comparisons ("what if you added containment on rows 3–5", "what if supply went 18→24 °C").

**Tier 3 — richer objectives & architectures (optional, later).**
Multi-objective (energy vs. capex vs. resilience) with a real Pareto front (`pymoo`); workload-placement as a decision variable; and the air-vs-liquid crossover (rear-door / direct-to-chip / immersion) as rack density climbs past what air can hold. Only worth doing once a real case motivates it.

## 3. The Tier-0 model (what's in this folder)

`dc_thermal.py` — the model:
- Rack air-side rise `d_i = Q_i/(ρ c_p q_i)` (server ΔT).
- Recirculation influence matrix `B` (row sum = recirculated fraction of inlet air), scaled by containment `η` and worsened under airflow under-provisioning. Inlet temps solve `(I−B)(T_in−T_s) = B d`. Perfect containment `η=1 ⇒ B=0 ⇒` every inlet = supply.
- Chiller `COP(T_s)` rising with supply temp; fan power ∝ provisioning³ (affinity law); PUE = 1 + cooling/IT.

`optimize_cooling.py` — SLSQP minimization of PUE over (supply temp, provisioning) subject to hottest-inlet ≤ ASHRAE limit and provisioning ≥ floor. Produces the inlet-temperature map and the containment sweep.

`test_dc_thermal.py` — analytic-limit + conservation validation (perfect-containment limit, energy-balance closure, COP/fan monotonicity, cubic fan scaling, optimizer feasibility, containment monotonicity, CFD containment floor). **11/11 pass.**

**Representative Tier-0 result** (generic 20-rack, 200 kW hall, moderate containment η=0.3): a naive "cold and over-blown" baseline of PUE 1.32 optimizes to **PUE 1.21 — a 33% cooling-energy cut** — by raising supply 18→23 °C until the hottest rack sits exactly on the 27 °C limit and trimming airflow. Across the containment sweep, PUE 1.227→1.205 and the affordable supply temperature climbs 20.6→25.6 °C. Those are realistic numbers, which is the point of a demonstrator.

> **Restated 2026-08-09 under the CFD default.** This paragraph previously read *PUE 1.22, 31% cut, supply 18→22 °C, sweep 1.24→1.20, supply climbing 18.6→27 °C*, and `fig1`/`fig2` (2026-07-23) were drawn from it. All of that was computed with `recirc_model="linear"`, the placeholder this document's own Tier-1 box says over-credited containment — the switch to `recirc_model="cfd"` landed 2026-07-25, two days later, and the figures were never redrawn. Under the calibrated curve the containment sweep buys **ΔPUE 0.022 instead of 0.039 (42% less)**, and the supply-temperature range both starts higher (20.6 vs 18.6 °C, because φ(0)=0.235 < 0.28) and ends lower (25.6 vs 27.0 °C, because φ(1)=0.071 ≠ 0). The η=0.3 optimum improves slightly (1.215 vs 1.220) for the same reason the η=0 end does. Figures regenerated 2026-08-09; the superseded pair is in `superseded_linear_2026-07-23/`.

## 4. Open-source simulation & optimization stack

| Layer | Tool | Role |
|---|---|---|
| Reduced-order thermal/energy | **NumPy/SciPy** (this folder) | the optimization engine; always the fast inner loop |
| Room CFD | **OpenFOAM** (`buoyantSimpleFoam`) | calibrate the recirculation matrix `B`; validate hotspots for a specific layout |
| Facility energy / PUE | **EnergyPlus** or **OpenModelica** | whole-plant energy, chiller/economizer curves, weather-driven free-cooling hours |
| Optimization | **SciPy** (single-objective) → **pymoo** (multi-objective/genetic) | setpoints now; Pareto fronts and layout search later |
| Geometry / meshing | **Gmsh**, **ParaView** | CFD pre/post |

Everything above is open-source. Commercial DC-CFD tools (6SigmaDCX, Future Facilities) are the industry default but closed; the open stack is a deliberate, defensible choice for a portfolio piece and avoids licensing when you share it.

## 5. Validation strategy (keep it honest)

- **Analytic limits** (built): perfect containment, zero recirculation, energy balance closes to machine precision.
- **Literature bands:** PUE, COP(T), and ASHRAE inlet limits kept in published ranges — flagged in code comments as "verify before formal use."
- **CFD cross-check (Tier 1):** the reduced-order hotspot ranking must match OpenFOAM's on the same geometry, or `B` is wrong.
- **Real-data back-test (Tier 2):** if a company shares even one operating point (measured supply temp, IT load, cooling power), the model must reproduce its PUE within a stated tolerance before any recommendation is trusted.

## 6. Company outreach & data roadmap

The model is built so that **each piece of real data upgrades a specific placeholder**, which is exactly what you want to be able to say when you ask for it:

| What to request | What it unlocks | Sensitivity |
|---|---|---|
| Rack coordinates + per-rack kW | replaces the generic load vector → a real room | high — this is the whole geometry |
| CRAH/CRAC positions, capacities, control | real supply/return + provisioning limits | high |
| Containment type (aisle/chimney/none) | pins `η` instead of sweeping it | medium |
| One measured operating point (T_s, IT kW, cooling kW) | back-test / calibrate COP & fan curves | high — turns a demo into a validated model |
| Chiller & economizer performance curves | free-cooling hours, weather coupling (Tier 1–2) | medium |

**Outreach framing (for when you're ready):** lead with the runnable demonstrator and a concrete offer — "I've built an open-source reduced-order + CFD optimization model that finds the PUE-optimal supply temperature and containment strategy under ASHRAE constraints; I'd like to calibrate it against a real layout. Even a single anonymized floor plan and one operating point would let me back-test it." Ask for the *minimum* that unlocks the most (rack loads + one operating point), offer an NDA, and offer to share the calibrated model back. Colo providers, hyperscaler sustainability teams, and CRAH vendors (Vertiv, STULZ) are the natural first contacts; university/research data centers are often the easiest first real layout.

## 7. Portfolio framing & immediate next steps

This sits in your **applied-optimization / fluids-adjacent** track (nearer the bidet wall-jet optimization than the soft-matter closures), and it's the most obviously industry-legible thing you'd have. Progress and next steps:
1. ~~**Tier 0** reduced-order model~~ — **DONE** (11/11 tests).
2. ~~**Tier 1 CFD calibration**, the credibility step~~ — **DONE (2D):** validated buoyant solver → measured φ(η) → fed back into the model, with a decision-relevant finding (containment ROI flattens; ~7% floor). This is the real technical content and it's now real.
3. **Next — 3D OpenFOAM** (`openfoam/` scaffold): rack-to-rack coupling + turbulence, to populate the full `B_ij` and confirm the 2D trend. Then **outreach**, using the validated demonstrator + the CFD calibration as the concrete artifact.
4. Optional Tier-0 polish: a second headline figure (PUE vs rack density / ASHRAE A1–A4 envelope).

## 8. Honest limitations of Tier 0 (state these plainly)

- `B`'s containment dependence is now **CFD-calibrated in 2D** (Tier 1 done); the remaining gap is the 3D rack-to-rack off-diagonal coupling (the `openfoam/` scaffold). Absolute hotspot temperatures remain indicative until the 3D run, but the containment ROI curve is now grounded.
- Single lumped supply plenum; no per-CRAH spatial supply, no transient/ride-through, no humidity/economizer, no raised-floor pressure field. All are deliberate Tier-0 omissions with a known home on the ladder.
- COP and fan curves are generic linear/cubic forms; real chiller curves replace them at Tier 2.
- **The airflow-provisioning axis is one-sided and its optimum is provisional.** Two gaps,
  both on the `p` leg only (the supply-temperature leg is unaffected):
  1. *Over-provisioning carries no thermal benefit.* `under = max(0, 1-p)` means `B` is
     identical at `p = 1` and `p = 1.4`, so blowing extra cold air costs `p^3` fan power and
     buys nothing. Physically it floods the cold aisle and suppresses recirculation. A
     two-sided `phi(p)` is the fix and needs the 3D run to calibrate.
  2. *Server-fan power is not charged.* Cutting CRAH airflow raises rack inlet temperatures,
     server fans ramp, and that power lands **inside** the IT denominator — so PUE improves
     while total hall energy worsens. This is the standard way a cooling model games PUE, and
     nothing in Tier 0 prevents it.
  Together these mean the optimizer can only ever want *less* airflow. The mass-balance floor
  added 2026-08-29 (`phi >= 1-p`) stops the worst of it and moves the optimum off the `p_min`
  bound to an interior `p = 0.68`, but until (1) and (2) are closed, **report the containment
  ROI at fixed `p = 1`** (`containment_roi()`), not the free `(T_s, p)` optimum. Under the
  free optimum the mass-balance floor binds at every `eta`, so the containment sweep is flat
  and measures the airflow axis rather than containment.
