# datacenter_cooling — reduced-order data-hall cooling optimization

Tier-0 of an industry-facing data-center cooling optimization project. A fast,
open-source reduced-order model that finds the **PUE-optimal supply temperature and
airflow** for an air-cooled hot-aisle/cold-aisle hall, subject to ASHRAE rack-inlet
limits — the object you optimize before any CFD, and the first artifact to show a company.

See `DESIGN.md` for the full problem framing, the tier ladder (reduced-order → OpenFOAM
calibration → real company layout), the open-source stack, and the outreach/data roadmap.

## The one tension it captures
Raise supply-air temperature → chiller COP rises → less cooling energy. But hot-air
**recirculation** then lifts the hottest rack's inlet toward its ASHRAE limit, and fighting
it with more airflow costs fan power ∝ flow³. The optimizer walks straight to that edge.

## What's here
```
dc_thermal.py        reduced-order model: rack heat balance, recirculation influence
                     matrix B, chiller/fan/pump power, PUE; generic_hall() builder
optimize_cooling.py  SLSQP minimize PUE over (supply temp, provisioning) s.t. ASHRAE
                     inlet limit; writes fig1_inlet_map.png, fig2_containment_sweep.png
test_dc_thermal.py   analytic-limit + conservation validation (16/16 pass)
cfd/                 Tier 1: de Vahl Davis-validated 2D buoyant solver that measures
                     the recirculation curve φ(η) the model now uses by default
openfoam/            Tier 2 scaffold: a v1 toolchain-validation run + read_openfoam.py
requirements.txt
```

## Run
```
pip install -r requirements.txt --break-system-packages
python test_dc_thermal.py      # 16/16 checks
python optimize_cooling.py     # prints the optimum, writes the two figures
```

## Representative result
Generic 20-rack, 200 kW hall, containment η=0.30: baseline PUE 1.320 → **optimized 1.235
(26% cooling-energy cut)** by raising supply 18→21 °C until the hottest rack sits on the
27 °C limit. **Containment ROI, at fixed design airflow: PUE 1.279→1.277 as η 0→1
(ΔPUE 0.0022), with the supply temperature it unlocks rising 26.4→27.0 °C.**

> **That ROI was ΔPUE 0.0126 until 2026-08-31, when the CFD gained an airflow axis**
> (`PROVISIONING_2026-08-31.md`). The 2D solver had no through-flow, so it had no
> provisioning ratio `p` in it and the eta-only curve was being applied at design
> airflow. With `p` measured, the cold aisle at design airflow is already flooded and
> containment has little left to remove — and below design airflow the ROI is *exactly
> zero* at every anchor, because mass balance pins `B`. Two more things came out of it:
> **the free optimum did not move at all** (21.02 °C, p = 0.685, PUE 1.2351, identical
> to five digits under all three recirculation models and across a 16× sweep of the p
> anchor), and **the "~7% recirculation floor" is thermal diffusion, not recirculated
> air** — sweeping the tracer's diffusivity at η = 0.95 takes the signal from 0.0025 at
> Le = 50 to 0.06945 at Le = 1, where it reproduces the temperature channel to five
> decimals, while the same sweep at η = 0 moves it 0.4%. This cell runs at Ra = 4×10³;
> a hall is 10⁹–10¹².

> **Airflow branch corrected 2026-08-29: the optimizer had been buying fan savings that
> violate mass conservation.** `provisioning` is defined as supplied/demanded airflow, and
> the row sum of `B` as the recirculated fraction of inlet air — so if the CRAHs deliver
> only *p* of what the servers ingest, at least (1−*p*) of every inlet stream *must* be
> recirculated exhaust. The model charged only a mild `φ(η)·(1+(1−p))` bump: at *p* = 0.6 it
> returned a recirculated fraction of **0.23 where conservation demands 0.40**. Because fan
> power goes as *p*³, under-blowing was close to a free lunch and the optimum sat exactly on
> the `p_min` = 0.6 bound. With `φ ≥ 1−p` enforced the optimum becomes *interior* (*p* =
> 0.68) and the headline moves **1.215 → 1.236, a 33% → 26% cooling-energy cut**.
>
> The containment number moved for a second reason. Under the free (T_s, *p*) optimum the
> mass-balance floor exceeds the geometric term at *every* η, so η drops out of `B` and that
> sweep is **flat** — it measures the airflow axis, not containment. The ROI quoted above is
> therefore taken at fixed *p* = 1, which is what an operator actually buys containment for.
> Its ΔPUE 0.0126 supersedes the 0.022 quoted below. Both curves are drawn in
> `fig2_containment_sweep.png`; the flat one is shown, not hidden.
>
> ⚠ **The *p* < 1 optimum still leans on a one-sided airflow axis** — see `DESIGN.md` §8.
> Over-provisioning carries no thermal benefit in this model, so the optimizer can only ever
> want less airflow.

> **Server-fan power is now charged, and PUE is no longer the only objective (2026-08-29).**
> Cutting CRAH airflow raises inlet temperatures, server fans ramp, and that power lands
> *inside* PUE's IT denominator — so a model that omits it can lower the ratio while raising
> the hall's bill. `dc_thermal` now adds a ramp of 0.4 % of rack IT power per K above a
> 25 °C knee (published measurements give 0.35–0.5 %/°C in the 23–30 °C band and 7–20 %
> total over 15→35 °C), charges the chiller for it, and exposes `total_kw()` alongside
> `pue()`. `optimize(..., objective="total")` ranks by energy instead of by ratio.
>
> **How much does the gaming matter here? Extra hall power from optimizing PUE rather than
> energy:**
>
> | ramp %/K | T_allow 27 | 32 | 35 | 40 |
> |---|---|---|---|---|
> | 0.35 | 0.00 kW | 0.00 | 0.07 | 0.45 |
> | **0.40** | **0.00 kW** | 0.15 | 0.54 | 1.34 |
> | 0.50 | 0.25 kW | 1.48 | 2.56 | 4.19 |
> | 0.75 | 1.44 kW | 5.57 | 8.34 | 12.03 |
> | 1.00 | 2.62 kW | 9.64 | 14.12 | 19.87 |
>
> At the mid-band coefficient and ASHRAE's 27 °C *recommended* limit the two objectives
> agree: **the inlet constraint, not the objective, is what closes the hole.** It opens on
> the *allowable* envelopes (A1 32 °C, A2 35 °C, A3 40 °C) and at the steep end of the
> published ramp, where the lowest-PUE design in the whole study — **PUE 1.197** at
> 1.0 %/K and T_allow 35 — gets there by moving **20 kW, 10 % of nominal IT**, into the
> denominator, and burns 14 kW more than the energy optimum.
>
> The practical reading: **PUE is safe as an objective for this hall only while the 27 °C
> recommended limit is held.** Relax to an allowable envelope and you must rank on kW.
> Which regime a real hall is in depends on its ramp coefficient, so that is the constant
> to pin first.

> **Numbers restated 2026-08-09 under the CFD-calibrated recirculation curve.** The
> paragraph above previously read *optimized 1.22, 31% cut, supply 18→22 °C, sweep
> 1.24→1.20*. Those came from the `"linear"` placeholder φ = 0.28(1−η), which
> `dc_thermal.py` stopped defaulting to on 2026-07-25 (`recirc_model="cfd"`) — two days
> *after* those figures and that paragraph were written. The placeholder drove recirculation
> to zero at η=1, so it over-credited containment: it showed the sweep buying ΔPUE = 0.039,
> where the calibrated curve buys **0.022 — 42% less**. The superseded figures are kept in
> `superseded_linear_2026-07-23/`. Regenerate with `python optimize_cooling.py`.

## Status & honesty
Tier 0 + Tier 1. The recirculation matrix `B` is no longer a pure placeholder: its
containment dependence φ(η) = 0.007 + 0.228·e^(−η/0.79) is measured by the validated 2D
solver in `cfd/`. What remains uncalibrated is the 3D rack-to-rack coupling (the
off-diagonal `B_ij`), so absolute hotspot temperatures are still indicative. The
`openfoam/` v1 run validates the toolchain only and yields no physical number — see
`DESIGN.md` §9 for exactly what a usable v2 needs. Full limitations in `DESIGN.md` §8.
