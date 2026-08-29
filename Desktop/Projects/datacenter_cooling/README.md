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
test_dc_thermal.py   analytic-limit + conservation validation (11/11 pass)
cfd/                 Tier 1: de Vahl Davis-validated 2D buoyant solver that measures
                     the recirculation curve φ(η) the model now uses by default
openfoam/            Tier 2 scaffold: a v1 toolchain-validation run + read_openfoam.py
requirements.txt
```

## Run
```
pip install -r requirements.txt --break-system-packages
python test_dc_thermal.py      # 11/11 checks
python optimize_cooling.py     # prints the optimum, writes the two figures
```

## Representative result
Generic 20-rack, 200 kW hall, containment η=0.30: baseline PUE 1.320 → **optimized 1.236
(26% cooling-energy cut)** by raising supply 18→21 °C until the hottest rack sits on the
27 °C limit. **Containment ROI, at fixed design airflow: PUE 1.294→1.281 as η 0→1
(ΔPUE 0.0126), with the supply temperature it unlocks rising 23.0→26.0 °C.**

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
> ⚠ **Do not read the *p* < 1 optimum as a recommendation yet** — see `DESIGN.md` §8. The
> airflow axis is still one-sided (over-provisioning carries no thermal benefit in this
> model) and server-fan power is not charged, which is the classic way a model games PUE.

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
