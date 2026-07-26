# P5 — Throwability Metric and Synthesis

Compresses P4's high-dimensional landing-statistics array into the
single-scalar throwability metric `Q(x, θ_wind)` defined in
`Project1_Frisbee_LoyolaWind_draft.tex` §10, and produces the paper figures.

## Throwability metric

For each launch position `x` and wind heading `θ_wind`:

```
Q(x, θ_wind) = wind_factor(x, θ_wind) * turbulence_factor(x, θ_wind) * success_rate(x, θ_wind)
```

with each factor normalized to [0, 1]:

- `wind_factor = 1 - |U(x, 2 m AGL, θ_wind)| / U_ref`, clipped to [0, 1]
- `turbulence_factor = 1 - I_2m(x, θ_wind)`, clipped to [0, 1]
- `success_rate` = averaged delivery-to-target rate across all throw types and bearings

## Files

- `throwability.py` — computes `Q(x, θ_wind)` from `flow_fields.npz` + `landing_stats.npz`.
- `plot_throwability_map.py` — saves contour plots per heading, summary plot of best/worst.
- `paper_figures.py` — builds the four figures the manuscript needs.

## Workflow

```bash
python throwability.py \
    --fields ../P3_CFD/flow_fields.npz \
    --landings ../P4_Coupling/landing_stats.npz \
    --u-ref 5.0 \
    --out Q.npz

python plot_throwability_map.py Q.npz --outdir figs/
python paper_figures.py Q.npz --outdir figs/
```

## Paper figures (per draft §3)

- **Fig. 1** — Site map with annotated buildings + IBL schematic.
- **Fig. 2** — 6-DOF trim-glide trajectory examples (3 wind regimes).
- **Fig. 3** — CFD mean-flow snapshots at 2 m AGL for 3 representative headings.
- **Fig. 4** — Throwability heatmap Q(x, y) at the best heading vs. worst heading.

Stub `paper_figures.py` generates the right axes / labels / panel layouts.
