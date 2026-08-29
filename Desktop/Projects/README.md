# Desktop/Projects — working-area index

> ## Current project state lives in the vault, not here
>
> **`C:\Users\compl\Documents\Obsidian Vault\Status\_Router.md`** has one card per active
> project: where it stands, the next action, the load-bearing numbers, the landmines.
> Start there. This file only says *what is in the folder and where*, never how a project
> is going.

Everything here is a working/scratch area. Nothing in it is a source of truth about
project status.

## Layout

Every project owns a directory. As of 2026-08-29 the top level holds no loose files
except this one: reports moved to `_reports/`, and the two manuscripts that used to sit
at the root (`chicago_run_bike_monitor`, `friday_experiments`) now have directories of
their own.

### `_reports/` — cross-project reports, audits and sweeps

Dated, point-in-time documents that span more than one project, so they never belonged
inside any single project directory. Overnight session reports, readiness sweeps, edit
passes, the limbo audit, the solver landscape, week plans. `_reports/README.md` indexes
them. They are historical records: read them for how something came to be, not for
where it stands now.

### Papers and manuscripts

| directory | what it is |
|---|---|
| `bidet/` | Minimum-water bidet nozzle optimization, CFD-anchored. Canonical: `bidet_optimization_1.JOURNAL.tex`. |
| `frisbee/` | Disc flight through urban CFD (Sean Earl Field). Pipeline in `SeanEarl_DiscFlight/P1..P5`. |
| `slime mold/` | Physarum-grown CTA 'L' expansion, plus `bike/` for bike-lane gap ranking. |
| `paper1/` | Drop impact, dense TPM vs cornstarch. See `PDF_PROVENANCE.md` for which PDF is current. |
| `friday_experiments/` | Five free-surface oobleck experiments. |
| `chicago_run_bike_monitor/` | DIY wearable air-quality sensing on Chicago running/cycling routes. |
| `LDV Paper/` | Laser Doppler Velocimetry reference library and slides. No manuscript source. |

### Code and models

| directory | what it is |
|---|---|
| `dst_fit/` | Constitutive identification for dense suspensions, with an identifiability verdict. |
| `datacenter_cooling/` | Reduced-order data-hall cooling optimizer, plus 2D CFD calibration. |
| `axidrop/` | Benchmarks and tolerances over the Basilisk suite. |
| `tier0_wyart_cates_closure/` | Tier-0 Wyart-Cates closure. |
| `impact_stress/` | Wagner / shear-stress toolkit for drop impact. |
| `hydrogels/` | Drying-fracture program: napari wrinkle pipeline plus the roadmap LaTeX. |
| `silver_to_color/` | Chunked video colourization pipeline. |
| `catfish school/` | Collective-motion numerics. |
| `uroflow_splashback/` | Uroflow splashback model. |
| `bloodstain_bpa/` | Bloodstain pattern analysis estimators. |
| `beach_sand_dilatancy/` | Dilatancy experiment design. |
| `rheo_products/` | Product spec for the computational-rheology line. |

### Earlier-stage and reference

`Chicago transit board/`, `Corn Starch Bulk Rheology AM/`, `ShearJamming_Theory/`,
`bio_fluids_bridge/`, `carbopol/`, `edm/`, `glymphatic_transport/`,
`thrombosis_jamming/`, and `ARCHIVE/` (old C. elegans / SCM work, correctly archived).

## Conventions

- **Build artifacts are not kept.** `.aux`, `.log`, `.out`, `.fls`, `.fdb_latexmk`,
  `.synctex.gz`, `.toc`, `.blg` and `__pycache__` are ignored by the repository's
  `.gitignore` and were cleared on 2026-08-29. Regenerate with `latexmk -pdf`.
- **Figures travel with a sibling `.md`.** A figure's caption, statistics and caveats go
  in a markdown file of the same name beside it, not printed on the axes.
- **A number in a manuscript should be generated, not typed.** Several projects here
  emit `*_generated.tex` macro files from their analysis code for exactly this reason,
  with a `\MISSINGNUM` fallback so an absent file renders visibly rather than silently.
