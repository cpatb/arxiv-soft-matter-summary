# Audit — Hydrogel camera-analysis scripts

**Date:** 2026-06-30 · **Scope:** `align.py`, `analyze.py`, `annotate.py`, `visualize.py`
**Reference program:** `Documents/hydrogel_drying_roadmap.tex` (crack/peel/curl mode-selection, camera-first)

---

## 1. What's there

Four Python scripts forming a **coherent, self-consistent pipeline** for **wrinkle / line-feature morphology** on TIFF image stacks:

| Script | Role | Key deps |
|---|---|---|
| `annotate.py` | napari GUI to hand-paint line features on each source `*.tif`; auto-saves `<stem>_roi.tif` (labeled) and `<stem>_skeleton.tif` on navigate/close. | napari, qtpy, tifffile, skimage, scipy |
| `align.py` | napari GUI to shift each timepoint's `*_roi.tif` layer into registration; writes `offsets.json` (dx, dy per file). | napari, tifffile, numpy |
| `analyze.py` | Batch: walk `*_roi.tif`, skeletonize, label 8-connected segments, measure per-segment length/orientation/centroid, count fork (branch) points; write `results.csv`. | numpy, pandas, tifffile, skimage, scipy |
| `visualize.py` | Read `results.csv` → 7 figure types (density+forks vs time, length histograms/violin, orientation roses, density-vs-forks, top/bottom counts, superimposed aligned ROIs). | pandas, numpy, matplotlib, tifffile |

**Intended workflow:** `annotate` → (`align`) → `analyze` → `visualize`. The data contract is consistent end-to-end (filename `<n>min`, folder layout `<root>/<condition>/<sample_N>/<face>/`, `offsets.json`, `results.csv` columns). This part is well done.

## 2. State / completeness

- **Syntax:** all four parse cleanly (the initial `annotate.py` failure was only a cp1252 read-encoding artifact from the `←`/`→` glyphs; it is valid UTF-8 Python).
- **Dependencies:** all resolve in the active environment: cv2 4.13 (unused by these scripts), numpy 2.4.2, scipy 1.17.1, skimage 0.26, pandas 3.0.1, tifffile 2026.3.3, napari 0.6.6, matplotlib 3.10.8, qtpy present.
- **Runnable:** core numeric paths exercised on synthetic data without error (`analyze_mask`, `count_forks`, `summarize_image`; `visualize.per_image` groupby.apply verified OK under pandas 3.0). The two napari scripts are GUI/interactive and were not launched (no footage present), but nothing blocks them.
- **Docstrings:** present in all four with usage lines. Bare `sys.argv` parsing (no argparse/`--help`), but a usage string prints on missing args.

## 3. Issues found (severity-ranked)

### S1 — CRITICAL (architectural, not a code bug): the pipeline measures the wrong physics for the stated program
These scripts are a **wrinkle / surface-pattern** toolkit (skeleton segments, orientation roses, branch-point "forks", areal density, top-vs-bottom counts). That is the vocabulary of the **swelling / surface-instability** work, **not** the camera-first **drying-fracture** program the roadmap describes. **None** of the roadmap's headline observables are implemented:
- mode classification (crack / peel / curl)
- onset time `t_c`
- peel thickness `h` (side-view / kymograph)
- crack spacing / polygon area `Λ, A_p`
- curl radius `ρ` (arc fit)
- crack velocity `v_c` (high-speed)
- mass-loss `m(t)`, side-view `H(t)`

There is also **no video ingestion** (works on TIFF stacks, not high-speed frames), **no side/kymograph view**, and **no px/mm calibration-from-fiducial step** even though the roadmap explicitly requires "fix px/mm per run." As-is, this pipeline does not serve the drying crack/peel/curl deliverable. It is reusable *scaffolding* (napari annotate/align + batch-measure + plot), but the measurement layer must be rewritten to the roadmap observables. **This needs author decision, not a mechanical fix.**

### S2 — MINOR: `count_forks` over-counts junctions
`skeletonize` frequently leaves a *cluster* of pixels with ≥3 neighbors at a single branch, so a single junction can report as several forks (a synthetic "+" reported 5). This inflates `n_forks` (used in `plot_timeseries`, `plot_density_vs_forks`). Fixable by connected-component-labeling the branch-point mask and counting components, but it changes scientific output; flagged, not silently changed.

### S3 — MINOR: fragile, undocumented folder/filename contract
`analyze.py` infers `condition/sample/face` purely from path depth (`<root>/<condition>/<sample_N>/<face>/`) and timepoint from a `(\d+)min` regex. A non-conforming layout silently yields `sample=None` or a mislabeled condition rather than an error. "top/bottom" faces are assumed by `plot_topbottom_counts`. Document or validate the expected tree.

### S4 — LOW: connectivity mismatch between annotate and analyze
`annotate._save` relabels with default 4-connectivity; `analyze` relabels the skeleton with 8-connectivity (`_EIGHT_CONN`). Harmless today (analyze only consumes `roi>0` as a mask, then re-skeletonizes/re-labels), but inconsistent and a latent trap if the stored `_roi` label IDs are ever trusted.

### S5 — LOW: minor robustness edges
- `summarize_image`/`per_image` `std` is `NaN` for single-segment images (errorbar tolerates it).
- `annotate._load_img` returns an RGB source array as-is; `viewer.add_image` is called without `rgb=True`, so a true-color source could be shown as a 3-slice stack. Only matters if inputs are RGB.

## 4. Positive findings (explicitly checked)

- **No hardcoded px/mm scale.** `px_per_um` is an *optional* CLI arg (`analyze.py [px_per_um]`), defaulting to `None` → lengths reported in px only. Compliant with the "calibrate per video, never hardcode" rule.
- **No hardcoded Shamah-derived numbers** anywhere (grep for `0.5`, `13.7`, `14×`, peel-plateau, velocity-ratio, swelling table → none). The quarantine is respected.
- **No hardcoded absolute paths.** All paths derive from the CLI `<root_folder>` argument.

## 5. Fixes applied

**None.** The scripts are syntactically valid, imports resolve, and they contain no hardcoded scale/Shamah numbers or broken imports, so there was no *safe, non-restructuring* fix to make. The one true defect (S2, fork over-counting) alters scientific output and the dominant problem (S1) is architectural; both are flagged for author judgment rather than changed unilaterally, per the conservative brief.

## 5b. Update 2026-07-06 — S1 addressed via new `drying_analysis/` package (option a)

Recommendation 1(a) taken: the legacy wrinkle scripts are left intact, and a **separate
`drying_analysis/`** package now implements the roadmap observables the wrinkle toolkit
lacked. It provides `calibration` (px/mm from a fiducial, Recommendation 2), `swelling_kinetics`
(D, k_m, Bi), `indentation_modulus` (E), `curl_radius` (ρ), `crack_metrics` (A_p, Λ),
`toughness` (Γ, l_fc, ℓ=Γ/E), and `kymograph` (side-view crack-vs-peel ordering,
Recommendation 3, the `t_c`/`ρ`/`h`-adjacent side path). **36/36 synthetic-ground-truth
checks pass** (`python -m drying_analysis.validate`), with convention cross-checks and a
hand-computed reference so a wrong constant can't pass silently. A science-code-review pass
found and fixed 5 silent-correctness bugs (Robin thickness convention, kymograph tip-tracking
+ classifier that were cancelling, t=0 normalization, pure-shear vs uniaxial f). Two MOC prose
discrepancies are flagged in `drying_analysis/NOTES_moc_flags.md`. **Still needs the user:** the
segmentation/tracing that produces crack masks / traced arcs / frame stacks from raw footage
(hand-trace via the existing `annotate.py` or threshold), and the real physical inputs. S2/S3/S4
in the legacy scripts remain open (unchanged).

## 6. Recommendations (author judgment required)

1. **Decide the scripts' role.** Either (a) keep them as the legacy wrinkle toolkit and build a *separate* drying-fracture pipeline for the roadmap observables, or (b) refactor the measurement layer onto crack/peel/curl metrics. Do not present the current outputs as drying-fracture results.
2. **Add per-run px/mm calibration** (fiducial/ruler in frame → scale stored per video), as the roadmap mandates, before any length is reported in physical units.
3. **Add video ingestion + a side/kymograph path** (`t_c`, `h`, `v_c`, `ρ`): the current top-view-only, static-TIFF design cannot produce the headline phase diagram.
4. If keeping `count_forks`, fix the junction over-count (label the branch-point mask, count components).
5. Validate/document the `condition/sample/face/…min` folder contract, or replace with an explicit manifest.
