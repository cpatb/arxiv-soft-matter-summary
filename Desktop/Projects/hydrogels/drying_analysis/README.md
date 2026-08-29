# drying_analysis — camera-first measurement layer for PEG drying-fracture

Created 2026-07-06. This package implements the **crack / peel / curl** roadmap
observables from `Documents/Claude/Projects/Hydrogels/hydrogel_MOC.md`, closing the
`AUDIT_camera.md` **S1** gap: the existing `align/analyze/annotate/visualize.py`
scripts measure surface-**wrinkle** morphology (skeleton segments, orientation
roses), which is the *swelling* program's vocabulary. **None** of the
drying-fracture observables were implemented. This package adds them.

## What each module measures (→ MOC observable)

| Module | Measures | MOC ref | Method |
|---|---|---|---|
| `calibration` | px/mm from a fiducial | §5B "fix px/mm per run" | ratio of pixel-dist to known length |
| `swelling_kinetics` | D, k_m, **Bi** from m(t) | §5A | Crank early-√t + long-time exp; Robin-BC (D,k_m) fit |
| `indentation_modulus` | **E** | §5 | Hertz + Dimitriadis finite-thickness correction |
| `curl_radius` | curl radius **ρ** | §6.1 | algebraic (Kåsa) circle fit to side-view arc |
| `crack_metrics` | polygon area **A_p**, spacing **Λ** | §6.1 | label complement of crack mask → regionprops |
| `toughness` | **Γ** (pure shear), **l_fc**, **ℓ=Γ/E** | §5 | ∫σ dλ; notch-series knee fit |
| `kymograph` | crack-vs-peel **ordering** | §5B | depth-time image + tip-trajectory classifier |

## Status & scope

- **All 36 validation checks pass** against synthetic ground truth
  (`python -m drying_analysis.validate`). Every module recovers a known input
  (Crank series → D; forward Hertz → E; noisy arc → ρ; grid → A_p; neo-Hookean →
  Γ; synthetic descent → crack→peel). The suite includes convention cross-checks
  (Robin sim → Crank limit), a hand-computed indentation reference, and
  discriminating cases (crack-only vs crack→peel, surface-at-top vs -bottom) so a
  wrong constant or a cancelling pair of bugs cannot pass silently. Two MOC prose
  discrepancies surfaced during the build are flagged in `NOTES_moc_flags.md`
  (not silently "fixed").
- These are the **pure-numeric measurement functions**: they run on real data
  the moment it exists. What they need from you (the physical-data frontier):
  - **Segmentation/tracing upstream:** `crack_metrics` takes a binary crack mask
    and `curl_radius`/`kymograph` take traced points / a frame stack. Producing
    those from raw footage (threshold, edge-detect, or hand-trace via the existing
    napari `annotate.py`) is the one interactive step not automated here.
  - **Real inputs:** gravimetric m(t) curves (balance), indentation δ + sphere
    spec, notch-series λ_c, pure-shear σ(λ) ladder, high-speed side-view stacks.
- **Not included** (need decisions/data, per MOC): the full mode *classifier*
  (crack/peel/curl label per run), which is a downstream rule over these outputs
  once real morphologies are in hand; and onset time t_c, which is just the first
  frame the relevant module fires (trivial once footage exists).

## Conventions (per lab CLAUDE.md)

- Pixel coordinates are `(row, col) = (y, x)` throughout (NumPy/skimage order).
- No hardcoded scale: every physical output goes through `px_per_mm` from
  `calibration`.
- SI units internally for physics (metres, Pa, J/m²); calibration converts px↔mm.
- Determinism: `validate.py` seeds `np.random.default_rng(20260706)`.

## Run

```
cd Desktop/Projects/hydrogels
python -m drying_analysis.validate      # 27/27 expected
```

## Key references (methods used)
- Crank, *The Mathematics of Diffusion* (1975); Tanaka & Fillmore, JCP 70, 1214 (1979).
- Dimitriadis et al., Biophys. J. 82, 2798 (2002) — finite-thickness indentation.
- Rivlin & Thomas, J. Polym. Sci. 10, 291 (1953) — pure-shear toughness.
- Chen, Wang & Suo, EML 10, 50 (2017) — fractocohesive length.
- Ma, Lowensohn & Burton, PRE 99, 012802 (2019) — polygon crack scaling.
- Chung, Regev & Mahadevan, Soft Matter 12, 7855 (2016) — exfoliation ordering.
