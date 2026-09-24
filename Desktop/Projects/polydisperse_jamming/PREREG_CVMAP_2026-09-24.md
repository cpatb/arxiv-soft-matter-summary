# PREREG — `cvmap`: the founding observation through the drop solver (CV at fixed weighed φ)

Written 2026-09-24, before any cell (commit time is the record). Inputs and blind numbers:
`predict_cvmap.py` → `out/predict_cvmap.txt`, `out/cvmap_cells.txt` (committed with this file; must
not be regenerated after data). Launcher `scripts/cvmap_2026-09-24.sh` (reuses the G0-gated jamwin
binary `db413e87dd21`, as kmermap). Solver `wc_drop_impact.c`, cs_fit base, L9, WC_CLOCK 1, θ₀ 5°,
TEND 6, MP_EXPO 2, default EPS_MIN 0.01, DNS PHI0 0.62 (so PHI = x·0.62, PHIM = χ·0.62, as kmermap).

## Why
The project exists because of one bench observation (author, 2026-08-31): **at similar particle size,
raising CV moves an impacting drop from purely dilatant to spreading and inertial-dominated.** The
DEM + closure say CV acts through x = φ/φ₀ — it raises φ₀ (0.6318 → 0.6441 → 0.6943) while χ is nearly
flat (0.9106 / 0.9121 / 0.9170, μ_m 0.3) — so at a FIXED WEIGHED φ raising CV walks the suspension down
and out of the window. That was COMPOSED in `jam_map.py` §4 and never run in the drop. The doublet axis
was (`kmermap`); the founding axis was not.

## Closed form (no DNS; `out/predict_cvmap.txt` §1)
A CV-0 suspension that flows at rest has x₀ < 1; the CV-0.5 suspension at the same weighed φ has
x = x₀·φ₀(0)/φ₀(0.5) = 0.91004·x₀, and can jam only if that is ≥ χ(0.5) = 0.91702. At μ_m 0.3 the
point estimate says NO weighed φ exists at which it can (margin +0.0070), **but the margin is 1.2×
its propagated PHI_SCALE uncertainty — NOT resolved** (φ₀(0.5) has n = 1 and a capped cell reads LOW,
which would widen the margin; that is a direction, not a measurement). What IS resolved: CV 0 → 0.5 moves
x by ≈ one full window log-width, and at μ_m 3 a band x₀ ∈ [0.969, 1) can still jam (friction-scoped).
The grid below sits at x₀ ≤ 0.95, where every CV-0.5 cell is outside its window by ≥ 0.05 in x.

## Inputs (all DEM, μ_p 0.3, derjaguin; from graded artifacts, parsed by the script)
| CV | χ | δ ratio | σ\* ratio (d10) | source |
|---|---|---|---|---|
| 0 | 0.91063 | 1 | 1 | `three_knobs` p̂ 3e-5 |
| 0.2 | 0.91215 | 0.938 | 0.957 | `actshape_cv02.txt` |
| 0.5 | 0.91702 | 0.815 | 0.735 (d32: 1.145) | `grade_actshape_fs.txt` F2/F3, `grade_actshape.txt` H2 |
DNS δ = 1.3169 × ratio, σ\* = 300 Pa × ratio (kmermap convention). ⚠ χ(CV 0) here is 0.9106
(the isotropic-ladder records), not kmer's 0.9083 (the primary-built singlets): same rung, different
campaigns; this campaign uses ONE source for all three CVs.

## Grid (27 cells)
x₀ ∈ {0.920, 0.932, 0.950} (x at CV 0) × CV ∈ {0, 0.2, 0.5} × U₀ ∈ {2, 4} m/s = 18 cells, plus a channel
split at U₀ 4 for every x₀ (9 cells):
- **W** (window only): x, χ from CV 0.5; δ, σ\* from CV 0.
- **A** (activation only): x, χ from CV 0; δ, σ\* (d10) from CV 0.5.
- **D** (full CV 0.5 with σ\* at fixed d32, 343.5 Pa): the size-convention bracket.
JAM_FHAT = f_c inside the window, 1.0 outside (J ≡ 0 there by construction; never reported as a jam).
In-window cells: CV 0 at all three x₀ (f_c 0.895 / 0.761 / 0.560); CV 0.2 at x₀ 0.932 (f_c 0.977) and
0.95 (0.776). All CV 0.5 cells are outside.

## Predictions
- **C1 (the switch is by construction; recorded, not graded as evidence).** Every CV-0.5 cell has J ≡ 0.
- **C2 (spreading — the author's observation, the substantive DNS gate).** At every (x₀, U₀) where the
  CV 0 and CV 0.5 cells are both healthy: max Rg_r(CV 0.5) > max Rg_r(CV 0) (strict, continuous
  channel), and β_max(CV 0.5) ≥ β_max(CV 0) − 1/64. Where CV 0.2 is also healthy, max Rg_r is
  non-decreasing across 0 → 0.2 → 0.5. Informative only if ≥ 2 (x₀, U₀) pairs are decided.
- **C3 (channel).** At each x₀ where c0.0_U4, c0.5_U4, W, A are all healthy, with
  Δ = ln Rg(c0.5) − ln Rg(c0.0):
  (a) window share ln Rg(W)/Rg(c0.0) ÷ Δ ∈ [0.8, 1.2];
  (b) |ln Rg(A)/Rg(c0.0)| ≤ 0.2 |Δ|;
  (c) convention: |ln Rg(D)/Rg(c0.5)| ≤ 0.2 |Δ|.
  PASS if every decided x₀ passes, and ≥ 1 x₀ is decided.
- **C4 (apex recoil out of sample, on a NEW knob; S2\* = 0.1, S2 as `grade_jamsig.sig`).** Every
  healthy J+ cell (J > 0.01) has S2 < 0.1; every healthy J0 cell S2 ≥ 0.1. Here σ\*, δ and χ change
  with the class, so the test also asks whether recoil tracks the jam state across a composition
  change. Informative only if ≥ 1 healthy cell of each class. The CV 0.2 x₀ 0.920 cell (outside,
  Θ 83.6, close to the kmermap doublet x 0.85 cell that died) is the hardest J0 case.
- **C5 (health).** rc 0 and t/τ ≥ 5.9, else EXCLUDED and printed (never a jam, never a FAIL of
  C2–C4). Expected deaths: CV 0 at x₀ 0.95 (f_c 0.56; kmermap died at f_c 0.58, U₀ 1/2).

## Scope
Stress-controlled activation (jamwin P2), L9, default cap (J magnitudes are cap numbers; classes and
orderings are what `epssig` showed robust), DEM → DNS transfer of χ, σ\*, δ, the DEM φ axis is not
the bench's (only x transfers).

## Amendment 1 — pre-data referee (registered 06:28 CDT 2026-09-24, before any cell landed; commit time is the record)
A physics-referee pass on this file and `grade_cvmap.py` returned no FATAL and 3 MATERIAL, all acted on
here. No bar is loosened; each change separates a question the original wording had merged.
1. **C3(c) is graded separately.** The registered "decided" list for C3 is (c0.0, c0.5, W, A) at U₀ 4,
   but the original grader counted a DEAD D arm as a FAIL of C3. Now: C3 = (a) AND (b) over decided x₀;
   **C3c** (d32 convention, bar unchanged ≤ 0.2) is graded over the x₀ where D is also healthy, and a
   dead D is "undecided", never a fail.
2. **A floor on Δ.** C3's shares divide by Δ = ln Rg(c0.5)/Rg(c0.0). An x₀ with |Δ| < 0.05 is
   UNDECIDED for C3 (printed with its numbers, not graded). 0.05 is set before data from the only
   prior scale on record: kmerdecomp's activation-only arms moved ln Rg by ≤ 0.03 (≤ 3 %), so a Δ
   below 0.05 cannot separate a window share of 1 from the activation channel's own size.
3. **C2 split into C2a and C2b.** C2a (the author's observation, the headline): CV 0 → 0.5 endpoints,
   Rg strictly larger and β not smaller by more than 1/64, at every decided (x₀, U₀), ≥ 2 decided.
   C2b (unchanged strict 0 → 0.2 → 0.5 monotonicity) is graded as registered and reported beside it.
   Reason: at x₀ 0.920 the CV-0.2 cell sits 0.0097 below its χ with Θ 83.6, near the closure's
   divergence, so C2b can fail on a near-singular thickening that is not evidence against the window
   mechanism — and C2b failing there would be a finding about thickening, not a reason to drop it.
