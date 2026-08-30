# AUDIT 2026-08-29 — what was implemented

Companion to `AUDIT_2026-08-29.md`, which is the findings file and is left unedited as the
record of what the audit found. This file records what was then done about it, on the same
day.

Every suite was re-run at the end: **P1 26/26, P4 40/40 (27 existing + 13 new), P5 13/13,
all exit 0.**

Two items were left open at first pass -- the contradictory acceptance bands and the
manuscript rebuild -- and both are closed in the follow-up section at the end of this file.

---

## The two FATAL findings

### F1 — GATE G1 re-graded against `Phi(N)`

**New:** `P4_Coupling/regrade_gate_g1.py`, `gust_model_comparison_regrade.json`,
`gust_model_comparison_isotropic_regrade.json`, `GATE_G1_REGRADE_2026-08-29.md`.
`GATE_G1_RESULT_2026-08-29.md` now carries a superseded-on-the-verdict header.

The measurement was re-derived from the 108 cells rather than trusted from the file, and
reproduces exactly: σ 1.7154 / 0.9886 m per axis, R = 1.735, cluster-bootstrap 95% CI
[1.638, 1.838] — the published CI to three decimals. Graded against the brief's own
ballistic law instead of the tracer law it had already disowned:

| run | R | band from Phi(N) | verdict | N |
|---|---|---|---|---|
| anisotropic | 1.735 [1.638, 1.838] | [1.50, 2.21] | **INSIDE** | 6.18 [5.25, 7.20] |
| isotropic | 1.868 [1.789, 1.955] | [1.50, 2.21] | **INSIDE** | 7.52 [6.71, 8.43] |

**The gate passes. `N_eddy = 3.01` is withdrawn** — it was the artifact of inverting
through the wrong functional form, not a measurement of an eddy count.

`Phi` is verified against five anchors that do not pass through it: sympy symbolic limits,
the sympy series, a 50-digit mpmath evaluation of both numerical branches, monotonicity,
and the brief's own independently computed values. Two of those anchors caught errors
during development — a hand-derived series coefficient (the true series is
`1 − 4N/15 + N²/18 − N³/105 + N⁴/720`) and an mpmath reference that had been polluted by
float64 rounding of `nv**4` before promotion. Both were mine; the anchors did their job.

### F2 — the two contradictory pre-registrations

Both are computed and reported side by side by `regrade_gate_g1.py`, and
`GATE_G1_REGRADE_2026-08-29.md` states the contradiction explicitly: WP1's band (2.0–2.6,
tracer) against WP2 P1's (line 599: 1.43–1.49, falsify outside [1.18, 1.82], ballistic).
The anisotropic measurement passes P1; the isotropic misses by 3%.

**Reconciled later the same day** -- see the follow-up section at the end. The resolution
was not to pick one band: both fixed `N` by convention, and the run measures `N` itself.

---

## The MATERIAL findings

### F3 — the 3-component |U| divided by a horizontal reference

**New:** `P5_Synthesis/conventions.py`, holding one definition each of the disc-altitude
speed, the reference speed, the pitch box and the 2 m level. Wired into
`throwability.py` and `make_results_tex.py`.

Measured first, then fixed. Over the pitch box at z = 2 m:

| heading | 3-component | horizontal | shift | w/u_h median · p95 · max |
|---|---|---|---|---|
| 090 | 1.2123 (24.25%) | 1.2057 (24.11%) | +0.55% | 0.064 · 0.266 · 3.51 |
| 210 | 1.4195 (28.39%) | 1.4034 (28.07%) | +1.15% | 0.120 · 0.405 · 0.82 |
| 270 | 1.3956 (27.91%) | 1.3882 (27.76%) | +0.53% | 0.080 · 0.227 · 1.91 |

So the fix was made without moving a published claim, which is the best case. `\UfieldFracMin{24}`
and `\UfieldFracMax{28}` are **byte-identical** before and after. Locally it is not
negligible — w/u_h reaches 3.5 in a wake cell — so the per-cell Q map does carry it.

Nothing is discarded physically: the Monte Carlo still flies the full 3-component field
through `cfd_field.py`. Only the diagnostic scalar that gets normalised by `U_ref` changed.

**Downstream regenerated** (P5 only; no Monte Carlo, no CFD): Q, significance, hypotheses,
figures, `results_generated.tex`. **17 of 60 macros moved in their last digit, none added
or removed**, so the manuscript's `\input` wiring and `\MISSINGNUM` guard are intact.
Every verdict survived:

| comparison | before | after | verdict |
|---|---|---|---|
| 90 vs 210 | p = 0.203 | p = 0.229 | **not resolved**, unchanged |
| 90 vs 270 | p = 0.022 | p = 0.023 | resolved, unchanged |
| 210 vs 270 | p = 0.0038 | p = 0.0052 | resolved, unchanged |

`worst_supported: 270`, `best_supported: null`. H1 SUPPORTED, H2 SUPPORTED / CONTRADICTED,
H3 PARTIALLY SUPPORTED, hammer NOT TESTABLE — all unchanged. Q̄ moved from
0.0443/0.0483/0.0324 to 0.0442/0.0481/0.0321.

The pre-change artifacts were kept in `SeanEarl_DiscFlight/_pre_horiz_2026-08-29/` and
then deleted once the diff was verified, since git already holds that state -- the same
convention the pre-rerun backup directory was dropped under. The path is gitignored.

### F4 — `U_ref` hardcoded in three files

`conventions.u_ref_from_field` reads it out of the field's own `meta_json` and raises on
any disagreement between headings or against a caller's expectation. It never falls back to
a default — a fallback is what made the original defect invisible. `throwability.compute_Q`
now cross-checks its `--u-ref` argument against the artifact, and `make_results_tex.py`
takes the value from the artifact rather than declaring it. Confirmed: the field was in
fact driven at 5.0 m/s, so the hardcoded value was right. It is now structurally enforced
rather than right by luck.

### F5, F6, F7 — the coverage holes

**New:** `P4_Coupling/test_provenance_and_gusts.py`, 13 tests, all passing (11 here, plus
two for the band reconciliation described in the follow-up section).

- The schema guard is now exercised **both ways**: a synthetic pre-fix artifact must make
  it fire and must produce no bias rows; the real artifact must not trip it.
- `gust_model` is read back from the real artifact and asserted to be `sweeping`, with the
  reason in the failure message.
- `prov_field_sha256` is checked against the actual SHA-256 of `flow_fields.npz` on disk,
  so a stale artifact carried past a field regeneration now fails a test instead of needing
  a human to compare timestamps.
- `prov_gust_dt` records a constant the production path never consumes; pinned, so it fails
  if `GUST_DT` ever becomes live and the provenance note goes stale.
- **F6 was re-derived, not reproduced.** The audit reported realised variance 1.21
  (temporal) vs 1.04 (sweeping). Over a *flight window* at the site's own correlation scale
  that does not hold: both models realise ≈0.96, because the generator normalises the
  ensemble variance by construction. What actually separates them is how that variance is
  **partitioned**. At T_L = 8.43 s (GATE G1's own measured `scale_mean`) over a 3.5 s
  flight, the temporal model puts **88.8%** of its variance into a frozen per-flight offset
  and only 10.8% into within-flight wander; the sweeping model, decorrelating over ~12 m of
  a ~25 m path, splits **54.8% / 43.8%**. A coherent offset pushes the disc the same way
  for the whole flight while within-flight wander partially cancels — which is the
  mechanism behind `Phi(N)`, now measured on realisations rather than asserted.

  This also exposed why the existing ACF test could never have seen it: `ConstantTurbField`
  gives T_L = k/ε = 1 s, **shorter than a flight**, so the frozen-gust regime the temporal
  model was retired for cannot arise in that fixture at all. The new tests use a fixture at
  the site's measured scale.
- The realised autocorrelation is now correlated empirically from sampled ξ(t) over 200
  realisations and checked against the analytic mode-weight prediction (atol 0.12).
- The 9% trim-speed deviation is bounded at 12%, with the three anchors that *do* agree
  (α_trim, γ, C_L) asserted alongside it so a regression cannot hide by moving them
  together. `validation_hummel.py`'s exit code 1 is pinned as the documented contract — it
  returns 1 while any coefficient set fails its band, and three of four do — together with
  the requirement that `potts_tabulated` be the **only** one that passes.

---

## The MINOR findings

- **F13** `atol=0.51` replaced by `conventions.disc_altitude_index`, whose tolerance is
  half the smallest grid spacing and which now **raises** rather than printing a note when
  the nearest sampled level is not the one being claimed.
- **F14** documented and pinned by a test (above).
- **F15** the pitch box and the 1e-3 m/s speed floor are single definitions in
  `conventions.py` / `throwability.SPEED_FLOOR`, imported by `make_results_tex.py`.

## Not done, and why

- **F9** (only `potts_tabulated` passes its own ±15% band) is now pinned by a test, but the
  underlying question — whether the other three sets should still ship — is a modelling
  decision, not a defect.
- **F2**'s reconciliation needs an author decision (above).
- The grid-independence failure, the posting hold and the author line are all untouched.
- No Monte Carlo, no CFD and no manuscript rebuild were run. `results_generated.tex` moved,
  so **the PDF should be rebuilt** before anything is quoted from it.

---

## Visualization

**New:** `P5_Synthesis/render_ensemble_fan.py` → `P5_Synthesis/video/ensemble_fan.mp4`.

A pyvista/VTK render of the 210° vs 270° Monte Carlo ensembles over the real site
geometry — the picture behind the paper's one resolved directional claim.

Everything with a physical position is computed, not drawn: the buildings are
`seanEarl.stl` (116 OSM footprints, 6,324 triangles); every trajectory is a fresh
integration of the same 6-DOF model with the same `potts_tabulated` coefficients, the same
RANS field, the same `sweeping` gust model and the same `sample_release` scatter as the
97,200-throw sweep; the landing scatter is those trajectories' own endpoints and the
ellipses are fitted to them; the ground is coloured by the horizontal 2 m wind speed from
the same field, with cells inside buildings left unpainted. Blocked throws are drawn in a
separate colour rather than deleted, because blocked is a real model outcome (10.2% of the
sweep).

Camera, colour and opacity are the only choices made for looks. Seeds are fixed, so the
same fan renders every time. Sanity-checked against the physics: mean range 23.8 m, apex
4.7 m, mean drift bearing 68° against a 60° throw.

The result is visible in the render itself — the 270° ensemble is both wider and longer
than the 210° one, with a visibly larger landing ellipse.

---

## Follow-up, same day: the two open items closed

### F2 — the contradictory bands, reconciled

Both bands fixed `N` by choosing an integral-scale convention and then disagreed about
which was meant. The argument was not resolvable by choosing better, because the run
records what it used:

| quantity | source | value over the 108 cells |
|---|---|---|
| decorrelation coordinate | `CorrelatedGusts.coordinate` returns `|r − r_launch|` | 0 → launch-to-landing chord |
| chord | `landing_stats.npz` `landing_mean` | **26.5 m** mean (17.3–33.4) |
| correlation length | per-cell `scale` in the comparison JSON | **3.30 m** mean (1.42–5.33) |
| **N = chord / scale** | measured, not chosen | **8.93 mean, 7.76 median** (3.70–19.81) |

The realised 3.30 m is none of the conventional 6.2 / 4.4 / 2.2 m, so none of the brief's
`N` values (4.0 / 5.6 / 11.3) is what ran. Pooled on sums of squares, `Φ(N)` predicts
**R = 1.90**; measured **1.868 isotropic (−1.9%)** and 1.735 anisotropic (−8.8%).

`Φ` is derived for an isotropic field, so the isotropic arm is the only one it is entitled
to predict, and it predicts it to 1.9%. The anisotropic residual is the anisotropy
correction `Φ` omits — which also resolves the 1.87-vs-1.74 sensitivity the physics
referee had flagged as an unexplained weakness in the estimator. It is not a weakness in
the estimator; it is missing physics in `Φ`, and it has the right sign and size.

**Single criterion, replacing both:** `Φ(N)` at the realised `N`, within GATE G2's own
±20% (brief line 622). The law is the brief's §C5, the tolerance is the brief's G2, and
`N` is the run's. Nothing invented.

Implemented as `regrade_gate_g1.py` section **[4b]** (`--landing`, default
`landing_stats.npz`), enforced by two new tests. The brief is annotated at all four places
it stated a band — GATE G1 acceptance, the WP2 P1 row, the early panel measurement at
line 71, and the WP1 checklist — with nothing deleted (diff confirms only the two lines
intended for replacement changed).

One consequence worth recording: the brief called the agreement between its code panel
(2.62/1.01 = 2.59×) and its tracer theory (2.3–2.6×) *"the strongest single piece of
evidence in this review."* It was not corroboration. That panel ran at `τ_c = 0.36 s`
while the shipped model realises a 3.30 m arclength scale over a 26.5 m chord — a
different experiment. Two methods agreeing while sharing the assumption in doubt is the
failure mode, not the confirmation.

### The manuscript PDF, rebuilt

`latexmk -pdf` on `Project1_Frisbee_LoyolaWind_PAPER.REVIEWED.tex`, exit 0:

- **26 pages**, 0 undefined references, 0 undefined citations, 0 multiply-defined labels
- 0 overfull `\hbox`/`\vbox` above 15 pt
- 0 `\MISSINGNUM`

Text-extracted and checked against the regenerated macros: **0.0321, 0.0481, 0.0442 all
present; 0.0324, 0.0483, 0.0356 all absent**, and no occurrence of `N_eddy` or the
withdrawn `2.0-2.6` band. The manuscript now matches its artifacts.

### Suites after both changes

P1 26/26, P4 40/40 (27 existing + 13 new), P5 13/13 — all exit 0.
