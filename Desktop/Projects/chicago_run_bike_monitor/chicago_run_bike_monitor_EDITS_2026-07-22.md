# Edits — chicago_run_bike_monitor (2026-07-22)

**Targets chosen because:** ⚠ **this project does not pass the edit brief's own selection gate**, and was edited anyway on explicit instruction. Recommended stage is *experimental validation* (score 3); preprint and journal both score 1. But the edit brief names the H2 σ_y → σ_z reframe by name in its Tier-B list. I treated the specific naming as overriding the general filter, because the correction has a hard deadline that the filter does not see: **fixing it after a field campaign means re-running 20+ transects in matched meteorology.** Say the word and I will revert.

**Files edited (original → copy):** `chicago_run_bike_monitor.TEX` → `chicago_run_bike_monitor.EDITED.TEX`

Original untouched (verified: 46,102 B, 2026-06-30, unchanged).

## Tier A (mechanical) — applied

| File | Location | Before → After | Why |
|---|---|---|---|
| `chicago_run_bike_monitor.EDITED.TEX` | l.6 header comment | `Compile with: pdflatex chicago_airquality_project.tex` → `… chicago_run_bike_monitor.EDITED.TEX` | `AUDIT_papers.md` Paper-3 issue #5. The named file does not exist in the tree. |

## Tier B (physics) — applied, needs your eyes

| File | Location | Before → After | Basis | Confidence |
|---|---|---|---|---|
| `.EDITED.TEX` | §H2 hypothesis, l.~190–208 | Crosswind Gaussian `C(y) = Q/(√(2π)σ_y u)·exp(−y²/2σ_y²)`, fitted width interpreted as σ_y and compared to PG curves → ground-level infinite-line-source relation `C(x) = √(2/π)·Q_L/(σ_z(x)·u)` with `σ_z(x) = a x^b`, exponent compared to PG class-D σ_z | `AUDIT_papers.md` Paper-3 issue #3. For an effectively infinite crosswind line source the concentration **does not fall off in the crosswind direction at all** — the source is translation-invariant along its own axis. What a perpendicular transect measures is decay with *downwind* distance, governed by *vertical* dispersion σ_z(x): the plume dilutes upward, not sideways. Standard line-source Gaussian plume theory (Turner, *Workbook of Atmospheric Dispersion Estimates*; PG σ_z parameterisation). | **High** on the diagnosis (σ_y is unambiguously the wrong coordinate); **Medium** on this specific replacement — see the flag below |
| `.EDITED.TEX` | §Level 3 analysis pipeline, l.~750–759 | Bin by "cross-wind distance from the highway centerline", fit `C(y)=A·exp(−y²/2σ_y²)` via `curve_fit`, extract σ_y → bin by **downwind** distance x from the road edge, fit `C(x) = C₀a⁻¹x^(−b)` in log–log space, extract a and b, convert prefactor to per-unit-length Q_L | Same basis; the pipeline had to follow the hypothesis. The old recipe conflated the two coordinates outright — distance from the centerline *along a perpendicular transect* **is** the downwind coordinate x, not the crosswind coordinate y. | High |
| `.EDITED.TEX` | §Level 3, added items | Added: restrict fits to near-perpendicular wind; fit in log–log space; **add a background term** `C(x) = C_bg + C₀a⁻¹x^(−b)` or subtract a control-site background | Two latent analysis bugs, not just a mislabelled axis. (i) Fitting `exp(−y²/2σ²)` to a monotonically decaying transect with no interior peak is poorly conditioned — `curve_fit` would routinely return an arbitrary width with a plausible-looking covariance. (ii) A power law fitted without a background term biases the exponent low, which would corrupt the PG comparison that is the whole point of H2. | High |

## NOT edited — needs author decision

- **⚠ Infinite-line vs finite-line idealisation.** The applied correction adopts the **infinite**-line form, valid when the road is long compared with the sampling distance (fine for I-90/94 out to a few hundred metres; questionable near interchanges or a road end). The alternative is to treat the highway as a **finite** line and integrate along it explicitly, which retains a genuine crosswind dependence and would justify a two-parameter (σ_y, σ_z) fit. The infinite-line form is the standard default and the document already calls the highway "a continuous line source", so I went with it — **but this is a modelling commitment you should make deliberately, not inherit from an edit pass.** Flagged inline in the file. This is why confidence is Medium rather than High.
- **The PG class-D σ_z exponent b ≈ 0.7–0.8** is written into the edited text as the expected range over 10–300 m, with an explicit "verify against the specific PG parameterisation used" instruction attached. PG σ_z curves differ between the Briggs, Turner, and ISC formulations; I did not want to pin a number without knowing which you will use.
- **H1, H3, H4 remain "measure and see" rather than falsifiable predictions.** Reframing them quantitatively is net-new scientific authoring, well outside an edit pass.
- **No co-location calibration protocol exists.** This is arguably the single biggest determinant of whether the open dataset has quantitative value, but writing one is net-new content. Raised in `READINESS.md`.
- **"Five EPA monitoring stations"** — unverified factual claim. Not edited; I cannot confirm the station count from here, and guessing would violate the no-fabrication guardrail.
- **"PM2.5 R² > 0.95 vs FEM" for the SEN5x** — optimistic and strongly condition-dependent. Tempering it requires citing specific evaluations I have not read. Left flagged.

## Compile / run check

| File | Passes | Exit | Pages (orig → edited) | Undefined refs/cites |
|---|---|---|---|---|
| `chicago_run_bike_monitor.EDITED.TEX` | 2 | 0, 0 | 20 → **20** | 0 |

## Remaining shortest path to the recommended stage (experimental validation)

1. **Decide infinite-line vs finite-line**, and confirm which PG σ_z parameterisation you will compare against. Everything else in H2 follows.
2. Write the co-location calibration protocol (park beside an EPA reference station for 24–48 h before and after the campaign; report the regression).
3. Verify the EPA station count; temper the SEN5x accuracy claim.
4. Then run the campaign — with the control site, since without it a fall-off cannot be distinguished from a background gradient or sensor drift.
