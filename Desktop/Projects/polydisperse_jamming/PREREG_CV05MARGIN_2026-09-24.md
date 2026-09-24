# PREREG — `cv05margin`: can a CV-0.5 suspension jam at a weighed φ where the CV-0 one still flows?

Written 2026-09-24, before any cell (commit time is the record). Script `margin_cv05.py` (extended to read the new
seeds explicitly); launcher `scripts/launch_cv05margin.sh`.

## The question, restated (supersedes `predict_cvmap.txt` §1's propagation)
At fixed weighed φ, CV 0.5 can shear-jam iff φ ≥ φ_m(CV 0.5, μ 0.3), and CV 0 flows at rest iff φ < φ₀(CV 0, μ 0).
Such a φ exists iff **D = φ_m(0.5, μ 0.3) − φ₀(0, μ 0) < 0**. φ₀(CV 0.5) cancels, so its single converged seed
never mattered. `predict_cvmap.py` propagated it anyway (the "1.2×" figure); that figure is withdrawn.
On disk (n = 3, seed-paired): D = +0.00485, per seed +0.0067 / +0.0081 / −0.0003. t = 1.88 against t_crit 4.30,
so NOT RESOLVED (`out/margin_cv05.txt`).

## Cells
Two (μ, CV) cells, lognormal, N = 1500, rung p̂ 3e-5, PJ_MAX_ITER 1e6 (the cap of the existing records):
(μ_p 0.3, CV 0.5) and (μ_p 0, CV 0), seeds 4–9 (NEW), plus seed 1 as a wiring control.
Output is pulled to `out/cv05margin/`, NOT to `out/`: `three_knobs.load()` globs `out/phase2*.jsonl`, so
dropping the files there would silently change every analysis that reads it.

## Gates and decision rule
- **W1 (wiring):** the seed-1 re-runs reproduce the existing φ of both cells to |Δφ| ≤ 1e-10 (same code era;
  the md5 of phase2_phiJ.py/packing.py/measure.py is identical repo vs WSL apart from line endings).
- Only converged records count. Every non-converged new seed is printed and excluded, never imputed.
- **Decision (n = all converged paired seeds, 1–9):** t = mean D / max(se_paired, PHI_SCALE·√(2/n)) against
  t_crit(df n−1).
  **RESOLVED-NO-OVERLAP** if t > t_crit: raising CV to 0.5 at fixed weighed φ removes shear-jamming for EVERY
  suspension that flows at rest at CV 0 (μ_m 0.3).
  **RESOLVED-OVERLAP** if t < −t_crit: a band [φ_m(0.5), φ₀(0)) exists where both hold.
  **UNDETERMINED** otherwise. Then |D| is below what this protocol resolves, and the honest statement is that
  CV 0.5 moves the frictional jamming point to within ±|D| of monodisperse frictionless jamming.
- Expected (not graded): D > 0 as a point estimate, from the n = 3 data.
- **Observable at a pressure, as `phase2_phiJ.py --phats` requires stating:** φ here is read AT p̂ 3e-5, not
  extrapolated to φ_J. That is the convention of the window endpoints in `jam_map.py` and `predict_cvmap.py`
  (§1 "window endpoints at the quoted rung"). The question is therefore about those rung-defined endpoints.
  The rung-convention caveat that already travels with every χ (ANSWER §6) applies unchanged.
