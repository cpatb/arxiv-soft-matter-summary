"""Fit the structural family to the Horner 2020 blood data: steady curve for the steady
parameters, the RISE of the 300 -> 1 1/s step for the one kinetic parameter the steady
curve cannot see.

Why only the rise: the thesis attributes the decay after the peak to syneresis (wall
depletion, ``data/horner2020/SOURCE.md``), so the bulk signal in the step sheet is the
samples from the first (1 s) to the maximum (4-6 s; 16 s on donor D). The steady curve
fixes ``eta_inf``, ``beta`` and the product ``a = alpha t_r`` (derivation fact 1) and
carries no information on the relaxation time (fact 3), so the joint problem factorises:
steady fit first, then a one-parameter fit of ``tau`` at 1 1/s to the rise, with

    sigma(t) = gd eta_inf (1 + beta lam(t)),
    lam(t)   = lam_ss(1) + (lam_ss(300) - lam_ss(1)) exp(-t / tau),        (fact 5)

and ``t_r = tau (1 + a)``, ``alpha = a / t_r``, ``t95 = tau ln 20``. The uncertainty of
the steady parameters is NOT propagated into ``tau`` (stated, not done). Time zero is
the step; the sheet's first sample is at 1.00037 s; whether the ARES-G2 had settled at
1 1/s by then is not stated in the thesis and is a caveat on the first sample.

Pre-registered gate (NEXT_STEPS.md item 2, written before this script ran):
  G1  t95 in [3, 8] s AND tau `measurable` or `bounded_only` on >= 16/20 donors;
  G2  rise-fit rmse <= 2x the sample noise (sd of the last 10 samples after a linear
      detrend) on >= 16/20 donors: the shape gate, since 3-6 points cannot carry a
      runs test.
``tau``'s 95 % interval is a one-parameter profile of the log-stress SSE on a log grid
with the Bates-Watts threshold ``SSE_min (1 + F(1, n-1; 0.95)/(n-1))``; the verdict is
``measurable`` if both ends are inside TAU_BOX, ``bounded_only`` if one end hits it,
``not_identifiable`` if both do. Writes ``data/horner2020/structural_step.md``.

Post-hoc diagnostics (added after the first run, which gave G1 20/20 and G2 6/20; the
gates above were not touched): (a) the rise refitted without the 1 s sample, to see
whether an instrument transient carries the failure; (b) the rise fitted between two
DATA anchors instead of the steady fit's values -- ``sigma(inf) = steady stress at
1 1/s`` (the peak itself, by protocol) and ``sigma(0) = steady stress at 350 1/s / 350``
(the sheet has 350, not 300: a 17 % rate mismatch, stated) -- so the steady fit's own
misfit (4-14 % at 1 1/s) is out of the residual. Both are reported beside G2, never in
its place. Each donor's residual sign pattern is printed because a shape failure shows
up there before it shows up in an rmse.

Run:  python scripts/structural_step_horner.py
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy import optimize, stats

from biorheo.constitutive.structural import STRUCTURAL_STEADY, lambda_steady
from biorheo.identify import fit_rate_curve

DATA = Path(__file__).resolve().parent.parent / "data" / "horner2020"
GD_HIGH, GD_LOW = 300.0, 1.0        # the step in the sheet's name
RATE_INDEX_1 = 4                    # the 1 1/s row of the steady sheet, asserted below
TAU_BOX = (1e-2, 1e2)               # s; profile box, two decades either side of the 1 s sampling
N_TAU_GRID = 801
NOISE_TAIL = 10                     # samples used for the noise estimate
T95_BAND = (3.0, 8.0)               # s, gate G1
SHAPE_FACTOR = 2.0                  # gate G2
GATE_MIN_DONORS = 16
ANCHOR_RATE = 350.0                 # 1/s; post-hoc (b) takes sigma(0) from the sheet's 350 1/s point, there is no 300
LN20 = float(np.log(20.0))


def read_matrix(name: str, prefix: str) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    with open(DATA / name, newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    x = np.array([r[0] for r in rows[1:]], dtype=float)
    cols = {h.removeprefix(prefix): np.array([np.nan if r[i] == "" else float(r[i]) for r in rows[1:]])
            for i, h in enumerate(rows[0][1:], start=1)}
    return x, cols


def stress_rise(t: np.ndarray, tau: float, eta_inf: float, beta: float, a: float) -> np.ndarray:
    lam_ss, lam0 = lambda_steady(GD_LOW, a), lambda_steady(GD_HIGH, a)
    lam = lam_ss + (lam0 - lam_ss) * np.exp(-t / tau)
    return GD_LOW * eta_inf * (1.0 + beta * lam)


def fit_tau(t: np.ndarray, sigma: np.ndarray, eta_inf: float, beta: float, a: float):
    """One-parameter log-stress least squares plus a profile interval. Returns
    (tau, lo, hi, verdict, rmse_log)."""
    def sse(log_tau: float) -> float:
        r = np.log(sigma) - np.log(stress_rise(t, np.exp(log_tau), eta_inf, beta, a))
        return float(r @ r)
    grid = np.linspace(np.log(TAU_BOX[0]), np.log(TAU_BOX[1]), N_TAU_GRID)
    s = np.array([sse(g) for g in grid])
    k = int(np.argmin(s))
    lo_b, hi_b = grid[max(k - 1, 0)], grid[min(k + 1, grid.size - 1)]
    res = optimize.minimize_scalar(sse, bounds=(lo_b, hi_b), method="bounded",
                                   options={"xatol": 1e-10})
    log_tau, s_min = float(res.x), float(res.fun)
    n = t.size
    thr = s_min * (1.0 + stats.f.ppf(0.95, 1, n - 1) / (n - 1))
    inside = s <= thr
    # the accepted region is the connected run of grid points around the minimum
    lo_i = k
    while lo_i > 0 and inside[lo_i - 1]:
        lo_i -= 1
    hi_i = k
    while hi_i < grid.size - 1 and inside[hi_i + 1]:
        hi_i += 1
    lo_at, hi_at = lo_i == 0, hi_i == grid.size - 1
    verdict = ("not_identifiable" if lo_at and hi_at else
               "bounded_only" if lo_at or hi_at else "measurable")
    return (np.exp(log_tau), np.exp(grid[lo_i]), np.exp(grid[hi_i]), verdict,
            float(np.sqrt(s_min / n)))


def fit_tau_anchored(t: np.ndarray, sigma: np.ndarray, s_0: float, s_inf: float):
    """Post-hoc (b): ``sigma(t) = s_inf + (s_0 - s_inf) exp(-t/tau)`` with both endpoints from
    the steady sheet; one parameter, same log-stress objective and grid as :func:`fit_tau`.
    Returns (tau, rmse, signs)."""
    def sse(log_tau: float) -> float:
        r = np.log(sigma) - np.log(s_inf + (s_0 - s_inf) * np.exp(-t / np.exp(log_tau)))
        return float(r @ r)
    grid = np.linspace(np.log(TAU_BOX[0]), np.log(TAU_BOX[1]), N_TAU_GRID)
    k = int(np.argmin([sse(g) for g in grid]))
    res = optimize.minimize_scalar(sse, bounds=(grid[max(k - 1, 0)], grid[min(k + 1, grid.size - 1)]),
                                   method="bounded", options={"xatol": 1e-10})
    tau = float(np.exp(res.x))
    pred = s_inf + (s_0 - s_inf) * np.exp(-t / tau)
    return tau, float(np.sqrt(np.mean((pred - sigma) ** 2))), "".join("+" if r > 0 else "-" for r in sigma - pred)


def main() -> None:
    t, step = read_matrix("step_300_to_1_human.csv", "stress_Pa_")
    rates, steady = read_matrix("steady_shear_human.csv", "stress_Pa_")
    assert rates[RATE_INDEX_1] == GD_LOW and set(step) == set(steady)
    print(f"step sheet: {t.size} samples, {len(step)} donors; steady sheet: {rates.size} rates")
    lines = ["# Structural family on blood: steady curve + step rise", "",
             "Generated by `scripts/structural_step_horner.py`. Steady columns: the "
             "`structural_steady` fit (Cross with m = 1) to the 16-point curve, log-eta residuals. "
             "Rise columns: the one-parameter fit of tau to the step samples up to the maximum "
             "(n_rise of them), holding the steady parameters. `t_r = tau (1 + a)`, "
             "`t95 = tau ln 20`. `noise` is the sd of the last 10 samples after a linear detrend; "
             "`rmse` is the rise fit's in mPa, `rmse/amp` the same over the model's rise amplitude "
             "sigma(inf) - sigma(0). `signs` is sign(data - model) per rise sample. "
             "G1: t95 in [3, 8] s and tau measurable or bounded_only. G2: rmse <= 2 noise.", "",
             "| donor | steady rmse_log | structured | eta_inf mPa s | beta | a s | n_rise | "
             "tau s [95 %] | verdict | t_r s | alpha | t95 s | rmse mPa | rmse/amp % | noise mPa | signs | G1 | G2 |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    n_g1 = n_g2 = n_struct = 0
    n_g2_drop1 = n_g2_anchored = n_anchored_pattern = n_pattern = 0
    t95s, taus, ramps, anchored_rows = [], [], [], []
    i350 = int(np.argmin(np.abs(rates - ANCHOR_RATE)))
    assert rates[i350] == ANCHOR_RATE
    for d in sorted(step):
        eta = steady[d] / rates
        st = fit_rate_curve(STRUCTURAL_STEADY, rates, eta)
        eta_inf, beta, a = st.values
        n_struct += st.residual_structure.structured
        v = step[d]
        ok = np.isfinite(v)
        vv, tt = v[ok], t[ok]
        imax = int(np.argmax(vv))
        t_rise, s_rise = tt[:imax + 1], vv[:imax + 1]
        tau, lo, hi, verdict, _ = fit_tau(t_rise, s_rise, eta_inf, beta, a)
        pred = stress_rise(t_rise, tau, eta_inf, beta, a)
        rmse = float(np.sqrt(np.mean((pred - s_rise) ** 2)))
        tail_t, tail_s = tt[-NOISE_TAIL:], vv[-NOISE_TAIL:]
        noise = float(np.std(tail_s - np.polyval(np.polyfit(tail_t, tail_s, 1), tail_t), ddof=2))
        t_r, alpha, t95 = tau * (1.0 + a), a / (tau * (1.0 + a)), tau * LN20
        g1 = T95_BAND[0] <= t95 <= T95_BAND[1] and verdict in ("measurable", "bounded_only")
        g2 = rmse <= SHAPE_FACTOR * noise
        n_g1 += g1
        n_g2 += g2
        amp = GD_LOW * eta_inf * beta * (lambda_steady(GD_LOW, a) - lambda_steady(GD_HIGH, a))
        signs = "".join("+" if r > 0 else "-" for r in s_rise - pred)
        n_pattern += signs[0] == "+" and signs[1:3] == "--"
        # post-hoc (a): without the 1 s sample
        tau1, *_ = fit_tau(t_rise[1:], s_rise[1:], eta_inf, beta, a)
        rmse1 = float(np.sqrt(np.mean((stress_rise(t_rise[1:], tau1, eta_inf, beta, a) - s_rise[1:]) ** 2)))
        n_g2_drop1 += rmse1 <= SHAPE_FACTOR * noise
        # post-hoc (b): both endpoints from the steady sheet, one-parameter exponential between them
        s_inf, s_0 = steady[d][RATE_INDEX_1], steady[d][i350] / ANCHOR_RATE * GD_LOW
        tau_b, rmse_b, signs_b = fit_tau_anchored(t_rise, s_rise, s_0, s_inf)
        g2_b = rmse_b <= SHAPE_FACTOR * noise
        n_g2_anchored += g2_b
        n_anchored_pattern += signs_b.startswith("+-") and signs_b.endswith("+")
        anchored_rows.append(f"| {d} | {s_0 * 1e3:.2f} | {s_inf * 1e3:.2f} | {tau_b:.2f} | {rmse_b * 1e3:.2f} | "
                             f"{signs_b} | {'pass' if g2_b else 'FAIL'} |")
        t95s.append(t95)
        ramps.append(rmse / amp * 100)
        taus.append(tau)
        lines.append(f"| {d} | {st.fit.rmse_log:.3f} | {'yes' if st.residual_structure.structured else 'no'} | "
                     f"{eta_inf * 1e3:.2f} | {beta:.1f} | {a:.2f} | {t_rise.size} | "
                     f"{tau:.2f} [{lo:.2f}, {hi:.2f}] | {verdict} | {t_r:.2f} | {alpha:.3f} | {t95:.1f} | "
                     f"{rmse * 1e3:.2f} | {rmse / amp * 100:.1f} | {noise * 1e3:.2f} | {signs} | {'pass' if g1 else 'FAIL'} | "
                     f"{'pass' if g2 else 'FAIL'} |")
    t95s, taus, ramps = np.array(t95s), np.array(taus), np.array(ramps)
    lines += ["", "## Counts", "",
              f"- Steady `structural_steady` fit flagged structured: {n_struct}/20 (Cross with m = 1 is a "
              f"misfit where free-m Cross was not; compare model_comparison.md).",
              f"- G1 (t95 in [{T95_BAND[0]:.0f}, {T95_BAND[1]:.0f}] s, tau determined): {n_g1}/20 against "
              f"a pre-registered {GATE_MIN_DONORS}/20 -> {'PASS' if n_g1 >= GATE_MIN_DONORS else 'FAIL'}.",
              f"- G2 (rise rmse <= {SHAPE_FACTOR:.0f}x sample noise): {n_g2}/20 against "
              f"{GATE_MIN_DONORS}/20 -> {'PASS' if n_g2 >= GATE_MIN_DONORS else 'FAIL'}. Residual sign "
              f"pattern `+--...` (data above the model at 1 s, below at 2 and 3 s) on {n_pattern}/20.",
              f"- tau: min {taus.min():.2f}, median {np.median(taus):.2f}, max {taus.max():.2f} s; "
              f"t95: min {t95s.min():.1f}, median {np.median(t95s):.1f}, max {t95s.max():.1f} s.", "",
              "## Post-hoc diagnostics (gates above unchanged)", "",
              f"- (a) Without the 1 s sample, G2 passes on {n_g2_drop1}/20.",
              f"- (b) Rise fitted between two data anchors, sigma(0) = steady stress at {ANCHOR_RATE:.0f} 1/s "
              f"/ {ANCHOR_RATE:.0f} (the sheet has no 300 1/s point) and sigma(inf) = the steady 1 1/s "
              f"value: G2 passes on {n_g2_anchored}/20; the residual sign pattern is `+-...+` (data above "
              f"the exponential at 1 s and at the peak, below between) on {n_anchored_pattern}/20.", "",
              "| donor | sigma(0) mPa | sigma(inf) mPa | tau s | rmse mPa | signs | G2 |",
              "|---|---|---|---|---|---|---|", *anchored_rows, "",
              "## Reading", "",
              "The steady curve alone leaves alpha and t_r degenerate (tests/test_structural.py); the "
              "rise samples carry tau on its own, so the pair is separated by construction whenever "
              "tau's profile closes inside the box. What the gate tests is whether the separated "
              "values are physically sensible (t95 at the observed peak time) and whether a single "
              "exponential is the shape of the rise at the sample noise. Caveats: the steady "
              "parameters' uncertainty is not propagated; the first sample at 1 s may include the "
              "instrument's own rate transient; the decay after the peak is excluded on the thesis's "
              "syneresis reading, which is not tested here.", "",
              "Result of the first run (2026-09-05): G1 passed on every donor and G2 failed. The "
              "timescale is right -- the model's approach to steady state lands where the peak is -- "
              "but a single exponential is not the shape of the rise at the sample noise: the fitted "
              f"rise misses by {ramps.min():.1f}-{ramps.max():.1f} % of its own amplitude (median "
              f"{np.median(ramps):.1f} %) where the tail noise is 0.04-1 mPa, "
              "and the sign pattern is systematic. Neither post-hoc variant changes that, and the "
              "data-anchored one removes the steady fit's misfit at 1 1/s from the residual and still "
              "fails, so the failure is in the kinetics, not inherited from the steady curve. The rise "
              "is less concave than an exponential between the same endpoints (nearer to linear in t), "
              "which first-order kinetics with one structure variable cannot produce. What would: a "
              "lag or a second timescale (two-structure or strain-dependent kinetics), or an instrument "
              "transient at 1 s that the drop-first variant does not fully remove. Untested here.", ""]
    (DATA / "structural_step.md").write_text("\n".join(lines), encoding="utf-8", newline="\n")
    print("\n".join(l for l in lines if l.startswith("- ")))


if __name__ == "__main__":
    main()
