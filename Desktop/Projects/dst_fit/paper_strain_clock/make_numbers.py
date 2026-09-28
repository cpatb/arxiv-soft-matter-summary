"""Write numbers.tex: every number the strain-clock note quotes, computed here.

Nothing in main.tex is typed by hand. Sources: the gate's own functions (imported from
scripts/gate_g11.py), the figure computations, and the direct-control probe.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
for p in (ROOT / "src", ROOT / "scripts", ROOT / "figures"):
    sys.path.insert(0, str(p))

import fig_clock as FC  # noqa: E402
import fig_reversal as FR  # noqa: E402
import gate_g11 as G  # noqa: E402
from dst_fit.transient import (  # noqa: E402
    admissible_window,
    d_reduced_jamming_distance,
    fit_relaxation,
    frictional_fraction_from_viscosity,
    log_bin,
)


def fmt(x, nd=2):
    return f"{x:.{nd}f}"


def main() -> None:
    m = {}
    rev = FR.compute()
    cl = FC.clocks()
    g, curves = FC.control()
    for panel, tag in (("Fig_2b", "Si"), ("Fig_2d", "Pm")):
        d = rev[panel]
        m[f"phi{tag}"] = fmt(d["phi"], 3)
        m[f"phim{tag}"] = fmt(d["phi_m"], 3)
        m[f"floor{tag}"] = fmt(d["floor"], 2)
        m[f"sweepMin{tag}"] = fmt(d["sweep_min"], 2)
        # depth of the SMOOTHED minimum below the floor, as the gate computes it
        depths = []
        for rate, gg, ee in d["records"]:
            try:
                w = admissible_window(gg, ee, d["phi"], d["phi_m"], G.PHI0_DEFAULT)
                depths.append(w.depth_below_floor)
            except ValueError:
                pass  # record entirely below the floor (PMMA 0.05)
        m[f"depthLo{tag}"] = fmt(min(depths), 2)
        m[f"depthHi{tag}"] = fmt(max(depths), 2)
        m[f"nDepth{tag}"] = str(len(depths))
        rows = [x for x in cl[panel] if x["included"]]
        full = np.array([x["full"] for x in rows])
        g50 = np.array([x["g50"] for x in rows])
        rates = np.array([x["rate"] for x in rows])
        m[f"fullLo{tag}"], m[f"fullHi{tag}"] = fmt(full.min()), fmt(full.max())
        m[f"gFiftyLo{tag}"], m[f"gFiftyHi{tag}"] = fmt(g50.min()), fmt(g50.max())
        t50 = g50 / rates
        m[f"gFiftySpread{tag}"] = fmt(g50.max() / g50.min())
        m[f"tFiftySpread{tag}"] = fmt(t50.max() / t50.min(), 1)
        m[f"rateSpan{tag}"] = f"{rates.max() / rates.min():.0f}"
        m[f"nRates{tag}"] = str(len(rows))
        allv = np.concatenate([full, [x["trunc"] for x in rows], g50])
        m[f"allLo{tag}"], m[f"allHi{tag}"] = fmt(allv.min()), fmt(allv.max())
    si = [x for x in cl["Fig_2b"] if x["included"]]
    m["fullFallSi"] = fmt(si[0]["full"] / si[-1]["full"])
    m["ctlDirect"] = fmt(curves["direct"][3], 3)
    m["ctlInverted"] = fmt(curves["inverted"][3], 3)
    m["ctlDiffPct"] = f"{100 * (curves['inverted'][3] - curves['direct'][3]) / curves['inverted'][3]:.0f}"
    # --- referee 2026-09-27, finding 2: depth below the steady sweep, from the SMOOTHED
    # minima the gate uses, as a percentage drop below the lowest steady viscosity.
    drops = []
    for panel in ("Fig_2b", "Fig_2d"):
        d = rev[panel]
        for rate, gg, ee in d["records"]:
            try:
                w = admissible_window(gg, ee, d["phi"], d["phi_m"], G.PHI0_DEFAULT)
            except ValueError:
                continue
            drops.append(1.0 - w.eta_min / d["sweep_min"])
    m["dropLoPct"] = f"{100 * min(drops):.0f}"
    m["dropHiPct"] = f"{100 * max(drops):.0f}"

    # --- referee finding 1: sensitivity of the clock to WHERE the fit starts. The start is
    # the binned reversal minimum of a broad trough (probe_rebuild_lag.py showed a delay
    # measured from it is ill-posed); move it anywhere within +SHIFT strain AFTER it and
    # refit. Starts BEFORE the minimum are not scanned into the quoted number: they include
    # the still-falling viscosity, which the relaxation law does not describe, and on
    # silica 5 1/s the fit then has no interior optimum (scripts/probe_start_sensitivity.py,
    # logs/probe_start_sensitivity.txt). That is stated in the text, not hidden here.
    SHIFT = 0.2  # strain; covers the ~0.15 trough ambiguity probe_rebuild_lag found
    worst = 1.0
    lam_range = {}  # (panel, rate) -> (min, max) over after-minimum starts
    for panel, steady in G.PAIRS:
        d = rev[panel]
        included = {x["rate"] for x in cl[panel] if x["included"]}
        for rate, gg, ee in d["records"]:
            if rate not in included:
                continue
            gb, eb, sb = log_bin(gg, ee)
            i0 = int(np.argmin(eb))

            def clock_from(j):
                f = frictional_fraction_from_viscosity(eb[j:], d["phi"], d["phi_m"], G.PHI0_DEFAULT)
                fe = (d["phi"] / (G.PHI0_DEFAULT - d["phi_m"])
                      * np.abs(d_reduced_jamming_distance(eb[j:])) * sb[j:])
                return fit_relaxation(gb[j:], f, fe).lam

            base = clock_from(i0)
            ref = [x["full"] for x in cl[panel] if x["rate"] == rate][0]
            assert abs(base - ref) < 1e-12, (panel, rate, base, ref)  # same code path as the gate
            lams = []
            for j in range(i0, gb.size - 12):
                if gb[j] - gb[i0] > SHIFT:
                    break
                lam = clock_from(j)
                lams.append(lam)
                r = max(lam / base, base / lam)
                if r > worst:
                    worst = r
            lam_range[(panel, rate)] = (min(lams), max(lams))
    m["startShiftRange"] = fmt(SHIFT)
    m["startShiftWorst"] = fmt(worst)
    # does the silica fall with rate survive the start ambiguity? lower bound on the fall
    # = smallest clock at the lowest rate over the largest at the highest rate
    si_rates = sorted(r for (p, r) in lam_range if p == "Fig_2b")
    m["fallLowerSi"] = fmt(lam_range[("Fig_2b", si_rates[0])][0]
                           / lam_range[("Fig_2b", si_rates[-1])][1])
    ordered = all(lam_range[("Fig_2b", a)][0] > lam_range[("Fig_2b", b)][1]
                  for a, b in zip(si_rates, si_rates[1:]))
    # rates whose ranges are disjoint from the next one: the part of the trend that survives
    sep = [f"{a:g}/{b:g}" for a, b in zip(si_rates, si_rates[1:])
           if lam_range[("Fig_2b", a)][0] > lam_range[("Fig_2b", b)][1]]
    print("silica after-min clock ranges:",
          {r: tuple(round(v, 3) for v in lam_range[("Fig_2b", r)]) for r in si_rates},
          "fully ordered:", ordered, "separated pairs:", sep)
    # main.tex (Results, "on silica the clock falls") states in PROSE that 0.5 > 1 > 2 stay
    # ordered and 2 vs 5 overlap. Fail loudly if the data ever stops saying that.
    assert sep == ["0.5/1", "1/2"] and not ordered, (
        f"separated pairs now {sep}: rewrite the ordering sentence in main.tex")

    # --- referee findings 3/4: the two-component numbers, parsed from the probe's output
    out = subprocess.run([sys.executable, str(ROOT / "scripts" / "probe_two_clock.py")],
                         capture_output=True, text=True, cwd=ROOT, check=True).stdout
    spreads = re.findall(r"spread ([0-9.]+)x", out)
    assert len(spreads) >= 2, out[:500]
    m["twoClockMOne"], m["twoClockMFour"] = spreads[0], spreads[1]
    # --- second protocol: constant-stress start-up on silica (Richards 2019 deposit);
    # PREREG_startup_clock_2026-09-27.md, P1/P3/P4 FAILED as registered.
    import probe_startup_p4 as S4  # noqa: E402
    jt = S4.jam_times()
    rows_j, onset_j = S4.jam_table()
    m["suOnset"] = f"{onset_j:.1f}"
    m["suNAbove"] = str(len(jt))
    m["suJamLo"], m["suJamHi"] = fmt(jt[:, 1].min()), fmt(jt[:, 1].max())
    m["suTjLo"], m["suTjHi"] = fmt(jt[:, 2].min()), fmt(jt[:, 2].max())
    m["suJamSpread"] = fmt(jt[:, 1].max() / jt[:, 1].min())
    m["suTjSpread"] = fmt(jt[:, 2].max() / jt[:, 2].min())
    m["suSigSpread"] = fmt(jt[:, 0].max() / jt[:, 0].min())
    ab_j = [r for r in rows_j if r["above"]]  # referee pass 3: scope to above onset
    m["suJamAboveR"] = f"{sum(r['gamma_j'] > r['R'] for r in ab_j)}/{len(ab_j)}"
    nz, nab, _ = S4.zero_crossings()
    m["suZeroCross"] = f"{nz}/{nab}"
    # fit-free landmark against fit-free landmark: smallest reversal gamma_50 on Lin silica
    # over the largest start-up jam strain on Richards silica (different samples)
    g50_si = [x["g50"] for x in cl["Fig_2b"] if x["included"]]
    m["reversalOverJam"] = f"{min(g50_si) / jt[:, 1].max():.1f}"
    m["solverGammaZero"] = f"{G.SOLVER_GAMMA_0:g}"
    allv = [float(v) for k in ("allLoSi", "allLoPm") for v in [m[k]]]
    m["minOverSolver"] = fmt(min(allv) / G.SOLVER_GAMMA_0, 1)

    lines = [f"% generated by {Path(__file__).name}; do not edit"]
    for k, v in sorted(m.items()):
        lines.append(f"\\newcommand{{\\{k}}}{{{v}}}")
    (HERE / "numbers.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")
    for k, v in sorted(m.items()):
        print(f"{k:16s} {v}")


if __name__ == "__main__":
    main()
