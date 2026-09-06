"""Can the Horner 2020 workbook's OSCILLATORY sheets discriminate a lag from a second
timescale in the blood step-shear rise? (NEXT_STEPS item 7.)

This is a capability measurement, not a fit. Item 2 fitted first-order structural
kinetics to the rise of the 300 -> 1 1/s step: the timescale gate passed 20/20 (t95 in
3-8 s, so tau ~ 1-2.7 s) and the SHAPE gate FAILED 6/20 -- the rise is less concave than
a single exponential, with one residual sign pattern on 15/20. Two readings survive that
and a 4-6 sample rise cannot separate them: a LAG before build-up starts, or a SECOND,
faster process in series. Item 7 asks whether the workbook's oscillatory sheets can.

NEXT_STEPS item 7 framed the question as "do the LAOS sheets have enough samples per
cycle". They do, by a wide margin, so that is NOT what decides it; three other things
do. All five criteria below were written before any of them was graded.

PRE-REGISTERED CRITERIA (2026-09-05; a failure is recorded, never softened):

  C1 POPULATION. A discrimination is only a corroboration of item 2 if it speaks for the
     same donors. Bar: >= 16 of the 20 human donors present on the sheets the
     discrimination would use -- the same 80 % threshold item 2's own G1/G2 used.
  C2 THE LINEAR CHANNEL IS A MEASUREMENT. On the Frequency Sweep, G' must be positive at
     every frequency for a donor to count (a negative storage modulus is not a
     measurement of stored elastic energy), on >= 16 of 20 donors.
  C3 LINEARITY IS DEMONSTRATED, NOT ASSUMED. The textbook signatures of a lag and of a
     second timescale are LINEAR-response statements, so they apply only if the
     frequency sweep sits inside the linear viscoelastic region. Bar: the Amplitude
     Sweep must reach a strain at most HALF the frequency sweep's own operating strain,
     so that a plateau containing that strain can actually be shown from below.
  C4 THE DISCRIMINATING CHANNEL EXISTS. A lag theta and a second timescale differ chiefly
     in PHASE: a lag contributes a phase -w*theta that is linear in w without a
     corresponding change in the modulus, a second Debye process contributes a bounded
     phase excursion with a matching dispersion in the modulus. Bar: the workbook must
     carry a loss modulus, a phase angle, or a raw stress waveform aligned to the strain
     from which one can be computed, for the sheets C1-C3 admit.
  C5 SAMPLING (item 7's own stated question). Bar: >= 10 samples per cycle on the LAOS
     sheets, which is the coarsest sampling at which a waveform's phase is recoverable
     at all.

VERDICT RULE, fixed in advance: the linear-response route is CLOSED if C2, C3 or C4
fails. The nonlinear (LAOS) route cannot substitute for it if C1 fails, because a
LAOS-only answer would speak for a different set of donors than item 2's finding does.
If the route is closed, item 7 is closed the way item 2b was -- as not runnable from this
data -- and what a discriminating measurement WOULD need is written down instead.

Source: `data/horner2020/Horner-Thesis-Data.xlsx` (gitignored; URL and md5 in
`data/horner2020/SOURCE.md`). Run: python scripts/oscillatory_capability_horner.py
Writes `data/horner2020/oscillatory_capability.md`. Pinned in
`tests/test_oscillatory_capability.py`.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import openpyxl

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data" / "horner2020"
WORKBOOK = DATA / "Horner-Thesis-Data.xlsx"

# The 20 human donors of the steady-shear sheet and of every blood record in the
# registry. E is absent from the source workbook; the letters are the source's own.
HUMAN_DONORS = tuple("ABCDFGHIJKLMNOPQRSTU")
COVERAGE_BAR = 16          # C1, C2: of 20, the threshold item 2's gates used
MIN_ROWS_FOR_A_SERIES = 50  # a donor column with fewer numeric cells is not a series
MIN_ROWS_SWEEP = 5          # the sweeps have 40-52 rows; 5 is "this column has data"
SAMPLES_PER_CYCLE_BAR = 10  # C5
LVE_HEADROOM = 0.5          # C3: the amplitude sweep must reach <= 0.5x the freq strain

# Read Me, verbatim: "The representative amplitude sweep was conducted at a frequency of
# 2 Hz and the frequency sweep was conducted at a strain amplitude of 10". The unit of
# that 10 is AMBIGUOUS in the source -- the Amplitude Sweep axis is in %, while the LAOS
# sheet labels are dimensionless by the Read Me's own statement -- so both readings are
# carried and C3 is graded under each. It is not resolved here and must not be guessed.
FREQ_SWEEP_STRAIN_READINGS = {"10 % (Amplitude Sweep's own unit)": 10.0,
                              "10 dimensionless = 1000 % (the LAOS label unit)": 1000.0}


def _sheet_rows(wb, name: str) -> list[tuple]:
    return list(wb[name].iter_rows(values_only=True))


def _donor_columns(header: tuple) -> dict[str, int]:
    """Column index per human donor. The header row labels donors by a single letter;
    animal columns ('Guinea Pig', 'Cow', ...) are longer and are excluded here, because
    every claim in this repo is about the 20 human donors."""
    return {c: i for i, c in enumerate(header)
            if isinstance(c, str) and c in HUMAN_DONORS}


def _numeric(rows: list[tuple], col: int) -> np.ndarray:
    return np.array([r[col] for r in rows if isinstance(r[col], (int, float))], dtype=float)


def _paired(rows: list[tuple], col: int) -> tuple[np.ndarray, np.ndarray]:
    """(x, y) keeping only rows where BOTH the axis and the donor column are numeric.
    Filtering the two columns independently and truncating to the shorter one would
    MISALIGN them the moment a donor column has a gap the axis does not."""
    both = [(r[0], r[col]) for r in rows
            if isinstance(r[0], (int, float)) and isinstance(r[col], (int, float))]
    if not both:
        return np.empty(0), np.empty(0)
    a = np.array(both, dtype=float)
    return a[:, 0], a[:, 1]


def _blocks(x: np.ndarray) -> list[slice]:
    """Both sweep sheets are TWO CONCATENATED REPEATS of the same ramp, not one sweep:
    the Amplitude Sweep is 26 + 26 points from 9.95 %, the Frequency Sweep 20 + 20 from
    0.5 rad/s. An argmax over the whole column therefore lands in whichever repeat
    happens to be larger, and mixes two runs. Split on every point where the axis stops
    increasing."""
    cuts = [0] + [j for j in range(1, x.size) if x[j] <= x[j - 1]] + [x.size]
    return [slice(cuts[i], cuts[i + 1]) for i in range(len(cuts) - 1)]


def laos_census(wb) -> list[dict]:
    """One row per LAOS / UD-LAOS sheet: amplitude, frequency, sampling, donors."""
    out = []
    for name in wb.sheetnames:
        m = re.match(r"(UD-)?LAOS\s*\(?\s*(\d+\.?\d*),\s*(\d+\.?\d*) rads\)", name)
        if m is None:
            continue
        rows = _sheet_rows(wb, name)
        header, body = rows[0], rows[3:]
        t = _numeric(body, 0)
        cols = _donor_columns(header)
        donors = "".join(sorted(d for d, i in cols.items()
                                if np.size(_numeric(body, i)) > MIN_ROWS_FOR_A_SERIES))
        dt = float(np.median(np.diff(t)))
        period = 2.0 * np.pi / float(m.group(3))
        out.append({"sheet": name, "unidirectional": m.group(1) is not None,
                    "strain_amplitude": float(m.group(2)),
                    "omega_rad_s": float(m.group(3)), "n_samples": int(t.size),
                    "dt_s": dt, "span_s": float(t.max()),
                    "cycles": float(t.max() / period),
                    "samples_per_cycle": float(period / dt),
                    "donors": donors, "n_donors": len(donors)})
    return out


def sweep_census(wb, name: str) -> dict:
    rows = _sheet_rows(wb, name)
    header, body = rows[0], rows[3:]
    x = _numeric(body, 0)
    cols = _donor_columns(header)
    per = {}
    for d, i in sorted(cols.items()):
        v = _numeric(body, i)
        if v.size < MIN_ROWS_SWEEP:
            continue
        per[d] = {"n": int(v.size), "min": float(v.min()), "max": float(v.max()),
                  "n_nonpositive": int((v <= 0).sum())}
    blocks = _blocks(x)
    return {"sheet": name, "x_unit": str(rows[2][0]), "n_x": int(x.size),
            "n_repeats": len(blocks),
            "points_per_repeat": [int(b.stop - b.start) for b in blocks],
            "x_min": float(x.min()), "x_max": float(x.max()),
            "per_donor": per, "n_donors": len(per),
            "donors": "".join(sorted(per)),
            "n_donors_all_positive": sum(1 for v in per.values()
                                         if v["n_nonpositive"] == 0)}


def lve_plateau(wb) -> dict:
    """C3, measured rather than argued: is the frequency sweep's operating strain on a
    PLATEAU in the amplitude sweep? A linear region is flat and extends to lower strain.

    Done PER RAMP, because the sheet is two concatenated ramps over the same axis (see
    `_blocks`): within each ramp, is G' at the lowest strain that ramp's own maximum, or
    is it still rising? The two ramps are NOT assumed to be repeats -- their spread at the
    same strain is measured and reported, and it is far too large for that word."""
    rows = _sheet_rows(wb, "Amplitude Sweep")
    header, body = rows[0], rows[3:]
    per = {}
    for d, i in sorted(_donor_columns(header).items()):
        x, v = _paired(body, i)
        if x.size < MIN_ROWS_SWEEP:
            continue
        blocks = _blocks(x)
        reps = []
        for b in blocks:
            xx, vv = x[b], v[b]
            j = int(np.argmax(vv))
            k = int(np.argmin(np.abs(xx - 1000.0)))
            reps.append({"n": int(xx.size), "strain_min_pct": float(xx[0]),
                         "G_at_min_strain": float(vv[0]), "G_max": float(vv[j]),
                         "strain_at_G_max_pct": float(xx[j]),
                         "first_point_is_the_max": bool(j == 0),
                         "ratio_at_first_point": float(vv[0] / vv[j]),
                         "strain_near_1000_pct": float(xx[k]),
                         "ratio_near_1000_pct": float(vv[k] / vv[j])})
        # the two repeats at their common lowest strain
        spread = (abs(reps[1]["G_at_min_strain"] / reps[0]["G_at_min_strain"] - 1.0)
                  if len(reps) > 1 and reps[0]["G_at_min_strain"] else float("nan"))
        per[d] = {"n_repeats": len(reps), "repeats": reps,
                  "repeat_spread_at_min_strain": float(spread)}
    n_rising = sum(1 for v in per.values()
                   if all(not r["first_point_is_the_max"] for r in v["repeats"]))
    ratios = [r["ratio_near_1000_pct"] for v in per.values() for r in v["repeats"]]
    spreads = [v["repeat_spread_at_min_strain"] for v in per.values()
               if np.isfinite(v["repeat_spread_at_min_strain"])]
    return {"per_donor": per, "n_donors": len(per),
            "n_rising_at_first_point_in_every_repeat": n_rising,
            "strain_min_pct": float(next(iter(per.values()))["repeats"][0]["strain_min_pct"]),
            "median_ratio_near_1000_pct": float(np.median(ratios)),
            "max_ratio_near_1000_pct": float(max(ratios)),
            "median_repeat_spread_at_min_strain": float(np.median(spreads)) if spreads else float("nan"),
            "max_repeat_spread_at_min_strain": float(max(spreads)) if spreads else float("nan")}


def moduli_channels(wb) -> dict:
    """C4: which oscillatory QUANTITIES the workbook stores. A lag and a second timescale
    are separated by PHASE, so the presence or absence of a loss modulus / phase angle
    decides whether the linear channel can discriminate at all."""
    found: dict[str, list[str]] = {}
    for name in wb.sheetnames:
        if name == "Read Me":
            continue
        rows = _sheet_rows(wb, name)
        if len(rows) < 2:
            continue
        labels = {str(c).strip() for c in rows[1] if isinstance(c, str) and c.strip()}
        for label in labels:
            found.setdefault(label, []).append(name)
    phase_like = {k: v for k, v in found.items()
                  if re.search(r"loss|phase|delta|tan|G''|viscous", k, re.I)}
    return {"labels": {k: len(v) for k, v in sorted(found.items())},
            "phase_like": {k: len(v) for k, v in sorted(phase_like.items())}}


def main() -> None:
    if not WORKBOOK.exists():
        raise SystemExit(f"missing {WORKBOOK}; see data/horner2020/SOURCE.md")
    wb = openpyxl.load_workbook(WORKBOOK, read_only=True)

    laos = laos_census(wb)
    freq = sweep_census(wb, "Frequency Sweep")
    amp = sweep_census(wb, "Amplitude Sweep")
    chan = moduli_channels(wb)
    lve = lve_plateau(wb)

    laos_union = set().union(*[set(r["donors"]) for r in laos])
    laos_every = set(laos[0]["donors"])
    for r in laos:
        laos_every &= set(r["donors"])

    # ---------------------------------------------------------------- the five criteria
    c1_n = len(laos_union)
    c1 = {"name": "C1 population", "bar": f">= {COVERAGE_BAR}/20 donors",
          "measured": f"{c1_n}/20 on at least one LAOS sheet, "
                      f"{len(laos_every)}/20 on every LAOS sheet, "
                      f"{freq['n_donors']}/20 on the Frequency Sweep",
          "pass": c1_n >= COVERAGE_BAR}
    c2 = {"name": "C2 linear channel is a measurement",
          "bar": f"G' > 0 at every frequency on >= {COVERAGE_BAR}/20 donors",
          "measured": f"{freq['n_donors_all_positive']}/20 donors have G' positive "
                      f"everywhere ({freq['n_donors']}/20 have data at all; "
                      f"{freq['n_donors'] - freq['n_donors_all_positive']} of those carry "
                      "a non-positive G')",
          "pass": freq["n_donors_all_positive"] >= COVERAGE_BAR}
    c3_readings = {}
    for label, strain in FREQ_SWEEP_STRAIN_READINGS.items():
        need = LVE_HEADROOM * strain
        c3_readings[label] = {"freq_sweep_strain_pct": strain,
                              "needed_amplitude_sweep_min_pct": need,
                              "actual_amplitude_sweep_min_pct": amp["x_min"],
                              "pass": amp["x_min"] <= need}
    c3 = {"name": "C3 linearity demonstrated",
          "bar": f"Amplitude Sweep reaches <= {LVE_HEADROOM:g}x the frequency sweep's strain",
          "measured": "; ".join(
              f"under '{k}': needs <= {v['needed_amplitude_sweep_min_pct']:g} %, "
              f"sweep stops at {v['actual_amplitude_sweep_min_pct']:.4g} % "
              f"-> {'PASS' if v['pass'] else 'FAILED'}" for k, v in c3_readings.items()),
          "pass": all(v["pass"] for v in c3_readings.values()),
          "readings": c3_readings,
          "plateau": (f"and measured directly, per repeat: G' is still RISING at the "
                      f"sweep's lowest strain ({lve['strain_min_pct']:.4g} %) in every "
                      f"ramp on {lve['n_rising_at_first_point_in_every_repeat']}/"
                      f"{lve['n_donors']} donors, so that point is on a flank and not a "
                      f"plateau; by ~1000 % strain G' has fallen to a median "
                      f"{lve['median_ratio_near_1000_pct']:.3f} (worst "
                      f"{lve['max_ratio_near_1000_pct']:.3f}) of its own maximum")}
    c4 = {"name": "C4 discriminating channel exists",
          "bar": "a loss modulus, phase angle, or strain-aligned stress waveform",
          "measured": ("phase-like column labels found: "
                       + (", ".join(f"{k} ({n} sheets)" for k, n in chan["phase_like"].items())
                          if chan["phase_like"] else "NONE")
                       + "; the workbook's oscillatory column labels are "
                       + ", ".join(sorted(chan["labels"]))),
          "pass": bool(chan["phase_like"])}
    spc_min = min(r["samples_per_cycle"] for r in laos)
    spc_max = max(r["samples_per_cycle"] for r in laos)
    c5 = {"name": "C5 sampling (item 7's stated question)",
          "bar": f">= {SAMPLES_PER_CYCLE_BAR} samples per cycle",
          "measured": f"{spc_min:.0f} to {spc_max:.0f} samples per cycle over "
                      f"{len(laos)} sheets, {min(r['cycles'] for r in laos):.1f} to "
                      f"{max(r['cycles'] for r in laos):.1f} cycles each",
          "pass": spc_min >= SAMPLES_PER_CYCLE_BAR}
    criteria = [c1, c2, c3, c4, c5]

    linear_route_open = c2["pass"] and c3["pass"] and c4["pass"]
    verdict = "OPEN" if (linear_route_open and c1["pass"]) else "CLOSED"

    # ------------------------------------------------------------------------- the page
    L = [f"# Can the oscillatory sheets say whether the blood step rise is a lag or a "
         f"second timescale? **{verdict}.**", "",
         "Generated by `scripts/oscillatory_capability_horner.py`; every number below is "
         "printed by it. NEXT_STEPS item 7. The five criteria and the verdict rule were "
         "written into that script's docstring before any of them was graded.", "",
         "## The five pre-registered criteria", "",
         "| criterion | bar | measured | verdict |", "|---|---|---|---|"]
    for c in criteria:
        L.append(f"| {c['name']} | {c['bar']} | {c['measured']} | "
                 f"{'PASS' if c['pass'] else '**FAILED**'} |")
    L += ["",
          f"Verdict rule, fixed in advance: the linear-response route is closed if C2, C3 "
          f"or C4 fails. C2 {'passed' if c2['pass'] else 'FAILED'}, "
          f"C3 {'passed' if c3['pass'] else 'FAILED'}, "
          f"C4 {'passed' if c4['pass'] else 'FAILED'}. "
          f"The route is **{'open' if linear_route_open else 'closed'}**.", "",
          "## What actually decides it, and what does not", "",
          "**Item 7's own question is not the binding one.** It asked whether the LAOS "
          f"sheets have enough samples per cycle; they have {spc_min:.0f} to {spc_max:.0f}, "
          f"over {min(r['cycles'] for r in laos):.1f} to {max(r['cycles'] for r in laos):.1f} "
          "cycles, against a 4-6 sample rise in the step sheet. Sampling is not the "
          "constraint and the item cannot be closed on that ground. Three other things "
          "close it.", "",
          "**C4: the workbook stores no phase.** "
          + ("A phase-like channel was found, so this paragraph is stale and the script "
             "should be re-read." if chan["phase_like"] else
             "The Frequency Sweep and Amplitude Sweep sheets carry a **storage modulus "
             "only** -- there is no loss modulus, no phase angle, and no tan(delta) "
             "anywhere in the file. A lag and a second timescale are separated chiefly by "
             "PHASE: a dead time theta contributes a phase -w*theta that grows without "
             "bound and leaves the modulus alone, while a second Debye process contributes "
             "a bounded phase excursion with a matching dispersion in the modulus. With "
             "only |G'| the two candidates are being compared on the channel where they "
             "agree most."), "",
          f"**C3: linearity is never demonstrated.** The Read Me says the frequency sweep "
          f"ran at \"a strain amplitude of 10\" and the Amplitude Sweep's own axis starts "
          f"at {amp['x_min']:.4g} % and runs to {amp['x_max']:.4g} %. The unit of that 10 "
          "is ambiguous in the source -- the Amplitude Sweep is labelled in %, while the "
          "Read Me says the LAOS sheet labels are dimensionless -- so both readings are "
          "carried here and **neither is guessed**:", ""]
    for label, v in c3_readings.items():
        L.append(f"- **{label}**: an LVE plateau containing that strain would have to be "
                 f"shown from below, i.e. the amplitude sweep would have to reach "
                 f"{v['needed_amplitude_sweep_min_pct']:g} % or lower. It stops at "
                 f"{v['actual_amplitude_sweep_min_pct']:.4g} %. "
                 f"{'PASS' if v['pass'] else '**FAILED**'}.")
    L += ["",
          "**Neither reading gives a linear measurement, and that is measured, not "
          "argued.** Under the first, the amplitude sweep's very first point IS the "
          f"frequency sweep's operating point, so nothing below it was measured -- and G' "
          f"is still RISING there, in both ramps, on "
          f"{lve['n_rising_at_first_point_in_every_repeat']} of {lve['n_donors']} donors, "
          "so that point is on a flank rather than on a plateau. Under the second, the "
          "frequency sweep "
          f"sits at ~1000 % strain, where G' has fallen to a median "
          f"{lve['median_ratio_near_1000_pct']:.3f} of its own maximum (worst case "
          f"{lve['max_ratio_near_1000_pct']:.3f}) -- roughly two decades in STRAIN past "
          "where G' peaks, with the modulus down by a factor of about five, so it is "
          "nonlinear by the source's own data. **The conclusion does not depend on which "
          "reading is right**, which is why the ambiguity is reported rather than "
          "resolved. Only the C3 pass/fail arithmetic differs between them, and the rule "
          "fixed in advance requires every reading to pass.", "",
          f"**C1: the population is not item 2's population.** The step-rise result is on "
          f"all 20 donors. The Frequency Sweep has {freq['n_donors']}/20, the LAOS sheets "
          f"have {c1_n}/20 on at least one sheet and only {len(laos_every)}/20 on every "
          f"sheet. Donors **{''.join(sorted(set(HUMAN_DONORS) - laos_union))}** have no "
          "LAOS data at all. An oscillatory answer would speak for a different set of "
          "people than the finding it is meant to explain.", "",
          "## What a discriminating measurement would need", "",
          "Not in this workbook, and worth writing down rather than working around:", "",
          "1. A **phase or loss channel** on the same donors -- G''(w) or tan(delta), or "
          "the raw stress waveform aligned to the strain so one can be computed.",
          "2. An **amplitude sweep that reaches into the linear region**, decades below "
          "the strain the frequency sweep uses, so that the linear signatures apply at all.",
          "3. Either that, or a **step experiment sampled fast enough on the rise**. Item "
          "2's rise carries 4-6 samples; a lag of order 0.1-1 s against a tau of 1-2.7 s "
          "needs roughly an order of magnitude more, which is a request to the instrument, "
          "not to the analysis.", "",
          "**Both sweep sheets are two concatenated ramps over the same axis**, not one "
          f"sweep: the Amplitude Sweep is {amp['points_per_repeat']} points from "
          f"{amp['x_min']:.4g} {amp['x_unit']} and the Frequency Sweep is "
          f"{freq['points_per_repeat']} from {freq['x_min']:g} {freq['x_unit']}. That "
          "matters mechanically -- a maximum taken over the whole column lands in "
          "whichever ramp happens to be larger and mixes two runs, so every number above "
          "is computed per ramp -- and it raises a question the workbook does not answer. "
          "At the lowest strain the two ramps of the amplitude sweep differ by a median "
          f"{lve['median_repeat_spread_at_min_strain']:.0%} (worst "
          f"{lve['max_repeat_spread_at_min_strain']:.0%}). **That is far too large to call "
          "them repeats**; they are more likely two conditions (the Read Me says the "
          "amplitude sweep shown is 'representative' and was run at 2 Hz, without saying "
          "what the second ramp is). Either way it does not move the C3 verdict, whose "
          f"decisive fact needs no ramp reasoning at all: nothing below {amp['x_min']:.4g} % "
          "strain was measured, on either ramp.", "",
          "## The LAOS census, for the record", "",
          "| sheet | strain | w (rad/s) | n | cycles | samples/cycle | donors |",
          "|---|---|---|---|---|---|---|"]
    for r in laos:
        L.append(f"| {r['sheet']} | {r['strain_amplitude']:g} | {r['omega_rad_s']:g} | "
                 f"{r['n_samples']} | {r['cycles']:.1f} | {r['samples_per_cycle']:.0f} | "
                 f"{r['n_donors']} ({r['donors']}) |")
    L += ["",
          f"Frequency Sweep: {freq['n_x']} points, {freq['x_min']:g} to {freq['x_max']:g} "
          f"{freq['x_unit']}, {freq['n_donors']}/20 human donors ({freq['donors']}). "
          f"Amplitude Sweep: {amp['n_x']} points, {amp['x_min']:.4g} to {amp['x_max']:.4g} "
          f"{amp['x_unit']}, {amp['n_donors']}/20 ({amp['donors']}).", ""]

    (DATA / "oscillatory_capability.md").write_text("\n".join(L), encoding="utf-8",
                                                    newline="\n")
    payload = {"verdict": verdict, "linear_route_open": linear_route_open,
               "criteria": criteria, "laos": laos, "frequency_sweep": freq,
               "amplitude_sweep": amp, "channels": chan, "lve_plateau": lve,
               "laos_donors_union": "".join(sorted(laos_union)),
               "laos_donors_every_sheet": "".join(sorted(laos_every)),
               "donors_with_no_laos": "".join(sorted(set(HUMAN_DONORS) - laos_union))}
    (DATA / "oscillatory_capability.json").write_text(json.dumps(payload, indent=2),
                                                      encoding="ascii")
    for c in criteria:
        print(f"{'PASS' if c['pass'] else 'FAILED'}  {c['name']}: {c['measured']}")
    print(f"\nlinear route {'OPEN' if linear_route_open else 'CLOSED'}; item 7 {verdict}")


if __name__ == "__main__":
    main()
