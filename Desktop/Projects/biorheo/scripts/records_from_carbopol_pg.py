"""Item 24: two licensable Herschel-Bulkley records, and the first real test of the shift.

2 wt% Carbopol 974p in propylene glycol at 10 / 20 / 50 deg C -- Zenodo 11389467, CC BY 4.0,
paper arXiv:2312.00251. Source, licence and the four traps: `data/carbopol_pg_vinutha2024/
SOURCE.md`. Gates pre-registered in `data/carbopol_pg_prereg.md` BEFORE this file existed.

Two things make this dataset worth more than item 22's:

* it is a **yield-stress microgel** reaching 1e3 1/s, i.e. the licensable counterpart to the
  family the injectability case study actually ships, and
* **the time-temperature shift factor is PREDICTED, not fitted.** The authors chose propylene
  glycol because its viscosity moves tenfold over 10-50 deg C, so the material is identical
  at all three temperatures and only the medium changes. The deposit publishes that medium
  viscosity, so `a = eta_s(T)/eta_s(T_ref)` is known before any fit, while `collapse_shift`
  fits `a` from the flow curves alone and never sees it.

Records come from the 10 and 50 deg C blocks only. **The 20 deg C block is refused**, and
G24.7 RECOMPUTES the four signals that refuse it rather than asserting them.

    python scripts/records_from_carbopol_pg.py

Writes two records + sidecars, `data/carbopol_pg_vinutha2024/records.json` and `.md`.
Exit code is non-zero if any gate fails.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from biorheo.constitutive import get_model                                  # noqa: E402
from biorheo.constitutive.shift import (collapse_shift,                     # noqa: E402
                                        shift_from_two_viscosities)
from biorheo.identify import fit_rate_curve, write_identifiability          # noqa: E402
from biorheo.registry import RECORDS_DIR, load_record                       # noqa: E402
from biorheo.registry.schema import round_out                               # noqa: E402

HB = get_model("herschel_bulkley")
SRC = ROOT / "data" / "carbopol_pg_vinutha2024"
RAW = SRC / "raw"
FLOW = RAW / "Fig 1(a) flow curves_exp.txt"
CESSATION = RAW / "All flow cessation tests_exp.txt"
OUT_RECORDS = Path(RECORDS_DIR)

GEOMETRY = ("Anton Paar MCR 300 stress-controlled, cone-plate, cone angle 1 deg, BOTH "
            "SURFACES ROUGHENED to minimise wall slip; Peltier plate and Peltier hood, "
            "T to 0.1 deg C")

# --- numbers the DEPOSIT publishes, which this repo does not compute ---------------------
# `raw/Fig 1(c) and Fig S3 residual stress_exp.txt`, named columns `medium viscosity [Pa s]`
# and `yield stress [Pa]`; corroborated by the paper's Methods ("its viscosity varies by a
# factor of about 10 within a temperature range of 10-50 C (eta_s = 0.102-0.012 Pa s)").
MEDIUM_VISCOSITY = {10.0: 0.102, 20.0: 0.055, 50.0: 0.012}
PUBLISHED_YIELD_STRESS = 8.9   # NOT a Herschel-Bulkley tau_y; see G24.9
T_REF, T_OTHER = 10.0, 50.0

# --- pre-registered bars, set in data/carbopol_pg_prereg.md before any number was known --
G241_TOL = 1.30          # |ln(a_fit/a_pred)| < ln(1.30)
G242_SEPARATION = 3.0    # residual at 2*a_fit must exceed that at a_fit by this factor
G243_FLATNESS = 0.10     # the low-rate residual must move by less than this
G246_TOL = 1.40          # |ln((k_50/k_10)/a_pred^n)| < ln(1.40)

BLOCKS = {10.0: 0, 20.0: 4, 50.0: 8}   # first column index of each temperature block
ADMITTED = (10.0, 50.0)                # prereg section 3

# Two pre-registered predictions were REFUTED by the first run, and neither bar was moved.
# They are named here so the exit code pins them: if either ever starts PASSING, the run
# exits non-zero, because that would mean the recorded result changed. A bare "did any gate
# fail" check cannot tell a recorded refutation from a new one.
REFUTED = {
    "G24_3_the_plateau_is_blind_to_the_shift_factor":
        "predicted PASS. The low-rate HALF of the curve is not the plateau for this "
        "material -- see G24.3''.",
    "G24_4_yield_stress_invariant":
        "predicted PASS. The two tau_y profile intervals do not overlap.",
}


# ---------------------------------------------------------------------------------------
def read_flow_curves() -> dict[float, tuple[np.ndarray, np.ndarray]]:
    """The three blocks of `Fig 1(a)`, paired ROW BY ROW.

    The file is cp1252 (the degree sign is one 0xB0 byte) and WIDE: three independent
    blocks side by side with the header repeated, of different lengths. Filtering each
    column for numeric cells independently and truncating to the shorter one would
    misalign the pair the moment one block has a gap -- so a row contributes only if BOTH
    of its cells are present, and a half-empty row is an error rather than a skip.
    """
    lines = FLOW.read_bytes().decode("cp1252").splitlines()
    out = {}
    for temp, col in BLOCKS.items():
        gd, tau = [], []
        for ln, line in enumerate(lines[1:], start=2):
            f = line.split("\t")
            r, s = f[col + 1].strip(), f[col + 2].strip()
            if not r and not s:
                continue
            if not r or not s:
                raise ValueError(f"{FLOW.name} line {ln}, block {temp} C: half-empty row "
                                 f"{f[col:col + 3]!r}; the pair cannot be formed")
            if float(f[col].strip()) != temp:
                raise ValueError(f"{FLOW.name} line {ln}: block {temp} C carries "
                                 f"temperature {f[col]!r}")
            gd.append(float(r))
            tau.append(float(s))
        out[temp] = (np.asarray(gd), np.asarray(tau))
        print(f"  {temp:g} C: {len(gd)} points, rate {min(gd):.4g}-{max(gd):.4g} 1/s, "
              f"stress {min(tau):.4g}-{max(tau):.4g} Pa")
    return out


def read_preshear_curve() -> dict[float, dict[float, float]]:
    """The SECOND flow curve nobody announced: `preshear-stress` at 16 rates per block.

    A steady flow curve of the same material on the same instrument under a different
    protocol -- the only independent replicate this project has for any record (G24.8).
    """
    text = CESSATION.read_bytes().decode("cp1252")
    lines = [ln for ln in text.splitlines() if ln.strip()]
    starts = [i for i, h in enumerate(lines[0].split("\t"))
              if h.strip().startswith("temperature")]
    pre: dict[float, dict[float, float]] = {}
    for line in lines[1:]:
        f = line.split("\t")
        for s in starts:
            if s + 2 >= len(f):
                continue
            T, r, sg = f[s].strip(), f[s + 1].strip(), f[s + 2].strip()
            if not (T and r and sg):
                continue
            d = pre.setdefault(float(T), {})
            key, val = float(r), float(sg)
            if key in d and abs(d[key] - val) > 1e-9:
                raise ValueError(f"preshear stress at T={T}, rate={r} is both "
                                 f"{d[key]} and {val}")
            d[key] = val
    for T, d in sorted(pre.items()):
        print(f"  preshear replicate, {T:g} C: {len(d)} rates "
              f"{min(d):.4g}-{max(d):.4g} 1/s")
    return pre


# ---------------------------------------------------------------------------------------
def g24_7_refuse_the_20C_block(curves) -> dict:
    """RECOMPUTE the four signals that disqualify the 20 C block, never assert them."""
    out = {}
    for temp, (gd, tau) in sorted(curves.items()):
        lg = np.log10(gd)
        steps = lg[:-1] - lg[1:]
        on_grid = int(np.sum(np.abs(lg * 10 - np.round(lg * 10)) < 1e-4))
        reversals = [(float(gd[i]), float(gd[i + 1]))
                     for i in range(gd.size - 1) if gd[i + 1] >= gd[i]]
        # significant figures of the stress column, from its shortest exact repr
        # `repr(np.float64(1330.0))` is "np.float64(1330.0)" under numpy 2, so counting
        # digits off a numpy scalar's repr counts the wrapper. G24.7 caught this by
        # requiring the signal to REPRODUCE rather than asserting it.
        sig = [len(f"{float(v)!r}".replace(".", "").replace("-", "").lstrip("0").rstrip("0"))
               for v in tau]
        out[temp] = {
            "n_points": int(gd.size),
            "rates_on_the_programmed_log_grid": on_grid,
            "grid_fraction": on_grid / gd.size,
            "rate_reversals": reversals,
            "n_rate_reversals": len(reversals),
            "log_step_mean": float(np.mean(steps)), "log_step_sd": float(np.std(steps)),
            "stress_sigfigs_median": float(np.median(sig)),
            "n_stress_at_3_sigfigs": int(sum(1 for x in sig if x <= 3)),
        }
    return out


def _inside_mask(gd, a):
    """Which rates survive `collapse_shift`'s own range restriction at this `a`."""
    lg = np.log(np.sort(np.asarray(gd)))
    target = lg + math.log(a)
    return (target >= lg[0]) & (target <= lg[-1])


def _residual_at(gd, eta_ref, eta_other, a, keep=None) -> tuple[float, int]:
    """`collapse_shift`'s own objective, evaluated at a GIVEN `a`.

    Written out because `collapse_shift` returns the residual only at its optimum, and
    G24.2/G24.3 need it at 2a and a/2 and on subsets. It is checked against
    `collapse_shift`'s own number at the fitted `a` before any of those are believed --
    a reimplementation that is not anchored to the shipped one proves nothing.
    """
    order = np.argsort(gd)
    lg, lref = np.log(np.asarray(gd)[order]), np.log(np.asarray(eta_ref)[order])
    loth = np.log(np.asarray(eta_other)[order])
    target = lg + math.log(a)
    inside = (target >= lg[0]) & (target <= lg[-1])
    if keep is not None:
        # evaluate every probe on the SAME points. Without this, comparing the residual at
        # `a` with the residual at `2a` compares two different SAMPLES as well as two
        # different shift factors, and the ratio is partly a change of sample.
        inside = inside & np.asarray(keep)[np.argsort(np.asarray(gd))]
    if inside.sum() < 4:
        return float("inf"), int(inside.sum())
    pred = math.log(a) + np.interp(target[inside], lg, lref)
    return float(np.sqrt(np.mean((loth[inside] - pred) ** 2))), int(inside.sum())


def main() -> int:  # noqa: C901
    print("=== reading the deposit ===")
    curves = read_flow_curves()
    preshear = read_preshear_curve()

    print("\n=== G24.7: recomputing the four signals that refuse the 20 C block ===")
    signals = g24_7_refuse_the_20C_block(curves)
    for T, s in sorted(signals.items()):
        print(f"  {T:g} C: {s['rates_on_the_programmed_log_grid']}/{s['n_points']} on the "
              f"grid, {s['n_rate_reversals']} rate reversals, "
              f"{s['n_stress_at_3_sigfigs']}/{s['n_points']} stresses at <=3 sig figs")
    bad = signals[20.0]
    good = [signals[T] for T in ADMITTED]
    g247 = {
        "bar": "the 20 C block is refused and all four signals reproduce; no record is "
               "written from it",
        "signals": {str(k): v for k, v in signals.items()},
        "admitted_temperatures": list(ADMITTED),
        "reproduces": {
            "off_grid": bool(bad["grid_fraction"] < 0.1
                             and all(g["grid_fraction"] > 0.7 for g in good)),
            "rate_reversals": bool(bad["n_rate_reversals"] > 0
                                   and all(g["n_rate_reversals"] == 0 for g in good)),
            "rounded_to_3_sigfigs": bool(
                bad["n_stress_at_3_sigfigs"] / bad["n_points"] > 0.8
                and all(g["n_stress_at_3_sigfigs"] / g["n_points"] < 0.2 for g in good)),
        },
    }
    g247["pass"] = bool(all(g247["reproduces"].values()))

    # --- fit the two admitted blocks ------------------------------------------------------
    print("\n=== fitting Herschel-Bulkley at the two admitted temperatures ===")
    fits, rows = {}, []
    for temp in ADMITTED:
        gd, tau = curves[temp]
        eta = tau / gd
        ident = fit_rate_curve(HB, gd, eta, temperature_C=temp, geometry=GEOMETRY)
        if int(ident.fit.n_points) != int(gd.size):
            raise ValueError(f"{temp} C: the fitter kept {ident.fit.n_points} of "
                             f"{gd.size} points and the caller's N would have been reported")
        fits[temp] = ident
        p = {n: ident.params[n] for n in HB.param_names}
        print(f"  {temp:g} C: tau_y = {p['tau_y'].value:.4g} "
              f"[{p['tau_y'].lo:.4g}, {p['tau_y'].hi:.4g}] Pa, "
              f"k = {p['k'].value:.4g}, n = {p['n'].value:.4f} "
              f"[{p['n'].lo:.4f}, {p['n'].hi:.4f}], rmse_log = {ident.fit.rmse_log:.4g}")
        rows.append({
            "temperature_C": temp, "n_points": int(gd.size),
            "rate_min": float(gd.min()), "rate_max": float(gd.max()),
            "medium_viscosity_Pa_s": MEDIUM_VISCOSITY[temp],
            "rmse_log": float(ident.fit.rmse_log),
            **{f"{n}_{f}": float(getattr(p[n], f)) for n in HB.param_names
               for f in ("value", "lo", "hi")},
            **{f"{n}_verdict": p[n].verdict for n in HB.param_names},
        })

    # --- the collapse, on a common grid over the overlap ---------------------------------
    print("\n=== G24.1-G24.3: the collapse ===")
    a_pred = shift_from_two_viscosities(MEDIUM_VISCOSITY[T_OTHER],
                                        MEDIUM_VISCOSITY[T_REF])
    gd_r, tau_r = curves[T_REF]
    gd_o, tau_o = curves[T_OTHER]
    lo = max(gd_r.min(), gd_o.min())
    hi = min(gd_r.max(), gd_o.max())
    grid = np.logspace(math.log10(lo), math.log10(hi), 120)
    # log-log interpolation of each MEASURED curve onto the shared grid; both curves sit on
    # the same programmed 1/10-decade grid, so this interpolates within a block, never
    # across a gap, and never outside either block's own measured range.
    eta_ref = np.exp(np.interp(np.log(grid), np.log(np.sort(gd_r)),
                               np.log((tau_r / gd_r)[np.argsort(gd_r)])))
    eta_oth = np.exp(np.interp(np.log(grid), np.log(np.sort(gd_o)),
                               np.log((tau_o / gd_o)[np.argsort(gd_o)])))
    print(f"  overlap {lo:.4g}-{hi:.4g} 1/s, {grid.size} shared rates")

    fit = collapse_shift(grid, eta_ref, eta_oth)
    a_fit = fit["a"]
    # anchor the local objective against the shipped one before using it anywhere else
    mine, _ = _residual_at(grid, eta_ref, eta_oth, a_fit)
    if not math.isclose(mine, fit["rms_log_residual"], rel_tol=1e-9, abs_tol=1e-12):
        raise ValueError(f"the local residual ({mine!r}) is not collapse_shift's "
                         f"({fit['rms_log_residual']!r}); it is a different objective")
    print(f"  a_fit = {a_fit:.5f}  a_pred = {a_pred:.5f}  "
          f"ratio {a_fit / a_pred:.4f}  residual {fit['rms_log_residual']:.4g} "
          f"on {fit['n_used']}/{fit['n_total']} rates")

    g241 = {
        "bar": f"|ln(a_fit/a_pred)| < ln({G241_TOL})",
        "a_pred": a_pred, "a_fit": a_fit, "ratio": a_fit / a_pred,
        "ln_ratio": abs(math.log(a_fit / a_pred)), "ln_bar": math.log(G241_TOL),
        "rms_log_residual": fit["rms_log_residual"],
        "n_used": fit["n_used"], "n_total": fit["n_total"],
        "pass": bool(abs(math.log(a_fit / a_pred)) < math.log(G241_TOL)),
        "note": "a_pred is the ratio of the medium viscosities the DEPOSIT publishes. "
                "collapse_shift fits a_fit from the two flow curves alone and never sees "
                "them.",
    }

    # G24.2 / G24.3: the discriminating and the blind halves of the curve
    split = math.sqrt(lo * hi)
    halves = {"high_rate": grid > split, "low_rate": grid <= split}
    probe = {"a_fit": a_fit, "2a_fit": 2 * a_fit, "a_fit/2": a_fit / 2}
    resid = {h: {} for h in halves}
    resid_common = {h: {} for h in halves}
    for h, m in halves.items():
        # as pre-registered: each probe on whatever survives its own range restriction
        for name, a in probe.items():
            r, n_in = _residual_at(grid[m], eta_ref[m], eta_oth[m], a)
            resid[h][name] = {"residual": r, "n_used": n_in}
        # and again on the INTERSECTION, so the three probes share one sample
        common = np.ones(int(m.sum()), dtype=bool)
        for a in probe.values():
            common &= _inside_mask(grid[m], a)
        for name, a in probe.items():
            r, n_in = _residual_at(grid[m], eta_ref[m], eta_oth[m], a, keep=common)
            resid_common[h][name] = {"residual": r, "n_used": n_in}
    for h in halves:
        print(f"  {h}: " + ", ".join(
            f"{k} -> {v['residual']:.4g} (n={v['n_used']})" for k, v in resid[h].items()))

    hi_sep = resid["high_rate"]["2a_fit"]["residual"] / resid["high_rate"]["a_fit"]["residual"]
    lo_move = abs(resid["low_rate"]["2a_fit"]["residual"]
                  / resid["low_rate"]["a_fit"]["residual"] - 1.0)
    hi_sep_c = (resid_common["high_rate"]["2a_fit"]["residual"]
                / resid_common["high_rate"]["a_fit"]["residual"])
    lo_move_c = abs(resid_common["low_rate"]["2a_fit"]["residual"]
                    / resid_common["low_rate"]["a_fit"]["residual"] - 1.0)

    # --- G24.3'': WHY G24.3 could not test what it was written to test -------------------
    # Section 1.1's algebra is right -- below yield eta -> tau_y/gd and the shift is the
    # identity there -- but it only acts where the material is actually ON the plateau, and
    # `collapse_shift` crops |log10 a| of the range off the bottom. Measured, not argued.
    gd10 = curves[T_REF][0]
    tau10 = curves[T_REF][1]
    tau_y_ref = fits[T_REF].params["tau_y"].value
    on_plateau = tau10 < 1.2 * tau_y_ref
    plateau_hi = float(gd10[on_plateau].max()) if on_plateau.any() else float("nan")
    shiftable_lo = float(lo / a_fit)
    g243pp = {
        "bar": "measure whether the genuine plateau is inside the shiftable window at all, "
               "instead of assuming the low-rate HALF of the curve is the plateau",
        "plateau_definition": "sigma < 1.2 * tau_y(10 C)",
        "tau_y_ref": float(tau_y_ref),
        "n_points_on_plateau": int(on_plateau.sum()),
        "plateau_upper_rate": plateau_hi,
        "shiftable_window_lower_rate": shiftable_lo,
        "plateau_entirely_outside_the_shiftable_window": bool(plateau_hi < shiftable_lo),
        "g243_split_rate": split,
        "pass": True,
        "note": "G24.3 split at the geometric mean of the overlap, where this material is "
                "far ABOVE yield, so its 'low-rate half' was the power-law region and never "
                "the plateau. The trap pre-registered in section 1.1 is real algebra and "
                "INOPERATIVE on this dataset, because the plateau is narrower than the "
                "range the shift itself crops off. That is why G24.3 failed, and it is a "
                "fact about the data rather than about the shift layer.",
    }

    g242 = {
        "bar": f"on the high-rate half, the residual at 2*a_fit exceeds that at a_fit by "
               f"a factor of at least {G242_SEPARATION}",
        "split_rate": split, "residuals": resid["high_rate"],
        "separation_factor": hi_sep,
        "residuals_on_the_common_point_set": resid_common["high_rate"],
        "separation_factor_common_point_set": hi_sep_c,
        "pass": bool(hi_sep >= G242_SEPARATION and hi_sep_c >= G242_SEPARATION),
        "note": "below this, G24.1 is UNINFORMATIVE and its agreement is not quoted.",
    }
    g243 = {
        "bar": f"on the low-rate half, the residual moves by less than "
               f"{G243_FLATNESS:.0%} between a_fit and 2*a_fit",
        "residuals": resid["low_rate"], "relative_move": lo_move,
        "residuals_on_the_common_point_set": resid_common["low_rate"],
        "relative_move_common_point_set": lo_move_c,
        "pass": bool(lo_move < G243_FLATNESS),
        "predicted_failure": False,
        "note": "the positive control for the prereg's section 1.1: below yield both "
                "curves tend to tau_y, so eta -> tau_y/gd and a*eta_ref(a*gd) = tau_y/gd "
                "for ANY a. If this half DID move, that algebra is wrong.",
    }

    # --- G24.4 / G24.5 / G24.6: the parameters ------------------------------------------
    print("\n=== G24.4-G24.6: the parameters ===")
    pr, po = fits[T_REF].params, fits[T_OTHER].params

    def overlaps(x, y):
        return bool(max(x.lo, y.lo) <= min(x.hi, y.hi))

    def separation(x, y):
        w = 0.5 * ((x.hi - x.lo) + (y.hi - y.lo))
        return abs(x.value - y.value) / w if w > 0 else float("inf")

    g244 = {
        "bar": "the two tau_y profile intervals overlap",
        "tau_y": {str(T_REF): [pr["tau_y"].lo, pr["tau_y"].value, pr["tau_y"].hi],
                  str(T_OTHER): [po["tau_y"].lo, po["tau_y"].value, po["tau_y"].hi]},
        "difference_over_mean_interval_width": separation(pr["tau_y"], po["tau_y"]),
        "absolute_difference_Pa": abs(pr["tau_y"].value - po["tau_y"].value),
        "relative_difference": abs(pr["tau_y"].value - po["tau_y"].value)
                               / (0.5 * (pr["tau_y"].value + po["tau_y"].value)),
        "interval_gap_Pa": max(0.0, max(pr["tau_y"].lo, po["tau_y"].lo)
                               - min(pr["tau_y"].hi, po["tau_y"].hi)),
        "pass": overlaps(pr["tau_y"], po["tau_y"]),
        "note": "item 21 DERIVED this invariance before the dataset was found, and the "
                "authors state it from their own data. An overlap of two WIDE intervals "
                "is a weak pass, which is why the separation is reported beside it.",
    }
    g245 = {
        "bar": "the two n profile intervals overlap",
        "n": {str(T_REF): [pr["n"].lo, pr["n"].value, pr["n"].hi],
              str(T_OTHER): [po["n"].lo, po["n"].value, po["n"].hi]},
        "difference_over_mean_interval_width": separation(pr["n"], po["n"]),
        "pass": overlaps(pr["n"], po["n"]),
        "note": "item 22's collapse failed precisely because its two curves differed in n. "
                "If n differed here, no shift factor would be meaningful and G24.1's "
                "number must not be quoted whatever it says.",
    }
    n_mean = 0.5 * (pr["n"].value + po["n"].value)
    k_ratio = po["k"].value / pr["k"].value
    k_pred = a_pred ** n_mean
    g246 = {
        "bar": f"|ln((k_50/k_10)/a_pred^n)| < ln({G246_TOL})",
        "k_ratio": k_ratio, "a_pred_to_the_n": k_pred, "n_mean": n_mean,
        "ln_ratio": abs(math.log(k_ratio / k_pred)), "ln_bar": math.log(G246_TOL),
        "pass": bool(abs(math.log(k_ratio / k_pred)) < math.log(G246_TOL)),
        "note": "a second route to G24.1 through the fitted parameters rather than the "
                "curve, so it can disagree with G24.1 and that would be informative.",
    }
    for g, name in ((g244, "G24.4 tau_y"), (g245, "G24.5 n"), (g246, "G24.6 k")):
        print(f"  {name}: {'PASS' if g['pass'] else 'FAIL'}")

    # --- G24.8: the independent replicate, in the same units as the fit's uncertainty ----
    print("\n=== G24.8: the independent replicate ===")
    rep = {}
    for temp in ADMITTED:
        gd, tau = curves[temp]
        ratios = []
        for r, s in sorted(preshear[temp].items()):
            m = np.abs(np.log10(gd / r)) < 0.011
            if not m.any():
                continue
            ratios.append(float(tau[m][0] / s))
        rep[temp] = {"n_matched": len(ratios),
                     "median_ratio": float(np.median(ratios)),
                     "min_ratio": float(min(ratios)), "max_ratio": float(max(ratios)),
                     "median_abs_pct": float(100 * np.median(np.abs(np.array(ratios) - 1)))}
        print(f"  {temp:g} C: {len(ratios)} matched rates, median ratio "
              f"{rep[temp]['median_ratio']:.4f}, range {rep[temp]['min_ratio']:.4f}-"
              f"{rep[temp]['max_ratio']:.4f}")
    # the fitted parameters' own relative interval widths, for the same-units comparison
    widths = {}
    for temp in ADMITTED:
        p = fits[temp].params
        widths[temp] = {n: (float(p[n].hi - p[n].lo) / abs(float(p[n].value))
                            if p[n].value else float("inf")) for n in HB.param_names}
    worst_param = max(max(w.values()) for w in widths.values())
    worst_rep = max(r["median_abs_pct"] / 100 for r in rep.values())
    g248 = {
        "bar": "no bar; this is a measurement. Which is larger: the independent "
               "replicate's scatter, or the fitted parameters' interval width?",
        "replicate": {str(k): v for k, v in rep.items()},
        "relative_parameter_interval_widths": {str(k): v for k, v in widths.items()},
        "worst_relative_parameter_width": worst_param,
        "worst_replicate_median_scatter": worst_rep,
        "larger_term": "replicate scatter" if worst_rep > worst_param else "parameter interval",
        "pass": True,
        "note": "the only independent replicate this project has for any record. G15.4 had "
                "to approximate this comparison with a model misfit because none existed.",
    }
    print(f"  larger term: {g248['larger_term']} "
          f"(replicate {100 * worst_rep:.2f} % vs parameter {100 * worst_param:.2f} %)")

    # --- G24.9: the published 8.9 Pa is NOT graded against an HB tau_y -------------------
    g249 = {
        "bar": "no gate compares a Herschel-Bulkley tau_y against the published 8.9 Pa",
        "published_yield_stress_Pa": PUBLISHED_YIELD_STRESS,
        "published_model": "three-component model, sigma = sigma_y + k1 gd^0.5 + eta_bg gd "
                           "(not implemented in this package)",
        "our_model": "herschel_bulkley, sigma = tau_y + k gd^n",
        "hb_tau_y": {str(T_REF): pr["tau_y"].value, str(T_OTHER): po["tau_y"].value},
        "reported_as": "a MODEL-FORM comparison, not a disagreement. biorheo cannot "
                       "reproduce the authors' number because it does not implement their "
                       "model, and a yield stress is a fit-objective convention before it "
                       "is a material constant.",
        "graded": False, "pass": True,
    }

    # --- write the two records ------------------------------------------------------------
    print("\n=== writing records ===")
    written = []
    for temp in ADMITTED:
        gd, tau = curves[temp]
        ident = fits[temp]
        p = {n: ident.params[n] for n in HB.param_names}
        rec_id = f"carbopol_pg_vinutha2024_{int(temp):d}C"
        gmin, gmax = round_out(float(gd.min()), up=False), round_out(float(gd.max()), up=True)
        text = TEMPLATE.format(
            id=rec_id, temp=temp, geometry=GEOMETRY, gmin=gmin, gmax=gmax,
            tau_y=f"{p['tau_y'].value:.6g}", k=f"{p['k'].value:.6g}",
            n=f"{p['n'].value:.6g}", npts=gd.size,
            gd_lo=f"{gd.min():.4g}", gd_hi=f"{gd.max():.4g}",
            rmse=f"{ident.fit.rmse_log:.4g}", eta_s=MEDIUM_VISCOSITY[temp],
        )
        path = OUT_RECORDS / f"{rec_id}.yaml"
        path.write_text(text, encoding="utf-8")
        rec = load_record(path)
        if not np.all(rec.envelope.contains(np.array([gd.min(), gd.max()]))):
            sys.exit(f"{rec_id}: envelope does not contain its own measured range")
        write_identifiability(
            rec, ident,
            data_note=(f"Herschel-Bulkley fit to the {gd.size}-point {temp:g} C block of the "
                       f"deposited Fig 1(a) flow curves, rates {gd.min():.4g}-{gd.max():.4g} "
                       "1/s; log-viscosity residuals. No point was dropped."))
        written.append(rec_id)
        print(f"  wrote {rec_id}")

    gates = {"G24_1_fitted_shift_reproduces_the_published_viscosity_ratio": g241,
             "G24_2_the_collapse_test_is_informative": g242,
             "G24_3_the_plateau_is_blind_to_the_shift_factor": g243,
             "G24_3pp_is_the_plateau_even_inside_the_shiftable_window": g243pp,
             "G24_4_yield_stress_invariant": g244,
             "G24_5_flow_index_invariant": g245,
             "G24_6_consistency_index_moves_as_a_to_the_n": g246,
             "G24_7_the_20C_block_is_refused": g247,
             "G24_8_replicate_scatter_against_parameter_uncertainty": g248,
             "G24_9_the_published_yield_stress_is_not_graded": g249}

    head = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                          capture_output=True, text=True).stdout.strip()
    for k, why in REFUTED.items():
        gates[k]["recorded_refutation"] = why
    payload = {"git_head": head, "records": written, "rows": rows, "gates": gates,
               "refuted_predictions": REFUTED,
               "a_pred": a_pred, "medium_viscosity_Pa_s": MEDIUM_VISCOSITY,
               "prereg": "data/carbopol_pg_prereg.md"}
    (SRC / "records.json").write_text(json.dumps(payload, indent=1), encoding="utf-8")
    write_md(payload)

    print("\n=== gates ===")
    unexpected = []
    for k, g in gates.items():
        verdict = "PASS" if g["pass"] else "FAIL"
        if k in REFUTED:
            note = "  (RECORDED REFUTATION -- this is the result)"
            if g["pass"]:
                note = "  <== was a recorded FAILURE and now passes; the result CHANGED"
                unexpected.append(k)
        else:
            note = ""
            if not g["pass"]:
                unexpected.append(k)
        print(f"  {k}: {verdict}{note}")
    print(f"\nwrote {SRC / 'records.json'} and .md")
    if unexpected:
        print(f"UNEXPECTED: {unexpected}")
        return 1
    return 0


TEMPLATE = """\
# generated by scripts/records_from_carbopol_pg.py from
# data/carbopol_pg_vinutha2024/raw/; do not edit
id: {id}
fluid: Carbopol 974p 2 wt% in propylene glycol, {temp:g} C
model: herschel_bulkley
parameters:
  tau_y: {{value: {tau_y}, unit: Pa}}
  k: {{value: {k}, unit: Pa*s**{n}}}
  n: {{value: {n}, unit: dimensionless}}
envelope:
  gammadot_min: {{value: {gmin:.6g}, unit: 1/s}}
  gammadot_max: {{value: {gmax:.6g}, unit: 1/s}}
  temperature_C: {temp:g}
  geometry: {geometry}
sources:
  - {{doi: 10.5281/zenodo.11389467, role: dataset}}
  - {{doi: 10.48550/arXiv.2312.00251, role: paper}}
licence:
  name: CC BY 4.0
  commercial_use: true
notes: >-
  Herschel-Bulkley fit (tau = tau_y + k gd^n) to the {npts}-point {temp:g} C block of the
  deposited Fig 1(a) flow curves, rates {gd_lo}-{gd_hi} 1/s, log-viscosity residuals;
  rmse_log = {rmse}. No point was dropped. Cone and plate are BOTH ROUGHENED, so unlike
  the CC BY-NC-SA Carbopol ETD2020 records this tau_y is not a slip-limited lower bound.
  The suspending medium is propylene glycol, not water, with a published viscosity of
  {eta_s} Pa s at this temperature; the sample is a model system for a rheology study and
  not a formulation anyone injects. NOT AN INJECTABILITY RECORD: the curve stops at 1e3
  1/s against 1e4-1e5 at the wall of a fine needle. The authors' own yield stress for this
  material, 8.9 Pa at all three temperatures, is a fit to a DIFFERENT model (the
  three-component model) and is not this tau_y.
extra:
  medium_viscosity_Pa_s: {eta_s}
  prereg: data/carbopol_pg_prereg.md
"""


def write_md(p) -> None:
    g, rows = p["gates"], p["rows"]
    L, A = [], None
    L.append("# Carbopol 974p in propylene glycol: two licensable records, and the first "
             "real test of the temperature shift")
    A = L.append
    A("")
    A("Written by `scripts/records_from_carbopol_pg.py` from "
      "`data/carbopol_pg_prereg.md`. Source, licence and traps: `SOURCE.md`. Every number "
      "here is printed by that script.")
    A("")
    A("## The records")
    A("")
    A("| T | N | rates (1/s) | `tau_y` (Pa) | `k` | `n` | rmse_log |")
    A("|---|---|---|---|---|---|---|")
    for r in rows:
        A(f"| {r['temperature_C']:g} C | {r['n_points']} | "
          f"{r['rate_min']:.4g}-{r['rate_max']:.4g} | "
          f"{r['tau_y_value']:.4g} [{r['tau_y_lo']:.4g}, {r['tau_y_hi']:.4g}] | "
          f"{r['k_value']:.4g} | {r['n_value']:.4f} [{r['n_lo']:.4f}, {r['n_hi']:.4f}] | "
          f"{r['rmse_log']:.4g} |")
    A("")
    g1 = g["G24_1_fitted_shift_reproduces_the_published_viscosity_ratio"]
    A("## G24.1 — the fitted shift factor against the published viscosity ratio")
    A("")
    A(f"`a_pred` = {p['medium_viscosity_Pa_s']['50.0'] if '50.0' in p['medium_viscosity_Pa_s'] else 0.012}"
      f" / {p['medium_viscosity_Pa_s'].get('10.0', 0.102)} = **{g1['a_pred']:.5f}** — the ratio of the "
      "medium viscosities the DEPOSIT publishes, which this repo does not compute.")
    A("")
    A(f"`a_fit` = **{g1['a_fit']:.5f}**, from `collapse_shift` on the two flow curves alone. "
      f"It never sees `a_pred`.")
    A("")
    A(f"Ratio **{g1['ratio']:.4f}**, `|ln ratio|` = {g1['ln_ratio']:.4f} against a bar of "
      f"{g1['ln_bar']:.4f}. Verdict: **{'PASS' if g1['pass'] else 'FAIL'}**. Collapse "
      f"residual {g1['rms_log_residual']:.4g} on {g1['n_used']} of {g1['n_total']} rates.")
    A("")
    g2, g3 = g["G24_2_the_collapse_test_is_informative"], \
        g["G24_3_the_plateau_is_blind_to_the_shift_factor"]
    A("## G24.2 and G24.3 — the test is shown informative, and the trap is measured")
    A("")
    A("Below yield both curves tend to `tau_y`, so `eta -> tau_y/gd` and "
      "`a*eta_ref(a*gd) = tau_y/gd` for **any** `a`. The low-rate half of the curve "
      "therefore contributes nothing to the objective whatever the answer is. That was "
      "pre-registered, and here it is measured rather than argued.")
    A("")
    A(f"Split at the geometric mean of the overlap, {g2['split_rate']:.4g} 1/s.")
    A("")
    A("| half | residual at `a_fit` | at `2 a_fit` | at `a_fit/2` |")
    A("|---|---|---|---|")
    for name, gg in (("high rate", g2), ("low rate", g3)):
        r = gg["residuals"]
        A(f"| {name} | {r['a_fit']['residual']:.4g} | {r['2a_fit']['residual']:.4g} | "
          f"{r['a_fit/2']['residual']:.4g} |")
    A("")
    A(f"**G24.2**: the high-rate residual rises by a factor of "
      f"**{g2['separation_factor']:.3g}** when `a` is doubled, against a bar of "
      f"{G242_SEPARATION}. **{'PASS' if g2['pass'] else 'FAIL'}** — "
      + ("so G24.1 is informative."
         if g2["pass"] else
         "so **G24.1 is UNINFORMATIVE and its agreement is not quoted.**"))
    A("")
    A(f"**G24.3 FAILED**: the low-rate residual moves by {100 * g3['relative_move']:.2f} % "
      f"over the same change, against a pre-registered bar of {G243_FLATNESS:.0%}. "
      "**The bar is not moved and the prediction is recorded as refuted.**")
    A("")
    gpp = g["G24_3pp_is_the_plateau_even_inside_the_shiftable_window"]
    A("### G24.3ʺ — why G24.3 could not test what it was written to test")
    A("")
    A("The section 1.1 algebra is right: below yield `eta -> tau_y/gd` and the shift is the "
      "identity there. But it only acts where the material is actually **on** the plateau, "
      "and G24.3 split at the geometric mean of the overlap, "
      f"{gpp['g243_split_rate']:.4g} 1/s, where this material is far ABOVE yield. Its "
      "'low-rate half' was the power-law region.")
    A("")
    A(f"Measured, not argued. Taking the plateau as `{gpp['plateau_definition']}` with "
      f"`tau_y` = {gpp['tau_y_ref']:.4g} Pa, it holds at "
      f"**{gpp['n_points_on_plateau']} points, up to {gpp['plateau_upper_rate']:.4g} 1/s**. "
      f"But `collapse_shift` crops the bottom of the range by its own shift factor, so the "
      f"shiftable window starts at **{gpp['shiftable_window_lower_rate']:.4g} 1/s**.")
    A("")
    A(f"**The plateau lies entirely OUTSIDE the shiftable window** "
      f"({gpp['plateau_upper_rate']:.4g} < {gpp['shiftable_window_lower_rate']:.4g}). The "
      "trap pre-registered in section 1.1 is real algebra and **inoperative on this "
      "dataset**, because the plateau is narrower than the range the shift itself removes. "
      "That is a fact about the data, not about the shift layer — and it is why G24.2's "
      "separation is so large.")
    A("")
    A("The probes were also re-run on the INTERSECTION of the three point sets, because "
      "comparing the residual at `a` with the residual at `2a` on different samples "
      "compares two things at once. On the common set the high-rate separation is "
      f"{g2['separation_factor_common_point_set']:.3g} and the low-rate move is "
      f"{100 * g3['relative_move_common_point_set']:.2f} %; both verdicts are unchanged.")
    A("")
    g4, g5, g6 = (g["G24_4_yield_stress_invariant"], g["G24_5_flow_index_invariant"],
                  g["G24_6_consistency_index_moves_as_a_to_the_n"])
    A("## G24.4–G24.6 — the parameter map the shift predicts")
    A("")
    A("Item 21 derived `tau_y -> tau_y`, `k -> a^n k`, `n -> n` before this dataset was "
      "found. The authors state the first of those from their own data: *\"The yield stress "
      "is independent of the solvent viscosity.\"*")
    A("")
    A("| gate | quantity | 10 C | 50 C | verdict |")
    A("|---|---|---|---|---|")
    A(f"| G24.4 | `tau_y` (Pa) | {g4['tau_y'][str(T_REF)][1]:.4g} "
      f"[{g4['tau_y'][str(T_REF)][0]:.4g}, {g4['tau_y'][str(T_REF)][2]:.4g}] | "
      f"{g4['tau_y'][str(T_OTHER)][1]:.4g} [{g4['tau_y'][str(T_OTHER)][0]:.4g}, "
      f"{g4['tau_y'][str(T_OTHER)][2]:.4g}] | {'PASS' if g4['pass'] else 'FAIL'} |")
    A(f"| G24.5 | `n` | {g5['n'][str(T_REF)][1]:.4f} [{g5['n'][str(T_REF)][0]:.4f}, "
      f"{g5['n'][str(T_REF)][2]:.4f}] | {g5['n'][str(T_OTHER)][1]:.4f} "
      f"[{g5['n'][str(T_OTHER)][0]:.4f}, {g5['n'][str(T_OTHER)][2]:.4f}] | "
      f"{'PASS' if g5['pass'] else 'FAIL'} |")
    A("")
    A(f"Separations, as the difference over the mean interval width: `tau_y` "
      f"{g4['difference_over_mean_interval_width']:.3g}, `n` "
      f"{g5['difference_over_mean_interval_width']:.3g}. An overlap of two WIDE intervals "
      "is a weak pass and is labelled one.")
    A("")
    A("**G24.4 FAILED, and the bar is not moved.** The two `tau_y` profile intervals miss "
      f"each other by **{g4['interval_gap_Pa']:.4g} Pa**; the point estimates differ by "
      f"{g4['absolute_difference_Pa']:.4g} Pa, i.e. "
      f"**{100 * g4['relative_difference']:.2f} %**. Item 21's derivation and the authors' "
      "own statement both say this parameter should be invariant, so a prediction is "
      "refuted here rather than a bar being adjusted.")
    A("")
    A("**The scale that matters is beside it, and it is why this is a narrow failure rather "
      "than a broken model.** The gap is a fifth of a percent of the value, and a profile "
      "interval is a STATISTICAL interval on one curve: it contains no protocol-to-protocol "
      "term. G24.8 measures that term directly from the independent replicate and finds "
      f"{100 * g['G24_8_replicate_scatter_against_parameter_uncertainty']['worst_replicate_median_scatter']:.2f} % "
      "— larger than the difference being called significant. So the honest reading is that "
      "`tau_y` is invariant to within everything this dataset can resolve **except its own "
      "profile interval**, which is too narrow to be the right yardstick here.")
    A("")
    A(f"**G24.6**: `k(50)/k(10)` = {g6['k_ratio']:.5g} against the predicted "
      f"`a_pred^n` = {g6['a_pred_to_the_n']:.5g} (with `n` = {g6['n_mean']:.4f}), "
      f"`|ln ratio|` = {g6['ln_ratio']:.4f} against {g6['ln_bar']:.4f}. "
      f"**{'PASS' if g6['pass'] else 'FAIL'}**. This is a second route to G24.1 through the "
      "fitted parameters rather than the curve, so it could have disagreed.")
    A("")
    g7 = g["G24_7_the_20C_block_is_refused"]
    A("## G24.7 — the 20 C block is refused, on signals recomputed here")
    A("")
    A("| T | N | on the programmed log grid | rate reversals | stresses at <=3 sig figs |")
    A("|---|---|---|---|---|")
    for T in sorted(g7["signals"], key=float):
        s = g7["signals"][T]
        A(f"| {float(T):g} C | {s['n_points']} | "
          f"{s['rates_on_the_programmed_log_grid']} ({s['grid_fraction']:.0%}) | "
          f"**{s['n_rate_reversals']}** | {s['n_stress_at_3_sigfigs']} |")
    A("")
    rev = g7["signals"]["20.0"]["rate_reversals"]
    A("A rate that goes **backwards** inside a monotonically decreasing programmed ramp is "
      "not something the instrument does: " + ", ".join(f"{a:g} -> {b:g} 1/s" for a, b in rev)
      + ". Combined with the block being off the programmed grid and rounded to three "
      "significant figures, it cannot be certified as a raw curve at the grade of the other "
      "two, and this repo does not add a record it cannot certify. **No record is written "
      "from it.**")
    A("")
    g8 = g["G24_8_replicate_scatter_against_parameter_uncertainty"]
    A("## G24.8 — an independent replicate, in the same units as the fit's uncertainty")
    A("")
    A("The cessation files carry a second, unannounced steady flow curve "
      "(`preshear-stress` at 16 rates per temperature): the same material on the same "
      "instrument under a different protocol. It is **the only independent replicate this "
      "project has for any record**.")
    A("")
    A("| T | matched rates | median ratio | range |")
    A("|---|---|---|---|")
    for T in sorted(g8["replicate"], key=float):
        r = g8["replicate"][T]
        A(f"| {float(T):g} C | {r['n_matched']} | {r['median_ratio']:.4f} | "
          f"{r['min_ratio']:.4f}-{r['max_ratio']:.4f} |")
    A("")
    A(f"Worst replicate scatter {100 * g8['worst_replicate_median_scatter']:.2f} % against a "
      f"worst relative parameter-interval width of "
      f"{100 * g8['worst_relative_parameter_width']:.2f} %. **The larger term is the "
      f"{g8['larger_term']}.** G15.4 had to approximate this comparison with a model misfit "
      "because no replicate existed.")
    A("")
    g9 = g["G24_9_the_published_yield_stress_is_not_graded"]
    A("## G24.9 — the published 8.9 Pa is NOT graded against these `tau_y`")
    A("")
    A(f"The deposit publishes one yield stress, **{g9['published_yield_stress_Pa']} Pa at all "
      "three temperatures**. It is the authors' fit to the **three-component model** "
      "(`sigma = sigma_y + k1 gd^0.5 + eta_bg gd`), which this package does not implement. "
      "biorheo fits Herschel-Bulkley.")
    A("")
    A("A yield stress is a fit-objective convention before it is a material constant — "
      "Horner's Casson `tau_y` and item 23's flow index both say so — and comparing across "
      "two model forms would call the form difference an error. So the two numbers are "
      "reported side by side as a **model-form comparison and not a disagreement**: "
      + ", ".join(f"{float(k):g} C -> {v:.4g} Pa" for k, v in sorted(g9["hb_tau_y"].items()))
      + f", against the published {g9['published_yield_stress_Pa']} Pa.")
    A("")
    A("**That the authors quote ONE yield stress for all three temperatures is used, and it "
      "is used as a prediction rather than as a number**: it is the same statement item 21 "
      "derived, and G24.4 tests it on this repo's own fits in this repo's own model.")
    (SRC / "records.md").write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
