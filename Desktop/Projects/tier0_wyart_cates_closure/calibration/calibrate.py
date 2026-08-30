"""Calibrate the tier-0 closure against three published flow-curve datasets.

This closes the README's open item, "calibrate ``(phi, phi_m, sigma_star,
sigma_yield, f_res, tau)`` to real flow curves", and it closes it with a
negative for three of the six: a steady flow curve contains no information about
``tau``, ``f_res`` or ``sigma_yield`` at all, and only two of the remaining
three transfer into a solver.

Reuse, not reimplementation
---------------------------
The readers, the fitter, the uncertainty machinery and the four per-dataset
fixed constants all come from ``dst_fit``, which already carries them with their
provenance and their tests. Nothing about a flow curve is re-derived here. What
is new is the *transfer* question: given a fit, what may a solver take from it.

The datasets
------------
Three Edinburgh DataShare deposits, chosen because in all three the volume (or
weight) fraction is known independently of any rheology:

- Guy, Hermes & Poon, PRL 115, 088304 (2015), PMMA at d = 404 nm and 3770 nm.
  Two sizes of the same chemistry, 9.3x apart, which is what makes the transfer
  map reproducible rather than anecdotal.
- Richards, Royer, Liebchen, Guy & Poon, PRL 123, 038004 (2019), corn starch in
  50 wt% glycerol-water. Different lab, different material, weight fractions.

Run:  python calibration/calibrate.py
Writes ``calibration/out/calibration.json`` and prints the tables.
"""
from __future__ import annotations

import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DST_FIT = Path(r"C:\Users\compl\Desktop\Projects\dst_fit")
sys.path.insert(0, str(DST_FIT / "src"))
sys.path.insert(0, str(HERE))

from dst_fit import fit_model, thickening_branch  # noqa: E402
from dst_fit.io_literature import (  # noqa: E402
    GUY_FIG1A_PHI,
    GUY_FIG1B_PHI,
    GUY_PMMA_404,
    GUY_PMMA_3770,
    RICHARDS_CORNSTARCH,
    read_guy_flow_curves,
    read_richards_series,
)
from dst_fit.models import WYART_CATES, WYART_CATES_KNOWN_PHI  # noqa: E402

from brownian import PE_NON_BROWNIAN, invert_plateau, plateau_admissible  # noqa: E402
from steady import is_steady  # noqa: E402
from transfer import (  # noqa: E402
    MaterialSpec,
    jamming_stress,
    jams_under_stress,
    transfer_map,
    verdict_band,
)

GUY_DIR = (DST_FIT / "Towards a unified description of the rheology of hard-particle suspensions"
           / "Fig 1" / "Fig 1")
RICHARDS_DIR = (DST_FIT / "Competing Timescales Lead to Oscillations in Shear-Thickening Suspensions"
                / "Figure1" / "Plotting")

# Particle radii. Guy et al. state d = 404 nm and d = 3770 nm in the caption of
# the figure these curves come from.
R_GUY_404 = 404e-9 / 2
R_GUY_3770 = 3770e-9 / 2

# Richards et al. do not state a granule size, and this calibration will not
# invent one. Their paper calls the system "non-Brownian (or granular)" in its
# first sentence and uses weight fractions "due to the porosity of our main
# model system (cornstarch)". A conservative lower bound on a corn-starch
# granule radius is 2.5 um (granules run roughly 5 to 20 um in diameter); the
# Peclet number scales as a^3, so a bound on the radius is a bound on Pe, and if
# the bound already clears the threshold the unknown size does not matter.
R_CORNSTARCH_LOWER_BOUND = 2.5e-6
# The deposit records no temperature either. 293 K is used only inside a Peclet
# bound that clears its threshold by many decades, where a few K cannot matter.
T_ROOM = 293.15

N_STARTS = 16


def _spec(cset, radius, temperature) -> MaterialSpec:
    """Lift a dst_fit ConstantSet into a MaterialSpec, adding size and temperature."""
    return MaterialSpec(
        name=cset.name, eta_s=cset.eta_s, phi_0=cset.phi_0, n=cset.n, beta=cset.beta,
        phi_units=cset.phi_units, radius=radius, temperature=temperature,
        phi_m_published=cset.phi_m_published, source=cset.source,
    )


SPECS = {
    "guy_404": _spec(GUY_PMMA_404, R_GUY_404, 292.15),
    "guy_3770": _spec(GUY_PMMA_3770, R_GUY_3770, 292.15),
    "richards_cornstarch": _spec(RICHARDS_CORNSTARCH, R_CORNSTARCH_LOWER_BOUND, T_ROOM),
}


@dataclass
class CurveCalibration:
    """All three protocols applied to one flow curve."""

    label: str
    phi_known: float
    # protocol C, fit-free
    pe_lowest: float
    plateau_admitted: bool
    plateau_reason: str
    eta_r_lowest: float
    phi_plateau: float
    # protocol A, free fit
    phi_free: float
    phi_free_err: float
    phi_m_free: float
    phi_m_free_err: float
    sigma_star_free: float
    sigma_star_free_err: float
    rmse_free: float
    at_bound_free: tuple[str, ...]
    # protocol B, phi pinned at the known value
    phi_m_pinned: float
    phi_m_pinned_err: float
    sigma_star_pinned: float
    sigma_star_pinned_err: float
    rmse_pinned: float
    at_bound_pinned: tuple[str, ...]
    # what a solver would conclude
    jams_free: bool
    jams_pinned: bool
    sigma_jam_free: float
    sigma_jam_pinned: float


def calibrate_curve(sweep, spec: MaterialSpec) -> CurveCalibration | None:
    """Run all three protocols on one flow curve. Returns None if it cannot be fit."""
    if sweep.phi is None:
        return None
    free_model = WYART_CATES.with_fixed(**spec.as_fixed())
    branch = thickening_branch(sweep.viscosity)
    try:
        free = fit_model(free_model, sweep.stress, sweep.viscosity,
                         branch=branch, n_starts=N_STARTS)
    except (RuntimeError, ValueError):
        return None

    pinned_model = WYART_CATES_KNOWN_PHI.with_fixed(phi=sweep.phi, **spec.as_fixed())
    try:
        pinned = fit_model(pinned_model, sweep.stress, sweep.viscosity,
                           branch=branch, n_starts=N_STARTS)
    except (RuntimeError, ValueError):
        return None

    verdict = plateau_admissible(
        sweep.stress, sweep.shear_rate, sweep.viscosity,
        radius=spec.radius, eta_s=spec.eta_s, temperature=spec.temperature,
        pe_min=PE_NON_BROWNIAN,
    )
    phi_plateau = (float(invert_plateau(np.array([verdict.eta_r]), spec.phi_0, spec.n)[0])
                   if verdict.admitted else float("nan"))

    p_free = dict(zip(free.model.param_names, free.params))
    e_free = dict(zip(free.model.param_names, free.stderr))
    p_pin = dict(zip(pinned.model.param_names, pinned.params))
    e_pin = dict(zip(pinned.model.param_names, pinned.stderr))

    return CurveCalibration(
        label=sweep.sheet,
        phi_known=float(sweep.phi),
        pe_lowest=verdict.pe,
        plateau_admitted=verdict.admitted,
        plateau_reason=verdict.reason,
        eta_r_lowest=verdict.eta_r,
        phi_plateau=phi_plateau,
        phi_free=p_free["phi"], phi_free_err=e_free["phi"],
        phi_m_free=p_free["phi_m"], phi_m_free_err=e_free["phi_m"],
        sigma_star_free=p_free["sigma_star"], sigma_star_free_err=e_free["sigma_star"],
        rmse_free=free.rmse_log, at_bound_free=free.parameters_at_bound(),
        phi_m_pinned=p_pin["phi_m"], phi_m_pinned_err=e_pin["phi_m"],
        sigma_star_pinned=p_pin["sigma_star"], sigma_star_pinned_err=e_pin["sigma_star"],
        rmse_pinned=pinned.rmse_log, at_bound_pinned=pinned.parameters_at_bound(),
        jams_free=jams_under_stress(p_free["phi"], p_free["phi_m"]),
        jams_pinned=jams_under_stress(sweep.phi, p_pin["phi_m"]),
        sigma_jam_free=jamming_stress(p_free["phi"], p_free["phi_m"], spec.phi_0,
                                      p_free["sigma_star"], spec.beta),
        sigma_jam_pinned=jamming_stress(sweep.phi, p_pin["phi_m"], spec.phi_0,
                                        p_pin["sigma_star"], spec.beta),
    )


def load_series(key: str):
    if key == "guy_404":
        return read_guy_flow_curves(GUY_DIR / "Fig 1a" / "Fig1a.xls", GUY_FIG1A_PHI)
    if key == "guy_3770":
        return read_guy_flow_curves(GUY_DIR / "Fig 1b" / "Fig1b.xlsx", GUY_FIG1B_PHI)
    if key == "richards_cornstarch":
        return read_richards_series(RICHARDS_DIR)
    raise KeyError(key)


def calibrate_series(key: str) -> tuple[MaterialSpec, list[CurveCalibration], list[str]]:
    """Fit every curve. Returns the spec, the successes, and the REFUSALS by name.

    The refusals are returned rather than filtered away because one of them --
    Guy 404's 40% sheet, whose thickening branch is a single point -- is the
    lowest-phi curve of a five-point regression, i.e. its highest-leverage
    point. A count of 5 with no note that a sixth was refused reads as a
    complete series.
    """
    spec = SPECS[key]
    ok, refused = [], []
    for sweep in load_series(key):
        row = calibrate_curve(sweep, spec)
        if row is None:
            refused.append(sweep.sheet)
        else:
            ok.append(row)
    return spec, ok, refused


def _usable(rows: list[CurveCalibration]) -> list[CurveCalibration]:
    """Curves whose free fit is a measurement, not a constrained optimum.

    This is protocol A's admission rule and it must not be applied to protocol
    C. ``at_bound_free`` is an *optimiser* output; using it to select the sample
    for the fit-free plateau route would make that route's independence a
    fiction. Guy's 19% curve is the case in point: it inverts perfectly well
    (Pe = 157, phi_plateau = 0.165) and is excluded only because its free fit
    railed. Including it moves that series' plateau slope from 1.152 to 1.020.
    """
    return [r for r in rows if not r.at_bound_free]


def _plateau_usable(rows: list[CurveCalibration]) -> list[CurveCalibration]:
    """Protocol C's own rule: non-Brownian at the lowest stress, and invertible.

    Nothing about the optimiser enters.
    """
    return [r for r in rows
            if r.plateau_admitted and np.isfinite(r.phi_plateau)]


def _steady(key: str, rows: list[CurveCalibration]) -> list[CurveCalibration]:
    """Usable curves the *source lab* calls steady. See calibration/steady.py."""
    return [r for r in _usable(rows) if is_steady(key, r.phi_known)]


def _band(rows: list[CurveCalibration], attr: str) -> dict | None:
    if len(rows) < 2:
        return None
    b = verdict_band([getattr(r, attr) for r in rows], [r.label for r in rows])
    return {"lo": b.lo, "hi": b.hi, "width": b.width, "n": len(rows)}


def summarise(key: str, spec: MaterialSpec, rows: list[CurveCalibration]) -> dict:
    usable = _usable(rows)
    out: dict = {
        "dataset": key,
        "material": spec.name,
        "source": spec.source,
        "phi_units": spec.phi_units,
        "constants": {"eta_s": spec.eta_s, "phi_0": spec.phi_0,
                      "n": spec.n, "beta": spec.beta},
        "phi_m_published": spec.phi_m_published,
        "n_curves": len(rows),
        "n_usable": len(usable),
        "curves": [asdict(r) for r in rows],
    }
    steady = _steady(key, rows)
    steady_all = [r for r in rows if is_steady(key, r.phi_known)]
    out["n_steady"] = len(steady)
    out["steady_labels"] = [r.label for r in steady]

    for tag, subset, plateau_pool in (("", usable, rows), ("_steady", steady, steady_all)):
        if len(subset) >= 3:
            out[f"transfer_free{tag}"] = asdict(transfer_map(
                [r.phi_known for r in subset], [r.phi_free for r in subset]))
        plateau = _plateau_usable(plateau_pool)
        out[f"n_plateau_admitted{tag}"] = len(plateau)
        out[f"plateau_labels{tag}"] = [r.label for r in plateau]
        if len(plateau) >= 3:
            out[f"transfer_plateau{tag}"] = asdict(transfer_map(
                [r.phi_known for r in plateau], [r.phi_plateau for r in plateau]))
        out[f"verdict_band_free{tag}"] = _band(subset, "phi_m_free")
        out[f"verdict_band_pinned{tag}"] = _band(subset, "phi_m_pinned")

    # Finding 9: the steady 4-subset is also the LOWEST-phi 4-subset, and phi_m
    # rises with phi, so it is the tightest block available by construction.
    # The null is what a random 4-subset of the same curves gives.
    if len(usable) >= 5:
        out["band_null_random4"] = _random_subset_band_null(usable, k=len(steady) or 4)
    return out


def _random_subset_band_null(rows: list[CurveCalibration], k: int) -> dict:
    """Mean and minimum phi_m range over every k-subset of these curves.

    The steady restriction removes the highest-phi curves, and phi_m rises with
    phi, so the steady block is the tightest k-subset available whatever
    steadiness has to do with it. Reporting the null says how much of the
    narrowing is the restriction and how much is arithmetic.
    """
    from itertools import combinations
    vals = np.array([r.phi_m_free for r in rows])
    ranges = [float(np.ptp(vals[list(c)])) for c in combinations(range(len(vals)), k)]
    return {"k": k, "n_subsets": len(ranges), "mean": float(np.mean(ranges)),
            "min": float(np.min(ranges)), "max": float(np.max(ranges))}


def _fmt(x, w=8, p=4):
    return f"{x:{w}.{p}f}" if np.isfinite(x) else f"{'--':>{w}s}"


def print_series(key: str, spec: MaterialSpec, rows: list[CurveCalibration]) -> None:
    print(f"\n{'=' * 100}\n{spec.name}")
    print(f"  eta_s = {spec.eta_s:g} Pa.s, phi_0 = {spec.phi_0:.4f}, n = {spec.n}, "
          f"beta = {spec.beta}, phi in {spec.phi_units} fraction")
    print(f"  {'curve':>7s} {'phi':>7s} | {'Pe(lo)':>9s} {'phi_plat':>8s} | "
          f"{'phi_free':>8s} {'phi_m':>8s} {'sigma*':>9s} | "
          f"{'phi_m|phi':>9s} {'sigma*':>9s} | {'jams?':>11s}")
    for r in rows:
        flag = " [bound]" if r.at_bound_free else ""
        jam = f"{'Y' if r.jams_free else 'n'}/{'Y' if r.jams_pinned else 'n'}"
        print(f"  {r.label:>7s} {r.phi_known:7.4f} | {r.pe_lowest:9.2e} "
              f"{_fmt(r.phi_plateau)} | {_fmt(r.phi_free)} {_fmt(r.phi_m_free)} "
              f"{r.sigma_star_free:9.2f} | {_fmt(r.phi_m_pinned)} "
              f"{r.sigma_star_pinned:9.2f} | {jam:>11s}{flag}")

    usable = _usable(rows)
    if len(usable) >= 3:
        tm = transfer_map([r.phi_known for r in usable], [r.phi_free for r in usable])
        print(f"\n  free-fit transfer map    : {tm}")
        print(f"    identity (slope 1, intercept 0) inside both intervals? "
              f"{'YES' if tm.is_identity else 'NO'}")
        plateau = [r for r in usable if r.plateau_admitted and np.isfinite(r.phi_plateau)]
        print(f"    plateau route admitted on {len(plateau)}/{len(usable)} usable curves "
              f"(Pe >= {PE_NON_BROWNIAN:g})")
        if len(plateau) >= 3:
            tp = transfer_map([r.phi_known for r in plateau],
                              [r.phi_plateau for r in plateau])
            print(f"  plateau transfer map     : {tp}")
    if len(usable) >= 2:
        vf = verdict_band([r.phi_m_free for r in usable])
        vp = verdict_band([r.phi_m_pinned for r in usable])
        print("\n  phi_m across this material's own curves")
        print(f"    free fit  : {vf}")
        print(f"    phi pinned: {vp}")
        if spec.phi_m_published is not None:
            inside_f = "inside" if vf.contains(spec.phi_m_published) else "OUTSIDE"
            inside_p = "inside" if vp.contains(spec.phi_m_published) else "OUTSIDE"
            print(f"    published phi_m = {spec.phi_m_published:.3f}: "
                  f"{inside_f} the free band, {inside_p} the pinned band")


def main() -> int:
    out = {}
    for key in SPECS:
        spec, rows, refused = calibrate_series(key)
        print_series(key, spec, rows)
        if refused:
            print(f"\n  REFUSED by the fitter (not counted anywhere else): {refused}")
        out[key] = summarise(key, spec, rows)
        out[key]["refused_by_fitter"] = refused
    dest = HERE / "out"
    dest.mkdir(exist_ok=True)
    (dest / "calibration.json").write_text(json.dumps(out, indent=2, default=str),
                                           encoding="utf-8")
    print(f"\nwrote {dest / 'calibration.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
