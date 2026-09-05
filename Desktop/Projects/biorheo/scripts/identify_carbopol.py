"""Fit the 15 deposited Carbopol ETD2020 curves with biorheo.identify and write the
identifiability sidecars beside the registry records.

Two things happen, and both are checks:

1. Independent reproduction. ``carbopol/fit_hb.py`` fitted ``tau = tau_y + k gd^n`` with
   ``curve_fit(sigma=tau)``; here ``dst_fit`` fits log-viscosity residuals with a
   multi-start optimiser and a different guess. The two agree to < 0.1 % on the 12
   curves with R2 >= 0.999 and to ~3.5 % on the three noisy 0.1 wt% curves (R2 0.93-0.96),
   which is the honest scale of "the same number from a different objective". The
   script refuses to write a sidecar whose interval does not contain the csv value
   (``write_identifiability`` enforces that).

2. Verdicts. Every parameter is profiled at 95 %; the sidecar carries value, interval,
   verdict and width in decades, and the loader attaches it to the record.

The reader is the same 10 lines as ``fit_hb.load`` (replicate mean of tau and of
tau/eta, latin-1, decimal comma); it is copied rather than imported because
``fit_hb.py`` runs its fit on import.

Run:  python scripts/identify_carbopol.py
"""

from __future__ import annotations

import glob
import io
import os
from pathlib import Path

import numpy as np

from biorheo.constitutive import get_model
from biorheo.identify import fit_rate_curve, write_identifiability
from biorheo.registry import load_all

HERE = Path(__file__).resolve().parent
DATA = HERE.parents[1] / "carbopol" / "data" / "etd2020_biermann2024"
HB = get_model("herschel_bulkley")


def load(path: str) -> tuple[np.ndarray, np.ndarray]:
    """(recovered rate, mean stress) for one condition, as fit_hb.load does it."""
    raw = io.open(path, encoding="latin-1").read().splitlines()
    rows = [r.split(";") for r in raw[2:] if r.strip()]
    num = lambda c: float(c.replace(",", ".")) if c.strip() else np.nan
    a = np.array([[num(c) for c in r] for r in rows])
    tau = a[:, [2, 4, 6]]                     # Pa
    eta = a[:, [3, 5, 7]] * 1e-3              # mPa.s -> Pa.s
    return np.nanmean(tau / eta, axis=1), np.nanmean(tau, axis=1)


def record_id(path: str) -> str:
    wt, tea = os.path.basename(path).replace("ETD2020-", "").replace(".csv", "").split("-")
    return f"carbopol_etd2020_{int(wt[:2]) / 10:g}wt_{tea[:2]}uL".replace(".", "p")


def main() -> None:
    records = load_all(allow_network=False)
    files = sorted(glob.glob(str(DATA / "*.csv")))
    assert len(files) == 15, files
    worst = 0.0
    print(f"{'record':36s} {'tau_y':>8s} {'k':>8s} {'n':>7s}  max |fit/csv - 1|  verdicts")
    for f in files:
        g, t = load(f)
        good = np.isfinite(g) & np.isfinite(t) & (g > 0) & (t > 0)
        print(f"  {os.path.basename(f)}: N = {g.size} -> {good.sum()} after finite/positive mask")
        g, t = g[good], t[good]
        rec = records[record_id(f)]
        ident = fit_rate_curve(HB, g, t / g, temperature_C=rec.envelope.temperature_C,
                               geometry=rec.envelope.geometry)
        dev = np.max(np.abs(ident.values / rec.params_si - 1.0))
        worst = max(worst, dev)
        write_identifiability(rec, ident, data_note=(
            "fit of the replicate-mean flow curve from the deposited csv, rate recovered "
            "as tau/eta; stored parameters come from carbopol/fit_hb.py (curve_fit, "
            f"sigma = tau) and differ from this fit by at most {dev:.2%}"))
        v = ident.values
        print(f"{rec.id:36s} {v[0]:8.3f} {v[1]:8.3f} {v[2]:7.4f}  {dev:16.2%}  "
              + ", ".join(f"{k} {e.verdict} ({e.decades:.2f} dec)" for k, e in ident.params.items()))
    print(f"worst deviation from the csv over 15 records: {worst:.2%}")


if __name__ == "__main__":
    main()
