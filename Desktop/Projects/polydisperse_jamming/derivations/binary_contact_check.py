"""binary_contact_check.py -- OUT OF SAMPLE test, on binary-mixture packs, of three contact-geometry results
from 2026-09-24 (sigma_star_geometry §5b, rattler_closure §3b-§4):

  T1 configuration model: with class degree sums D_S = n_S z_S, D_L = n_L z_L and M = D_S + D_L contact
     ends, random pairing of ends gives nc_SS = D_S^2/(2M), nc_SL = D_S D_L/M, nc_LL = D_L^2/(2M).
     Printed: measured/predicted for each pair type.
  T2 mean-degree exponent: p = ln(z_L/z_S)/ln(d_L/d_S) per record (isotropic per-grain dumps gave 1.55-1.70).
  T3 rattling ~ surface area: Poisson(Z a^2/<a^2>) per class, P(k < 2), against the per-class iterated
     rattler fractions rat_S_k2, rat_L_k2 (and the size-blind null, Poisson(Z)).

Records: out/phase12*.jsonl with a `by_class` block and both classes populated (0 < x_small_vol < 1),
converged, deduplicated by key. Diameters from psd.BinaryPSD(d10=1, ratio, x_small).

    python derivations/binary_contact_check.py  ->  out/binary_contact_check.txt
"""
from __future__ import annotations

import glob
import json
import math
import os
import sys

import numpy as np
from scipy.stats import poisson

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from psd import BinaryPSD  # noqa: E402

lines: list[str] = []


def say(s: str = "") -> None:
    print(s)
    lines.append(s)


def ternary_block() -> None:
    """T3 on TERNARY packs (phase12_tern*: `by_design_class` c0..c2 with d, n, z_mean, rattler_frac); only
    records with all three classes populated. No pair counts are recorded there, so T1 cannot be run."""
    recs = {}
    for f in sorted(glob.glob(os.path.join(ROOT, "out", "phase12_tern*.jsonl"))):
        for line in open(f, encoding="utf-8"):
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            bd = r.get("by_design_class")
            if not r.get("converged") or not bd or len(bd) != 3 or any(not bd[c]["n"] for c in bd):
                continue
            recs[r["key"]] = r
    say(f"\nTERNARY (T3 only): {len(recs)} converged records with three populated classes")
    err, err0 = [], []
    for r in recs.values():
        cl = list(r["by_design_class"].values())
        n = np.array([c["n"] for c in cl], float)
        d = np.array([c["d"] for c in cl], float)
        z = np.array([c["z_mean"] for c in cl], float)
        fr = np.array([c["rattler_frac"] for c in cl], float)
        Z = np.sum(n * z) / n.sum()
        a2 = np.sum(n * d**2) / n.sum()
        err += list(poisson.cdf(1, Z * d**2 / a2) - fr)
        err0 += list(poisson.cdf(1, Z) - fr)
    e, e0 = np.array(err), np.array(err0)
    say(f"  per-class rattler fraction: a^2 rms {np.sqrt(np.mean(e**2)):.4f} (mean {e.mean():+.4f}); size-blind rms "
        f"{np.sqrt(np.mean(e0**2)):.4f} (mean {e0.mean():+.4f}); n {e.size} class-records")


def main() -> int:
    recs = {}
    for f in sorted(glob.glob(os.path.join(ROOT, "out", "phase12*.jsonl"))):
        for line in open(f, encoding="utf-8"):
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            bc = r.get("by_class")
            if (not bc or not r.get("converged") or "mix_ratio" not in r or "x_small_vol" not in r
                    or not (0.0 < r["x_small_vol"] < 1.0) or r["mix_ratio"] <= 1.0
                    or not bc.get("n_S") or not bc.get("n_L")):
                continue
            recs[r["key"]] = r
    say(f"binary records (converged, both classes, dedup by key): {len(recs)}")
    rows = []
    for r in recs.values():
        bc = r["by_class"]
        psd = BinaryPSD(d10=1.0, ratio=r["mix_ratio"], x_small=r["x_small_vol"])
        dS, dL = psd.d_small, psd.d_large
        nS, nL = bc["n_S"], bc["n_L"]
        zS, zL = bc["z_S"], bc["z_L"]
        DS, DL = nS * zS, nL * zL
        M = DS + DL
        pred = {"SS": DS**2 / (2 * M), "SL": DS * DL / M, "LL": DL**2 / (2 * M)}
        meas = {"SS": bc["nc_SS"], "SL": bc["nc_SL"], "LL": bc["nc_LL"]}
        cons = abs((2 * meas["SS"] + meas["SL"]) - DS) + abs((2 * meas["LL"] + meas["SL"]) - DL)
        p_exp = math.log(zL / zS) / math.log(dL / dS) if zS > 0 and zL > 0 else float("nan")
        Z = (DS + DL) / (nS + nL)
        a2 = (nS * dS**2 + nL * dL**2) / (nS + nL)
        lamS, lamL = Z * dS**2 / a2, Z * dL**2 / a2
        rows.append(dict(ratio=r["mix_ratio"], x=r["x_small_vol"], ns=nS / (nS + nL),
                         rSL=meas["SL"] / pred["SL"] if pred["SL"] else np.nan,
                         rSS=meas["SS"] / pred["SS"] if pred["SS"] else np.nan,
                         rLL=meas["LL"] / pred["LL"] if pred["LL"] else np.nan,
                         p=p_exp, cons=cons,
                         ratS=bc["rat_S_k2"], ratL=bc["rat_L_k2"],
                         predS=poisson.cdf(1, lamS), predL=poisson.cdf(1, lamL),
                         nullS=poisson.cdf(1, Z), nullL=poisson.cdf(1, Z)))
    worst_cons = max(x["cons"] for x in rows)
    say(f"bookkeeping: |2 nc_SS + nc_SL - n_S z_S| + |2 nc_LL + nc_SL - n_L z_L| worst over records = {worst_cons:.3g}")
    say("\nT1 configuration model, measured/predicted contact counts by pair type, grouped by size ratio")
    say("  ratio   n    SL: median [min,max]        SS: median          LL: median")
    for ratio in sorted({x["ratio"] for x in rows}):
        g = [x for x in rows if x["ratio"] == ratio]
        sl = np.array([x["rSL"] for x in g])
        say(f"  {ratio:4.2f}  {len(g):4d}   {np.median(sl):.3f} [{sl.min():.3f}, {sl.max():.3f}]     "
            f"{np.nanmedian([x['rSS'] for x in g]):.3f}            {np.nanmedian([x['rLL'] for x in g]):.3f}")
    sl = np.array([x["rSL"] for x in rows])
    say(f"  ALL: SL measured/predicted median {np.median(sl):.3f}, IQR [{np.percentile(sl, 25):.3f}, {np.percentile(sl, 75):.3f}]")
    say("\nT2 mean-degree exponent p = ln(z_L/z_S)/ln(d_L/d_S), grouped by size ratio (median [IQR])")
    for ratio in sorted({x["ratio"] for x in rows}):
        pp = np.array([x["p"] for x in rows if x["ratio"] == ratio])
        say(f"  ratio {ratio:4.2f}: p {np.median(pp):.3f} [{np.percentile(pp, 25):.3f}, {np.percentile(pp, 75):.3f}]  n {pp.size}")
    say("\nT3 per-class rattler fraction: Poisson with a^2 mean vs size-blind, rms error over records")
    for cls in ("S", "L"):
        m = np.array([x["rat" + cls] for x in rows])
        pa = np.array([x["pred" + cls] for x in rows])
        pn = np.array([x["null" + cls] for x in rows])
        ok = np.isfinite(m)
        say(f"  class {cls}: a^2 rms {np.sqrt(np.mean((pa[ok] - m[ok]) ** 2)):.4f} (mean {np.mean(pa[ok] - m[ok]):+.4f}); "
            f"size-blind rms {np.sqrt(np.mean((pn[ok] - m[ok]) ** 2)):.4f} (mean {np.mean(pn[ok] - m[ok]):+.4f}); n {ok.sum()}")
    say("  by size ratio (a^2 rms S / L):")
    for ratio in sorted({x["ratio"] for x in rows}):
        g = [x for x in rows if x["ratio"] == ratio]
        eS = np.array([x["predS"] - x["ratS"] for x in g])
        eL = np.array([x["predL"] - x["ratL"] for x in g])
        say(f"    ratio {ratio:4.2f}: S {np.sqrt(np.mean(eS**2)):.4f} (mean {eS.mean():+.4f})   "
            f"L {np.sqrt(np.mean(eL**2)):.4f} (mean {eL.mean():+.4f})   n {len(g)}")
    ternary_block()
    with open(os.path.join(ROOT, "out", "binary_contact_check.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
