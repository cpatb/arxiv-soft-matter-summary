"""firstpassage.py — the first-passage clock, and the one place it is implemented.

PRE-REGISTRATION: PREREG_FIRSTPASSAGE_2026-09-13.md.

WHY THIS OBSERVABLE
-------------------
`grade_strainclock.clock` measures the strain to reach **half of the cell's own `p_ss`**.
That is not well posed near the critical line: `p_ss` is itself collapsing as
`phi -> phi_c`, so the target moves, and the strain-clock campaign found 0 of 16 jammed
cells readable at `x` <= +0.10 (`PREREG_STRAINCLOCK_2026-09-13.md` §7, G2).

`t_fp` fixes the target instead. It is the strain at which `p_hat` first reaches an
ABSOLUTE threshold, so it needs no `p_ss`, no monotonicity assumption, and no averaging
window. What it needs instead is a check that the passage is not a transient excursion,
and that is `frac_above`: the fraction of post-passage samples still at or above the
threshold. A cell that pokes through the threshold and falls back is refused, with its
number printed.

⚠ The threshold is a CONVENTION and the exponent it returns depends on it (0.08-0.37 in
`a` between 1e-04 and 3e-04 on the exploratory data). Both are computed everywhere and
the spread between them is the floor below which no exponent difference is meaningful.
This module never picks one.
"""

from __future__ import annotations

import glob
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

# The two thresholds, both reported everywhere. 1e-04 is the primary (§3: it is the one
# where every x = +0.10 cell already satisfies frac_above > 0.9); 3e-04 is `p*`, the
# pressure the constant-pressure ladder held and `phi_c` was defined at.
THRESHOLDS = (1e-4, 3e-4)
P_PRIMARY = 1e-4

FRAC_ABOVE_MIN = 0.90        # F1
P_PREP_MAX = 1e-6            # F1
FIT_RUNGS = (0.35, 0.60, 0.85)      # F2: the fit sees only these
PREDICT_RUNGS = (0.05, 0.10)        # F2/F4: these are predicted, never fitted

PHI_MU = {0.0: 0.591127, 0.2: 0.601989, 0.5: 0.644447}


def first_passage(gamma, p, thr: float) -> tuple[float, float]:
    """(strain at first passage to `thr`, fraction of post-passage samples >= `thr`).

    Linear interpolation between the two bracketing samples. Returns (nan, nan) when the
    threshold is never reached -- a pack that never builds that much pressure is a
    result, not an error.
    """
    g = np.asarray(gamma, float)
    y = np.asarray(p, float)
    if g.size < 3 or not np.isfinite(thr) or thr <= 0:
        return float("nan"), float("nan")
    hit = y >= thr
    if not hit.any():
        return float("nan"), float("nan")
    i = int(np.argmax(hit))
    frac = float(np.mean(hit[i:]))
    if i == 0:
        return float(g[0]), frac
    y0, y1 = y[i - 1], y[i]
    t = float(g[i]) if y1 == y0 else float(g[i - 1] + (thr - y0) * (g[i] - g[i - 1])
                                           / (y1 - y0))
    return t, frac


def load_records(patterns=("phase10_sc_shard*.jsonl", "phase10_fp*_shard*.jsonl")):
    """Every strain-clock and first-passage record, deduplicated by cell key.

    The two campaigns share a driver, a protocol string, a rate, `n_samples` and N, and
    their keys carry `phi` and `seed`, so they pool by construction. `dedupe` is not
    cosmetic: a restarted shard re-writes nothing, but a re-LAUNCHED one would.
    """
    seen: set[str] = set()
    recs, files = [], []
    for pat in patterns:
        files += sorted(glob.glob(os.path.join(OUT, pat)))
    dupes = 0
    for f in files:
        for line in open(f, encoding="utf-8"):
            d = json.loads(line)
            if d["key"] in seen:
                dupes += 1
                continue
            seen.add(d["key"])
            recs.append(d)
    return recs, files, dupes


def cell_table(recs) -> list[dict]:
    """One row per cell, with `t_fp` and `frac_above` at every threshold."""
    rows = []
    for d in recs:
        if d.get("error"):
            continue
        r = {
            "cv": d["cv"], "seed": d["seed"], "x": d.get("x_rung", float("nan")),
            "phi_target": d["phi_target"],
            "dphi": d["phi_target"] - PHI_MU[d["cv"]],
            "p_prep": float(d.get("p_hat_prep", float("nan"))),
            "p_ss": float(d.get("p_ss_hat", float("nan"))),
            "z_nr_ss": d.get("z_nr_ss", float("nan")),
            "gdot": d.get("gdot", float("nan")),
            "protocol": d.get("protocol", ""), "mu_prep": d.get("mu_prep", float("nan")),
            "fswitch": bool(d.get("friction_switch_ok", False)),
            "n_samples": d.get("n_samples", 0), "n_requested": d.get("n_requested", 0),
            "key": d["key"],
        }
        for thr in THRESHOLDS:
            t, frac = first_passage(d.get("gamma") or [], d.get("p_hat_series") or [], thr)
            r[f"t_{thr:g}"] = t
            r[f"frac_{thr:g}"] = frac
        rows.append(r)
    return rows


def admitted(r: dict, thr: float) -> bool:
    """F1, and it is the whole admission rule."""
    return (math.isfinite(r[f"t_{thr:g}"])
            and r[f"frac_{thr:g}"] >= FRAC_ABOVE_MIN
            and r["p_prep"] < P_PREP_MAX)


def rung_means(rows, cv: float, thr: float) -> dict[float, tuple[float, float, int]]:
    """{x: (geometric mean t_fp, sd of ln t_fp, n)} over the admitted cells."""
    out = {}
    for x in sorted({round(r["x"], 4) for r in rows if r["cv"] == cv}):
        v = [math.log(r[f"t_{thr:g}"]) for r in rows
             if r["cv"] == cv and abs(r["x"] - x) < 1e-9 and admitted(r, thr)]
        if v:
            out[x] = (float(np.exp(np.mean(v))),
                      float(np.std(v, ddof=1)) if len(v) > 1 else float("nan"), len(v))
    return out


def powerlaw_fit(rows, cv: float, thr: float, rungs=FIT_RUNGS):
    """ln t_fp = ln A - a ln(dphi), fitted on `rungs` only, per CELL not per rung mean.

    Per cell, so every admitted seed carries its own weight and a rung with fewer
    admitted cells does not count the same as a full one. Returns (a, lnA, n_cells).
    """
    xs, ys = [], []
    for r in rows:
        if r["cv"] != cv or not admitted(r, thr):
            continue
        if not any(abs(r["x"] - x) < 1e-9 for x in rungs):
            continue
        if r["dphi"] <= 0:
            continue
        xs.append(math.log(r["dphi"]))
        ys.append(math.log(r[f"t_{thr:g}"]))
    if len(xs) < 3:
        return float("nan"), float("nan"), len(xs)
    b, c = np.polyfit(xs, ys, 1)
    return float(-b), float(c), len(xs)


def bootstrap_exponent(rows, cv: float, thr: float, rungs=FIT_RUNGS,
                       n_boot: int = 2000, seed: int = 20260913):
    """se of the exponent by resampling SEEDS within each rung.

    Resampling seeds and not cells: the seeds are the independent replicates and the
    rungs are the design, so a bootstrap over all cells pooled would let a rung vanish
    and change the leverage rather than the noise.
    """
    rng = np.random.default_rng(seed)
    by_rung: dict[float, list[dict]] = {}
    for r in rows:
        if r["cv"] == cv and admitted(r, thr) and r["dphi"] > 0:
            if any(abs(r["x"] - x) < 1e-9 for x in rungs):
                by_rung.setdefault(round(r["x"], 4), []).append(r)
    if len(by_rung) < 2:
        return float("nan")
    out = []
    for _ in range(n_boot):
        xs, ys = [], []
        for rr in by_rung.values():
            pick = rng.integers(0, len(rr), len(rr))
            for i in pick:
                xs.append(math.log(rr[i]["dphi"]))
                ys.append(math.log(rr[i][f"t_{thr:g}"]))
        if len(set(xs)) < 2:
            continue
        b, _ = np.polyfit(xs, ys, 1)
        out.append(-b)
    return float(np.std(out, ddof=1)) if len(out) > 10 else float("nan")


def pss_exponent(rows, cv: float, thr: float, rungs=FIT_RUNGS):
    """`b` in `p_ss ~ dphi^b`, fitted on the SAME admitted cells as `powerlaw_fit`.

    Same sample by construction, because F6's statistic is the ratio a/b and a ratio of
    two exponents fitted on different cells is not a ratio of anything. `thr` enters only
    through the admission rule.
    """
    xs, ys = [], []
    for r in rows:
        if r["cv"] != cv or not admitted(r, thr) or r["dphi"] <= 0:
            continue
        if not any(abs(r["x"] - x) < 1e-9 for x in rungs):
            continue
        if not (math.isfinite(r["p_ss"]) and r["p_ss"] > 0):
            continue
        xs.append(math.log(r["dphi"]))
        ys.append(math.log(r["p_ss"]))
    if len(xs) < 3:
        return float("nan"), len(xs)
    b, _ = np.polyfit(xs, ys, 1)
    return float(b), len(xs)


def bootstrap_ratio(rows, cv: float, thr: float, rungs=FIT_RUNGS,
                    n_boot: int = 2000, seed: int = 20260913):
    """se of a/b by resampling SEEDS within each rung, refitting BOTH exponents on the
    SAME resample -- so the numerator's and denominator's shared sampling noise cancels
    exactly as it does in the point estimate."""
    rng = np.random.default_rng(seed)
    by_rung: dict[float, list[dict]] = {}
    for r in rows:
        if (r["cv"] == cv and admitted(r, thr) and r["dphi"] > 0
                and math.isfinite(r["p_ss"]) and r["p_ss"] > 0
                and any(abs(r["x"] - x) < 1e-9 for x in rungs)):
            by_rung.setdefault(round(r["x"], 4), []).append(r)
    if len(by_rung) < 2:
        return float("nan")
    out = []
    for _ in range(n_boot):
        xs, yt, yp = [], [], []
        for rr in by_rung.values():
            for i in rng.integers(0, len(rr), len(rr)):
                xs.append(math.log(rr[i]["dphi"]))
                yt.append(math.log(rr[i][f"t_{thr:g}"]))
                yp.append(math.log(rr[i]["p_ss"]))
        if len(set(xs)) < 2:
            continue
        a = -np.polyfit(xs, yt, 1)[0]
        b = np.polyfit(xs, yp, 1)[0]
        if b != 0:
            out.append(a / b)
    return float(np.std(out, ddof=1)) if len(out) > 10 else float("nan")
