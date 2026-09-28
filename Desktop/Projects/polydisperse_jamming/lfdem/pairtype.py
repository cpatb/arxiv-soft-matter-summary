"""pairtype.py -- grade PREREG_LFSIG2 Amendment 7: does a pair's activation stress scale with its reduced radius as
R*^-1 (Derjaguin threshold, load ~ sigma R*^2) or R*^-2 (size-independent critical load, Guy et al. 2020)?

    python lfdem/pairtype.py   ->  out/grade_lfsig2bi.txt

Binary suspension, size ratio 2 (radii 1 and 2), small species 50 % of the solid volume, phi 0.55, mu 0.3, 3 seeds, stress
ladder {1, 1.5, 2, 3, 4.5, 7, 10, 15, 100}. Pair types SS (R* 0.5), SL (R* 2/3), LL (R* 1).
Per cell, interaction snapshots in strain [2, 4]: for each pair type, the fraction of LISTED pairs of that type (every
pair inside lub_max_gap) that touch (state >= 1). f_type(sigma) normalised by its own sigma-100 plateau. sigma_half(type)
= the stress where the normalised fraction crosses 1/2, log-linear between bracketing rungs (NaN if not monotone or
not bracketed). Per seed: ln[sigma_half(SS)/sigma_half(LL)] (predicted ln 2 = 0.693 Derjaguin, ln 4 = 1.386 constant
force) and ln[sigma_half(SL)/sigma_half(LL)] (ln 1.5 = 0.405 vs ln 2.25 = 0.811).
Registered verdict B1 (95 %, df 2): INFORMATIVE iff se < (ln 4 - ln 2)/(2 x 4.303) = 0.0805; then DERJAGUIN if the mean is
within 4.303 se of ln 2 and not of ln 4, CONSTANT-FORCE for the reverse, BOTH/NEITHER otherwise (NEITHER printed with the
measured exponent m_eff = mean/ln 2, where Derjaguin is m 1 and constant force m 2).
Health: rc 0 (or 124 after the window), strain >= 3.95, overlap < 0.05 (lfsig2's rule).
"""
from __future__ import annotations

import glob
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from mech_lfsig import snapshots  # noqa: E402

OUT = os.path.join(os.path.dirname(HERE), "out")
W0, W1 = 2.0, 4.0
SIGMAS = (1.0, 1.5, 2.0, 3.0, 4.5, 7.0, 10.0, 15.0, 100.0)
SEEDS = (1, 2, 3)
T_CRIT = 4.303
SE_MAX = (math.log(4) - math.log(2)) / (2 * T_CRIT)
lines: list[str] = []


def say(s: str = "") -> None:
    print(s)
    lines.append(s)


def cell(d: str) -> dict | None:
    meta = open(os.path.join(d, "META"), encoding="utf-8").read()
    kv = dict(re.findall(r"(\w+)=(\S+)", meta))
    data = glob.glob(os.path.join(d, "data_*.dat"))
    ints = glob.glob(os.path.join(d, "int_*.dat"))
    if not data or not ints or not os.path.exists(os.path.join(d, "conf.dat")):
        return None
    a = np.loadtxt(data[0], comments="#", ndmin=2)
    rc = int(kv["rc"]) if "rc" in kv else 124
    w = (a[:, 1] >= W0) & (a[:, 1] <= W1)
    overlap = float(-a[w, 15].min()) if w.any() else float("nan")
    healthy = rc in (0, 124) and a[:, 1].max() >= W1 - 0.05 and overlap < 0.05
    r = np.loadtxt(os.path.join(d, "conf.dat"), comments="#", ndmin=2)[:, 3]
    small = r < 1.5
    tot = {"SS": 0, "SL": 0, "LL": 0}
    hit = {"SS": 0, "SL": 0, "LL": 0}
    for strain, _rate, s in snapshots(ints[0]):
        if not (W0 <= strain <= W1):
            continue
        i, j = s[:, 0].astype(int), s[:, 1].astype(int)
        ns = small[i].astype(int) + small[j].astype(int)
        t = s[:, 2] >= 1
        for key, k in (("SS", 2), ("SL", 1), ("LL", 0)):
            m = ns == k
            tot[key] += int(m.sum())
            hit[key] += int((m & t).sum())
    frac = {k: hit[k] / tot[k] if tot[k] else float("nan") for k in tot}
    return {"label": kv["LABEL"], "sig": float(kv["SIG"]), "seed": int(kv["SEED"]), "rc": rc, "healthy": healthy,
            "overlap": overlap, "frac": frac, "tot": tot}


def shalf(sig: np.ndarray, f: np.ndarray) -> float:
    if not np.all(np.diff(f) >= -1e-12):
        return float("nan")
    for k in range(1, len(f)):
        if f[k - 1] < 0.5 <= f[k]:
            t = (0.5 - f[k - 1]) / (f[k] - f[k - 1])
            return float(math.exp(math.log(sig[k - 1]) + t * (math.log(sig[k]) - math.log(sig[k - 1]))))
    return float("nan")


def shalf_first(sig: np.ndarray, f: np.ndarray) -> float:
    """POST HOC (added 2026-09-28 04:45, after seed 1 showed the LL series dipping by ~1 binomial s.e. AFTER it crossed):
    the FIRST crossing of 1/2, log-linear, without the registered monotonicity requirement."""
    for k in range(1, len(f)):
        if f[k - 1] < 0.5 <= f[k]:
            t = (0.5 - f[k - 1]) / (f[k] - f[k - 1])
            return float(math.exp(math.log(sig[k - 1]) + t * (math.log(sig[k]) - math.log(sig[k - 1]))))
    return float("nan")


def post_hoc(by: dict) -> None:
    say("\n  POST HOC (not registered): first-crossing sigma_half per type (no monotonicity requirement); binomial s.e. of the")
    say("  plateau-normalised LL fraction shown so the size of the dips that void the registered readout can be judged")
    vals = []
    for sd in SEEDS:
        if not all((sd, s) in by for s in SIGMAS):
            continue
        s = np.array(SIGMAS)
        sh = {}
        for key in ("SS", "SL", "LL"):
            f = np.array([by[(sd, x)]["frac"][key] for x in SIGMAS])
            sh[key] = shalf_first(s, f / f[-1])
        n_ll = np.mean([by[(sd, x)]["tot"]["LL"] for x in SIGMAS])
        f_ll = by[(sd, 100.0)]["frac"]["LL"]
        se_ll = math.sqrt(f_ll * (1 - f_ll) / n_ll) / f_ll
        say(f"    seed {sd}: SS {sh['SS']:.3f}  SL {sh['SL']:.3f}  LL {sh['LL']:.3f}  -> SS/LL {sh['SS'] / sh['LL']:.3f}, "
            f"SL/LL {sh['SL'] / sh['LL']:.3f}  (LL pairs per dump ~{n_ll:.0f}, normalised binomial s.e. ~{se_ll:.3f})")
        vals.append(math.log(sh["SS"] / sh["LL"]))
    if len(vals) >= 2:
        v = np.array(vals); m = float(v.mean()); se = float(v.std(ddof=1) / math.sqrt(len(v)))
        say(f"    SS/LL first-crossing: ratio {math.exp(m):.3f}, ln se {se:.4f} (n {len(v)}); m_eff = ln ratio / ln 2 = {m / math.log(2):.2f} "
            f"(Derjaguin 1, constant force 2)")


def main() -> int:
    say("pairtype -- PREREG_LFSIG2 Amendment 7: per-pair-type activation in a binary suspension (ratio 2, x_small 0.5, phi 0.55)")
    dirs = sorted(glob.glob(os.path.join(OUT, "lfsig2bi", "bi_*")))
    cells = [c for c in (cell(d) for d in dirs) if c]
    say(f"  cells with dumps: {len(cells)} of {len(SIGMAS) * len(SEEDS)}")
    say(f"  {'label':14s} rc ok {'ovlp':>6s}  f_SS    f_SL    f_LL     (listed pairs SS/SL/LL)")
    for c in sorted(cells, key=lambda c: (c["seed"], c["sig"])):
        f = c["frac"]
        say(f"  {c['label']:14s} {c['rc']:2d} {int(c['healthy']):2d} {c['overlap']:6.4f}  {f['SS']:.4f}  {f['SL']:.4f}  {f['LL']:.4f}   "
            f"({c['tot']['SS']}/{c['tot']['SL']}/{c['tot']['LL']})")
    by = {(c["seed"], c["sig"]): c for c in cells if c["healthy"]}
    res = {"SS/LL": [], "SL/LL": []}
    for sd in SEEDS:
        if not all((sd, s) in by for s in SIGMAS):
            say(f"  seed {sd}: INCOMPLETE ({sum((sd, s) in by for s in SIGMAS)}/{len(SIGMAS)} healthy)")
            continue
        s = np.array(SIGMAS)
        sh = {}
        for key in ("SS", "SL", "LL"):
            f = np.array([by[(sd, x)]["frac"][key] for x in SIGMAS])
            sh[key] = shalf(s, f / f[-1])
        say(f"  seed {sd}: sigma_half SS {sh['SS']:.3f}  SL {sh['SL']:.3f}  LL {sh['LL']:.3f}")
        if all(np.isfinite(v) for v in sh.values()):
            res["SS/LL"].append(math.log(sh["SS"] / sh["LL"]))
            res["SL/LL"].append(math.log(sh["SL"] / sh["LL"]))
    say("")
    for key, (pd, pc) in (("SS/LL", (math.log(2), math.log(4))), ("SL/LL", (math.log(1.5), math.log(2.25)))):
        v = np.array(res[key])
        if len(v) < 3:
            say(f"  {key}: NOT GRADED ({len(v)}/3 seeds)")
            continue
        m, se = float(v.mean()), float(v.std(ddof=1) / math.sqrt(3))
        say(f"  {key}: ln ratio mean {m:.4f} se {se:.4f} (ratio {math.exp(m):.3f}); per seed "
            + " ".join(f"{x:.3f}" for x in v) + f"; Derjaguin {pd:.3f}, constant force {pc:.3f}; m_eff {m / pd:.2f}")
        if key == "SS/LL":
            if se >= SE_MAX:
                say(f"  B1: NOT INFORMATIVE (se {se:.4f} >= {SE_MAX:.4f})")
            else:
                d = abs(m - pd) <= T_CRIT * se
                c = abs(m - pc) <= T_CRIT * se
                v1 = "BOTH" if d and c else "DERJAGUIN" if d else "CONSTANT-FORCE" if c else "NEITHER"
                say(f"  B1 (registered, 95 %): {v1}  (|m - ln2| = {abs(m - pd) / se:.1f} se, |m - ln4| = {abs(m - pc) / se:.1f} se)")
    post_hoc(by)
    outp = os.path.join(OUT, "grade_lfsig2bi.txt")
    with open(outp, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwrote {outp}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
